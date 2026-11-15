/* ********************************************************************** **
**  Ladruno recorder — modular sibling of MPCORecorder (frozen).     **
**  Ladruno_Sinks.cpp — StreamingSink + EnvelopeSink (Phase 2).             **
**                                                                        **
**  All HDF5 *writes* (group create, dataset create+write, attribute      **
**  write) go through ladruno::h5::*. The only raw HDF5 C-API calls are      **
**  *navigation* — testing existence (H5Lexists) and opening an already-   **
**  existing parent group (H5Gopen2) — and deleting a stale envelope       **
**  dataset on in-place rewrite (H5Ldelete). Those are not result writes;  **
**  they resolve where in the tree the wrapper writes.                     **
** ********************************************************************** */

// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

#include "Ladruno_Sinks.h"
#include "Ladruno_Hdf5.h" // pulls Ladruno_Types.h; provides ladruno::h5::*

#include "OPS_Globals.h"   // opserr (WP-163 R5 failure reports)

#include <sstream>
#include <cmath>

namespace ladruno {

	/* ===================================================================== */
	/* ResultFamily                                                          */
	/* ===================================================================== */

	const char* ResultFamily::groupName(ResultFamily::Enum f)
	{
		switch (f) {
		case ResultFamily::OnNodes:    return "ON_NODES";
		case ResultFamily::OnElements: return "ON_ELEMENTS";
		case ResultFamily::OnDomain:   return "ON_DOMAIN";
		case ResultFamily::OnRegions:  return "ON_REGIONS";
		default:                       return "ON_NODES";
		}
	}

	/* ===================================================================== */
	/* tree navigation helpers (local to this TU)                            */
	/*                                                                       */
	/* The recorder (LadrunoRecorder::writeModel) created the stage as   */
	/* "MODEL_STAGE[<current_model_stage_id>]" off the file root and an      */
	/* (empty) RESULTS child. The sinks resolve that same path here. We      */
	/* open-or-create each link in the chain so a sink works whether the     */
	/* recorder pre-created RESULTS or not, and whether this is the first    */
	/* result of the stage or a sibling.                                     */
	/* ===================================================================== */

	namespace {

		// MODEL_STAGE group name, matching LadrunoRecorder::writeModel exactly.
		std::string stageGroupName(const detail::ProcessInfo& info)
		{
			std::stringstream ss;
			ss << "MODEL_STAGE[" << info.current_model_stage_id << "]";
			return ss.str();
		}

		// Open child `name` under `loc` if it exists, else create it (creation
		// goes through the wrapper; opening an existing group is plain navigation).
		hid_t openOrCreateGroup(hid_t loc, const char* name, hid_t gplist)
		{
			htri_t exists = H5Lexists(loc, name, H5P_DEFAULT);
			if (exists > 0) {
				return H5Gopen2(loc, name, H5P_DEFAULT);
			}
			return h5::group::create(loc, name, H5P_DEFAULT, gplist, H5P_DEFAULT);
		}

		// Resolve MODEL_STAGE[..]/RESULTS/<family> and return its open handle.
		// Caller closes it. Returns HID_INVALID on failure.
		hid_t openFamilyGroup(detail::ProcessInfo& info, ResultFamily::Enum family)
		{
			hid_t h_stage = openOrCreateGroup(info.h_file_id,
				stageGroupName(info).c_str(), info.h_group_proplist);
			if (h_stage == HID_INVALID || h_stage < 0)
				return HID_INVALID;

			hid_t h_results = openOrCreateGroup(h_stage, "RESULTS", info.h_group_proplist);
			h5::group::close(h_stage);
			if (h_results == HID_INVALID || h_results < 0)
				return HID_INVALID;

			hid_t h_family = openOrCreateGroup(h_results,
				ResultFamily::groupName(family), info.h_group_proplist);
			h5::group::close(h_results);
			return h_family;
		}

	} // anonymous namespace

	/* ===================================================================== */
	/* StreamingSink                                                         */
	/* ===================================================================== */

	void StreamingSink::fail(const ResultSchema& schema, const char* why)
	{
		m_name = schema.name;
		failNamed(why);
	}

	void StreamingSink::failNamed(const char* why)
	{
		opserr << "LadrunoRecorder error: result \"" << m_name.c_str()
		       << "\" " << why << "; it will NOT be recorded from here on\n";
		m_dead = true;
		m_initialized = true;   // never re-run begin() (no per-step retry/spam)
		m_has_pending = false;
		m_pending.clear();
	}

	void StreamingSink::begin(detail::ProcessInfo& info, const ResultSource& src)
	{
		const ResultSchema& schema = src.schema();
		if (schema.num_components < 1)
			return;
		// An empty channel (no ids on this process) has nothing to write. Its
		// [T x 0 x C] DATA could never be created (chunk dim 1 > max dim 0) and
		// used to leave a DATA-less result group behind; write nothing instead.
		if (src.ids().empty()) {
			m_dead = true;
			m_initialized = true;
			return;
		}

		// WP-233: objects are created only on a layout step (SWMR forbids it
		// once armed). The recorder makes a stage's first recorded step a layout
		// step, so reaching here while armed is an internal error: refuse loudly.
		if (!info.layout_open) {
			fail(schema, "would be created while the file is in SWMR-write mode "
			             "(internal error: not on a layout step)");
			return;
		}

		hid_t h_family = openFamilyGroup(info, m_family);
		if (h_family == HID_INVALID || h_family < 0) {
			fail(schema, "could not open its RESULTS family group");
			return;
		}

		// WP-163 R5/R4: result group names are unique per channel within a
		// MODEL_STAGE (element buckets embed classTag + rule + header index), so
		// a pre-existing group means this result was requested twice (e.g.
		// `-N displacement displacement`, or the alias pair tieForce /
		// constraintTieForce) or two recorders target one file. Appending into
		// the existing group doubled its rows silently; refuse instead.
		if (H5Lexists(h_family, schema.name.c_str(), H5P_DEFAULT) > 0) {
			h5::group::close(h_family);
			fail(schema, "already exists in this MODEL_STAGE (duplicate request, "
			             "or two recorders writing the same file)");
			return;
		}

		// Result group with the self-describing attrs (schema §7.1).
		hid_t h_gp_result = h5::group::createResultGroup(
			h_family, info.h_group_proplist,
			schema.name, schema.display_name,
			schema.components_csv, schema.num_components,
			schema.dimension, schema.description,
			(int)schema.result_type, (int)schema.data_type);
		if (h_gp_result < 0) {
			h5::group::close(h_family);
			fail(schema, "could not create its result group");
			return;
		}
		// WP-126: how a reader combines this result across partition files (§7.1).
		h5::attribute::write(h_gp_result, "PARTITION_REDUCTION",
			std::string(src.partitionReduction()));

		// ID dataset [nIds x 1] (same shape as the frozen recorder).
		const std::vector<int>& ids = src.ids();
		const size_t n_ids = ids.size();
		const size_t n_comp = (size_t)schema.num_components;
		hid_t h_dset_id = h5::dataset::createAndWrite(
			h_gp_result, "ID", ids, n_ids, 1);

		// Chunked time-series layout (schema D3): DATA is now a single
		// [T x nIds x nComp] extensible dataset that grows one slab per step,
		// with matching TIME[T] (double) / STEP[T] (int) axes — replacing the
		// old DATA group of per-step STEP_<k> datasets. accept() appends.
		// On-disk DATA type is f64 (default, lossless) or f32 (opt-in
		// `-precision f32` lossy mode); TIME/STEP stay f64/int either way.
		hid_t data_disk_type = info.store_data_f32 ? H5T_IEEE_F32LE : H5T_IEEE_F64LE;
		// WP-233: remember where and how, to reopen after a layout step.
		{
			std::stringstream ss_path;
			ss_path << stageGroupName(info) << "/RESULTS/"
			        << ResultFamily::groupName(m_family) << "/" << schema.name;
			m_path = ss_path.str();
		}
		m_name = schema.name;
		m_disk_type = data_disk_type;
		m_deflate = info.deflate_level;
		// WP-164/WP-233: the SWMR-safe chunk plan + sized chunk cache, deflate
		// level from `-compress` (0 = unfiltered 64 KiB chunks; > 0 = one deflated
		// step per chunk), and the handles are KEPT OPEN until the sink is
		// destroyed or a layout step releases them.
		m_data = h5::dataset::createTimeSeries3dOpen(
			h_gp_result, "DATA", (hsize_t)n_ids, (hsize_t)n_comp, data_disk_type,
			info.deflate_level);
		m_time = h5::dataset::createTimeAxis1d(h_gp_result, "TIME", H5T_IEEE_F64LE);
		m_step = h5::dataset::createTimeAxis1d(h_gp_result, "STEP", H5T_STD_I32LE);
		hsize_t mdims[3] = { 1, (hsize_t)n_ids, (hsize_t)n_comp };
		m_mspace = H5Screate_simple(3, mdims, NULL);
		hsize_t mdims1[1] = { 1 };
		m_mspace1 = H5Screate_simple(1, mdims1, NULL);
		m_n_ids = n_ids;
		m_n_comp = n_comp;
		m_t = 0;
		const bool ok = (h_dset_id >= 0 && m_data >= 0 && m_time >= 0 && m_step >= 0 &&
		                 m_mspace >= 0 && m_mspace1 >= 0);

		if (h_dset_id >= 0) h5::dataset::close(h_dset_id);
		h5::group::close(h_gp_result);
		h5::group::close(h_family);

		if (!ok) {
			closeHandles();
			// WP-163 R5: previously ignored — every later accept() then failed
			// H5Dopen2("DATA") and returned silently, dropping the whole result.
			fail(schema, "could not create its ID/DATA/TIME/STEP datasets "
			             "(disk full or quota, or HDF5 refused the dataset)");
			return;
		}
		m_initialized = true;
	}

	void StreamingSink::accept(detail::ProcessInfo& info, const ResultSource& src,
	                           const std::vector<double>& buffer)
	{
		const ResultSchema& schema = src.schema();
		if (schema.num_components < 1)
			return;

		// Defensive: if accept() is reached without a begin() for this stage
		// (e.g. a re-begin path missed it), create the group lazily so a step
		// is never silently dropped.
		if (!m_initialized)
			begin(info, src);
		if (m_dead)
			return;

		const std::vector<int>& ids = src.ids();
		const size_t n_ids = ids.size();
		const size_t n_comp = (size_t)schema.num_components;

		// WP-163 R5: a short buffer skips the WHOLE step (DATA and TIME/STEP
		// together) so the three axes stay aligned; previously DATA was skipped
		// while TIME/STEP were still appended, shifting every later row.
		if (buffer.size() < n_ids * n_comp) {
			if (!m_warned_short) {
				opserr << "LadrunoRecorder warning: result \"" << schema.name.c_str()
				       << "\" got " << (int)buffer.size() << " values, expected "
				       << (int)(n_ids * n_comp) << "; step(s) skipped\n";
				m_warned_short = true;
			}
			return;
		}

		const double t_val = info.current_time_step;
		const int s_val = info.current_time_step_id;
		// WP-233: on a layout step the file is not in SWMR mode yet; hold the row
		// and append it in commitPending(), after the recorder arms SWMR, so no
		// streaming row is ever written outside SWMR.
		if (info.layout_open) {
			m_pending.assign(buffer.begin(), buffer.begin() + n_ids * n_comp);
			m_pending_time = t_val;
			m_pending_step = s_val;
			m_has_pending = true;
			return;
		}
		writeRow(&buffer[0], t_val, s_val);
	}

	void StreamingSink::writeRow(const double* row, double t_val, int s_val)
	{
		// WP-164 (P2): write into the handles held open since begin(). The slab
		// lands in the dataset's chunk cache and reaches the file at the next
		// flush or eviction. Under SWMR (WP-233) each dataset's own writes reach
		// the disk in order, but DATA, TIME and STEP are separate datasets: a kill
		// can leave them at different lengths, and a reader cuts to the shortest
		// (schema §7b.2).
		herr_t st = h5::dataset::writeSlab3dAt(m_data, m_t, row,
			(hsize_t)m_n_ids, (hsize_t)m_n_comp, m_mspace);
		if (st < 0) {
			// WP-163 R5: if DATA was extended before the write failed, that row
			// exists (fill value) — give it its TIME/STEP so the axes stay aligned,
			// then stop the channel.
			const bool grew = h5::dataset::extent0(m_data) > m_t;
			failNamed("failed to append a step to DATA (disk full or quota?)");
			if (!grew)
				return;
		}
		h5::dataset::writeScalar1dAt(m_time, m_t, H5T_NATIVE_DOUBLE, &t_val, m_mspace1);
		h5::dataset::writeScalar1dAt(m_step, m_t, H5T_NATIVE_INT, &s_val, m_mspace1);
		++m_t;
	}

	void StreamingSink::commitPending(detail::ProcessInfo& /*info*/)
	{
		if (!m_has_pending)
			return;
		m_has_pending = false;
		if (!m_dead && m_data >= 0)
			writeRow(&m_pending[0], m_pending_time, m_pending_step);
		m_pending.clear();
	}

	void StreamingSink::releaseHandles()
	{
		// WP-233: the datasets only; the memory spaces are not file objects.
		if (m_step >= 0) H5Dclose(m_step);
		if (m_time >= 0) H5Dclose(m_time);
		if (m_data >= 0) H5Dclose(m_data);
		m_data = m_time = m_step = HID_INVALID;
	}

	bool StreamingSink::reacquireHandles(detail::ProcessInfo& info)
	{
		if (m_dead || !m_initialized || m_path.empty())
			return true;
		hid_t h_grp = H5Gopen2(info.h_file_id, m_path.c_str(), H5P_DEFAULT);
		if (h_grp >= 0) {
			hid_t dapl = h5::dataset::timeSeriesAccessList(
				(hsize_t)m_n_ids, (hsize_t)m_n_comp, m_disk_type, m_deflate);
			m_data = H5Dopen2(h_grp, "DATA", dapl >= 0 ? dapl : H5P_DEFAULT);
			if (dapl >= 0) H5Pclose(dapl);
			m_time = H5Dopen2(h_grp, "TIME", H5P_DEFAULT);
			m_step = H5Dopen2(h_grp, "STEP", H5P_DEFAULT);
			H5Gclose(h_grp);
		}
		if (h_grp < 0 || m_data < 0 || m_time < 0 || m_step < 0) {
			releaseHandles();
			failNamed("could not reopen its DATA/TIME/STEP after a layout step");
			return false;
		}
		return true;
	}

	void StreamingSink::closeHandles()
	{
		if (m_mspace1 >= 0) H5Sclose(m_mspace1);
		if (m_mspace >= 0) H5Sclose(m_mspace);
		if (m_step >= 0) H5Dclose(m_step);
		if (m_time >= 0) H5Dclose(m_time);
		if (m_data >= 0) H5Dclose(m_data);
		m_data = m_time = m_step = m_mspace = m_mspace1 = HID_INVALID;
	}

	StreamingSink::~StreamingSink()
	{
		closeHandles();
	}

	void StreamingSink::finalize(detail::ProcessInfo& /*info*/)
	{
		// Every step is handed to HDF5 in accept(); partially filled chunks sit in
		// the open dataset's chunk cache until the recorder's H5Fflush (its flush
		// cadence) or the handles close in the destructor. Nothing to do here.
	}

	/* ===================================================================== */
	/* EnvelopeSink                                                          */
	/* ===================================================================== */

	void EnvelopeSink::begin(detail::ProcessInfo& /*info*/, const ResultSource& src)
	{
		// Cache identity from the source (stable for its lifetime). Accumulators
		// are not allocated here — they are seeded on the first accept() so MIN/MAX
		// begin at the real first sample rather than +/-inf sentinels.
		const ResultSchema& schema = src.schema();
		m_name = schema.name;
		m_display_name = schema.display_name;
		m_components_csv = schema.components_csv;
		m_dimension = schema.dimension;
		m_description = schema.description;
		m_result_type = (int)schema.result_type;
		m_data_type = (int)schema.data_type;
		m_partition_reduction = src.partitionReduction();   // WP-126
		m_n_comp = (size_t)(schema.num_components < 0 ? 0 : schema.num_components);
		m_ids = src.ids();
		m_n_ids = m_ids.size();
	}

	void EnvelopeSink::accept(detail::ProcessInfo& info, const ResultSource& src,
	                          const std::vector<double>& buffer)
	{
		const ResultSchema& schema = src.schema();
		if (schema.num_components < 1)
			return;

		// Ensure identity + self-describing metadata are cached even if begin() was
		// skipped (the recorder drives accept() directly, no begin() call).
		if (m_name.empty()) {
			m_name = schema.name;
			m_display_name = schema.display_name;
			m_components_csv = schema.components_csv;
			m_dimension = schema.dimension;
			m_description = schema.description;
			m_result_type = (int)schema.result_type;
			m_data_type = (int)schema.data_type;
			m_partition_reduction = src.partitionReduction();   // WP-126 (begin() is skipped here)
		}
		m_n_comp = (size_t)schema.num_components;
		if (m_ids.empty()) {
			m_ids = src.ids();
			m_n_ids = m_ids.size();
		}

		const size_t n = m_n_ids * m_n_comp;
		if (n == 0 || buffer.size() < n)
			return;

		const int step = info.current_time_step_id;

		if (!m_seeded) {
			// First accept: seed every accumulator from this step's values.
			m_min.assign(buffer.begin(), buffer.begin() + n);
			m_max.assign(buffer.begin(), buffer.begin() + n);
			m_absmax.resize(n);
			m_arg_step.assign(n, step);
			for (size_t i = 0; i < n; ++i)
				m_absmax[i] = std::abs(buffer[i]);
			m_seeded = true;
			return;
		}

		// Subsequent accepts: element-wise componentwise update (ADR D7 — each
		// component tracked independently; ARG_STEP records the step at which the
		// absolute extreme for THAT component was last attained).
		for (size_t i = 0; i < n; ++i) {
			const double v = buffer[i];
			// WP-163 (ROB-7): NaN is STICKY. Ordered comparisons with NaN are
			// false, so a run that went NaN at step k used to keep its finite
			// pre-divergence extremes and the envelope looked healthy. Now the
			// first NaN poisons MIN/MAX/ABSMAX of that component and ARG_STEP
			// records the step it appeared (a NaN first sample already stuck).
			if (v != v) {
				if (m_absmax[i] == m_absmax[i]) {
					m_min[i] = v; m_max[i] = v; m_absmax[i] = v;
					m_arg_step[i] = step;
				}
				continue;
			}
			if (v < m_min[i]) m_min[i] = v;
			if (v > m_max[i]) m_max[i] = v;
			const double a = std::abs(v);
			if (a > m_absmax[i]) {
				m_absmax[i] = a;
				m_arg_step[i] = step;
			}
		}
	}

	void EnvelopeSink::writeEnvelope(detail::ProcessInfo& info)
	{
		if (!m_seeded || m_name.empty() || m_n_ids == 0 || m_n_comp == 0)
			return;

		// WP-164 (P1): after the first write, overwrite the four accumulator
		// datasets in place — no group walk, no delete/recreate, no attributes.
		if (m_created) {
			if (m_dmin >= 0) H5Dwrite(m_dmin, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_min[0]);
			if (m_dmax >= 0) H5Dwrite(m_dmax, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_max[0]);
			if (m_dabs >= 0) H5Dwrite(m_dabs, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_absmax[0]);
			if (m_darg >= 0) H5Dwrite(m_darg, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &m_arg_step[0]);
			return;
		}

		// WP-233: the group is created on a layout step only (SWMR forbids it
		// once armed). An envelope first seeded on an armed step waits: the
		// recorder sees needsLayout() and makes the next recorded step a layout
		// step, which creates it holding the extremes so far.
		if (!info.layout_open)
			return;

		// ENVELOPES/<family>/<name> lives under the stage's RESULTS tree.
		hid_t h_stage = openOrCreateGroup(info.h_file_id,
			stageGroupName(info).c_str(), info.h_group_proplist);
		if (h_stage == HID_INVALID || h_stage < 0)
			return;
		hid_t h_results = openOrCreateGroup(h_stage, "RESULTS", info.h_group_proplist);
		h5::group::close(h_stage);
		if (h_results == HID_INVALID || h_results < 0)
			return;
		hid_t h_env = openOrCreateGroup(h_results, "ENVELOPES", info.h_group_proplist);
		h5::group::close(h_results);
		if (h_env == HID_INVALID || h_env < 0)
			return;
		hid_t h_family = openOrCreateGroup(h_env,
			ResultFamily::groupName(m_family), info.h_group_proplist);
		h5::group::close(h_env);
		if (h_family == HID_INVALID || h_family < 0)
			return;

		// Element result names are nested ("<display>/<bucket>", e.g.
		// "stress/204-FourNodeQuad[201:0:0]"); node/domain names are flat. The group
		// proplist carries no create-intermediate-group link property, so H5Gcreate
		// (in createResultGroup below) will NOT auto-create the "<display>" parent —
		// pre-create each ancestor prefix once. The delete-recreate just below only
		// removes the leaf link, so these intermediates persist across flushes.
		for (size_t slash = m_name.find('/'); slash != std::string::npos;
		     slash = m_name.find('/', slash + 1)) {
			std::string prefix = m_name.substr(0, slash);
			if (H5Lexists(h_family, prefix.c_str(), H5P_DEFAULT) <= 0) {
				hid_t h_mid = h5::group::create(h_family, prefix.c_str(),
					H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
				if (h_mid >= 0)
					h5::group::close(h_mid);
			}
		}

		// WP-164 (P1): first write of this stage — create the group once. (A
		// pre-existing group can only be a duplicate request, which the parser
		// now drops; it is not deleted and recreated here any more.)
		if (H5Lexists(h_family, m_name.c_str(), H5P_DEFAULT) > 0) {
			opserr << "LadrunoRecorder error: envelope \"" << m_name.c_str()
			       << "\" already exists in this MODEL_STAGE; it will NOT be recorded\n";
			m_created = true;   // never retry (handles stay invalid)
			h5::group::close(h_family);
			return;
		}

		// Self-describing result group (same COMPONENTS/DISPLAY_NAME/DIMENSION attrs
		// as the time-series StreamingSink writes) so envelope output carries
		// authoritative component names in-file (finding (h)), not just bare arrays.
		hid_t h_name = h5::group::createResultGroup(
			h_family, info.h_group_proplist, m_name, m_display_name,
			m_components_csv, (int)m_n_comp, m_dimension, m_description,
			m_result_type, m_data_type);
		// WP-126: same attribute as the streaming group. In a partitioned run the
		// recorder refuses to envelope anything but "NONE", so a written envelope is
		// always "NONE" there; serial files keep the source's own value.
		h5::attribute::write(h_name, "PARTITION_REDUCTION", m_partition_reduction);

		// ID [nIds x 1], and the four [nIds x nComp] accumulators (schema §7.4).
		// WP-164: the four accumulator handles are kept for in-place rewrites.
		hid_t d_id = h5::dataset::createAndWrite(h_name, "ID", m_ids, m_n_ids, 1);
		m_dmin = h5::dataset::createAndWrite(h_name, "MIN", m_min, m_n_ids, m_n_comp);
		m_dmax = h5::dataset::createAndWrite(h_name, "MAX", m_max, m_n_ids, m_n_comp);
		m_dabs = h5::dataset::createAndWrite(h_name, "ABSMAX", m_absmax, m_n_ids, m_n_comp);
		m_darg = h5::dataset::createAndWrite(h_name, "ARG_STEP", m_arg_step, m_n_ids, m_n_comp);
		m_created = true;
		{
			std::stringstream ss_path;
			ss_path << stageGroupName(info) << "/RESULTS/ENVELOPES/"
			        << ResultFamily::groupName(m_family) << "/" << m_name;
			m_path = ss_path.str();
		}

		h5::dataset::close(d_id);
		h5::group::close(h_name);
		h5::group::close(h_family);
	}

	void EnvelopeSink::releaseHandles()
	{
		if (m_darg >= 0) H5Dclose(m_darg);
		if (m_dabs >= 0) H5Dclose(m_dabs);
		if (m_dmax >= 0) H5Dclose(m_dmax);
		if (m_dmin >= 0) H5Dclose(m_dmin);
		m_dmin = m_dmax = m_dabs = m_darg = HID_INVALID;
	}

	bool EnvelopeSink::reacquireHandles(detail::ProcessInfo& info)
	{
		if (!m_created || m_path.empty())
			return true;
		hid_t h_grp = H5Gopen2(info.h_file_id, m_path.c_str(), H5P_DEFAULT);
		if (h_grp >= 0) {
			m_dmin = H5Dopen2(h_grp, "MIN", H5P_DEFAULT);
			m_dmax = H5Dopen2(h_grp, "MAX", H5P_DEFAULT);
			m_dabs = H5Dopen2(h_grp, "ABSMAX", H5P_DEFAULT);
			m_darg = H5Dopen2(h_grp, "ARG_STEP", H5P_DEFAULT);
			H5Gclose(h_grp);
		}
		if (h_grp < 0 || m_dmin < 0 || m_dmax < 0 || m_dabs < 0 || m_darg < 0) {
			releaseHandles();
			opserr << "LadrunoRecorder error: envelope \"" << m_name.c_str()
			       << "\" could not be reopened after a layout step; its last "
			          "written extremes stay in the file\n";
			m_path.clear();
			return false;
		}
		return true;
	}

	EnvelopeSink::~EnvelopeSink()
	{
		if (m_darg >= 0) H5Dclose(m_darg);
		if (m_dabs >= 0) H5Dclose(m_dabs);
		if (m_dmax >= 0) H5Dclose(m_dmax);
		if (m_dmin >= 0) H5Dclose(m_dmin);
	}

	void EnvelopeSink::flush(detail::ProcessInfo& info)
	{
		writeEnvelope(info);
	}

	void EnvelopeSink::finalize(detail::ProcessInfo& info)
	{
		writeEnvelope(info);
	}

	void EnvelopeSink::reset()
	{
		m_seeded = false;
		m_n_ids = 0;
		m_n_comp = 0;
		m_name.clear();
		m_ids.clear();
		m_min.clear();
		m_max.clear();
		m_absmax.clear();
		m_arg_step.clear();
		if (m_darg >= 0) H5Dclose(m_darg);
		if (m_dabs >= 0) H5Dclose(m_dabs);
		if (m_dmax >= 0) H5Dclose(m_dmax);
		if (m_dmin >= 0) H5Dclose(m_dmin);
		m_dmin = m_dmax = m_dabs = m_darg = HID_INVALID;
		m_created = false;
	}

} // namespace ladruno
