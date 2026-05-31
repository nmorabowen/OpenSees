/* ********************************************************************** **
**  MPCO_Ladruno recorder — modular sibling of MPCORecorder (frozen).     **
**  MPCOL_Sinks.cpp — StreamingSink + EnvelopeSink (Phase 2).             **
**                                                                        **
**  All HDF5 *writes* (group create, dataset create+write, attribute      **
**  write) go through mpcol::h5::*. The only raw HDF5 C-API calls are      **
**  *navigation* — testing existence (H5Lexists) and opening an already-   **
**  existing parent group (H5Gopen2) — and deleting a stale envelope       **
**  dataset on in-place rewrite (H5Ldelete). Those are not result writes;  **
**  they resolve where in the tree the wrapper writes.                     **
** ********************************************************************** */

#include "MPCOL_Sinks.h"
#include "MPCOL_Hdf5.h" // pulls MPCOL_Types.h; provides mpcol::h5::*

#include <sstream>
#include <cmath>

namespace mpcol {

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
	/* The recorder (MPCORecorderLadruno::writeModel) created the stage as   */
	/* "MODEL_STAGE[<current_model_stage_id>]" off the file root and an      */
	/* (empty) RESULTS child. The sinks resolve that same path here. We      */
	/* open-or-create each link in the chain so a sink works whether the     */
	/* recorder pre-created RESULTS or not, and whether this is the first    */
	/* result of the stage or a sibling.                                     */
	/* ===================================================================== */

	namespace {

		// MODEL_STAGE group name, matching MPCORecorderLadruno::writeModel exactly.
		std::string stageGroupName(const mpco::ProcessInfo& info)
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
		hid_t openFamilyGroup(mpco::ProcessInfo& info, ResultFamily::Enum family)
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

	void StreamingSink::begin(mpco::ProcessInfo& info, const ResultSource& src)
	{
		const ResultSchema& schema = src.schema();
		if (schema.num_components < 1)
			return;

		hid_t h_family = openFamilyGroup(info, m_family);
		if (h_family == HID_INVALID)
			return;

		// Result group with the self-describing attrs (schema §7.1).
		hid_t h_gp_result = h5::group::createResultGroup(
			h_family, info.h_group_proplist,
			schema.name, schema.display_name,
			schema.components_csv, schema.num_components,
			schema.dimension, schema.description,
			(int)schema.result_type, (int)schema.data_type);

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
		hid_t h_data = h5::dataset::createTimeSeries3d(
			h_gp_result, "DATA", (hsize_t)n_ids, (hsize_t)n_comp);
		hid_t h_time = h5::dataset::createTimeAxis1d(h_gp_result, "TIME", H5T_IEEE_F64LE);
		hid_t h_step = h5::dataset::createTimeAxis1d(h_gp_result, "STEP", H5T_STD_I32LE);

		h5::dataset::close(h_step);
		h5::dataset::close(h_time);
		h5::dataset::close(h_data);
		h5::dataset::close(h_dset_id);
		h5::group::close(h_gp_result);
		h5::group::close(h_family);

		m_initialized = true;
	}

	void StreamingSink::accept(mpco::ProcessInfo& info, const ResultSource& src,
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

		const std::vector<int>& ids = src.ids();
		const size_t n_ids = ids.size();
		const size_t n_comp = (size_t)schema.num_components;

		hid_t h_family = openFamilyGroup(info, m_family);
		if (h_family == HID_INVALID)
			return;

		// Open the (already-created) result group, then append this step's slab
		// to DATA[T x nIds x nComp] and the matching TIME/STEP axes (schema D3).
		hid_t h_gp_result = H5Gopen2(h_family, schema.name.c_str(), H5P_DEFAULT);
		if (h_gp_result < 0) {
			h5::group::close(h_family);
			return;
		}
		hid_t h_data = H5Dopen2(h_gp_result, "DATA", H5P_DEFAULT);
		if (h_data < 0) {
			h5::group::close(h_gp_result);
			h5::group::close(h_family);
			return;
		}
		if (buffer.size() >= n_ids * n_comp)
			h5::dataset::appendSlab3d(h_data, &buffer[0], (hsize_t)n_ids, (hsize_t)n_comp);
		h5::dataset::close(h_data);

		hid_t h_time = H5Dopen2(h_gp_result, "TIME", H5P_DEFAULT);
		if (h_time >= 0) {
			h5::dataset::appendDouble1d(h_time, info.current_time_step);
			h5::dataset::close(h_time);
		}
		hid_t h_step = H5Dopen2(h_gp_result, "STEP", H5P_DEFAULT);
		if (h_step >= 0) {
			h5::dataset::appendInt1d(h_step, info.current_time_step_id);
			h5::dataset::close(h_step);
		}

		h5::group::close(h_gp_result);
		h5::group::close(h_family);
	}

	void StreamingSink::finalize(mpco::ProcessInfo& /*info*/)
	{
		// Streaming has no deferred state: every step is fully written in accept().
		// Handles are opened/closed per call, so there is nothing to flush.
	}

	/* ===================================================================== */
	/* EnvelopeSink                                                          */
	/* ===================================================================== */

	void EnvelopeSink::begin(mpco::ProcessInfo& /*info*/, const ResultSource& src)
	{
		// Cache identity from the source (stable for its lifetime). Accumulators
		// are not allocated here — they are seeded on the first accept() so MIN/MAX
		// begin at the real first sample rather than +/-inf sentinels.
		const ResultSchema& schema = src.schema();
		m_name = schema.name;
		m_n_comp = (size_t)(schema.num_components < 0 ? 0 : schema.num_components);
		m_ids = src.ids();
		m_n_ids = m_ids.size();
	}

	void EnvelopeSink::accept(mpco::ProcessInfo& info, const ResultSource& src,
	                          const std::vector<double>& buffer)
	{
		const ResultSchema& schema = src.schema();
		if (schema.num_components < 1)
			return;

		// Ensure identity is cached even if begin() was skipped.
		if (m_name.empty())
			m_name = schema.name;
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
			if (v < m_min[i]) m_min[i] = v;
			if (v > m_max[i]) m_max[i] = v;
			const double a = std::abs(v);
			if (a > m_absmax[i]) {
				m_absmax[i] = a;
				m_arg_step[i] = step;
			}
		}
	}

	void EnvelopeSink::writeEnvelope(mpco::ProcessInfo& info)
	{
		if (!m_seeded || m_name.empty() || m_n_ids == 0 || m_n_comp == 0)
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

		// Delete a stale <name> group from a prior flush (in-place rewrite), then
		// recreate it fresh. Delete-and-recreate keeps the writer trivial.
		if (H5Lexists(h_family, m_name.c_str(), H5P_DEFAULT) > 0)
			H5Ldelete(h_family, m_name.c_str(), H5P_DEFAULT);

		hid_t h_name = h5::group::create(
			h_family, m_name.c_str(), H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);

		// ID [nIds x 1], and the four [nIds x nComp] accumulators (schema §7.4).
		hid_t d_id  = h5::dataset::createAndWrite(h_name, "ID", m_ids, m_n_ids, 1);
		hid_t d_min = h5::dataset::createAndWrite(h_name, "MIN", m_min, m_n_ids, m_n_comp);
		hid_t d_max = h5::dataset::createAndWrite(h_name, "MAX", m_max, m_n_ids, m_n_comp);
		hid_t d_abs = h5::dataset::createAndWrite(h_name, "ABSMAX", m_absmax, m_n_ids, m_n_comp);
		hid_t d_arg = h5::dataset::createAndWrite(h_name, "ARG_STEP", m_arg_step, m_n_ids, m_n_comp);

		h5::dataset::close(d_id);
		h5::dataset::close(d_min);
		h5::dataset::close(d_max);
		h5::dataset::close(d_abs);
		h5::dataset::close(d_arg);
		h5::group::close(h_name);
		h5::group::close(h_family);
	}

	void EnvelopeSink::flush(mpco::ProcessInfo& info)
	{
		writeEnvelope(info);
	}

	void EnvelopeSink::finalize(mpco::ProcessInfo& info)
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
	}

} // namespace mpcol
