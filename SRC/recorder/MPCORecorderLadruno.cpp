/* ********************************************************************** **
**  MPCO_Ladruno recorder — Phase-1 skeleton implementation.              **
**                                                                        **
**  Goal of Phase 1 (ADR 03_mpco_ladruno, execution plan):               **
**    - prove the `namespace mpcol` modules build AND link next to the    **
**      frozen MPCORecorder translation unit (the ODR-safety gate);       **
**    - register `recorder mpcoLadruno`;                                  **
**    - emit a schema-v1-valid .ladruno file (INFO + one MODEL_STAGE      **
**      with MODEL/NODES) that passes ladruno_recorder_tests/validate.    **
**                                                                        **
**  Results (node/element/domain), envelopes, parallel part-files and the **
**  ResultSource/ResultSink wiring arrive in Phases 2-3.                  **
** ********************************************************************** */

#include "MPCORecorderLadruno.h"

// mpcol modules (all symbols namespaced to avoid ODR clash with the frozen file)
#include "MPCOL_Hdf5.h"     // pulls MPCOL_Types.h
#include "MPCOL_ResultIO.h"

// OpenSees
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <classTags.h>

#include <string>
#include <vector>
#include <sstream>

namespace mpcolns = mpcol; // local alias

/* ===================================================================== */
/* private_data                                                          */
/* ===================================================================== */

class MPCORecorderLadruno::private_data
{
public:
	private_data()
		: filename()
		, initialized(false)
		, info()
	{}

	std::string filename;
	bool initialized;
	mpcolns::mpco::ProcessInfo info;
};

/* ===================================================================== */
/* ctor / dtor                                                           */
/* ===================================================================== */

MPCORecorderLadruno::MPCORecorderLadruno()
	: Recorder(RECORDER_TAGS_MPCOLadrunoRecorder)
	, m_data(new private_data())
{
}

MPCORecorderLadruno::~MPCORecorderLadruno()
{
	if (m_data) {
		if (m_data->initialized && m_data->info.h_file_id != mpcolns::HID_INVALID) {
			mpcolns::h5::file::flush(m_data->info.h_file_id);
			mpcolns::h5::file::close(m_data->info.h_file_id);
			mpcolns::h5::plist::close(m_data->info.h_group_proplist);
			mpcolns::h5::plist::close(m_data->info.h_file_proplist);
		}
		delete m_data;
	}
}

/* ===================================================================== */
/* lifecycle                                                             */
/* ===================================================================== */

int MPCORecorderLadruno::setDomain(Domain& theDomain)
{
	m_data->info.domain = &theDomain;
	return 0;
}

int MPCORecorderLadruno::restart(void)
{
	return 0;
}

int MPCORecorderLadruno::domainChanged(void)
{
	return 0;
}

int MPCORecorderLadruno::record(int commitTag, double timeStamp)
{
	m_data->info.current_time_step_id = commitTag;
	m_data->info.current_time_step = timeStamp;

	if (!m_data->initialized) {
		if (initialize() != 0)
			return -1;
		m_data->initialized = true;
	}

	// Phase 1 writes the model once on first record. Result streaming arrives in
	// Phase 3; for now keep the file flushed so it is always readable.
	if (m_data->info.h_file_id != mpcolns::HID_INVALID)
		mpcolns::h5::file::flush(m_data->info.h_file_id);

	return 0;
}

/* ===================================================================== */
/* file creation                                                         */
/* ===================================================================== */

int MPCORecorderLadruno::initialize()
{
	mpcolns::mpco::ProcessInfo& info = m_data->info;

	if (info.domain == 0) {
		opserr << "MPCORecorderLadruno error: domain is not defined\n";
		return -1;
	}

	// spatial dimension from the first node (enforce >0 and <=3)
	if (info.num_dimensions < 1) {
		NodeIter& nit = info.domain->getNodes();
		Node* nptr = 0;
		while ((nptr = nit()) != 0) {
			info.num_dimensions = nptr->getCrds().Size();
			break;
		}
		if (info.num_dimensions < 1 || info.num_dimensions > 3) {
			opserr << "MPCORecorderLadruno error: invalid spatial dimension "
			       << info.num_dimensions << "\n";
			return -1;
		}
	}

	// property lists (link creation order tracked/indexed, as the frozen recorder)
	info.h_file_proplist = mpcolns::h5::plist::crate(mpcolns::h5::plist::FileCreate);
	mpcolns::h5::plist::setLinkCreationOrder(
		info.h_file_proplist, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
	info.h_group_proplist = mpcolns::h5::plist::crate(mpcolns::h5::plist::GroupCreate);
	mpcolns::h5::plist::setLinkCreationOrder(
		info.h_group_proplist, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);

	// create the file (Phase 1: no SWMR; H5P_DEFAULT file access)
	info.h_file_id = mpcolns::h5::file::create(
		m_data->filename.c_str(), info.h_file_proplist, H5P_DEFAULT);
	if (info.h_file_id == mpcolns::HID_INVALID) {
		opserr << "MPCORecorderLadruno error: cannot create file \""
		       << m_data->filename.c_str() << "\"\n";
		return -1;
	}

	// INFO group (schema v1 §2)
	hid_t h_info = mpcolns::h5::group::create(
		info.h_file_id, "INFO", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	mpcolns::h5::attribute::write(h_info, "GENERATOR", std::string("MPCO_Ladruno"));
	mpcolns::h5::attribute::write(h_info, "FORMAT_VERSION", (int)1);
	mpcolns::h5::attribute::write(h_info, "SOLVER_NAME", std::string("OpenSees"));
	std::vector<int> solver_version;
	solver_version.push_back(3);
	solver_version.push_back(5);
	solver_version.push_back(1);
	mpcolns::h5::attribute::write(h_info, "SOLVER_VERSION", solver_version);
	mpcolns::h5::attribute::write(h_info, "SPATIAL_DIM", info.num_dimensions);
	// PROVENANCE group (lineage slots; empty in Phase 1)
	hid_t h_prov = mpcolns::h5::group::create(
		h_info, "PROVENANCE", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	mpcolns::h5::group::close(h_prov);
	mpcolns::h5::group::close(h_info);

	return writeModel();
}

int MPCORecorderLadruno::writeModel()
{
	mpcolns::mpco::ProcessInfo& info = m_data->info;
	const int ndim = info.num_dimensions;

	info.current_model_stage_id = info.domain->hasDomainChanged();

	// collect nodes
	std::vector<int> node_ids;
	std::vector<double> node_coords; // row-major [nNodes x ndim]
	{
		NodeIter& nit = info.domain->getNodes();
		Node* nptr = 0;
		while ((nptr = nit()) != 0) {
			node_ids.push_back(nptr->getTag());
			const Vector& crds = nptr->getCrds();
			for (int d = 0; d < ndim; ++d)
				node_coords.push_back(d < crds.Size() ? crds(d) : 0.0);
		}
	}
	const size_t nnodes = node_ids.size();

	// MODEL_STAGE[<stamp>]
	std::stringstream ss_stage;
	ss_stage << "MODEL_STAGE[" << info.current_model_stage_id << "]";
	hid_t h_stage = mpcolns::h5::group::create(
		info.h_file_id, ss_stage.str().c_str(), H5P_DEFAULT,
		info.h_group_proplist, H5P_DEFAULT);
	mpcolns::h5::attribute::write(h_stage, "STEP", info.current_time_step_id);
	mpcolns::h5::attribute::write(h_stage, "TIME", info.current_time_step);
	mpcolns::h5::attribute::write(h_stage, "KIND", std::string("static"));

	// MODEL
	hid_t h_model = mpcolns::h5::group::create(
		h_stage, "MODEL", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);

	// MODEL/NODES
	hid_t h_nodes = mpcolns::h5::group::create(
		h_model, "NODES", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	hid_t d_id = mpcolns::h5::dataset::createAndWrite(h_nodes, "ID", node_ids);
	mpcolns::h5::dataset::close(d_id);
	hid_t d_coord = mpcolns::h5::dataset::createAndWrite(
		h_nodes, "COORDINATES", node_coords, nnodes, (size_t)ndim);
	mpcolns::h5::dataset::close(d_coord);
	mpcolns::h5::group::close(h_nodes);

	// MODEL/ELEMENTS (empty in Phase 1; element groups arrive in Phase 2/3)
	hid_t h_elems = mpcolns::h5::group::create(
		h_model, "ELEMENTS", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	mpcolns::h5::group::close(h_elems);

	mpcolns::h5::group::close(h_model);

	// RESULTS (empty in Phase 1)
	hid_t h_results = mpcolns::h5::group::create(
		h_stage, "RESULTS", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	mpcolns::h5::group::close(h_results);

	mpcolns::h5::group::close(h_stage);
	return 0;
}

/* ===================================================================== */
/* parallel (Phase-1 stubs; real send/recv in Phase 3)                   */
/* ===================================================================== */

int MPCORecorderLadruno::sendSelf(int /*commitTag*/, Channel& /*theChannel*/)
{
	return 0;
}

int MPCORecorderLadruno::recvSelf(int /*commitTag*/, Channel& /*theChannel*/,
                                  FEM_ObjectBroker& /*theBroker*/)
{
	return 0;
}

/* ===================================================================== */
/* command parser  (recorder mpcoLadruno <file> ...)                     */
/* ===================================================================== */

void* OPS_MPCOLadrunoRecorder()
{
	int num_remaining = OPS_GetNumRemainingInputArgs();
	if (num_remaining < 1) {
		opserr << "MPCORecorderLadruno error: insufficient args; expected a filename\n";
		return 0;
	}

	const char* filename = OPS_GetString();
	num_remaining--;

	Domain* domain = OPS_GetDomain();
	if (domain == 0) {
		opserr << "MPCORecorderLadruno error: domain is not defined\n";
		return 0;
	}

	// Phase 1: consume and ignore any remaining options (-N/-E/-T/-R parsing
	// arrives with the result engine in Phase 3).
	while (num_remaining > 0) {
		OPS_GetString();
		num_remaining--;
	}

	MPCORecorderLadruno* recorder = new MPCORecorderLadruno();
	recorder->m_data->filename = filename;
	return recorder;
}
