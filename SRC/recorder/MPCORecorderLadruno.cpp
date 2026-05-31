/* ********************************************************************** **
**  MPCO_Ladruno recorder — Phase-3 integration.                          **
**                                                                        **
**  Wires the namespace-mpcol Source/Sink modules into a working          **
**  recorder. Each -N / -NS request becomes a NodeResultSource +          **
**  StreamingSink(OnNodes); each -E request becomes a set of              **
**  ElementResultSource + StreamingSink(OnElements) buckets, discovered    **
**  via the ported element engine (ElementCollection / OutputDescriptor). **
**  record() drives sink.accept(source.evaluate()) per step, honouring    **
**  the frozen subtleties: reaction-flag grouping, eigen/mode detection,  **
**  and the -T output-frequency gate. The Source/Sink split changes the   **
**  STRUCTURE, not the math (1e-12 parity gate vs frozen `recorder mpco`). **
** ********************************************************************** */

#include "MPCORecorderLadruno.h"

// mpcol modules (all symbols namespaced to avoid ODR clash with the frozen file)
#include "MPCOL_Hdf5.h"     // pulls MPCOL_Types.h
#include "MPCOL_ResultIO.h"
#include "MPCOL_NodeResults.h"
#include "MPCOL_ElementResults.h"
#include "MPCOL_DomainResults.h"
#include "MPCOL_Sinks.h"

// OpenSees
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <classTags.h>
#include <Response.h>
#include <CompositeResponse.h>
#include <MeshRegion.h>
#include <Pressure_ConstraintIter.h>
#include <Pressure_Constraint.h>
#include "section/SectionForceDeformation.h"

#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <cmath>
#include <cstring>

namespace mpcolns = mpcol; // local alias

// max number of iterations to guess the number of fibers (frozen MPCO_MAX_TRIAL_NFIB)
#ifndef MPCO_LADRUNO_MAX_TRIAL_NFIB
#define MPCO_LADRUNO_MAX_TRIAL_NFIB 100000
#endif

/* ===================================================================== */
/* private_data                                                          */
/* ===================================================================== */

class MPCORecorderLadruno::private_data
{
public:
	private_data()
		: filename()
		, initialized(false)
		, first_domain_changed_done(false)
		, info()
		, output_freq()
		, has_region(false)
		, region_tag(0)
		, node_set()
		, elem_set()
		, energy_requested(false)
		, energy_region_tags()
		, nodal_results_requests()
		, sens_grad_indices()
		, elemental_results_requests()
		, nodes()
		, elements()
		, elem_ngauss_nfiber_info()
		, elemental_responses()
	{}

	std::string filename;
	bool initialized;
	// first_domain_changed_done — false until the first record() resolves the
	// model-stage stamp; thereafter a stamp change triggers a model rebuild
	// (frozen MPCORecorder::record() multi-stage block).
	bool first_domain_changed_done;
	mpcolns::mpco::ProcessInfo info;

	// -T output frequency
	mpcolns::mpco::OutputFrequency output_freq;

	// -R region restriction
	bool has_region;
	int region_tag;              // the -R region tag (0 if none); for MODEL/SETS
	std::vector<int> node_set;
	std::vector<int> elem_set;

	// -G energy balance (ADR D8): whole-model and/or per-region energy.
	bool energy_requested;               // -G energy   -> whole-model ON_DOMAIN
	std::vector<int> energy_region_tags; // -G <tag...> -> per-region ON_REGIONS

	// parsed requests
	std::vector<mpcolns::mpco::NodalResultType::Enum> nodal_results_requests;
	std::vector<int> sens_grad_indices;
	std::vector<std::vector<std::string> > elemental_results_requests;

	// the node set (Node* in write order) used by all node sources
	std::vector<Node*> nodes;

	// element engine
	mpcolns::mpco::element::ElementCollection elements;
	// for each element tag: per gauss point (fiber_base_index, num_fibers)
	std::map<int, std::vector<std::pair<int, int> > > elem_ngauss_nfiber_info;
	// all element responses, owned here (released in clearSources)
	std::vector<Response*> elemental_responses;

	/* --- node sources/sinks --- */
	struct NodeChannel {
		mpcolns::ResultSource* source;   // owned
		mpcolns::StreamingSink* sink;    // owned
		int reaction_flag;               // -1 if not a reaction source
		bool is_modes;                   // ModesOfVibration(Rotational) special path
		NodeChannel() : source(0), sink(0), reaction_flag(-1), is_modes(false) {}
	};
	std::vector<NodeChannel> node_channels;

	/* --- element sources/sinks --- */
	struct ElemChannel {
		mpcolns::ElementResultSource* source;             // owned
		mpcolns::StreamingSink* sink;                     // owned
		const mpcolns::mpco::element::OutputDescriptorHeader* header; // not owned
		bool column_map_written;
		ElemChannel() : source(0), sink(0), header(0), column_map_written(false) {}
	};
	std::vector<ElemChannel> elem_channels;
	// the per-request response-collection trees (own the OutputResponseCollection
	// buckets referenced by ElementResultSource); kept alive for source lifetime.
	mpcolns::mpco::element::ResultRecorderCollection elemental_recorders;

	/* --- domain/region sources/sinks (ADR D8 energy balance) --- */
	struct DomainChannel {
		mpcolns::ResultSource* source;   // owned
		mpcolns::StreamingSink* sink;    // owned
		DomainChannel() : source(0), sink(0) {}
	};
	std::vector<DomainChannel> domain_channels;
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
		clearSources();
		if (m_data->initialized && m_data->info.h_file_id != mpcolns::HID_INVALID) {
			mpcolns::h5::file::flush(m_data->info.h_file_id);
			mpcolns::h5::file::close(m_data->info.h_file_id);
			mpcolns::h5::plist::close(m_data->info.h_group_proplist);
			mpcolns::h5::plist::close(m_data->info.h_file_proplist);
		}
		delete m_data;
	}
}

int MPCORecorderLadruno::clearSources()
{
	for (size_t i = 0; i < m_data->node_channels.size(); ++i) {
		delete m_data->node_channels[i].source;
		delete m_data->node_channels[i].sink;
	}
	m_data->node_channels.clear();
	for (size_t i = 0; i < m_data->elem_channels.size(); ++i) {
		delete m_data->elem_channels[i].source;
		delete m_data->elem_channels[i].sink;
	}
	m_data->elem_channels.clear();
	for (size_t i = 0; i < m_data->domain_channels.size(); ++i) {
		delete m_data->domain_channels[i].source;
		delete m_data->domain_channels[i].sink;
	}
	m_data->domain_channels.clear();
	// element responses are owned here (frozen clearElementRecorders)
	for (size_t i = 0; i < m_data->elemental_responses.size(); ++i) {
		if (m_data->elemental_responses[i])
			delete m_data->elemental_responses[i];
	}
	m_data->elemental_responses.clear();
	m_data->elemental_recorders.clear();
	return 0;
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
	mpcolns::mpco::ProcessInfo& info = m_data->info;
	info.current_time_step_id = commitTag;
	info.current_time_step = timeStamp;

	// -T output-frequency gate (frozen record() do_record logic).
	bool do_record = false;
	if (!m_data->initialized) {
		do_record = true; // always record the first time
	}
	else {
		if (m_data->output_freq.type == mpcolns::mpco::OutputFrequency::NumberOfSteps) {
			if ((info.current_time_step_id - m_data->output_freq.last_step) >= m_data->output_freq.nsteps)
				do_record = true;
		}
		else if (m_data->output_freq.type == mpcolns::mpco::OutputFrequency::DeltaTime) {
			if (std::abs(info.current_time_step - m_data->output_freq.last_time) >= m_data->output_freq.dt)
				do_record = true;
		}
	}
	if (do_record) {
		m_data->output_freq.last_step = info.current_time_step_id;
		m_data->output_freq.last_time = info.current_time_step;
	}
	else {
		return 0;
	}

	if (!m_data->initialized) {
		if (initialize() != 0)
			return -1;
		m_data->initialized = true;
	}

	// Multi-stage detection (frozen MPCORecorder::record() rebuild_model block).
	// On the first record, and again whenever the domain-change stamp moves
	// (a new MODEL_STAGE — wipeAnalysis + new elements/nodes, staged analysis),
	// (re)write the model and rebuild every node/element source+sink. This is
	// what re-acquires fresh Element*/Response* pointers; without it the cached
	// pointers from the prior stage dangle and the new stage is never written.
	// NOTE: the parallel per-process stamp Allreduce (frozen lambdaHasDomainChanged)
	// is deferred with the rest of the parallel channel; this recorder is
	// single-process / per-partition for now.
	bool rebuild_model = false;
	int new_stamp = info.domain->hasDomainChanged();
	if (!m_data->first_domain_changed_done) {
		info.current_model_stage_id = new_stamp;
		m_data->first_domain_changed_done = true;
		rebuild_model = true;
	}
	else if (new_stamp != info.current_model_stage_id) {
		info.current_model_stage_id = new_stamp;
		rebuild_model = true;
	}
	if (rebuild_model) {
		if (writeModel() != 0)
			return -1;
	}

	if (recordResultsOnNodes() != 0)
		return -1;
	if (recordResultsOnElements() != 0)
		return -1;
	if (recordResultsOnDomain() != 0)
		return -1;

	if (info.h_file_id != mpcolns::HID_INVALID)
		mpcolns::h5::file::flush(info.h_file_id);

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
	hid_t h_prov = mpcolns::h5::group::create(
		h_info, "PROVENANCE", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	mpcolns::h5::group::close(h_prov);
	mpcolns::h5::group::close(h_info);

	// NOTE: model writing + source/sink building are NOT done here. They live in
	// writeModel(), driven by record()'s multi-stage rebuild_model block, so they
	// re-run on every MODEL_STAGE change (mirrors the frozen recorder, where
	// initialize() only sets up the file and record() owns the model rebuild).

	return 0;
}

/* ===================================================================== */
/* model writing                                                         */
/* ===================================================================== */

int MPCORecorderLadruno::writeModel()
{
	mpcolns::mpco::ProcessInfo& info = m_data->info;

	// info.current_model_stage_id is set by record() before calling writeModel()
	// (the multi-stage rebuild_model block), so the stamp is consistent across
	// the MODEL_STAGE group and every result/source path resolved from it.

	// MODEL_STAGE[<stamp>] + MODEL + RESULTS skeleton (mirror frozen writeModel)
	std::stringstream ss_stage;
	ss_stage << "MODEL_STAGE[" << info.current_model_stage_id << "]";
	hid_t h_stage = mpcolns::h5::group::create(
		info.h_file_id, ss_stage.str().c_str(), H5P_DEFAULT,
		info.h_group_proplist, H5P_DEFAULT);
	mpcolns::h5::attribute::write(h_stage, "STEP", info.current_time_step_id);
	mpcolns::h5::attribute::write(h_stage, "TIME", info.current_time_step);
	mpcolns::h5::attribute::write(h_stage, "KIND", std::string("static"));

	hid_t h_model = mpcolns::h5::group::create(
		h_stage, "MODEL", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	mpcolns::h5::group::close(h_model);

	hid_t h_results = mpcolns::h5::group::create(
		h_stage, "RESULTS", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	hid_t h_on_nodes = mpcolns::h5::group::create(
		h_results, "ON_NODES", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	mpcolns::h5::group::close(h_on_nodes);
	hid_t h_on_elems = mpcolns::h5::group::create(
		h_results, "ON_ELEMENTS", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	mpcolns::h5::group::close(h_on_elems);
	// ON_DOMAIN / ON_REGIONS skeleton (ADR D8 energy balance lands here)
	hid_t h_on_domain = mpcolns::h5::group::create(
		h_results, "ON_DOMAIN", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	mpcolns::h5::group::close(h_on_domain);
	hid_t h_on_regions = mpcolns::h5::group::create(
		h_results, "ON_REGIONS", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	mpcolns::h5::group::close(h_on_regions);
	mpcolns::h5::group::close(h_results);

	mpcolns::h5::group::close(h_stage);

	if (writeModelNodes() != 0)
		return -1;
	if (writeModelElements() != 0)
		return -1;
	if (writeModelSets() != 0)
		return -1;
	if (writeSections() != 0)
		return -1;

	// (Re)build the source/sink channels for this stage. clearSources() releases
	// the prior stage's sources, sinks, and cached element Response* (which would
	// otherwise dangle after the domain rebuild); the fresh sinks are not yet
	// initialized, so they re-create their result groups under the new MODEL_STAGE.
	// Mirrors frozen writeModel()'s trailing initNodeRecorders()/initElementRecorders().
	if (clearSources() != 0)
		return -1;
	if (initNodeSources() != 0)
		return -1;
	if (initElementSources() != 0)
		return -1;
	if (initDomainSources() != 0)
		return -1;

	return 0;
}

int MPCORecorderLadruno::writeModelNodes()
{
	mpcolns::mpco::ProcessInfo& info = m_data->info;
	const int ndim = info.num_dimensions;

	// gather the node set (skip pressure nodes), frozen writeModelNodes()
	m_data->nodes.clear();
	ID pressure_node_tags(0, info.domain->getNumPCs());
	Pressure_ConstraintIter& pc_iter = info.domain->getPCs();
	Pressure_Constraint* pc = 0;
	while ((pc = pc_iter()) != 0) {
		Node* pn = pc->getPressureNode();
		if (pn)
			pressure_node_tags.insert(pn->getTag());
	}
	if (m_data->has_region) {
		for (size_t i = 0; i < m_data->node_set.size(); ++i) {
			int nd = m_data->node_set[i];
			if (pressure_node_tags.getLocationOrdered(nd) < 0) {
				Node* cn = info.domain->getNode(nd);
				if (cn)
					m_data->nodes.push_back(cn);
			}
		}
	}
	else {
		NodeIter& nit = info.domain->getNodes();
		Node* nptr = 0;
		while ((nptr = nit()) != 0) {
			int nd = nptr->getTag();
			if (pressure_node_tags.getLocationOrdered(nd) < 0)
				m_data->nodes.push_back(nptr);
		}
	}

	const size_t nnodes = m_data->nodes.size();
	if (nnodes == 0) {
		opserr << "MPCORecorderLadruno Error: no nodes to write\n";
		return -1;
	}

	std::vector<int> node_ids(nnodes);
	std::vector<double> node_coords(nnodes * (size_t)ndim);
	for (size_t i = 0; i < nnodes; ++i) {
		Node* inode = m_data->nodes[i];
		node_ids[i] = inode->getTag();
		const Vector& crds = inode->getCrds();
		size_t j = i * (size_t)ndim;
		node_coords[j] = crds[0];
		if (ndim > 1) {
			node_coords[j + 1] = crds[1];
			if (ndim > 2)
				node_coords[j + 2] = crds[2];
		}
	}

	std::stringstream ss_nodes;
	ss_nodes << "MODEL_STAGE[" << info.current_model_stage_id << "]/MODEL/NODES";
	hid_t h_nodes = mpcolns::h5::group::create(
		info.h_file_id, ss_nodes.str().c_str(), H5P_DEFAULT,
		info.h_group_proplist, H5P_DEFAULT);
	hid_t d_id = mpcolns::h5::dataset::createAndWrite(h_nodes, "ID", node_ids);
	mpcolns::h5::dataset::close(d_id);
	hid_t d_coord = mpcolns::h5::dataset::createAndWrite(
		h_nodes, "COORDINATES", node_coords, nnodes, (size_t)ndim);
	mpcolns::h5::dataset::close(d_coord);
	mpcolns::h5::group::close(h_nodes);

	return 0;
}

int MPCORecorderLadruno::writeModelElements()
{
	mpcolns::mpco::ProcessInfo& info = m_data->info;

	std::stringstream ss_elems;
	ss_elems << "MODEL_STAGE[" << info.current_model_stage_id << "]/MODEL/ELEMENTS";
	hid_t h_gp_elements = mpcolns::h5::group::create(
		info.h_file_id, ss_elems.str().c_str(), H5P_DEFAULT,
		info.h_group_proplist, H5P_DEFAULT);

	// map all elements (also populates registered custom rules)
	m_data->elements.mapElements(info.domain, m_data->has_region, m_data->elem_set);

	namespace me = mpcolns::mpco::element;
	for (me::ElementCollection::submap_type::iterator it1 = m_data->elements.items.begin();
		it1 != m_data->elements.items.end(); ++it1) {
		me::ElementWithSameClassTagCollection& elem_by_tag = it1->second;
		for (me::ElementWithSameClassTagCollection::submap_type::iterator it2 = elem_by_tag.items.begin();
			it2 != elem_by_tag.items.end(); ++it2) {
			me::ElementWithSameIntRuleCollection& elem_by_rule = it2->second;
			for (me::ElementWithSameIntRuleCollection::submap_type::iterator it3 = elem_by_rule.items.begin();
				it3 != elem_by_rule.items.end(); ++it3) {
				me::ElementWithSameCustomIntRuleCollection& elem_by_custom_rule = it3->second;
				if (elem_by_custom_rule.items.empty())
					continue;
				// name <classTag>-<className>[<int_rule>:<custom_rule_index>]
				std::stringstream ss_dset;
				ss_dset << elem_by_tag.class_tag << "-" << elem_by_tag.class_name
					<< "[" << elem_by_rule.int_rule_type << ":"
					<< elem_by_custom_rule.custom_int_rule_index << "]";
				elem_by_custom_rule.name = ss_dset.str();

				// CONNECTIVITY [nElem x (1 + num_nodes)]: (elemTag, c1..cK)
				std::vector<int> buffer(elem_by_custom_rule.items.size() * (size_t)(1 + elem_by_tag.num_nodes));
				size_t offset = 0;
				for (me::ElementWithSameCustomIntRuleCollection::collection_type::iterator it4 = elem_by_custom_rule.items.begin();
					it4 != elem_by_custom_rule.items.end(); ++it4) {
					Element* elem = *it4;
					buffer[offset] = elem->getTag();
					for (int j = 0; j < elem_by_tag.num_nodes; ++j)
						buffer[(size_t)j + 1 + offset] = elem->getExternalNodes()(j);
					offset += (size_t)(1 + elem_by_tag.num_nodes);
				}
				// Element group <name> is a GROUP (schema v1) holding CONNECTIVITY +
				// BASIS attrs + (for custom rules) a QUADRATURE child. A dataset cannot
				// parent the QUADRATURE group, which silently failed before this.
				hid_t h_eg = mpcolns::h5::group::create(
					h_gp_elements, elem_by_custom_rule.name.c_str(), H5P_DEFAULT,
					info.h_group_proplist, H5P_DEFAULT);
				hid_t d_conn = mpcolns::h5::dataset::createAndWrite(
					h_eg, "CONNECTIVITY", buffer,
					elem_by_custom_rule.items.size(), (size_t)(1 + elem_by_tag.num_nodes));
				mpcolns::h5::dataset::close(d_conn);

				// BASIS / QUADRATURE descriptors (schema §3) DERIVED from the
				// legacy (geometry, integration rule) pair. We keep the legacy
				// GEOMETRY/INTEGRATION_RULE attrs (for tooling that still reads them)
				// AND add the derived self-describing descriptors.
				mpcolns::mpco::ElementGeometryType::Enum geom = elem_by_tag.geom_type;
				mpcolns::mpco::ElementIntegrationRuleType::Enum irule = elem_by_rule.int_rule_type;
				mpcolns::h5::attribute::write(h_eg, "GEOMETRY", (int)geom);
				mpcolns::h5::attribute::write(h_eg, "INTEGRATION_RULE", (int)irule);

				// derive TOPOLOGY / FAMILY / ORDER / PARAM_DOMAIN / NUM_CTRL
				std::string topology = "custom";
				std::string param_domain = "[-1,1]";
				std::vector<int> order;
				switch (geom) {
				case mpcolns::mpco::ElementGeometryType::Line_2N:
				case mpcolns::mpco::ElementGeometryType::Line_3N:
					topology = "line"; param_domain = "[-1,1]"; order.push_back(1); break;
				case mpcolns::mpco::ElementGeometryType::Triangle_3N:
				case mpcolns::mpco::ElementGeometryType::Triangle_6N:
					topology = "tri"; param_domain = "bary"; order.push_back(1); break;
				case mpcolns::mpco::ElementGeometryType::Quadrilateral_4N:
				case mpcolns::mpco::ElementGeometryType::Quadrilateral_8N:
				case mpcolns::mpco::ElementGeometryType::Quadrilateral_9N:
				case mpcolns::mpco::ElementGeometryType::Quadrilateral_CohesiveBand_4N:
					topology = "quad"; param_domain = "[-1,1]"; order.push_back(1); order.push_back(1); break;
				case mpcolns::mpco::ElementGeometryType::Tetrahedron_4N:
				case mpcolns::mpco::ElementGeometryType::Tetrahedron_10N:
					topology = "tet"; param_domain = "bary"; order.push_back(1); break;
				case mpcolns::mpco::ElementGeometryType::Hexahedron_8N:
				case mpcolns::mpco::ElementGeometryType::Hexahedron_20N:
				case mpcolns::mpco::ElementGeometryType::Hexahedron_27N:
					topology = "hex"; param_domain = "[-1,1]"; order.push_back(1); order.push_back(1); order.push_back(1); break;
				default:
					topology = "custom"; break;
				}
				// A self-describing element (Element-contract Part A) overrides the
				// geometry-derived guesses above. Lagrange elements don't answer the
				// "basisInfo" probe (basis_info.valid == false) so the derived defaults
				// are kept; a Bernstein element (BezierTri6) instead records
				// FAMILY=bernstein, ORDER=2.
				const me::BasisInfo& binfo = elem_by_tag.basis_info;
				std::string family = "lagrange";
				int rational = 0;
				int num_ctrl = elem_by_tag.num_nodes;
				if (binfo.valid) {
					if (!binfo.family.empty())       family = binfo.family;
					if (!binfo.topology.empty())     topology = binfo.topology;
					if (!binfo.param_domain.empty()) param_domain = binfo.param_domain;
					if (binfo.rational >= 0)         rational = binfo.rational;
					if (binfo.num_ctrl > 0)          num_ctrl = binfo.num_ctrl;
					if (binfo.order > 0) {
						order.clear();
						// ORDER carries one entry per GP_PARAM column (the validator equates
						// len(ORDER) with the parametric dimension). For a multi-dim custom
						// rule (BezierTri6 barycentric: 2 free area coords) replicate the
						// declared total degree across the parametric directions.
						int ndir = 1;
						if (irule == mpcolns::mpco::ElementIntegrationRuleType::CustomIntegrationRule
							&& elem_by_custom_rule.custom_int_rule_index != 0) {
							int d = m_data->elements.registered_custom_rules[
								elem_by_custom_rule.custom_int_rule_index].custom_rule_dimension;
							if (d > 1) ndir = d;
						}
						for (int k = 0; k < ndir; ++k) order.push_back(binfo.order);
					}
				}
				mpcolns::h5::attribute::write(h_eg, "TOPOLOGY", topology);
				mpcolns::h5::attribute::write(h_eg, "FAMILY", family);
				if (order.size() > 0)
					mpcolns::h5::attribute::write(h_eg, "ORDER", order);
				mpcolns::h5::attribute::write(h_eg, "PARAM_DOMAIN", param_domain);
				mpcolns::h5::attribute::write(h_eg, "RATIONAL", (int)rational);
				mpcolns::h5::attribute::write(h_eg, "NUM_CTRL", (int)num_ctrl);

				// ---- Resolve the quadrature rule: custom -> standard table -> none ----
				// QUADRATURE/{GP_PARAM[NUM_GP x NDIR], GP_WEIGHT[NUM_GP]} + NDIR + NUM_GP
				// are now written for BOTH custom-rule elements (BezierTri6 etc.) AND
				// standard Gauss-Legendre elements (quad/hex/tri/tet/line). For standard
				// rules the GP coords + weights come from the built-in getStandardQuadrature
				// table ("parity by derivation" — STKO held this reader-side; we move it
				// write-side so the reader never KeyErrors on a legacy element). The bucket
				// <name> is a GROUP (see above), so the QUADRATURE child group is valid.
				std::vector<double> q_gp_param, q_gp_weight;
				int q_num_gp = 0, q_ndir = 0;
				bool have_rule = false;
				if (irule == mpcolns::mpco::ElementIntegrationRuleType::CustomIntegrationRule) {
					mpcolns::h5::attribute::write(h_eg, "CUSTOM_INTEGRATION_RULE", elem_by_custom_rule.custom_int_rule_index);
					if (elem_by_custom_rule.custom_int_rule_index != 0) {
						me::ElementIntegrationRule& custom_rule =
							m_data->elements.registered_custom_rules[elem_by_custom_rule.custom_int_rule_index];
						// legacy GP_X attribute (flat coords) kept for tooling that reads it.
						mpcolns::h5::attribute::write(h_eg, "GP_X", custom_rule.x);
						mpcolns::h5::attribute::write(h_eg, "CUSTOM_INTEGRATION_RULE_DIMENSION", custom_rule.custom_rule_dimension);
						// reshape: legacy line rules store 1 coord/point (num_gp==0 => cols=1);
						// multi-dim parametric rules (e.g. BezierTri6 barycentric) store num_gp
						// rows of (x.size()/num_gp) coords each.
						q_num_gp = custom_rule.num_gp > 0 ? custom_rule.num_gp : (int)custom_rule.x.size();
						q_ndir = (custom_rule.num_gp > 0 && custom_rule.num_gp <= (int)custom_rule.x.size())
							? (int)custom_rule.x.size() / custom_rule.num_gp : 1;
						q_gp_param = custom_rule.x;
						q_gp_weight = custom_rule.w;
						have_rule = (q_num_gp > 0 && q_ndir > 0);
					}
				}
				else {
					// standard Gauss rule (returns false for untabulated rules ->
					// no QUADRATURE written; the reader tolerates its absence).
					have_rule = me::getStandardQuadrature(irule, q_gp_param, q_gp_weight, q_num_gp, q_ndir);
				}

				if (have_rule) {
					mpcolns::h5::attribute::write(h_eg, "NDIR", q_ndir);
					mpcolns::h5::attribute::write(h_eg, "NUM_GP", q_num_gp);
					// QUADRATURE child group: GP_PARAM (NUM_GP x NDIR) + optional GP_WEIGHT.
					hid_t h_q = mpcolns::h5::group::create(
						h_eg, "QUADRATURE", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
					hid_t d_gpp = mpcolns::h5::dataset::createAndWrite(
						h_q, "GP_PARAM", q_gp_param, (size_t)q_num_gp, (size_t)q_ndir);
					mpcolns::h5::dataset::close(d_gpp);
					if (!q_gp_weight.empty()) {
						hid_t d_gpw = mpcolns::h5::dataset::createAndWrite(
							h_q, "GP_WEIGHT", q_gp_weight);  // 1-D [nGP] per schema v1
						mpcolns::h5::dataset::close(d_gpw);
					}
					mpcolns::h5::group::close(h_q);

					// ---- belt-and-suspenders GLOBAL_GP_COORDS (ADR D2) -------------
					// Static, reference-config, C++-computed x(ξ_g)=ΣN_i(ξ_g)X_i per
					// element via the write-side basis evaluator. Stored 2-D
					// [nElem x (NUM_GP*ndim)] because MPCOL_Hdf5 has no rank-3 writer;
					// reader reshapes to [nElem, NUM_GP, ndim]. Skipped (graceful) when
					// computeGlobalGP doesn't yet have the topology's basis. Reuses the
					// CONNECTIVITY node access pattern (domain->getNode->getCrds).
					int ndim = 0;
					{
						const ID& nd0 = elem_by_custom_rule.items.front()->getExternalNodes();
						Node* n0 = (nd0.Size() > 0) ? info.domain->getNode(nd0(0)) : 0;
						if (n0) ndim = n0->getCrds().Size();
					}
					if (ndim > 0) {
						size_t nElem = elem_by_custom_rule.items.size();
						std::vector<double> ggp;
						ggp.reserve(nElem * (size_t)q_num_gp * (size_t)ndim);
						bool ok = true;
						for (me::ElementWithSameCustomIntRuleCollection::collection_type::iterator it5 = elem_by_custom_rule.items.begin();
							it5 != elem_by_custom_rule.items.end() && ok; ++it5) {
							Element* elem = *it5;
							const ID& nd = elem->getExternalNodes();
							int nn = elem_by_tag.num_nodes;
							std::vector<std::vector<double> > node_xyz(
								(size_t)nn, std::vector<double>((size_t)ndim, 0.0));
							for (int j = 0; j < nn && ok; ++j) {
								Node* nj = (j < nd.Size()) ? info.domain->getNode(nd(j)) : 0;
								if (nj == 0) { ok = false; break; }
								const Vector& c = nj->getCrds();
								for (int d = 0; d < ndim && d < c.Size(); ++d)
									node_xyz[(size_t)j][(size_t)d] = c(d);
							}
							if (!ok) break;
							std::vector<double> ex;
							if (!me::computeGlobalGP(geom, nn, node_xyz, ndim,
								q_gp_param, q_num_gp, q_ndir, ex)) { ok = false; break; }
							ggp.insert(ggp.end(), ex.begin(), ex.end());
						}
						if (ok && !ggp.empty()) {
							hid_t d_g = mpcolns::h5::dataset::createAndWrite(
								h_eg, "GLOBAL_GP_COORDS", ggp,
								nElem, (size_t)q_num_gp * (size_t)ndim);
							mpcolns::h5::dataset::close(d_g);
						}
					}
				}
				mpcolns::h5::group::close(h_eg);
			}
		}
	}

	mpcolns::h5::group::close(h_gp_elements);
	return 0;
}

int MPCORecorderLadruno::writeModelSets()
{
	// MODEL/SETS — one SET_<tag> per region whose identity the file must be
	// self-describing about (ADR D8 / H): the per-region energy tags plus the
	// -R restriction region (if any). Each SET_<tag> carries NODES + ELEMENTS
	// tag datasets and a TAG attr. If no region is referenced, the (empty) SETS
	// group is still created for a consistent layout.
	mpcolns::mpco::ProcessInfo& info = m_data->info;

	// union of region tags to describe (preserve order, drop duplicates)
	std::vector<int> set_tags;
	{
		std::set<int> seen;
		if (m_data->has_region && m_data->region_tag != 0) {
			if (seen.insert(m_data->region_tag).second)
				set_tags.push_back(m_data->region_tag);
		}
		for (size_t i = 0; i < m_data->energy_region_tags.size(); ++i) {
			int t = m_data->energy_region_tags[i];
			if (seen.insert(t).second)
				set_tags.push_back(t);
		}
	}

	// Nothing to describe -> don't create an empty MODEL/SETS group (it adds
	// benign HDF5 open-probe noise on region-free models for no benefit).
	if (set_tags.empty())
		return 0;

	std::stringstream ss_sets;
	ss_sets << "MODEL_STAGE[" << info.current_model_stage_id << "]/MODEL/SETS";
	hid_t h_sets = mpcolns::h5::group::create(
		info.h_file_id, ss_sets.str().c_str(), H5P_DEFAULT,
		info.h_group_proplist, H5P_DEFAULT);

	for (size_t i = 0; i < set_tags.size(); ++i) {
		int tag = set_tags[i];
		MeshRegion* reg = info.domain->getRegion(tag);
		if (reg == 0) {
			opserr << "MPCORecorderLadruno - WARNING region " << tag
			       << " referenced by the recorder was not found; its MODEL/SETS "
			          "entry is skipped\n";
			continue;
		}
		std::stringstream ss;
		ss << "MODEL_STAGE[" << info.current_model_stage_id
		   << "]/MODEL/SETS/SET_" << tag;
		hid_t h = mpcolns::h5::group::create(
			info.h_file_id, ss.str().c_str(), H5P_DEFAULT,
			info.h_group_proplist, H5P_DEFAULT);
		mpcolns::h5::attribute::write(h, "TAG", tag);

		const ID& rnodes = reg->getNodes();
		std::vector<int> node_ids((size_t)rnodes.Size());
		for (int k = 0; k < rnodes.Size(); ++k)
			node_ids[(size_t)k] = rnodes(k);
		const ID& relems = reg->getElements();
		std::vector<int> elem_ids((size_t)relems.Size());
		for (int k = 0; k < relems.Size(); ++k)
			elem_ids[(size_t)k] = relems(k);

		hid_t d_n = mpcolns::h5::dataset::createAndWrite(h, "NODES", node_ids);
		if (d_n != mpcolns::HID_INVALID)
			mpcolns::h5::dataset::close(d_n);
		hid_t d_e = mpcolns::h5::dataset::createAndWrite(h, "ELEMENTS", elem_ids);
		if (d_e != mpcolns::HID_INVALID)
			mpcolns::h5::dataset::close(d_e);

		mpcolns::h5::group::close(h);
	}

	mpcolns::h5::group::close(h_sets);
	return 0;
}

int MPCORecorderLadruno::writeSections()
{
	mpcolns::mpco::ProcessInfo& info = m_data->info;
	namespace me = mpcolns::mpco::element;

	std::map<int, me::SectionAssignment> sec_assignments;
	std::map<me::FiberSectionData, me::SectionAssignment> aux_sec_assignments;

	m_data->elem_ngauss_nfiber_info.clear();

	std::stringstream ss_sec_dir;
	ss_sec_dir << "MODEL_STAGE[" << info.current_model_stage_id << "]/MODEL/SECTION_ASSIGNMENTS";
	hid_t h_gp_sec = mpcolns::h5::group::create(
		info.h_file_id, ss_sec_dir.str().c_str(), H5P_DEFAULT,
		info.h_group_proplist, H5P_DEFAULT);

	// loop all element groups (frozen writeSections section/fiber discovery)
	for (me::ElementCollection::submap_type::iterator it1 = m_data->elements.items.begin();
		it1 != m_data->elements.items.end(); ++it1) {
		me::ElementWithSameClassTagCollection& elem_by_tag = it1->second;
		for (me::ElementWithSameClassTagCollection::submap_type::iterator it2 = elem_by_tag.items.begin();
			it2 != elem_by_tag.items.end(); ++it2) {
			me::ElementWithSameIntRuleCollection& elem_by_rule = it2->second;
			for (me::ElementWithSameIntRuleCollection::submap_type::iterator it3 = elem_by_rule.items.begin();
				it3 != elem_by_rule.items.end(); ++it3) {
				me::ElementWithSameCustomIntRuleCollection& elem_by_custom_rule = it3->second;
				for (std::vector<Element*>::iterator it4 = elem_by_custom_rule.items.begin();
					it4 != elem_by_custom_rule.items.end(); ++it4) {
					Element* elem = *it4;

					std::string request1 = mpcolns::utils::shell::isShellElementTag(elem->getClassTag()) ? "material" : "section";
					std::string request2 = "fiber";
					std::string request3 = "stress";
					int argc = 5;
					const char** argv = new const char*[argc];
					argv[0] = request1.c_str();
					argv[2] = request2.c_str();
					argv[4] = request3.c_str();
					int woagg_agrc = 6;
					const char** woagg_argv = new const char*[woagg_agrc];
					woagg_argv[0] = request1.c_str();
					woagg_argv[2] = request1.c_str();
					woagg_argv[3] = request2.c_str();
					woagg_argv[5] = request3.c_str();

					std::vector<int> elem_gauss_id;
					std::vector<int> elem_sec_id;
					std::vector<bool> elem_dummy_sec_flags;
					std::vector<me::FiberSectionData> elem_sections;
					std::vector<int> elem_fiber_base_index;

					int trial_num_sec = 0;
					while (true) {
						trial_num_sec++;
						if (trial_num_sec > MPCO_MAX_TRIAL_NSEC)
							break;
						std::stringstream ss_ts; ss_ts << trial_num_sec;
						std::string s_ts = ss_ts.str();
						argv[1] = s_ts.c_str();
						woagg_argv[1] = s_ts.c_str();

						int trial_num_fib = -1;
						bool break_sec_loop = false;
						bool first_fiber_done = false;
						bool do_workaround_for_aggregator = false;
						while (true) {
							trial_num_fib++;
							if (trial_num_fib > MPCO_LADRUNO_MAX_TRIAL_NFIB)
								break;
							std::stringstream ss_tf; ss_tf << trial_num_fib;
							std::string s_tf = ss_tf.str();
							argv[3] = s_tf.c_str();
							woagg_argv[4] = s_tf.c_str();

							me::OutputDescriptor eo_descriptor;
							me::OutputDescriptorStream eo_stream(&eo_descriptor);
							Response* eo_response = do_workaround_for_aggregator ?
								elem->setResponse(woagg_argv, woagg_agrc, eo_stream) :
								elem->setResponse(argv, argc, eo_stream);
							eo_stream.finalizeSetResponse();
							if (eo_response)
								delete eo_response;

							std::vector<me::FiberData> trial_fiberdata;
							std::vector<int> trial_fiber_material_id;
							std::vector<int> trial_sec_id;
							std::vector<int> trial_gp_id;
							std::vector<bool> trial_dummy_flag;
							eo_descriptor.getFiberData(trial_fiberdata, trial_fiber_material_id,
								trial_sec_id, trial_gp_id, trial_dummy_flag);

							if (trial_gp_id.size() == 0) {
								break_sec_loop = true;
								break;
							}
							if (trial_sec_id.size() == 0)
								break;

							if (trial_sec_id.size() != trial_fiberdata.size()) trial_sec_id.resize(trial_fiberdata.size());
							if (trial_gp_id.size() != trial_fiberdata.size()) trial_gp_id.resize(trial_fiberdata.size());
							if (trial_dummy_flag.size() != trial_fiberdata.size()) trial_dummy_flag.resize(trial_fiberdata.size());

							if (trial_fiberdata.size() == 0 && trial_num_fib == 0)
								continue;

							if (!first_fiber_done) {
								elem_gauss_id.push_back(trial_gp_id[0]);
								elem_sec_id.push_back(trial_sec_id[0]);
								elem_dummy_sec_flags.push_back(trial_dummy_flag[0]);
								elem_sections.push_back(me::FiberSectionData());
								elem_fiber_base_index.push_back(trial_num_fib);
								first_fiber_done = true;
							}

							if (trial_fiberdata.size() == 0) {
								if (!do_workaround_for_aggregator && (trial_num_fib == 1)) {
									do_workaround_for_aggregator = true;
									if (first_fiber_done) {
										elem_gauss_id.pop_back();
										elem_sec_id.pop_back();
										elem_dummy_sec_flags.pop_back();
										elem_sections.pop_back();
										elem_fiber_base_index.pop_back();
									}
									trial_num_fib = -1;
									break_sec_loop = false;
									first_fiber_done = false;
									continue;
								}
								break;
							}

							if (trial_fiberdata.size() == 1) {
								me::FiberSectionData& curr = elem_sections.back();
								curr.fibers.push_back(trial_fiberdata[0]);
								curr.materials.push_back(trial_fiber_material_id[0]);
							}
						} // fiber loop
						if (break_sec_loop)
							break;
					} // section loop

					delete[] argv;
					delete[] woagg_argv;

					// store ngauss/nfiber info for the element-result section path
					std::vector<std::pair<int, int> >& nfpg = m_data->elem_ngauss_nfiber_info[elem->getTag()];
					nfpg.resize(elem_gauss_id.size());
					for (size_t igp = 0; igp < elem_gauss_id.size(); ++igp) {
						nfpg[igp].first = elem_fiber_base_index[igp];
						nfpg[igp].second = (int)elem_sections[igp].fibers.size();
					}

					size_t ngauss = elem_gauss_id.size();
					for (size_t igp = 0; igp < ngauss; ++igp) {
						int igp_sec_id = elem_sec_id[igp];
						bool igp_dummy = elem_dummy_sec_flags[igp];
						me::FiberSectionData& igp_data = elem_sections[igp];
						if (igp_dummy) {
							me::SectionAssignment& ca = aux_sec_assignments[igp_data];
							if (ca.is_new) { ca.fiber_section_data = igp_data; ca.is_new = false; }
							ca.assignments.push_back(me::ElemGaussPair(elem->getTag(), (int)igp));
						}
						else {
							me::SectionAssignment& ca = sec_assignments[igp_sec_id];
							if (ca.is_new) { ca.fiber_section_data = igp_data; ca.is_new = false; }
							ca.assignments.push_back(me::ElemGaussPair(elem->getTag(), (int)igp));
						}
					}
				} // each element
			}
		}
	}

	// resolve section class names
	for (std::map<int, me::SectionAssignment>::iterator it = sec_assignments.begin();
		it != sec_assignments.end(); ++it) {
		SectionForceDeformation* sfd = OPS_getSectionForceDeformation(it->first);
		if (sfd)
			it->second.name = sfd->getClassType();
	}

	// move dummy (floating-fiber) sections to tags above max
	int max_sec_tag = sec_assignments.size() > 0 ? (--sec_assignments.end())->first : 0;
	for (std::map<me::FiberSectionData, me::SectionAssignment>::iterator it = aux_sec_assignments.begin();
		it != aux_sec_assignments.end(); ++it) {
		sec_assignments[++max_sec_tag] = it->second;
	}

	// write each assignment (schema §6: KIND=fiber)
	for (std::map<int, me::SectionAssignment>::iterator it = sec_assignments.begin();
		it != sec_assignments.end(); ++it) {
		me::SectionAssignment& sec_asn = it->second;
		if (sec_asn.assignments.size() == 0)
			continue;
		std::stringstream ss_name;
		ss_name << "SECTION_" << it->first << "[" << sec_asn.name << "]";
		hid_t h_isec = mpcolns::h5::group::create(
			h_gp_sec, ss_name.str().c_str(), H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
		mpcolns::h5::attribute::write(h_isec, "ID", it->first);
		mpcolns::h5::attribute::write(h_isec, "NAME", sec_asn.name);
		mpcolns::h5::attribute::write(h_isec, "KIND", std::string("fiber"));
		{
			hid_t dset = mpcolns::h5::dataset::createAndWrite(h_isec, "ASSIGNMENT", sec_asn.assignments);
			mpcolns::h5::dataset::close(dset);
		}
		{
			hid_t dset = mpcolns::h5::dataset::createAndWrite(h_isec, "FIBER_DATA", sec_asn.fiber_section_data.fibers);
			mpcolns::h5::dataset::close(dset);
		}
		{
			hid_t dset = mpcolns::h5::dataset::createAndWrite(h_isec, "FIBER_MATERIALS", sec_asn.fiber_section_data.materials);
			mpcolns::h5::dataset::close(dset);
		}
		mpcolns::h5::group::close(h_isec);
	}

	mpcolns::h5::group::close(h_gp_sec);
	return 0;
}

/* ===================================================================== */
/* node sources                                                          */
/* ===================================================================== */

int MPCORecorderLadruno::initNodeSources()
{
	mpcolns::mpco::ProcessInfo& info = m_data->info;
	namespace NR = mpcolns::mpco; // for NodalResultType

	// node-tag set passed to each source (matches writeModelNodes order)
	std::vector<int> node_ids(m_data->nodes.size());
	for (size_t i = 0; i < m_data->nodes.size(); ++i)
		node_ids[i] = m_data->nodes[i]->getTag();

	for (size_t i = 0; i < m_data->nodal_results_requests.size(); ++i) {
		NR::NodalResultType::Enum rtype = m_data->nodal_results_requests[i];
		int grad = m_data->sens_grad_indices[i];

		private_data::NodeChannel ch;
		mpcolns::ResultSource* src = 0;
		int reac = -1;
		bool is_modes = false;

		switch (rtype) {
		case NR::NodalResultType::Displacement:
			src = new mpcolns::DisplacementSource(info, node_ids); break;
		case NR::NodalResultType::Rotation:
			src = new mpcolns::RotationSource(info, node_ids); break;
		case NR::NodalResultType::Velocity:
			src = new mpcolns::VelocitySource(info, node_ids); break;
		case NR::NodalResultType::AngularVelocity:
			src = new mpcolns::AngularVelocitySource(info, node_ids); break;
		case NR::NodalResultType::Acceleration:
			src = new mpcolns::AccelerationSource(info, node_ids); break;
		case NR::NodalResultType::AngularAcceleration:
			src = new mpcolns::AngularAccelerationSource(info, node_ids); break;
		case NR::NodalResultType::Pressure:
			src = new mpcolns::PressureSource(info, node_ids); break;
		case NR::NodalResultType::ReactionForce: {
			mpcolns::ReactionForceSource* s = new mpcolns::ReactionForceSource(info, node_ids);
			reac = s->reactionFlag(); src = s; break; }
		case NR::NodalResultType::ReactionMoment: {
			mpcolns::ReactionMomentSource* s = new mpcolns::ReactionMomentSource(info, node_ids);
			reac = s->reactionFlag(); src = s; break; }
		case NR::NodalResultType::ReactionForceIncludingInertia: {
			mpcolns::ReactionForceIncInertiaSource* s = new mpcolns::ReactionForceIncInertiaSource(info, node_ids);
			reac = s->reactionFlag(); src = s; break; }
		case NR::NodalResultType::ReactionMomentIncludingInertia: {
			mpcolns::ReactionMomentIncInertiaSource* s = new mpcolns::ReactionMomentIncInertiaSource(info, node_ids);
			reac = s->reactionFlag(); src = s; break; }
		case NR::NodalResultType::RayleighForce: {
			mpcolns::RayleighForceSource* s = new mpcolns::RayleighForceSource(info, node_ids);
			reac = s->reactionFlag(); src = s; break; }
		case NR::NodalResultType::RayleighMoment: {
			mpcolns::RayleighMomentSource* s = new mpcolns::RayleighMomentSource(info, node_ids);
			reac = s->reactionFlag(); src = s; break; }
		case NR::NodalResultType::UnbalancedForce:
			src = new mpcolns::UnbalancedForceSource(info, node_ids); break;
		case NR::NodalResultType::UnbalancedMoment:
			src = new mpcolns::UnbalancedMomentSource(info, node_ids); break;
		case NR::NodalResultType::UnbalancedForceIncludingInertia:
			src = new mpcolns::UnbalancedForceIncInertiaSource(info, node_ids); break;
		case NR::NodalResultType::UnbalancedMomentIncludingInertia:
			src = new mpcolns::UnbalancedMomentIncInertiaSource(info, node_ids); break;
		case NR::NodalResultType::ModesOfVibration:
			src = new mpcolns::ModesOfVibrationSource(info, node_ids); is_modes = true; break;
		case NR::NodalResultType::ModesOfVibrationRotational:
			src = new mpcolns::ModesOfVibrationRotationalSource(info, node_ids); is_modes = true; break;
		case NR::NodalResultType::DisplacementSensitivity:
			src = new mpcolns::DisplacementSensitivitySource(info, node_ids, grad); break;
		case NR::NodalResultType::RotationSensitivity:
			src = new mpcolns::RotationSensitivitySource(info, node_ids, grad); break;
		case NR::NodalResultType::VelocitySensitivity:
			src = new mpcolns::VelocitySensitivitySource(info, node_ids, grad); break;
		case NR::NodalResultType::AngularVelocitySensitivity:
			src = new mpcolns::AngularVelocitySensitivitySource(info, node_ids, grad); break;
		case NR::NodalResultType::AccelerationSensitivity:
			src = new mpcolns::AccelerationSensitivitySource(info, node_ids, grad); break;
		case NR::NodalResultType::AngularAccelerationSensitivity:
			src = new mpcolns::AngularAccelerationSensitivitySource(info, node_ids, grad); break;
		default:
			break;
		}

		if (src == 0)
			continue;

		ch.source = src;
		ch.sink = new mpcolns::StreamingSink(mpcolns::ResultFamily::OnNodes);
		ch.reaction_flag = reac;
		ch.is_modes = is_modes;
		m_data->node_channels.push_back(ch);
	}

	return 0;
}

/* ===================================================================== */
/* element sources                                                       */
/* ===================================================================== */

int MPCORecorderLadruno::initElementSources()
{
	mpcolns::mpco::ProcessInfo& info = m_data->info;
	namespace me = mpcolns::mpco::element;

	if (m_data->elemental_results_requests.empty())
		return 0;

	// Build the per-request response-collection trees (port of initElementRecorders).
	m_data->elemental_recorders.resize(m_data->elemental_results_requests.size());

	for (size_t cnt_request = 0; cnt_request < m_data->elemental_results_requests.size(); ++cnt_request) {
		const std::vector<std::string>& request = m_data->elemental_results_requests[cnt_request];

		std::vector<std::string> request_mod;
		bool do_all_materials = false;
		bool do_all_sections = false;
		bool do_all_fibers = false;
		size_t material_id_placeholder_index = 0;
		size_t section_id_placeholder_index = 0;
		size_t fiber_id_placeholder_index = 0;
		for (size_t i = 0; i < request.size(); ++i) {
			const std::string& cr = request[i];
			request_mod.push_back(cr);
			if (i == 0) {
				if (cr == "section") {
					request_mod.push_back("");
					section_id_placeholder_index = request_mod.size() - 1;
					do_all_sections = true;
				}
				else if (cr == "material") {
					request_mod.push_back("");
					material_id_placeholder_index = request_mod.size() - 1;
					do_all_materials = true;
				}
			}
			else {
				if (cr == "fiber" && do_all_sections) {
					request_mod.push_back("");
					fiber_id_placeholder_index = request_mod.size() - 1;
					do_all_fibers = true;
				}
			}
		}
		int argc = (int)request_mod.size();
		const char** argv = new const char*[argc];
		for (size_t i = 0; i < request_mod.size(); ++i)
			argv[i] = request_mod[i].c_str();

		std::string aux_section_keyword_standard = "section";
		std::string aux_section_keyword_for_shells = "material";

		me::ResultRecorder& recorder = m_data->elemental_recorders[cnt_request];
		recorder.result_request = request;

		for (me::ElementCollection::submap_type::iterator it1 = m_data->elements.items.begin();
			it1 != m_data->elements.items.end(); ++it1) {
			me::ElementWithSameClassTagCollection& elem_by_tag = it1->second;
			me::OutputWithSameClassTagCollection& eo_by_tag = recorder.response_map[it1->first];

			bool standard_section_keyword_modified = false;
			if (do_all_sections && mpcolns::utils::shell::isShellElementTag(it1->first)) {
				argv[0] = aux_section_keyword_for_shells.c_str();
				standard_section_keyword_modified = true;
			}

			for (me::ElementWithSameClassTagCollection::submap_type::iterator it2 = elem_by_tag.items.begin();
				it2 != elem_by_tag.items.end(); ++it2) {
				me::ElementWithSameIntRuleCollection& elem_by_rule = it2->second;
				me::OutputWithSameIntRuleCollection& eo_by_rule = eo_by_tag.items[it2->first];
				for (me::ElementWithSameIntRuleCollection::submap_type::iterator it3 = elem_by_rule.items.begin();
					it3 != elem_by_rule.items.end(); ++it3) {
					me::ElementWithSameCustomIntRuleCollection& elem_by_custom_rule = it3->second;
					me::OutputWithSameCustomIntRuleCollection& eo_by_custom_rule = eo_by_rule.items[it3->first];
					int header_local_index = 0;
					for (me::ElementWithSameCustomIntRuleCollection::collection_type::iterator it4 = elem_by_custom_rule.items.begin();
						it4 != elem_by_custom_rule.items.end(); ++it4) {
						Element* elem = *it4;
						me::OutputDescriptor eo_descriptor;
						Response* eo_response = 0;
						me::OutputDescriptorStream eo_stream(&eo_descriptor);

						if (do_all_sections) {
							std::map<int, std::vector<std::pair<int, int> > >::const_iterator it_info =
								m_data->elem_ngauss_nfiber_info.find(elem->getTag());
							if (it_info != m_data->elem_ngauss_nfiber_info.end()) {
								CompositeResponse* sec_comp_response = new CompositeResponse();
								int num_sec_responses = 0;
								for (size_t section_id = 0; section_id < it_info->second.size(); ++section_id) {
									std::stringstream ss_sid; ss_sid << section_id + 1;
									std::string s_sid = ss_sid.str();
									argv[section_id_placeholder_index] = s_sid.c_str();
									Response* sec_response = 0;
									if (do_all_fibers) {
										CompositeResponse* fib_comp_response = new CompositeResponse();
										int num_fib_responses = 0;
										int fiber_base_id = it_info->second[section_id].first;
										int num_fibers = it_info->second[section_id].second;
										for (int fiber_id = 0; fiber_id < num_fibers; ++fiber_id) {
											std::string s_fid = mpcolns::utils::strings::to_string(fiber_id + fiber_base_id);
											argv[fiber_id_placeholder_index] = s_fid.c_str();
											bool was_valid_before = (eo_stream.error_code == me::OutputDescriptorStream::ERROR_CODE_OK);
											Response* fib_response = elem->setResponse(argv, argc, eo_stream);
											if (fib_response) {
												num_fib_responses = fib_comp_response->addResponse(fib_response);
											}
											else {
												if (was_valid_before && (eo_stream.error_code == me::OutputDescriptorStream::ERROR_CODE_SECTION_AFTER_FIBER)) {
													eo_descriptor.fixSectionAfterFiberDueToFiberOutputFail();
													eo_stream.error_code = me::OutputDescriptorStream::ERROR_CODE_OK;
												}
											}
										}
										if (num_fib_responses == 0)
											delete fib_comp_response;
										else
											sec_response = fib_comp_response;
									}
									else {
										sec_response = elem->setResponse(argv, argc, eo_stream);
									}
									if (sec_response)
										num_sec_responses = sec_comp_response->addResponse(sec_response);
								}
								if (num_sec_responses == 0)
									delete sec_comp_response;
								else
									eo_response = sec_comp_response;
							}
						}
						else if (do_all_materials) {
							CompositeResponse* mat_comp_response = new CompositeResponse();
							int num_mat_responses = 0;
							int material_id = 0;
							while (true) {
								material_id++;
								if (material_id > MPCO_MAX_TRIAL_NSEC)
									break;
								std::stringstream ss_mid; ss_mid << material_id;
								std::string s_mid = ss_mid.str();
								argv[material_id_placeholder_index] = s_mid.c_str();
								Response* mat_response = elem->setResponse(argv, argc, eo_stream);
								if (mat_response)
									num_mat_responses = mat_comp_response->addResponse(mat_response);
								else
									break;
							}
							if (num_mat_responses == 0)
								delete mat_comp_response;
							else
								eo_response = mat_comp_response;
						}
						else {
							eo_response = elem->setResponse(argv, argc, eo_stream);
						}
						eo_stream.finalizeSetResponse();
						if (do_all_fibers)
							eo_descriptor.purge();

						if (eo_response && (eo_stream.error_code == me::OutputDescriptorStream::ERROR_CODE_OK)) {
							me::OutputDescriptorHeader header(eo_descriptor.makeHeader());
							header.workaroundForSizeInconsistency(eo_response->getInformation().getData().Size());
							header.workaroundForDuplicatedComponents();
							me::OutputResponseCollection& eo_by_header = eo_by_custom_rule.items[header];
							if (eo_by_header.is_new) {
								std::stringstream ss_buf;
								ss_buf << elem_by_tag.class_tag << "-" << elem_by_tag.class_name
									<< "[" << elem_by_rule.int_rule_type << ":"
									<< elem_by_custom_rule.custom_int_rule_index << ":"
									<< header_local_index++ << "]";
								eo_by_header.dir_name = ss_buf.str();
								eo_by_header.is_new = false;
							}
							eo_by_header.items.push_back(me::OutputResponse(elem, eo_response));
							m_data->elemental_responses.push_back(eo_response);
						}
					}
				}
			}
			if (standard_section_keyword_modified) {
				argv[0] = aux_section_keyword_standard.c_str();
				standard_section_keyword_modified = false;
			}
		}
		delete[] argv;
	}

	// Now build one ElementResultSource + StreamingSink per (request, class, rule,
	// custom-rule, header) bucket. The bucket references live in
	// m_data->elemental_recorders, kept alive for the source lifetime.
	for (size_t cnt_request = 0; cnt_request < m_data->elemental_recorders.size(); ++cnt_request) {
		me::ResultRecorder& recorder = m_data->elemental_recorders[cnt_request];
		const std::vector<std::string>& request = recorder.result_request;
		for (me::ResultRecorder::collection_type::iterator it1 = recorder.response_map.begin();
			it1 != recorder.response_map.end(); ++it1) {
			me::OutputWithSameClassTagCollection& eo_by_tag = it1->second;
			for (me::OutputWithSameClassTagCollection::collection_type::iterator it2 = eo_by_tag.items.begin();
				it2 != eo_by_tag.items.end(); ++it2) {
				me::OutputWithSameIntRuleCollection& eo_by_rule = it2->second;
				for (me::OutputWithSameIntRuleCollection::collection_type::iterator it3 = eo_by_rule.items.begin();
					it3 != eo_by_rule.items.end(); ++it3) {
					me::OutputWithSameCustomIntRuleCollection& eo_by_custom_rule = it3->second;
					for (me::OutputWithSameCustomIntRuleCollection::collection_type::iterator it4 = eo_by_custom_rule.items.begin();
						it4 != eo_by_custom_rule.items.end(); ++it4) {
						const me::OutputDescriptorHeader& header = it4->first;
						me::OutputResponseCollection& bucket = it4->second;
						if (bucket.items.empty())
							continue;
						private_data::ElemChannel ech;
						ech.source = new mpcolns::ElementResultSource(request, header, bucket);
						ech.sink = new mpcolns::StreamingSink(mpcolns::ResultFamily::OnElements);
						ech.header = &header;
						ech.column_map_written = false;
						m_data->elem_channels.push_back(ech);
						// The element source name is "<display>/<bucket>" (nested).
						// The StreamingSink creates that as a 2-level path under
						// ON_ELEMENTS, so the intermediate <display> group ("force",
						// "stress", ...) must already exist. Pre-create it here once.
						{
							std::stringstream ss_disp_grp;
							ss_disp_grp << "MODEL_STAGE[" << info.current_model_stage_id
								<< "]/RESULTS/ON_ELEMENTS/" << m_data->elem_channels.back().source->schema().display_name;
							if (H5Lexists(info.h_file_id, ss_disp_grp.str().c_str(), H5P_DEFAULT) <= 0) {
								hid_t h_disp = mpcolns::h5::group::create(
									info.h_file_id, ss_disp_grp.str().c_str(), H5P_DEFAULT,
									info.h_group_proplist, H5P_DEFAULT);
								mpcolns::h5::group::close(h_disp);
							}
						}
					}
				}
			}
		}
	}

	return 0;
}

/* ===================================================================== */
/* per-step recording                                                    */
/* ===================================================================== */

int MPCORecorderLadruno::recordResultsOnNodes()
{
	mpcolns::mpco::ProcessInfo& info = m_data->info;

	if (m_data->nodes.empty() || m_data->node_channels.empty())
		return 0;

	// EIGEN/MODES detection (frozen recordResultsOnNodes preamble).
	info.record_eigen_on_this_step = false;
	int num_eigen = *OPS_GetNumEigen();
	if (num_eigen > 0) {
		bool eigen_requested = false;
		for (size_t i = 0; i < m_data->node_channels.size(); ++i) {
			if (m_data->node_channels[i].is_modes) { eigen_requested = true; break; }
		}
		if (eigen_requested) {
			double eigen_set_time = info.domain->getTimeEigenvaluesSet();
			const Vector& new_eigenvalues = info.domain->getEigenvalues();
			if (info.eigen_first_initialization_done) {
				if (std::abs(eigen_set_time - info.eigen_last_time_set) > 1.0e-10) {
					info.record_eigen_on_this_step = true;
				}
				else if (num_eigen != info.eigen_last_values.Size()) {
					info.record_eigen_on_this_step = true;
				}
				else {
					for (int i = 0; i < num_eigen; ++i) {
						if (std::abs(new_eigenvalues[i] - info.eigen_last_values[i]) > 1.0e-10) {
							info.record_eigen_on_this_step = true;
							break;
						}
					}
				}
			}
			else {
				info.record_eigen_on_this_step = true;
				info.eigen_first_initialization_done = true;
			}
			if (info.record_eigen_on_this_step) {
				info.eigen_last_time_set = eigen_set_time;
				info.eigen_last_values = new_eigenvalues;
			}
		}
	}

	// Drive each channel. Group by reaction flag (frozen previous_reac_type),
	// priming Domain::calculateNodalReactions(flag) only when the flag changes.
	int previous_reac_type = -1;
	std::vector<double> buffer;
	for (size_t i = 0; i < m_data->node_channels.size(); ++i) {
		private_data::NodeChannel& ch = m_data->node_channels[i];

		int curr_reac_type = ch.reaction_flag;
		if (curr_reac_type != previous_reac_type) {
			if (curr_reac_type > -1 && curr_reac_type < 3)
				info.domain->calculateNodalReactions(curr_reac_type);
			previous_reac_type = curr_reac_type;
		}

		if (ch.is_modes) {
			// Modes special path: gated on record_eigen_on_this_step; drive the
			// mode loop, writing one MODE_<k> per mode with LAMBDA/OMEGA/...
			if (!info.record_eigen_on_this_step)
				continue;
			recordModeChannel((int)i);
		}
		else {
			ch.source->evaluate(info, buffer);
			ch.sink->accept(info, *ch.source, buffer);
		}
	}

	return 0;
}

/* Drive the mode loop for a ModesOfVibration(Rotational) source and write one
   MODE_<k> dataset per mode under STEP_<commitTag>, mirroring the frozen
   ResultRecorderModesOfVibration::record(). The StreamingSink writes one DATA
   group per source; modes need a STEP group containing MODE_<k> datasets, so we
   write that structure directly through the h5 wrapper here. */
void MPCORecorderLadruno::recordModeChannel(int node_channel_index)
{
	mpcolns::mpco::ProcessInfo& info = m_data->info;
	private_data::NodeChannel& ch = m_data->node_channels[(size_t)node_channel_index];
	mpcolns::ModesOfVibrationSource* msrc =
		dynamic_cast<mpcolns::ModesOfVibrationSource*>(ch.source);
	if (msrc == 0)
		return;

	const mpcolns::ResultSchema& schema = msrc->schema();
	if (schema.num_components < 1)
		return;

	// Ensure the result group + ID + DATA exist (StreamingSink::begin path).
	ch.sink->begin(info, *msrc);

	// Resolve the result + DATA group: MODEL_STAGE/RESULTS/ON_NODES/<name>/DATA
	std::stringstream ss_data;
	ss_data << "MODEL_STAGE[" << info.current_model_stage_id
		<< "]/RESULTS/ON_NODES/" << schema.name << "/DATA";
	// create the STEP group under DATA
	std::stringstream ss_step;
	ss_step << ss_data.str() << "/STEP_" << info.current_time_step_id;
	hid_t h_gp_step = mpcolns::h5::group::create(
		info.h_file_id, ss_step.str().c_str(), H5P_DEFAULT,
		info.h_group_proplist, H5P_DEFAULT);
	mpcolns::h5::attribute::write(h_gp_step, "STEP", info.current_time_step_id);
	mpcolns::h5::attribute::write(h_gp_step, "TIME", info.current_time_step);

	const std::vector<int>& ids = msrc->ids();
	const size_t n_ids = ids.size();
	const size_t n_comp = (size_t)schema.num_components;

	std::vector<double> buffer;
	int nmodes = msrc->numModes(info);
	for (int k = 0; k < nmodes; ++k) {
		msrc->setCurrentMode(k);
		msrc->evaluate(info, buffer);
		double lambda, omega, freq, period;
		msrc->modeInfo(info, k, lambda, omega, freq, period);
		std::stringstream ss_mode; ss_mode << "MODE_" << k;
		hid_t h_dset = mpcolns::h5::dataset::createAndWrite(
			h_gp_step, ss_mode.str().c_str(), buffer, n_ids, n_comp);
		mpcolns::h5::attribute::write(h_dset, "MODE", k);
		mpcolns::h5::attribute::write(h_dset, "LAMBDA", lambda);
		mpcolns::h5::attribute::write(h_dset, "OMEGA", omega);
		mpcolns::h5::attribute::write(h_dset, "FREQUENCY", freq);
		mpcolns::h5::attribute::write(h_dset, "PERIOD", period);
		mpcolns::h5::dataset::close(h_dset);
	}

	mpcolns::h5::group::close(h_gp_step);
}

int MPCORecorderLadruno::recordResultsOnElements()
{
	mpcolns::mpco::ProcessInfo& info = m_data->info;
	if (m_data->elem_channels.empty())
		return 0;

	std::vector<double> buffer;
	for (size_t i = 0; i < m_data->elem_channels.size(); ++i) {
		private_data::ElemChannel& ech = m_data->elem_channels[i];
		ech.source->evaluate(info, buffer);
		ech.sink->accept(info, *ech.source, buffer);
		// Write the COLUMN_MAP metadata once (after the sink created the result
		// group on first accept()).
		if (!ech.column_map_written) {
			writeElementColumnMap((int)i);
			ech.column_map_written = true;
		}
	}
	return 0;
}

/* ===================================================================== */
/* domain / region sources (ADR D8 energy balance)                       */
/* ===================================================================== */

int MPCORecorderLadruno::initDomainSources()
{
	mpcolns::mpco::ProcessInfo& info = m_data->info;

	// whole-model energy -> ON_DOMAIN/energyBalance
	if (m_data->energy_requested) {
		private_data::DomainChannel ch;
		ch.source = new mpcolns::EnergyBalanceSource(info);
		ch.sink = new mpcolns::StreamingSink(mpcolns::ResultFamily::OnDomain);
		m_data->domain_channels.push_back(ch);
	}

	// per-region energy -> ON_REGIONS/energyBalance
	if (!m_data->energy_region_tags.empty()) {
		private_data::DomainChannel ch;
		ch.source = new mpcolns::EnergyBalanceSource(info, m_data->energy_region_tags);
		ch.sink = new mpcolns::StreamingSink(mpcolns::ResultFamily::OnRegions);
		m_data->domain_channels.push_back(ch);
	}

	return 0;
}

int MPCORecorderLadruno::recordResultsOnDomain()
{
	mpcolns::mpco::ProcessInfo& info = m_data->info;

	std::vector<double> buffer;
	for (size_t i = 0; i < m_data->domain_channels.size(); ++i) {
		private_data::DomainChannel& ch = m_data->domain_channels[i];
		ch.source->evaluate(info, buffer);
		ch.sink->accept(info, *ch.source, buffer);
	}
	return 0;
}

/* Write the COLUMN_MAP child group + SECTION_MAP under the element result group
   (schema §7.2). Called once, after the StreamingSink has created
   ON_ELEMENTS/<display>/<bucket>. */
void MPCORecorderLadruno::writeElementColumnMap(int elem_channel_index)
{
	mpcolns::mpco::ProcessInfo& info = m_data->info;
	private_data::ElemChannel& ech = m_data->elem_channels[(size_t)elem_channel_index];
	const mpcolns::mpco::element::OutputDescriptorHeader& header = *ech.header;
	const mpcolns::ResultSchema& schema = ech.source->schema();

	// Open the bucket result group: ON_ELEMENTS/<display>/<bucket>
	std::stringstream ss_grp;
	ss_grp << "MODEL_STAGE[" << info.current_model_stage_id
		<< "]/RESULTS/ON_ELEMENTS/" << schema.name;
	hid_t h_grp = H5Gopen2(info.h_file_id, ss_grp.str().c_str(), H5P_DEFAULT);
	if (h_grp < 0)
		return;

	mpcolns::h5::attribute::write(h_grp, "NUM_COLUMNS", header.num_columns);

	// COLUMN_MAP structured metadata (schema §7.2). One row per descriptor block.
	const size_t k = header.num_components.size();
	std::vector<int> levels(k, 0);   // block path depth code = last element of path
	std::vector<int> gauss_id(k, -1);
	std::vector<int> section_tag(k, -1); // section tag not tracked at header level: -1
	std::vector<int> fiber_id(k, -1);    // fiber index not tracked at header level: -1
	std::vector<int> num_comp(k, 0);
	std::vector<int> multiplicity(k, 0);
	std::stringstream ss_comp_names;
	for (size_t i = 0; i < k; ++i) {
		// LEVELS: deepest path code in this block's components_path
		int depth_code = 0;
		if (i < header.components_path.size() && !header.components_path[i].empty())
			depth_code = header.components_path[i].back();
		levels[i] = depth_code;
		gauss_id[i] = (i < header.gauss_id.size()) ? header.gauss_id[i] : -1;
		num_comp[i] = header.num_components[i];
		multiplicity[i] = (i < header.multiplicity.size()) ? header.multiplicity[i] : 1;
		// COMP_NAMES: one line per block, CSV of component names
		if (i > 0) ss_comp_names << "\n";
		if (i < header.components.size()) {
			const std::vector<std::string>& cn = header.components[i];
			for (size_t j = 0; j < cn.size(); ++j) {
				if (j > 0) ss_comp_names << ",";
				ss_comp_names << cn[j];
			}
		}
	}

	hid_t h_cm = mpcolns::h5::group::create(
		h_grp, "COLUMN_MAP", H5P_DEFAULT, info.h_group_proplist, H5P_DEFAULT);
	if (k > 0) {
		hid_t d;
		d = mpcolns::h5::dataset::createAndWrite(h_cm, "LEVELS", levels, k, 1); mpcolns::h5::dataset::close(d);
		d = mpcolns::h5::dataset::createAndWrite(h_cm, "GAUSS_ID", gauss_id, k, 1); mpcolns::h5::dataset::close(d);
		d = mpcolns::h5::dataset::createAndWrite(h_cm, "SECTION_TAG", section_tag, k, 1); mpcolns::h5::dataset::close(d);
		d = mpcolns::h5::dataset::createAndWrite(h_cm, "FIBER_ID", fiber_id, k, 1); mpcolns::h5::dataset::close(d);
		d = mpcolns::h5::dataset::createAndWrite(h_cm, "NUM_COMP", num_comp, k, 1); mpcolns::h5::dataset::close(d);
		d = mpcolns::h5::dataset::createAndWrite(h_cm, "MULTIPLICITY", multiplicity, k, 1); mpcolns::h5::dataset::close(d);
	}
	mpcolns::h5::attribute::write(h_cm, "COMP_NAMES", ss_comp_names.str());
	mpcolns::h5::group::close(h_cm);

	mpcolns::h5::group::close(h_grp);
}

/* ===================================================================== */
/* parallel (Phase-3 stubs; real send/recv deferred)                     */
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
	int numdata = OPS_GetNumRemainingInputArgs();
	if (numdata < 1) {
		opserr << "MPCORecorderLadruno error: insufficient args; expected a filename\n"
		       << "  usage: recorder mpcoLadruno <file> [-N <res...>] [-NS <res...> <grad>] "
		          "[-E <res...>] [-R <regionTag>] [-G energy <regionTag...>] "
		          "[-T dt|nsteps <v>]\n";
		return 0;
	}

	const char* filename = OPS_GetString();
	numdata--;

	Domain* domain = OPS_GetDomain();
	if (domain == 0) {
		opserr << "MPCORecorderLadruno error: domain is not defined\n";
		return 0;
	}

	// parse optional arguments (mirror frozen OPS_MPCORecorder)
	mpcolns::utils::parsing::option_type curr_opt = mpcolns::utils::parsing::opt_none;
	std::vector<mpcolns::mpco::NodalResultType::Enum> nodal_results_requests;
	std::vector<int> sens_grad_indices;
	std::vector<std::vector<std::string> > elemental_results_requests;
	std::vector<std::string> tokens;
	mpcolns::mpco::OutputFrequency output_freq;
	bool has_region = false;
	int region_tag = 0;                  // -R region tag (for MODEL/SETS)
	std::set<int> node_set;
	std::set<int> elem_set;
	bool energy_requested = false;       // -G energy
	std::vector<int> energy_region_tags; // -G <regionTag...>
	int one_item = 1;

	while (numdata > 0) {
		const char* data = OPS_GetString();
		numdata--;
		if (strcmp(data, "-N") == 0) {
			curr_opt = mpcolns::utils::parsing::opt_result_on_nodes;
		}
		else if (strcmp(data, "-NS") == 0) {
			curr_opt = mpcolns::utils::parsing::opt_result_on_nodes_sens;
		}
		else if (strcmp(data, "-E") == 0) {
			curr_opt = mpcolns::utils::parsing::opt_result_on_elements;
		}
		else if (strcmp(data, "-T") == 0) {
			curr_opt = mpcolns::utils::parsing::opt_time;
			output_freq.reset();
		}
		else if (strcmp(data, "-R") == 0) {
			curr_opt = mpcolns::utils::parsing::opt_region;
			if (numdata > 0) {
				int rtag = 0;
				if (OPS_GetInt(&one_item, &rtag) != 0) {
					opserr << "MPCORecorderLadruno error: option -R (region) requires an int region tag\n";
					return 0;
				}
				MeshRegion* region = domain->getRegion(rtag);
				if (region == 0) {
					opserr << "MPCORecorderLadruno error: region " << rtag << " is null\n";
					return 0;
				}
				const ID& node_ids = region->getNodes();
				const ID& elem_ids = region->getElements();
				for (int i = 0; i < node_ids.Size(); ++i)
					node_set.insert(node_ids(i));
				for (int i = 0; i < elem_ids.Size(); ++i)
					elem_set.insert(elem_ids(i));
				has_region = true;
				region_tag = rtag; // retained for MODEL/SETS self-description
				numdata--;
			}
			else {
				opserr << "MPCORecorderLadruno error: option -R (region) requires an int region tag\n";
				return 0;
			}
			has_region = true;
		}
		else if (strcmp(data, "-G") == 0) {
			// -G energy <regionTag...> : whole-model energy balance (ON_DOMAIN)
			// + optional per-region balance (ON_REGIONS). ADR D8.
			// Consume eagerly: the "energy" keyword is a real string arg, but the
			// region tags are integers -> they MUST be read with OPS_GetIntInput,
			// not OPS_GetString (which returns "" for a numeric OpenSeesPy arg,
			// so atoi() would silently yield 0). Same idiom as -node.
			curr_opt = mpcolns::utils::parsing::opt_global;
			if (numdata > 0) {
				const char* gkind = OPS_GetString();
				numdata--;
				if (strcmp(gkind, "energy") == 0) {
					energy_requested = true;
				}
				else {
					// tolerate a tag passed as a string (Tcl path)
					int t = atoi(gkind);
					if (t != 0 && domain->getRegion(t) != 0)
						energy_region_tags.push_back(t);
					else
						opserr << "MPCORecorderLadruno warning: -G expects 'energy' "
						          "[regionTag...]; ignoring token (" << gkind << ")\n";
				}
			}
			// 0+ integer region tags follow; stop at the first non-int (an option
			// like -T) and hand it back to the outer parse loop.
			while (numdata > 0) {
				int tag = 0; int n = 1;
				if (OPS_GetIntInput(&n, &tag) < 0) {
					OPS_ResetCurrentInputArg(-1);
					break;
				}
				numdata--;
				MeshRegion* gregion = domain->getRegion(tag);
				if (gregion == 0)
					opserr << "MPCORecorderLadruno warning: -G region " << tag
					       << " not found in the domain; skipping its energy block\n";
				else
					energy_region_tags.push_back(tag);
			}
		}
		else {
			switch (curr_opt) {
			case mpcolns::utils::parsing::opt_result_on_nodes: {
				if (strcmp(data, "displacement") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::Displacement);
				else if (strcmp(data, "rotation") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::Rotation);
				else if (strcmp(data, "velocity") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::Velocity);
				else if (strcmp(data, "angularVelocity") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::AngularVelocity);
				else if (strcmp(data, "acceleration") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::Acceleration);
				else if (strcmp(data, "angularAcceleration") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::AngularAcceleration);
				else if (strcmp(data, "reactionForce") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::ReactionForce);
				else if (strcmp(data, "reactionMoment") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::ReactionMoment);
				else if (strcmp(data, "reactionForceIncludingInertia") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::ReactionForceIncludingInertia);
				else if (strcmp(data, "reactionMomentIncludingInertia") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::ReactionMomentIncludingInertia);
				else if (strcmp(data, "rayleighForce") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::RayleighForce);
				else if (strcmp(data, "rayleighMoment") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::RayleighMoment);
				else if (strcmp(data, "unbalancedForce") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::UnbalancedForce);
				else if (strcmp(data, "unbalancedMoment") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::UnbalancedMoment);
				else if (strcmp(data, "unbalancedForceIncludingInertia") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::UnbalancedForceIncludingInertia);
				else if (strcmp(data, "unbalancedMomentIncludingInertia") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::UnbalancedMomentIncludingInertia);
				else if (strcmp(data, "pressure") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::Pressure);
				else if (strcmp(data, "modesOfVibration") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::ModesOfVibration);
				else if (strcmp(data, "modesOfVibrationRotational") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::ModesOfVibrationRotational);
				else {
					opserr << "MPCORecorderLadruno error: option -N with unknown result type (" << data << ")\n";
					return 0;
				}
				sens_grad_indices.push_back(0); // placeholder
				break;
			}
			case mpcolns::utils::parsing::opt_result_on_nodes_sens: {
				if (strcmp(data, "displacementSensitivity") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::DisplacementSensitivity);
				else if (strcmp(data, "rotationSensitivity") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::RotationSensitivity);
				else if (strcmp(data, "velocitySensitivity") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::VelocitySensitivity);
				else if (strcmp(data, "angularVelocitySensitivity") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::AngularVelocitySensitivity);
				else if (strcmp(data, "accelerationSensitivity") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::AccelerationSensitivity);
				else if (strcmp(data, "angularAccelerationSensitivity") == 0)
					nodal_results_requests.push_back(mpcolns::mpco::NodalResultType::AngularAccelerationSensitivity);
				else {
					opserr << "MPCORecorderLadruno error: option -NS with unknown result type (" << data << ")\n";
					return 0;
				}
				if (numdata > 0) {
					int grad_index;
					if (OPS_GetInt(&one_item, &grad_index) != 0) {
						opserr << "MPCORecorderLadruno error: option -NS requires an int sensitivity parameter index\n";
						return 0;
					}
					numdata--;
					sens_grad_indices.push_back(grad_index);
				}
				else {
					opserr << "MPCORecorderLadruno error: option -NS requires a sensitivity parameter index\n";
					return 0;
				}
				break;
			}
			case mpcolns::utils::parsing::opt_result_on_elements: {
				std::string temp(data);
				mpcolns::utils::strings::split(temp, '.', tokens, true);
				if (tokens.size() > 0)
					elemental_results_requests.push_back(tokens);
				break;
			}
			case mpcolns::utils::parsing::opt_time: {
				if (strcmp(data, "dt") == 0) {
					output_freq.type = mpcolns::mpco::OutputFrequency::DeltaTime;
					output_freq.nsteps = 1;
					if (numdata > 0) {
						if (OPS_GetDouble(&one_item, &output_freq.dt) != 0) {
							opserr << "MPCORecorderLadruno error: invalid double argument for the delta time\n";
							return 0;
						}
						if (output_freq.dt < 0.0) output_freq.dt = 0.0;
						numdata--;
					}
					else {
						opserr << "MPCORecorderLadruno error: option -T dt requires a delta time argument\n";
						return 0;
					}
				}
				else if (strcmp(data, "nsteps") == 0) {
					output_freq.type = mpcolns::mpco::OutputFrequency::NumberOfSteps;
					output_freq.dt = 0.0;
					if (numdata > 0) {
						if (OPS_GetInt(&one_item, &output_freq.nsteps) != 0) {
							opserr << "MPCORecorderLadruno error: invalid int argument for the number of steps\n";
							return 0;
						}
						if (output_freq.nsteps < 1) output_freq.nsteps = 1;
						numdata--;
					}
					else {
						opserr << "MPCORecorderLadruno error: option -T nsteps requires a number-of-steps argument\n";
						return 0;
					}
				}
				else {
					opserr << "MPCORecorderLadruno error: option -T with unknown frequency type (" << data << ")\n";
					return 0;
				}
				break;
			}
			case mpcolns::utils::parsing::opt_global: {
				// -G energy           -> whole-model energy balance (ON_DOMAIN)
				// -G energy <tag...>  -> additionally per-region (ON_REGIONS)
				if (strcmp(data, "energy") == 0) {
					energy_requested = true;
				}
				else {
					// interpret as a region tag (same int-parse style as -R)
					int tag = atoi(data);
					MeshRegion* region = domain->getRegion(tag);
					if (region == 0) {
						opserr << "MPCORecorderLadruno warning: option -G region " << tag
						       << " not found in the domain; skipping its energy block\n";
					}
					else {
						energy_region_tags.push_back(tag);
					}
				}
				break;
			}
			default: {
				opserr << "MPCORecorderLadruno error: unknown arg with option none " << data << "\n";
				return 0;
			}
			}
		}
	}

	MPCORecorderLadruno* recorder = new MPCORecorderLadruno();
	recorder->m_data->filename = filename;
	recorder->m_data->output_freq = output_freq;
	recorder->m_data->nodal_results_requests.swap(nodal_results_requests);
	recorder->m_data->sens_grad_indices.swap(sens_grad_indices);
	recorder->m_data->elemental_results_requests.swap(elemental_results_requests);
	recorder->m_data->has_region = has_region;
	recorder->m_data->region_tag = region_tag;
	for (std::set<int>::const_iterator it = node_set.begin(); it != node_set.end(); ++it)
		recorder->m_data->node_set.push_back(*it);
	for (std::set<int>::const_iterator it = elem_set.begin(); it != elem_set.end(); ++it)
		recorder->m_data->elem_set.push_back(*it);
	recorder->m_data->energy_requested = energy_requested;
	recorder->m_data->energy_region_tags.swap(energy_region_tags);
	return recorder;
}
