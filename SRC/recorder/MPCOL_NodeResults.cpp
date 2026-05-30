/* ********************************************************************** **
**  MPCO_Ladruno recorder — modular sibling of MPCORecorder (frozen).     **
**  MPCOL_NodeResults.cpp — NODAL ResultSource implementations (Phase 2). **
**                                                                        **
**  Faithful port of the frozen `mpco::node` ResultRecorder subclasses    **
**  (MPCORecorder.cpp lines 1478-2594). Every evaluate() reproduces the   **
**  corresponding frozen bufferResponse() byte-for-byte (1e-12 parity);   **
**  every schema() reproduces the frozen ctor metadata. The group name    **
**  string built in each schema is identical to the frozen m_result_name  **
**  (MODEL_STAGE[<id>]/RESULTS/ON_NODES/<GROUP>).                          **
** ********************************************************************** */

#include "MPCOL_NodeResults.h"

// OPS_GetNumEigen() lives in the element API (frozen recorder uses it directly).
#include "elementAPI.h"
#include "OPS_Globals.h"   // opserr
#include "Matrix.h"        // Node::getEigenvectors() returns const Matrix&

#include <cmath>
#include <sstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace mpcol {

	// Build the frozen result-group name: MODEL_STAGE[<id>]/RESULTS/ON_NODES/<group>
	static std::string makeResultName(const mpco::ProcessInfo& info, const std::string& group)
	{
		std::stringstream ss;
		ss << "MODEL_STAGE[" << info.current_model_stage_id << "]/RESULTS/ON_NODES/" << group;
		return ss.str();
	}

	/* ===================================================================== */
	/* NodeResultSource base                                                 */
	/* ===================================================================== */

	NodeResultSource::NodeResultSource(const mpco::ProcessInfo& info)
		: m_ndim(info.num_dimensions)
	{
		// Gather all nodes from the domain (frozen recorder default node set,
		// matching MPCORecorderLadruno::writeModel()'s NodeIter walk).
		if (info.domain != 0) {
			NodeIter& nit = info.domain->getNodes();
			Node* nptr = 0;
			while ((nptr = nit()) != 0)
				m_ids.push_back(nptr->getTag());
		}
	}

	NodeResultSource::NodeResultSource(const mpco::ProcessInfo& info, const std::vector<int>& node_ids)
		: m_ndim(info.num_dimensions)
		, m_ids(node_ids)
	{
	}

	Node* NodeResultSource::getNode(const mpco::ProcessInfo& info, int tag) const
	{
		if (info.domain == 0)
			return 0;
		return info.domain->getNode(tag);
	}

	/* ===================================================================== */
	/* DisplacementSource                                                    */
	/* ===================================================================== */

	DisplacementSource::DisplacementSource(const mpco::ProcessInfo& info)
		: NodeResultSource(info) { build(); }
	DisplacementSource::DisplacementSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids) { build(); }

	void DisplacementSource::build()
	{
		m_schema.name = "DISPLACEMENT";
		m_schema.display_name = "Displacement";
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = "Ux";
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = "Ux,Uy";
			m_schema.num_components = 2;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = "Ux,Uy,Uz";
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "L";
		m_schema.description = "Nodal displacement field";
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void DisplacementSource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = getNode(info, m_ids[i]);
			if (n) utils::misc::bufferNodeResponseVec3u(i, m_ndim, n->getTrialDisp(), buffer);
		}
	}

	/* ===================================================================== */
	/* RotationSource                                                        */
	/* ===================================================================== */

	RotationSource::RotationSource(const mpco::ProcessInfo& info)
		: NodeResultSource(info) { build(); }
	RotationSource::RotationSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids) { build(); }

	void RotationSource::build()
	{
		m_schema.name = "ROTATION";
		m_schema.display_name = "Rotation";
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = "Rz";
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = "Rx,Ry,Rz";
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "A";
		m_schema.description = "Nodal rotation field";
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void RotationSource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = getNode(info, m_ids[i]);
			if (n) utils::misc::bufferNodeResponseVec3r(i, m_ndim, n->getTrialDisp(), buffer);
		}
	}

	/* ===================================================================== */
	/* VelocitySource                                                        */
	/* ===================================================================== */

	VelocitySource::VelocitySource(const mpco::ProcessInfo& info)
		: NodeResultSource(info) { build(); }
	VelocitySource::VelocitySource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids) { build(); }

	void VelocitySource::build()
	{
		m_schema.name = "VELOCITY";
		m_schema.display_name = "Velocity";
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = "Vx";
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = "Vx,Vy";
			m_schema.num_components = 2;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = "Vx,Vy,Vz";
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "L/t";
		m_schema.description = "Nodal velocity field";
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void VelocitySource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = getNode(info, m_ids[i]);
			if (n) utils::misc::bufferNodeResponseVec3u(i, m_ndim, n->getTrialVel(), buffer);
		}
	}

	/* ===================================================================== */
	/* AngularVelocitySource                                                 */
	/* ===================================================================== */

	AngularVelocitySource::AngularVelocitySource(const mpco::ProcessInfo& info)
		: NodeResultSource(info) { build(); }
	AngularVelocitySource::AngularVelocitySource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids) { build(); }

	void AngularVelocitySource::build()
	{
		m_schema.name = "ANGULAR_VELOCITY";
		m_schema.display_name = "Angular Velocity";
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = "RVz";
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = "RVx,RVy,RVz";
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "L/t*A";
		m_schema.description = "Nodal angular velocity field";
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void AngularVelocitySource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = getNode(info, m_ids[i]);
			if (n) utils::misc::bufferNodeResponseVec3r(i, m_ndim, n->getTrialVel(), buffer);
		}
	}

	/* ===================================================================== */
	/* AccelerationSource                                                    */
	/* ===================================================================== */

	AccelerationSource::AccelerationSource(const mpco::ProcessInfo& info)
		: NodeResultSource(info) { build(); }
	AccelerationSource::AccelerationSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids) { build(); }

	void AccelerationSource::build()
	{
		m_schema.name = "ACCELERATION";
		m_schema.display_name = "Acceleration";
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = "Ax";
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = "Ax,Ay";
			m_schema.num_components = 2;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = "Ax,Ay,Az";
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "L/t^2";
		m_schema.description = "Nodal acceleration field";
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void AccelerationSource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = getNode(info, m_ids[i]);
			if (n) utils::misc::bufferNodeResponseVec3u(i, m_ndim, n->getTrialAccel(), buffer);
		}
	}

	/* ===================================================================== */
	/* AngularAccelerationSource                                             */
	/* ===================================================================== */

	AngularAccelerationSource::AngularAccelerationSource(const mpco::ProcessInfo& info)
		: NodeResultSource(info) { build(); }
	AngularAccelerationSource::AngularAccelerationSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids) { build(); }

	void AngularAccelerationSource::build()
	{
		m_schema.name = "ANGULAR_ACCELERATION";
		m_schema.display_name = "Angular Acceleration";
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = "RAz";
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = "RAx,RAy,RAz";
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "L/t^2*A";
		m_schema.description = "Nodal angular acceleration field";
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void AngularAccelerationSource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = getNode(info, m_ids[i]);
			if (n) utils::misc::bufferNodeResponseVec3r(i, m_ndim, n->getTrialAccel(), buffer);
		}
	}

	/* ===================================================================== */
	/* PressureSource                                                        */
	/* ===================================================================== */

	PressureSource::PressureSource(const mpco::ProcessInfo& info)
		: NodeResultSource(info) { build(); }
	PressureSource::PressureSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids) { build(); }

	void PressureSource::build()
	{
		m_schema.name = "PRESSURE";
		m_schema.display_name = "Pressure";
		m_schema.components_csv = "p";
		m_schema.num_components = 1;
		m_schema.data_type = mpco::ResultDataType::Scalar;
		m_schema.dimension = "F/L^2";
		m_schema.description = "Nodal pressure field";
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void PressureSource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		// Pressure is the extra velocity DOF: idx 2 (2D u-p) or idx 3 (3D u-p).
		for (size_t i = 0; i < m_ids.size(); i++) {
			double pressure = 0.0;
			Node* inode = getNode(info, m_ids[i]);
			if (inode) {
				const Vector& vel = inode->getTrialVel();
				if (inode->getCrds().Size() == 2) {
					if (vel.Size() == 3) {
						pressure = vel[2];
					}
				}
				else if (inode->getCrds().Size() == 3) {
					if (vel.Size() == 4) {
						pressure = vel[3];
					}
				}
			}
			buffer[i] = pressure;
		}
	}

	/* ===================================================================== */
	/* ReactionForceSource (+ inertia/Rayleigh variants)                     */
	/* ===================================================================== */

	ReactionForceSource::ReactionForceSource(const mpco::ProcessInfo& info)
		: NodeResultSource(info)
	{ build("REACTION_FORCE", "Reaction Force", "Nodal reaction force field"); }
	ReactionForceSource::ReactionForceSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids)
	{ build("REACTION_FORCE", "Reaction Force", "Nodal reaction force field"); }

	void ReactionForceSource::build(const char* group, const char* display, const char* description)
	{
		m_schema.name = group;
		m_schema.display_name = display;
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = "RFx";
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = "RFx,RFy";
			m_schema.num_components = 2;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = "RFx,RFy,RFz";
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "F";
		m_schema.description = description;
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void ReactionForceSource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = getNode(info, m_ids[i]);
			if (n) utils::misc::bufferNodeResponseVec3u(i, m_ndim, n->getReaction(), buffer);
		}
	}

	ReactionForceIncInertiaSource::ReactionForceIncInertiaSource(const mpco::ProcessInfo& info)
		: ReactionForceSource(info)
	{ build("REACTION_FORCE_INCLUDING_INERTIA", "Reaction Force Including Inertia",
	        "Nodal reaction force field including inertia"); }
	ReactionForceIncInertiaSource::ReactionForceIncInertiaSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: ReactionForceSource(info, ids)
	{ build("REACTION_FORCE_INCLUDING_INERTIA", "Reaction Force Including Inertia",
	        "Nodal reaction force field including inertia"); }

	RayleighForceSource::RayleighForceSource(const mpco::ProcessInfo& info)
		: ReactionForceSource(info)
	{ build("RAYLEIGH_FORCE", "Rayleigh Force", "Nodal Rayleigh force field"); }
	RayleighForceSource::RayleighForceSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: ReactionForceSource(info, ids)
	{ build("RAYLEIGH_FORCE", "Rayleigh Force", "Nodal Rayleigh force field"); }

	/* ===================================================================== */
	/* ReactionMomentSource (+ inertia/Rayleigh variants)                    */
	/* ===================================================================== */

	ReactionMomentSource::ReactionMomentSource(const mpco::ProcessInfo& info)
		: NodeResultSource(info)
	{ build("REACTION_MOMENT", "Reaction Moment", "Nodal reaction moment field"); }
	ReactionMomentSource::ReactionMomentSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids)
	{ build("REACTION_MOMENT", "Reaction Moment", "Nodal reaction moment field"); }

	void ReactionMomentSource::build(const char* group, const char* display, const char* description)
	{
		m_schema.name = group;
		m_schema.display_name = display;
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = "RMz";
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = "RMx,RMy,RMz";
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "F*L";
		m_schema.description = description;
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void ReactionMomentSource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = getNode(info, m_ids[i]);
			if (n) utils::misc::bufferNodeResponseVec3r(i, m_ndim, n->getReaction(), buffer);
		}
	}

	ReactionMomentIncInertiaSource::ReactionMomentIncInertiaSource(const mpco::ProcessInfo& info)
		: ReactionMomentSource(info)
	{ build("REACTION_MOMENT_INCLUDING_INERTIA", "Reaction Moment Including Inertia",
	        "Nodal reaction moment field including inertia"); }
	ReactionMomentIncInertiaSource::ReactionMomentIncInertiaSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: ReactionMomentSource(info, ids)
	{ build("REACTION_MOMENT_INCLUDING_INERTIA", "Reaction Moment Including Inertia",
	        "Nodal reaction moment field including inertia"); }

	RayleighMomentSource::RayleighMomentSource(const mpco::ProcessInfo& info)
		: ReactionMomentSource(info)
	{ build("RAYLEIGH_MOMENT", "Rayleigh Moment", "Nodal Rayleigh moment field"); }
	RayleighMomentSource::RayleighMomentSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: ReactionMomentSource(info, ids)
	{ build("RAYLEIGH_MOMENT", "Rayleigh Moment", "Nodal Rayleigh moment field"); }

	/* ===================================================================== */
	/* UnbalancedForceSource (+ IncInertia variant)                          */
	/* ===================================================================== */

	UnbalancedForceSource::UnbalancedForceSource(const mpco::ProcessInfo& info)
		: NodeResultSource(info)
	{ build("UNBALANCED_FORCE", "Unbalanced Force", "Nodal unbalanced force field"); }
	UnbalancedForceSource::UnbalancedForceSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids)
	{ build("UNBALANCED_FORCE", "Unbalanced Force", "Nodal unbalanced force field"); }

	void UnbalancedForceSource::build(const char* group, const char* display, const char* description)
	{
		m_schema.name = group;
		m_schema.display_name = display;
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = "Fx";
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = "Fx,Fy";
			m_schema.num_components = 2;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = "Fx,Fy,Fz";
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "F";
		m_schema.description = description;
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void UnbalancedForceSource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = getNode(info, m_ids[i]);
			if (n) utils::misc::bufferNodeResponseVec3u(i, m_ndim, n->getUnbalancedLoad(), buffer);
		}
	}

	UnbalancedForceIncInertiaSource::UnbalancedForceIncInertiaSource(const mpco::ProcessInfo& info)
		: UnbalancedForceSource(info)
	{ build("UNBALANCED_FORCE_INCLUDING_INERTIA", "Unbalanced Force Including Inertia",
	        "Nodal unbalanced force field including inertia"); }
	UnbalancedForceIncInertiaSource::UnbalancedForceIncInertiaSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: UnbalancedForceSource(info, ids)
	{ build("UNBALANCED_FORCE_INCLUDING_INERTIA", "Unbalanced Force Including Inertia",
	        "Nodal unbalanced force field including inertia"); }

	void UnbalancedForceIncInertiaSource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = getNode(info, m_ids[i]);
			if (n) utils::misc::bufferNodeResponseVec3u(i, m_ndim, n->getUnbalancedLoadIncInertia(), buffer);
		}
	}

	/* ===================================================================== */
	/* UnbalancedMomentSource (+ IncInertia variant)                         */
	/* ===================================================================== */

	UnbalancedMomentSource::UnbalancedMomentSource(const mpco::ProcessInfo& info)
		: NodeResultSource(info)
	{ build("UNBALANCED_MOMENT", "Unbalanced Moment", "Nodal unbalanced moment field"); }
	UnbalancedMomentSource::UnbalancedMomentSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids)
	{ build("UNBALANCED_MOMENT", "Unbalanced Moment", "Nodal unbalanced moment field"); }

	void UnbalancedMomentSource::build(const char* group, const char* display, const char* description)
	{
		m_schema.name = group;
		m_schema.display_name = display;
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = "Mz";
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = "Mx,My,Mz";
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "F*L";
		m_schema.description = description;
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void UnbalancedMomentSource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = getNode(info, m_ids[i]);
			if (n) utils::misc::bufferNodeResponseVec3r(i, m_ndim, n->getUnbalancedLoad(), buffer);
		}
	}

	UnbalancedMomentIncInertiaSource::UnbalancedMomentIncInertiaSource(const mpco::ProcessInfo& info)
		: UnbalancedMomentSource(info)
	{ build("UNBALANCED_MOMENT_INCLUDING_INERTIA", "Unbalanced Moment Including Inertia",
	        "Nodal unbalanced moment field including inertia"); }
	UnbalancedMomentIncInertiaSource::UnbalancedMomentIncInertiaSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: UnbalancedMomentSource(info, ids)
	{ build("UNBALANCED_MOMENT_INCLUDING_INERTIA", "Unbalanced Moment Including Inertia",
	        "Nodal unbalanced moment field including inertia"); }

	void UnbalancedMomentIncInertiaSource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = getNode(info, m_ids[i]);
			if (n) utils::misc::bufferNodeResponseVec3r(i, m_ndim, n->getUnbalancedLoadIncInertia(), buffer);
		}
	}

	/* ===================================================================== */
	/* ModesOfVibrationSource (translational) — SPECIAL multi-mode case.     */
	/* ===================================================================== */

	ModesOfVibrationSource::ModesOfVibrationSource(const mpco::ProcessInfo& info)
		: NodeResultSource(info), m_current_mode(0) { build(); }
	ModesOfVibrationSource::ModesOfVibrationSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: NodeResultSource(info, ids), m_current_mode(0) { build(); }

	void ModesOfVibrationSource::build()
	{
		m_schema.name = "MODES_OF_VIBRATION(U)";
		m_schema.display_name = "Modes of Vibration (U)";
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = "Ux";
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = "Ux,Uy";
			m_schema.num_components = 2;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = "Ux,Uy,Uz";
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = "Eigenvector (translational components)";
		m_schema.result_type = mpco::ResultType::Modal;
	}

	int ModesOfVibrationSource::numModes(const mpco::ProcessInfo& /*info*/) const
	{
		return *OPS_GetNumEigen();
	}

	void ModesOfVibrationSource::modeInfo(const mpco::ProcessInfo& info, int k,
	                                      double& lambda, double& omega,
	                                      double& freq, double& period) const
	{
		// Frozen ModesOfVibration::record() per-mode scalars.
		lambda = k < info.eigen_last_values.Size() ? info.eigen_last_values[k] : 0.0;
		omega = std::sqrt(lambda);
		freq = omega / (2.0 * M_PI);
		period = 1.0 / freq;
	}

	void ModesOfVibrationSource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		// Translational eigenvector components for the CURRENT mode.
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = getNode(info, m_ids[i]);
			if (!n) continue;
			size_t j = i * m_ndim;
			const Matrix& inode_eigenvec = n->getEigenvectors();
			buffer[j] = inode_eigenvec(0, m_current_mode);
			if (m_ndim > 1) {
				buffer[j + 1] = inode_eigenvec(1, m_current_mode);
				if (m_ndim > 2) {
					buffer[j + 2] = inode_eigenvec(2, m_current_mode);
				}
			}
		}
	}

	/* ===================================================================== */
	/* ModesOfVibrationRotationalSource                                      */
	/* ===================================================================== */

	ModesOfVibrationRotationalSource::ModesOfVibrationRotationalSource(const mpco::ProcessInfo& info)
		: ModesOfVibrationSource(info) { build(); }
	ModesOfVibrationRotationalSource::ModesOfVibrationRotationalSource(const mpco::ProcessInfo& info, const std::vector<int>& ids)
		: ModesOfVibrationSource(info, ids) { build(); }

	void ModesOfVibrationRotationalSource::build()
	{
		m_schema.name = "MODES_OF_VIBRATION(R)";
		m_schema.display_name = "Modes of Vibration (R)";
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = "Rz";
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = "Rx,Ry,Rz";
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = "Eigenvector (rotational components)";
		m_schema.result_type = mpco::ResultType::Modal;
	}

	void ModesOfVibrationRotationalSource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		// Rotational eigenvector slice depends on the node's DOF count (rows>2 in
		// 2D, rows>5 in 3D); nodes without rotational DOFs stay zero.
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* n = getNode(info, m_ids[i]);
			if (!n) continue;
			const Matrix& inode_eigenvec = n->getEigenvectors();
			if (m_ndim == 2 && inode_eigenvec.noRows() > 2) {
				size_t j = i;
				buffer[j] = inode_eigenvec(2, m_current_mode);
			}
			else if (m_ndim == 3 && inode_eigenvec.noRows() > 5) {
				size_t j = i * 3;
				buffer[j] = inode_eigenvec(3, m_current_mode);
				buffer[j + 1] = inode_eigenvec(4, m_current_mode);
				buffer[j + 2] = inode_eigenvec(5, m_current_mode);
			}
		}
	}

	/* ===================================================================== */
	/* Sensitivity sources                                                   */
	/* ===================================================================== */

	DisplacementSensitivitySource::DisplacementSensitivitySource(const mpco::ProcessInfo& info, int grad)
		: NodeResultSource(info), m_grad(grad) { build(); }
	DisplacementSensitivitySource::DisplacementSensitivitySource(const mpco::ProcessInfo& info, const std::vector<int>& ids, int grad)
		: NodeResultSource(info, ids), m_grad(grad) { build(); }

	void DisplacementSensitivitySource::build()
	{
		m_schema.name = MPCO_MAKE_STRING("DISPLACEMENT_SENSITIVITY_" << m_grad);
		m_schema.display_name = MPCO_MAKE_STRING("Displacement Sensitivity d_U/d_q" << m_grad);
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = MPCO_MAKE_STRING("d_Ux/d_q" << m_grad);
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = MPCO_MAKE_STRING("d_Ux/d_q" << m_grad << ",d_Uy/d_q" << m_grad);
			m_schema.num_components = 2;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = MPCO_MAKE_STRING("d_Ux/d_q" << m_grad << ",d_Uy/d_q" << m_grad << ",d_Uz/d_q" << m_grad);
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = MPCO_MAKE_STRING("Nodal displacement sensitivity field d_U/d_q" << m_grad);
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void DisplacementSensitivitySource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
			opserr << "MPCORecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
			return;
		}
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* inode = getNode(info, m_ids[i]);
			if (!inode) continue;
			size_t j = i * m_ndim;
			buffer[j] = inode->getDispSensitivity(1, m_grad - 1);
			if (m_ndim > 1) {
				buffer[j + 1] = inode->getDispSensitivity(2, m_grad - 1);
				if (m_ndim > 2) {
					buffer[j + 2] = inode->getDispSensitivity(3, m_grad - 1);
				}
			}
		}
	}

	RotationSensitivitySource::RotationSensitivitySource(const mpco::ProcessInfo& info, int grad)
		: NodeResultSource(info), m_grad(grad) { build(); }
	RotationSensitivitySource::RotationSensitivitySource(const mpco::ProcessInfo& info, const std::vector<int>& ids, int grad)
		: NodeResultSource(info, ids), m_grad(grad) { build(); }

	void RotationSensitivitySource::build()
	{
		m_schema.name = MPCO_MAKE_STRING("ROTATION_SENSITIVITY_" << m_grad);
		m_schema.display_name = MPCO_MAKE_STRING("Rotation Sensitivity d_R/d_q" << m_grad);
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = MPCO_MAKE_STRING("d_Rz/d_q" << m_grad);
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = MPCO_MAKE_STRING("d_Rx/d_q" << m_grad << ",d_Ry/d_q" << m_grad << ",d_Rz/d_q" << m_grad);
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = MPCO_MAKE_STRING("Nodal rotation sensitivity field d_R/d_q" << m_grad);
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void RotationSensitivitySource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
			opserr << "MPCORecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
			return;
		}
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* inode = getNode(info, m_ids[i]);
			if (!inode) continue;
			const Vector& iresponse = inode->getTrialDisp(); // aux, just to check for rotational components
			if (m_ndim == 2) {
				size_t j = i; // node_counter * 1 (1 rotation component in 2D)
				if (iresponse.Size() > 2)
					buffer[j] = inode->getDispSensitivity(3, m_grad - 1);
				else
					buffer[j] = 0.0;
			}
			else if (m_ndim == 3) {
				size_t j = i * 3; // node_counter * 3 (3 rotation components in 3D)
				if (iresponse.Size() > 5) {
					buffer[j] = inode->getDispSensitivity(4, m_grad - 1);
					buffer[j + 1] = inode->getDispSensitivity(5, m_grad - 1);
					buffer[j + 2] = inode->getDispSensitivity(6, m_grad - 1);
				}
				else {
					buffer[j] = 0.0;
					buffer[j + 1] = 0.0;
					buffer[j + 2] = 0.0;
				}
			}
		}
	}

	VelocitySensitivitySource::VelocitySensitivitySource(const mpco::ProcessInfo& info, int grad)
		: NodeResultSource(info), m_grad(grad) { build(); }
	VelocitySensitivitySource::VelocitySensitivitySource(const mpco::ProcessInfo& info, const std::vector<int>& ids, int grad)
		: NodeResultSource(info, ids), m_grad(grad) { build(); }

	void VelocitySensitivitySource::build()
	{
		m_schema.name = MPCO_MAKE_STRING("VELOCITY_SENSITIVITY_" << m_grad);
		m_schema.display_name = MPCO_MAKE_STRING("Velocity Sensitivity d_V/d_q" << m_grad);
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = MPCO_MAKE_STRING("d_Vx/d_q" << m_grad);
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = MPCO_MAKE_STRING("d_Vx/d_q" << m_grad << ",d_Vy/d_q" << m_grad);
			m_schema.num_components = 2;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = MPCO_MAKE_STRING("d_Vx/d_q" << m_grad << ",d_Vy/d_q" << m_grad << ",d_Vz/d_q" << m_grad);
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = MPCO_MAKE_STRING("Nodal velocity sensitivity field d_V/d_q" << m_grad);
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void VelocitySensitivitySource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
			opserr << "MPCORecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
			return;
		}
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* inode = getNode(info, m_ids[i]);
			if (!inode) continue;
			size_t j = i * m_ndim;
			buffer[j] = inode->getVelSensitivity(1, m_grad - 1);
			if (m_ndim > 1) {
				buffer[j + 1] = inode->getVelSensitivity(2, m_grad - 1);
				if (m_ndim > 2) {
					buffer[j + 2] = inode->getVelSensitivity(3, m_grad - 1);
				}
			}
		}
	}

	AngularVelocitySensitivitySource::AngularVelocitySensitivitySource(const mpco::ProcessInfo& info, int grad)
		: NodeResultSource(info), m_grad(grad) { build(); }
	AngularVelocitySensitivitySource::AngularVelocitySensitivitySource(const mpco::ProcessInfo& info, const std::vector<int>& ids, int grad)
		: NodeResultSource(info, ids), m_grad(grad) { build(); }

	void AngularVelocitySensitivitySource::build()
	{
		m_schema.name = MPCO_MAKE_STRING("ANGULAR_VELOCITY_SENSITIVITY_" << m_grad);
		m_schema.display_name = MPCO_MAKE_STRING("Angular Velocity Sensitivity d_RV/d_q" << m_grad);
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = MPCO_MAKE_STRING("d_RVz/d_q" << m_grad);
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = MPCO_MAKE_STRING("d_RVx/d_q" << m_grad << ",d_RVy/d_q" << m_grad << ",d_RVz/d_q" << m_grad);
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = MPCO_MAKE_STRING("Nodal angular velocity sensitivity field d_RV/d_q" << m_grad);
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void AngularVelocitySensitivitySource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
			opserr << "MPCORecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
			return;
		}
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* inode = getNode(info, m_ids[i]);
			if (!inode) continue;
			const Vector& iresponse = inode->getTrialVel(); // aux, just to check for angular velocityal components
			if (m_ndim == 2) {
				size_t j = i; // node_counter * 1 (1 angular velocity component in 2D)
				if (iresponse.Size() > 2)
					buffer[j] = inode->getVelSensitivity(3, m_grad - 1);
				else
					buffer[j] = 0.0;
			}
			else if (m_ndim == 3) {
				size_t j = i * 3; // node_counter * 3 (3 angular velocity components in 3D)
				if (iresponse.Size() > 5) {
					buffer[j] = inode->getVelSensitivity(4, m_grad - 1);
					buffer[j + 1] = inode->getVelSensitivity(5, m_grad - 1);
					buffer[j + 2] = inode->getVelSensitivity(6, m_grad - 1);
				}
				else {
					buffer[j] = 0.0;
					buffer[j + 1] = 0.0;
					buffer[j + 2] = 0.0;
				}
			}
		}
	}

	AccelerationSensitivitySource::AccelerationSensitivitySource(const mpco::ProcessInfo& info, int grad)
		: NodeResultSource(info), m_grad(grad) { build(); }
	AccelerationSensitivitySource::AccelerationSensitivitySource(const mpco::ProcessInfo& info, const std::vector<int>& ids, int grad)
		: NodeResultSource(info, ids), m_grad(grad) { build(); }

	void AccelerationSensitivitySource::build()
	{
		m_schema.name = MPCO_MAKE_STRING("ACCELERATION_SENSITIVITY_" << m_grad);
		m_schema.display_name = MPCO_MAKE_STRING("Acceleration Sensitivity d_A/d_q" << m_grad);
		m_schema.num_components = 0;
		if (m_ndim == 1) {
			m_schema.components_csv = MPCO_MAKE_STRING("d_Ax/d_q" << m_grad);
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else if (m_ndim == 2) {
			m_schema.components_csv = MPCO_MAKE_STRING("d_Ax/d_q" << m_grad << ",d_Ay/d_q" << m_grad);
			m_schema.num_components = 2;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		else if (m_ndim == 3) {
			m_schema.components_csv = MPCO_MAKE_STRING("d_Ax/d_q" << m_grad << ",d_Ay/d_q" << m_grad << ",d_Az/d_q" << m_grad);
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = MPCO_MAKE_STRING("Nodal acceleration sensitivity field d_A/d_q" << m_grad);
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void AccelerationSensitivitySource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
			opserr << "MPCORecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
			return;
		}
		// NOTE: byte-faithful to frozen ResultRecorderAccelerationSensitivity —
		// the frozen bufferResponse uses getVelSensitivity() here (NOT
		// getAccSensitivity); preserved verbatim for the parity gate.
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* inode = getNode(info, m_ids[i]);
			if (!inode) continue;
			size_t j = i * m_ndim;
			buffer[j] = inode->getVelSensitivity(1, m_grad - 1);
			if (m_ndim > 1) {
				buffer[j + 1] = inode->getVelSensitivity(2, m_grad - 1);
				if (m_ndim > 2) {
					buffer[j + 2] = inode->getVelSensitivity(3, m_grad - 1);
				}
			}
		}
	}

	AngularAccelerationSensitivitySource::AngularAccelerationSensitivitySource(const mpco::ProcessInfo& info, int grad)
		: NodeResultSource(info), m_grad(grad) { build(); }
	AngularAccelerationSensitivitySource::AngularAccelerationSensitivitySource(const mpco::ProcessInfo& info, const std::vector<int>& ids, int grad)
		: NodeResultSource(info, ids), m_grad(grad) { build(); }

	void AngularAccelerationSensitivitySource::build()
	{
		m_schema.name = MPCO_MAKE_STRING("ANGULAR_ACCELERATION_SENSITIVITY_" << m_grad);
		m_schema.display_name = MPCO_MAKE_STRING("Angular Acceleration Sensitivity d_RA/d_q" << m_grad);
		m_schema.num_components = 0;
		if (m_ndim == 2) {
			m_schema.components_csv = MPCO_MAKE_STRING("d_RAz/d_q" << m_grad);
			m_schema.num_components = 1;
			m_schema.data_type = mpco::ResultDataType::Scalar;
		}
		else {
			m_schema.components_csv = MPCO_MAKE_STRING("d_RAx/d_q" << m_grad << ",d_RAy/d_q" << m_grad << ",d_RAz/d_q" << m_grad);
			m_schema.num_components = 3;
			m_schema.data_type = mpco::ResultDataType::Vectorial;
		}
		m_schema.dimension = "";
		m_schema.description = MPCO_MAKE_STRING("Nodal angular acceleration sensitivity field dRA/d_q" << m_grad);
		m_schema.result_type = mpco::ResultType::Generic;
	}

	void AngularAccelerationSensitivitySource::evaluate(const mpco::ProcessInfo& info, std::vector<double>& buffer)
	{
		buffer.assign(m_ids.size() * m_schema.num_components, 0.0);
		if (m_grad < 1 || m_grad > info.domain->getNumParameters()) {
			opserr << "MPCORecorder sensitivity parameter is out of range: grad index = " << m_grad << ", domain->getNumParameters() = " << info.domain->getNumParameters() << "\n";
			return;
		}
		// NOTE: byte-faithful to frozen ResultRecorderAngularAccelerationSensitivity
		// — uses getVelSensitivity() (NOT acceleration); preserved verbatim.
		for (size_t i = 0; i < m_ids.size(); i++) {
			Node* inode = getNode(info, m_ids[i]);
			if (!inode) continue;
			const Vector& iresponse = inode->getTrialVel(); // aux, just to check for angular accelerational components
			if (m_ndim == 2) {
				size_t j = i; // node_counter * 1 (1 angular acceleration component in 2D)
				if (iresponse.Size() > 2)
					buffer[j] = inode->getVelSensitivity(3, m_grad - 1);
				else
					buffer[j] = 0.0;
			}
			else if (m_ndim == 3) {
				size_t j = i * 3; // node_counter * 3 (3 angular acceleration components in 3D)
				if (iresponse.Size() > 5) {
					buffer[j] = inode->getVelSensitivity(4, m_grad - 1);
					buffer[j + 1] = inode->getVelSensitivity(5, m_grad - 1);
					buffer[j + 2] = inode->getVelSensitivity(6, m_grad - 1);
				}
				else {
					buffer[j] = 0.0;
					buffer[j + 1] = 0.0;
					buffer[j + 2] = 0.0;
				}
			}
		}
	}

} // namespace mpcol
