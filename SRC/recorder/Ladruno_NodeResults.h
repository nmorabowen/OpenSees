/* ********************************************************************** **
**  Ladruno recorder — modular sibling of MPCORecorder (frozen).     **
**  Ladruno_NodeResults.h — NODAL ResultSource implementations (Phase 2).   **
**                                                                        **
**  Each class here is a faithful port of one frozen `detail::node`         **
**  ResultRecorder subclass (MPCORecorder.cpp lines 1478-2594) onto the   **
**  ladruno::ResultSource contract (Ladruno_ResultIO.h):                      **
**                                                                        **
**    frozen ctor metadata        -> ResultSchema (schema())             **
**    frozen node gather          -> ids()                               **
**    frozen bufferResponse(...)   -> evaluate(info, buffer)             **
**    frozen getReactionFlag()>=0  -> requiresPartitionReduction()==true **
**                                                                        **
**  The numerical computation inside evaluate() is byte-faithful to the   **
**  frozen recorder (it feeds a 1e-12 parity gate). These sources do NOT  **
**  touch the HDF5 layer; persistence is the sink's job (Agent C).        **
** ********************************************************************** */

#ifndef Ladruno_NodeResults_h
#define Ladruno_NodeResults_h

#include "Ladruno_ResultIO.h"  // ladruno::ResultSource, ladruno::ResultSchema
#include "Ladruno_Types.h"     // ladruno::detail::ProcessInfo, enums, utils

// OpenSees (only what the node sources need; NO h5 layer here)
#include "Node.h"
#include "NodeIter.h"
#include "Domain.h"
#include "Vector.h"

#include <string>
#include <vector>

namespace ladruno {

	/*
	NodeResultSource — common base for all nodal sources.

	Mirrors the frozen `detail::node::ResultRecorder` base: it holds the schema and
	the spatial dimension (m_ndim), plus the row-identity vector of node tags.
	The node set is gathered exactly like the frozen recorder's `record()` /
	`writeModel()` paths (NodeIter over the ProcessInfo domain), either at
	construction (via the convenience ctor) or supplied explicitly.

	Subclasses fill in the schema in their ctor and override evaluate().
	*/
	class NodeResultSource : public ResultSource {
	public:
		// Gather all nodes from the domain (frozen recorder default node set).
		explicit NodeResultSource(const detail::ProcessInfo& info);
		// Or accept an explicit node-tag set (mirrors a region-restricted recorder).
		NodeResultSource(const detail::ProcessInfo& info, const std::vector<int>& node_ids);

		virtual ~NodeResultSource() = default;

		const ResultSchema& schema() const override { return m_schema; }
		const std::vector<int>& ids() const override { return m_ids; }

	protected:
		// Resolve a node tag -> Node* via the domain held in `info`.
		Node* getNode(const detail::ProcessInfo& info, int tag) const;

	protected:
		int m_ndim;                 // info.num_dimensions, frozen m_ndim
		std::vector<int> m_ids;     // node tags this source covers
		ResultSchema m_schema;      // constant for the lifetime of the source
	};

	/* ===================================================================== */
	/* Kinematics — consistent quantities, requiresPartitionReduction()==false */
	/* ===================================================================== */

	class DisplacementSource : public NodeResultSource {
	public:
		explicit DisplacementSource(const detail::ProcessInfo& info);
		DisplacementSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	class RotationSource : public NodeResultSource {
	public:
		explicit RotationSource(const detail::ProcessInfo& info);
		RotationSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	class VelocitySource : public NodeResultSource {
	public:
		explicit VelocitySource(const detail::ProcessInfo& info);
		VelocitySource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	class AngularVelocitySource : public NodeResultSource {
	public:
		explicit AngularVelocitySource(const detail::ProcessInfo& info);
		AngularVelocitySource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	class AccelerationSource : public NodeResultSource {
	public:
		explicit AccelerationSource(const detail::ProcessInfo& info);
		AccelerationSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	class AngularAccelerationSource : public NodeResultSource {
	public:
		explicit AngularAccelerationSource(const detail::ProcessInfo& info);
		AngularAccelerationSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	class PressureSource : public NodeResultSource {
	public:
		explicit PressureSource(const detail::ProcessInfo& info);
		PressureSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	/* ===================================================================== */
	/* Reactions — additive at partition-boundary nodes.                     */
	/* requiresPartitionReduction()==true (ADR decision D6): the reaction at  */
	/* a node shared between partitions is the SUM of each partition's        */
	/* contribution, so a per-step partition reduction must precede sink      */
	/* accumulation. The frozen recorder signalled this via                   */
	/* getReactionFlag()>=0; here the flag also selects which reaction the     */
	/* Node exposes (0=reaction, 1=incl. inertia, 2=Rayleigh) — the frozen     */
	/* recorder sets that flag on the Node externally before record(), so we   */
	/* expose it via reactionFlag() for the recorder driver to honour.         */
	/* ===================================================================== */

	class ReactionForceSource : public NodeResultSource {
	public:
		explicit ReactionForceSource(const detail::ProcessInfo& info);
		ReactionForceSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
		// D6: reactions are additive across partition boundaries.
		bool requiresPartitionReduction() const override { return true; }
		// Frozen getReactionFlag(): 0=reaction, 1=incl. inertia, 2=Rayleigh.
		virtual int reactionFlag() const { return 0; }
	protected:
		void build(const char* group, const char* display, const char* description);
	};

	class ReactionMomentSource : public NodeResultSource {
	public:
		explicit ReactionMomentSource(const detail::ProcessInfo& info);
		ReactionMomentSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
		bool requiresPartitionReduction() const override { return true; }
		virtual int reactionFlag() const { return 0; }
	protected:
		void build(const char* group, const char* display, const char* description);
	};

	// Inertia / Rayleigh variants only change schema name + reaction flag; the
	// numerical getReaction() path is identical (the Node returns the variant the
	// recorder driver primed via setReactionFlag-equivalent before record()).
	class ReactionForceIncInertiaSource : public ReactionForceSource {
	public:
		explicit ReactionForceIncInertiaSource(const detail::ProcessInfo& info);
		ReactionForceIncInertiaSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		int reactionFlag() const override { return 1; }
	};

	class ReactionMomentIncInertiaSource : public ReactionMomentSource {
	public:
		explicit ReactionMomentIncInertiaSource(const detail::ProcessInfo& info);
		ReactionMomentIncInertiaSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		int reactionFlag() const override { return 1; }
	};

	class RayleighForceSource : public ReactionForceSource {
	public:
		explicit RayleighForceSource(const detail::ProcessInfo& info);
		RayleighForceSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		int reactionFlag() const override { return 2; }
	};

	class RayleighMomentSource : public ReactionMomentSource {
	public:
		explicit RayleighMomentSource(const detail::ProcessInfo& info);
		RayleighMomentSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		int reactionFlag() const override { return 2; }
	};

	/* ===================================================================== */
	/* Unbalanced forces/moments — also additive at partition boundaries.    */
	/* requiresPartitionReduction()==true (ADR D6).                          */
	/* ===================================================================== */

	class UnbalancedForceSource : public NodeResultSource {
	public:
		explicit UnbalancedForceSource(const detail::ProcessInfo& info);
		UnbalancedForceSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
		// D6: unbalanced load is additive across partition boundaries.
		bool requiresPartitionReduction() const override { return true; }
	protected:
		void build(const char* group, const char* display, const char* description);
	};

	class UnbalancedMomentSource : public NodeResultSource {
	public:
		explicit UnbalancedMomentSource(const detail::ProcessInfo& info);
		UnbalancedMomentSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
		bool requiresPartitionReduction() const override { return true; }
	protected:
		void build(const char* group, const char* display, const char* description);
	};

	// IncInertia variants override evaluate() to use getUnbalancedLoadIncInertia().
	class UnbalancedForceIncInertiaSource : public UnbalancedForceSource {
	public:
		explicit UnbalancedForceIncInertiaSource(const detail::ProcessInfo& info);
		UnbalancedForceIncInertiaSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	};

	class UnbalancedMomentIncInertiaSource : public UnbalancedMomentSource {
	public:
		explicit UnbalancedMomentIncInertiaSource(const detail::ProcessInfo& info);
		UnbalancedMomentIncInertiaSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	};

	/* ===================================================================== */
	/* Modes of vibration — SPECIAL CASE (multi-mode per step).              */
	/*                                                                       */
	/* In the frozen recorder ModesOfVibration overrides record() to loop    */
	/* over OPS_GetNumEigen() modes, writing one MODE_<k> dataset per mode    */
	/* under a per-step group, each carrying MODE/LAMBDA/OMEGA/FREQUENCY/     */
	/* PERIOD attributes. The ResultSource contract yields ONE flat buffer    */
	/* per evaluate(), so we expose:                                          */
	/*   - evaluate(info, buf): the current mode (set via setCurrentMode())   */
	/*   - numModes(info), setCurrentMode(k), modeInfo(info, k, ...)          */
	/* The StreamingSink (Agent C) MUST drive the mode loop: for k in         */
	/* [0,numModes): setCurrentMode(k); evaluate(...); write MODE_<k> + attrs.*/
	/* requiresPartitionReduction()==false (eigenvectors are consistent).     */
	/* Recording is also gated on info.record_eigen_on_this_step by the sink. */
	/* ===================================================================== */

	class ModesOfVibrationSource : public NodeResultSource {
	public:
		explicit ModesOfVibrationSource(const detail::ProcessInfo& info);
		ModesOfVibrationSource(const detail::ProcessInfo& info, const std::vector<int>& ids);

		// Fill buffer for the CURRENT mode (translational components).
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;

		// Mode-aware API for the mode-loop driver (StreamingSink).
		int numModes(const detail::ProcessInfo& info) const;  // *OPS_GetNumEigen()
		void setCurrentMode(int k) { m_current_mode = k; }
		int currentMode() const { return m_current_mode; }
		// Per-mode scalars (frozen ModesOfVibration::record math).
		void modeInfo(const detail::ProcessInfo& info, int k,
		              double& lambda, double& omega,
		              double& freq, double& period) const;

	protected:
		void build();
		int m_current_mode;  // frozen m_current_mode_to_buffer
	};

	class ModesOfVibrationRotationalSource : public ModesOfVibrationSource {
	public:
		explicit ModesOfVibrationRotationalSource(const detail::ProcessInfo& info);
		ModesOfVibrationRotationalSource(const detail::ProcessInfo& info, const std::vector<int>& ids);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
	};

	/* ===================================================================== */
	/* Sensitivity sources — consistent, requiresPartitionReduction()==false. */
	/* Each carries a gradient index (m_grad); evaluate() guards m_grad        */
	/* against domain->getNumParameters() exactly as the frozen recorder.      */
	/* ===================================================================== */

	class DisplacementSensitivitySource : public NodeResultSource {
	public:
		DisplacementSensitivitySource(const detail::ProcessInfo& info, int grad);
		DisplacementSensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
		int m_grad;
	};

	class RotationSensitivitySource : public NodeResultSource {
	public:
		RotationSensitivitySource(const detail::ProcessInfo& info, int grad);
		RotationSensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
		int m_grad;
	};

	class VelocitySensitivitySource : public NodeResultSource {
	public:
		VelocitySensitivitySource(const detail::ProcessInfo& info, int grad);
		VelocitySensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
		int m_grad;
	};

	class AngularVelocitySensitivitySource : public NodeResultSource {
	public:
		AngularVelocitySensitivitySource(const detail::ProcessInfo& info, int grad);
		AngularVelocitySensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
		int m_grad;
	};

	class AccelerationSensitivitySource : public NodeResultSource {
	public:
		AccelerationSensitivitySource(const detail::ProcessInfo& info, int grad);
		AccelerationSensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
		int m_grad;
	};

	class AngularAccelerationSensitivitySource : public NodeResultSource {
	public:
		AngularAccelerationSensitivitySource(const detail::ProcessInfo& info, int grad);
		AngularAccelerationSensitivitySource(const detail::ProcessInfo& info, const std::vector<int>& ids, int grad);
		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
	private:
		void build();
		int m_grad;
	};

} // namespace ladruno

#endif // Ladruno_NodeResults_h
