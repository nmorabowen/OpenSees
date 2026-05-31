/* ********************************************************************** **
**  Ladruno recorder — modular sibling of MPCORecorder (frozen).     **
**  Ladruno_DomainResults.h — DOMAIN / REGION ResultSource implementations.  **
**                                                                        **
**  EnergyBalanceSource (ADR D8): whole-model OR per-region structural-    **
**  dynamics energy balance (KE/IE/DW/ULW/RES/ERR). The energy math is     **
**  shared with EnergyBalanceRecorder through ebkernel (EnergyBalanceKernel**
**  .h) — ONE definition. The whole-model source feeds RESULTS/ON_DOMAIN/  **
**  energyBalance; the per-region source feeds RESULTS/ON_REGIONS/         **
**  energyBalance (one row per region tag).                                **
** ********************************************************************** */

#ifndef Ladruno_DomainResults_h
#define Ladruno_DomainResults_h

#include "Ladruno_ResultIO.h"      // ladruno::ResultSource, ladruno::ResultSchema
#include "Ladruno_Types.h"         // ladruno::detail::ProcessInfo, enums

#include "EnergyBalanceKernel.h" // ebkernel::EnergyAccumulator + sweep helpers

#include "Vector.h"

#include <string>
#include <vector>

namespace ladruno {

	/*
	EnergyBalanceSource — whole-model OR per-region energy balance.

	Two construction modes:
	  - Whole-model: ids() == {0}; evaluate() runs one accumulator over the
	    whole domain -> 6 doubles.
	  - Per-region : ids() == [regionTags]; evaluate() sweeps each region (in
	    the given tag order), stepping its own accumulator -> 6 doubles per
	    region, row-major [nRegions x 6].

	Unlike the stateless node/element sources, this source carries cross-step
	state (the trapezoidal work integrals + closure live in the accumulators,
	plus the previous time stamp). evaluate() therefore advances that state on
	every call — it must be invoked exactly once per recorded step.
	*/
	class EnergyBalanceSource : public ResultSource {
	public:
		// Whole-model (ON_DOMAIN). ids() == {0}.
		explicit EnergyBalanceSource(const detail::ProcessInfo& info);
		// Per-region (ON_REGIONS). ids() == region_tags.
		EnergyBalanceSource(const detail::ProcessInfo& info,
		                    const std::vector<int>& region_tags);

		~EnergyBalanceSource() override = default;

		const ResultSchema& schema() const override { return m_schema; }
		const std::vector<int>& ids() const override { return m_ids; }

		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;

		// ADR D6/D8: energy is an additive/global quantity that needs a per-step
		// partition reduction (Allreduce) before a sink may accumulate it. v2
		// ships serial-correct; the actual MPI reduce is v3b (NOT implemented
		// here — this only sets the routing flag).
		bool requiresPartitionReduction() const override { return true; }

	private:
		void buildSchema();

	private:
		bool m_per_region;                 // false => whole-model {0}
		std::vector<int> m_ids;            // {0} or the region tags
		ResultSchema m_schema;

		// cross-step accumulator state (ebkernel = single source of truth)
		ebkernel::EnergyAccumulator m_model_acc;              // whole-model
		std::vector<ebkernel::EnergyAccumulator> m_region_acc; // per-region

		double m_prev_time;                // previous recorded time (for dT)
		bool m_first;                      // first evaluate() only seeds rates
		Vector m_velScratch;               // caller-owned element-DOF scratch
	};

} // namespace ladruno

#endif // Ladruno_DomainResults_h
