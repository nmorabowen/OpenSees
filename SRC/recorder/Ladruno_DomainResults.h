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
	plus the previous time stamp).

	WP-165 (R3): that state lives in an EnergyState OWNED BY THE RECORDER, so it
	survives the source being rebuilt at a MODEL_STAGE change (it used to restart
	at zero and take a rate x t_total jump on the first step of every stage), and
	advance() — idempotent per commit tag — is called by the recorder on EVERY
	commit, so the work integrals are not sampled only at the -T recorded steps.
	evaluate() = advance() + copy the latest values.
	*/
	struct EnergyState {
		ebkernel::EnergyAccumulator model;                // whole-model
		std::vector<ebkernel::EnergyAccumulator> regions; // per-region, tag order
		double prev_time = 0.0;      // time of the last integrated commit
		bool first = true;           // the first integrated commit seeds the rates
		bool has_last = false;
		int last_tag = 0;            // commit tag last integrated (idempotency key)
		std::vector<double> last_out;   // [nIds x 6] of the last integrated commit
		Vector velScratch;           // element-DOF scratch
	};

	class EnergyBalanceSource : public ResultSource {
	public:
		// Whole-model (ON_DOMAIN). ids() == {0}.
		// `state` (optional) is recorder-owned and outlives the source (WP-165 R3);
		// without it the source owns a private one (stage-local, the old behaviour).
		explicit EnergyBalanceSource(const detail::ProcessInfo& info, EnergyState* state = 0);
		// Per-region (ON_REGIONS). ids() == region_tags.
		EnergyBalanceSource(const detail::ProcessInfo& info,
		                    const std::vector<int>& region_tags, EnergyState* state = 0);

		~EnergyBalanceSource() override = default;

		const ResultSchema& schema() const override { return m_schema; }
		const std::vector<int>& ids() const override { return m_ids; }

		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;
		// WP-165 (R3): integrate this commit (no-op if this commit tag was already
		// integrated). The recorder calls it on every commit.
		void advance(const detail::ProcessInfo& info);

		// ADR D6/D8: energy is an additive/global quantity that needs a per-step
		// partition reduction (Allreduce) before a sink may accumulate it. v2
		// ships serial-correct; the actual MPI reduce is v3b (NOT implemented
		// here — this only sets the routing flag).
		bool requiresPartitionReduction() const override { return true; }
		// WP-126: NOT "SUM" -- KE at a node shared by two partitions is counted in both,
		// and RES/ERR are derived from the other components, so a componentwise sum of
		// the partition files is wrong. Only the v3b Allreduce + recompute can merge it.
		const char* partitionReduction() const override { return "UNSUPPORTED"; }

	private:
		void buildSchema();

	private:
		bool m_per_region;                 // false => whole-model {0}
		std::vector<int> m_ids;            // {0} or the region tags
		ResultSchema m_schema;

		// cross-step state (WP-165 R3: recorder-owned when given)
		EnergyState m_own;                 // used only when no state was passed
		EnergyState* m_state;              // -> recorder's state, or &m_own
	};

} // namespace ladruno

#endif // Ladruno_DomainResults_h
