/* ********************************************************************** **
**  Ladruno recorder — modular sibling of MPCORecorder (frozen).     **
**  Ladruno_OverlayResults.h — OVERLAY pore-pressure ResultSource (ADR 73). **
**                                                                        **
**  OverlayPressureSource (ADR-73 P4): the committed pore-pressure field   **
**  pⁿ⁺¹ of one LadrunoPorousOverlay (PATTERN_TAG 33022), exposed as a      **
**  per-region-node nodal scalar. ids() == the overlay's region-node order **
**  (getRegionNodeTags — THE canonical id order); evaluate() copies        **
**  getCommittedP() (aligned to that order). One component "p". The field  **
**  is per-process serial in v1, so no partition reduction. The source     **
**  touches no HDF5; persistence is the generic StreamingSink/EnvelopeSink **
**  (OnNodes), so this drops in with ZERO sink changes.                    **
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

#ifndef Ladruno_OverlayResults_h
#define Ladruno_OverlayResults_h

#include "Ladruno_ResultIO.h"      // ladruno::ResultSource, ladruno::ResultSchema
#include "Ladruno_Types.h"         // ladruno::detail::ProcessInfo

#include <string>
#include <vector>

class LadrunoPorousOverlay;

namespace ladruno {

	/*
	OverlayPressureSource — one overlay's committed pore-pressure field.

	Constructed from a resolved (snapshot-ready) LadrunoPorousOverlay*. The region
	node tags are COPIED into m_ids at construction so schema()/ids() stay constant
	for the source lifetime (the recorder rebuilds every source on a stage change,
	so a copy never goes stale within a stage). evaluate() copies the overlay's
	committed p (freshly committed at Domain::commit BEFORE recorders run — the
	4.3 ordering pin), row-major over ids() x 1 component.

	Removal robustness (P4 panel F-1): `remove loadPattern` DELETES the pattern
	without firing a domain change (no owned SPs), so a cached pointer would
	dangle for the rest of the stage. evaluate() therefore re-resolves the
	overlay BY TAG through info.domain each step (the Monitor idiom) and
	zero-fills when it is gone — never a stale deref.
	*/
	class OverlayPressureSource : public ResultSource {
	public:
		// overlay MUST be non-null and snapshotReady() (the recorder guards this).
		OverlayPressureSource(const detail::ProcessInfo& info,
		                      LadrunoPorousOverlay* overlay);
		~OverlayPressureSource() override = default;

		const ResultSchema& schema() const override { return m_schema; }
		const std::vector<int>& ids() const override { return m_ids; }

		void evaluate(const detail::ProcessInfo& info, std::vector<double>& buffer) override;

		// Overlay is per-process serial (v1); its field never needs a cross-
		// partition reduction before a sink accumulates it.
		bool requiresPartitionReduction() const override { return false; }

	private:
		void buildSchema(int patternTag);

	private:
		int m_patternTag;                  // re-resolved by tag in evaluate() (F-1)
		std::vector<int> m_ids;            // region-node tags (copy, construction-time)
		ResultSchema m_schema;
	};

} // namespace ladruno

#endif // Ladruno_OverlayResults_h
