---
wp: ADR-78
title: "ADR-78 unified the mortar-tie guards for ALL facet orders — linear decks with cancelling gap fields now refuse (by design)"
legacy_seq: 265
---
### ADR-78 unified the mortar-tie guards for ALL facet orders — linear decks with cancelling gap fields now refuse (by design)
- **Bites:** the pre-ADR-78 conforming-gap guard tested the per-node SIGNED weighted gap `|Σ∫N_I g_N|/cover`, which a gap field that cancels inside a node's support could slip through (a warp/antisymmetric offset reading as "conforming"). The ADR-78 per-facet `∫|g_N|/area` L1 guard has no cancellation blind spot and — per the OQ-2 "unify" sign-off — applies to tri3/quad4 decks too. A linear deck that previously built its tie may now refuse with the conforming-gap message. The emitted P (weights) is byte-identical for linear inputs; only refusal behaviour changed.
- **Rule:** if a formerly-working linear tie now refuses on the gap guard, the geometry genuinely is off the master surface somewhere — fix the interface or consciously relax `-tol`.
- **Workaround/status:** ✅ intentional behaviour change, documented here + ADR-78 D3/BLOCKER-3; regression-tested (`test_refuse_accordion_gap_L1`). Two adversarial-gate amendments same day: (a) the threshold scale is `0.5·sqrt(areaFull)` — a bare `sqrt(A)` was 2× LOOSER than the shipped per-node tributary scale (review MINOR); (b) the per-facet area-coverage sum counts MULTIPLICITY, so a self-overlapping / doubly-listed MASTER surface could exactly mask an uncovered slave strip — now a named refusal (coincident master-master overlap, mean-gap-gated so curved masters stay legal; review MAJOR, `test_refuse_duplicate_master_facets`). A linear deck with duplicated master facets that previously "worked" now refuses. *2026-08-04.*
