---
wp: LEGACY
title: "MPIDiagonalSOE::assembleSharedSum reuses structures built by the FIRST solve() — calling it pre-solve SILENTLY drops the cross-rank sum (now tripwired)."
legacy_seq: 114
---
### `MPIDiagonalSOE::assembleSharedSum` reuses structures built by the FIRST `solve()` — calling it pre-solve SILENTLY drops the cross-rank sum (now tripwired).
- The neighbour exchange (`myActualNeighborsBsToSend`/`myNeighborsSizes`/posloc) and the factored GLOBAL `getScalingDiagonalA()` (= `1/mass` summed across ranks) are produced by `MPIDiagonalSolver::solve()`'s first (`notSet`) pass. The distributed consistent PCG's correctness therefore depends on the implicit invariant **first `solve()` precedes first `refineAccel()` on every rank** — true for all 3 shipped consistent integrators (CDL starter + ExplicitBathe/LNVD both sub-steps all do `solve()`→`refineAccel`). A new explicit integrator using the consistent path MUST preserve that order. A one-time tripwire warning now fires if `assembleSharedSum` runs with neighbours but un-built buffers (was a silent no-op → wrong answer). This was the unanimous residual flag from the 4-lens adversarial review (which found ZERO actual bugs). Learned 2026-06-21, [[38_ladruno_consistent_mass_scaling_adr]] V5.
