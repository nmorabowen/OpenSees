---
wp: WP-152
title: "InitialStateAnalysis off ends with a domain update on the ZEROED displacements: a SANISAND point's trial jumps by −ε_n (pre-existing, found in WP-152)"
legacy_seq: 538
---
### `InitialStateAnalysis off` ends with a domain update on the ZEROED displacements: a SANISAND point's trial jumps by −ε_n (pre-existing, found in WP-152)
- **Bites:** reading a material's trial state right after `InitialStateAnalysis off` (or taking a step from it) with a strain-driven material that keeps its committed strain under ISA (ManzariDafalias / LadrunoSANISAND).
  - `OPS_InitialStateAnalysis` "off" calls `Domain::revertToStart`, which zeroes the displacements and ENDS with `this->update()`; the material's `revertToStart` keeps σ and ε_n under ISA, so the update integrates an increment of −ε_n.
  - Measured: a NORMAL point at p ≈ 1.77 kPa: trial p 1.769 → 1.731 kPa after ISA off; a separated point reads the same jump as closing (trial re-contact at 2.66 kPa).
- **Rule:** Under ISA, trust the COMMITTED state (and `sepActive`), not the trial read right after "off"; a deck that continues from ISA with SANISAND inherits the strain-frame jump.
- **Workaround/status:** open, not WP-152's (it predates it). WP-152 only makes the separation state survive the ISA revert (review #5).
