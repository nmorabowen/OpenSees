---
wp: LEGACY
title: "updateMaterialStage CANNOT reach a LoadPattern (overlay moduli must be re-set via the parameter ... loadPattern route after a stage flip)"
legacy_seq: 187
---
### `updateMaterialStage` CANNOT reach a `LoadPattern` (overlay moduli must be re-set via the `parameter ... loadPattern` route after a stage flip)
- **Bites:** you flip a soil constitutive stage mid-analysis with `updateMaterialStage $matTag $stage` (PDMY/PM4Sand elastic→plastic), expecting the `LadrunoPorousOverlay`'s fixed-stress `L` factor (which uses the drained skeleton moduli) to follow the new stage. It does not — the overlay keeps its stage-0 moduli, so `L` is stale (stable, but convergence-degraded — the fixed-stress split's L is a preconditioner-like term, not a physics error).
- **Why:** `MaterialStageParameter` registers only the FIRST accepting ELEMENT in its domain scan and never scans load patterns (the ADR-71 sibling-broadcast trap, family-documented). A `LoadPattern` subclass like the overlay is simply not on the path `updateMaterialStage` walks.
- **Workaround/status (2026-07-17, ADR-73 P2):** the transport contract is explicit — after a stage flip the USER re-sets overlay moduli through the EXISTING parameter route: `parameter $p loadPattern $overlayTag E $newE` (or `nu`, `layerE $i`, `layerNu $i`), which marks the overlay `moduliDirty_` and lazily rebuilds `aS_`/`aL_` at the next fluid use. A flip without a re-set keeps stage-0 `L`. The PDMY staged-liquefaction battery + the P4 guide pin the recipe.
