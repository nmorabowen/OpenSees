---
wp: LEGACY
title: "LadrunoBrick::update() DISCARDED setTrialStrain's return code on four of its five formulation paths — so no material could ever fail a step through it"
legacy_seq: 373
---
## `LadrunoBrick::update()` DISCARDED `setTrialStrain`'s return code on four of its five formulation paths — so no material could ever fail a step through it

**Found 2026-09-05, WP-86b (ADR-86b), while wiring the SANISAND substep cap. Fixed in the same PR.**

The `stdBrick`-swallows-return-codes trap is already on this ledger. What was not known is that the
fork's own brick had the same hole on most of its paths, and that the two paths which DID propagate
are the newest ones — so it reads, from the code, as if the element already had the contract.

`LadrunoBrick::update()` dispatches by formulation. Before WP-86b:

| path | propagated a material failure? |
|---|---|
| `-geom finite` -> `updateFinite()` | yes |
| `-geom hypo` -> `updateHypo()` (`:1720`, `< 0`) | yes |
| `Formulation::EAS` -> `formEAStrue()` | yes |
| `Formulation::SSP` (centroid slot 0) | **no** — bare call, `return 0` |
| `Formulation::URI` + `Hourglass::PHYSICAL` (8 GPs) | **no** |
| `Formulation::URI` (perturbation, centroid) | **no** |
| std / **b-bar** (the 8-GP default loop) | **no** |

So on `-formulation bbar` — the formulation ADR-90 S3 freezes for the whole SANISAND campaign — a
material that refused the increment was indistinguishable from one that integrated it. All four now
test for a refusal, name the element and Gauss point through one throttled reporter, and return -1.

- **Why it stayed invisible:** almost no OpenSees `NDMaterial` returns non-zero from
  `setTrialStrain` — the whole UW family hardcodes `return 0`, and so did both `LadrunoSANISAND`
  wrappers until ADR-86b. A contract nothing exercises is a contract nobody notices is missing.
- **The tell, if you are looking for it:** an element `update()` that ends in a bare
  `materialPointers[i]->setTrialStrain(strain);` with no `if`. Grep for that shape before assuming
  any element propagates material failure.
- **Upstream `stdBrick` still does not propagate.** Any gate on a material's failure return must use
  `LadrunoBrick` (or another element you have checked), or it silently tests nothing.
