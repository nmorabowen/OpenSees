---
wp: PR-877
title: "A failed implicit return REFUSES under -implex and under explicit integrators (WP concrete3d-hang-diagnosis, #877)"
pr: "#877"
date: 2026-10-10
---
### A failed implicit return REFUSES under `-implex` and under explicit integrators (WP concrete3d-hang-diagnosis, [#877](https://github.com/nmorabowen/OpenSees/pull/877))
- **Bites:** the local return map used to fall back silently to the elastic trial and commit it (`Kernel.h` returnMap tier-2 status), so an IMPL-EX or explicit run "converged" on garbage (one committed wild step left hundreds of MPa of nominal stress). `setTrialStrain` now returns `LADRUNO_MATERIAL_REFUSED` and `commitState` declares the refusal to `Domain::commit()` (WP-99 counter). Hosts that forward the code cut the step at the trial (`LadrunoBrick`, `LadrunoQuad`, `LadrunoSolidShell`, `LogStrainNDMaterial`); a host that discards it (vanilla `stdBrick`) is only ever stopped at COMMIT, after the nodes and sibling integration points have committed - the model is then inconsistent and must be restarted (apeGmsh `AnalysisAbortedError`). `-implexControl` is NOT an option of this material (the parser rejects it): under `-implex` the implicit companion of the IMPL-EX pair runs at TRIAL time inside `returnMap`, so a trial refuses whenever that implicit solve fails, and the recovery is the global step cut. A commit-time implicit companion does not exist (register 32's trial-vs-commit rule would need the implicit solve moved into `commitState`).
