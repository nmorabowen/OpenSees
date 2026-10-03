---
wp: LEGACY
title: "system(\"FullGeneral\") hard-crashes the whole process on a model with zero free equations"
legacy_seq: 297
---
### `system("FullGeneral")` hard-crashes the whole process on a model with zero free equations
- **Bites:** any model where every DOF ends up fixed or sp-prescribed — e.g. a fully strain-driven single-brick material-point driver under `constraints("Transformation")`, the standard ASDPlasticMaterial3D unit-test rig — dies the instant `system("FullGeneral")` is selected. It is material-independent: `ElasticIsotropic`, plain `MohrCoulomb_YF`, and the new `MohrCoulombTensionCutoff_YF` all die identically at analysis step 1. There is no Python traceback, `faulthandler` prints nothing, and the process exits 255/-1 — it looks exactly like a fresh material bug in whatever is under test, and cost about an hour of bisection before the culprit turned out to be the solver, not the material. `UmfPack` returns `rc=0` on the byte-identical model.
- **Why:** the `FullGeneral` SOE/solver path does not guard `N==0` free equations; a fully-prescribed system legitimately has none.
- **Workaround/status (2026-08-12, found during PR #741 / ADR-84 P0 verification):** use `UmfPack` (or another solver with an N=0 guard) for fully-prescribed material-point drivers. Root fix (a `FullGenLinSOE` N=0 guard) is filed as its own task, not part of this PR — pre-existing core defect, out of scope for the MCTC feature.
