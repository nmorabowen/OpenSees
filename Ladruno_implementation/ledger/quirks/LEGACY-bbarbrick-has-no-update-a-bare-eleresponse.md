---
wp: LEGACY
title: "BbarBrick has no update() — a bare eleResponse(\"stresses\") reads the *predictor* (u=0) state"
legacy_seq: 18
---
### `BbarBrick` has no `update()` — a bare `eleResponse("stresses")` reads the *predictor* (u=0) state
- **Bites:** after a static/linear solve, `ops.eleResponse(tag, "stresses")` on an
  upstream `bbarBrick` returns **all zeros** even though the displacements are
  correct. A regression that compares element stresses against `bbarBrick` then
  fails with `ours=<real> vs ref=0`. (`stdBrick` does **not** show this — it reads
  back real stresses.)
- **Why:** `Brick` overrides `Element::update()` to push the material trial strain
  every step, so its committed `getStress()` reflects the solved displacement.
  `BbarBrick` has **no `update()`** — it calls `setTrialStrain` lazily, only inside
  `formResidAndTangent`. The `"stresses"` response (responseID 3) just returns
  `materialPointers[i]->getStress()` *without* recomputing, so it reflects whatever
  the last `formResidAndTangent` saw — which, after a linear step, is the predictor
  state (u=0 → zero strain → zero stress).
- **Workaround:** read `"forces"` (responseID 1 → `getResistingForce` →
  `formResidAndTangent` → `setTrialStrain` at the committed u) **before** reading
  `"stresses"`. Then both lazy- and eager-strain elements report the solved state.
  `LadrunoBrick` implements `update()` (like `Brick`), so it is order-insensitive —
  this only matters when comparing against `bbarBrick`.
- **How it surfaced:** `tests/test_ladrunoBrick_element.py::test_bbar_matches_bbarBrick`
  (LadrunoBrick `bbar`↔`bbarBrick` regression, PR [#65](https://github.com/nmorabowen/OpenSees/pull/65)).
  Displacements matched to 1e-9 (kernel-equivalent); only the stress readback timing
  differed. Date learned: 2026-06-01.
