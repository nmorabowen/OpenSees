---
wp: WP-219
title: "LadrunoDispBeamColumn -spanLoad bubble diverged on cracked sections and reported success"
pr: "#1056"
date: 2026-11-01
---
### `LadrunoDispBeamColumn -spanLoad` push stops on a fine mesh (Newton cycles at a large norm), `-noSpanLoad` completes, force changes after a converged `analyze`

- **Bites:** `-spanLoad` with fiber RC sections on 4+ elements per member. frame-ladruno-validation
  04a: 80/141 and 21/141 monotonic steps at 4 and 8 elements. 1 and 2 elements happened to be fine.
- **Why:** a cracked section with one bar layer yielded leaves the bubble stiffness nearly rank one,
  and softening fibers make it indefinite (scaled `sPP*sZZ - sPZ^2` about -0.004). The inner Newton
  diverged (`rp` 1e6 N against a test of 1e-3), hit its cap and returned success, 3356 times in the
  4-element run, storing a beta that did not match the section state. The force was then not a
  function of the deformation, which also explains the "same residual under both tests, force moves
  after analyze": `ops.reactions()` re-applies loads, and the re-solve warm-starts from the new
  committed beta. A line search on `|r|` is not enough: it stalls at the kinks of `|r|`.
- **Fix (WP-219):** with the committed history fixed each fiber's stress is a function of its strain,
  so the residual is the gradient of a potential in beta. Use a positive definite Newton direction and
  a line search on the directional derivative of that potential. On exhausting the budget, stand on the
  best iterate, accept it within 1e-8 of the force scale, otherwise return -2. Unchanged where the old
  solve converged: 1 and 2 elements, same peak to 6 digits.
- **Result:** 4 and 8 elements complete monotonic (141/141) and cyclic (276/276). 8 elements: peak
  483.96 kN against 483.90 with -noSpanLoad and vanilla; cyclic 471.71 against 471.57.
- **Open:** under `constraints Transformation`, `enforceSPs` still drops the -2 (vanilla handler, see
  the WP-219 recommendation); LadrunoFBC has its own `trialFailed` guard, this element does not.
- **Gate:** `tests/test_ladrunoDispBeamColumn_spanload_softening.py` (2D and 3D).
