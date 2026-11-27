---
wp: WP-261
title: "LadrunoIMKBeam safeguarded local hinge solve (unloading 2-cycle fix)"
pr: "#1168"
date: 2026-11-26
status: "fixed on the branch"
section: "table"
---
| **IMK beam safeguarded local hinge solve (WP-261)**: the local hinge Newton of `ladrunoIMKSolveAxis` / `ladrunoIMKSolveBiaxial` 2-cycled when a hinge unloaded from a plastic state, and on exit it returned q from the last θ with the materials at the previous iterate. A line-search recovery could then commit a hinge off its own law. Now the plain loop runs first, unchanged (bit-identical wherever it converged). If it fails, `ladrunoIMKSafeguardedSolve` restarts it: Newton with a bracketed line search on the directional residual (descent on the convex series potential), law-agnostic. If that also fails, `update()` returns −1 so the step is rejected. **Behaviour change** only on states the old loop did not converge. Quirk: "`LadrunoIMKBeam` frame needs KrylovNewton / line-search recovery at hinge unloading". | bug fix | — (existing IMK tags) | `SRC/element/ladrunoIMKBeam/LadrunoIMKSafeSolve.h` (new), `LadrunoIMKHinge.h`, `LadrunoIMKBeam{,2d}.cpp`, `tests/test_ladrunoIMK_safe_solve_cpp.py`, `tests/ladrunoIMK_safe_solve_check.cpp`, `tests/test_ladrunoIMKBeam_unload.py` | **fixed on the branch** | #1168 |
