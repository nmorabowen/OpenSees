---
wp: WP-200
title: "Force-based element energy test |dW| < tol passes on a flat section plateau whatever the section unbalance"
date: 2026-10-29
---
### Force-based column "converges" with tip axial displacement ~1e12 mm, sections out of equilibrium with the element force

- **Bites:** a crushed or fully yielded force-based column (every fibre on a
  plateau with tangent ~1e-10) commits a state with the element force equal to
  the applied P but section axial forces off by 1-400 % and section deformations
  unrelated to the nodal ones. frame-ladruno-validation recipe 02 (2026-10-29),
  `LadrunoForceBeamColumn -iter 10 1e-4`: Gill S1 increment 344 (tip u_z 8.7e11
  mm, every IP at ε = -2.6e4, section N = -1837 kN against q0 = -1815 kN),
  Takizawa & Aoyama 3D step 467 (u_z 1.3e12 mm, section N from -329 to +685 kN
  against q0 = -157 kN). The next increment fails on every rung.
- **Why:** `evalBulkState` (shared with upstream `ForceBeamColumn`) accepts on
  `|dW| < tol`, `dW = dv . kv . dv`, `dv = v - vr`, where `vr` integrates
  `vs + fs (Ss - Ssr)`. On a plateau `fs` is huge, so `kv = f^-1 -> 0` (2.5e-8
  N/mm on Gill) and the linearized correction `fs (Ss - Ssr)` absorbs ANY target
  `v`. The energy test is then blind to the unbalance itself. The global
  `NormDispIncr` test accepts too: with an axial stiffness of 2.5e-8 N/mm the
  corrections are tiny once the element force equals P.
- **Measure, not dW:** the compatibility error of the section deformations
  actually set, `e = v - (vr - integral(b^T fs (Ss - Ssr)))`. It is round-off on
  the tangent branch (`l = 0`), up to ~1e-2 relative on the initial-flexibility
  branches (`l = 1, 2`, measured over the beam/hinge suites and the three
  specimens), and about 1 in every garbage state.
- **Fix/status:** fixed 2026-10-29 in `LadrunoForceBeamColumn2d/3d` (vanilla
  `ForceBeamColumn` keeps the blind test). A bulk state is accepted only if also
  `max|e| <= 1e-2 x deformation scale` (axial component / L); otherwise the
  element keeps iterating and, after the subdivisions, fails with
  `failed to get compatible element forces & deformations for element: <tag>
  (dW, tol, compatibility error)`. Histories that never hit a garbage state are
  unchanged: Gill to 1.5e-9 of the peak force, LM1 and Takizawa bit-identical
  up to the old garbage step. Regression:
  `tests/test_ladruno_force_beam_column_garbage_state.py`
  (`test_lfbc2d/3d_plateau_jump_is_solution_or_failure`).
