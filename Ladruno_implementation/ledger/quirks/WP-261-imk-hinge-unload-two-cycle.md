---
wp: WP-261
title: "LadrunoIMKBeam local hinge Newton 2-cycles on unloading and committed the unconverged iterate silently"
date: 2026-11-26
status: "fixed-by WP-261"
---
### `LadrunoIMKBeam` frame needs KrylovNewton / line-search recovery at hinge unloading, or a hinge commits a moment off its own material law

- **Bites:** steel moment frames with bilinear `Steel01` hinges (b = 1e-4, k0 = 100 × 6EI/L,
  so k0 ≈ 150 a with a = 4EI/L). The global Newton cycles to its iteration limit at steps where
  hinges unload, so the driver falls back to KrylovNewton or NewtonLineSearch (20-story frame:
  420 recovered steps in 21 records). Worse, the state can be wrong without any failure. In
  20-story NST-N R = 2, step 5644 (recovered by line search) committed two roof-beam hinges with
  M ≈ 2 My away from the Steel01 clip of their own rotation, and `NormDispIncr` passed because the
  line-search increment was small. A single element under rotation control shows it too
  (`tests/test_ladrunoIMKBeam_unload.py`). Unloading by 0.002 rad from 0.02 rad "converged" in one
  global iteration onto θj = 0.0177, Mj = 1.4e3. The clip of that rotation is 1.0e5, and the
  exact answer is θj = 1.19e-3, Mj = 8.95e4.
- **Why:** the element-internal Newton of `ladrunoIMKSolveAxis` / `ladrunoIMKSolveBiaxial`
  (`LadrunoIMKHinge.h`) solves g(θ) = a(v − θ) − M(θ) = 0 per hinge. With a bilinear law g is
  monotone but S-shaped: slope ≈ a on the plastic branches, a + k0 on the narrow elastic band.
  From a plastic committed state, a Newton step with the flat plastic slope overshoots the band
  by ~150× onto the opposite plastic branch, and the flat slope there sends it back: the arctan
  2-cycle. Loading from the elastic side converges in two iterations, so only unloading bites.
  Then, after 50 iterations, the loop exited and returned anyway. Its last update had moved θ
  AFTER the materials were evaluated, so `q` came from the last θ while the material's trial
  stress, tangent and (at commit) state came from the previous iterate. `update()` returned 0.
  A transcription of the local Newton reproduces 266 of the 268 recovered steps of a 69-run suite.
  Given the exact element state, the global Newton never needs more than 4 iterations.
- **Fix (WP-261):** the plain loop still runs first, unchanged. A state it solves is
  bit-identical. When it exits unconverged, `ladrunoIMKSafeguardedSolve`
  (`LadrunoIMKSafeSolve.h`) restarts from the warm start. It does Newton with a bracketed line
  search on h(t) = G(x + t d)·d, which is −dΠ/dt of the convex series potential while
  K_el + K_H is positive definite (any non-softening law), so it is nonincreasing in t. When the
  full step passes the line minimum, the root of h in (0, 1) is found by Newton-in-t /
  Illinois / bisection. If the tangent is indefinite, it falls back to steepest descent. The
  fallback uses only material trial calls, so it is law-agnostic (UniaxialMaterial and the biaxial
  NDMaterial path). Its last material call is at the returned point. If it still fails,
  `update()` returns −1 with a capped warning, so `Domain::update` → integrator → algorithm
  rejects the step (`ladrunoKeepWorstCode`; every Newton and line-search `update()` check is
  `< 0`). The materials, moments and q stay consistent with the returned θ.
- **Tolerance trap found on the way:** the plain test `|g| ≤ 1e-10 (|q| + |M|) + 1e-12` is below
  round-off whenever a hinge holds a plastic rotation at M ≈ 0. q = a(v − θ) is a difference
  of large terms. The old loop silently accepted those states. Failing `update()` there
  would reject good steps. The fallback adds a round-off floor `1e-12 × Σ|terms|` and
  reports failure only above `1e-8 × (|q| + |M| + Σ|terms|)`.
- **Second trap: Steel01's absolute dead zone.** `Steel01::setTrialStrain` ignores
  |Δε| ≤ `DBL_EPSILON` (absolute, not relative) and returns the committed stress. So near a
  committed point, M(θ) jumps by k0·2.2e-16. For ε < v < (1 + k0/a)·ε, the first steps of
  every ground motion, the hinge equation has **no root**. The first WP-261 build reported that
  as a failure and rejected step 0 of 13 of 26 frame9 runs. The fallback now also accepts
  residuals below `k_ref × 1e-14 rad` (tight) / `k_ref × 1e-12 rad` (reporting), with k_ref =
  a + b + the hinges' initial stiffnesses. A real 2-cycle misses by ~My, ten orders above. Any
  law with an absolute strain threshold behaves the same way.
- **Gates:** `tests/test_ladrunoIMK_safe_solve_cpp.py` (g++ harness: the plain loop fails the
  unloading case and 555 of 4000 sweep cases. The fallback solves all of them to round-off against
  Steel01 region enumeration, and returns false on a no-root law with its last evaluation at x).
  `tests/test_ladrunoIMKBeam_unload.py` (2D and 3D element vs an exact series oracle through a
  cyclic protocol).
