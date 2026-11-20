---
wp: WP-144
title: "Kernel-vs-oracle TANGENT parity at an exact Lode corner is limited by round-off (~1e-8): theta comes from arccos of a value within ulps of 1 (WP-144)"
date: 2026-11-06
---
### Kernel-vs-oracle TANGENT parity at an exact Lode corner is limited by round-off (~1e-8): theta comes from arccos of a value within ulps of 1 (WP-144)
- **Bites:** on an exactly axisymmetric path (TXC; any three-invariant model whose zeta_y has a corner branch), theta = acos(sqrt(6) y)/3 with sqrt(6) y a few ulps from 1. So |sin 3theta| takes the quantised values 0, 1.5e-8, 2.1e-8, 2.6e-8, ... (sqrt(2k eps)).
  - zeta_y switches between the corner limit (S.9) and the csc 3theta quotient (S.8). They agree only to O(|theta - theta_c|), and the quotient itself cancels catastrophically there (sheet 144a §3.1).
  - The STATE is insensitive (y_a = O(theta - theta_c), so q_a and the converged point do not see it). The consistent TANGENT is not: zeta_y multiplies y_ab ∝ 1/R in q_ab.
  - Measured in the O2 oracle itself: its tangent moves by 1e-9 to 1.8e-8 when its committed elastic strain moves by ONE ulp. The same holds near the vertex (R < 1e-3 |p|: n_hat_ab ∝ 1/R, y_ab ∝ 1/R^2).
  - So a C++ kernel that reproduces the oracle's algorithm exactly, but sums in another order (BLAS vs hand loops, LAPACK eigh vs Jacobi), matches the state to 1e-14 and the off-corner tangent to 1e-12, and still differs by 1.5e-8 in the corner tangent. A uniform 1e-10 tangent gate fails for a non-bug.
- **Rule:**
  - Gate the tangent inside the corner / vertex bands separately (WP-144: |sin 3theta| < 1e-7 or R < 1e-3 |p| → 1e-7). Print the oracle's own 1-ulp self-sensitivity next to the kernel error as the evidence. Keep every STATE quantity at the full gate inside the bands, and keep a generic-theta path (non-coaxial) at the full tangent gate, so the formula is still verified.
  - Measure a per-step INCREMENT (the step's dissipation D) against the path's scale, not its own value: on a near-neutral plastic step dlambda is set by F_tr ≈ 1e-7 kPa, a difference of O(100 kPa) terms.
  - Gate the local / nested iteration counts too: a contract constant changed without changing the answer (the pi_i scan step) shows only there.
- **Workaround/status:** ✅ In force in `Ladruno_files/testbed/norsand_oracle/kernel_parity/test_kernel_parity.py` (WP-144 P1a, 2026-10-01). Nine mutants of the kernel (spin ½, shear-column sum, v vs v0, cap w_eta term, scan step, corner branch, trial threshold, pi* Omega term, substep depth) each fail it.
