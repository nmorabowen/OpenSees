---
wp: WP-183
title: "eigen under Penalty / Auto / Lagrange constraints returns ~1e-3 or 1.0 eigenvalues (fixed by WP-183)"
date: 2026-10-11
status: "fixed-by WP-183"
---
### `eigen` after gravity (or after a Linear analysis) returns eigenvalues ~ +-1e-3 under `constraints Penalty` / `Auto` / `Lagrange`; `-fullGenLapack` / `-feast` return lambda ~ 1 even at K0; `modalDamping` built on such an eigen is ~50 % off (fixed by WP-183)
- **Bites:** any model with equalDOF / rigidLink / fix under a Penalty or Lagrange handler, and `constraints Auto` (it makes a PenaltyMP_FE for every MP). Transformation and Plain are unaffected. ARPACK at K0 often looks right, so the error appears only after the first analysis step that evaluates the residual at a deformed state (Newton; Linear after its first step).
- **Why (two defects, upstream code identical on master):** (A) the constraint FEs' `getM_Force` returned `*resid`, the last constraint residual, so ARPACK's `M x` (used when M is not diagonal) and `IncrementalIntegrator::doMv` (modalDamping) added `-alpha C^T C (u - u0)`; (B) the eigen mass loops added every FE's `getTangent`, and a constraint FE ignores the integrator and returns its penalty stiffness, so `alpha C^T C` / `alpha` went into M (that is what made M non-diagonal and fed (A), and what explicit-M solvers saw directly).
- **Fix:** WP-183 (`ledger/vanilla/WP-183-*`): zero `getM_Force` / `getC_Force` on Penalty/Lagrange MP/EQ/SP FEs; the mass loops skip FEs with no Element. Rayleigh transients are bitwise unchanged.
- **Separate, not fixed:** a penalty value must be compared with the largest stiffness it fights, rotational ones included. In N-mm a 2D column's base `4EI/L` is ~1.7e11, so `Penalty 1e12 1e12` on the `fix` is a semi-rigid base (10 % eigenvalue shift, 30 % transient shift). Test: `tests/test_wp183_constraint_eigen_mass.py`.
