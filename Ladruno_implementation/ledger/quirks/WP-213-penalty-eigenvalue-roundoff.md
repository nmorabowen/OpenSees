---
wp: WP-213
title: "Under Penalty the eigenvalue carries round-off ~ eps*alpha/(lambda*m); the Rayleigh quotient does not"
date: 2026-10-26
pr: "#1023"
---
### Penalty model: eigenvalue and Rayleigh quotient of the same mode differ by 1e-4 on an unchanged model

- **Bites:** any check that compares an eigenvalue from `eigen` with φᵀKφ/φᵀMφ
  recomputed on the same tangent. Penalty constraints with α = 1e14 (two
  nodes per floor tied by `equalDOF`, kip-in units) showed relative gaps of
  1.6e-4 (ARPACK) to 2.7e-4 (`-fullGenLapack`); in kN-m about 1e-6.
  WP-213's first state check refused such unchanged models as "changed after
  eigen".
- **Why:** the eigen solve works on a matrix whose entries span α to k, so
  its eigenvalue carries round-off of order ε·|φ|ᵀ|K||φ| / φᵀKφ ≈ ε·α/(λ m).
  The Rayleigh quotient of the returned vector is second-order in the vector
  error, so it is exact.
- **Fix/status:** scale any such tolerance by that floor. `residualModes`
  refuses above max(1e-4, 100 f) and warns above max(1e-7, 10 f), with f
  accumulated in the K·Φ element sweep (`LadrunoRitzKernel::applyK` absQuad).
  Learned 2026-10-26 (WP-213, review round 2).

#### Follow-up (WP-214)

The 100 f / 10 f margins were loose: on the tied-floor model (kip-in, 10
storeys) f₁ = 4.0e-3 while the measured gap is only 2.5e-4, so after `eigen 1`
a 20 % stiffness loss only warned and a 2 % loss passed silently. WP-214 sets
refuse above max(1e-4, 10 f) and warn above max(1e-7, f); the unchanged
models (kip-in and kN-m, `-fullGenLapack` and ARPACK; reference LAPACK's gap
≈ 2.1e-6) stay well under f.
