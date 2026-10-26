---
wp: WP-213
title: "A static-correction vector alone does not reach 100 % modal mass; -period fixes lambda, not only Sa"
date: 2026-10-26
---
### Static-correction (missing-mass) mode shows < 100 % in `modalProperties`; `-period 0.05` changes the amplitude

- **Bites:** a single residual vector x = K⁻¹Mι deflated against the modes
  (the textbook static correction, NASTRAN RESVEC, the WP-213 guide's first
  design) reproduces the static response exactly but its effective mass is
  (xᵀMr)²/(xᵀMx·rᵀMr) of the missing mass r — 93.8 % on the V3 podium model.
  Separately, giving a pseudo-mode a fixed period through λ = (2π/T)² scales
  the RSA displacement φΓSa/λ by λ_Rayleigh/λ_T, not only the Sa lookup:
  V3 with `-period 0.05` gives 0.085 × the reference base shear.
- **Why:** Cauchy–Schwarz (equality only if x ∥ r); and RSA uses the same λ
  for the Sa lookup and the 1/ω² amplitude.
- **Fix/status:** `residualModes` uses span{r, K⁻¹Mr} (100 % mass and the
  static response, ADR 213 D3). `-period $T` is documented as conservative
  only when the Rayleigh periods are already below T; the ASCE rigid-mode
  rule proper needs an RSA-side period override (ADR 213 follow-up 1).
  Learned 2026-10-26 (WP-213).
