---
wp: WP-213
title: "A static-correction vector alone does not reach 100 % modal mass; a fixed pseudo-mode period must not go through lambda"
date: 2026-10-26
pr: "#1023"
---
### Static-correction (missing-mass) mode shows < 100 % in `modalProperties`; a fixed period through λ rescales the RSA amplitude

- **Bites:** a single residual vector x = K⁻¹Mι deflated against the modes
  (the textbook static correction, NASTRAN RESVEC, the WP-213 guide's first
  design) reproduces the static response exactly, but its effective mass is
  (xᵀMr)²/(xᵀMx·rᵀMr) of the missing mass r: 93.8 % on the V3 podium model.
  Separately, giving a pseudo-mode a fixed period by setting λ = (2π/T)²
  scales the RSA displacement φΓSa/λ by λ_Rayleigh/λ_T as well as changing the
  Sa lookup: V3 with T = 0.05 s gave 0.085 × the reference base shear.
- **Why:** Cauchy–Schwarz (equality only if x ∥ r); and stock RSA uses the
  mode's λ both for the Sa lookup and for the 1/ω² amplitude.
- **Fix/status:** `residualModes` uses span{r, K⁻¹Mr}, which gives 100 % mass
  and the exact static response (ADR 213 D3). `-period $T` keeps the Rayleigh
  λ and records T as a per-mode Sa-period override that
  `responseSpectrumAnalysis` reads on every path, with ω = 2π/T in CQC too so
  the overridden modes act as one rigid mode (ADR 213 D3b; keeping their own
  ω made them combine as SRSS, 1.278 × instead of 1.651 × on the review's
  stiff podium). Even so, `-period 0.05` is only conservative when the
  pseudo-modes really are rigid (Rayleigh periods below 0.05 s): on V3 as
  specified (0.31 s, 0.16 s) it gives 0.662 × the reference, and
  `residualModes` warns. Learned 2026-10-26
  (WP-213).
