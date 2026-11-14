---
wp: WP-232b
title: "uGIMP with fixed particle half-lengths does not converge under large compression"
date: 2026-11-13
---
### uGIMP with fixed particle half-lengths does not converge under large compression

- **Bites:** the explicit MPM self-weight column (Charlton et al. 2017 closed form
  σ = −ρ0·g·(l0 − Y)), compressed 20% (`ρ0·g·l0/E` = 0.2): the oracle's worst particle error is
  0.072 of `ρ0·g·l0` at h = 0.5. Refining makes it **worse**: 0.067 at h = 0.25, and 0.106 at
  h = 0.125. At small strain (`ρ0·g·l0/E` = 0.002) the same column gives 0.013. Found by the WP-232b
  numpy oracle (`tests/mpm_reference.py`, case `column_selfweight_2d_gimp`).
- **Why:** cell-crossing noise. The particle domains of uniform GIMP (uGIMP) keep their initial
  half-length `lp`, so they stop tracking the deformation. Compressed particles overlap, and the
  basis's integration error at every cell crossing does not shrink with h. The linear hat is
  worse still (0.47 at h = 0.5).
- **Workaround/status:** the slice-1 gate for the closed form is the small-deformation column
  (`column_selfweight_small_2d_gimp`, gated at 0.03). The 20% column stays a strict xfail in
  `tests/test_mpm_reference.py`, as the acceptance test for phase 3's remedy: a domain-tracking
  GIMP (cpGIMP / CPDI) or B-spline basis (ADR-232 §9, phase 3). Interface Amendment 4
  (`tests/mpm_interface.md`).
