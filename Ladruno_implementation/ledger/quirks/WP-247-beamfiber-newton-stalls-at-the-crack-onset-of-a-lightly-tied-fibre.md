---
wp: WP-247
title: "LadrunoRCConcrete BeamFiber view: the lateral Newton stalls at the crack onset of a lightly tied fibre (rho fyt < f_cr), though a root sits next to the committed strain"
pr: "TBD"
date: 2026-11-22
---
### LadrunoRCConcrete BeamFiber view: the lateral Newton stalls at the crack onset of a lightly tied fibre (ρ fyt < f_cr), though a root sits next to the committed strain (WP-247 adds Model C's root scan)
- **Bites:** a force-based fibre column of the program's lab set (Sezen–Moehle 2CLD12, `LadrunoNDFiberSection3d` with parabolic shear, core fibres with `-ties 0.0017 0 469 2e5`) stopped at 0.22 % drift. The step where its first core fibre cracked in shear was refused: "BeamFiber lateral balance against the ties did not converge", εxx = −1.2e-4, γxy = −2.1e-4. A single fibre reproduces it on a monotonic ramp, but only for some step sizes: 60 steps to that strain refuse on the last step, and 300 steps to a larger strain never do.
- **Why:** past the tension peak, the concrete's `σyy` softens faster than `ρ Es` stiffens it, so `r_y(εyy) = σyy + ρ f_s` is not monotone. A full Newton step from the uncracked side lands on the falling branch, where `∂r/∂εyy < 0`, and cycles there. Backtracking on `|r|` only shortens the step along the same bad direction. Model C finds the root, `εyy = 3.5e-6`, by scanning, not by Newton.
- **Fix/status (WP-247):** when Newton fails, `driveBeamFiber` applies Model C's rule:
  - It samples the committed lateral normal ± geometric offsets, along y and then z.
  - It bisects the increasing sign change nearest the committed value.
  - It balances the other two lateral unknowns by a 2×2 Newton at each sample.

  The result is accepted on the stalled-Newton balance test. Newton-reached states are unchanged.
- **Second symptom, same root cause:** a cover fibre that is fully cracked in y has a nearly flat `σyy`. Its lateral Newton "converged" to `εyy ≈ 3e4` with a small residual (Sezen–Moehle 2CHD12, high axial load) and refused 250 steps later. WP-247 projects the lateral unknowns into Model C's box `[−0.01, 0.15]`.
- **Gates:** `tests/_testbed/rc_rootscan_gpp.cpp` R1–R4 (the scan) and R5 (the box). See [[LadrunoRCConcrete_guide]] §12.1.
