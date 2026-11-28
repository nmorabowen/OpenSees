---
wp: WP-267
title: "A DDM run is only as exact as the tangent the material reports: LadrunoBondSlip reports k0, not 0, on its plateau and slider"
pr: "#1176"
date: 2026-11-28
status: "open"
---
### `sensNodeDisp` disagrees with the finite difference although the material's own DDM check passes: the step sat on a branch whose reported tangent is not the consistent one

- **Bites:** a DDM run (`parameter` + `sensitivityAlgorithm`) through a material whose
  `getTangent` is deliberately not d(stress)/d(strain) on some branch. The structure
  solves `K dU/dθ = dP/dθ − dF/dθ|cond` with the tangent the materials report, so on
  that branch `dU/dθ` is wrong even when `getStressSensitivity` and
  `commitSensitivity` are exact. In `LadrunoBondSlip` (WP-267) that is the plateau
  (`s1 < |s| <= s2`) and the residual branch (true slope 0, reported `k0` and
  `−1e-12·k0` "to avoid a hard 0 pivot"), and friction sliding at the cap (true slope
  0, reported `k0`). The material's committed history sensitivities stay exact; the
  error is in `dU/dθ`, and through the next `commitSensitivity` it enters the history.
- **Why:** `LadrunoBondSlip::setTrialStrain` (now `LadrunoBondSlipKernel.h::trialState`)
  replaces a zero tangent so the global matrix never gets a zero pivot. Newton still
  converges with that tangent (as a modified Newton). The DDM solve does not iterate,
  so it keeps the wrong tangent.
- **What to do:** keep a DDM-checked run off those branches (the WP-267 end-to-end
  test does: `tests/test_ladrunoBondSlip_ddm.py`), or make the tangent consistent
  first. That is a production tangent change (worker rule 7: behind an option, with
  numbers), not part of WP-267. The kernel check pins the guard per branch, so a
  change to it is visible: `tests/test_ladrunoBondSlip_dual_cpp.py::test_tangent_is_consistent_off_the_plateau_and_the_slider`.
  `LadrunoCohesiveHinge` has no such branch: its loading tangent is the envelope
  slope, its unload tangent the secant, and the IMPL-EX tangent the frozen secant,
  each the derivative of the stress it returns.
