---
wp: ADR-31
title: "LadrunoConcrete3D accepted a failed return map: setTrialStrain returned 0, so nothing cut the step"
date: 2026-10-29
---
### LadrunoConcrete3D's "step-cut" warning never cut a step: a failed return map returned 0

- **Bites:** a softening `LadrunoConcrete3D` run prints millions of
  `return map did not converge ... step-cut` warnings, every step reports
  converged, and the peak load moves with the BLAS rounding. apeGmsh #2602:
  the Cook-Mitchell double corbel peaked at 287, 342 or 379 kN per side
  depending on `MKL_CBWR`, and a 1e-12 change in `ft` on one MKL path moved
  it from 318 to 414 kN.
- **Why:** when the kernel's return map fails it hands back the elastic
  predictor with status 2 (`LadrunoConcrete3DKernel.h`, the PR #249 safe
  fallback), and the undamaged elastic tangent. `setTrialStrain` threw the
  status away and returned 0 on all three paths (3D / plane strain, the
  `sigma_22 = 0` condensation loop, the confined BeamFiber view). The residual
  then jumps wherever a Gauss point flips between a converged return and the
  trial stress, Newton stalls, and a fallback accepts a state that a
  last-bit change in strain decides. Passing the +2 through would not have
  helped: `DisplacementControl` and `LoadControl` refuse only `< 0`, and
  `LadrunoBrick` only the exact sentinel.
- **The failure set is fractal.** From the virgin state, 13% of the random
  plane-strain increments of size 1e-4 that go plastic fail (4.5% of all of
  them; the battery's material), and a 1e-6 relative change in one component
  flips some of them. A test that drives a failure must pick an
  input whose outcome survives every perturbation it can see (the #2617 tests
  check 400 random ones from 1e-12 to 1e-4).
- **Fix (apeGmsh #2617):** `setTrialStrain` returns `LADRUNO_MATERIAL_REFUSED`
  when any `integrate()` in the call failed, on all three paths. In the
  condensation loop the failure accumulates over the passes, per the
  maintainer ruling on #2602: a pass that fell back fed its elastic
  `Dtan6[2][2]` into the Newton on `eps_22`. Measured on the kernel at
  increments of 3e-4: 10.5% of plane-stress trials fail on the last pass and
  a further 36.8% only on an earlier pass, so accumulation refuses about four
  times as many plane-stress trials as a last-pass check would.
- **This departs from the ADR-33/34 rule** that a material's best-effort
  miss returns 0 and lets the global Newton arbitrate. That rule holds for a
  best state; LadrunoConcrete3D's fallback is not one (it did not integrate
  the increment, and its committed state is untouched), which is the case the
  sentinel exists for (LadrunoMaterialStatus.h, (b)).
- **Hosts:** the code reaches the analysis only through a FORWARD or SENTINEL
  element of the WP-99 refusal roster. On a DISCARD host (`stdBrick`,
  `SSPquad`, `FourNodeTetrahedron`, `LadrunoSolidShell`, ...) the failed
  return is still accepted silently.
