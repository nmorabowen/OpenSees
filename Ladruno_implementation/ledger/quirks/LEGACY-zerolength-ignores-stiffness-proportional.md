---
wp: LEGACY
title: "zeroLength ignores stiffness-proportional Rayleigh unless -doRayleigh 1"
legacy_seq: 15
---
### zeroLength ignores stiffness-proportional Rayleigh unless `-doRayleigh 1`
- **Bites:** a `zeroLength` / `zeroLengthSection` element contributes **zero**
  stiffness-proportional Rayleigh damping (`betaK`, `betaKinit`/`betaK0`,
  `betaKcomm`/`betaKc`) by default. You set `rayleigh 0 0 0.0159 0`, expect
  ζ≈0.05, and measure ζ≈0. Mass-proportional `alphaM` (which lives on the node,
  not the element) works regardless, which masks the problem.
- **Why:** the element carries an internal `doRayleigh` flag, default **0**, that
  gates whether `getDamp()`/`getResistingForceIncInertia()` include the element's
  stiffness term. The `-doRayleigh` option flips it: `element zeroLength … -dir 1
  -doRayleigh 1`. Most other elements default the flag on; zeroLength does not.
- **Workaround/status (2026-06-01):** pass `-doRayleigh 1` whenever you want
  stiffness-proportional Rayleigh on a zeroLength; or (better) model the damping
  physically with a `Viscous`/`ViscousDamper` uniaxial material on the DOF — that
  enters R(u̇) directly and is explicit-safe. Pinned by
  `tests/test_damping_channels.py::test_zeroLength_doRayleigh_default_off`
  (ζ≈0 with default) and `::test_betaK0_realises_target_zeta` (ζ=0.05 with flag).
  Full map: [[12_damping_channels]].
