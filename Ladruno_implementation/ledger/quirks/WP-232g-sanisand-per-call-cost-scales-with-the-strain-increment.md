---
wp: WP-232g
title: "LadrunoSANISAND's per-call cost scales with the strain increment"
date: 2026-11-15
---
### LadrunoSANISAND's per-call cost scales with the strain increment

- **Bites:** any cost figure for SANISAND per material call (a cost model, a throughput
  projection, a "SANISAND is N× DruckerPrager" ratio). The cost depends on the strain increment
  per call, not just on the material. Driver: a zero-free-DOF stdBrick at p ≈ 10 kPa, isochoric
  steps, default options (WP-232g M2 smoke, build `9db5803f6`). Plastic per-call cost at ×0.5,
  ×1 and ×2 of a 3.4e-5 strain-increment norm:

  | Material | ×0.5 | ×1 | ×2 |
  |---|---|---|---|
  | LadrunoSANISAND | 8.0e4 ns | 1.45e5 ns | 2.67e5 ns (about linear) |
  | DruckerPrager | 1.85 µs | 1.77 µs | 1.78 µs (flat) |

  The SANISAND / DP ratio therefore runs 43 → 82 → 150 across the band. The explicit integrator
  sub-steps in proportion to the increment.
- **Workaround:** quote a SANISAND cost only together with its strain increment per call. To
  project an MPM or FE run, match the increment to the target run's per-step strain increment
  (`Ladruno_scripts/mpm_m0b/m2_explicit.py` measures the column's and drives the point at it).
  A smaller explicit Δt makes each call cheaper, but not the run, which takes more calls.
