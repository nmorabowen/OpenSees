---
wp: WP-115
title: "Vanilla ElasticBeam2d subtracts the ground-motion load Q TWICE with element -mass — UniformExcitation drives its element mass at 2·a_g (found WP-115)"
legacy_seq: 474
---
### Vanilla `ElasticBeam2d` subtracts the ground-motion load Q TWICE with element `-mass` — `UniformExcitation` drives its element mass at 2·a_g (found WP-115)
- **Bites:** a 2D `elasticBeamColumn ... -mass rho` under `pattern UniformExcitation` responds to twice the ground acceleration. Measured WP-115: a 4-element cantilever, constant a_g in y — the tip history with element `-mass` is exactly 2.000× the same beam with the identical lumped masses given as nodal `mass` (ratio 2.0000 at every step). `ElasticBeam3d` is correct; nodal masses are correct.
- **Why:** `ElasticBeam2d::addInertiaLoadToUnbalance` accumulates `Q -= m·a_g` (lumped) or `Q -= M·a_g` (consistent). `getResistingForce()` then does `if (rho != 0) P.addVector(1.0, Q, -1.0);` (`ElasticBeam2d.cpp` ~:1089) AND `getResistingForceIncInertia()` calls `getResistingForce()` and subtracts Q again (~:1010). `ElasticBeam3d::getResistingForceIncInertia` has no second subtraction. Upstream code, present at the fork's import (`30cc727df`); not a Ladruno edit.
- **Workaround/status:** put 2D beam mass on the nodes (`ops.mass`) for ground-motion runs, or use `ElasticBeam3d`/a fork beam. Do NOT use a 2D `elasticBeamColumn` with element `-mass` as a -Q oracle — `tests/test_rayleigh_inertia_bezier_imk.py` uses nodal masses for exactly this reason. **FIXED by WP-119 (#854)**, together with the same defect in `ElasticTimoshenkoBeam2d/3d` — see the entry "subtracted the ground-motion load TWICE with element `-mass`" (WP-119). On builds before WP-119, the nodal-mass workaround above applies.
