---
wp: WP-119
title: "Vanilla ElasticBeam2d / ElasticTimoshenkoBeam2d / ElasticTimoshenkoBeam3d subtracted the ground-motion load TWICE with element -mass — UniformExcitation drove…"
legacy_seq: 476
---
### Vanilla `ElasticBeam2d` / `ElasticTimoshenkoBeam2d` / `ElasticTimoshenkoBeam3d` subtracted the ground-motion load TWICE with element `-mass` — `UniformExcitation` drove them at 2·a_g (FIXED WP-119)
- **Bites:** a beam with element `-mass` under `pattern UniformExcitation` responds to twice the ground acceleration. Rigid-body probe (beam free only in the excitation direction, a_g = 2.0): relative acceleration −4.000 for these three where every other element gives −2.000; a cantilever with element `-mass` moves exactly 2× the same cantilever with the identical lumped masses given as nodal `mass`. `ElasticBeam3d` and nodal masses were always correct.
- **Why:** the load vector (`Q` / `theLoad`) holds only the inertia load `−M·R·a_g`. `getResistingForce()` subtracts it when `rho != 0`, and `getResistingForceIncInertia()` calls `getResistingForce()` and then subtracted it again. Upstream code, still in OpenSees master (checked 2026-09-24).
- **Workaround/status:** FIXED WP-119 (#854): second subtraction commented out, marked `// Ladruno WP-119` (LEDGER_vanilla_files; upstreamable). Audit of every element's `getResistingForce` / `getResistingForceIncInertia` pair found no other. General check for any element: the load vector must be subtracted ONCE — if `IncInertia` calls `getResistingForce()`, it must not subtract the load again (LadrunoBrick/SolidShell avoid it by rebuilding the residual without calling `getResistingForce()`). Gate `tests/test_elastic_beam_ground_motion.py`. On builds before WP-119, put 2D beam mass on the nodes for ground-motion runs.
