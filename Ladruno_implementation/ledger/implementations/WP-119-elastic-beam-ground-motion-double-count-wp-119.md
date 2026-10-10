---
wp: WP-119
title: "Elastic-beam ground-motion double count (WP-119)"
pr: "#854"
status: "fixed on the branch"
section: "table"
legacy_seq: 21
---
| **Elastic-beam ground-motion double count (WP-119)** — vanilla `ElasticBeam2d`, `ElasticTimoshenkoBeam2d`, `ElasticTimoshenkoBeam3d` subtracted the `UniformExcitation` inertia load in both `getResistingForce()` and `getResistingForceIncInertia()`: element mass driven at 2·a_g (rigid-body probe −4.0 vs −2.0). Second subtraction removed (marked vanilla edit, see LEDGER_vanilla_files); `ElasticBeam3d` was already correct; audit found no other element with the double subtraction. **Behaviour change** for these beams with element `-mass` under ground motion. Gate `tests/test_elastic_beam_ground_motion.py` (rigid-body probe + element-vs-nodal-mass differential, ElasticBeam3d as control). | upstream bug fix | — (vanilla elements) | `SRC/element/elasticBeamColumn/{ElasticBeam2d,ElasticTimoshenkoBeam2d,ElasticTimoshenkoBeam3d}.cpp`, `tests/test_elastic_beam_ground_motion.py` | **fixed on the branch** | #854 |
