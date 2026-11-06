---
wp: WP-117
title: "Bezier ground-motion sign fix (WP-117)"
pr: "#852"
status: "fixed on the branch"
section: "table"
legacy_seq: 22
---
| **Bezier ground-motion sign fix (WP-117)** — `BezierTri6` / `BezierTet10` `addInertiaLoadToUnbalance` built `Q += +M·R·a_g` while `getResistingForce()` subtracts `Q`, so `UniformExcitation` drove the Bezier element mass with −a_g (mesh shaken the wrong way; mixed Bezier-soil + nodal-mass-structure models shaken in opposite directions). Factor `+1.0 → −1.0`, the vanilla and fork-wide convention; audit of all 14 fork `addInertiaLoadToUnbalance` implementations found no other. **Behaviour change** for every Bezier ground-motion run. Gate `tests/test_bezier_ground_motion.py`: rigid-body relative accel == −a_g over {lumped, -cMass} × {element, material ρ}, and UniformExcitation == equivalent −m·a_g nodal loads; all 10 cases fail on the unfixed binary. | bug fix | — (existing BezierTri6 / BezierTet10 tags) | `SRC/element/bezierTriangle/BezierTri6.cpp`, `SRC/element/bezierTetrahedron/BezierTet10.cpp`, `tests/test_bezier_ground_motion.py` | **fixed on the branch** | #852 |
