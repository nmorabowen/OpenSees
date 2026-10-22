---
wp: WP-204
title: "A diaphragm connected at points: the displacement at a point load or point support never converges with mesh refinement"
pr: "#1013"
date: 2026-10-22
---
### Diaphragm (membrane) displacement at a column or a point load changes with the mesh and never converges — in ETABS too

- **What bites.** Validating `LadrunoMembrane` against a refined ASDShellQ4 membrane mesh: the displacement AT a loaded corner of an 8 x 6 panel went 4.21e-5 -> 4.85e-5 -> 5.48e-5 as the reference spacing halved (0.1, 0.05, 0.025); a "1-2 % vs the refined mesh" target on such a displacement is unreachable for ANY model.
- **Why.** In 2-D elasticity the displacement under a point force is logarithmically singular (and a point support becomes infinitely flexible as h -> 0). Every model that connects a slab to a column at one node (ETABS, OpenSees, this element) has a point-connection stiffness that depends on its mesh size.
- **What converges.** The reaction SPLIT between supports (0.4-2.6 % over `-meshSize` 1.0 .. 0.25 on the deep-beam case), resultants of line supports (walls with many nodes), and drifts dominated by panel deformation. Compare displacements only at matched resolution.
- **How to apply.** To reproduce an ETABS diaphragm, set `-meshSize` to the ETABS slab auto-mesh size; connect walls through their mesh nodes (hole loops or `-diaphragmNodes`), not one node. Validation tests: compare reaction splits and line resultants, not point displacements (`tests/test_wp204_membrane_diaphragm.py`).
- **Related trap.** `build_pslg_cdt_mesh` caps the vertex count (default 1200) and silently coarsens beyond it; the diaphragm mesh now passes 6000. A too-small cap showed up as an asymmetric, wrong answer on a symmetric 12 x 12 slab at `-meshSize 0.25`.
