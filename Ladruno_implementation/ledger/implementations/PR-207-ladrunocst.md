---
wp: PR-207
title: "LadrunoCST"
pr: "#207"
status: "shipped — P1 (std, -geom linear); mirrors Tri31 (covered in the 13/13 plane bat…"
section: "table"
legacy_seq: 80
---
| **LadrunoCST** — unified 3-node constant-strain triangle, the thin 2D sibling of `LadrunoQuad`: `-formulation std` only (1-pt triangle is rank-sufficient — no hourglass; bbar/ssp/eas have nothing to average) + `-geom {linear\|corot\|finite}` + crack-band `lch=√(2A)`. std↔upstream `Tri31`. Ships as the trivial baseline / triangular-mesh fallback / future E-FEM carrier ([[53_ladruno_embedded_discontinuity_adr]]) — **honestly low-value** (plain CST volumetric-locks and mesh-biases localization; see [[26_ladruno_plane_frontier_adr]] §CST). ADR [[25_ladruno_plane_elements_adr]]. | Element | **ELE_TAG 33008** (after LadrunoQuad=33007; renumbered from the ADR-reserved 33007) | `SRC/element/ladrunoPlane/{LadrunoCST.{cpp,h},OPS_LadrunoCST.cpp}`, `tests/test_ladrunoPlane_element.py` | shipped — **P1** (std, `-geom linear`); mirrors `Tri31` (covered in the 13/13 plane battery) | [#207](https://github.com/nmorabowen/OpenSees/pull/207) |
