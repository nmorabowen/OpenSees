---
wp: ADR-25
title: "LadrunoLST"
pr: "#554"
status: "shipped — P3, std-only (std × -geom {linear|finite}; reduce-to SixNodeTri + ran…"
class_tags: ["ELE_TAG_LadrunoLST", "ND_TAG_LogStrain2D"]
section: "table"
legacy_seq: 83
---
| **LadrunoLST** — 6-node linear-strain triangle (T6), the F-bar-friendly finite-strain triangle ([[70_ladruno_plane_finite_triangles_adr]]): std × `-geom {linear\|finite}` + 3-point interior integration (matches upstream `SixNodeTri` = the reduce-to gate, mirroring Quad↔FourNodeQuad / CST↔Tri31). The linear strain field carries an inclined shear band (unlike constant-strain `LadrunoCST`, which snaps bands to mesh lines). The ADR's "element-level F-bar can average the non-constant J" premise was REFUTED at P3 (rank-deficient — see status cell); the usable near-incompressible triangle still awaits P4. Consumes the shared `LadrunoFiniteStrain2D` kernel. corot stays out of scope (ADR-25 P4 `SolidTransformation2DCorot`). | Element | **`ELE_TAG_LadrunoLST` 33016** (ELE registry; per-registry, distinct from `ND_TAG_LogStrain2D`=33016) | `SRC/element/ladrunoPlane/{LadrunoLST.{cpp,h},OPS_LadrunoLST.cpp}` | shipped — **P3, std-only** (std × `-geom {linear\|finite}`; reduce-to `SixNodeTri` + rank/RBM + linear-strain-field patch + finite oracle gates in `tests/test_ladrunolst_element.py` + `tests/test_ladrunolst_finite.py`). **bbar/F-bar REFUTED at P3**: constant element-mean dilatation is rank-deficient on the T6 — the 2 quadratic conformal modes (Re/Im z²) carry zero dev strain and zero MEAN dilatation ⇒ 5 zero-energy modes on a free element (stacked-B̄ rank 7 of 9; caught by the locked T0 zero-energy gate, confirmed by numpy rank + compiled eigen). `-formulation bbar` parser-refused; triangle volumetric cure = P4 (nodal F-bar-Patch / P1-projected dilatation) | [#554](https://github.com/nmorabowen/OpenSees/pull/554) |
