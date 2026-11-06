---
wp: PR-6
title: "BezierTri6"
pr: "#6, #10, #709, #848"
status: "shipped — unsym-tangent fix (getTangentStiff/getInitialStiff triu-mirror → full…"
section: "table"
legacy_seq: 50
---
| **BezierTri6** — quadratic Bézier triangle element (Kadapa 2018) | Element | 33000 | `SRC/element/bezierTriangle/BezierTri6.{cpp,h}` | shipped — **unsym-tangent fix** (`getTangentStiff`/`getInitialStiff` triu-mirror → full BᵀDB; same defect as BezierTet10, see its row); **WP-114 plane-strain B-bar fix** (TIMs F17): `computeBBarMatrix` used the 3D ÷3 split on the 3-row plane-strain B (material trace (θ+2θ̄)/3, 3 volumetric constraints/element = no relief) → now the 2D ½ split (LadrunoQuad's, 1 constraint/element); ψ=0 punch no longer walls on apex GPs; `tests/test_beziertri6_bbar_plane_strain.py` | [#6](https://github.com/nmorabowen/OpenSees/pull/6), [#10](https://github.com/nmorabowen/OpenSees/pull/10), [#709](https://github.com/nmorabowen/OpenSees/pull/709), [#848](https://github.com/nmorabowen/OpenSees/pull/848) |
