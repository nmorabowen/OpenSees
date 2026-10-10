---
wp: PR-264
title: "LadrunoCohesiveHinge"
pr: "#264"
status: "shipped"
section: "table"
legacy_seq: 94
---
| **LadrunoCohesiveHinge** ([[32_ladruno_dispbeamcolumn_regularization_adr]]) — discrete cohesive moment–rotation `UniaxialMaterial` `M([[θ]])` carrying `Gf` per hinge (rigid-softening penalty w/ guarded floor `Mc²/2Gf`; exp/linear envelope calibrated to `∫M d[[θ]]==Gf`; irreversible secant damage; `getEnergy()`) — the cohesive law consumed by the DispBeamColumn `-hinge`. | Uniaxial material | **MAT_TAG 33003** | `SRC/material/uniaxial/LadrunoCohesiveHinge.{h,cpp}`, `tests/test_ladrunoCohesiveHinge_material.py` | shipped | [#264](https://github.com/nmorabowen/OpenSees/pull/264) |
