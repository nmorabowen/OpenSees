---
wp: PR-264
title: "Stage 2 (cohesive material SHIPPED #264): LadrunoCohesiveHinge MAT_TAG 33003 (uniaxial) BUILT — discrete cohesive momen…"
section: "history"
legacy_seq: 119
---
- **Stage 2 (cohesive material SHIPPED [#264](https://github.com/nmorabowen/OpenSees/pull/264)):** `LadrunoCohesiveHinge` `MAT_TAG 33003` (uniaxial) **BUILT** — discrete cohesive moment–rotation law `M([[θ]])` carrying `Gf` per hinge (rigid-softening penalty w/ guarded floor `Mc²/2Gf`; exp/linear envelope calibrated to `∫M d[[θ]]==Gf`; irreversible secant damage; `getEnergy()`). `SRC/material/uniaxial/LadrunoCohesiveHinge.{h,cpp}`, `tests/test_ladrunoCohesiveHinge_material.py` (10/10, energy gate exact 1e-9 LINEAR).
