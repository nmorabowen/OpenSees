---
wp: PR-275
title: "LadrunoCohesiveHingeBiaxial"
pr: "#275, #302"
status: "shipped"
section: "table"
legacy_seq: 95
---
| **LadrunoCohesiveHingeBiaxial** ([[34_ladruno_cohesive_hinge_biaxial_adr]]) — coupled Mz–My cohesive interaction-surface `nDMaterial` (order 2): elliptical onset `√((Mz/Mcz)²+(My/Mcy)²)=1` + isotropic secant damage on the normalized norm, per-axis `Gf`, Benzeggagh-Kenane mode-mix (`-bk η`); reduces EXACTLY to the scalar `LadrunoCohesiveHinge` on each pure axis. Consumed by `LadrunoDispBeamColumn3d -hingeBiaxial`. | nD material | **ND_TAG 33004** | `SRC/material/nD/LadrunoCohesiveHingeBiaxial.{h,cpp}`, `tests/test_ladrunoCohesiveHingeBiaxial_material.py` | shipped | [#275](https://github.com/nmorabowen/OpenSees/pull/275), [#302](https://github.com/nmorabowen/OpenSees/pull/302) |
