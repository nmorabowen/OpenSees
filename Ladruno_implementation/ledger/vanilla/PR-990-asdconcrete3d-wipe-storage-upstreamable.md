---
wp: PR-990
title: "990 -- upstreamable-table row(s)"
date: 2026-10-15
pr: "#990"
files: ["`SRC/material/nD/ASDConcrete3DMaterial.cpp`", "`SRC/material/uniaxial/ASDConcrete1DMaterial.cpp`"]
table: "upstreamable"
---
| `SRC/material/nD/ASDConcrete3DMaterial.cpp` | `// Ladruno #2044` — **`HardeningLawStorage` survived `wipe`.** The registry is process-wide, keyed only by material tag, and store-if-absent; `HardeningLaw::regularize()` → `deRegularize()` copies the stored law over the instance. So under `-autoRegularization` with `lch != lch_ref`, a same-tag material defined after a `wipe` silently took the FIRST model's backbones, elastic modulus included (E=3e6 then E=30e6: sigma_zz = -720 for -7200, Newton diverges in an elastic step). Now cleared from `OPS_clearAllNDMaterial()` (`NDMaterial.cpp` hook; upstream would call it from the same place). Gated by `tests/test_asdconcrete3d_wipe_storage.py` + `tests/tcl/asdconcrete3d_wipe_storage.tcl` (both fail on 0c4c6a1d). | [#990](https://github.com/nmorabowen/OpenSees/pull/990) |
| `SRC/material/uniaxial/ASDConcrete1DMaterial.cpp` | `// Ladruno #2044` — **the identical store-if-absent `HardeningLawStorage` survived `wipe`** (E=3e6 then E=30e6 under one tag: singular truss, u = 0, stress 0). Cleared from `OPS_clearAllUniaxialMaterial()`. Gated by `tests/test_asdconcrete1d_wipe_storage.py` + `tests/tcl/asdconcrete1d_wipe_storage.tcl` (both fail on 0c4c6a1d). | [#990](https://github.com/nmorabowen/OpenSees/pull/990) |
