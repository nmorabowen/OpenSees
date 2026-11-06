---
wp: PR-1005
title: "1005 -- upstreamable-table row(s)"
date: 2026-10-22
pr: "#1005"
files: ["`SRC/material/nD/UWmaterials/ManzariDafalias.cpp`"]
table: "upstreamable"
---
| `SRC/material/nD/UWmaterials/ManzariDafalias.cpp` | `// Ladruno #2444` — **vanilla's classTag constructor resets the class-wide `mElastFlag = 0`, and `ManzariDafalias3D/PlaneStrain::getCopy()` runs it** (`new ManzariDafalias3D(); *clone = *this;`), so any element that copies an already-copied MD material after `updateMaterialStage -stage 1` sends every live MD instance back to elastic. Vanilla's own full constructors don't reset (so the usual prototype -> `getCopy(type)` path is safe on stock 3.8.0); the fork had made it worse by adding the reset to them too (TIMs item 7), which is what apeGmsh hit: activating MD elements in a later stage returned uy -1.4923e-3 (elastic) instead of -1.6597e-3. Fix: no constructor writes the static; recvSelf already restores it from `data(25)`. Gated by `tests/test_manzari_stage_static_activation.py` + `tests/tcl/manzari_stage_static_activation.tcl`. | [#1005](https://github.com/nmorabowen/OpenSees/pull/1005) |
