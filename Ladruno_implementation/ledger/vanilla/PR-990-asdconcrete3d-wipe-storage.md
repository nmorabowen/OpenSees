---
wp: PR-990
title: "ASDConcrete3D: wipe clears the per-tag hardening-law registry (apeGmsh #2044)"
date: 2026-10-15
pr: "#990"
files: ["`SRC/material/nD/ASDConcrete3DMaterial.h`", "`SRC/material/nD/ASDConcrete3DMaterial.cpp`", "`SRC/material/nD/NDMaterial.cpp`"]
---
| `SRC/material/nD/ASDConcrete3DMaterial.h` | `// Ladruno #2044`: declares `clear()` on `HardeningLawStorage` and `CrackPlanesStorage`. | [#990](https://github.com/nmorabowen/OpenSees/pull/990) |
| `SRC/material/nD/ASDConcrete3DMaterial.cpp` | `// Ladruno #2044`: the two `clear()` bodies and the free function `ladrunoAsdConcrete3DResetStorage()` that clears both registries on `wipe`. | [#990](https://github.com/nmorabowen/OpenSees/pull/990) |
| `SRC/material/nD/NDMaterial.cpp` | `// Ladruno #2044`: `OPS_clearAllNDMaterial()` calls `ladrunoAsdConcrete3DResetStorage()` (local extern, same idiom as WP-104). | [#990](https://github.com/nmorabowen/OpenSees/pull/990) |
