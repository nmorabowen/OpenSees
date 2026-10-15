---
wp: PR-990
title: "ASDConcrete3D + ASDConcrete1D: wipe clears the per-tag hardening-law registry (apeGmsh #2044)"
date: 2026-10-15
pr: "#990"
files: ["`SRC/material/nD/ASDConcrete3DMaterial.h`", "`SRC/material/nD/ASDConcrete3DMaterial.cpp`", "`SRC/material/nD/NDMaterial.cpp`", "`SRC/material/uniaxial/ASDConcrete1DMaterial.h`", "`SRC/material/uniaxial/ASDConcrete1DMaterial.cpp`", "`SRC/material/uniaxial/UniaxialMaterial.cpp`"]
---
| `SRC/material/nD/ASDConcrete3DMaterial.h` | `// Ladruno #2044`: declares `clear()` on `HardeningLawStorage` and `CrackPlanesStorage`. | [#990](https://github.com/nmorabowen/OpenSees/pull/990) |
| `SRC/material/nD/ASDConcrete3DMaterial.cpp` | `// Ladruno #2044`: the two `clear()` bodies and the free function `ladrunoAsdConcrete3DResetStorage()` that clears both registries on `wipe`. | [#990](https://github.com/nmorabowen/OpenSees/pull/990) |
| `SRC/material/nD/NDMaterial.cpp` | `// Ladruno #2044`: `OPS_clearAllNDMaterial()` calls `ladrunoAsdConcrete3DResetStorage()` (local extern, same idiom as WP-104). | [#990](https://github.com/nmorabowen/OpenSees/pull/990) |
| `SRC/material/uniaxial/ASDConcrete1DMaterial.h` | `// Ladruno #2044`: declares `clear()` on `HardeningLawStorage`. | [#990](https://github.com/nmorabowen/OpenSees/pull/990) |
| `SRC/material/uniaxial/ASDConcrete1DMaterial.cpp` | `// Ladruno #2044`: the `clear()` body and the free function `ladrunoAsdConcrete1DResetStorage()`. | [#990](https://github.com/nmorabowen/OpenSees/pull/990) |
| `SRC/material/uniaxial/UniaxialMaterial.cpp` | `// Ladruno #2044`: `OPS_clearAllUniaxialMaterial()` calls `ladrunoAsdConcrete1DResetStorage()` (local extern). | [#990](https://github.com/nmorabowen/OpenSees/pull/990) |
