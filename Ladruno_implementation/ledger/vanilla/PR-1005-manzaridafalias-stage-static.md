---
wp: PR-1005
title: "ManzariDafalias: no constructor resets the class-wide stage flag; wipe does (apeGmsh #2444)"
date: 2026-10-22
pr: "#1005"
files: ["`SRC/material/nD/UWmaterials/ManzariDafalias.cpp`", "`SRC/material/nD/UWmaterials/ManzariDafalias.h`", "`SRC/material/nD/NDMaterial.cpp`"]
---
| `SRC/material/nD/UWmaterials/ManzariDafalias.cpp` | `// Ladruno #2444`: file-scope `mElastFlag` initialiser `1` -> `0` (elastic default); the `mElastFlag = 0` lines removed from both full constructors (fork TIMs item 7, 41cfb81de) and from vanilla's classTag constructor (the one `getCopy()` runs); free function `ladrunoManzariDafaliasResetStage()` for the wipe hook. | [#1005](https://github.com/nmorabowen/OpenSees/pull/1005) |
| `SRC/material/nD/UWmaterials/ManzariDafalias.h` | `// Ladruno #2444`: public `static void resetStageOnWipe()` (sets `mElastFlag = 0`); comment on `mElastFlag` corrected (0 = elastic). | [#1005](https://github.com/nmorabowen/OpenSees/pull/1005) |
| `SRC/material/nD/NDMaterial.cpp` | `// Ladruno #2444`: `OPS_clearAllNDMaterial()` calls `ladrunoManzariDafaliasResetStage()` (local extern, same idiom as WP-104 / #2044). | [#1005](https://github.com/nmorabowen/OpenSees/pull/1005) |
