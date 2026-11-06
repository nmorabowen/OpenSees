---
wp: PR-97
title: "97 -- 5 vanilla row(s)"
pr: "#97"
files: ["`SRC/classTags.h`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/interpreter/OpenSeesNDMaterialCommands.cpp`", "`SRC/material/nD/CMakeLists.txt`", "`SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}`"]
table: "main"
legacy_seq: [94, 95, 102, 103, 113]
---
| `SRC/classTags.h` | Register `ND_TAG_LadrunoJ2Finite`=33012 (finite-strain-native combined-hardening J2 with co-rotating backstress) | [#97](https://github.com/nmorabowen/OpenSees/pull/97) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno`: `#include "LadrunoJ2Finite.h"` + `case ND_TAG_LadrunoJ2Finite: return new LadrunoJ2Finite();` so parallel/database `recvSelf` can reconstruct it | [#97](https://github.com/nmorabowen/OpenSees/pull/97) |
| `SRC/interpreter/OpenSeesNDMaterialCommands.cpp` | `// Ladruno`: fwd-decl `OPS_LadrunoJ2Finite()` + `nDMaterialsMap["LadrunoJ2Finite"]` (shared OpenSeesPy/openseesmp/interpreter-Tcl registry) | [#97](https://github.com/nmorabowen/OpenSees/pull/97) |
| `SRC/material/nD/CMakeLists.txt` | Add `LadrunoJ2Finite.cpp` to `OPS_Material` sources + `LadrunoJ2Finite.h` to headers | [#97](https://github.com/nmorabowen/OpenSees/pull/97) |
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature list regen (`FEATURES-START/END`) via `patch_banner.py` — add LadrunoJ2Finite line | [#97](https://github.com/nmorabowen/OpenSees/pull/97) |
