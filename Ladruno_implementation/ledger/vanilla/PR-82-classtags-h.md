---
wp: PR-82
title: "82 -- 5 vanilla row(s)"
pr: "#82"
files: ["`SRC/classTags.h`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/interpreter/OpenSeesNDMaterialCommands.cpp`", "`SRC/material/nD/CMakeLists.txt`", "`SRC/interpreter/PythonModule.cpp`"]
table: "main"
legacy_seq: [89, 90, 91, 92, 93]
---
| `SRC/classTags.h` | Register `ND_TAG_LadrunoJ2`=33011 (combined iso + Chaboche AF kinematic J2) | [#82](https://github.com/nmorabowen/OpenSees/pull/82) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno`: `#include "LadrunoJ2.h"` + `case ND_TAG_LadrunoJ2: return new LadrunoJ2();` so parallel/database `recvSelf` can reconstruct it | [#82](https://github.com/nmorabowen/OpenSees/pull/82) |
| `SRC/interpreter/OpenSeesNDMaterialCommands.cpp` | `// Ladruno`: fwd-decl `OPS_LadrunoJ2()` + `nDMaterialsMap["LadrunoJ2"]` (shared OpenSeesPy/openseesmp/interpreter-Tcl registry) | [#82](https://github.com/nmorabowen/OpenSees/pull/82) |
| `SRC/material/nD/CMakeLists.txt` | Add `LadrunoJ2.cpp` to `OPS_Material` sources + `LadrunoJ2.h` to headers | [#82](https://github.com/nmorabowen/OpenSees/pull/82) |
| `SRC/interpreter/PythonModule.cpp` | Splash-banner feature list regen (`FEATURES-START/END`) via `patch_banner.py` — add LadrunoJ2 line | [#82](https://github.com/nmorabowen/OpenSees/pull/82) |
