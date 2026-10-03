---
wp: PR-264
title: "264 -- 6 vanilla row(s)"
pr: "#264"
files: ["`SRC/classTags.h`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/interpreter/OpenSeesUniaxialMaterialCommands.cpp`", "`SRC/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp`", "`SRC/material/uniaxial/CMakeLists.txt`", "`SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}`"]
table: "main"
legacy_seq: [174, 175, 176, 177, 178, 179]
---
| `SRC/classTags.h` | Register `MAT_TAG_LadrunoCohesiveHinge`=33003 (discrete cohesive moment-rotation hinge for the LadrunoDispBeamColumn Tier-2 embedded strong-discontinuity hinge; uniaxial band, after LadrunoBondSlip=33002) | [#264](https://github.com/nmorabowen/OpenSees/pull/264) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno`: `#include "LadrunoCohesiveHinge.h"` + `case MAT_TAG_LadrunoCohesiveHinge: return new LadrunoCohesiveHinge();` so parallel/database `recvSelf` can reconstruct it | [#264](https://github.com/nmorabowen/OpenSees/pull/264) |
| `SRC/interpreter/OpenSeesUniaxialMaterialCommands.cpp` | `// Ladruno`: fwd-decl `OPS_LadrunoCohesiveHinge()` + `uniaxialMaterialsMap["LadrunoCohesiveHinge"]` (OpenSeesPy/openseesmp/interpreter-Tcl registry) | [#264](https://github.com/nmorabowen/OpenSees/pull/264) |
| `SRC/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp` | `// Ladruno`: extern `OPS_LadrunoCohesiveHinge()` + `strcmp(argv[1],"LadrunoCohesiveHinge")` dispatch block (classic-Tcl `OpenSees.exe` path) | [#264](https://github.com/nmorabowen/OpenSees/pull/264) |
| `SRC/material/uniaxial/CMakeLists.txt` | Add `LadrunoCohesiveHinge.cpp` to `OPS_Material` sources + `LadrunoCohesiveHinge.h` to headers | [#264](https://github.com/nmorabowen/OpenSees/pull/264) |
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature list regen (`FEATURES-START/END`) via `patch_banner.py` — add LadrunoCohesiveHinge line | [#264](https://github.com/nmorabowen/OpenSees/pull/264) |
