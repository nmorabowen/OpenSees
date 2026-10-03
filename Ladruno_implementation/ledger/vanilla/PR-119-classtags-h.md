---
wp: PR-119
title: "119 -- 5 vanilla row(s)"
pr: "#119"
files: ["`SRC/classTags.h`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/interpreter/OpenSeesUniaxialMaterialCommands.cpp`", "`SRC/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp`", "`SRC/material/uniaxial/CMakeLists.txt`"]
table: "main"
legacy_seq: [130, 131, 132, 133, 134]
---
| `SRC/classTags.h` | Register `MAT_TAG_LadrunoRebarBuckling`=33001 (rebar-buckling wrapper, Dhakal-Maekawa; sibling of LadrunoUniaxialJ2=33000 in the uniaxial band) | [#119](https://github.com/nmorabowen/OpenSees/pull/119) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno`: `#include "LadrunoRebarBuckling.h"` + `case MAT_TAG_LadrunoRebarBuckling: return new LadrunoRebarBuckling();` (nested-material `recvSelf` reconstruction) | [#119](https://github.com/nmorabowen/OpenSees/pull/119) |
| `SRC/interpreter/OpenSeesUniaxialMaterialCommands.cpp` | `// Ladruno`: fwd-decl `OPS_LadrunoRebarBuckling()` + `uniaxialMaterialsMap["LadrunoRebarBuckling"]` (shared OpenSeesPy/openseesmp/interpreter-Tcl registry) | [#119](https://github.com/nmorabowen/OpenSees/pull/119) |
| `SRC/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp` | `// Ladruno`: extern `OPS_LadrunoRebarBuckling()` + `strcmp(argv[1],"LadrunoRebarBuckling")` dispatch block (classic-Tcl `OpenSees.exe` path) | [#119](https://github.com/nmorabowen/OpenSees/pull/119) |
| `SRC/material/uniaxial/CMakeLists.txt` | Add `LadrunoRebarBuckling.cpp` to `OPS_Material` sources + `LadrunoRebarBuckling.h` to headers | [#119](https://github.com/nmorabowen/OpenSees/pull/119) |
