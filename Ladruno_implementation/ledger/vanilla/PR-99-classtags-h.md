---
wp: PR-99
title: "99 -- 6 vanilla row(s)"
pr: "#99"
files: ["`SRC/classTags.h`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/interpreter/OpenSeesUniaxialMaterialCommands.cpp`", "`SRC/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp`", "`SRC/material/uniaxial/CMakeLists.txt`", "`SRC/material/nD/CMakeLists.txt`"]
table: "main"
legacy_seq: [124, 125, 126, 127, 128, 129]
---
| `SRC/classTags.h` | Register `MAT_TAG_LadrunoUniaxialJ2`=33000 (uniaxial combined iso + Chaboche AF kinematic J2; first Ladruno *uniaxial* tag, 33000-band is per-registry) | [#99](https://github.com/nmorabowen/OpenSees/pull/99) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno`: `#include "LadrunoUniaxialJ2.h"` + `case MAT_TAG_LadrunoUniaxialJ2: return new LadrunoUniaxialJ2();` so parallel/database `recvSelf` can reconstruct it | [#99](https://github.com/nmorabowen/OpenSees/pull/99) |
| `SRC/interpreter/OpenSeesUniaxialMaterialCommands.cpp` | `// Ladruno`: fwd-decl `OPS_LadrunoUniaxialJ2()` + `uniaxialMaterialsMap["LadrunoUniaxialJ2"]` (shared OpenSeesPy/openseesmp/interpreter-Tcl registry) | [#99](https://github.com/nmorabowen/OpenSees/pull/99) |
| `SRC/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp` | `// Ladruno`: extern `OPS_LadrunoUniaxialJ2()` + `strcmp(argv[1],"LadrunoUniaxialJ2")` dispatch block (classic-Tcl `OpenSees.exe` path) | [#99](https://github.com/nmorabowen/OpenSees/pull/99) |
| `SRC/material/uniaxial/CMakeLists.txt` | Add `LadrunoUniaxialJ2.cpp` to `OPS_Material` sources + `LadrunoUniaxialJ2.h` to headers | [#99](https://github.com/nmorabowen/OpenSees/pull/99) |
| `SRC/material/nD/CMakeLists.txt` | Add header-only `LadrunoHardening.h` (shared isotropic law, consumed by both `LadrunoJ2` and uniaxial `LadrunoUniaxialJ2` — the oracle contract) to `OPS_Material` PUBLIC headers | [#99](https://github.com/nmorabowen/OpenSees/pull/99) |
