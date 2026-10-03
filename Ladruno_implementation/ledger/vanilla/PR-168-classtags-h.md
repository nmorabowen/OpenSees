---
wp: PR-168
title: "168 -- 5 vanilla row(s)"
pr: "#168"
files: ["`SRC/classTags.h`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/interpreter/OpenSeesUniaxialMaterialCommands.cpp`", "`SRC/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp`", "`SRC/material/uniaxial/CMakeLists.txt`"]
table: "main"
legacy_seq: [207, 208, 209, 210, 211]
---
| `SRC/classTags.h` | `// Ladruno`: register `MAT_TAG_LadrunoBondSlip`=33002 (uniaxial band) | [#168](https://github.com/nmorabowen/OpenSees/pull/168) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno`: `#include "LadrunoBondSlip.h"` + `case MAT_TAG_LadrunoBondSlip: return new LadrunoBondSlip();` | [#168](https://github.com/nmorabowen/OpenSees/pull/168) |
| `SRC/interpreter/OpenSeesUniaxialMaterialCommands.cpp` | `// Ladruno`: fwd-decl `OPS_LadrunoBondSlip()` + `uniaxialMaterialsMap["LadrunoBondSlip"]` | [#168](https://github.com/nmorabowen/OpenSees/pull/168) |
| `SRC/material/uniaxial/TclModelBuilderUniaxialMaterialCommand.cpp` | `// Ladruno`: extern `OPS_LadrunoBondSlip()` + `strcmp(argv[1],"LadrunoBondSlip")` classic-Tcl dispatch | [#168](https://github.com/nmorabowen/OpenSees/pull/168) |
| `SRC/material/uniaxial/CMakeLists.txt` | `// Ladruno`: add `LadrunoBondSlip.cpp`/`.h` to `OPS_Material` | [#168](https://github.com/nmorabowen/OpenSees/pull/168) |
