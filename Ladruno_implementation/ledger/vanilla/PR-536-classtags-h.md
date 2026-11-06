---
wp: PR-536
title: "536 -- 4 vanilla row(s)"
pr: "#536"
files: ["`SRC/classTags.h`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/interpreter/OpenSeesNDMaterialCommands.cpp`", "`SRC/material/nD/CMakeLists.txt`"]
table: "main"
legacy_seq: [104, 110, 111, 112]
---
| `SRC/classTags.h` | `// Ladruno` (ADR-25 P5): register `ND_TAG_LogStrain2D`=33016 (2D plane Hencky finite-strain adaptor); replaces the prior "33016 reserved" placeholder comment | [#536](https://github.com/nmorabowen/OpenSees/pull/536) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno` (ADR-25 P5): `#include "LogStrain2D.h"` + `case ND_TAG_LogStrain2D: return new LogStrain2D();` so parallel/database `recvSelf` can reconstruct it | [#536](https://github.com/nmorabowen/OpenSees/pull/536) |
| `SRC/interpreter/OpenSeesNDMaterialCommands.cpp` | `// Ladruno` (ADR-25 P5): fwd-decl `OPS_LogStrain2D()` + `nDMaterialsMap["LogStrain2D"]` (shared OpenSeesPy/openseesmp/interpreter-Tcl registry) | [#536](https://github.com/nmorabowen/OpenSees/pull/536) |
| `SRC/material/nD/CMakeLists.txt` | `// Ladruno` (ADR-25 P5): add `LogStrain2D.cpp` to `OPS_Material` sources + `LogStrain2D.h` / `FiniteStrainND2DMaterial.h` to headers | [#536](https://github.com/nmorabowen/OpenSees/pull/536) |
