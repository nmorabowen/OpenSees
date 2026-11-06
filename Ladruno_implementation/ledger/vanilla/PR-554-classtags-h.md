---
wp: PR-554
title: "554 -- 4 vanilla row(s)"
pr: "#554"
files: ["`SRC/classTags.h`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/interpreter/OpenSeesElementCommands.cpp`", "`SRC/element/TclElementCommands.cpp`"]
table: "main"
legacy_seq: [106, 107, 108, 109]
---
| `SRC/classTags.h` | `// Ladruno` (ADR-70 P3): register `ELE_TAG_LadrunoLST`=33016 (6-node linear-strain triangle, ELE registry — distinct from `ND_TAG_LogStrain2D`=33016); replaces the P0 "33016 RESERVED" placeholder comment | [#554](https://github.com/nmorabowen/OpenSees/pull/554) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno` (ADR-70 P3): `#include "ladrunoPlane/LadrunoLST.h"` + `case ELE_TAG_LadrunoLST: return new LadrunoLST();` so parallel/database `recvSelf` can reconstruct it | [#554](https://github.com/nmorabowen/OpenSees/pull/554) |
| `SRC/interpreter/OpenSeesElementCommands.cpp` | `// Ladruno` (ADR-70 P3): fwd-decl `OPS_LadrunoLST()` + `functionMap` entries `LadrunoLST`/`ladrunoLST` | [#554](https://github.com/nmorabowen/OpenSees/pull/554) |
| `SRC/element/TclElementCommands.cpp` | `// Ladruno` (ADR-70 P3): extern `OPS_LadrunoLST` + classic-Tcl dispatch row `LadrunoLST`/`ladrunoLST` | [#554](https://github.com/nmorabowen/OpenSees/pull/554) |
