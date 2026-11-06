---
wp: PR-221
title: "221 -- 5 vanilla row(s)"
pr: "#221"
files: ["`SRC/classTags.h`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/interpreter/OpenSeesElementCommands.cpp`", "`SRC/element/CMakeLists.txt`", "`SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}`"]
table: "main"
legacy_seq: [168, 169, 170, 171, 172]
---
| `SRC/classTags.h` | `// Ladruno` (ADR 29): register `ELE_TAG_LadrunoKinematicCoupling`=33012 (RBE2 / kinematic coupling), after RBE3=33011 | [#221](https://github.com/nmorabowen/OpenSees/pull/221) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno`: `#include "ladrunoKinematicCoupling/LadrunoKinematicCoupling.h"` + `case ELE_TAG_LadrunoKinematicCoupling: return new LadrunoKinematicCoupling();` for parallel/database `recvSelf` | [#221](https://github.com/nmorabowen/OpenSees/pull/221) |
| `SRC/interpreter/OpenSeesElementCommands.cpp` | `// Ladruno`: register `element` dispatch for `LadrunoKinematicCoupling`/`ladrunoKinematicCoupling` (fwd-decl + 2 `functionMap` inserts) | [#221](https://github.com/nmorabowen/OpenSees/pull/221) |
| `SRC/element/CMakeLists.txt` | `// Ladruno`: `add_subdirectory(ladrunoKinematicCoupling)` (classTag 33012) | [#221](https://github.com/nmorabowen/OpenSees/pull/221) |
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature list regen (`FEATURES-START/END`) via `patch_banner.py` — add LadrunoKinematicCoupling line | [#221](https://github.com/nmorabowen/OpenSees/pull/221) |
