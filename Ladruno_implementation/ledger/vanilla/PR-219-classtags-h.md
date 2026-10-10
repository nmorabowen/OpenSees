---
wp: PR-219
title: "219 -- 5 vanilla row(s)"
pr: "#219"
files: ["`SRC/classTags.h`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/interpreter/OpenSeesElementCommands.cpp`", "`SRC/element/CMakeLists.txt`", "`SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}`"]
table: "main"
legacy_seq: [163, 164, 165, 166, 167]
---
| `SRC/classTags.h` | `// Ladruno` (ADR 28): register `ELE_TAG_LadrunoDistributingCoupling`=33011 (RBE3 / distributing coupling); 33009/33010 noted reserved for the VEM/SBFEM frontier elements | [#219](https://github.com/nmorabowen/OpenSees/pull/219) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno`: `#include "ladrunoDistributingCoupling/LadrunoDistributingCoupling.h"` + `case ELE_TAG_LadrunoDistributingCoupling: return new LadrunoDistributingCoupling();` for parallel/database `recvSelf` | [#219](https://github.com/nmorabowen/OpenSees/pull/219) |
| `SRC/interpreter/OpenSeesElementCommands.cpp` | `// Ladruno`: register `element` dispatch for `LadrunoDistributingCoupling`/`ladrunoDistributingCoupling` (fwd-decl + 2 `functionMap` inserts) | [#219](https://github.com/nmorabowen/OpenSees/pull/219) |
| `SRC/element/CMakeLists.txt` | `// Ladruno`: `add_subdirectory(ladrunoDistributingCoupling)` (classTag 33011) | [#219](https://github.com/nmorabowen/OpenSees/pull/219) |
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature list regen (`FEATURES-START/END`) via `patch_banner.py` — add LadrunoDistributingCoupling line | [#219](https://github.com/nmorabowen/OpenSees/pull/219) |
