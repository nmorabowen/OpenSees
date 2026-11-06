---
wp: PR-561
title: "561 -- 7 vanilla row(s)"
pr: "#561"
files: ["`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/interpreter/OpenSeesElementCommands.cpp`", "`SRC/element/TclElementCommands.cpp`", "`SRC/recorder/VTK_Recorder.cpp`", "`SRC/recorder/VTKHDF_Recorder.cpp`", "`SRC/recorder/PVDRecorder.cpp`", "`SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}`"]
table: "main"
legacy_seq: [279, 280, 281, 282, 283, 284, 285]
---
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno` ADR-72 P1: `#include "ladrunoBrick/LadrunoBrick20.h"` + `getNewElement` case `ELE_TAG_LadrunoBrick20` → `new LadrunoBrick20()` so DB/MPI `recvSelf` can reconstruct the 20-node hex | [#561](https://github.com/nmorabowen/OpenSees/pull/561) |
| `SRC/interpreter/OpenSeesElementCommands.cpp` | `// Ladruno` ADR-72 P1: `element` dispatch for `LadrunoBrick20`/`ladrunoBrick20` (fwd-decl + 2 `functionMap` inserts), serving both Tcl and openseespy | [#561](https://github.com/nmorabowen/OpenSees/pull/561) |
| `SRC/element/TclElementCommands.cpp` | `// Ladruno` ADR-72 P1: classic-Tcl dispatch row `{"LadrunoBrick20","ladrunoBrick20",OPS_LadrunoBrick20}` in the `ladrunoElementTable` + extern fwd-decl | [#561](https://github.com/nmorabowen/OpenSees/pull/561) |
| `SRC/recorder/VTK_Recorder.cpp` | `// Ladruno` ADR-72 P1: `vtktypes[ELE_TAG_LadrunoBrick20] = VTK_QUADRATIC_HEXAHEDRON` (one-liner beside `Twenty_Node_Brick`) | [#561](https://github.com/nmorabowen/OpenSees/pull/561) |
| `SRC/recorder/VTKHDF_Recorder.cpp` | `// Ladruno` ADR-72 P1: `vtktypes[ELE_TAG_LadrunoBrick20] = VTK_QUADRATIC_HEXAHEDRON` | [#561](https://github.com/nmorabowen/OpenSees/pull/561) |
| `SRC/recorder/PVDRecorder.cpp` | `// Ladruno` ADR-72 P1: `vtktypes[ELE_TAG_LadrunoBrick20] = VTK_QUADRATIC_HEXAHEDRON` | [#561](https://github.com/nmorabowen/OpenSees/pull/561) |
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature list regen (`FEATURES-START/END`) via `patch_banner.py` — add LadrunoBrick20 line | [#561](https://github.com/nmorabowen/OpenSees/pull/561) |
