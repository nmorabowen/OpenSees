---
wp: PR-482
title: "482 -- 6 vanilla row(s)"
pr: "#482"
files: ["`SRC/classTags.h`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/interpreter/OpenSeesElementCommands.cpp`", "`SRC/element/TclElementCommands.cpp`", "`SRC/element/CMakeLists.txt`", "`SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}`"]
table: "main"
legacy_seq: [218, 219, 220, 221, 223, 224]
---
| `SRC/classTags.h` | `// Ladruno` ADR-66: `ELE_TAG_LadrunoSolidShell`=33020 (8-node ANS/EAS solid-shell; honored as reserved by ADR 19 — 33016-33019 remain free ELE slots) | [#482](https://github.com/nmorabowen/OpenSees/pull/482) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno` ADR-66: `getNewElement` case `ELE_TAG_LadrunoSolidShell` → `new LadrunoSolidShell()` (+include) so DB/MPI `recvSelf` can reconstruct it | [#482](https://github.com/nmorabowen/OpenSees/pull/482) |
| `SRC/interpreter/OpenSeesElementCommands.cpp` | `// Ladruno` ADR-66: `element` dispatch for `LadrunoSolidShell`/`ladrunoSolidShell` (fwd-decl + `functionMap`), serving both Tcl and openseespy | [#482](https://github.com/nmorabowen/OpenSees/pull/482) |
| `SRC/element/TclElementCommands.cpp` | `// Ladruno` ADR-66: classic-Tcl element table entry `{"LadrunoSolidShell","ladrunoSolidShell",OPS_LadrunoSolidShell}` (+extern decl) | [#482](https://github.com/nmorabowen/OpenSees/pull/482) |
| `SRC/element/CMakeLists.txt` | `// Ladruno` ADR-66: `add_subdirectory(ladrunoSolidShell)` (classTag 33020) | [#482](https://github.com/nmorabowen/OpenSees/pull/482) |
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature regen via `patch_banner.py` — add the `LadrunoSolidShell` line (ADR-66 P5.1). Banner strings only; the element sources are fork-authored (see [[LEDGER_implementations]]). | [#482](https://github.com/nmorabowen/OpenSees/pull/482) |
