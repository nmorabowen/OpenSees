---
wp: PR-558
title: "558 -- upstreamable-table row(s)"
pr: "#558"
files: ["`SRC/classTags.h`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/interpreter/OpenSeesElementCommands.cpp`", "`SRC/element/TclElementCommands.cpp`", "`SRC/tcl/tclMain.cpp` + `SRC/interpreter/PythonModule.cpp`"]
table: "upstreamable"
legacy_seq: [317, 318, 319, 320, 321]
---
| `SRC/classTags.h` | `// Ladruno` (ADR-70 P4a): register `ELE_TAG_LadrunoCSTPair`=33021 (disjoint 2-triangle F-bar-Patch macro-element, dSNPO 15.1.9) after the SolidShell 33020 row | [#558](https://github.com/nmorabowen/OpenSees/pull/558) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno` (ADR-70 P4a): include + `case ELE_TAG_LadrunoCSTPair` broker constructor | [#558](https://github.com/nmorabowen/OpenSees/pull/558) |
| `SRC/interpreter/OpenSeesElementCommands.cpp` | `// Ladruno` (ADR-70 P4a): fwd-decl + `functionMap` entries `LadrunoCSTPair`/`ladrunoCSTPair` | [#558](https://github.com/nmorabowen/OpenSees/pull/558) |
| `SRC/element/TclElementCommands.cpp` | `// Ladruno` (ADR-70 P4a): extern + classic-Tcl dispatch row `LadrunoCSTPair`/`ladrunoCSTPair` | [#558](https://github.com/nmorabowen/OpenSees/pull/558) |
| `SRC/tcl/tclMain.cpp` + `SRC/interpreter/PythonModule.cpp` | banner FEATURES block regen (`patch_banner.py`) — plane-family line now lists `LadrunoCSTPair` (F-bar-Patch macro) | [#558](https://github.com/nmorabowen/OpenSees/pull/558) |
