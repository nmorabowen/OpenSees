---
wp: PR-718
title: "718 -- upstreamable-table row(s)"
pr: "#718"
files: ["`SRC/interpreter/{OpenSeesOutputCommands.cpp,OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}`", "`SRC/tcl/commands.cpp`"]
table: "upstreamable"
legacy_seq: [414, 419]
---
| `SRC/interpreter/{OpenSeesOutputCommands.cpp,OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}` | `// Ladruno` build-stamp query: `ladrunoBuild` command — `OPS_LadrunoBuild()` returns the `OPENSEES_VERSION` compile define (the CMake-stamped git hash the banner prints) via `OPS_SetString`; decl + Py + Tcl wrappers. Machine-readable engine provenance after the TIMs T1 wrong-build incident. | [#718](https://github.com/nmorabowen/OpenSees/pull/718) |
| `SRC/tcl/commands.cpp` | `// Ladruno` build-stamp query: classic-Tcl `ladrunoBuild` command (the `version` pattern — `Tcl_SetResult` of `OPENSEES_VERSION`) so `OpenSees.exe` answers without banner scraping. | [#718](https://github.com/nmorabowen/OpenSees/pull/718) |
