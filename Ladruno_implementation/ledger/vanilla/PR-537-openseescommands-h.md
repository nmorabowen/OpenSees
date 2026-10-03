---
wp: PR-537
title: "537 -- 4 vanilla row(s)"
pr: "#537"
files: ["`SRC/interpreter/OpenSeesCommands.h`", "`SRC/interpreter/PythonWrapper.cpp`", "`SRC/interpreter/TclWrapper.cpp`", "`SRC/analysis/analysis/CMakeLists.txt`"]
table: "main"
legacy_seq: [247, 248, 249, 250]
---
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno` ADR44 P1a: declare `int OPS_LadrunoModalResponseHistory()` (the `modalResponseHistory` command entry; the body lives in `LadrunoModalResponse.cpp`, not `OpenSeesCommands.cpp`). Additive. | [#537](https://github.com/nmorabowen/OpenSees/pull/537) |
| `SRC/interpreter/PythonWrapper.cpp` | `// Ladruno` ADR44 P1a: `Py_ops_modalResponseHistory` wrapper + `addCommand("modalResponseHistory", ...)` (openseespy registration, mirroring the complexEigen precedent). | [#537](https://github.com/nmorabowen/OpenSees/pull/537) |
| `SRC/interpreter/TclWrapper.cpp` | `// Ladruno` ADR44 P1a: `Tcl_ops_modalResponseHistory` wrapper + `addCommand(interp, "modalResponseHistory", ...)` (Tcl registration; dual-wired w/ Python). | [#537](https://github.com/nmorabowen/OpenSees/pull/537) |
| `SRC/analysis/analysis/CMakeLists.txt` | ADR44 P1a: add `LadrunoModalResponse.{cpp,h}` to `OPS_Analysis` target_sources (build wiring only). | [#537](https://github.com/nmorabowen/OpenSees/pull/537) |
