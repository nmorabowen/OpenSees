---
wp: PR-506
title: "506 -- 4 vanilla row(s)"
pr: "#506"
files: ["`SRC/interpreter/OpenSeesCommands.{h,cpp}`", "`SRC/interpreter/PythonWrapper.cpp`", "`SRC/interpreter/TclWrapper.cpp`", "`SRC/analysis/analysis/CMakeLists.txt`"]
table: "main"
legacy_seq: [243, 244, 245, 246]
---
| `SRC/interpreter/OpenSeesCommands.{h,cpp}` | `// Ladruno` ADR46 P0: `OPS_complexEigen()` — complex/state-space modal command (P0: only the `-qz p Mt Ct Kt` reduced-pencil kernel entry for oracle verification, 7 doubles/mode incl. residual; the domain-coupled projection is P1) + declaration + `#include <LadrunoComplexEigen.h>`. | [#506](https://github.com/nmorabowen/OpenSees/pull/506) |
| `SRC/interpreter/PythonWrapper.cpp` | `// Ladruno` ADR46 P0: `Py_ops_complexEigen` wrapper + `addCommand("complexEigen", ...)` (openseespy registration). | [#506](https://github.com/nmorabowen/OpenSees/pull/506) |
| `SRC/interpreter/TclWrapper.cpp` | `// Ladruno` ADR46 P0: `Tcl_ops_complexEigen` wrapper + `addCommand(interp, "complexEigen", ...)` (Tcl registration). | [#506](https://github.com/nmorabowen/OpenSees/pull/506) |
| `SRC/analysis/analysis/CMakeLists.txt` | ADR46 P0: add `LadrunoComplexEigen.{cpp,h}` to `OPS_Analysis` target_sources (build wiring only, no upstream logic touched). | [#506](https://github.com/nmorabowen/OpenSees/pull/506) |
