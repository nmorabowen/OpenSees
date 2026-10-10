---
wp: PR-454
title: "454 -- 3 vanilla row(s)"
pr: "#454"
files: ["`SRC/domain/constraints/CMakeLists.txt`", "`SRC/interpreter/{OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}`", "`SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}`"]
table: "main"
legacy_seq: [212, 213, 214]
---
| `SRC/domain/constraints/CMakeLists.txt` | `// Ladruno` ADR-62: add `LadrunoTie.{cpp,h}` to the `OPS_Domain` sources (the kinematic mesh-tie generator) | [#454](https://github.com/nmorabowen/OpenSees/pull/454) |
| `SRC/interpreter/{OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}` | `// Ladruno` ADR-62: add the `LadrunoTie` command — `OPS_LadrunoTie()` decl + `Py_ops_LadrunoTie`/`Tcl_ops_LadrunoTie` wrappers + dual `addCommand` (mirrors the `equationConstraint` registration). Additive; no existing command touched. | [#454](https://github.com/nmorabowen/OpenSees/pull/454) |
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature regen via `patch_banner.py` — add the `LadrunoTie` line | [#454](https://github.com/nmorabowen/OpenSees/pull/454) |
