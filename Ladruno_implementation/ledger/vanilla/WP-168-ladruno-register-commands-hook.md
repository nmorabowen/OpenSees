---
wp: WP-168
pr: "#922"
title: "WP-168 vanilla rows"
date: 2026-10-05
files: ["`SRC/tcl/commands.cpp`", "`SRC/interpreter/TclWrapper.cpp`", "`SRC/interpreter/PythonWrapper.cpp`"]
table: "main"
---
| `SRC/tcl/commands.cpp` | `// Ladruno WP-168`: the 23 hand-written fork `Tcl_CreateCommand` registrations (contact family, ADR-44 modal family, `ladrunoArcLength`/`ladrunoDR`, `ladrunoBuild`/`ladrunoThreads`/`ladrunoMutation`/`ladrunoSANISANDReplay`, `profiler`, `ladrunoNumbering`, `LadrunoStaggeredAnalyze`) replaced by `#include "LadrunoCommandsClassicTcl.h"` before `OpenSeesAppInit` and one `Ladruno_registerCommands(interp)` call in it. The classic bridges themselves are unchanged. | [#922](https://github.com/nmorabowen/OpenSees/pull/922) |
| `SRC/interpreter/TclWrapper.cpp` | `// Ladruno WP-168`: the 29 fork `addCommand(interp, …)` registrations and their 29 hand-written `Tcl_ops_*` bridges removed (net fork footprint shrinks); `#include "LadrunoCommandsTclWrapper.h"` + one `Ladruno_registerCommands(this, interp)` call, which generates the same bridges from the table. | [#922](https://github.com/nmorabowen/OpenSees/pull/922) |
| `SRC/interpreter/PythonWrapper.cpp` | `// Ladruno WP-168`: the 29 fork `addCommand(…)` registrations and their 29 hand-written `Py_ops_*` bridges removed; `#include "LadrunoCommandsPython.h"` + one `Ladruno_registerCommands(this)` call, which generates the same bridges from the table. | [#922](https://github.com/nmorabowen/OpenSees/pull/922) |
