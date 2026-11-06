---
wp: PR-176
title: "176 -- 4 vanilla row(s)"
pr: "#176, #178"
files: ["`SRC/interpreter/OpenSeesCommands.cpp`", "`SRC/interpreter/OpenSeesCommands.h`", "`SRC/interpreter/PythonWrapper.cpp`", "`SRC/interpreter/TclWrapper.cpp`"]
table: "main"
legacy_seq: [136, 137, 138, 139]
---
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno`: (a) backfill the `integrator LadrunoArcLength` dispatch branch → `OPS_LadrunoArcLength()` (originally added in #157, un-ledgered); (b) new `OPS_LadrunoArcLengthCmd()` runtime command (Layer-B, ADR §4.3) — reads `cmds->getStaticIntegrator()`, guards on `INTEGRATOR_TAGS_LadrunoArcLength`, dispatches `reduceStep`/`increaseStep`/`setArcLength` mutators + `arcLength`/`deltaLambdaStep`/`currentLambda`/`sign`/`deltaUstepNorm` queries; +`#include <LadrunoArcLength.h>`, `<classTags.h>`. Additive. | [#176](https://github.com/nmorabowen/OpenSees/pull/178) |
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno`: declare `int OPS_LadrunoArcLengthCmd()` | [#176](https://github.com/nmorabowen/OpenSees/pull/178) |
| `SRC/interpreter/PythonWrapper.cpp` | `// Ladruno`: register the `ladrunoArcLength` openseespy command (`Py_ops_ladrunoArcLength` → `OPS_LadrunoArcLengthCmd`) — Layer-B script-driven cut-and-retry exposure | [#176](https://github.com/nmorabowen/OpenSees/pull/178) |
| `SRC/interpreter/TclWrapper.cpp` | `// Ladruno`: register the `ladrunoArcLength` Tcl command (modern interpreter path, `Tcl_ops_ladrunoArcLength` → `OPS_LadrunoArcLengthCmd`) | [#176](https://github.com/nmorabowen/OpenSees/pull/178) |
