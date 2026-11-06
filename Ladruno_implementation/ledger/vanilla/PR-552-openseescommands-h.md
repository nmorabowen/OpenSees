---
wp: PR-552
title: "552 -- 6 vanilla row(s)"
pr: "#552"
files: ["`SRC/interpreter/OpenSeesCommands.h`", "`SRC/interpreter/PythonWrapper.cpp`", "`SRC/interpreter/TclWrapper.cpp`", "`SRC/tcl/commands.cpp`", "`SRC/tcl/commands.h`", "`SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}`"]
table: "main"
legacy_seq: [291, 292, 293, 294, 295, 296]
---
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno` ADR44 P3: declare `int OPS_LadrunoRandomResponse()` (the `randomResponse` command entry; body lives in `LadrunoModalResponse.cpp`). Additive. | [#552](https://github.com/nmorabowen/OpenSees/pull/552) |
| `SRC/interpreter/PythonWrapper.cpp` | `// Ladruno` ADR44 P3: `Py_ops_randomResponse` wrapper + `addCommand("randomResponse", ...)` (openseespy registration, mirroring the P2 precedent). | [#552](https://github.com/nmorabowen/OpenSees/pull/552) |
| `SRC/interpreter/TclWrapper.cpp` | `// Ladruno` ADR44 P3: `Tcl_ops_randomResponse` wrapper + `addCommand(interp, "randomResponse", ...)` (dual-wired w/ Python). | [#552](https://github.com/nmorabowen/OpenSees/pull/552) |
| `SRC/tcl/commands.cpp` | `// Ladruno` ADR44 P3: classic-Tcl bridge for `randomResponse` — `extern OPS_LadrunoRandomResponse` + handler + `Tcl_CreateCommand`, mirroring the #546 pattern (scalar/stats return via `OPS_SetDoubleOutput`→`Tcl_SetObjResult`). Additive. | [#552](https://github.com/nmorabowen/OpenSees/pull/552) |
| `SRC/tcl/commands.h` | `// Ladruno` ADR44 P3: forward declaration for the `randomResponse` handler. Additive. | [#552](https://github.com/nmorabowen/OpenSees/pull/552) |
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature regen via `patch_banner.py` — extend the `modal frequency domain` line with `randomResponse` (ADR44 P3). Banner strings only. | [#552](https://github.com/nmorabowen/OpenSees/pull/552) |
