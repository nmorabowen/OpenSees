---
wp: PR-544
title: "544 -- 4 vanilla row(s)"
pr: "#544"
files: ["`SRC/interpreter/OpenSeesCommands.h`", "`SRC/interpreter/PythonWrapper.cpp`", "`SRC/interpreter/TclWrapper.cpp`", "`SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}`"]
table: "main"
legacy_seq: [253, 254, 255, 256]
---
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno` ADR44 P2: declare `int OPS_LadrunoFrequencyResponse()` + `int OPS_LadrunoSteadyStateDynamics()` (the `frequencyResponse` / `steadyStateDynamics` command entries; bodies live in `LadrunoModalResponse.cpp`). Additive. | [#544](https://github.com/nmorabowen/OpenSees/pull/544) |
| `SRC/interpreter/PythonWrapper.cpp` | `// Ladruno` ADR44 P2: `Py_ops_frequencyResponse` / `Py_ops_steadyStateDynamics` wrappers + `addCommand(...)` (openseespy registration, mirroring the P1a modalResponseHistory precedent). | [#544](https://github.com/nmorabowen/OpenSees/pull/544) |
| `SRC/interpreter/TclWrapper.cpp` | `// Ladruno` ADR44 P2: `Tcl_ops_frequencyResponse` / `Tcl_ops_steadyStateDynamics` wrappers + `addCommand(interp, ...)` (dual-wired w/ Python). | [#544](https://github.com/nmorabowen/OpenSees/pull/544) |
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature regen via `patch_banner.py` — add the `modal frequency domain — frequencyResponse / steadyStateDynamics` line (ADR44 P2). Banner strings only; `LadrunoModalResponse.{cpp,h}` are fork-authored. | [#544](https://github.com/nmorabowen/OpenSees/pull/544) |
