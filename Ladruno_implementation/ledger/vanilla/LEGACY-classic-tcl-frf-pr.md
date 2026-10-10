---
wp: LEGACY
title: "classic-Tcl FRF PR -- 2 vanilla row(s)"
files: ["`SRC/tcl/commands.cpp`", "`SRC/tcl/commands.h`"]
table: "main"
legacy_seq: [277, 278]
---
| `SRC/tcl/commands.cpp` | `// Ladruno` ADR-44: wire the modal-response family into classic-Tcl dispatch — `extern OPS_LadrunoModalResponseHistory/OPS_LadrunoFrequencyResponse/OPS_LadrunoSteadyStateDynamics` + three handler fns (`modalResponseHistory`/`frequencyResponse`/`steadyStateDynamics`, each `OPS_ResetInputNoBuilder`→`OPS_*`, mirroring `responseSpectrumAnalysis`) + three `Tcl_CreateCommand` registrations. The splash banner (`tclMain.cpp`) already advertised these but only `TclWrapper.cpp` (unbuilt `TclInterpreter`) wired them, so classic `.tcl` users hit "invalid command name". The `{f,…}` return list is delivered via `OPS_SetDoubleListsOutput`→`Tcl_SetObjResult` (works in classic mode). Additive. | classic-Tcl FRF PR |
| `SRC/tcl/commands.h` | `// Ladruno` ADR-44: forward declarations for the three modal-response handlers (`modalResponseHistory`/`frequencyResponse`/`steadyStateDynamics`) so the `Tcl_CreateCommand` block resolves them ahead of their definitions. | classic-Tcl FRF PR |
