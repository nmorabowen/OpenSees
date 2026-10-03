---
wp: PR-340
title: "340 -- 1 vanilla row(s)"
pr: "#340"
files: ["`SRC/tcl/commands.cpp`"]
table: "main"
legacy_seq: [29]
---
| `SRC/tcl/commands.cpp` | `// Ladruno` (SMS): register the 6 selective-mass-scaling integrators in the legacy `specifyIntegrator()` Tcl parser (extern decls + `else if` branches for `CentralDifferenceSMS`/`...SMSConsistent`, `ExplicitBatheSMS`/`...SMSConsistent`, `ExplicitBatheLNVDSMS`/`...SMSConsistent`, each `→ OPS_*()` with a null-guard before `setIntegrator`). They were openseespy-only (`OpenSeesCommands.cpp`) despite the Tcl splash banner advertising them → `integrator CentralDifferenceSMS …` errored "No Integrator type exists" in classic `OpenSees.exe`/`OpenSeesMP.exe`. Mirrors the existing `CentralDifferenceLadruno` Tcl wiring; additive. | [#340](https://github.com/nmorabowen/OpenSees/pull/340) |
