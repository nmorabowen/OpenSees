---
wp: WP-171
title: "Refused integrator is an error (openseespy + classic Tcl)"
date: 2026-10-06
files: ["`SRC/interpreter/OpenSeesCommands.cpp`", "`SRC/tcl/commands.cpp`"]
table: "upstreamable"
---
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno WP-171`: `OPS_Integrator()` — **a null factory result is now `return -1`** (one guard after the type ladder, covering every factory and the unknown-type branch). It used to fall through to `return 0`, so `Py_ops_integrator` raised nothing for a refused integrator (`ExplicitBathe 0.54 -lnvd 1.5`, `LoadControl` with no args, an unknown type) and the PREVIOUS integrator silently stayed in force. The ADR-76 fix for `OPS_Algorithm` in the same file, applied to integrators. Inert for every accepted integrator. | WP-171 |
| `SRC/tcl/commands.cpp` | `// Ladruno WP-171`: `specifyIntegrator` — **a refused integrator no longer segfaults or silently swaps.** The stock body is renamed `ladrunoSpecifyIntegratorImpl` (static) and unchanged except that its 66 `theXAnalysis->setIntegrator(*theXIntegrator)` calls (13 static, 53 transient; ~55 of them had no null check) go through the null-safe `ladrunoSetIntegratorIfAny`. Measured pre-fix on `OpenSees.exe`: with an analysis present, `integrator ExplicitBathe 0.54 -lnvd 1.5` **segfaulted** (`*null`); without one it returned `TCL_OK` and nulled the global, so `analysis Transient` fell back to the default Newmark. The new `specifyIntegrator` wrapper saves both global integrator pointers, calls the impl, and when no NEW integrator was created restores them and returns `TCL_ERROR` with `previous integrator left unchanged`. A branch that already failed correctly (e.g. stock `LoadControl` arg errors) keeps its own message plus that line. Inert for every accepted integrator. | WP-171 |
