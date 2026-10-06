---
wp: WP-171
title: "A refused integrator raised nothing in openseespy and segfaulted classic Tcl"
date: 2026-10-06
---
### A refused `integrator` raised nothing in openseespy (previous one silently kept) and SEGFAULTED classic Tcl when an analysis existed (WP-171, 2026-10-06)
- **openseespy:** `OPS_Integrator()` returned `0` whatever the factory returned; `Py_ops_integrator` raises only on `< 0`. A refused integrator printed its WARNING and the run carried on with the previous integrator. ADR-76 had closed this hole for `OPS_Algorithm` only. Test code written before WP-171 that "checks a refusal" must have asserted on stderr, never on an exception (e.g. `tests/test_wp170_typed_arg_peeks.py` `_integrator_stderr`).
- **classic Tcl:** almost every `specifyIntegrator` branch is `theXIntegrator = OPS_X(); if (theXAnalysis != 0) theXAnalysis->setIntegrator(*theXIntegrator);` with no null check. With an analysis present, a refused integrator dereferenced null (exit 0xC0000005 / 139). Without one, it returned `TCL_OK` with the global set to null, so the next `analysis` printed `no Integrator specified, ... default will be used`. A parse error in a Ladruno factory could therefore change the integrator to a default with only that one line to show for it.
- **General shape — check it for every factory command:** the factory may return null; the dispatcher must (a) fail with an error code the interpreter turns into an exception or `TCL_ERROR`, and (b) leave the previous object in place. `OPS_Algorithm` (ADR-76) and `OPS_Integrator`/`specifyIntegrator` (WP-171) now do. Not audited: `OPS_System`, `OPS_Numberer`, `OPS_ConstraintHandler`, `OPS_Test` and their classic-Tcl twins.
- **Fix + gate:** see `ledger/vanilla/WP-171-*`; `tests/test_wp171_refused_integrator.py` (openseespy static + transient, with/without analysis, unknown type; classic-Tcl decks via `dist/bin/OpenSees.exe`, one of which segfaulted pre-fix).
