---
wp: PR-381
title: "381 -- 2 vanilla row(s)"
pr: "#381"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`", "`SRC/interpreter/{OpenSeesOutputCommands.cpp,OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}`"]
table: "main"
legacy_seq: [59, 60]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno` ADR-41 C4: the `-mortar` contact parse gains `-tie` (a PERMANENT mesh-tie bond — the zero-gap limit) + `-epsTie auto\|<v>` (an alias for the `-epsN` penalty slot; a tie has one penalty) → threaded into `addMortarContact(..., isTie)`. `-tie` is REFUSED with `-mu/-cohesion/-tauMax` (a tie has no friction cone) and REQUIRES `-mortar` (validated after the option loop). | [#381](https://github.com/nmorabowen/OpenSees/pull/381) |
| `SRC/interpreter/{OpenSeesOutputCommands.cpp,OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}` | `// Ladruno` ADR-41 C4: add the `ladrunoMortarTieResidual` query → `OPS_LadrunoMortarTieResidual` returns the max tie bond `‖r̄‖_∞` over tie slave nodes (`LadrunoContactDomain::getMaxMortarTieResidual`; 0 if no engine/tie) — the convergence measure `analyze_augmented(query=ops.ladrunoMortarTieResidual)` reads to stop the tie Uzawa. OPS_ body + decl + Py (`Py_ops_LadrunoMortarTieResidual`) + Tcl (`Tcl_ops_LadrunoMortarTieResidual`) wrappers + dual `addCommand`. Additive; no existing command touched. | [#381](https://github.com/nmorabowen/OpenSees/pull/381) |
