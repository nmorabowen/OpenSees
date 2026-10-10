---
wp: PR-375
title: "375 -- 1 vanilla row(s)"
pr: "#375"
files: ["`SRC/interpreter/{OpenSeesOutputCommands.cpp,OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}`"]
table: "main"
legacy_seq: [58]
---
| `SRC/interpreter/{OpenSeesOutputCommands.cpp,OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}` | `// Ladruno` ADR-41 C2.2: add the `ladrunoMortarPenetration` query → `OPS_LadrunoMortarPenetration` returns the max KKT-active mortar penetration `‖ḡ‖_∞` (`LadrunoContactDomain::getMaxMortarPenetration`; 0 if no engine/contact) — the convergence measure the held-load `analyze_augmented` proc reads to stop Uzawa augmenting. OPS_ body + decl + Py (`Py_ops_LadrunoMortarPenetration`) + Tcl (`Tcl_ops_LadrunoMortarPenetration`) wrappers + dual `addCommand`. Additive; no existing command touched. | [#375](https://github.com/nmorabowen/OpenSees/pull/375) |
