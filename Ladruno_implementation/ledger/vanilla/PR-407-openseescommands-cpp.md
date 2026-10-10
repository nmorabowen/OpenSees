---
wp: PR-407
title: "407 -- 2 vanilla row(s)"
pr: "#407"
files: ["`SRC/interpreter/{OpenSeesCommands.cpp,OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}`", "`SRC/interpreter/{OpenSeesMiscCommands.cpp,OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}`"]
table: "main"
legacy_seq: [186, 187]
---
| `SRC/interpreter/{OpenSeesCommands.cpp,OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}` | `// Ladruno` ADR-52 W1-I1b: add the read-only `ladrunoTrialResidualNorm <loadTime>` query → `OPS_LadrunoTrialResidualNorm` (`Domain::update()` to refresh element internal state to the node TRIAL disp, then the active integrator's `formUnbalance()`; returns the inf-norm of the free-DOF dynamic unbalance; optional `loadTime` re-applies loads at the step midpoint via `applyLoad`). Supplies the half-increment-residual primitive OpenSeesPy lacks (no commit; committed state untouched). OPS_ body + decl + Py/Tcl wrappers + dual `addCommand`. Additive; no existing command touched. | [#407](https://github.com/nmorabowen/OpenSees/pull/407) |
| `SRC/interpreter/{OpenSeesMiscCommands.cpp,OpenSeesCommands.h,PythonWrapper.cpp,TclWrapper.cpp}` | `// Ladruno` ADR-52 W1-I1b: add the `ladrunoSetNodeTrial nodeTag <ndof disp><ndof vel><ndof accel>` full-vector trial setter → `OPS_LadrunoSetNodeTrial` (sets the complete trial disp/vel/accel of a node in one call; the per-dof `setNodeDisp/Vel/Accel` each restart from the COMMITTED vector so repeated calls cannot build a multi-dof trial state). No commit. OPS_ body + decl + Py/Tcl wrappers + dual `addCommand`. Additive; no existing command touched. | [#407](https://github.com/nmorabowen/OpenSees/pull/407) |
