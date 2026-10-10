---
wp: PR-184
title: "184 -- 7 vanilla row(s)"
pr: "#184"
files: ["`SRC/classTags.h`", "`SRC/analysis/integrator/CMakeLists.txt`", "`SRC/convergenceTest/CMakeLists.txt`", "`SRC/interpreter/OpenSeesCommands.{cpp,h}`", "`SRC/tcl/commands.cpp`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/runtime/runtime/TclPackageClassBroker.cpp`"]
table: "main"
legacy_seq: [153, 154, 155, 156, 157, 158, 159]
---
| `SRC/classTags.h` | `// Ladruno` (ADR-20 §8 #3/#4): register `INTEGRATOR_TAGS_LadrunoIndirectControl`=33006 (indirect/CMOD control) + `CONVERGENCE_TEST_LadrunoStabilizedUnbalance`=33000 (true-equilibrium NormUnbalance; per-registry band ≥33000) | [#184](https://github.com/nmorabowen/OpenSees/pull/184) |
| `SRC/analysis/integrator/CMakeLists.txt` | `// Ladruno`: add `LadrunoIndirectControl.cpp`/`.h` to `OPS_Analysis` integrator sources/headers | [#184](https://github.com/nmorabowen/OpenSees/pull/184) |
| `SRC/convergenceTest/CMakeLists.txt` | `// Ladruno`: add `LadrunoStabilizedUnbalance.cpp`/`.h` to `OPS_ConvergenceTest` sources/headers | [#184](https://github.com/nmorabowen/OpenSees/pull/184) |
| `SRC/interpreter/OpenSeesCommands.{cpp,h}` | `// Ladruno`: integrator dispatch `LadrunoIndirectControl` → `OPS_LadrunoIndirectControl()` + test dispatch `LadrunoStabilizedUnbalance` → `OPS_LadrunoStabilizedUnbalance()` (+ the two `OPS_*` decls in the header; shared OpenSeesPy/openseesmp/interpreter-Tcl registry) | [#184](https://github.com/nmorabowen/OpenSees/pull/184) |
| `SRC/tcl/commands.cpp` | `// Ladruno`: classic-Tcl `integrator LadrunoIndirectControl` branch (extern + `strcmp` dispatch → `OPS_LadrunoIndirectControl`) **and** the `test LadrunoStabilizedUnbalance` branch in `specifyCTest` (`#include <LadrunoStabilizedUnbalance.h>` + `new LadrunoStabilizedUnbalance(tol,numIter,printIt,normType,maxTol)`) so the convergence test is reachable from the classic `OpenSees` Tcl binary, not just openseespy (review finding) | [#184](https://github.com/nmorabowen/OpenSees/pull/184) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno`: `getNewStaticIntegrator` case `INTEGRATOR_TAGS_LadrunoIndirectControl` → `new LadrunoIndirectControl()` + `getNewConvergenceTest` case `CONVERGENCE_TEST_LadrunoStabilizedUnbalance` → `new LadrunoStabilizedUnbalance()` (+ includes), for parallel/database `recvSelf` | [#184](https://github.com/nmorabowen/OpenSees/pull/184) |
| `SRC/runtime/runtime/TclPackageClassBroker.cpp` | `// Ladruno`: same two broker cases (integrator + convergence test) + includes, for the Tcl-package class broker | [#184](https://github.com/nmorabowen/OpenSees/pull/184) |
