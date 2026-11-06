---
wp: PR-305
title: "305 -- 7 vanilla row(s)"
pr: "#305"
files: ["`SRC/classTags.h`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/system_of_eqn/linearSOE/diagonal/DiagonalSOE.h`", "`SRC/interpreter/OpenSeesCommands.cpp`", "`SRC/tcl/commands.cpp`", "`SRC/analysis/handler/CMakeLists.txt`", "`SRC/analysis/analysis/DirectIntegrationAnalysis.cpp`"]
table: "main"
legacy_seq: [22, 23, 24, 27, 28, 30, 69]
---
| `SRC/classTags.h` | `// Ladruno` ADR-30: `HANDLER_TAG_LadrunoProjectionHandler`=33001 (first fork handler, private band) | [#305](https://github.com/nmorabowen/OpenSees/pull/305) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno` ADR-30: `getNewConstraintHandler` case `HANDLER_TAG_LadrunoProjectionHandler` → `new LadrunoProjectionHandler()` (+include) so DB/MPI restore reconstructs the handler | [#305](https://github.com/nmorabowen/OpenSees/pull/305) |
| `SRC/system_of_eqn/linearSOE/diagonal/DiagonalSOE.h` | `// Ladruno` ADR-30: add `getDiagonalA()` const accessor — the projector reads the assembled lumped-mass diagonal (the exact M the integrator inverts) before the solver factors it in place | [#305](https://github.com/nmorabowen/OpenSees/pull/305) |
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno` ADR-30: `OPS_ConstraintHandler()` branch `constraints LadrunoProjection` → `OPS_LadrunoProjectionHandler()` (openseespy path) | [#305](https://github.com/nmorabowen/OpenSees/pull/305) |
| `SRC/tcl/commands.cpp` | `// Ladruno` ADR-30: `specifyConstraintHandler()` branch `constraints LadrunoProjection` (classic-Tcl path) | [#305](https://github.com/nmorabowen/OpenSees/pull/305) |
| `SRC/analysis/handler/CMakeLists.txt` | Add `LadrunoProjectionHandler.{cpp,h}`, `LadrunoConstraintProjector.{cpp,h}`, `LadrunoProjectionConsumer.h` to `OPS_Analysis` sources | [#305](https://github.com/nmorabowen/OpenSees/pull/305) |
| `SRC/analysis/analysis/DirectIntegrationAnalysis.cpp` | `// Ladruno` ADR-30: `domainChanged()` now honors the `ConstraintHandler::handle()` / `doneNumberingDOF()` / `Integrator::domainChanged()` error contract (`<0` ⇒ return −1); same guard added to the `setAlgorithm()` / `setIntegrator()` setter paths (mid-session swap). Upstream ignored these returns, so a handler/integrator that DETECTS a bad setup (the projection handler's chain / double-constraint / SP-on-slave / massless / IC-compliance diagnostics, non-projection-aware-integrator refusal) printed its named error but the analysis ran on regardless. Additive — `≥0` on success, `<0` only on a genuine error (verified vs full Zone-A, 0 regressions). | [#305](https://github.com/nmorabowen/OpenSees/pull/305) |
