---
wp: PR-345
title: "345 -- 11 vanilla row(s)"
pr: "#345"
files: ["`SRC/classTags.h`", "`SRC/interpreter/OpenSeesCommands.cpp`", "`SRC/tcl/commands.cpp`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/analysis/handler/CMakeLists.txt`", "`SRC/domain/domain/Domain.{h,cpp}`", "`SRC/domain/CMakeLists.txt`", "`SRC/interpreter/OpenSeesOutputCommands.cpp`", "`SRC/interpreter/OpenSeesCommands.h`", "`SRC/interpreter/PythonWrapper.cpp`", "`SRC/interpreter/TclWrapper.cpp`"]
table: "main"
legacy_seq: [31, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51]
---
| `SRC/classTags.h` | `// Ladruno` ADR-39: `HANDLER_TAG_LadrunoContactHandler`=33002 (ContactDomain handler, fork private band, after 33001) | [#345](https://github.com/nmorabowen/OpenSees/pull/345) |
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno` ADR-39: `OPS_ConstraintHandler()` branch `constraints LadrunoContact` → `OPS_LadrunoContactHandler()` (openseespy path) | [#345](https://github.com/nmorabowen/OpenSees/pull/345) |
| `SRC/tcl/commands.cpp` | `// Ladruno` ADR-39: `specifyConstraintHandler()` branch `constraints LadrunoContact` (classic-Tcl path) | [#345](https://github.com/nmorabowen/OpenSees/pull/345) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno` ADR-39: `getNewConstraintHandler` case `HANDLER_TAG_LadrunoContactHandler` → `new LadrunoContactHandler()` (+include) for DB/MPI restore | [#345](https://github.com/nmorabowen/OpenSees/pull/345) |
| `SRC/analysis/handler/CMakeLists.txt` | `// Ladruno` ADR-39: add `LadrunoContactHandler.{cpp,h}` + `LadrunoContactFE.{cpp,h}` to `OPS_Analysis` sources | [#345](https://github.com/nmorabowen/OpenSees/pull/345) |
| `SRC/domain/domain/Domain.{h,cpp}` | `// Ladruno` ADR-39 P1b: own an optional `LadrunoContactDomain *theContactDomain` (declared LAST → no `-Wreorder`; init `0` in all 4 ctors; `delete` in `~Domain` AND `Domain::clearAll()` (the wipe path — mirrors the ADR-30 `theEQs->clearAll()` leak-fix; domainChanged runs `AnalysisModel::clearAll` so the engine SURVIVES re-analysis); `set/getLadrunoContactDomain`); + `Domain::commit()` → `theContactDomain->commit()` and `Domain::revertToLastCommit()` → `theContactDomain->revertToLastCommit()` (the integrator-agnostic contact-state commit/revert choke point — design-gate BLOCKER-1/2) | [#345](https://github.com/nmorabowen/OpenSees/pull/345) |
| `SRC/domain/CMakeLists.txt` | `// Ladruno` ADR-39 P1b: `add_subdirectory(contact)` (new `OPS_Domain` ContactDomain subsystem) | [#345](https://github.com/nmorabowen/OpenSees/pull/345) |
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno` ADR-39 P1b: `OPS_LadrunoContactSurface` / `OPS_LadrunoContact` / `OPS_LadrunoContactInfo` command bodies (lazily create + populate the Domain-owned `LadrunoContactDomain`) | [#345](https://github.com/nmorabowen/OpenSees/pull/345) |
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno` ADR-39 P1b: declare the 3 contact OPS_ command fns | [#345](https://github.com/nmorabowen/OpenSees/pull/345) |
| `SRC/interpreter/PythonWrapper.cpp` | `// Ladruno` ADR-39 P1b: `Py_ops_*` wrappers + `addCommand` for `contactSurface`/`contact`/`ladrunoContactInfo` | [#345](https://github.com/nmorabowen/OpenSees/pull/345) |
| `SRC/interpreter/TclWrapper.cpp` | `// Ladruno` ADR-39 P1b: `Tcl_ops_*` wrappers + `addCommand` for `contactSurface`/`contact`/`ladrunoContactInfo` (dual-wired w/ Python per the P1a code-gate lesson) | [#345](https://github.com/nmorabowen/OpenSees/pull/345) |
