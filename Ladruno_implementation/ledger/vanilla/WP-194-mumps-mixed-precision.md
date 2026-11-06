---
wp: WP-194
title: "Mixed-precision MUMPS hooks in MumpsParallelSolver, the two system parsers and CMake"
date: 2026-10-14
pr: "#982"
files: ["`SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.h`", "`SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.cpp`", "`SRC/tcl/commands.cpp`", "`SRC/interpreter/OpenSeesCommands.cpp`", "`CMakeLists.txt`"]
---
| `SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.h` | `// Ladruno WP-194`: `setPrecision()`, static `singlePrecisionAvailable()`, `solveMixed()`; members `precision`, `refineSteps`, `LadrunoMumpsMixed *theMixed` (null = off). Not filed upstream. | [#982](https://github.com/nmorabowen/OpenSees/pull/982) |
| `SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.cpp` | `// Ladruno WP-194`: `solve()` diverts to `solveMixed()` when single (the DMUMPS path below is unchanged and never initialized then); destructor and `setCommunicator()` also tear down the SMUMPS instance; `sendSelf`/`recvSelf` append an `ID(2)` {precision, refine} after the existing controls. Not filed upstream. | [#982](https://github.com/nmorabowen/OpenSees/pull/982) |
| `SRC/tcl/commands.cpp` | `// Ladruno WP-194`: `system Mumps` parses `-precision single|double`, `-mixed`, `-refine k`; refuses bad values, `-refine` without single, a missing value, single on the serial path or on a build without SMUMPS (TCL_ERROR). The parser OpenSeesSP.exe and OpenSeesMP.exe use. Not filed upstream. | [#982](https://github.com/nmorabowen/OpenSees/pull/982) |
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno WP-194`: the same options and refusals in `OPS_MumpsSolver` (openseespy MP). Not filed upstream. | [#982](https://github.com/nmorabowen/OpenSees/pull/982) |
| `CMakeLists.txt` | `# Ladruno WP-194`: `LADRUNO_MUMPS_SINGLE` (AUTO/ON/OFF) detects libsmumps + `smumps_c.h`, prepends the library to `MUMPS_LIBRARIES` and appends `_LADRUNO_SMUMPS` to `MUMPS_FLAG`; `LadrunoMumpsMixed.cpp` in the OpenSeesSP, OpenSeesMP and OpenSeesPyMP source lists. | [#982](https://github.com/nmorabowen/OpenSees/pull/982) |
