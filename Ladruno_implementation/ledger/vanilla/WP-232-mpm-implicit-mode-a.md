---
wp: WP-232
title: "WP-232f: MPM mode A -- vanilla rows"
pr: "#1083"
date: 2026-11-14
table: "main"
files: ["`SRC/classTags.h`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/element/CMakeLists.txt`", "`SRC/interpreter/OpenSeesCommands.cpp`", "`SRC/interpreter/OpenSeesCommands.h`", "`SRC/tcl/commands.cpp`", "`SRC/analysis/integrator/LoadControl.h`"]
---
| `SRC/classTags.h` | `ELE_TAG_LadrunoMPMPoint 33026` (Ladruno comment on the define): the mode-A MPM particle element's ELE-registry tag. | WP-232f [#1083](https://github.com/nmorabowen/OpenSees/pull/1083) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno WP-232f`: `getNewElement` builds `LadrunoMPMPoint` and `getNewLoadPattern` builds `LadrunoMPMPattern` (with their includes), so a database restore rebuilds a mode-A model. | WP-232f [#1083](https://github.com/nmorabowen/OpenSees/pull/1083) |
| `SRC/element/CMakeLists.txt` | `# Ladruno WP-232f`: `add_subdirectory(ladrunoMPMPoint)`. | WP-232f [#1083](https://github.com/nmorabowen/OpenSees/pull/1083) |
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno WP-232f`: `OPS_LadrunoMPMDL()`, the DL engines' (Python, DL Tcl) `mpm` entry: hands the `cmds` singleton's static / transient analysis, integrators and handler to the shared `OPS_LadrunoMPMWithAnalysis` (mode A's `mpm analyze`), the `OPS_LadrunoStaggeredAnalyze` pattern. | WP-232f [#1083](https://github.com/nmorabowen/OpenSees/pull/1083) |
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno WP-232f`: declares `OPS_LadrunoMPMDL()` for the generated DL / Python bridges of the `mpm` table row. | WP-232f [#1083](https://github.com/nmorabowen/OpenSees/pull/1083) |
| `SRC/tcl/commands.cpp` | `// Ladruno WP-232f`: the classic-Tcl `ladrunoMPM` bridge passes `theStaticAnalysis`, `theTransientAnalysis`, `theStaticIntegrator`, `theTransientIntegrator` and `theHandler` to `OPS_LadrunoMPMWithAnalysis` (classic Tcl keeps its analysis outside `cmds`). | WP-232f [#1083](https://github.com/nmorabowen/OpenSees/pull/1083) |
| `SRC/analysis/integrator/LoadControl.h` | `// Ladruno WP-232f`: inline `ladrunoGetIncrement` / `ladrunoSetIncrement` (Δλ, min, max): the MPM quasi-static driver steps with Δλ = dt and puts the user's values back. `setDeltaLambda` alone is not enough: `newStep` clamps Δλ to [min, max], which default to Δλ. | WP-232f [#1083](https://github.com/nmorabowen/OpenSees/pull/1083) |
