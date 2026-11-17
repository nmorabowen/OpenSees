---
wp: WP-232j
title: "WP-232j: MPM anti-locking -- vanilla rows"
pr: "#1105"
date: 2026-11-17
table: "main"
files: ["`SRC/tcl/commands.cpp`", "`SRC/interpreter/OpenSeesCommands.cpp`"]
---
| `SRC/tcl/commands.cpp` | `// Ladruno WP-232j`: the classic-Tcl `ladrunoMPM` bridge also passes `theSOE` and, under `_PARDISO`, the PARDISO SOE's `getMatType()` to `OPS_LadrunoMPMWithAnalysis`, so mode A refuses a symmetric-storage system with `mpm antiLocking fbar` (its tangent is unsymmetric). | WP-232j [#1105](https://github.com/nmorabowen/OpenSees/pull/1105) |
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno WP-232j`: `OPS_LadrunoMPMDL()` passes `cmds->getSOE()` and, under `_PARDISO`, its `getMatType()` the same way (the DL / Python bridge). | WP-232j [#1105](https://github.com/nmorabowen/OpenSees/pull/1105) |
