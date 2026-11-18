---
wp: WP-232j
title: "WP-232j: MPM anti-locking -- vanilla rows"
pr: "#1105"
date: 2026-11-17
table: "main"
files: ["`SRC/tcl/commands.cpp`", "`SRC/interpreter/OpenSeesCommands.cpp`", "`SRC/system_of_eqn/linearSOE/mumps/MumpsSOE.h`"]
---
| `SRC/tcl/commands.cpp` | `// Ladruno WP-232j`: the classic-Tcl `ladrunoMPM` bridge also passes `theSOE` and, under `_PARDISO` / `_MUMPS`, the PARDISO or MUMPS SOE's matrix type to `OPS_LadrunoMPMWithAnalysis`, so mode A refuses a symmetric-storage system with `mpm antiLocking fbar` (its tangent is unsymmetric). | WP-232j [#1105](https://github.com/nmorabowen/OpenSees/pull/1105) |
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno WP-232j`: `OPS_LadrunoMPMDL()` passes `cmds->getSOE()` and, under `_PARDISO` / `_MUMPS`, its matrix type the same way (the DL / Python bridge). | WP-232j [#1105](https://github.com/nmorabowen/OpenSees/pull/1105) |
| `SRC/system_of_eqn/linearSOE/mumps/MumpsSOE.h` | `// Ladruno WP-232j`: inline `ladrunoGetMatType()` (the protected `matType`: 0 unsymmetric, 1 SPD, 2 symmetric), so the MPM driver refuses `system Mumps -matrixType 1\|2` with `mpm antiLocking fbar`. | WP-232j [#1105](https://github.com/nmorabowen/OpenSees/pull/1105) |
