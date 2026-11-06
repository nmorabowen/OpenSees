---
wp: LEGACY
title: "The serial MumpsSolver is never compiled in this fork (TIMs F5, 2026-09-07)"
date: 2026-09-07
legacy_seq: 407
---
## The serial `MumpsSolver` is never compiled in this fork (TIMs F5, 2026-09-07)

`CMakeLists.txt` defines `_MUMPS` only for the parallel targets (`OpenSeesSP`,
`OpenSeesMP`, `OpenSeesPyMP`: lines ~996/1217/1292/1408); the serial `OpenSees` /
`OpenSeesPy` targets get no MUMPS at all — ADR-75 P1b kept MUMPS as the CLUSTER
solver and made PARDISO the desktop one (`CMakeLists.txt:588-594`). Measured on
the F4 build: `system Mumps -stats` on `OpenSees.exe` and the pyd both answer
"unknown system type". So the "silent serial `MumpsSolver.cpp` path" the TIMs
note cites (`MumpsSolver.cpp:150-200`, no `printStats`, both parsers construct it
2-arg and `commands.cpp:4341-4345` warns `-stats` is ignored) is real in the
source but UNREACHABLE in any shipped serial binary. Wiring `-stats` there would
be dead, unverifiable code; MUMPS statistics on the desktop come from a one-rank
`mpiexec -n 1 OpenSeesMP` / `openseesmp` run (rank 0 prints them,
`MumpsParallelSolver.cpp:295-319`) — which needs the packaged `dist\openseesmp`
runtime (a no-arg `build.bat`), not the 4-target build. Linking MUMPS into the
serial targets is an ADR-75 policy reversal for the owner, not a "small" WP.
