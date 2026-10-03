---
wp: ADR-43
title: "ADR43 P3a -- 2 vanilla row(s)"
files: ["`SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.{h,cpp}`", "`SRC/interpreter/OpenSeesCommands.cpp`"]
table: "main"
legacy_seq: [261, 262]
---
| `SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.{h,cpp}` | `// Ladruno` ADR43 P3a: make the solver honor an arbitrary MPI communicator — the upstream `mpi_comm` constructor argument was accepted and **silently discarded** while `initializeMumps()` hardcoded `id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD)` and probed rank/size on WORLD (D2 spike Finding 1, [[LEDGER_quirks]]). Adds an `MPI_Comm theComm` member (defaults WORLD → default behavior byte-identical) + `setCommunicator(MPI_Comm)` (tears down any live MUMPS instance, forces re-init); `initializeMumps()` now uses `theComm` for `comm_fortran` and rank/size. The legacy `_OPENMPI` `comm_fortran=0` magic ("MUMPS's own COMM_WORLD") is kept only when `theComm==MPI_COMM_WORLD`. Additive. | ADR43 P3a |
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno` ADR43 P3a: (1) `system('Mumps', '-commSplit', color)` — `MPI_Comm_split(WORLD, color, worldRank)`, `solver->setCommunicator(subComm)`, and per-group SOE channel wiring (group-local rank 0 = MUMPS host + channel hub, `MPI_Channel` built from `MPI_Allgather`-ed world ranks; channels process-lifetime since the SOE frees only its pointer array). Collective: every rank must call `system(...)` with some color. (2) `OPS_MumpsSolver` option-loop guard `> 2` → `> 1` — the old bound silently skipped parsing when exactly one `-opt value` pair remained, i.e. the common `system('Mumps','-ICNTL14',n)` spelling never parsed its option. +`#include <MPI_Channel.h>`. | ADR43 P3a |
