---
wp: LEGACY
title: "MumpsParallelSolver accepts an mpi_comm constructor argument and silently ignores it — always binds to MPI_COMM_WORLD"
legacy_seq: 161
---
### `MumpsParallelSolver` accepts an `mpi_comm` constructor argument and silently ignores it — always binds to `MPI_COMM_WORLD`

- **Bites:** ADR 43 P3 (FEAST-over-sub-communicator, D2 spike, 2026-07-07). Anyone
  assuming `MumpsParallelSOE`/`MumpsParallelSolver` can be pointed at an
  `MPI_Comm_split` sub-communicator today by passing a comm handle — it compiles and
  runs, but silently uses `MPI_COMM_WORLD` for the factorization anyway (wrong ranks
  participate, or it deadlocks/hangs on a partial-world sub-comm).
- **Why:** `MumpsParallelSolver::MumpsParallelSolver(int mpi_comm, int ICNTL7, int
  ICNTL14)` (`MumpsParallelSolver.cpp:54-64`) takes the parameter but never stores it —
  no member is set. `initializeMumps()` (`:93-105`) hardcodes
  `id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD)` on the Intel-MPI path (this fork's
  Windows/oneAPI build), `0` (MUMPS's own WORLD) under `_OPENMPI`; the rank/size probe
  two lines later also reads `MPI_COMM_WORLD` directly. Dead parameter, not a config
  toggle.
- **Fix (SHIPPED, ADR 43 P3a):** `setCommunicator(MPI_Comm)` on the solver stores the
  comm, `MPI_Comm_c2f`s *it* (not WORLD) into `id.comm_fortran`, uses it for the
  rank/size probe, tears down any live MUMPS instance, and clears the SOE's
  `factored` flag. Test hook: `system('Mumps', '-commSplit', color)` (collective —
  every rank must call it). Gate: `feast_d2_spike/p3a_commsplit_gate.py` (4 ranks,
  2 concurrent groups vs serial oracles). **Residual subtlety:** `MPI_Channel`
  hardcodes WORLD/tag-0, so only the MUMPS factor/solve is comm-isolated — the SOE's
  B/X exchange rides WORLD envelopes, safe today via disjoint (src,dst) pairs +
  MPI non-overtaking with phase ordering; true envelope isolation is a P3c-MPI item.
  **Two deliberate scope notes** (both differ from this bullet's pre-ship draft, so
  don't read them as omissions): (1) the comm lives on the SOLVER — it is NOT threaded
  through `MumpsParallelSOE`'s constructor, which still takes
  `(MumpsParallelSolver &theSolvr, int matType)` (`MumpsParallelSOE.cpp:43`); the
  wiring runs `OpenSeesCommands.cpp:5086` → `solver->setCommunicator(subComm)`.
  (2) `-commSplit` is **openseespy / interpreter-Tcl only** — classic Tcl deliberately
  refuses it (`SRC/tcl/commands.cpp:4296-4312`, ADR-75 P2h).
