---
wp: WP-221
title: "Half-storage asymmetry guard: ladrunoSymmetryGuard query + -strictSymmetry"
date: 2026-11-02
status: "draft"
pr: "#1059"
---
| **WP-221 — `ladrunoSymmetryGuard` + `-strictSymmetry`** (`system Pardiso -matrixType 1\|2`, `system ProfileSPD`): the ADR-75 P1d guard's state is returned to openseespy / DL Tcl / classic Tcl; the strict flag checks every assembly and refuses the solve once latched, so `analyze`/`eigen` fail instead of solving the reflected system; ArpackSOE checks every eigen K under ANY storage (`ladrunoSymmetryGuard -eigen`; ARPACK is symmetric, so full CSR is no remedy for eigen) | solver diagnostic | — | `SRC/system_of_eqn/linearSOE/LadrunoSymmetryGuard.h` (new), `LadrunoSymmetryGuardQuery.h` (new), `PARDISOGenLinSOE.*`, `ProfileSPDLinSOE.*`, `ProfileSPDLinDirectSolver.cpp`, `ArpackSolver.cpp`, `ArpackSOE.*`, `LinearSOE.h`, `LadrunoCommandTable.h`, `tests/test_wp221_symmetry_guard*.py`, `tests/tcl/wp221_symmetry_guard.tcl` | draft | #1059 |
