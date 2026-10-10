---
wp: WP-132
title: "WP-132 (TIMs F22) — deterministic mode for the implicit static path"
pr: "#864"
status: "draft"
section: "table"
legacy_seq: 9
---
| **WP-132 (TIMs F22) — deterministic mode for the implicit static path** ([[75c_pardiso_solver_recipe]] Trap 7 → "The deterministic mode") — `system Pardiso -deterministic` (MKL CNR, AUTO branch; keeps an `MKL_CBWR` the launcher fixed) and `-cbwr <AUTO\|COMPATIBLE\|AVX2\|AVX512\|…[,STRICT]>` (explicit branch for mixed-CPU cross-node runs, TIMs §1.6; implies -deterministic), same spelling in Tcl and Python. `mkl_cbwr_set` at `system` time + PARDISO `iparm(34)` = MKL thread count per symbolic phase + a one-time runtime notice read back from MKL (the splash line is compile-time). Guide lists what stays order-dependent (WP-107 loop = integer-only reduction, bit-identical; SANISAND refused from it; `-implex` process-wide counters; PFEM omp compiled out; MUMPS/MPI not covered). | solver option (vanilla solver) | — (no class tag) | `SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSolver.{h,cpp}`, `SRC/tcl/commands.cpp`, `SRC/interpreter/OpenSeesCommands.cpp`, banner, `Ladruno_implementation/75c_pardiso_solver_recipe.md`, `tests/test_wp132_deterministic_pardiso.py` | **draft** | [#864](https://github.com/nmorabowen/OpenSees/pull/864) |
