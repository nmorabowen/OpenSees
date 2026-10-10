---
wp: LEGACY
title: "Quirk: bundled OTHER/LAPACK ships dsygvx but NOT dsygv (Linux link break)"
legacy_seq: 36
---
## Quirk: bundled `OTHER/LAPACK` ships `dsygvx` but NOT `dsygv` (Linux link break)

The fork's in-tree reference LAPACK (`OTHER/LAPACK/`, statically linked by the
Unix/Ubuntu Zone-A build) is a *curated subset* of LAPACK, not the full library.
It provides the **expert/range** symmetric-definite generalized eigensolver
`dsygvx.f` (line ~147 of `OTHER/LAPACK/CMakeLists.txt`) but does **not** provide
the simple driver `dsygv.f`. It also ships `dsyevx` but not `dsyev`, `dggev` but
the `x`-suffixed expert drivers are the ones actually bundled for the symmetric
path.

Consequence: fork code that calls `dsygv_` compiles per-TU fine and links fine
on **Windows** (full MKL resolves `DSYGV`), but fails the **Linux** link with
`undefined reference to 'dsygv_'`. This silently blocked the entire ladruno
Zone-A CI Build step (and thus every downstream runtime test battery) after
`CriticalTimeStep.cpp` (CentralDifferenceLadruno, PR #22) introduced a `dsygv_`
call — even though the recorder/MPCO work was unrelated.

**Rule for fork code:** for a symmetric-definite generalized eigenproblem, call
the bundled **`dsygvx_`** (expert driver), not `dsygv_`. To get just the largest
eigenvalue use `RANGE='I'` with `IL=IU=N`; the value lands in `W[0]`. Mirror the
existing, proven usage in
`SRC/system_of_eqn/eigenSOE/SymmGeneralizedEigenSolver.cpp`. Adding `dsygv.f` to
the bundle instead is the wrong fix: it cascades (dsygv calls `dsyev`, also not
bundled) and edits the upstream LAPACK source list.
*(Fixed in CriticalTimeStep.cpp via `fix/criticaltimestep-dsygvx-link`.)*
