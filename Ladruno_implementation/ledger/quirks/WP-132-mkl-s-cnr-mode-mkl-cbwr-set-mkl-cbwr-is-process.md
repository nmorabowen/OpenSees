---
wp: WP-132
title: "MKL's CNR mode (mkl_cbwr_set, MKL_CBWR) is PROCESS-WIDE and sticky — -deterministic on one model is on for every later model in the interpreter, and there is n…"
legacy_seq: 502
---
### MKL's CNR mode (`mkl_cbwr_set`, `MKL_CBWR`) is PROCESS-WIDE and sticky — `-deterministic` on one model is on for every later model in the interpreter, and there is no "off" (WP-132)
- **Bites:** `system Pardiso -deterministic` in one model of a pytest session / apeGmsh parameter study / Jupyter kernel leaves every LATER solve in that process on the CNR branch, including models that never asked for it (and all other MKL calls: BLAS/LAPACK too). Setting it inside a test pollutes the rest of the suite.
- **Why:** CNR is MKL library state, not solver state; the solver only requests it. There is no per-call CNR.
- **Workaround/status:** tests run every CNR case in a child process (`tests/test_wp132_deterministic_pardiso.py`). The solver never calls `mkl_cbwr_set` when the flag is absent (byte-identical default), and re-requesting the branch already in force is a no-op. Documented in [[75c_pardiso_solver_recipe]]. See [[LEDGER_vanilla_files]] WP-132 rows. [#864](https://github.com/nmorabowen/OpenSees/pull/864)
