---
wp: WP-132
title: "-cbwr AVX2 (any instruction-set CNR branch) is REFUSED on an AMD CPU — only AUTO and COMPATIBLE exist there; pin cross-node runs to COMPATIBLE (WP-132)"
legacy_seq: 505
---
### `-cbwr AVX2` (any instruction-set CNR branch) is REFUSED on an AMD CPU — only `AUTO` and `COMPATIBLE` exist there; pin cross-node runs to `COMPATIBLE` (WP-132)
- **Bites:** the obvious cross-node recipe for TIMs §1.6 ("pin every node to AVX2") fails on this box: on an AMD Ryzen AI 7 PRO 350, `mkl_cbwr_set` returned `MKL_CBWR_ERR_UNSUPPORTED_BRANCH` (-3) for SSE4_2, AVX, AVX2, AVX512 and AVX10 (measured 2026-09-26). `AUTO`, `COMPATIBLE` and `AUTO,STRICT` were accepted. `mkl_cbwr_get_auto_branch()` also returns `AUTO` itself there, not a concrete branch.
- **Why:** MKL's instruction-set CNR branches are defined for Intel CPUs.
- **Workaround/status:** use `-cbwr COMPATIBLE` for a branch every x86 node can run. Keep `-cbwr AVX2` only for all-Intel clusters. Tcl stops on the refusal and Python warns; both name AUTO/COMPATIBLE. The notice omits the `(-> branch)` suffix when MKL reports only AUTO. Documented in [[75c_pardiso_solver_recipe]]. [#864](https://github.com/nmorabowen/OpenSees/pull/864)
