---
wp: WP-260
title: "PARDISO iparm[24]=0 is already parallel and iparm[27]=1 requires float* buffers"
date: 2026-11-25
---
### PARDISO `iparm[24]=0` is already parallel and `iparm[27]=1` requires `float*` buffers

- **Bites:** two separate PARDISO assumptions in `PARDISOGenLinSolver.cpp`: (1) the comment at `:367` states `iparm[24] = 0` is a sequential forward/backward solve, and (2) attempting single-precision PARDISO by setting `iparm[27] = 1` (Fortran `IPARM(28)`) while still passing `double*` `a`, `b`, `x` arrays corrupts memory / produces garbage.
- **Why:** per the Intel oneMKL specification, `iparm[24]` (`iparm(25)`) = `0` **already runs parallel** triangular solves (`1` forces sequential; `2` forces matrix-partitioned parallel solve), while `iparm[27]` (`iparm(28)`) = `1` instructs PARDISO to read/write `a`, `b`, and `x` strictly as `float*` and perform all internal operations (including internal refinement) in FP32.
- **Workaround / status (2026-11-25, WP-260 — disclosed):** fix the `:367` comment, target `iparm[23] = 1|10` (two-level parallel factorization) and `iparm[1] = 3` (parallel METIS) for multicore scaling, and implement mixed-precision desktop PARDISO only with a `std::vector<float>` shadow of `A` + `float` `b_f`/`x_f` buffers wrapped in an outer FP64 GMRES-IR (`mkl_rci.h` `dfgmres`) loop.
