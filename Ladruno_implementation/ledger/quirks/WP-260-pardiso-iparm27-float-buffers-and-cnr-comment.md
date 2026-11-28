---
wp: WP-260
title: "PARDISO iparm[27]=1 requires float* buffers; the iparm[24] comment sits inside the CNR block"
date: 2026-11-25
---
### PARDISO `iparm[27]=1` requires `float*` buffers; the `iparm[24]` comment sits inside the `-deterministic` (CNR) block

- **Bites:** (1) setting `iparm[27] = 1` (Fortran `IPARM(28)`, single precision) while still passing the `double*` `a`, `b`, `x` arrays makes PARDISO read and write them as `float`: garbage results or memory corruption. (2) `PARDISOGenLinSolver.cpp:366-368` describes `iparm[23]/[24] = 0` as "classic factorization / sequential forward-backward solve" and as "already the deterministic choices"; the Intel oneMKL `iparm` table describes `iparm[24] = 0` as the parallel solve (`1` sequential, `2` the partitioned parallel solve). Because the text is inside the WP-132 CNR block, "fixing the comment" is not cosmetic: the determinism argument it makes has to be re-checked.
- **Why:** `iparm[27] = 1` switches all input and internal arrays to single precision. The CNR comment was written for the WP-132 `-deterministic` path, where `iparm[33] > 0` turns on deterministic scheduling; whether that covers a parallel solve step was argued from the comment's reading of `iparm[24]`, not measured.
- **Workaround / status (2026-11-25, WP-260 — disclosed, not yet acted on):** quote the oneMKL **2024.2** table text (the fork's oldest MKL) when changing the comment, and re-prove `-deterministic` bit identity run-to-run at several thread counts (see "Threaded MKL PARDISO is NOT byte-reproducible run-to-run"). New parallel options (`iparm[23] = 1|10`, `iparm[24] = 2`, `iparm[1] = 3`) are opt-in, `G-OPT`, and refused under `-deterministic` unless proved there. Mixed-precision PARDISO uses a `std::vector<float>` shadow of `A` plus `float` `b_f`/`x_f`, wrapped in an outer FP64 GMRES-IR loop (`mkl_rci.h` `dfgmres`), and must pass the CI job `MKL 2024.2 compat (WP-148)`.
