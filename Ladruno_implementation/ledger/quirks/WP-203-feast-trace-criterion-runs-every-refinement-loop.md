---
wp: WP-203
title: "eigen -feast on a building model ran every refinement loop under MKL's default trace stop criterion"
pr: "#1010"
date: 2026-10-22
---
### `eigen -feast` takes minutes on a building model: every solve runs all `maxRefine` refinement loops (MKL's default trace stop never fires)

- **What bites.** YOO SLE port: `eigen -feast 0 fmax` took 406 s for 17 modes on 150k DOF vs 12.4 s for 30 modes with ARPACK + Pardiso. Reproduced on the `bldg.py` frame + ASDShellQ4 floor model (38k DOF, 6 modes in band): 127 s and 42 printed refinement loops (two FEAST calls, the stochastic m0 seed and the solve, 20 + 1 loops each), eigenvalues correct.
- **Why (observed, oneAPI 2025.1 MKL).** `feastinit` leaves `fpm(6)=0`, the relative-trace stop test. On this model MKL kept refining to `fpm(4)` = 20 loops although the trace error it prints is 1e-13..1e-15 from loop 2 (below the 10^-12 of the default `-tol 12`). With `fpm(6)=1` (relative residual) the same solve stops after 3 loops: 21 s, identical eigenvalues. `-tol 10` with the trace test also stops early (15.5 s, 6 loops), so the trace test is simply not meeting 10^-12 on this model, whatever MKL compares internally.
- **Fix (WP-203).** `FeastEigenSolver` sets `fpm[5] = 1` for the seed and the solve; `eigen -feast ... -traceStop` restores MKL's default. Gate B3 (info=2 accepted only if every MKL residual meets the tolerance) is unchanged.
- **Also (same WP).** The command-layer reconcile trimmed only `Domain::setEigenvalues`, not the node eigenvectors, so `modalProperties` after `eigen -feast` aborted ("inconsistent eigenvector matrix and ID", every constraint handler); `FeastEigenSOE::trimDomainModes` trims both.
