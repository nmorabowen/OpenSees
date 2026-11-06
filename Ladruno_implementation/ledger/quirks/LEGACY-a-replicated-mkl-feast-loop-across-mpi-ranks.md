---
wp: LEGACY
title: "A REPLICATED MKL FEAST loop across MPI ranks DEADLOCKS unless MKL is pinned single-threaded"
legacy_seq: 166
---
### A REPLICATED MKL FEAST loop across MPI ranks DEADLOCKS unless MKL is pinned single-threaded

- **Bites:** any design that runs an MKL routine (FEAST `dfeast_srci`/`dfeast_scsrgv`,
  or any threaded LAPACK) REPLICATED on every MPI rank in lockstep, with
  collective work between the rank-local MKL calls. This is exactly the ADR 43
  P3c-MPI (L3-only) model: every rank runs the same `dfeast_srci` outer loop and
  they cooperate on ONE distributed MUMPS inner solve at each `ijob=10/11`.
- **Why:** multi-threaded MKL LAPACK is **not bitwise-reproducible across
  processes** (dynamic work scheduling, thread count, CPU affinity all vary the
  reduction order → last-ULP differences). Two failure modes, both silent until
  they aren't: (1) the reduced m0×m0 eig inside `dfeast_srci` rounds differently
  per rank → the refinement-loop **count** diverges → one rank returns `ijob=0`
  (done) while another issues another collective inner solve → **hard deadlock**;
  (2) the stochastic auto-seed (`dfeast_scsrgv` with `fpm[13]=2`, taken when
  `-m0` is unset) rounds `mEst` differently per rank → different `m0` → different
  `nrhs` at the first solve → `MPI_Bcast` length mismatch → MPI abort. MKL often
  runs small dense solves single-threaded (below an internal size threshold), so
  it **usually works and intermittently hangs** — the worst profile. Caught by
  the P3c-MPI adversarial gate (the differential test masked it by pinning
  `MKL_NUM_THREADS=1` in its env; the shipped code enforced nothing).
- **Fix (what actually shipped, after a false start):** the obvious fix —
  `mkl_set_num_threads_local(1)` around the replicated span — **SEGFAULTS at
  runtime in this MP MKL DLL layout** (c0000005 on the first call; the global
  `mkl_set_num_threads` too). Bisected: guard disabled ⇒ distributed solve runs
  clean; guard enabled ⇒ crash. The symbol links (import lib has it) but the
  call dies — never exercised in the serial build because the guard is inactive
  there. So do NOT reach for MKL thread control here. Instead, attack the
  **catastrophic** desync directly: broadcast the stochastic auto-seed `m0` from
  rank 0 (`LadrunoFeastInnerSolve::agreeInt` → `MPI_Bcast`, serial default =
  identity) right after the estimate, so the FIRST distributed solve's block
  width (hence its `MPI_Bcast` length) is identical on every rank. The *other*
  desync source — the reduced m0×m0 eig inside `dfeast_srci` rounding
  differently under threading — does **not** bite in practice: m0 is small
  (tens), below MKL's LAPACK threading threshold, so it runs single-threaded and
  stays bit-identical; the enlargement decisions key off the collective solve's
  `m`/`info` (identical on all ranks post-broadcast). **PROVEN**: the P3c-MPI
  gate passes at `MKL_NUM_THREADS=4` (multi-threaded) with only the broadcast,
  no thread pin. For pathologically large bands (m0 in the hundreds, where MKL
  might thread the reduced eig) set `MKL_NUM_THREADS=1` — documented, not
  enforced (the enforcement API is the thing that crashes). See `agreeInt` +
  the m0-broadcast in `FeastEigenSolver::solve()`.
