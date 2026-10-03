---
wp: LEGACY
title: "MUMPS Error −9 (\"work array too small\") is RANK-COUNT-dependent and non-monotonic — and ICNTL14 strongly perturbs wall time, so hold it uniform across any np s…"
legacy_seq: 5
---
### MUMPS Error −9 ("work array too small") is RANK-COUNT-dependent and non-monotonic — and ICNTL14 strongly perturbs wall time, so hold it uniform across any np sweep
- **Bites:** parallel implicit runs (`system Mumps`) that work at np=1/2/8 can **fail at step 0 with INFO(1)=−9 at np=4** (or any specific rank count) — the workspace under-prediction depends on how the partition shapes the distributed factorization, not on problem size alone. Worse for benchmarking: the `-ICNTL14` value needed to survive also changes performance dramatically (measured: np=8 wall 90.5 s at `-ICNTL14 200` vs 33.0 s at `-ICNTL14 2000` on the same model), so a sweep with per-np ICNTL14 values produces meaningless comparisons.
- **Why:** ICNTL(14) is MUMPS's % workspace relaxation over its analysis-phase estimate; the estimate quality varies with the partition/ordering interaction. OpenSees defaults it to 20 (`OpenSeesCommands.cpp:4252`); the −9 handler prints the "make ICNTL14 larger" hint (`MumpsParallelSolver.cpp:168`). Default-20 failed at np=4 on a plain 5488-hex/18.9k-DOF block; 200 did NOT rescue it; 2000 did.
- **Workaround/status (2026-07-06, ADR-40 rank-3 MUMPS scaling measurement):** for any np sweep set `-ICNTL14` high (≥200, be ready for 2000) and **uniform across all rank counts**; verify every config actually completed (a DNF leaves stale output/h5 files from prior runs — delete artifacts before a fresh sweep). Measured sweep + table in [[40b_phase0_dominance_report]] §MUMPS addendum.
