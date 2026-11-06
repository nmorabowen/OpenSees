---
wp: LEGACY
title: "linearSolve is NOT the solver cost when the integrator solves too — sum soe.factor/soe.trisolve over EVERY site"
legacy_seq: 218
---
### `linearSolve` is NOT the solver cost when the integrator solves too — sum `soe.factor`/`soe.trisolve` over EVERY site
- **Bites.** You read the `linearSolve` phase as "time in the solver" and understate it. Under `DisplacementControl` the integrator calls `setB(phat); solve()` for the reference displacement `dUhat` in **both** `update()` and `newStep()` — a full extra factorization of the same `K`, booked under those phases, not under `linearSolve`. Measured on lane A (ADR-75b L3-0): `soe.factor` summed over all sites = **505.7 ms > `linearSolve` = 449.4 ms**, so true solver work is **7.86% of step, not the 5.65% `linearSolve` reports**. ADR-40b found the same thing on lane E at far greater severity (59% of step in factorization, two-thirds of it outside `linearSolve`) — and a later report still walked into it, which is why this is its own row.
- **Workaround/status:** to get true solver cost, sum `soe.factor + soe.trisolve` over **every** site in the tree, not the `linearSolve` phase. Affected integrators are any that solve inside `update`/`newStep` — `DisplacementControl` confirmed; `ArcLength`/`MinUnbalDispNorm` share the shape. *2026-07-25 (ADR-75b L3-0 adversarial review pass 2).*
- ✅ **The PARDISO half of this row is RESOLVED** — see the next row. It used to read "you try the same cross-check on a PARDISO run, get `soe.factor = 0.00%`, and conclude no factorization cost", because the scopes lived in `UmfpackGenLinSolver.cpp` only. PARDISO has had its own brackets since 2026-07-27.
