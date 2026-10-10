---
wp: WP-132
title: "PARDISO's iparm is ZEROED at every symbolic phase in PARDISOGenLinSolver::solve() — any iparm option (CNR iparm[33], …) must be re-applied inside that block, n…"
legacy_seq: 504
---
### PARDISO's `iparm` is ZEROED at every symbolic phase in `PARDISOGenLinSolver::solve()` — any iparm option (CNR `iparm[33]`, …) must be re-applied inside that block, not once in a setter (WP-132)
- **Bites:** the symbolic block (`needsSymbolic || !init`) runs `for (i<64) iparm[i]=0` and rebuilds the control array on every new sparsity pattern (`setSize` → domain change, re-emission, element removal). An option written into `iparm` from a setter survives only until the first pattern change — here it would silently drop CNR (`iparm(34)` back to 0) mid-run.
- **Workaround/status:** WP-132 stores the request in a member (`cnrBranch`) and writes `iparm[33]` inside the symbolic block, like `-stats` does for `iparm[17]/[18]`. Same rule for any future iparm lever. [#864](https://github.com/nmorabowen/OpenSees/pull/864)
