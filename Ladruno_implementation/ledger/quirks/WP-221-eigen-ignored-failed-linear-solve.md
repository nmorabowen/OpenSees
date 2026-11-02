---
wp: WP-221
title: "eigen ignored a failed linear solve, and its K assembly is only sampled by the asymmetry guard"
pr: "TBD"
date: 2026-11-02
---
### `eigen` returns 0 although the linear system failed to solve; `system Pardiso -matrixType 2` eigen of an unsymmetric K gives no error

- **ArpackSolver ignored `LinearSOE::solve()`.** In the shift-invert loop the ido -1 branch stored the
  return in `ierr` and never read it; the ido 1 branch did not store it. A refused or failed
  factorization left ARPACK iterating on whatever X held, and `eigen` could return 0 with modes of a
  system that was never solved. The analysis-level eigen loops do the same with `addA`: they set
  `result = -2` on an addA failure and then return 0. WP-221 makes ArpackSolver return -1 on a failed
  solve; it does not change the analysis loops. **A failure meant to stop `eigen` must therefore come
  from the SOE's `solve()`, not from `addA`.** That is where `-strictSymmetry` refuses.
- **`eigen` does go through the active system's guard:** ArpackSOE forwards `zeroA`/`addA` to the
  `system` SOE, so `system Pardiso -matrixType 2` checks the eigen's K. But without
  `-strictSymmetry` the guard samples (3 of every 64 assemblies). It is armed after `setSize`, which
  an eigen on a new or changed domain triggers, and an eigen after a long analysis on an unchanged
  pattern can fall in a sleep window and go unchecked. `ladrunoSymmetryGuard checkedAssemblies`
  tells which; `-strictSymmetry` checks every assembly.
- **Not covered:** `-fullGenLapack`, `-symmBandLapack` and `-feast` assemble into their own storage,
  not the `system` SOE, so no guard applies. `DistributedProfileSPDLinSOE` (SP/MP) overrides
  `solve()`, so its strict refusal is not wired (classic `system ProfileSPD -strictSymmetry` warns
  there).
