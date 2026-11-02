---
wp: WP-221
title: "eigen ignored a failed linear solve; an unsymmetric K is wrong for eigen under any storage"
pr: "#1059"
date: 2026-11-02
---
### `eigen` returns 0 although the linear system failed to solve; eigenvalues of a loaded corotational / non-associated model change from run to run, even on `system Pardiso` (full CSR)

- **An unsymmetric K is wrong for `eigen` whatever the `system` stores.** The default eigen solver is
  ARPACK's symmetric Lanczos (`dsaupd`). Switching `-matrixType 2` to full storage fixes `analyze`,
  not `eigen`: on a 6-element corotational cantilever under a tip load (YOO probe, 2026-11-02),
  mode 1 came out 143.5 / 152.3 / 149.3 / 144.0 on four identical runs under `-matrixType 0`.
  WP-221 makes ArpackSOE check every element K of every eigen and warn once;
  `ladrunoSymmetryGuard -eigen` reads it, and the system's `-strictSymmetry` fails the eigen.

- **ArpackSolver ignored `LinearSOE::solve()`.** In the shift-invert loop the ido -1 branch stored the
  return in `ierr` and never read it; the ido 1 branch did not store it. A refused or failed
  factorization left ARPACK iterating on whatever X held, and `eigen` could return 0 with modes of a
  system that was never solved. The analysis-level eigen loops do the same with `addA`: they set
  `result = -2` on an addA failure and then return 0. WP-221 makes ArpackSolver return -1 on a failed
  solve; it does not change the analysis loops. **A failure meant to stop `eigen` must therefore come
  from the SOE's `solve()`, not from `addA`.** That is where `-strictSymmetry` refuses.
- **`eigen` also goes through the active system's half-storage guard:** ArpackSOE forwards `zeroA`/`addA` to the
  `system` SOE, so `system Pardiso -matrixType 2` checks the eigen's K. But without
  `-strictSymmetry` the guard samples (3 of every 64 assemblies). It is armed after `setSize`, which
  an eigen on a new or changed domain triggers, and an eigen after a long analysis on an unchanged
  pattern can fall in a sleep window and go unchecked. `ladrunoSymmetryGuard checkedAssemblies`
  tells which; `-strictSymmetry` checks every assembly.
- **Not covered:** `-fullGenLapack`, `-symmBandLapack` and `-feast` assemble into their own storage,
  not the `system` SOE, so no guard applies. `DistributedProfileSPDLinSOE` (SP/MP) overrides
  `solve()`, so its strict refusal is not wired (classic `system ProfileSPD -strictSymmetry` warns
  there).
