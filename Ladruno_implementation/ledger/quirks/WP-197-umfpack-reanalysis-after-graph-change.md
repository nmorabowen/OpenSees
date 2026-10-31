---
wp: WP-197
title: "UmfPack after an sp / element change mid-analysis: DISCARDED entries, then setSize has not been called"
date: 2026-10-15
pr: "#988"
status: "fixed-by WP-197"
---
### `system UmfPack` after an `sp` (or element) is added or removed between `analyze` calls: "DISCARDED" from `addA`, then "setSize has not been called" on every solve

- **What bit.** `UmfpackGenLinSOE::setSize` appended the new CSC pattern to the old one
  (`reserve` + `push_back`, never `clear`). The first `setSize` is fine; the second, after the
  DOF graph changes (`sp` on a free DOF with `constraints Transformation`, a new or removed
  element), left `Ap[0..n]` pointing at the OLD pattern. Symptoms: the once-per-SOE `addA`
  warning "an element entry has no slot in the sparsity pattern; its contribution was
  DISCARDED", `symbolic analysis returns -8` and then "setSize has not been called" on every
  later solve; with a shrinking nnz, `addA` writes past the end of `Ax`. Reported on the
  WP-138 footing replica (issue #981). Upstream master has the same code.
- **Status.** Fixed by WP-197 (`Ap.clear(); Ai.clear(); Ax.assign(nnz, 0.0)`). On a binary built
  before WP-197, call `wipeAnalysis` after every graph change and rebuild the analysis, and
  treat any "DISCARDED" line as a failed run. Pardiso, MUMPS, SparseGEN, BandGen, ProfileSPD,
  FullGeneral and PFEMGeneral rebuild their patterns by index or clear first and are not affected.
