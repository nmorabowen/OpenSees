---
wp: PR-877
title: "LadrunoConcrete3D: the deterministic sub-incrementation is DISCONTINUOUS (n boundaries, ladder switches) and its tangent is the CHAIN's, not the last piece's (…"
pr: "#877"
date: 2026-10-10
---
### LadrunoConcrete3D: the deterministic sub-incrementation is DISCONTINUOUS (n boundaries, ladder switches) and its tangent is the CHAIN's, not the last piece's (WP concrete3d-hang-diagnosis, [#877](https://github.com/nmorabowen/OpenSees/pull/877))
- **Bites:** `-subIncr deterministic` (the default) integrates every plastic trial in `n = ceil(f_tr/0.3)` pieces. A chain of `n` and one of `n + 1` pieces (or two ladder levels) are two different, each consistent, integrations of the same increment: `sigma_eff` jumps 0.06-0.40 MPa on 15-41 MPa (0.4-1.5 %) at `n` boundaries and 1-2 % with kappa_p jumps of 1-2.5 at ladder switches near first cracking. The ladder flickers because the first plastic tensile piece (sigma_xx ~ 0.92 ft, kappa_p ~ 0.05) fails its direct return in ~30 % of 1e-12 perturbations. Consequences: nominally identical Gauss points differ (hosting parity ~1e-3, EAS alpha ~1e-11), `F` vs `Q@F` objectivity holds only on one branch, and a finite-difference tangent reference that lets a leg straddle a boundary differentiates the jump.
- **The tangent:** the chain is differentiated forward through the pieces with the kappa-history coupling (`PieceSens`); ignoring that coupling left 6e-2..7e-1 error at n = 2..64. Pin `n` on both FD legs (`Params::subIncrForceN`, `returnMapPieces()`, oracle `subincr_force_n`).
- **Workaround/status:** documented, tests derive tolerances from the measured jump; the root cause (direct-return convergence at the tension onset) is not fixed. `substeps` counts forced chains too.
