---
wp: WP-210
title: "A radial-return (p, q) tangent that special-cases q_tr = 0 by dropping the deviatoric block gives ZERO shear stiffness under isotropic loading"
date: 2026-11-06
---
### Isotropic consolidation on a cap model gives a singular element matrix (or Newton stalls) at the very first plastic step

- **Bites:** an invariant-space return map reconstructs `sigma = p I + (q/q_tr) s_tr`
  and builds the tangent from `n_s = s_tr / q_tr`. At `q_tr = 0` the direction is
  undefined, so the natural guard `if (q_tr > tol) {...} else { C = m (x) grad_p : C_e }`
  keeps only the volumetric part. The material tangent then has rank 1: C44 = C55 = C66 = 0
  and C11 = C12. Every isotropically consolidated CIU/CID test, and every near-isotropic
  Gauss point in a BVP, hits it. LadrunoHS shipped this (`LadrunoHSKernel.h`, regime-3
  reconstruction) until WP-210 round 3.
- **Why the FD gate missed it:** the cap FD case used a strain with shear components, so
  `q_tr > 0` and the q = 0 branch never ran. A tangent test needs a case AT the special
  point (purely isotropic strain from an isotropic state) where FD perturbs only the
  shear strains off zero.
- **Fix:** the limit of `q/q_tr` as `q_tr -> 0` is `dq/dq_tr` (the regime's `sens[1][1]`;
  for the HS cap, `p_hat / (p_hat + c_q deps_v^pc)`), so the deviatoric block is
  `sens[1][1] * I_dev`, not zero. Gate: `ladruno_hs_kernel_check.cpp` case
  `cap_isotropic` and the `test_t2c_consistent_tangent_fd_all_regimes` isotropic cap row.
- **Generalises to:** any (p, q) return map (cap models, Cam-clay family, cone apex
  handling). Look for an `else` branch under `q > tol` that builds the tangent from
  `grad_p` alone.
