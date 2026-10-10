---
wp: WP-130
title: "The sign-corrected CPPM tangent is still NOT the consistent tangent -- four separate error sources (WP-130 review r1)"
legacy_seq: 515
---
### The sign-corrected CPPM tangent is still NOT the consistent tangent -- four separate error sources (WP-130 review r1)
- **Bites:** "`-cppmTangent fixed` = the algorithmic tangent, one iterate stale" (the WP-130 first
  claim) is wrong in degree. A reviewer's FD campaign (`review130/p1_tangent.py`) separates:
  (a) **staleness**, and the local convergence test `||R|| < TolR (1 + ||R0||)` mixes strain,
  back-stress and stress units, so R1 converges only to ~TolR*|F|: at the default TolR 1e-7 the
  g12 column is off by 0.27-0.53 relative (the first "1.2e-3" was one lucky state);
  (b) the **void-ratio dependence** eps -> e (inVar(37)) -> psi -> M^b, M^d, b0 is absent from
  dR/deps: 1.3-2.3e-4 on the volumetric column at p ~ 170, 1.24e-3 at p ~ 1.9, TolR-independent;
  (c) the **low-p D_factor derivative** in `NewtonSol` had the wrong sign and dropped the D < 0
  branch (`- one3*Macauley(D)*be*t/(1+t)^2` vs the correct `+ one3*(D/D_factor)*be*t/(1+t)^2`);
  it enters the local Jacobian too, i.e. Newton's rate;
  (d) after a SUCCESSFUL halving the tangent handed out is the SECOND half-increment's (O(1)).
- **Status (WP-130, #868):** (c) FIXED under `-cppmTangent fixed` (the LadrunoSANISAND default;
  vanilla ManzariDafalias and `-cppmTangent vanilla` keep the vanilla expression, byte-identical);
  it moves the local iterates, not the root, where p < 0.05 P_atm -- the fixed-default baseline
  was re-pinned. (a), (b), (d) documented, not fixed: (a) needs a unit-consistent local norm (a
  change to the return map's convergence semantics), (b) a d/de block in NewtonSol, (d) a
  chain-rule product across the halving tree. So the order-1.24 global convergence on the bearing
  deck has several causes, not one.
