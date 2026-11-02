---
wp: WP-219
title: "Concrete02 returned the previous trial's stress at its committed strain: the force-based 'class E' stall"
pr: "#1056"
date: 2026-11-01
---
### Force-based element stops with `dW` just above its tolerance, at any tolerance, chaotically in step size and solver (cyclic RC columns with Concrete02)

- **Bites:** `LadrunoForceBeamColumn` and vanilla `forceBeamColumn` with `Concrete02` fibers, under
  a constant axial load. frame-ladruno-validation blind recipe: U4 stopped at increment 297, L407
  at 791. `-iter 10 1e-4` only moved the stops. A displacement-based element hid it.
- **Why:** `Concrete02::setTrialStrain` returned early on `|deps| < DBL_EPSILON` with `sig` and `e`
  still holding the previous TRIAL (commit state was restored for `ecmin` and `dept` only). The stress
  therefore jumps by a whole iterate at the committed strain. A section whose strain barely moves
  (IP at the near-zero-moment end, constant N) lands there, and the element's internal iteration
  cycles between the two values forever. Debug trace on U4: the same fibers at IP 4 flip tangent
  17726 and 18473 (envelope and reload line) in an exact two-cycle, `dW` 5.9388e-12 / 5.9142e-12.
  `Concrete01` and `Concrete04` reset their trial stress first and do not have this.
- **The misleading signature:** the cycle amplitude follows the size of the step, so failures are
  always "just above" the tolerance: dW/W = 4e-18 at `tol` 1e-12 and 1.7e-14 at relative 1e-14.
  A convergence test relative to the element's work scale (W = sum |Se·v|) was tried at 1e-18 to 1e-12:
  histories were unchanged (to 5e-6 of peak) but stops only moved. Halving the reversing
  correction gave a three-cycle, and freezing the section flexibilities still cycled at dv = 1.5e-14,
  because the jump is in the stress, not the tangent. The energy test was never the defect.
- **Fix (WP-219):** return the committed response there (`sig = sigP; e = eP`). The element is not
  changed, and neither is its tolerance or the 1 % cap of its compatibility check.
- **Result, same recipe:** U4, L407, L815, L1015 now complete. Vanilla forceBeamColumn + Steel02: U4
  297 to done, L407 710 to done, U3 555 to done. Common increments agree to 3.7e-4 of peak or better.
  Remaining stops (U3 1014, U1 1014, L415 5291) are not this: fibers at strain below -1 or
  dW/W of about 1e-7 (crushed sections). L430 stops with no element failure at all.
- **Gate:** `tests/test_concrete02_trial_at_committed_strain.py`.
