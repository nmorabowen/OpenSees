---
wp: WP-216
title: "ASDShellQ4 forces depend at round-off on the element's call history; modalProperties -return re-normalizes"
date: 2026-10-27
pr: "#TBD"
---
### Per-mode ASDShellQ4 results are not bit-reproducible across call histories

- **Bites:** two per-mode RSA loops that read the same quantities with a
  different sequence of calls (an extra `eleForce`, another mode order, a
  native command that updates only some elements) differ at round-off. On a
  sum that cancels (the vertical reactions of a core wall under overturning)
  that is ~1e-12 of the RESULT, and ~1e-12 of the term scale on a stiff
  residual pseudo-mode that follows a mode of much larger amplitude. A 1e-12
  gate relative to the result fails for reasons unrelated to the code under test.
- **Why:** `ASDShellQ4::AGQIupdate` advances the condensed internal DOFs from
  `KQU` and `Q_residual` as left by the element's LAST `calculateAll` call
  (any of update / getResistingForce / getTangentStiff); `revertToLastCommit`
  restores `U` and `Q`, not those. Separately, `modalProperties -return`
  without `-unorm` re-stores the modal properties with another eigenvector
  scaling: φ s Γ is unchanged mathematically, but later RSA calls move by
  ~1e-14.
- **Fix/status:** not a bug to fix in the element (internal-DOF Newton by
  design). Gate such comparisons relative to the scale of the summed terms,
  and measure the oracle's own spread (run it twice with different histories);
  `tests/test_wp216_rsa_per_mode.py` does both. Read modal properties back
  with the same flags (`modalProperties -unorm -return`). Learned 2026-10-27 (WP-216).
