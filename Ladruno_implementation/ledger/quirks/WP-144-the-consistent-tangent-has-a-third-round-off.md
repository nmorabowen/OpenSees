---
wp: WP-144
title: "The consistent tangent has a THIRD round-off band: NEAR-COALESCENT principal elastic strains (WP-144; companion of the Lode-corner / vertex entry above)"
date: 2026-11-06
---
### The consistent tangent has a THIRD round-off band: NEAR-COALESCENT principal elastic strains (WP-144; companion of the Lode-corner / vertex entry above)
- **Bites:** the spectral tangent divides by the eigenvalue gap, g_ab = (sigma_a - sigma_b)/(eps_a - eps_b). Inside the switch tolerance (|eps_a - eps_b| < 1e-10, `REPEATED_EIG_TOL`) the kernel and O2 use the limit a~_aa - a~_ab, and just outside it the quotient cancels.
  - The (a,b) and (b,a) rows of the limit agree only for EXACTLY repeated eigenvalues (swap symmetry of the isotropic map); inside the tolerance they differ by O(|eps_a - eps_b|). Compressing C4 to 6x6 by the `i <= j` rows is a CONTRACT (sheet 144a 9.4), so kernel and O2 agree on it, but any other summation order moves the result.
  - Measured by the P1 parity drills (sheet 144a 9.4; scratch drills, not tracked): ~1e-8 relative on C on near-coalescent states (a NON-ZERO gap < 1e-6 max|eps~|; an EXACTLY zero gap, axisymmetric / isotropic, is gated at the full 1e-10, the row-convention effect needing a non-zero gap), while the STATE is unaffected.
- **Rule:** gate the tangent in this band separately, like the corner and vertex bands (a looser tangent gate; state quantities stay at the full gate), print the oracle's own 1-ulp self-sensitivity next to the kernel error, and keep a generic (well-separated eigenvalue) path at the full tangent gate so the formula is still verified. Never symmetrise g_ab or average the two rows to "fix" it: that moves parity at the 1e-8 level and is not the contract.
- **Workaround/status:** ✅ Contract in sheet 144a 9.4; the band's threshold lives in `test_kernel_parity.py` beside the corner and vertex bands (WP-144, 2026-10-01).
