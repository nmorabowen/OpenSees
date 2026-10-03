---
wp: ADR-70
title: "T6 quirks: constant-mean B-bar/F-bar is RANK-DEFICIENT; nodal-lumped corner masses are ZERO (ADR 70 P3)"
legacy_seq: 174
---
## T6 quirks: constant-mean B-bar/F-bar is RANK-DEFICIENT; nodal-lumped corner masses are ZERO (ADR 70 P3)

Two hard-won T6 (quadratic triangle) facts from building `LadrunoLST`:

1. **Constant element-mean dilatation (B-bar / centroid-sampled F-bar) loses
   rank on the T6.** The two quadratic CONFORMAL displacement modes (Re/Im of
   z²: u=(x²−y², −2xy)-type) have identically zero deviatoric strain, and
   their linear dilatation has zero element mean (mean over the 3-interior-
   point rule = centroid value for a linear integrand). Any constant-mean
   averaging therefore assigns them ZERO energy: a free element shows **5**
   zero-energy modes (3 RBM + 2 spurious; stacked-B̄ rank 7 of the required 9).
   The Q4 is immune only because z² fields are not in its bilinear space —
   "J varies so F-bar can average" is NOT sufficient for rank. Caught by the
   locked single-free-element `assert_zero_energy` T0 gate; confirmed by numpy
   rank and compiled `eigen`. Cure DECIDED by the ADR-70 P4 spike (2026-07-11):
   disjoint 2-triangle patch macro-element — patch-constant J̄ for T3 (dSNPO
   §15.1.9), patch-P1 projected dilatation for T6 (see the [[#Volumetric-
   projection traps on triangles (ADR 70 P4 spike)|P4 quirk entry]]). Until it
   ships, triangles are `std` and near-incompressible plane problems use the
   quad `bbar`.

2. **SixNodeTri-style plain N-lumping gives EXACTLY ZERO corner masses on the
   T6** — ∫N_corner over the 3-interior-point rule integrates to 0, and goes
   NEGATIVE on distorted elements (adversarial gate verified both). Unusable
   for explicit dynamics (diagonal M⁻¹). `LadrunoLST::getMass` therefore
   deliberately DIVERGES from upstream and uses **HRZ lumping** (∫ρN²_a dV
   rescaled to the exact total — strictly positive), the same call ADR-72 made
   for the H20 (row-sum corners were −M/8 there). The reduce-to-`SixNodeTri`
   gate is static, so the anchor is unaffected. Documented in `LadrunoLST.h`;
   gated by a free-free finite-eigenvalue test.

Also: upstream `SixNodeTri` registers in the shared functionMap as **`tri6n`**
(not "SixNodeTri"), and its `setPressureLoadAtNodes` carries a side-61
`dx61 = x4-x6` copy-paste typo that breaks closed-contour equilibrium of the
consistent pressure load — `LadrunoLST` fixes it (`x1-x6`) and pins the fix
with a zero-net-resultant test.
