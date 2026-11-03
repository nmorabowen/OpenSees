---
wp: WP-223
title: "PDelta geomTransf is not linear about a preloaded state: the increment carries dP x dDelta"
date: 2026-11-02
pr: "TBD"
---
### `geomTransf PDelta` is not exactly linear about its committed state: an imposed increment carries a dP·dΔ term

- **Bites:** anything that assumes an elastic P-Δ frame responds linearly to a
  trial displacement imposed about a preloaded committed state (scaling one
  RSA evaluation to another amplitude, `responseSpectrumPerMode -linear`,
  finite-difference tangents). The response to u and to 2u differ by more than
  round-off, and the gap grows with the amplitude.
- **Why:** `PDeltaCrdTransf{2d,3d}::getGlobalResistingForce` forms the P-Δ
  shear from the CURRENT basic axial force `pb(0)` (committed P plus the
  trial increment) times the trial chord drift. The increment over the
  committed state is P₀·dΔ + dP·Δ₀ + dP·dΔ: linear in the first two, second
  order in the last. The tangent drops dP·dΔ, so `eigen` is unaffected.
  `ASDShellQ4 -pDelta` is different: after `loadConst` it freezes the membrane
  forces (WP-209) and is exactly linear.
- **Fix/status:** by design (consistent second-order equilibrium). Measured
  on a WP-223 elastic frame under RSA amplitudes: see ADR 214 D10.
  `responseSpectrumPerMode -linear` measures this with its self-check and
  refuses above `-linearTol`. Learned 2026-11-02 (WP-223).
