---
wp: LEGACY
title: "section LayeredShell bending is exact only on its own midpoint rule — a predictable Σ h³/12 stiffness deficit vs the continuum (≈2% at 5 uniform layers)"
legacy_seq: 156
---
### `section LayeredShell` bending is exact only on its own midpoint rule — a predictable Σ h³/12 stiffness deficit vs the continuum (≈2% at 5 uniform layers)

- **Bites:** elastic cylindrical bending of a LayeredShell with n uniform layers
  undershoots `E·t³/12(1−ν²)` by `Σ E_i·h_i³/12` (each layer is ONE fiber at its
  centroid: the midpoint rule loses the layer's self-inertia). Compared against
  LadrunoSolidShell — whose `-nz` gauss/lobatto rule integrates z² EXACTLY — this
  reads as "the solid-shell is too stiff". It is the layered quadrature, on both
  the concrete AND the displaced-by-rebar bookkeeping.
- **Workaround/status (measured 2026-07-07, ADR-66 G7):** predict it (the G7
  elastic anchor asserts the layered arm to 1e-4 against the midpoint-rule closed
  form, deficit 2.07% at 5 core layers) or halve the layer thickness (error ∝ h²;
  0.5% at the G7 production layering of 3+10+3).
