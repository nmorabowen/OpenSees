---
wp: LEGACY
title: "Consistent/Olovsson SMS + LadrunoPorousOverlay: correct undrained PRICING but measured UNDER-DELIVERY on the coupling mode (use lumped SMS)"
legacy_seq: 195
---
### Consistent/Olovsson SMS + LadrunoPorousOverlay: correct undrained PRICING but measured UNDER-DELIVERY on the coupling mode (use lumped SMS)
- **Bites:** `CentralDifferenceSMSConsistent` on an overlay model, post-P3b: the sizing report correctly prices the undrained pencil (INFO line prints), yet the certified `dtTarget` march can still diverge — measured uniform ~×1.83/step growth from step 1 at dtTarget = 3× the unscaled pencil on the e72 column (battery gate (d) EXPECTED-LIMITED record), while the LUMPED builder's certified march on the same model is stable 4000+ steps.
- **Why (ADR §12 P3b item 5, panel-checked):** the Olovsson centroid-preserving M̄ blocks add inertia only to the NON-RIGID element modes (that is their design — element mass distribution preserved). The overlay's undrained volumetric coupling mode carries a large rigid-translation component that stays UNSCALED, so the coupled frequency scales by less than √s and the certified step over-promises for overlay-owned cells. This is a scheme interaction, not a wiring bug — the lumped builder injects real nodal mass and delivers.
- **Workaround/status (2026-07-19, ADR-73 P3b):** use lumped `CentralDifferenceSMS` with overlays, or size dt from the overlay-aware `criticalTimeStep()` report. The consistent builder keeps the undrained pricing (honest report) and prints a loud one-time warning when it scales overlay-owned elements ("the Olovsson centroid-preserving blocks under-scale the undrained COUPLING mode").
