---
wp: LEGACY
title: "Isotropic sqrt(A) getCharacteristicLength regularizes a mesh-objectivity band ONLY on in-plane-SQUARE elements (dx == dy)"
legacy_seq: 150
---
### Isotropic sqrt(A) `getCharacteristicLength` regularizes a mesh-objectivity band ONLY on in-plane-SQUARE elements (dx == dy)
- **Bites:** any crack-band/energy mesh-objectivity study (or user model) with `LadrunoSolidShell` (and any element whose lch is the isotropic sqrt of the midsurface/element area) meshed with in-plane rectangles. The crack band localizes in ONE element column, so the physical band width is the element size ALONG the band normal (dx) — but the material regularizes with lch = sqrt(dx*dy). For dx != dy the dissipated energy is off by sqrt(dy/dx), and a "refinement" that changes the aspect ratio reads as spurious energy drift even with `-autoRegularization` working perfectly.
- **Why:** the scalar lch has no direction; sqrt(A) == dx only when dx == dy. The through-thickness projection is already excluded by design (ADR 66 D6), but the in-plane anisotropy is not.
- **Workaround/status (2026-07-06, ADR 66 P5.2 G5):** keep localization-band meshes in-plane square (the G5 gate enforces dx == dy at all three densities: spread 0.9% across a 4x size range, vs the fixed-lch control at ~4x energy error). A directional `lch(n)` API is the ADR 66 O4 backlog item (shared with ADR 19's sqrt(2)-strut residual).
