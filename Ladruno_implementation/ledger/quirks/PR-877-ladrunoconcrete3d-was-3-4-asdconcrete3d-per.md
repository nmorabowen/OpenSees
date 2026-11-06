---
wp: PR-877
title: "LadrunoConcrete3D was 3–4× ASDConcrete3D per ELASTIC iteration — the damaged tangent finite-differenced damage gradients that are zero for every undamaged point"
pr: "#877"
date: 2026-10-10
---
### LadrunoConcrete3D was 3–4× ASDConcrete3D per ELASTIC iteration — the damaged tangent finite-differenced damage gradients that are zero for every undamaged point
- **Bites:** L&W beam coarse, 864 bbar bricks: 0.71–0.84 s vs 0.20 s per Newton iteration in the elastic range (README D2 saw 2.6–3.1 vs 0.6–0.8 s on a loaded machine). Looks like the 3⁶ spectral-tangent loops or the eigen solves.
- **Why:** `damagedTangent` always evaluated `dscalarDsig` three times (6-component central FD of `eig3sym` + `damageDrivers` = 36 eigendecompositions) and the plastic-strain gradient, although they feed only `∂ω/∂ε`, which is identically zero unless `0 < ω < 1`.
- **Fix (2026-09-26, Ladruno C3c):** skip that block when no ω is interior; bit-identical; 1.3× ASD on L&W. Profile before optimizing the "obvious" loops.
