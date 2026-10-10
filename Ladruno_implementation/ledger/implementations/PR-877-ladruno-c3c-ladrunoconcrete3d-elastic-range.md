---
wp: PR-877
title: "Ladruno C3c — LadrunoConcrete3D elastic-range cost: skip the ∂ω/∂ε block when no ω is interior"
pr: "#877"
date: 2026-10-10
status: "shipped on branch, unmerged"
section: "table"
---
| **Ladruno C3c — LadrunoConcrete3D elastic-range cost: skip the `∂ω/∂ε` block when no ω is interior** ([[LadrunoConcrete3D_guide]] §13) — `damagedTangent` computed three micro-FD scalar gradients (36 eigendecompositions per Gauss point) plus the `‖Δε_p‖` gradient for EVERY point, although they only feed `dω/dε`, which is zero unless `0 < ω < 1`. Guarded by `needDw`; `dwt = dwc = 0` otherwise, same assembly ⇒ bit-identical stress and tangent (dump of 320 steps × std/bbar × implicit/IMPL-EX, every stress + 24×24 stiffness entry equal to the pre-change build; kernel g++ check + LadrunoConcrete3D suites green). No change to the elastic predictor / Lode / eigen path (kept for byte identity; the eigen work was not the cost). **Gate:** per elastic iteration vs ASDConcrete3D on the L&W coarse model 3.6–4.2× → 1.3× (gate ≤ 1.5×); 512-brick block 4.7× → 1.2×. | performance (kernel) | — (33017 unchanged) | `SRC/material/nD/LadrunoConcrete3DKernel.h`, `Ladruno_implementation/LadrunoConcrete3D_guide.md` | shipped on branch, unmerged | [#877](https://github.com/nmorabowen/OpenSees/pull/877) |
