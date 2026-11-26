---
wp: WP-260
title: "HHTGeneralizedExplicit adds damping matrix C to the LHS tangent"
date: 2026-11-25
---
### `HHTGeneralizedExplicit` adds damping matrix `C` to the LHS tangent

- **Bites:** `HHTGeneralizedExplicit` (`SRC/analysis/integrator/HHTGeneralizedExplicit.cpp:238-257`) looks like a single-step dissipative explicit integrator, so pairing it with `system Diagonal` when Rayleigh stiffness-proportional damping (`betaK`) is active breaks the diagonal-mass assumption and forces a sparse linear solve (or fails diagonal assumptions).
- **Why:** `HHTGeneralizedExplicit::formEleTangent()` explicitly calls `theEle->addCtoTang(alphaF * c2);` alongside `addMtoTang(alphaI * c3);`, injecting off-diagonal `C` entries into the global LHS matrix. It also lacks `CentralDifferenceLadruno`'s mass-tangent caching, `-cfl`, `-commitSolveState`, and SMS hooks.
- **Workaround / status (2026-11-25, WP-260 — disclosed):** do not treat `HHTGeneralizedExplicit` as a pure-$M$ explicit wave integrator; build single-step dissipative explicit wave schemes (Hulbert–Chung / Chung–Lee / Soares) on `CentralDifferenceLadruno`'s cached-diagonal-mass chassis.
