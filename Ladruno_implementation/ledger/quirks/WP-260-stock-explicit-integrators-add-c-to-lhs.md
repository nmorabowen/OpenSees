---
wp: WP-260
title: "Stock explicit integrators add the damping matrix C to the LHS tangent"
date: 2026-11-25
---
### Stock explicit integrators (`CentralDifference`, `HHTExplicit`, `HHTGeneralizedExplicit`) add `C` to the LHS tangent

- **Bites:** with Rayleigh stiffness-proportional damping (`betaK`, or any element `C` with off-diagonal terms) the LHS of a stock explicit integrator is no longer diagonal. Under `system Diagonal` the off-diagonal part of `C` is then **silently dropped** (`DiagonalSOE::addA` keeps `m(i,i)` only, or row-sums with `-lump`): the step runs and reports numbers with the wrong damping. Under a sparse SOE every step pays a full factorization, so the integrator is "explicit" in name only.
- **Why:** `formEleTangent`/`formNodTangent` call `addCtoTang` next to `addMtoTang`: `CentralDifference.cpp:155-156,166-167` (`c2`), `HHTExplicit.cpp:237-238,248-249` (`alpha*c2`), `HHTGeneralizedExplicit.cpp:238-257` (`alphaF*c2`). `CentralDifferenceLadruno` adds only `M` (`CentralDifferenceLadruno.cpp:285-297`) and treats damping on the right-hand side; the stock family also lacks its mass-tangent cache, `-cfl`, `-commitSolveState` and SMS hooks. Same silent-drop mechanism as the quirk "`system Diagonal` silently drops the off-diagonal of a condensed mass".
- **Workaround / status (2026-11-25, WP-260 — disclosed):** for diagonal-mass explicit wave runs with damping use `CentralDifferenceLadruno`; build single-step dissipative explicit schemes (Hulbert–Chung / Chung–Lee / Soares) on its cached-diagonal-mass chassis, never on the stock HHT explicit family.
