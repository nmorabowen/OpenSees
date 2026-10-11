---
wp: WP-183
title: "Constraint FEs carry no mass: getM_Force/getC_Force return zero, eigen mass loops skip element-less FEs"
date: 2026-10-11
files: ["`SRC/analysis/fe_ele/penalty/PenaltyMP_FE.cpp`", "`SRC/analysis/fe_ele/penalty/PenaltyEQ_FE.cpp`", "`SRC/analysis/fe_ele/lagrange/LagrangeMP_FE.cpp`", "`SRC/analysis/fe_ele/lagrange/LagrangeSP_FE.cpp`", "`SRC/analysis/fe_ele/lagrange/LagrangeEQ_FE.cpp`", "`SRC/analysis/integrator/EigenIntegrator.cpp`", "`SRC/analysis/analysis/StaticAnalysis.cpp`", "`SRC/analysis/analysis/DirectIntegrationAnalysis.cpp`", "`SRC/analysis/analysis/StaticDomainDecompositionAnalysis.cpp`", "`SRC/analysis/analysis/TransientDomainDecompositionAnalysis.cpp`"]
table: "upstreamable"
---
| `SRC/analysis/fe_ele/penalty/PenaltyMP_FE.cpp` | `// Ladruno WP-183`: `getM_Force` / `getC_Force` zero `resid` and return it. Upstream returned `*resid`, the last constraint residual `-alpha C^T C (u - u0)`, which `ArpackSolver::myMv` (non-diagonal M) and `IncrementalIntegrator::doMv` (modalDamping setup) summed into `M x`: after a Newton or multi-step Linear analysis the eigenvalues came out ~ +-1e-3 and a modalDamping transient was off by ~50 %. Identical on OpenSees/OpenSees master; no upstream issue filed. | WP-183 |
| `SRC/analysis/fe_ele/penalty/PenaltyEQ_FE.cpp` | `// Ladruno WP-183`: same as PenaltyMP_FE (and the "not yet implemented" warnings go). | WP-183 |
| `SRC/analysis/fe_ele/lagrange/LagrangeMP_FE.cpp` | `// Ladruno WP-183`: same; the multiplier rows have no mass either (`LagrangeDOF_Group` already returned zero). | WP-183 |
| `SRC/analysis/fe_ele/lagrange/LagrangeSP_FE.cpp` | `// Ladruno WP-183`: same (its warning text also named PenaltySP_FE). | WP-183 |
| `SRC/analysis/fe_ele/lagrange/LagrangeEQ_FE.cpp` | `// Ladruno WP-183`: same. | WP-183 |
| `SRC/analysis/integrator/EigenIntegrator.cpp` | `// Ladruno WP-183`: `formM` skips an FE with `getElement() == 0` (Penalty/Lagrange MP/SP/EQ, LadrunoContactFE). Those FEs' `getTangent` ignores the integrator and returns the constraint stiffness (`alpha C^T C`, `alpha`), which upstream added to M: M non-diagonal (ARPACK then takes the `getM_Force` path above), and explicit-M solvers (`-fullGenLapack`, `-feast`) report lambda ~ 1 at every constrained dof even at K0. | WP-183 |
| `SRC/analysis/analysis/StaticAnalysis.cpp` | `// Ladruno WP-183`: same skip in the `eigen()` mass loop. This is the path classic-Tcl `eigen` takes once a static analysis exists (e.g. after gravity). | WP-183 |
| `SRC/analysis/analysis/DirectIntegrationAnalysis.cpp` | `// Ladruno WP-183`: same skip in the `eigen()` mass loop (transient analysis present, and the ephemeral analysis openseespy builds). | WP-183 |
| `SRC/analysis/analysis/StaticDomainDecompositionAnalysis.cpp` | `// Ladruno WP-183`: same skip (OpenSeesSP/MP subdomains). Built, not exercised by a test. | WP-183 |
| `SRC/analysis/analysis/TransientDomainDecompositionAnalysis.cpp` | `// Ladruno WP-183`: same skip (OpenSeesSP/MP subdomains). Built, not exercised by a test. | WP-183 |
