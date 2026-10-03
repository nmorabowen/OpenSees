---
wp: ADR-43
title: "ADR43 P1, P2 -- 1 vanilla row(s)"
files: ["`SRC/interpreter/OpenSeesCommands.cpp`"]
table: "main"
legacy_seq: [259]
---
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno` ADR43 P1: (1) `eigen -feast fmin fmax [-m0][-nq][-tol][-maxiter][-verbose]` — `OPS_eigenFeast()` dedicated flow (a frequency band instead of a mode count) building a configured `FeastEigenSOE`/`FeastEigenSolver` pair passed through the existing `providedEigenSOE` seam; (2) the found-mode reconcile in `OpenSeesCommands::eigen` — the band defines the count, so post-solve the domain eigenvalue list is trimmed to the solver's found m (the analysis loop runs with a cap; beyond-found getters return silent zeros). +2 includes. **P2:** `-certify` flag (Sturm/inertia completeness certificate — PARDISO negative-pivot counts at the band edges, refuse-on-mismatch AND refuse-on-perturbed-pivots per the adversarial gate). | ADR43 P1, P2 |
