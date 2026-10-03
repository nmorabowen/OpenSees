---
wp: ADR-84
title: "ASDPlasticMaterial3D's Backward_Euler ACCEPTS a non-converged return map — silently committing f > 0 as success (FIXED opt-in, ADR-84 P2a)"
legacy_seq: 301
---
### ASDPlasticMaterial3D's `Backward_Euler` ACCEPTS a non-converged return map — silently committing f > 0 as success (FIXED opt-in, ADR-84 P2a)

- **Bites:** every ASDP material on the default `Backward_Euler` integrator
  (VonMises, DruckerPrager, MohrCoulomb, the new MohrCoulombTensionCutoff, ...).
  The scalar-Newton consistency loop is written
  `for (int iter = 0; iter < max_iter; ++iter) { ... if (|Phi| < tol_yf) break; ... }`
  and **falls out of `max_iter` with no convergence check at all**, dropping
  straight through to `ComputeTangentStiffness(); return 0;`. A stalled or
  slowly-converging Gauss point therefore reports SUCCESS to the element, the
  element reports success to the algorithm, and the global Newton converges on
  a residual assembled from an inadmissible stress. Nothing anywhere in the
  output says a return map failed — `n_max_iterations` is not a budget, it is a
  silent truncation. This is the persistence mechanism behind the Cerro Lindo
  ADR-0005 M3 finding: 20 Gauss points sitting measurably OUTSIDE the yield
  surface (`f/(2c·cosφ) = +0.0299`) in a model whose every analysis step
  "converged".
- **Why:** convergence is signalled only by `break`, and C++ gives you no way to
  distinguish "broke out early" from "ran out of iterations" without a flag —
  so an author who forgets the flag gets the accepting behaviour by DEFAULT.
  The neighbouring paths are not written this way: `Modified_Euler_Error_Control`
  has an explicit `if (niter > max_iterations) { ...; return -1; }` and
  `Backward_Euler_LineSearch` tracks a `newton_ok` flag and returns -1 once
  substepping is exhausted. Only the plain BE — the DEFAULT integrator — accepts.
- **Workaround/status (2026-08-13, ADR-84 P2a, PR):** fixed **opt-in** via a new
  integration option `strict_convergence` (int, default 0; parsed in the
  `Begin_Integration_Options` block of `OPS_AllASDPlasticMaterial3Ds.cpp`, stored
  in the per-tag static map `INT_OPT_strict_convergence`). With
  `strict_convergence 1`, loop exhaustion with `|Phi| >= tol_yf` prints an
  `opserr` line naming the material tag, the final `|Phi|` and the tolerance,
  and returns -1 so the element reports the failure upward and the algorithm can
  cut back. Default 0 is byte-identical to upstream — deliberately, because
  turning it on changes convergence behaviour for every existing ASDP user.
  **If you are chasing "f > 0 at committed states" in an ASDP model, set
  `strict_convergence 1` before you suspect anything else**; a clean run under
  the flag rules this defect out in one shot.
