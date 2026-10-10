---
wp: ADR-94
title: "The ASD Drucker-Prager apex region test was EUCLIDEAN; zero-dilatancy decks walled at the footing edge (ADR-94 wp/94f)"
legacy_seq: 428
---
## The ASD Drucker-Prager apex region test was EUCLIDEAN; zero-dilatancy decks walled at the footing edge (ADR-94 wp/94f)

`DruckerPrager_YF::CHECK_APEX_REGION` classifies a trial state as "apex" when
`(p - p_apex) >= eta*q` — the normal-cone test drawn with a Euclidean protractor
in the `(p, sqrt(J2))` half-plane. The exact test lives in the ELASTIC metric:

```
apex  <=>  (p - p_apex) >= (K * etabar / G) * q
```

The two coincide only when `K*etabar/G == eta`. On a NON-ASSOCIATED deck with
zero dilatancy (`DP_etabar = 0`, the standard geotechnical choice, and the one
the ADR-95 Prandtl-Reissner footing uses) the exact test degenerates to
`p >= p_apex`, because a non-dilatant flank return has a traceless flow
direction and **cannot move the mean stress at all**. So the Euclidean test is
strictly too NARROW: every over-apex trial with `q > eta*(p - p_apex)` was sent
to the flank map, where `f = q + eta*p - xi_c >= eta*p - xi_c > 0` for every
`q >= 0` — no solution exists. The scalar Newton then exhausted its iterations
and, under `strict_convergence`, refused the step.

What made this expensive to find: **nothing printed the word apex**. The
`|f(sigma_apex)| <= tol_yf` guard #815 added only catches an apex that is not on
its own surface; a *misclassified* state never reaches the guard. The visible
symptom was 435 `scalar Newton exhausted` refusals and a footing that stopped at
s/B 0.011 with a perfectly healthy global tangent — which reads as a mesh or an
element problem, not a return-map problem.

Rules that generalize:

1. **A region test that needs the elastic tangent cannot live in the yield
   function.** The YF signature sees stress, internal variables and parameters —
   not `K`, not `G`, and not the plastic potential's dilatancy. Any "which
   branch of the return map" decision that depends on the metric belongs in the
   integrator. ADR-97's `Closest_Point` already knew this (`cp_apex_region`);
   `Backward_Euler` did not, and the two integrators disagreed about the same
   material for two PRs.
2. **Order the returns so the cheap map cannot be the last word.** "Flank first,
   apex if the flank fails and the trial is past the vertex" is robust to a
   wrong classification in a way that "classify, then commit" is not.
3. When a return map exhausts, ask whether the state it was handed has a
   solution at all before tuning `n_max_iterations` or the tolerance.
