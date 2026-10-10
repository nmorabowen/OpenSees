---
wp: LEGACY
title: "The apex region is in the RELATIVE deviator r = s - alpha; two of the three places that decide it were not (one fixed, one pinned)"
legacy_seq: 445
---
## The apex region is in the RELATIVE deviator `r = s - alpha`; two of the three places that decide it were not (one fixed, one pinned)

Drucker-Prager's surface is written in `r = dev(sigma) - alpha`, so every test
about where a trial sits relative to the vertex has to be written in `r` too.
Three places decide that, and they did not agree:

| | measured in | status |
|---|---|---|
| `DruckerPrager_YF::check_apex_region` (Euclidean) | **`r`** — correct variable, wrong metric | unchanged |
| `cp_apex_region` (elastic metric, the one the integrator uses) | `dev(sigma)` | **FIXED, F8 round 2** |
| `DruckerPrager_YF::apex_stress()` | ignores `alpha` entirely | **pinned, not fixed** |

The consequence, measured with `alpha = (0.02, 0.02, -0.04, 0, 0, 0)`, `Ht = 0`,
associated, trial `(q_tr, p_tr) = (0.02, 0.3810)`: the raw-deviator flip test
classifies it APEX although its exact return is the cone point
`sqrt(J2(r)) = 0.017321`, `p = 0.220190`. `be_apex_project` is then asked for a
vertex the yield function cannot supply — `apex_stress()` returns `p_apex*I`,
whose `|f|` is exactly `sqrt(J2(alpha))` = **0.034641**, not zero — so its own
`|f(sigma_apex)| <= tol_yf` guard fires and the step is **refused** (`rc = -3`,
both `strict_convergence` settings). An admissible return, reported impossible.

This is **older than F8** and is not the union's doing: the Euclidean test also
over-classifies this trial (`p - p_apex = 0.1220 >= eta*q_rel = 0.0244`), so
every arrangement since wp/94c made the apex projection live — Euclidean alone,
Euclidean OR elastic-metric, or elastic-metric alone — reached the same refusal.
What F8 round 2 changed is the one of the three that the integrator now relies
on: `cp_apex_region` measures in `r`, the trial classifies CONE, and the state
returns to the closed-form cone point (verified to 1e-6, `rc = 0`, both strict
settings).

**Still pinned:** `apex_stress()`. The vertex of this surface is at
`sigma = alpha + p_apex*I`, not `p_apex*I`, so a state that genuinely IS in the
apex region with `alpha != 0` still hits the `|f(sigma_apex)|` guard and is
refused rather than projected. Refusing is the safe half of wrong — the guard
exists exactly so a bad vertex is never committed — but the fix is one line in a
vanilla yield function that also changes `Closest_Point`'s apex return, so it
belongs in its own WP with its own gate. Both tests here use zero back stress in
every other case, which is why nothing else moved.

It also falsifies a comment that stood in `ASDPlasticMaterial3D.h`'s
`be_apex_project`: "every yield function that opts into `yf_has_apex` today is
perfectly plastic (Null hardening)". `DruckerPrager_YF` opts in with
`AlphaHardeningType` / `CohesionHardeningType` template parameters and the
registered specializations include both linear tensor and linear scalar
hardening. The comment is corrected in the same PR.

It also falsifies a comment that stood in `ASDPlasticMaterial3D.h`'s
`be_apex_project`: "every yield function that opts into `yf_has_apex` today is
perfectly plastic (Null hardening)". `DruckerPrager_YF` opts in with
`AlphaHardeningType` / `CohesionHardeningType` template parameters and the
registered specializations include both linear tensor and linear scalar
hardening. The comment is corrected in the same PR.
