---
wp: ADR-94
title: "...and the fix for it was a UNION, which is only safe while K*etabar/G <= eta (ADR-94 addendum, F8)"
legacy_seq: 443
---
## ...and the fix for it was a UNION, which is only safe while `K*etabar/G <= eta` (ADR-94 addendum, F8)

wp/94f (above) repaired the Euclidean apex test by **unioning** the yield
function's answer with the exact elastic-metric one inside `Backward_Euler`:

```cpp
bool be_in_apex = yf.check_apex_region(...);      // (p - p_apex) >= eta*q
if (!be_in_apex)
    be_in_apex = cp_apex_region(...);             // (p - p_apex) >= (K*etabar/G)*q
```

A union takes the **wider** of the two regions, and **which one is wider depends
on the flow rule**, because the two slopes are `eta` and `K*etabar/G`:

* `etabar = 0` (the deck wp/94f measured): exact slope 0, region `p >= p_apex`,
  which strictly CONTAINS the Euclidean cone. Union == replace. Correct.
* `etabar > 0`: the exact slope `K*etabar/G` **overtakes** the Euclidean one as
  soon as `etabar > eta*G/K`. That threshold is tiny: on this deck
  `G/K = 0.10345` and `eta = 0.4457`, so it is `etabar = 0.0461`, i.e.
  **psi ~ 2.3 deg** — not "associated", but essentially **every dilatant sand**.
  At `etabar = eta/2` (psi ~ phi/2) the exact slope is 2.1545 against the
  Euclidean 0.4457: a 4.8x-wide wedge. At `etabar = eta` (associated) it is
  4.3089, about **ten times** the Euclidean cone, and the union keeps the
  Euclidean answer at every one of these.

So a dilatant leg of the very same deck apex-projected every trial in the wedge
`eta*q <= p - p_apex < (K*etabar/G)*q`, whose correct return is to the cone
FLANK. Symptoms, all of them misleading:

* **no refusal, no message, no NaN** — the material reports success at every
  Gauss point, because an apex projection always "works";
* the committed stress at those Gauss points is `sigma_apex` with **zero
  deviator** (on the ADR-95 cone that is 0.259 kPa hydrostatic in a 200 kPa
  field), and under `tangent_type Continuum` the Gauss point also reports a
  **zero tangent** — by design, since the apex stress cannot move;
* the analysis therefore presents as "the element/mesh walls while still
  hardening": the step size collapses, the load is still climbing, and every
  refusal counter reads zero. Measured on the ADR-95 R3 gate deck at `h0 = 1.0`:
  ASD associated crawled to s/B 0.016 in 773 s where the vanilla UW material
  reached s/B 0.15 in 63 s with zero failed attempts.

Rules that generalize:

4. **Never union two classifications of the same region unless you have shown
   which one is wider, in every regime the code will see.** Two tests of "am I
   in region X" are not interchangeable with "am I in region X *or* Y". If one
   of them is the exact test, it should REPLACE the other, not widen it — that
   is what an opt-in trait like `yf_apex_elastic_metric` is for.
5. **A silent wrong answer is harder to find than a loud refusal, and fixing a
   refusal can create one.** wp/94f converted a refusal (435 `scalar Newton
   exhausted` lines) into silence; the leg that then walled had *zero*
   diagnostics. When a fix removes a failure mode, ask what it turns that
   failure into on the decks it was not measured on.
6. **A cone's apex test is flow-rule-dependent; its yield surface is not.** Two
   legs that differ only in `DP_etabar` exercise genuinely different branches of
   the return map, so a zero-dilatancy acceptance leg is not evidence about a
   dilatant one. Run both — and note how low the crossover is (psi ~ 2.3 deg
   here): "we only tested psi = 0" and "we tested the usual case" are not the
   same sentence.
