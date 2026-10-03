---
wp: ADR-80
title: "An integrator remedy's measured gain is a property of the ITERATION BUDGET, not of the idiom — ADR-80 P3's x18.7 does not transfer to a deck whose ladder opens…"
legacy_seq: 362
---
## An integrator remedy's measured gain is a property of the ITERATION BUDGET, not of the idiom — ADR-80 P3's x18.7 does not transfer to a deck whose ladder opens at 25 iterations

**Measured 2026-09-05, ADR-92 CP1 (`_adr92_cp1_surcharge_results` §7).**

`LadrunoLoadControl -tangentPredictor` (ADR-80 P3, #786) fixes a real defect: stock
`LoadControl` has nothing in front of `applyLoadDomain`, so with a non-homogeneous `sp` the
driven layer's first constitutive evaluation of each increment is overstrained by `L/h`. On
ADR-80's own gate that cost **23 cutbacks and x18.7 iterations**, which the flag took to
**0 and x1.00**.

The GATE U footing deck is the same idiom — prescribed settlement by `sp` + `LoadControl`,
adaptive march with cutbacks, and a maximally path-dependent material (`LadrunoSANISAND`).
A one-flag A/B on it measured **+7.6 % in settlement rate, the load unchanged to < 0.3 %,
and MORE ladder failures (37 -> 44)**, with both legs spending 2 of 80 subdivisions.

- **Why:** ADR-80's gate ran `maxIter = 5`, where the artificial residual is enough to FAIL
  the step and force a cutback — the cost was cutbacks, not arithmetic. This deck's ladder
  opens with `Newton` at **25** iterations on a graded mesh whose driven layer is one element
  deep, so the same residual is absorbed in roughly one extra iteration and never reaches the
  cutback path.
- **Rule:** before adopting a solver remedy on the strength of another lane's headline
  number, check whether *your* march can even reach the failure mode that number was measured
  on. A remedy whose benefit is "avoided cutbacks" is worth nothing to a march that was not
  cutting back. Keep the flag anyway when it is free and non-negative — but do not budget for
  the other lane's factor.
- **Not a defect in the flag.** It engaged (banner, `predictor: true` in the leg payload, no
  integrator warning) and it did not move the answer, which is what a correct predictor
  should do.
