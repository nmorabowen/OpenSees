---
wp: WP-217
title: "IMPL-EX hinges read ops_Dt live: under DisplacementControl the extrapolation moved with the Newton iterate"
date: 2026-10-29
---
### `-implex` hinge or IMK/Bilin spring gives a different answer under DisplacementControl than under LoadControl on the same path, or loses first-order convergence

- **Bites:** `LadrunoCohesiveHinge`, `LadrunoCohesiveHingeBiaxial`, `ModIMKPeakOriented`
  and `Bilin` with `-implex`, driven by `DisplacementControl` (or arc-length).
  frame-ladruno-validation 00a (B3): with Mc = 120, Gf = 1.5, Kpen = 48000 and uniform
  steps h = 1.40625e-4, step 18 (theta = 1.0125 kappa0) returned M = 121.5. That is the
  undamaged Kpen·theta. The published scheme (Prazeres et al. 2016 Eq. 17) and the same
  material under sp + LoadControl(1) give 119.6671. The max error against the converged
  reference went 2.28e-2, 1.53e-2, 1.53e-2, 2.87e-3 Mc for N = 40..320, which is not first order.
- **Why:** the extrapolation factor was `f = ops_Dt / dt_n`, with `ops_Dt` read on every
  trial. In a static analysis `ops_Dt` is the load-factor increment set by
  `Domain::update`. Under DisplacementControl that increment is part of the solution: it
  changes with each Newton iterate, flips sign on softening (then `f` fell back to 1), and
  depends on the material's own tangent. So the "explicit" step was neither explicit nor
  iterate-independent. The `LEGACY` dual-damage quirk already named the trap ("derive the
  IMPL-EX dt from a monotone control parameter instead of λ").
- **Fix (WP-217):** `SRC/material/uniaxial/LadrunoImplexClock.h`, shared by all four materials.
  The default is the `step` clock: pseudo-time advances by one per committed step, so `f = 1`.
  It is exact for uniform steps of any control variable and independent of the integrator.
  `-implexDt pseudo` keeps the old `ops_Dt` ratio for variable LoadControl/transient steps,
  latched at the first trial of the step and re-latched after a `revertToLastCommit`.
  Under DisplacementControl the latch is the predictor's dλ: explicit, but not a meaningful
  clock. `ModIMKPeakOriented` and `Bilin` call `revertToLastCommit()` at the top of
  `setTrialStrain`, so they save and restore the latch around that call; a reopen there
  would re-latch on every iterate.
- **Behaviour change:** decks with `-implex` and NON-uniform LoadControl/transient steps
  now extrapolate with `f = 1` instead of the step ratio. Add `-implexDt pseudo` for the old rule.
- **Not covered:** `LadrunoIMKBeam` latches `ops_Dt` at its first trial (`stepDtLatched`).
  That is explicit, but under DisplacementControl it is the predictor's dλ.
  `LadrunoSANISAND` has its own `-implexDt pseudo|strain|user` (ADR 92).
- **Gate:** `tests/test_ladruno_cohesive_hinge_implex_bk.py` (DisplacementControl ==
  sp + LoadControl for all four materials; the step-18 repro; the pseudo ratio).
