---
wp: WP-105
title: "IntScheme 2 (BackwardEuler_CPPM) on LadrunoSANISAND/ManzariDafalias is QUALIFIED per prescribed increment and REFUTED as a BVP's primary integrator — WP-105 (F…"
legacy_seq: 461
---
### `IntScheme 2` (BackwardEuler_CPPM) on `LadrunoSANISAND`/`ManzariDafalias` is QUALIFIED per prescribed increment and REFUTED as a BVP's primary integrator — WP-105 (F12)
- **Bites:** scheme 2 looks like a strictly-better companion to scheme 1 (`ModifiedEuler`) if you
  only run it at a given strain increment — and it IS, there. On a replayed strain path (zero free
  DOF, the increment supplied rather than proposed) it integrates the **same model to the same
  limit** as scheme 1: `1.3e-3` / `2.9e-3` maximum relative stress deviation over the path at
  `dEz = 1e-5` (`p0 = 100` / `20 kPa`), terminal `eta` within `4.2e-4` / `2.0e-4`; at the campaign's
  own increment (`dEz = 1e-4`) it is **3.7-4.3x more accurate and 4.2-7.6x cheaper**, and at
  `4.6e-4`, **7-30x more accurate and 10-13x cheaper**, with scheme 1 the one leaving its own
  bounding surface at `p0 = 20 kPa` (`eta/M^b = 1.056`, 160% wrong in stress norm). **Put it under a
  global Newton instead — the increment now PROPOSED, not given — and it inverts.** Free-standing
  drained-triaxial at `dEz >= 1e-4`: scheme 2 stalls in **8 of 8** arms (scheme 1 stalls in 1 of 8);
  loosening the global tolerance `1e-9 -> 1e-7` does not rescue it. Each failing step costs
  **12-134 s** against a 30 ms normal step (up to 4400x) grinding the recursive-halving ladder. On
  the real CP1/ADR-95 bearing leg (`x10z8`, `h1.0_e0.6944`, 1200 s budget each) scheme 2 committed
  11 steps to `s/B = 4e-5` against the scheme-1 baseline's 51 steps to `s/B = 0.019` — **475x
  shallower for the same wall clock** — `ds` pinned at 25x the subdivision floor, 100% of its
  committed steps on the relaxed rung 3.
- **Why:** the failure is on off-path trial iterates a global Newton proposes, not on the solution
  path — the replay arm (same path, same `dEz = 4.6e-4`, increment given) completes 40 of 40 in
  1.7 ms/step on the identical model. A coarse-step trial iterate is a strain increment
  `BackwardEuler_CPPM`'s 19-unknown Newton cannot return; it recurses (see the next two entries)
  through up to 512 half-increments before giving up, and every one of those half-increments is
  another 19-unknown Newton of up to 30 iterations. ADR-92 D3's own stated reason for keeping
  scheme 1 as the default ("58-74% of scheme 2's calls take the low-`p` branch and integrate
  explicitly") does NOT reproduce here: measured **0 of 1820** steps on every replayed
  drained-triaxial path at `p0 = 100` and `20 kPa`, and **0 of 80** on the descent of a `p -> p_min`
  path; the 58-74% figure reproduces (53%, `85/160` steps) only once the point is already pinned at
  `p_min` with a zero deviator, i.e. on steps where nothing is being integrated.
- **Workaround/status (2026-09-16):** D3's *conclusion* survives (scheme 1 + `-maxSubsteps` stays
  the companion default) but its *stated reason* does not — see the amendment in
  `92_ladruno_sanisand_implex_adr.md` and the new subsection in
  `LadrunoSANISAND_implex_guide.md` §9. Use scheme 2 only where the increment is already given (a
  prescribed-strain material-point study); do not make it the primary integrator of a load- or
  displacement-controlled BVP without a cap on the ladder — none was tried. `-implex` was OFF in
  every WP-105 arm, so whether scheme 2 is the better commit-time IMPL-EX companion (which runs at
  `commitState` on an increment nothing proposes off-path — exactly the regime where scheme 2 wins)
  is measured nowhere and remains open. Full tables and scripts:
  `Ladruno_files/testbed/hypo_bearing/adr92_f12/F12_intscheme2_verdict.md`.
- **Status (WP-130, #868):** the ladder is bounded on request: `-cppmHalvings n` (0..9, vanilla 9 =
  up to 2^9 half-increments) with `-cppmOnFail refuse`. On a one-quad free-DOF deck (100 kPa, one
  2 kPa push step the CPPM cannot return) the default grinds 31 global iterations, 224 local-Newton
  failures, 440 half-increments and 4 silent explicit fallbacks in 5.6 s before analyze returns
  -3; `-cppmOnFail refuse -cppmHalvings 0` returns -3 in 8-22 ms (pinned,
  `tests/test_ladruno_sanisand_cppm_newton.py`). Refusing fast ALONE does not carry the bearing
  leg (vanilla tangent: s/B 0.00000-0.00017); with the sign-corrected tangent it does -- see
  "IntScheme 2's TanType-2 tangent is MINUS" below.
