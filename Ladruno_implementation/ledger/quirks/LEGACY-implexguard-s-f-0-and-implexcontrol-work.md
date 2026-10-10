---
wp: LEGACY
title: "-implexGuard's f = 0 and -implexControl work against each other — the guard makes the prediction inaccurate and the control then refuses the step for being ina…"
legacy_seq: 436
---
### `-implexGuard`'s `f = 0` and `-implexControl` work against each other — the guard makes the prediction inaccurate and the control then refuses the step for being inaccurate
- **Bites:** the P2-2 guard (`:2359`) forces `f = 0` on any step whose committed
  predecessor showed a loading reversal or `Kp <= 0`. The guide says this trades
  "the prediction's accuracy, not the step". That is true with `-implexControl`
  **off**. With it on, `sigma~` is then a pure elastic predictor, `implexError`
  measures the whole plastic correction, and the step is refused — and P2-6's
  trial-time fallback cannot help, because it only runs when
  `mImplexFactor != 0.0` (`:3055`).
- **Measured:** **30 of leg B's 49 throttled refusal lines report `f = 0`**, and
  so do all four of the over-tolerance Gauss points in the fixed-`ds` census.
  `-implexGuard off` buys +85 % reach (0.0085 -> 0.0157); both guards off +117 %
  (0.0185) at 12x the refusal count. `implexGuards[1]` ran to 31 491 on leg B,
  about 3.8 % of all Gauss-point-steps.
- **Workaround/status (2026-09-14):** do NOT reach for `-implexGuard off` as a
  fix — the guard exists for ADR-93's softening seat, which that deck does not
  test. Turn the control off instead. The design question — whether a guarded
  point should be exempt from the control's tolerance for that step, the same
  shape as the un-primed exemption — is recorded in
  [[92b_implex_selfweight_wall_note]] section 10 for ADR 92 to settle.
