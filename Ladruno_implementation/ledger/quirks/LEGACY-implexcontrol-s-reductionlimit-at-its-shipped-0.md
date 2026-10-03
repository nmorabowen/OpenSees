---
wp: LEGACY
title: "-implexControl's reductionLimit at its shipped 0.01 puts the material's floor exactly ON the harness's own DS_MIN — it can never fire, and raising it is \"turn…"
legacy_seq: 437
---
### `-implexControl`'s `reductionLimit` at its shipped `0.01` puts the material's floor exactly ON the harness's own `DS_MIN` — it can never fire, and raising it is "turn the control off after one halving"
- **Bites:** `implexGuards[0]` (floor fallbacks) reads 0 and the P2-1 floor
  policy looks dead. Raising `reductionLimit` then looks like a free
  tolerance-preserving win (+71 % reach measured).
- **Why, arithmetically:** the floor is `reductionLimit * |dt0|`, and on the
  fork's R3-derived controller `|dt0| = DS_BASE = 2e-5` m while `DS_MIN = 2e-7`
  m. At the shipped `0.01` the material's floor is `0.01 x 2e-5 = 2e-7` — the
  harness's own floor to the digit, so the driver declares `FLOOR` at exactly the
  step where the branch would first become reachable. (The guide section 7
  records `reductionLimit` "measured inert"; this is why.)
- **What raising it actually does:** at `0.5` the floor sits at `1e-5`, one
  halving below the base step, so the P2-1 floor branch — which DELIVERS the
  companion and does NOT refuse — takes over almost immediately. That is not a
  gentler tolerance; it is the control switching itself off after one halving.
  Measured: 151 floor fallbacks, +71 % reach, `tol` untouched.
- **Workaround/status (2026-09-14):** if you want that behaviour, prefer saying
  so — run with the control off (above) rather than with a floor set so high the
  tolerance branch cannot reach.
