---
wp: ADR-93
title: "-implexControl under DisplacementControl walls on trial strains that do not scale with the step (ADR 93, 2026-09-07)"
date: 2026-09-07
legacy_seq: 384
---
### `-implexControl` under `DisplacementControl` walls on trial strains that do not scale with the step (ADR 93, 2026-09-07)
- **Symptom:** every IMPL-EX leg driven by `DisplacementControl` walls at s/B ~0.001–0.002 on control refusals at iteration 1, with the committed `implexError` at the free-surface ring RISING (×30) while the step SHRINKS (÷16); cap hits zero, D2 zero.
- **Cause (inferred from the refusal patterns and error scalings, not measured per iterate):** `DisplacementControl` predicts the load factor from the frozen elastic tangent, which lands an O(1) trial strain on a near-zero-stiffness ring regardless of ds; the control correctly refuses a trial whose error no step size cures.
- **Fix:** use the fork's push idiom — `LoadControl(-ds)` on a prescribed-settlement `sp` (Linear series, Transformation) — under which the error is O(ds), the registered arm walks at full step, and the curve overlays the implicit twin (Esmeralda jobs 146451/146453 vs 146438). The guide now says so. Recorded so nobody re-derives it from a "tolerance too tight" reading, which is what it looks like.
