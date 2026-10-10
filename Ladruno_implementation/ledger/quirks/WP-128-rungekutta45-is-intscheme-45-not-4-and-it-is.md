---
wp: WP-128
title: "RungeKutta45 is IntScheme 45, not 4 — and it is not a usable reference integrator: dT_min 1e-3 hard-coded, force-accept with the Mc clamp, no Stress_Correction…"
legacy_seq: 498
---
### `RungeKutta45` is IntScheme **45**, not 4 — and it is not a usable reference integrator: `dT_min` 1e-3 hard-coded, force-accept with the Mc clamp, no `Stress_Correction` (WP-128)
- **Bites:** `#define INT_RungeKutta45 45`. IntScheme 4 is `INT_MAXENE_FE`, so "IntScheme 4 tight" silently runs a different scheme. RK45's `dT_min = 1e-3` (`TolE = mTolR` is honoured) means at most 1000 substeps before the same forced accept ModifiedEuler has: radial clamp to `Mc` and α re-derived. At the TIMs worst point it "returns" η 1.33 from η 12.87 on a `1e-7` increment. Its `Stress_Correction` calls are commented out, so it committed `f > 1e-6` on 512/640 ring replays at TolR 1e-4. On one vertical-unload cycle chain it ran away to α/α^b 9541. It is also not counted by `substepStats` (ModifiedEuler only).
- **Workaround/status:** use ModifiedEuler at a tight `-honorTolR 1` TolR as a reference and exclude cases where it force-accepted. Remember it is α-blind too (row above). WP-128 used the α-aware port as the check.
