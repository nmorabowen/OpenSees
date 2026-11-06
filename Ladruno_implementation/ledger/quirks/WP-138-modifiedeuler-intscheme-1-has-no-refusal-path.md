---
wp: WP-138
title: "ModifiedEuler (IntScheme 1) has no refusal path — at a wall it force-accepts at dt_min, commits ρ_α ≫ 1, and its q–s bends spuriously UP (WP-138)"
legacy_seq: 512
---
### ModifiedEuler (IntScheme 1) has no refusal path — at a wall it force-accepts at dt_min, commits ρ_α ≫ 1, and its q–s bends spuriously UP (WP-138)
- **Bites:** reading a ModifiedEuler SANISAND load–settlement curve near its wall (TIMs' own curves included) as if every committed state were admissible. On the WP-138 footing, E_A forced 3 786 acceptances at dt_min over s/B 0.026–0.0293. It committed ρ_α up to 13.09, where SAS-ME had max 1.004 over the same window. Its q ended +6.4 % above SAS-ME at the same s/B, with an end slope of 1.37 × the initial slope. The counters (`forcedAtDTmin`, `capHits` in `substepStats`) are the only sign of it; the run itself "converges".
- **Workaround:** use IntScheme 129 (SAS-ME), which refuses instead of committing. With ModifiedEuler, gate every reported point on `forcedAtDTmin == 0` for that step and on max ρ_α ≤ 1 + tol. *2026-09-28.*
