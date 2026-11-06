---
wp: LEGACY
title: "-fsL zero under a quasi-static implicit fs1 march = the naive drained split -- diverges in ~4 steps at soil coupling (explicit-lane-only setting)"
legacy_seq: 193
---
### `-fsL zero` under a quasi-static implicit fs1 march = the naive drained split -- diverges in ~4 steps at soil coupling (explicit-lane-only setting)
- **Bites:** an overlay built with `-fsL zero` (the ADR-73 P3 explicit-lane setting) but driven by a quasi-static / implicit fs1 march (plain `analyze` under Newmark at consolidation time scales) diverges within ~4 steps, ~10 orders of magnitude, at realistic soil coupling strength tau = (alpha^2/K_dr)/storage ~ 1e3 (measured, ADR-73 SS3.2) -- loudly, not wrongly.
- **Why:** L = 0 removes the fixed-stress relaxation entirely, so the split IS the naive drained split the ADR bans for the implicit lane. The setting exists ONLY for the explicit lane, where stability is governed by dt <= 0.5x the discrete undrained pencil (E7.2: the L=0 implicit-fluid boundary = exactly 1.000x the pencil) and no iteration/L is needed.
- **Workaround/status (2026-07-18, ADR-73 P3):** by design. The parser prints a one-time loud advisory at `-fsL zero`; `LadrunoStaggeredAnalyze` refuses FSL_ZERO overlays with a loud fatal (iterating with L = 0 is the same drained split). Use `-fsL classic|oedometric` for implicit/driver lanes; reserve `zero` for `CentralDifferenceLadruno`/explicit runs at dt <= 0.5x the (overlay-aware) `criticalTimeStep` pencil.
