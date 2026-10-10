---
wp: LEGACY
title: "The P1 settlement-column test deck hits the companion cap at 20000 in 19 of 30 steps at ds 1e-3 (2026-09-07)"
date: 2026-09-07
legacy_seq: 386
---
### The P1 settlement-column test deck hits the companion cap at 20000 in 19 of 30 steps at ds 1e-3 (2026-09-07)
- **Symptom:** `_build_settlement_column` (tests/test_ladruno_sanisand_implex.py) pushed at ds = 1e-3 under `-implex` records `implexRefusals[3]` (companion cap at commit) on 19/30 steps even at `-maxSubsteps 20000`, silently force-accepted (`Domain::commit()` discards the refusal), and the implicit twin cannot converge that deck at ds 1e-3 or 1e-4 at all.
- **Consequence:** that column is a correctness deck for the operator's bookkeeping, NOT a deck on which `_CAP_ADEQUATE = 20000` is adequate; any test that reads a *history* off it must assert `implexRefusals[3]` unchanged first. The campaign deck (`adr92_bvp_fix/`) had zero cap hits at the same cap.
