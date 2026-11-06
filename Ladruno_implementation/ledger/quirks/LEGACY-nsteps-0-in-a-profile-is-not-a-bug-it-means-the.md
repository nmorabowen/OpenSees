---
wp: LEGACY
title: "nSteps=0 in a profile is NOT a bug — it means the run had no -perStep"
legacy_seq: 222
---
### `nSteps=0` in a profile is NOT a bug — it means the run had no `-perStep`
- **Bites:** you see `nSteps=0` next to a healthy 15-step rollup, add it to a bug list, and go looking for a counter that was never broken. (Done — an early ADR-75 P1h draft published it as a defect alongside the two real `threads`/`nElem` bugs, and it had to be retracted.)
- **Why:** `nSteps` is derived from the per-step **series**, exactly like `dt_min`/`dt_max` — `buildMeta()` only fills it under `if (config_.perStep)`. A coarse run has no series to count, so 0 is the correct answer. The rollup's `root/step` scope still carries the true `calls` count if you need it.
- **Workaround/status:** by design. Read step count from `root/step` `calls`, or run with `-perStep`. Distinguish this from the genuinely-broken attributes in the row above, which were wrong *regardless* of flags. *2026-07-27 (ADR-75 P1h/P1i).*
