---
wp: ADR-92
title: "implexRefusals is a process-wide singleton: read it as DELTAS, never as per-arm totals (ADR 92/93, 2026-09-07)"
date: 2026-09-07
legacy_seq: 385
---
### `implexRefusals` is a process-wide singleton: read it as DELTAS, never as per-arm totals (ADR 92/93, 2026-09-07)
- **Symptom:** running two IMPL-EX arms in one Python process (a probe script, a pytest module) shows the second arm's `implexRefusals` starting from the first arm's totals; `ops.wipe()` does not reset it.
- **Cause:** the four refusal counters live in a process-wide static (the same design as `avgImplexError`'s accumulator), by intent -- the campaign driver reads one element's response as the leg's census.
- **Status (2026-09-15, WP-104): `wipe` now ZEROES the singleton.** The rule below still holds for arms of ONE model (`ops.reset()` / `revertToStart()` deliberately do NOT reset the ledger), but two MODELS in one process no longer share totals -- see the WP-104 entry at the end of this file for the measured symptom (`[9,0,0,9,0,9]` on a fresh material) and what exactly resets.
- **Rule:** snapshot before, subtract after; one arm per process for a reported census. Also: a zero-increment trial (`LoadControl 0.0` hold) sets `f = 0` and measures NO `implexError` -- an `implexDetail[0] == 0` on a hold is by construction, not evidence of a consistent committed state. And on a FREE-DOF deck a hold is not a zero-strain step: Newton moves the nodes to close the committed state's equilibrium gap (the committed stress is the companion's, equilibrium was found on `sigma~`), which is IMPL-EX's defining property; a true zero-increment material probe needs the zero-free-DOF deck (`sani._build`).
