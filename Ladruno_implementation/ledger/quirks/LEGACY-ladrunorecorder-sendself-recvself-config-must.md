---
wp: LEGACY
title: "LadrunoRecorder sendSelf/recvSelf config must stay in LOCKSTEP — two fields were missing (FIXED) + the rule"
legacy_seq: 48
---
### LadrunoRecorder `sendSelf`/`recvSelf` config must stay in LOCKSTEP — two fields were missing (FIXED) + the rule
- **Bites:** any recorder config field set by the OPS_ parser must be (a) `ser.put_*`'d in
  `sendSelf` AND (b) `de.get_*`'d in `recvSelf` **in the same order**, or OpenSeesMP worker
  ranks (built via `recvSelf`) silently diverge from P0 (built via the parser). Two fields
  were configured but NOT transmitted: `m_data->envelope_mode` (`-envelope`) and
  `m_data->info.store_data_f32` (`-precision f32`). Consequence: with `-envelope` or
  `-precision f32` under MP, P0 wrote ENVELOPES/f32 while every worker wrote full
  time-series/f64 → `.part-N` files with a *different schema* than `.part-0`, breaking the
  apeGmsh stitch-on-read. **Fix:** append both to `sendSelf` (after the elemental-results
  block) and `recvSelf` (after the same block), same order. **Rule for the next agent: when
  you add ANY field to the recorder config, grep `sendSelf`/`recvSelf` and add it to both.**
  (Found by the 2026-06-03 adversarial review. The MP round-trip itself was verified by
  construction — symmetric put/get — not by a live 2-rank run, since the worktree had only
  the OpenSeesPy build; the `mp_parallel` gate needs OpenSeesPyMP.)
