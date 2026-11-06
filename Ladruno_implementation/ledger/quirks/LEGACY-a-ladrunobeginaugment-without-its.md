---
wp: LEGACY
title: "A ladrunoBeginAugment without its ladrunoEndAugment does not fail — it silently VOIDS every later recorder sample"
legacy_seq: 442
---
### A `ladrunoBeginAugment` without its `ladrunoEndAugment` does not fail — it silently VOIDS every later recorder sample
- **Bites:** the ADR-41 D1 held-load augmentation sweep suppresses recorders and the commitTag
  bump so its passes leave no trace in the output stream (`Domain::commit()`'s
  `if (!contactAugmenting)` guard). That is correct while the sweep is open — and catastrophic
  if it is never closed: the next ordinary step returns `ok = 0`, the time advances `1 → 2`
  exactly as expected, and the recorder file is **empty**. Every downstream check (peak
  displacement, time history, energy balance) reads a truncated or empty stream, and nothing in
  the run says why.
- **Why:** the flag lives on the `Domain` (`Domain::contactAugmenting`), not on the analysis or
  a scope guard, and both commands were bare idempotent setters.
- **Workaround/status (2026-09-14, WP-101):** a second `ladrunoBeginAugment` without an
  intervening `End` now **warns** (it is the signature of exactly this mistake), and the flag is
  cleared by `Domain::clearAll()` (`wipe`) and by `OpenSeesCommands::wipeAnalysis()`. Inside one
  analysis nothing else will catch it, so write the sweep as `try: … finally:
  ops.ladrunoEndAugment()`. Marked `// Ladruno` in `OpenSeesOutputCommands.cpp`,
  `OpenSeesCommands.cpp` and `Domain.cpp` — see [[LEDGER_vanilla_files]].
