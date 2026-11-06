---
wp: LEGACY
title: "timeSeries Path -time/-values -useLast silently dropped the flag — Path-driven sp snapped to zero at the FINAL analysis step"
legacy_seq: 249
---
### `timeSeries Path -time/-values -useLast` silently dropped the flag — Path-driven sp snapped to zero at the FINAL analysis step
- **Bites:** any Python/interpreter model driving a `Plain` pattern's `sp` with `timeSeries("Path", tag, "-time", ..., "-values", ..., "-useLast")`. The interpreter parse (`OPS_PathSeries`, `SRC/domain/pattern/PathSeries.cpp`) READ `-useLast` but then constructed the `-time` variant as `PathTimeSeries(tag, path, time, factor)` — dropping the flag (5th ctor arg defaults false). Because the DOMAIN time accumulates dt in floating point, the final step's pseudo-time overshoots the last path point by a few ulps, `getFactor` takes the beyond-the-end branch, and with `useLast == false` returns **factor 0**: the constrained DOFs snap back to zero exactly at the last step.
- **Tell:** a monotone Path-driven quantity that grows correctly for N−1 steps and collapses at step N (the ADR-79 P2 undrained gate saw p: 6.16e5 → 5.2e4 at the final step, top displacement 0.0975 → 0.0000). Invisible when the final state happens to satisfy the assertion anyway — the ADR-78 corot gate-3 test asserts p ≈ 0 under rigid rotation, so its final-step snap-back to u = 0 ALSO read p ≈ 0 and passed.
- **Workaround/status:** ✅ FIXED (ADR-79 P2 PR) — the parsed `useLast` is forwarded to the `PathTimeSeries` ctor (`// Ladruno` mark; [[LEDGER_vanilla_files]] row). Belt-and-braces for test authors: give the Path an extra terminal point beyond the last analysis time so no route depends on the beyond-the-end branch. *2026-07-28 (ADR-79 P2).*
