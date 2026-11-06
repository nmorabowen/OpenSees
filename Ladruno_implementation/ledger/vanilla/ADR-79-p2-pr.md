---
wp: ADR-79
title: "ADR-79 P2 PR -- 1 vanilla row(s)"
files: ["`SRC/domain/pattern/PathSeries.cpp`"]
table: "main"
legacy_seq: [11]
---
| `SRC/domain/pattern/PathSeries.cpp` | `// Ladruno (ADR 79 P2)`: the interpreter `-time`/`-values` Path route parsed `-useLast` but DROPPED it at the `PathTimeSeries` construction (5th ctor arg defaulted false) — a Path-driven sp snapped to factor 0 when the accumulated domain time overshot the last path point by ulps at the FINAL analysis step (found by the P2 undrained-compression gate; also silently affected the ADR-78 corot gate-3 test, whose final-step snap-back was invisible to its rigid/p≈0 assertions). One-argument forward fix. See [[LEDGER_quirks]]. | ADR-79 P2 PR |
