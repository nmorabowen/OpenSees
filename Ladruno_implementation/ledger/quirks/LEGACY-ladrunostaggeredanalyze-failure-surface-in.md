---
wp: LEGACY
title: "LadrunoStaggeredAnalyze failure surface in openseespy is split: parse/analysis-setup fatals raise, run-time aborts return a negative int"
legacy_seq: 189
---
### `LadrunoStaggeredAnalyze` failure surface in openseespy is split: parse/analysis-setup fatals raise, run-time aborts return a negative int
- **Bites:** `ops.LadrunoStaggeredAnalyze(...)` raises `OpenSeesError` for a static-analysis-active or no-transient-analysis misuse (detected in the command body), but returns a plain negative int (no exception) for run-time aborts — empty driven set, bad args caught in the core, solve/fluid failures, maxIter (−1/−2/−3/−6/−7). `assert rc == 0` catches the second class only.
- **Why:** the run path deliberately mirrors classic `analyze` (negative return, "failed, returned: N error flag" print) so scripted retry logic works the same for both commands; the command-body fatals happen before a result exists, so the wrapper turns them into exceptions (`Py_ops_...` NULL-return convention).
- **Workaround/status (2026-07-17, ADR-73 P2):** by design; check the integer return like you would for `analyze`, and wrap in try/except only for setup misuse. Panel 2.D robustness-9 recorded the split; the P2 battery handles both surfaces.
