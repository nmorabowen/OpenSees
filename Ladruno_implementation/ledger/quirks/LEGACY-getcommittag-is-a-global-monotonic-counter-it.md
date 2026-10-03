---
wp: LEGACY
title: "getCommitTag() is a GLOBAL monotonic counter — it does NOT reset on wipe()"
legacy_seq: 42
---
### `getCommitTag()` is a GLOBAL monotonic counter — it does NOT reset on `wipe()`
- **Bites:** any per-step recorder/series that uses `Domain::getCommitTag()` as its
  step axis (the analysis monitor's `STEP`, the profiler per-step series). Across
  several `analyze()` runs in one interpreter session — even with `wipe()` between
  them — the commitTag keeps climbing (run 1 → 0..199, run 2 → 200..399, ...). It is
  NOT a within-analysis 0-based step index.
- **Implication:** don't assert absolute step values or compare step arrays across
  runs by value; compare the *stride* (`np.diff(step) == every`) instead. For a live
  viewer, treat `STEP` as a monotonic id, not "step N of this analysis."
- Learned 2026-05-31 building the analysis monitor (`08_analysis_monitor.md`).
