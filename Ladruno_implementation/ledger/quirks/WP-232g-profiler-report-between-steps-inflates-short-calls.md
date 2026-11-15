---
wp: WP-232g
title: "A profiler report between steps inflates the per-call cost of short calls"
date: 2026-11-15
---
### A profiler report between steps inflates the per-call cost of short calls

- **Bites:** per-step profiling done as `profiler start` / one step / `stop` / `report` /
  `reset`, repeated. On a zero-free-DOF stdBrick driving ElasticIsotropic, the per-point cost
  (`elem.update` / 8 + `commit.elements` / 8, deep) read 85 ns with the 200 steps under one
  profiler run and 162 ns with each step its own run (WP-232g, build `84ff95616`). The report's
  HDF5 write and the tree reset evict what the next step needs, and a sub-microsecond call pays
  the cold start every time. A mode-A MPM step (≥ 0.1 s) does not notice it.
- **Also:** the profiler keeps one tree per OS thread, each rooted at `root`, and the report
  merges them by name. A scope opened on a worker thread is NOT nested under the scope that was
  open on the master: it appears at the top level, its wall summed over the threads. Bracket a
  threaded region on the master thread, never inside the worker body.
- **Workaround:** profile a contiguous block of steps as one run and take the per-step rows from
  `profiler start -perStep` (the `series/wall_ms` dataset); keep per-step runs for long steps
  only. `Ladruno_scripts/mpm_m0b/` does both (M2 blocks, M1 per step).
