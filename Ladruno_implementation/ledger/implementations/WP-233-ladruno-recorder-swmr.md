---
wp: WP-233
title: "WP-233 — LadrunoRecorder survives a hard kill (HDF5 SWMR)"
pr: "#1084"
date: 2026-11-14
status: "draft"
section: "table"
banner: "Ladruno recorder — SWMR writes, survives a hard kill (WP-233)"
---
| **WP-233 — LadrunoRecorder survives a hard kill** (#1082; schema §7b.2). The file uses the latest HDF5 format and is in SWMR-write mode except during a *layout step*: a stage's first recorded step, the `KIND` rewrite, a modes snapshot, or a late envelope. A layout step closes the file, reopens it read-write, creates the objects, arms `H5Fstart_swmr_write` and only then appends the step's rows, so no streaming row is written outside SWMR. Groups that gain links during a run (the root; a modes result's `DATA`) are presized and compact, and the provisional `KIND` is rewritten in place, so a kill in a layout step can damage only the stage being laid out. Streaming `DATA` layout that SWMR does not bloat: `-compress N > 0` (default 1) one step per chunk, written once; `-compress 0` unfiltered ~64 KiB chunks. Reader rule (plain, SWMR-read, repair a copy) and crash rule (cut to the shortest `DATA`/`TIME`/`STEP`) in `ladruno_format.py`, shared with the MPM particle reader. | recorder robustness | **none — no new classes** | `SRC/recorder/{LadrunoRecorder.{h,cpp}, Ladruno_Sinks.{h,cpp}, Ladruno_Hdf5.h, Ladruno_Types.h, Ladruno_ResultIO.h}`, `Ladruno_scripts/ladruno_recorder_tests/{ladruno_format.py, recorder_crash_replay.py, recorder_kill_stress.py, recorder_binary_kill_stress.py, recorder_swmr_layout_bench.py}`, `tests/mpm_contracts.py`, `Ladruno_implementation/{ladruno_schema_v1.md §7b.2, 164_ladruno_recorder_write_path.md, 03_ladruno_recorder.md}` | draft | [#1084](https://github.com/nmorabowen/OpenSees/pull/1084) |
