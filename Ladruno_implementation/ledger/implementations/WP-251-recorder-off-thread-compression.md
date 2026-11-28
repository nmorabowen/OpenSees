---
wp: WP-251
title: "WP-251 — recorder ladruno: off-thread chunk compression and threaded reactions"
pr: "#1146"
date: 2026-11-23
status: "draft"
section: "table"
banner: "Ladruno recorder — off-thread chunk compression, threaded reactions (WP-251)"
---
| **WP-251 — recorder ladruno: off-thread compression + threaded reactions.** A compressed streaming `DATA` (one step per chunk, WP-233) is shuffled + deflated on worker threads and written by `H5Dwrite_chunk` on the solver thread in step order with its `TIME`/`STEP` rows — same bytes on disk; only the solver thread calls HDF5. `-compressThreads <n>` (`-1` auto, `0` = in-line filter). Rows in flight are written before every flush, layout step, close and process exit (an at-exit hook that compresses on the exiting thread). `Domain::ladrunoThreadedReactions`: under `ladrunoThreads > 1`, `calculateNodalReactions` (flag 0/1) computes `Element::ladrunoThreadSafeReaction()` forces in parallel and scatters in serial order (bit-identical); quirk rule `reaction-scatter`. Recorder 124 → 44 ms/step serial, 99 → 19 ms/step at 8 threads (8192 stdBricks, all six results). | recorder performance | **none — no new classes** | `SRC/recorder/{Ladruno_ChunkCompress.{h,cpp}, Ladruno_Sinks.{h,cpp}, LadrunoRecorder.cpp, Ladruno_ResultIO.h, Ladruno_Types.h, CMakeLists.txt}`, `SRC/domain/domain/Domain.{h,cpp}`, `SRC/element/Element.h`, `SRC/element/absorbentBoundaries/ASDAbsorbingBoundary{2D,3D}.h`, `ci/{check_quirk_patterns.py, test_quirk_reaction_scatter.py, README.md}`, `tests/test_ladruno_recorder_compress_threads.py`, `Ladruno_implementation/{251_ladruno_recorder_off_thread_compression.md, 03_ladruno_recorder.md}` | draft | [#1146](https://github.com/nmorabowen/OpenSees/pull/1146) |
