---
wp: WP-251
title: "recorder ladruno per-step cost is deflate + reactions, not the collect walk; worker threads are dead at .pyd atexit"
date: 2026-11-23
---
### recorder ladruno slows every step: it is deflate and reactionForce, not getResponse

Profiled with `lr.collect` / `lr.write` / `lr.reactions` (stack profiler, WP-251) on
8192 stdBricks: collecting results is ~7 % of the recorder; HDF5's shuffle + deflate,
run on the solver thread inside `H5Dwrite`, is ~60 %; `-N reactionForce` re-runs
`Domain::calculateNodalReactions` — a full resisting-force sweep, about one residual
assembly — every recorded step (~28 %). Measure before moving data into flat arrays
or threading `getResponse` (many elements fill class-static buffers there). WP-251 moved
the deflate to worker threads (same bytes, `H5Dwrite_chunk`) and threads the reaction
sweep under `ladrunoThreads`.

#### Worker threads at process exit (WP-251)

In `opensees.pyd` the CRT runs atexit handlers at `DLL_PROCESS_DETACH`, after Windows
has terminated every other thread. A handler that waits on a worker pool (or locks a
mutex a worker may hold) hangs the exit forever, and a static `std::thread` destroyed
joinable calls `std::terminate`. The recorder's at-exit drain therefore redoes queued
work on the exiting thread and never touches the pool; the pool is leaked and its
threads detached. The hook must be registered AFTER HDF5 initialises (atexit is LIFO)
so it runs before HDF5 closes every ID.
