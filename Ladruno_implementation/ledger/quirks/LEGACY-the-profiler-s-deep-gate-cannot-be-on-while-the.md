---
wp: LEGACY
title: "The profiler's deep gate cannot be on while the loop is threaded"
legacy_seq: 455
---
### The profiler's deep gate cannot be on while the loop is threaded

`OPS_PROFILE_SCOPE_DEEP_NAMED` is gated on `enabled() && deep()`
(`ProfilerMacros.h:90-93`), and every `~ElemScope` does a lazy `std::map` insert
plus counter read-modify-writes on the node the master thread built
(`Profiler.cpp:115-124`). Concurrent `std::map` insertion is undefined behaviour,
and `Profiler.h:59-63` states the precondition ("each thread owns its own tree")
that this violates. So `Domain::ladrunoThreadedUpdate()` **refuses** to thread
while the deep gate is armed.

The practical consequence for anyone benchmarking: the per-loop *fraction* is a
property of the SERIAL baseline, so measure it with the deep profiler at 1
thread, and measure the threaded runs on wall time with the profiler off. Trying
to profile a threaded run deeply gets you a serial run and a confusing table.
