---
wp: WP-258
title: "opserr from an OpenMP worker under openseespy segfaults: PythonStream called CPython without a thread state (the WP-107 ManzariDafalias crash)"
date: 2026-11-25
---
### `opserr` from an OpenMP worker under openseespy segfaults — PythonStream called CPython without a thread state (WP-258)

- **Symptom.** WP-107: ManzariDafalias IntScheme 1 segfaulted (`0xC0000005`) on a worker as
  soon as the plastic branch ran in volume. The elastic branch was clean. Serialising the whole
  `update()` under `omp critical` changed nothing, nor did larger stacks or another solver.
- **Cause.** The plastic branch prints: the ModifiedEuler substep-cap and low-p clamp warnings.
  Under openseespy `opserr` is a `PythonStream`, and `err_out` called `PySys_FormatStderr`. On
  a thread with no Python thread state that dereferences null. A worker has no thread state
  whether or not it runs under a critical section, which is why serialising did not help. The
  WP-107 logs show it: the serial twin prints the warning, and the threaded run dies right
  after the `THREADED` line.
- **What does NOT fix it.** `LadrunoOutputGuard` (a mutex). Taking the GIL from the worker
  (`PyGILState_Ensure`) deadlocks: the main thread is thread 0 of the parallel region and holds
  the GIL until the region ends.
- **The fix (`PythonStream.h`, WP-258).** `err_out` writes a local string to C `stderr` when
  `!PyGILState_Check()`. The GIL-holding thread keeps the Python path. Tcl (`StandardStream`)
  was never affected.
- **Rule.** A threaded material or element may print, but only through `opserr` under
  `LadrunoOutputGuard` (so lines do not interleave). Never call any `PyObject` / `Py*` API from
  a body of `ladruno_parallelFor`. A test that threads a printing path must run under
  openseespy, not only Tcl: `test_wp258_soil_threads.py::test_worker_prints_under_openseespy`.
- **Corollary for the "a static grep is an audit" quirk (WP-107).** That entry's conclusion,
  that the next tool is TSan, was right in spirit. The defect, though, was not a race: it was a
  thread-affinity bug in a shared stream. A race detector would not have flagged it.
