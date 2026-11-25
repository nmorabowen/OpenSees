---
wp: WP-253
title: "ElasticIsotropicPlaneStrain2D is Update-safe but not Form-safe: getStress writes a class-wide static"
pr: "#1147"
date: 2026-11-23
---
### A material on the Update allowlist can still race in a host that calls getStress / getInitialTangent

`ElasticIsotropicPlaneStrain2D::ladrunoThreadSafeUpdate()` is `true`, and `ladrunoThreadSafeForm()` is
`false`, for a reason: `getStress()` and `getInitialTangent()` return the class-wide statics `sigma` and
`D`. The Update contract (WP-107) covers `setTrialStrain` alone. A host that reads `getStress()` right
after `setTrialStrain()` on worker threads, as the MPM particle update does, races on that static
under the Update answer alone. Mode F's threaded material loop (WP-253) therefore requires
`ladrunoThreadSafeUpdate() && ladrunoThreadSafeForm()` on every particle material; the 2D elastic
particles run their material calls serially in particle order while the transfer loops thread.
Check the Form answer too whenever a threaded caller does more than `setTrialStrain`.

#### Also learned (WP-253): a threaded commit cannot name its refuser from the shared counter

The WP-99 refusal counter is process-wide (atomic since WP-224). The serial MPM commit names the
first refusing particle from a before/after read of that counter (WP-232k D13); under threads,
another worker can bump it between the two reads. WP-253 adds a `thread_local` tally next to it
(`ladrunoThreadCommitRefusals()`, bumped in the same place, never reset), whose before/after
difference on the committing thread is exact. The global counter's meaning is unchanged.

#### Also learned (WP-253): running `ci_run_pytest.py` on nmb

`python3.12 -S Ladruno_scripts/ci_run_pytest.py` imports pytest from the interpreter's `purelib`
only; on nmb pytest lives in the **user** site (`%APPDATA%\Python\Python312\site-packages`), so set
`PYTHONPATH` to it or the runner dies on `import pytest`. The runner `chdir`s into `tests/`: pass
test paths relative to `tests/` (`test_x.py::case`). From Git Bash a `file.py::id` argument with a
leading path gets mangled into a `/c/...` path that pytest cannot find; run it from PowerShell.
`LADRUNO_TCL_EXE` points the Tcl-parity tests at an `OpenSees.exe` when the worktree has no `dist`.

#### Also learned (WP-253): an exception inside `ladruno_parallelFor` is undefined behaviour

A `std::bad_alloc` (or any C++ exception) escaping a body run under MSVC's OpenMP terminates the
process; there is no propagation to the caller. The threaded MPM step sizes most buffers before
its parallel regions, but the transpose's per-block discovery lists still `push_back` inside a body
(`LadrunoMPMThreads.h`, `transposeThreaded`), so an out-of-memory there ends the run without the
serial path's error. WP-192's loops carry the same exposure. Accepted for WP-253 (desktop memory).

#### Also learned (WP-253): libiomp's KMP_BLOCKTIME spin starves a serial phase between regions

The build links Intel `libiomp5md` (default `KMP_BLOCKTIME` 200 ms). After each parallel region the
workers spin for that long, so a serial loop that follows (MPM's serial material calls when the
materials are not thread-safe) shares its core with a spinning hyperthread: measured 2D SANISAND
S(8) = 0.89 and 2D DP 1.18 / 1.23, and 1.13 / 1.30 / 1.33 with `KMP_BLOCKTIME=0`. The fix sits in
the code, not the environment: `ladruno_ompQuietBegin/End` (`LadrunoThreads.{h,cpp}`) resolve
`kmp_set_blocktime` / `kmp_get_blocktime` at run time and set the master's blocktime to 0 around
the serial phase. **Timing trap:** libiomp copies the master's ICVs to the team at the FORK, and each
worker takes its post-region wait from them at the JOIN barrier, so the drop must come before the
fork of the region that PRECEDES the serial phase. Set after that join, it governs only the next
region's tail. `get` returns the value in `KMP_BLOCKTIME`'s own units, so `set(get())` restores
the VALUE exactly (checked: 200 ms, 500 us). It does not restore the "user-set" flag: libiomp's
`kmp_set_blocktime` sets `bt_set = TRUE` in the master's ICVs, and from the first span on that
switches off the runtime's oversubscription auto-zero (`KMP_ADJUST_BLOCKTIME`: blocktime 0 while
more threads than processors exist and `bt_set` is false). If a later phase oversubscribes (more
OpenMP and MKL threads than cores), the workers spin out the full blocktime there instead of
sleeping at once. This affects timing only, never a result. With `KMP_BLOCKTIME=infinite`, `get`
reports INT_MAX whatever is set, and `set(INT_MAX)` prints an `OMP: Info #109` line per call, so
the span changes nothing. The entry points come from the module that owns `omp_get_max_threads`
(an MSVC import thunk is followed to the IAT slot; checked against oneAPI 2026.1's
`libiomp5md.dll`), never from the first DLL found by name. `mpm debug quietSpan` says whether they
resolved.

#### Also learned (WP-253): `profiler stop` leaves `-deep` armed; Domain's threaded element loop then stays serial

`profiler start -deep` sets the profiler's deep flag. Neither `profiler stop` nor `profiler reset`
clears it, and Domain's WP-192 gate tests `deep()` alone. After one deep-profiled run, every later
threaded element loop in the same process runs SERIAL ("the DEEP profiler gate is armed"). A test
that deep-profiles must end with `profiler start` (no `-deep`), then `stop` and `reset`, or the
WP-192 tests that run after it in the same pytest process (Zone-A runs them all in one) fail.
WP-253's deep case did exactly this and broke four WP-192 tests. The gate fix is #1150.

#### Also learned (WP-253): an atomic claim that stores an owner id must use compare-exchange, not exchange

The threaded transpose claims each grid node with an atomic mark. `load() == 0 && exchange(1) == 0`
is correct while the mark is only a flag: a losing `exchange` writes 1 over 1. When the transpose
began storing the CLAIMING BLOCK in the mark (`b + 1`, for the windowed histogram), the same idiom
broke. The losing block's `exchange` overwrote the winner's id. Its own return value was still
right, but any block that loaded the mark later read the wrong owner. Its histogram window then
missed the node, and at 8 threads the result was heap corruption (`0xc0000374`) and non-identical
particles in about one run in five. 2 and 4 threads passed. Store an id with
`compare_exchange_strong(expected = 0, b + 1)`: a loser never writes, and gets the winner's id
back in `expected`. The windowed-histogram self-check that found it compared each entry's mark
with its slot's owner (`mark == owner + 1`).
