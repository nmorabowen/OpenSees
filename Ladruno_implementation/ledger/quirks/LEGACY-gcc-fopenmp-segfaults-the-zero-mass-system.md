---
wp: LEGACY
title: "gcc + -fopenmp SEGFAULTS the zero-mass system Diagonal path — the fork cannot be built with OpenMP on Linux"
legacy_seq: 469
---
## gcc + `-fopenmp` SEGFAULTS the zero-mass `system Diagonal` path — the fork cannot be built with OpenMP on Linux

- **Symptom.** Build the fork with `-DLADRUNO_OPENMP=ON` on gcc/Linux and
  `tests/test_adr30_projection_p0.py::test_massless_dof_is_not_policeable_by_the_soe_layer`
  dies with `Fatal Python error: Segmentation fault`, taking the whole pytest process with it
  (`Segmentation fault (core dumped)`, exit **139**) at the second test in the suite. Measured
  twice on the same commit, byte-identical traceback both times:
  [run 35164371356](https://github.com/nmorabowen/OpenSees/actions/runs/35164371356) (PR #843).
  The Python frame is `ops.analyze(1, 0.001)` inside `_run_massless(("Diagonal",), 0.0)`.
- **It is NOT the threaded loop.** The runtime thread count defaults to 1, at which
  `Domain::ladrunoThreadedUpdate()` returns `false` before touching anything and
  `Domain::update()` runs the unchanged serial loop. No `ladrunoThreads` call appears anywhere
  in the crashing file. Nor is it a dormant-pragma activation (`#pragma omp` / `_OPENMP` exist
  only in PFEM — `OPS_Element`, deliberately un-flagged — the interpreter, and WP-107's own
  code), nor an ODR/ABI split (no class layout is `#ifdef`-conditional).
- **What was left was "codegen/link" — and it was LINK (see RESOLVED below; this bullet is the pre-fix reasoning).** `-fopenmp` on `OPS_Domain` + `OPS_Utilities`, plus libgomp
  and `-pthread` on the link line, changes optimization and the glibc allocator's threading
  path — enough to turn a **latent defect in the singular-mass failure path** into a hard crash.
  That path is already on record two entries' worth: a free DOF with **zero lumped mass** makes
  the assembled `M` singular, `Diagonal` aborts with `aii = 0`, and Full/Band **return success
  with garbage**. The failure route also re-enters `Domain::update()` from
  `Domain::revertToLastCommit()` on the shared element iterator — the recorded reentrancy trap.
- **Does NOT reproduce on MSVC.** The same source with `LADRUNO_OPENMP=ON` passes that file
  locally (3/3), and the whole WP-107 file passes 18/18.
- **Status 2026-09-16 (PR #843):** banked, not fixed — `option(LADRUNO_OPENMP … OFF)`, build.bat ON, Zone-A
  skipping the WP-107 file. Superseded the next day:
- **RESOLVED 2026-09-17 ([WP-109](https://github.com/nmorabowen/OpenSees/pull/846)) — root cause is this fork's `CMakeLists.txt`, not gcc.** The runner
  itself ran Zone-A under gdb ([run 35171085324](https://github.com/nmorabowen/OpenSees/actions/runs/35171085324)): the fault is in
  `std::codecvt<char16_t>::do_unshift` called from `std::ostream::_M_insert<long>` from
  `PythonStream::err_out<int>` from `DiagonalDirectSolver::solve`'s `opserr << i` — a stream asked its
  locale for `num_put` and got a `codecvt` facet, i.e. **two libstdc++ runtimes in one process**. The GNU
  branch sets `CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++"`; `FindOpenMP` probes with an
  executable `try_compile`, sees `-Bstatic -lstdc++`, and records `OpenMP_CXX_LIB_NAMES =
  libstdc++;gomp;pthread` with `OpenMP_libstdc++_LIBRARY = …/libstdc++.a`; WP-107 appended
  `${OpenMP_CXX_LIBRARIES}` to the SHARED Python modules, which also `DT_NEED` `libstdc++.so.6`.
  esmeralda's module exported 179 libstdc++ internals; `LD_DEBUG=bindings` bound `num_put::id` to both
  copies. gcc 13 / 24.04 crashed, gcc 11 / 22.04 did not — load-order luck, not a compiler difference;
  esmeralda (Release, and ASAN) never reproduced, which is why the runner had to be the debugger.
  **Fix:** probe `FindOpenMP` with the EXE flags cleared AND link only `LADRUNO_OPENMP_LINK_LIBS` (the
  OpenMP runtime, everything else dropped loudly, `FATAL_ERROR` if a static C++/pthread runtime survives);
  `option(LADRUNO_OPENMP … ON)`. **Pinned:** `tests/test_wp109_module_single_libstdcxx.py` (pure-Python ELF
  `.dynsym` reader: no libstdc++ symbol defined, `libstdc++.so.6` NEEDED; red on the unfixed module, green
  after). Zone-A now runs the WP-107 file on Linux. See ADR-75b §14.5, BUILD_GOTCHAS §16,
  WORKFLOW_GOTCHAS §10.
- **Generalisation worth keeping:** any `find_package` that probes by linking an executable can return a
  static system runtime when `CMAKE_EXE_LINKER_FLAGS` carries `-static-*`; never append its `_LIBRARIES`
  to a `SHARED` target unread, and a Python extension must never define `std::locale`/`std::ios_base`
  symbols. `nm -D --defined-only module.so | grep -c _ZNSt6locale` must print 0.
- **Workaround/status (2026-09-16, PR #843, owner decision — not fixed).** `option(LADRUNO_OPENMP … OFF)`
  in `CMakeLists.txt`; `Ladruno_scripts\build.bat` turns it ON, so the ON path is the Windows/MSVC
  canonical build and nothing else. Consequence to keep in view: **Zone-A does not exercise the
  threaded element loop at all** — `tests/test_wp107_threaded_update.py` probes the binary and
  skips with a reason. Fixing this needs a Linux build under **ASAN/valgrind + gdb** as its own
  work package; when it lands, flip the CMake default to ON so CI gates the feature.
  See ADR-75b §14.4, [[107_ladruno_openmp_element_loop]] §3.1, BUILD_GOTCHAS §15.
