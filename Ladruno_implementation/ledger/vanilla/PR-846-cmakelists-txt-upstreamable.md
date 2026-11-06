---
wp: PR-846
title: "846 -- upstreamable-table row(s)"
pr: "#846"
files: ["`CMakeLists.txt` (root)"]
table: "upstreamable"
legacy_seq: [638]
---
| `CMakeLists.txt` (root) | `# Ladruno WP-109`: the gcc `-fopenmp` "segfault" of #843 was this file linking a SECOND, static libstdc++ into the SHARED Python modules — `FindOpenMP`'s executable probe inherited the GNU-branch `CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++"` and reported `libstdc++.a` (+ `libpthread.a`) as OpenMP implicit libraries, and WP-107 appended `${OpenMP_CXX_LIBRARIES}` to all five targets. Now: `find_package(OpenMP)` runs with the EXE linker flags temporarily cleared; the five link sites take `${LADRUNO_OPENMP_LINK_LIBS}` = `OpenMP_CXX_LIBRARIES` filtered to the OpenMP runtime (`gomp`/`omp`/`iomp5`), dropped entries announced, `FATAL_ERROR` if a static C++/pthread runtime survives (a stale cache still carries the bad `OpenMP_CXX_LIB_NAMES`); **`option(LADRUNO_OPENMP … ON)`** — the flip #843 deferred, so Zone-A gates WP-107. The 60-line option comment rewritten from "OFF because gcc" to the root cause. No `SRC/` change. Pinned by `tests/test_wp109_module_single_libstdcxx.py`. ADR-75b §14.5, BUILD_GOTCHAS §16. | [#846](https://github.com/nmorabowen/OpenSees/pull/846) |
