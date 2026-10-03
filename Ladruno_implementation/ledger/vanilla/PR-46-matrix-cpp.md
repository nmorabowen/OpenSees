---
wp: PR-46
title: "46 -- 2 vanilla row(s)"
pr: "#46"
files: ["`SRC/matrix/Matrix.cpp`", "`SRC/matrix/Vector.cpp`"]
table: "main"
legacy_seq: [114, 115]
---
| `SRC/matrix/Matrix.cpp` | `// Ladruno P4`: profiler memory hooks — `OPS_PROFILE_COUNT_ALLOC` after each per-object `data` `new[]`, `OPS_PROFILE_COUNT_FREE` before each `delete[]` (under the `fromFree==0` ownership guard), tagged `ops_profiler::ALLOC_MATRIX`. Runtime-gated on `mem()`; the static `matrixWork`/`intWork` scratch is deliberately not counted (never-freed whitelist). Additive. | [#46](https://github.com/nmorabowen/OpenSees/pull/46) |
| `SRC/matrix/Vector.cpp` | `// Ladruno P4`: same profiler memory hooks on the per-object `theData` `new[]`/`delete[]` sites (ctor/copy/dtor/setData/resize/`operator[]`/`operator=`/move-`operator=`), tagged `ops_profiler::ALLOC_VECTOR`. Runtime-gated on `mem()`. Additive. | [#46](https://github.com/nmorabowen/OpenSees/pull/46) |
