---
wp: PR-47
title: "47 -- 1 vanilla row(s)"
pr: "#47"
files: ["`SRC/matrix/ID.cpp`"]
table: "main"
legacy_seq: [116]
---
| `SRC/matrix/ID.cpp` | `// Ladruno P4`: same profiler memory hooks on the per-object `data` `new[]`/`delete[]` sites (ctor×2 + `malloc` branch / copy / dtor / setData / `unique` / `operator[]`-grow / `resize`-grow / `operator=` / `insert`-grow), tagged `ops_profiler::ALLOC_ID`; `arraySize` is the allocated capacity (exact frees). Bare-`if` deletes brace-converted; `operator=` captures the old `arraySize` in a temp before its overwrite. Additive. | [#47](https://github.com/nmorabowen/OpenSees/pull/47) |
