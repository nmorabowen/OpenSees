---
wp: PR-6
title: "6, #65, #108 -- 2 vanilla row(s)"
pr: "#6, #65, #108"
files: ["`SRC/interpreter/OpenSeesElementCommands.cpp`", "`SRC/element/CMakeLists.txt`"]
table: "main"
legacy_seq: [20, 21]
---
| `SRC/interpreter/OpenSeesElementCommands.cpp` | Register `element` dispatch for `BezierTri6`, `BezierTet10`, and `LadrunoBrick`/`ladrunoBrick` (fwd-decl + `functionMap`) | [#6](https://github.com/nmorabowen/OpenSees/pull/6), [#65](https://github.com/nmorabowen/OpenSees/pull/65), [#108](https://github.com/nmorabowen/OpenSees/pull/108) |
| `SRC/element/CMakeLists.txt` | `add_subdirectory(bezierTriangle)` + `add_subdirectory(bezierTetrahedron)` + `add_subdirectory(solidTransformation)` + `add_subdirectory(ladrunoBrick)` | [#6](https://github.com/nmorabowen/OpenSees/pull/6), [#65](https://github.com/nmorabowen/OpenSees/pull/65), [#108](https://github.com/nmorabowen/OpenSees/pull/108) |
