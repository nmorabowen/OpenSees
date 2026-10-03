---
wp: PR-1
title: "1 -- 5 vanilla row(s)"
pr: "#1"
files: ["`SRC/interpreter/PythonMPIModule.cpp`", "`SRC/interpreter/OpenSeesMiscCommands.cpp`", "`SRC/interpreter/CMakeLists.txt`", "`SRC/actor/channel/CMakeLists.txt`", "`SRC/actor/address/CMakeLists.txt`"]
table: "main"
legacy_seq: [16, 17, 71, 86, 87]
---
| `SRC/interpreter/PythonMPIModule.cpp` | `openseesmp` MPI module entry point (Patch 9) — re-includes `PythonModule.cpp` with `OPS_PY_MODULE_NAME` redefined | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/interpreter/OpenSeesMiscCommands.cpp` | `partition` command METIS 5 API guard / error message (Patch 9) | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/interpreter/CMakeLists.txt` | Build wiring for `openseesmp` and new interpreter sources | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/actor/channel/CMakeLists.txt` | MP build wiring | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
| `SRC/actor/address/CMakeLists.txt` | MP build wiring | [#1](https://github.com/nmorabowen/OpenSees/pull/1) |
