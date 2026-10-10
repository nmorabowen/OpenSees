---
wp: PR-2
title: "2, #8 -- 2 vanilla row(s)"
pr: "#2, #8"
files: ["`SRC/interpreter/OpenSeesOutputCommands.cpp`", "`SRC/recorder/CMakeLists.txt`"]
table: "main"
legacy_seq: [18, 72]
---
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | Register `recorder` keywords for EnergyBalance and Ladruno (`recorder ladruno`) — the shared map used by OpenSeesPy/openseesmp **and** the interpreter-based Tcl (`TclWrapper`) | [#2](https://github.com/nmorabowen/OpenSees/pull/2), [#8](https://github.com/nmorabowen/OpenSees/pull/8) |
| `SRC/recorder/CMakeLists.txt` | Add `EnergyBalanceRecorder` and `LadrunoRecorder` + `Ladruno_*` to the recorder target | [#2](https://github.com/nmorabowen/OpenSees/pull/2), [#8](https://github.com/nmorabowen/OpenSees/pull/8) |
