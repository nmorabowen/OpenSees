---
wp: PR-64
title: "64 -- 3 vanilla row(s)"
pr: "#64"
files: ["`SRC/classTags.h`", "`SRC/interpreter/OpenSeesOutputCommands.cpp`", "`SRC/recorder/CMakeLists.txt`"]
table: "main"
legacy_seq: [73, 74, 75]
---
| `SRC/classTags.h` | Register `RECORDER_TAGS_LadrunoMonitorRecorder`=33002 (live analysis-monitor recorder) | [#64](https://github.com/nmorabowen/OpenSees/pull/64) |
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | Register `recorder Monitor` keyword → `OPS_LadrunoMonitorRecorder` (HDF5 block) | [#64](https://github.com/nmorabowen/OpenSees/pull/64) |
| `SRC/recorder/CMakeLists.txt` | Add `LadrunoMonitorRecorder` + `LadrunoMonitorSink` to the recorder target (HDF5≥1.12 block, SWMR) | [#64](https://github.com/nmorabowen/OpenSees/pull/64) |
