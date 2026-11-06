---
wp: PR-555
title: "555 -- 1 vanilla row(s)"
pr: "#555"
files: ["`SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}`"]
table: "main"
legacy_seq: [299]
---
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature regen via `patch_banner.py` — note the `-load`/`-series` transient channel on the `modalResponseHistory` line (ADR44 transient -load). Banner strings only; the code lives in fork-authored `LadrunoModalResponse.{h,cpp}` (zero other vanilla edits). | [#555](https://github.com/nmorabowen/OpenSees/pull/555) |
