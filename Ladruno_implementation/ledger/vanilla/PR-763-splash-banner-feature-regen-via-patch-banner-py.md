---
wp: PR-763
title: "Splash-banner feature regen via patch_banner.py -- the LadrunoContact line's trailing 2D clause gains ; T4: radial end-cap open-terminal vertices (C0, replaces…"
pr: "#763"
files: ["`SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}`"]
table: "upstreamable"
legacy_seq: [478]
---
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature regen via `patch_banner.py` -- the `LadrunoContact` line's trailing 2D clause gains `; T4: radial end-cap open-terminal vertices (C0, replaces the NTS2D_END_SLACK window), Hertz cylinder-on-plane benchmark, node-union NTS/mortar routing match -- 2D lane complete)`. Source of truth `Ladruno_scripts/banner_features.txt`; never hand-edited. **No other vanilla file is touched by T4** -- the two C++ files it edits (`LadrunoContactFE.cpp`, `LadrunoContactHandler.cpp`) are both Ladruno-original (see [[LEDGER_implementations]]'s T4 row). PR: [#763](https://github.com/nmorabowen/OpenSees/pull/763) -- ADR-85 T4. |
