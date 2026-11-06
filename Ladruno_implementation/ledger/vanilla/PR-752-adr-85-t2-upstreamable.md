---
wp: PR-752
title: "752 -- ADR-85 T2 -- upstreamable-table row(s)"
pr: "#752"
files: ["`SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}`"]
table: "upstreamable"
legacy_seq: [475]
---
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature regen via `patch_banner.py` -- the `LadrunoContact` line's 2D clause now reads `T1b: frictionless NTS pairs + concave-vertex corners live (...); T2: friction live (scalar 1D return map, unified Coulomb/Tresca cone, -consistanttan, softKt, implicit-transient + -initial arms, removal incl. the vertex case); 2D mortar/tie lands in T3 -- refused by name` (honest scoping: mortar/tie is still T3). Banner strings only, regenerated from `Ladruno_scripts/banner_features.txt` (repo rule -- never hand-edited). **No other vanilla file is touched by T2** -- the four C++ files it edits (`LadrunoContactFE.{h,cpp}`, `LadrunoContactHandler.cpp`, `LadrunoFrictionKernel.h`) are all Ladruno-original. | [#752](https://github.com/nmorabowen/OpenSees/pull/752) -- ADR-85 T2 |
