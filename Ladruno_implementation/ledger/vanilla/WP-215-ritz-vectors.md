---
wp: WP-215
title: "ritzVectors wiring"
date: 2026-10-26
pr: "#TBD"
table: "main"
files: ["`SRC/tcl/commands.cpp`", "`SRC/tcl/commands.h`", "`SRC/interpreter/OpenSeesCommands.h`", "`SRC/analysis/analysis/CMakeLists.txt`"]
---
| `SRC/tcl/commands.cpp` | `// Ladruno WP-215`: `extern int OPS_LadrunoRitzVectors(void)` and the `ritzVectors` classic-Tcl bridge (`OPS_ResetInputNoBuilder` + call, the `residualModes` shape). Registered by the `LadrunoCommandTable.h` row. | #TBD |
| `SRC/tcl/commands.h` | `// Ladruno WP-215`: declaration of the `ritzVectors` bridge next to `residualModes`. | #TBD |
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno WP-215`: `int OPS_LadrunoRitzVectors();` for the generated Python / DL-Tcl bridges. | #TBD |
| `SRC/analysis/analysis/CMakeLists.txt` | Adds `LadrunoRitzVectors.cpp` (WP-215) to `OPS_Analysis`. | #TBD |
