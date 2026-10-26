---
wp: WP-213
title: "residualModes wiring (classic-Tcl bridge, DL declaration, CMake)"
date: 2026-10-26
table: "main"
files: ["`SRC/tcl/commands.cpp`", "`SRC/tcl/commands.h`", "`SRC/interpreter/OpenSeesCommands.h`", "`SRC/analysis/analysis/CMakeLists.txt`"]
---
| `SRC/tcl/commands.cpp` | `// Ladruno WP-213`: `extern int OPS_LadrunoResidualModes(void)` and the `residualModes` classic-Tcl bridge (`OPS_ResetInputNoBuilder` + call, the ADR-44 modal-bridge shape). Registered by the `LadrunoCommandTable.h` row, not by a `Tcl_CreateCommand` here. | — |
| `SRC/tcl/commands.h` | `// Ladruno WP-213`: declaration of the `residualModes` bridge next to the ADR-44 modal family. | — |
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno WP-213`: `int OPS_LadrunoResidualModes();` for the generated Python / DL-Tcl bridges. | — |
| `SRC/analysis/analysis/CMakeLists.txt` | Adds `LadrunoRitzKernel.{h,cpp}` and `LadrunoResidualModes.cpp` (WP-213) to `OPS_Analysis`. | — |
