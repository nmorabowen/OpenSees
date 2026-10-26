---
wp: WP-213
title: "residualModes wiring + RSA Sa-period override for residual modes"
date: 2026-10-26
pr: "#1023"
table: "main"
files: ["`SRC/tcl/commands.cpp`", "`SRC/tcl/commands.h`", "`SRC/interpreter/OpenSeesCommands.h`", "`SRC/analysis/analysis/CMakeLists.txt`", "`SRC/analysis/analysis/ResponseSpectrumAnalysis.cpp`"]
---
| `SRC/tcl/commands.cpp` | `// Ladruno WP-213`: `extern int OPS_LadrunoResidualModes(void)` and the `residualModes` classic-Tcl bridge (`OPS_ResetInputNoBuilder` + call, the ADR-44 modal-bridge shape). Registered by the `LadrunoCommandTable.h` row, not by a `Tcl_CreateCommand` here. | #1023 |
| `SRC/tcl/commands.h` | `// Ladruno WP-213`: declaration of the `residualModes` bridge next to the ADR-44 modal family. | #1023 |
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno WP-213`: `int OPS_LadrunoResidualModes();` for the generated Python / DL-Tcl bridges. | #1023 |
| `SRC/analysis/analysis/CMakeLists.txt` | Adds `LadrunoRitzKernel.{h,cpp}` and `LadrunoResidualModes.cpp` (WP-213) to `OPS_Analysis`. | #1023 |
| `SRC/analysis/analysis/ResponseSpectrumAnalysis.cpp` | `// Ladruno WP-213`: per-mode Sa-period override read from `LadrunoModeKinds` (set by `residualModes -period $T`). `solveMode()` (per mode, full loop, `-increment`) and `analyzeCombined()` (`-combine`) look Sa up at the override, and `analyzeCombined()` also uses ω = 2π/T in the CQC ρ and T in TenPercent (one rigid mode, ASCE 7-22 §12.9.1.1); λ in φΓSa/λ stays the mode's own. No override recorded => byte-identical to before. | #1023 |
