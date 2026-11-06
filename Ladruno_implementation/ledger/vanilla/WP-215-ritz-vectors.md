---
wp: WP-215
title: "ritzVectors wiring; analysis stamp accessors; Linear -factorOnce re-arms on domainChanged"
date: 2026-10-26
pr: "#1027"
table: "main"
files: ["`SRC/tcl/commands.cpp`", "`SRC/tcl/commands.h`", "`SRC/interpreter/OpenSeesCommands.h`", "`SRC/analysis/analysis/CMakeLists.txt`", "`SRC/analysis/analysis/StaticAnalysis.h`", "`SRC/analysis/analysis/DirectIntegrationAnalysis.h`", "`SRC/analysis/algorithm/equiSolnAlgo/Linear.h`"]
---
| `SRC/tcl/commands.cpp` | `// Ladruno WP-215`: `extern int OPS_LadrunoRitzVectors(void)` and the `ritzVectors` classic-Tcl bridge (`OPS_ResetInputNoBuilder` + call, the `residualModes` shape). Registered by the `LadrunoCommandTable.h` row. | #1027 |
| `SRC/tcl/commands.h` | `// Ladruno WP-215`: declaration of the `ritzVectors` bridge next to `residualModes`. | #1027 |
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno WP-215`: `int OPS_LadrunoRitzVectors();` for the generated Python / DL-Tcl bridges. | #1027 |
| `SRC/analysis/analysis/CMakeLists.txt` | Adds `LadrunoRitzVectors.cpp` (WP-215) to `OPS_Analysis`. | #1027 |
| `SRC/analysis/analysis/StaticAnalysis.h` | `// Ladruno WP-215`: inline `int ladrunoDomainStamp() const` returning the private `domainStamp`, so `ritzVectors` renumbers only when the domain stamp moved (eigen's rule). Additive. | #1027 |
| `SRC/analysis/analysis/DirectIntegrationAnalysis.h` | `// Ladruno WP-215`: the same `ladrunoDomainStamp()` accessor; mid-transient an unconditional `domainChanged()` reset `CentralDifference`'s U(t − Δt). Additive. | #1027 |
| `SRC/analysis/algorithm/equiSolnAlgo/Linear.h` | `// Ladruno WP-215`: `Linear::domainChanged()` re-arms `-factorOnce` (2 → 1). Before, a renumbered model or a command that factors another matrix into the analysis' LinearSOE (`ritzVectors`, `residualModes`) left the next step solving with a stale factor. Behaviour change only for `-factorOnce` after a `domainChanged()`; upstream candidate. | #1027 |
