---
wp: WP-216
title: "responseSpectrumPerMode wiring; RSA check()/getSa()/checkMass() public; TransformationFE dofData sized by node count"
date: 2026-10-27
pr: "#1032"
table: "main"
files: ["`SRC/analysis/analysis/ResponseSpectrumAnalysis.h`", "`SRC/analysis/fe_ele/transformation/TransformationFE.cpp`", "`SRC/tcl/commands.cpp`", "`SRC/tcl/commands.h`", "`SRC/interpreter/OpenSeesCommands.h`", "`SRC/analysis/analysis/CMakeLists.txt`"]
---
| `SRC/analysis/analysis/ResponseSpectrumAnalysis.h` | `// Ladruno WP-216`: `check()`, `getSa()` and the WP-212 `checkMass()` moved from private to public so `responseSpectrumPerMode` applies the same modalProperties consistency check, the same per-direction mass guard and the same spectral ordinate (`-scale` included). No behaviour change. | #1032 |
| `SRC/analysis/fe_ele/transformation/TransformationFE.cpp` | `// Ladruno WP-216`: the static `dofData` holds one DOF count per NODE (`numDOFs.setData(dofData, numGroups)`) but was allocated with `MAX_NUM_DOF` = 64 entries; the constructor now grows it to the element's node count (file-static `ladrunoSizeDofData`), and the WP-192 thread-scratch copy is sized the same way. Before, an element with > 64 nodes (a LadrunoMembrane floor) wrote past it on every assembly under `constraints Transformation`: heap corruption (0xc0000374 at a later `wipe`). Upstream candidate (upstream has the same buffer). | #1032 |
| `SRC/tcl/commands.cpp` | `// Ladruno WP-216`: `extern int OPS_LadrunoResponseSpectrumPerMode(void)` and the `responseSpectrumPerMode` classic-Tcl bridge (the `ritzVectors` shape). Registered by the `LadrunoCommandTable.h` row. | #1032 |
| `SRC/tcl/commands.h` | `// Ladruno WP-216`: declaration of the `responseSpectrumPerMode` bridge. | #1032 |
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno WP-216`: `int OPS_LadrunoResponseSpectrumPerMode();` for the generated Python / DL-Tcl bridges. | #1032 |
| `SRC/analysis/analysis/CMakeLists.txt` | Adds `LadrunoRSAPerMode.cpp` (WP-216) to `OPS_Analysis`. | #1032 |
