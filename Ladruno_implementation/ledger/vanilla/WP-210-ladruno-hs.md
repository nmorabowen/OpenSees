---
wp: WP-210
title: "LadrunoHS / LadrunoHSSmall: register ND_TAG_LadrunoHS (33029), OPS_LadrunoHS / OPS_LadrunoHSSmall, CMake, object broker, and splash banner"
date: 2026-10-24
pr: "—"
table: "main"
files: ["`SRC/classTags.h`", "`SRC/material/nD/CMakeLists.txt`", "`SRC/interpreter/OpenSeesNDMaterialCommands.cpp`", "`SRC/material/nD/TclModelBuilderNDMaterialCommand.cpp`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/tcl/tclMain.cpp`", "`SRC/interpreter/PythonModule.cpp`"]
---
| `SRC/classTags.h` | `// Ladruno (WP-210)`: allocate `ND_TAG_LadrunoHS 33029` in the Ladruno private nD material band. | — |
| `SRC/material/nD/CMakeLists.txt` | `// Ladruno (WP-210)`: add `LadrunoHS.cpp`, `LadrunoHS.h`, and `LadrunoHSKernel.h` to `OPS_Material`. | — |
| `SRC/interpreter/OpenSeesNDMaterialCommands.cpp` | `// Ladruno (WP-210)`: declare and register `OPS_LadrunoHS` (`"LadrunoHS"`) and `OPS_LadrunoHSSmall` (`"LadrunoHSSmall"`) in `nDMaterialsMap`. | — |
| `SRC/material/nD/TclModelBuilderNDMaterialCommand.cpp` | `// Ladruno (WP-210)`: declare and dispatch `"LadrunoHS"` and `"LadrunoHSSmall"` in the Tcl `nDMaterial` `strcmp` ladder for Tcl/Python parity. | — |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno (WP-210)`: include `LadrunoHS.h` and instantiate `new LadrunoHS()` on `ND_TAG_LadrunoHS` in `getNewNDMaterial`. | — |
| `SRC/tcl/tclMain.cpp` | Regenerated `FEATURES-START/END` block via `patch_banner.py` for WP-210 (`LadrunoHS/HSSmall`). | — |
| `SRC/interpreter/PythonModule.cpp` | Regenerated `FEATURES-START/END` block via `patch_banner.py` for WP-210 (`LadrunoHS/HSSmall`). | — |
