---
wp: WP-199
title: "LadrunoCPML element and CPML mutation-family registrations"
date: 2026-10-16
files: ["`CMakeLists.txt`", "`SRC/classTags.h`", "`DEVELOPER/core/classTags.h`", "`SRC/element/CMakeLists.txt`", "`SRC/element/TclElementCommands.cpp`", "`SRC/interpreter/OpenSeesElementCommands.cpp`", "`SRC/interpreter/OpenSeesOutputCommands.cpp`", "`SRC/interpreter/PythonModule.cpp`", "`SRC/tcl/commands.cpp`", "`SRC/tcl/tclMain.cpp`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`"]
---
| `CMakeLists.txt` | `// Ladruno (WP-199)`: added `CPML` to `_ladruno_valid_families` and the `LADRUNO_MUTATE_FAMILY` cache `STRINGS` property for the ADR-87 D2 mutation gate. | — |
| `SRC/classTags.h` | `// Ladruno (WP-199)`: defined `ELE_TAG_LadrunoCPML 33023` in the Ladruno element tag band. | — |
| `DEVELOPER/core/classTags.h` | `// Ladruno (WP-199)`: mirrored `ELE_TAG_LadrunoCPML 33023` to keep `DEVELOPER/core/classTags.h` in sync with `SRC/classTags.h`. | — |
| `SRC/element/CMakeLists.txt` | `// Ladruno (WP-199)`: added `add_subdirectory(ladrunoCPML)` to compile `LadrunoCPML` into `OPS_Element`. | — |
| `SRC/element/TclElementCommands.cpp` | `// Ladruno (WP-199)`: declared `OPS_LadrunoCPML` and registered `element LadrunoCPML` in classic Tcl (`TclElementCommands`). | — |
| `SRC/interpreter/OpenSeesElementCommands.cpp` | `// Ladruno (WP-199)`: declared `OPS_LadrunoCPML` and registered `LadrunoCPML` in the interpreter `functionMap` (`OpenSeesPy` and interpreter Tcl). | — |
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno (WP-199)`: added `{"CPML", LADRUNO_MUTATE_CPML}` to the `OPS_ladrunoMutation` reporting table. | — |
| `SRC/interpreter/PythonModule.cpp` | Regenerated `FEATURES-START/END` splash banner via `Ladruno_scripts/patch_banner.py` to include `LadrunoCPML`. | — |
| `SRC/tcl/commands.cpp` | `// Ladruno (WP-199)`: added `{"CPML", LADRUNO_MUTATE_CPML}` to the Tcl `ladrunoMutation` reporting table. | — |
| `SRC/tcl/tclMain.cpp` | Regenerated `FEATURES-START/END` splash banner via `Ladruno_scripts/patch_banner.py` to include `LadrunoCPML`. | — |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno (WP-199)`: `#include <LadrunoCPML.h>` and `case ELE_TAG_LadrunoCPML: return new LadrunoCPML();` in `getNewElement` for `sendSelf`/`recvSelf` checkpointing and parallel `OpenSeesSP`/`OpenSeesMP` serialization. | — |
