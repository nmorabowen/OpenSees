---
wp: PR-793
title: "793 -- upstreamable-table row(s)"
pr: "#793"
files: ["`SRC/classTags.h`", "`SRC/material/section/CMakeLists.txt`", "`SRC/material/section/Makefile`", "`SRC/Makefile`", "`SRC/interpreter/OpenSeesSectionCommands.cpp`", "`SRC/material/section/TclModelBuilderSectionCommand.cpp`", "`SRC/runtime/commands/modeling/section.cpp`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/runtime/runtime/TclPackageClassBroker.cpp`"]
table: "upstreamable"
legacy_seq: [504, 505, 506, 507, 508, 509, 510, 511, 512]
---
| `SRC/classTags.h` | `// Ladruno` (ADR 91): register `SEC_TAG_LadrunoShellModifier` = 33000 (SEC_TAG registry; per-registry band, not a collision with the 33000 used in NUMBERER_TAG/RECORDER_TAGS/INTEGRATOR_TAGS). | 793 |
| `SRC/material/section/CMakeLists.txt` | `# Ladruno` (ADR 91): add `LadrunoShellModifierSection.{cpp,h}` to `OPS_Material` sources. | 793 |
| `SRC/material/section/Makefile` | `# Ladruno` (ADR 91): add `LadrunoShellModifierSection.o` to `OBJS`. | 793 |
| `SRC/Makefile` | `# Ladruno` (ADR 91): add `$(FE)/material/section/LadrunoShellModifierSection.o`. | 793 |
| `SRC/interpreter/OpenSeesSectionCommands.cpp` | `// Ladruno` (ADR 91): extern `void* OPS_LadrunoShellModifierSection();` + `functionMap.insert(..."LadrunoShellModifier"...)` (Python/Tcl-in-interpreter `section` verb dispatch). | 793 |
| `SRC/material/section/TclModelBuilderSectionCommand.cpp` | `// Ladruno` (ADR 91): extern declaration + `else if (strcmp(argv[1],"LadrunoShellModifier")==0)` dispatch arm (classic-Tcl `section` command, Tcl path 1). | 793 |
| `SRC/runtime/commands/modeling/section.cpp` | `// Ladruno` (ADR 91): `extern OPS_Routine OPS_LadrunoShellModifierSection;` + dispatch arm (Xara-derived `section` command, Tcl path 2). | 793 |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno` (ADR 91): `#include "LadrunoShellModifierSection.h"` + `case SEC_TAG_LadrunoShellModifier: return new LadrunoShellModifierSection();` (OpenSeesSP/MP parallel broker). | 793 |
| `SRC/runtime/runtime/TclPackageClassBroker.cpp` | `// Ladruno` (ADR 91): `#include "LadrunoShellModifierSection.h"` + `case SEC_TAG_LadrunoShellModifier: return new LadrunoShellModifierSection();` (Xara-derived Tcl package broker, used by sendSelf/recvSelf and database restore). | 793 |
