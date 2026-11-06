---
wp: PR-295
title: "295 -- 6 vanilla row(s)"
pr: "#295"
files: ["`SRC/classTags.h`", "`SRC/interpreter/OpenSeesCommands.{h,cpp}`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/runtime/runtime/TclPackageClassBroker.cpp`", "`SRC/analysis/integrator/CMakeLists.txt`", "`SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}`"]
table: "main"
legacy_seq: [180, 181, 182, 183, 184, 185]
---
| `SRC/classTags.h` | Register `INTEGRATOR_TAGS_CentralDifferenceSMS`=33007 (selective mass-scaling explicit integrator; per-registry band, distinct from `ELE_TAG_LadrunoQuad`=33007) | [#295](https://github.com/nmorabowen/OpenSees/pull/295) |
| `SRC/interpreter/OpenSeesCommands.{h,cpp}` | `// Ladruno`: declare `OPS_CentralDifferenceSMS()` + `strcmp(type,"CentralDifferenceSMS")` dispatch in the integrator command | [#295](https://github.com/nmorabowen/OpenSees/pull/295) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno`: `#include "CentralDifferenceSMS.h"` + `case INTEGRATOR_TAGS_CentralDifferenceSMS: return new CentralDifferenceSMS();` (parallel/database `recvSelf` reconstruction) | [#295](https://github.com/nmorabowen/OpenSees/pull/295) |
| `SRC/runtime/runtime/TclPackageClassBroker.cpp` | `// Ladruno`: `#include "CentralDifferenceSMS.h"` + `case INTEGRATOR_TAGS_CentralDifferenceSMS: return new CentralDifferenceSMS();` | [#295](https://github.com/nmorabowen/OpenSees/pull/295) |
| `SRC/analysis/integrator/CMakeLists.txt` | Add `CentralDifferenceSMS.cpp` to sources + `CentralDifferenceSMS.h`/`LadrunoMassLumping.h`/`LadrunoMassScaling.h` to headers | [#295](https://github.com/nmorabowen/OpenSees/pull/295) |
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature regen via `patch_banner.py` — add CentralDifferenceSMS + HRZ-lumping lines | [#295](https://github.com/nmorabowen/OpenSees/pull/295) |
