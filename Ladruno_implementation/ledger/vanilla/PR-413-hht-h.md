---
wp: PR-413
title: "413 -- 8 vanilla row(s)"
pr: "#413"
files: ["`SRC/analysis/integrator/HHT.h`", "`SRC/classTags.h`", "`SRC/tcl/commands.cpp`", "`SRC/interpreter/OpenSeesCommands.{h,cpp}`", "`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`", "`SRC/runtime/runtime/TclPackageClassBroker.cpp`", "`SRC/analysis/integrator/{CMakeLists.txt,Makefile}`", "`SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}`"]
table: "main"
legacy_seq: [188, 189, 190, 191, 192, 193, 194, 195]
---
| `SRC/analysis/integrator/HHT.h` | **`// Ladruno` ADR-52 W3-I2 (the key vanilla edit):** (a) promote the integrator state members (`alpha/beta/gamma/deltaT/c1-c3` + the `Ut/U/Ualpha` response vectors) from `private:` → `protected:` so the sensitivity subclass `LadrunoHHT` can reach them — pure access-level change, no algorithm edit; (b) add ONE protected **inline** ctor `HHT(int classTag,double,double,double)` so the subclass registers under `INTEGRATOR_TAGS_LadrunoHHT` while reusing the full HHT algorithm (mirrors the classTag-param ctor already in `Newmark.h`). Inline ⇒ **`HHT.cpp` stays byte-identical** (header-only vanilla edit). | [#413](https://github.com/nmorabowen/OpenSees/pull/413) |
| `SRC/classTags.h` | `// Ladruno` ADR-52 W3-I2: register `INTEGRATOR_TAGS_LadrunoHHT`=33013 (lowest free integrator tag; W3-I3's reserved 33013 was never built — NO-GO #410) | [#413](https://github.com/nmorabowen/OpenSees/pull/413) |
| `SRC/tcl/commands.cpp` | `// Ladruno` ADR-52 W3-I2: `extern OPS_LadrunoHHT()` + `integrator LadrunoHHT` strcmp branch in the legacy `specifyIntegrator()` Tcl parser (null-guard before `setIntegrator`); mirrors the `HHT` wiring | [#413](https://github.com/nmorabowen/OpenSees/pull/413) |
| `SRC/interpreter/OpenSeesCommands.{h,cpp}` | `// Ladruno` ADR-52 W3-I2: declare `OPS_LadrunoHHT()` + `strcmp(type,"LadrunoHHT")` dispatch in the openseespy/interpreter integrator command | [#413](https://github.com/nmorabowen/OpenSees/pull/413) |
| `SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp` | `// Ladruno` ADR-52 W3-I2: `#include "LadrunoHHT.h"` + `case INTEGRATOR_TAGS_LadrunoHHT: return new LadrunoHHT();` in BOTH `getNewTransientIntegrator` and `getNewIncrementalIntegrator` (mirrors Newmark, which carries DDM too) for parallel/database `recvSelf` reconstruction | [#413](https://github.com/nmorabowen/OpenSees/pull/413) |
| `SRC/runtime/runtime/TclPackageClassBroker.cpp` | `// Ladruno` ADR-52 W3-I2: `#include "LadrunoHHT.h"` + `case INTEGRATOR_TAGS_LadrunoHHT: return new LadrunoHHT();` in both integrator broker switches | [#413](https://github.com/nmorabowen/OpenSees/pull/413) |
| `SRC/analysis/integrator/{CMakeLists.txt,Makefile}` | `// Ladruno` ADR-52 W3-I2: add `LadrunoHHT.cpp`/`.h` to the integrator sources/headers (CMake) + `LadrunoHHT.o` (Makefile) | [#413](https://github.com/nmorabowen/OpenSees/pull/413) |
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature regen via `patch_banner.py` — add the `LadrunoHHT — sensitivity/DDM HHT` line | [#413](https://github.com/nmorabowen/OpenSees/pull/413) |
