---
wp: PR-415
title: "415 -- 2 vanilla row(s)"
pr: "#415"
files: ["`SRC/analysis/integrator/GeneralizedAlpha.h`", "`SRC/{classTags.h, tcl/commands.cpp, interpreter/OpenSeesCommands.{cpp,h}, actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp, runtime/runtime/TclPackageClassBroker.cpp, analysis/integrator/CMakeLists.txt, analysis/integrator/Makefile, tcl/tclMain.cpp, interpreter/PythonModule.cpp}`"]
table: "main"
legacy_seq: [196, 197]
---
| `SRC/analysis/integrator/GeneralizedAlpha.h` | **`// Ladruno` ADR-52 W3-I2 PR2 (the key vanilla edit):** same header-only pattern as the `HHT.h` row above — promote the integrator state (`alphaM/alphaF/beta/gamma/deltaT/c1-c3` + the `Ut/U/Ualpha…` vectors incl. `Ualphadotdot`) `private:`→`protected:`, and add ONE protected **inline** classTag ctor `GeneralizedAlpha(int,double,double,double,double)` so `LadrunoGeneralizedAlpha` registers under `INTEGRATOR_TAGS_LadrunoGeneralizedAlpha` while reusing the full algorithm. Inline ⇒ **`GeneralizedAlpha.cpp` stays byte-identical**. | [#415](https://github.com/nmorabowen/OpenSees/pull/415) |
| `SRC/{classTags.h, tcl/commands.cpp, interpreter/OpenSeesCommands.{cpp,h}, actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp, runtime/runtime/TclPackageClassBroker.cpp, analysis/integrator/CMakeLists.txt, analysis/integrator/Makefile, tcl/tclMain.cpp, interpreter/PythonModule.cpp}` | `// Ladruno` ADR-52 W3-I2 PR2: extend the same integrator-registration set used for `LadrunoHHT` with the `LadrunoGeneralizedAlpha` sibling — classTag `INTEGRATOR_TAGS_LadrunoGeneralizedAlpha`=33014; `extern`+`strcmp` dispatch (Tcl `commands.cpp` + Py/interpreter `OpenSeesCommands.{cpp,h}`); `#include`+`case` in BOTH switches of each broker; sources/headers in CMake+Makefile; banner regen line. Each edit marked `// Ladruno`. | [#415](https://github.com/nmorabowen/OpenSees/pull/415) |
