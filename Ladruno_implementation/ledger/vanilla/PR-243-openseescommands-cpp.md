---
wp: PR-243
title: "243 -- 4 vanilla row(s)"
pr: "#243"
files: ["`SRC/interpreter/OpenSeesCommands.cpp`", "`SRC/interpreter/OpenSeesCommands.h`", "`SRC/interpreter/PythonWrapper.cpp`", "`SRC/interpreter/TclWrapper.cpp`"]
table: "main"
legacy_seq: [140, 141, 142, 143]
---
| `SRC/interpreter/OpenSeesCommands.cpp` | `// Ladruno`: (a) `ladrunoArcLength` Layer-1.5 stabilization-energy subcommands `dissipationRatio`/`dissipatedEnergy`/`referenceEnergy`/`scaleCVisc` (ADR-31 rung-4 gate + R-RAMPDOWN actuator); (b) new `OPS_LadrunoDRCmd()` runtime query — reads `cmds->getTransientIntegrator()`, guards on `INTEGRATOR_TAGS_LadrunoDynamicRelaxation`, dispatches `residualNorm`/`kineticEnergy` (ADR-31 rung-5 settling signal); +`#include <LadrunoDynamicRelaxation.h>`. Additive. | [#243](https://github.com/nmorabowen/OpenSees/pull/243) |
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno`: declare `int OPS_LadrunoDRCmd()` | [#243](https://github.com/nmorabowen/OpenSees/pull/243) |
| `SRC/interpreter/PythonWrapper.cpp` | `// Ladruno`: register the `ladrunoDR` openseespy command (`Py_ops_ladrunoDR` → `OPS_LadrunoDRCmd`) — rung-5 DR settling/micro-burst query | [#243](https://github.com/nmorabowen/OpenSees/pull/243) |
| `SRC/interpreter/TclWrapper.cpp` | `// Ladruno`: register the `ladrunoDR` Tcl command (`Tcl_ops_ladrunoDR` → `OPS_LadrunoDRCmd`) | [#243](https://github.com/nmorabowen/OpenSees/pull/243) |
