---
wp: WP-182
title: "relativeMassDamping hooks in Node, DOF_Group, Domain and the three command registries"
date: 2026-10-10
files: ["`SRC/domain/node/Node.h`", "`SRC/domain/node/Node.cpp`", "`SRC/domain/node/CMakeLists.txt`", "`SRC/domain/domain/Domain.cpp`", "`SRC/analysis/dof_grp/DOF_Group.cpp`", "`SRC/tcl/commands.cpp`", "`SRC/tcl/commands.h`", "`SRC/interpreter/TclWrapper.cpp`", "`SRC/interpreter/PythonWrapper.cpp`"]
pr: "#968"
---
| `SRC/domain/node/Node.h` | `// Ladruno (relativeMassDamping)`: forward decl, `setMassDampingReference` / `getMassDampingReference` / `getRayleighDampingFactor` / `getMassDampingRelativeVel`, two in-class-initialized members (non-owning reference pointer, scratch vector). | [#968](https://github.com/nmorabowen/OpenSees/pull/968) |
| `SRC/domain/node/Node.cpp` | `// Ladruno (relativeMassDamping)`: guarded relative branch (`alphaM M (v - v_ref)`) in `getUnbalancedLoadIncInertia`, `resetReactionForce` (flag 0/2 branch) and `getResponse(RayleighForces)`; setter and `getMassDampingRelativeVel` (domain time); scratch freed in the destructor. Original statements kept verbatim in the unguarded branch. | [#968](https://github.com/nmorabowen/OpenSees/pull/968) |
| `SRC/domain/node/CMakeLists.txt` | `RelativeMassDamping.{cpp,h}` added to `OPS_Domain`. | [#968](https://github.com/nmorabowen/OpenSees/pull/968) |
| `SRC/domain/domain/Domain.cpp` | `// Ladruno (relativeMassDamping)`: `clearAll()` empties the reference registry (`OPS_clearAllRelativeMassDampings`) before the nodes go. | [#968](https://github.com/nmorabowen/OpenSees/pull/968) |
| `SRC/analysis/dof_grp/DOF_Group.cpp` | `// Ladruno (relativeMassDamping)`: `addD_Force` refuses (FATAL + `_Exit(-1)`) a node carrying a reference with non-zero alphaM — the `*_TP` / sensitivity path would drop the reference silently. | [#968](https://github.com/nmorabowen/OpenSees/pull/968) |
| `SRC/tcl/commands.cpp` | `// Ladruno (relativeMassDamping)`: classic-Tcl bridge `relativeMassDampingCmd` (OPS_ResetInputNoBuilder → `OPS_RelativeMassDamping`) + registration under `OPS_RelativeMassDampingCmdName`. Follow-up (WP-183): the `#include <RelativeMassDamping.h>` sat inside the `#else` of `_PARALLEL_PROCESSING`, so OpenSeesSP (and the full `build.bat`) failed with C2065 `OPS_RelativeMassDampingCmdName`; it now follows the `#endif`. | [#968](https://github.com/nmorabowen/OpenSees/pull/968), WP-183 |
| `SRC/tcl/commands.h` | `// Ladruno (relativeMassDamping)`: declaration of `relativeMassDampingCmd`. | [#968](https://github.com/nmorabowen/OpenSees/pull/968) |
| `SRC/interpreter/TclWrapper.cpp` | `// Ladruno (relativeMassDamping)`: interpreter-Tcl registration of `relativeMassDamping`. | [#968](https://github.com/nmorabowen/OpenSees/pull/968) |
| `SRC/interpreter/PythonWrapper.cpp` | `// Ladruno (relativeMassDamping)`: openseespy registration of `relativeMassDamping`. | [#968](https://github.com/nmorabowen/OpenSees/pull/968) |
