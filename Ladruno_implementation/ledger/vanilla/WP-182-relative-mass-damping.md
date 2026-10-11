---
wp: WP-182
title: "relativeMassDamping hooks in Node, DOF_Group, Domain and the three command registries"
date: 2026-10-10
files: ["`SRC/domain/node/Node.h`", "`SRC/domain/node/Node.cpp`", "`SRC/domain/node/CMakeLists.txt`", "`SRC/domain/domain/Domain.cpp`", "`SRC/analysis/dof_grp/DOF_Group.cpp`", "`SRC/tcl/commands.cpp`", "`SRC/tcl/commands.h`", "`SRC/interpreter/TclWrapper.cpp`", "`SRC/interpreter/PythonWrapper.cpp`"]
---
| `SRC/domain/node/Node.h` | `// Ladruno (relativeMassDamping)`: forward decl, `setMassDampingReference` / `getMassDampingReference` / `getRayleighDampingFactor` / `getMassDampingRelativeVel`, two in-class-initialized members (non-owning reference pointer, scratch vector). | PR pending |
| `SRC/domain/node/Node.cpp` | `// Ladruno (relativeMassDamping)`: guarded relative branch (`alphaM M (v - v_ref)`) in `getUnbalancedLoadIncInertia`, `resetReactionForce` (flag 0/2 branch) and `getResponse(RayleighForces)`; setter and `getMassDampingRelativeVel` (domain time); scratch freed in the destructor. Original statements kept verbatim in the unguarded branch. | PR pending |
| `SRC/domain/node/CMakeLists.txt` | `RelativeMassDamping.{cpp,h}` added to `OPS_Domain`. | PR pending |
| `SRC/domain/domain/Domain.cpp` | `// Ladruno (relativeMassDamping)`: `clearAll()` empties the reference registry (`OPS_clearAllRelativeMassDampings`) before the nodes go. | PR pending |
| `SRC/analysis/dof_grp/DOF_Group.cpp` | `// Ladruno (relativeMassDamping)`: `addD_Force` refuses (FATAL + `_Exit(-1)`) a node carrying a reference with non-zero alphaM — the `*_TP` / sensitivity path would drop the reference silently. | PR pending |
| `SRC/tcl/commands.cpp` | `// Ladruno (relativeMassDamping)`: classic-Tcl bridge `relativeMassDampingCmd` (OPS_ResetInputNoBuilder → `OPS_RelativeMassDamping`) + registration under `OPS_RelativeMassDampingCmdName`. | PR pending |
| `SRC/tcl/commands.h` | `// Ladruno (relativeMassDamping)`: declaration of `relativeMassDampingCmd`. | PR pending |
| `SRC/interpreter/TclWrapper.cpp` | `// Ladruno (relativeMassDamping)`: interpreter-Tcl registration of `relativeMassDamping`. | PR pending |
| `SRC/interpreter/PythonWrapper.cpp` | `// Ladruno (relativeMassDamping)`: openseespy registration of `relativeMassDamping`. | PR pending |
