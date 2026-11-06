---
wp: PR-312
title: "312 -- 5 vanilla row(s)"
pr: "#312"
files: ["`SRC/interpreter/OpenSeesCommands.h`", "`SRC/interpreter/OpenSeesOutputCommands.cpp`", "`SRC/interpreter/PythonWrapper.cpp`", "`SRC/interpreter/TclWrapper.cpp`", "`SRC/domain/domain/Domain.cpp`"]
table: "main"
legacy_seq: [144, 145, 146, 147, 148]
---
| `SRC/interpreter/OpenSeesCommands.h` | `// Ladruno` (ADR-30 P3): declare `int OPS_LadrunoProjectionTieForce()` | [#312](https://github.com/nmorabowen/OpenSees/pull/312) |
| `SRC/interpreter/OpenSeesOutputCommands.cpp` | `// Ladruno` (ADR-30 P3): `OPS_LadrunoProjectionTieForce()` — the tie-force query `M(a_raw−a_proj)`; gets the active handler via `OPS_GetHandler()`, dynamic_casts to `LadrunoProjectionHandler`, returns `getTieForce(node,dof)` via `OPS_SetDoubleOutput` | [#312](https://github.com/nmorabowen/OpenSees/pull/312) |
| `SRC/interpreter/PythonWrapper.cpp` | `// Ladruno` (ADR-30 P3): register the `ladrunoProjectionTieForce` openseespy command (`Py_ops_LadrunoProjectionTieForce` → `OPS_LadrunoProjectionTieForce`) | [#312](https://github.com/nmorabowen/OpenSees/pull/312) |
| `SRC/interpreter/TclWrapper.cpp` | `// Ladruno` (ADR-30 P3): register the `ladrunoProjectionTieForce` Tcl command (`Tcl_ops_LadrunoProjectionTieForce` → `OPS_LadrunoProjectionTieForce`) | [#312](https://github.com/nmorabowen/OpenSees/pull/312) |
| `SRC/domain/domain/Domain.cpp` | `// Ladruno` (ADR-30 P3): `Domain::clearAll()` omitted `theEQs->clearAll()` — upstream added `EQ_Constraint` but never wired it into `clearAll`, so `wipe` LEAKED equation constraints into the next model. Invisible until a handler reads `getEQs()` (LadrunoProjection) → stale EQ mis-assembled. One-line additive fix mirroring `theMPs->clearAll()`. **Upstreamable bug.** | [#312](https://github.com/nmorabowen/OpenSees/pull/312) |
