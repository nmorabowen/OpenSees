---
wp: PR-623
title: "623 -- upstreamable-table row(s)"
pr: "#623"
files: ["`CMakeLists.txt`, `SRC/system_of_eqn/linearSOE/CMakeLists.txt`, `SRC/system_of_eqn/linearSOE/pardiso/CMakeLists.txt`", "`SRC/interpreter/OpenSeesCommands.{h,cpp}`"]
table: "upstreamable"
legacy_seq: [354, 355]
---
| `CMakeLists.txt`, `SRC/system_of_eqn/linearSOE/CMakeLists.txt`, `SRC/system_of_eqn/linearSOE/pardiso/CMakeLists.txt` | `// Ladruno ADR-75 P1b`: build+link MKL PARDISO for the **serial/desktop** targets. `add_subdirectory(pardiso)` gated on `MKL_FOUND` (non-MKL builds, e.g. Zone-A Ubuntu, never compile it and simply have no `system Pardiso`); new `PARDISO_FLAG`/`PARDISO_MKL_LIBRARIES` linked to `OpenSees` + `OpenSeesPy` only. **Threading layer is the crux:** parallel targets keep `mkl_sequential` (correct — their parallelism is MPI, and it is what ScaLAPACK/MUMPS use) while serial takes `mkl_intel_thread` + `libiomp5md`, because desktop PARDISO's entire value is shared-memory threading. Mixing layers in one process is unsupported by MKL, so they are kept on disjoint targets — **the MP/MUMPS lane is untouched**. | [#623](https://github.com/nmorabowen/OpenSees/pull/623) |
| `SRC/interpreter/OpenSeesCommands.{h,cpp}` | `// Ladruno ADR-75 P1b`: register `system Pardiso` (and `PARDISO`) + the `OPS_PARDISOGenLinSolver()` factory, under `#ifdef _PARDISO`. Note the factory uses a LOCAL SOE pointer: `theSOE` is a **member of `OpenSeesCommands`** (`OpenSeesCommands.h:178`), not a file-scope global, so a free factory cannot assign it. The serial `#else` branch of `OPS_MumpsSolver()` still does `theSOE = new MumpsSOE(...)` and carries that same undeclared-identifier bug — it survives only because `_MUMPS` is never defined without `_PARALLEL_INTERPRETERS`, i.e. a third independent proof (with the missing `libseq` and the CMake gating) that the **serial MUMPS path has never been compiled**. Left as-is: ADR-75 descopes serial MUMPS. | [#623](https://github.com/nmorabowen/OpenSees/pull/623) |
