---
wp: LEGACY
title: "pardiso-linux -- upstreamable-table row(s)"
files: ["`CMakeLists.txt`, `SRC/system_of_eqn/linearSOE/CMakeLists.txt`"]
table: "upstreamable"
legacy_seq: [497]
---
| `CMakeLists.txt`, `SRC/system_of_eqn/linearSOE/CMakeLists.txt` | `// Ladruno (TIMs, 2026-09-03)` **`system Pardiso` was Windows-only by accident of the build, not of the solver.** ADR-75 P1b gates `PARDISO_FLAG` on `MKL_LPATH`, which is assigned only in the Windows branch from `.lib` names, so no MKL on Linux could enable it; the serial Linux build fell back to UmfPack, whose 32-bit index path fails (`numeric analysis returns -1`, 28 GB free) between 49 626 and 93 246 DOF on the TIMs `LadrunoUP` footing, whose 64-bit path never runs the symbolic analysis in this build, and where SuperLU took >18x UmfPack's whole gravity stage on the smaller mesh without finishing a solve. New opt-in `LADRUNO_MKL_PARDISO_LINUX` (default OFF, byte-identical otherwise) mirrors `LADRUNO_MKL_FEAST_LINUX`: explicit LP64/**sequential**/core 3-layer through `MKL_RT_HINT` plus `mkl_pardiso.h`, attached to `OPS_SysOfEqn` PUBLIC; the `pardiso/` subdirectory gate accepts it beside `MKL_FOUND`. Sequential layer on purpose — the Linux targets already link sequential MKL and one process carries one threading layer. Gate: on esmeralda, `mid` (49 626 DOF) `Pardiso` vs `UmfPack` ten push steps to 1e-6, then the 93 246-DOF wide box's elastic gravity stage solving. | pardiso-linux |
