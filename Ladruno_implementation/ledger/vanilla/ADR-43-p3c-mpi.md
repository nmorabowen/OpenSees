---
wp: ADR-43
title: "ADR43 P3c-MPI -- 2 vanilla row(s)"
files: ["`SRC/system_of_eqn/eigenSOE/CMakeLists.txt`", "`CMakeLists.txt` (root)"]
table: "main"
legacy_seq: [266, 267]
---
| `SRC/system_of_eqn/eigenSOE/CMakeLists.txt` | `// Ladruno` ADR43 P3c: add `LadrunoFeastInnerSolve.{cpp,h}` (the RCI inner-solve seam + factory global) to `OPS_SysOfEqn` target_sources (build wiring only). | ADR43 P3c-MPI |
| `CMakeLists.txt` (root) | `# Ladruno` ADR43 P3c-MPI: add `LadrunoDistBlockZKernel.cpp` to the EXTRA_SOURCES of the two `_PARALLEL_INTERPRETERS` targets (`OpenSeesMP`, `OpenSeesPyMP`), mirroring `MumpsParallelSOE.cpp` — the distributed FEAST inner solve needs `_MUMPS`+`mpi.h`+`dmumps` on the link line, which `OPS_SysOfEqn` does NOT carry; deliberately NOT added to `OpenSeesSP` (`_PARALLEL_PROCESSING`, partitioned-domain model). Build wiring only. | ADR43 P3c-MPI |
