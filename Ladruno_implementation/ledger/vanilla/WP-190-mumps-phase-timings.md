---
wp: WP-190
title: "MumpsParallelSolver phase timings"
date: 2026-10-14
pr: "#TBD"
files: ["`SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.h`", "`SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.cpp`"]
---
| `SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.cpp` | `// WP-190`: with `-stats`, MPI_Wtime around JOB=1 / JOB=2 / JOB=3 (JOB=5 split into 2+3 only when stats are on), max-over-ranks `t_analysis`, `t_factor`, `t_solve (n)` lines appended to the existing `MUMPS stats:` block; `printStats` added to `sendSelf`/`recvSelf` (ID 3 -> 4). Stats off: unchanged JOB=5 path. Not filed upstream. | #TBD |
| `SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.h` | `// WP-190`: `tAnalysis`, `tSolveSum`, `nSolveSum`, `reportSolveTail()`. Not filed upstream. | #TBD |
