---
wp: WP-181
title: "PARDISOGenLinSolver.cpp: legacy CBWR branches behind #ifdef for oneMKL 2026.1"
date: 2026-10-09
files: ["`SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSolver.cpp`"]
table: "upstreamable"
---
| `SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSolver.cpp` | `// Ladruno WP-181`: the WP-132 `-cbwr` table entries `SSE3`, `SSSE3`, `SSE4_1`, `AVX`, `AVX512_MIC` and `AVX512_MIC_E1` each go behind `#ifdef MKL_CBWR_<name>`. oneMKL 2026.1 removed those six macros from `mkl_types.h` (it keeps `OFF`, `BRANCH_OFF`, `AUTO`, `COMPATIBLE`, `SSE2`, `SSE4_2`, `AVX2`, `AVX512`, `AVX512_E1`, `AVX10`), so the Windows build stopped compiling once `mkl\latest` pointed at 2026.1. Same pattern as #886's `AVX10` guard, from the other end of the version range. With 2026.1 those six names report an unknown branch; no behaviour change where the macros exist. | WP-181 |
