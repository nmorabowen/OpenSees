---
wp: WP-132
title: "886 -- upstreamable-table row(s)"
pr: "#886"
files: ["`SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSolver.cpp`"]
table: "upstreamable"
legacy_seq: [673]
---
| `SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSolver.cpp` | `// Ladruno (pardiso-linux)`: the WP-132 `-cbwr` table entry `{"AVX10", MKL_CBWR_AVX10}` goes behind `#ifdef MKL_CBWR_AVX10`. The macro is oneMKL 2025.0+, so the Linux opt-in (`-DLADRUNO_MKL_PARDISO_LINUX=ON`, esmeralda oneMKL 2024.2) failed to compile from #864 on; Windows (2025.1) and Zone-A (no MKL) could not see it. The WP-132 header comment now says the "an earlier PARDISO solve does not block CNR" measurement is Windows-only. Fork-authored lines in a vanilla file; no behaviour change where the macro exists. | [#886](https://github.com/nmorabowen/OpenSees/pull/886) |
