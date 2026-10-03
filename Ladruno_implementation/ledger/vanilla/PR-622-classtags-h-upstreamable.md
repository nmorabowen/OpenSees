---
wp: PR-622
title: "622 -- upstreamable-table row(s)"
pr: "#622"
files: ["`SRC/classTags.h`", "`SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSolver.{h,cpp}`"]
table: "upstreamable"
legacy_seq: [352, 353]
---
| `SRC/classTags.h` | `// Ladruno ADR-75 P1`: define `SOLVER_TAGS_PARDISOGenLinSolver`=33000 (solver registry band). **Upstream bug fix** — the tag is *used* at `PARDISOGenLinSolver.cpp:27` but was never defined anywhere upstream (upstream added `LinSOE_TAGS_PARDISOGenLinSOE`=99990 but not the solver tag), which is why that contributed 2019 MKL-PARDISO prototype has never compiled. Placed in the fork's private band rather than next to 99990 so a future upstream definition cannot collide; upstream's SOE tag left untouched. | [#622](https://github.com/nmorabowen/OpenSees/pull/622) |
| `SRC/system_of_eqn/linearSOE/pardiso/PARDISOGenLinSolver.{h,cpp}` | `// Ladruno ADR-75 P1`: **factorization reuse**. Stock ran PARDISO phases 11→22→33→−1 on EVERY solve with `pt[]`/`iparm` as function locals (redoing the METIS reorder *and* the numeric factorization per solve, then freeing them; `iparm` leaked). Now: handle/control array are members; phase 11 once per sparsity pattern (`setSize`), phase 22 only when `A` changed (gated on the SOE `factored` flag, mirroring `MumpsSolver` job=5/job=3), phase 33 per call, phase −1 once in the dtor. Also `iparm[2]=1`→0 (does not pin threads; `MKL_NUM_THREADS` does), `iparm[18]`→0 (Intel documents −1 as increasing reordering time; never printed at `msglvl=0`), MKL error-code decoding. Adversarial review fixed 3 defects: dtor use-after-free (`~LinearSOE` deletes the solver *after* the SOE freed its arrays ⇒ cache `cachedN`, pass dummies), a virgin-handle path (`needsSymbolic \|\| !init`), and the iparm[18] cost. Compile-verified; not yet built/registered. | [#622](https://github.com/nmorabowen/OpenSees/pull/622) |
