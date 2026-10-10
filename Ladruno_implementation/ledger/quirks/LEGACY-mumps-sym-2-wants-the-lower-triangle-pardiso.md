---
wp: LEGACY
title: "MUMPS SYM=2 wants the LOWER triangle; PARDISO mtype −2 wants the UPPER"
legacy_seq: 167
---
### MUMPS SYM=2 wants the LOWER triangle; PARDISO mtype −2 wants the UPPER

- When porting the block-real `(zM−K)` solve from serial PARDISO
  (`LadrunoBlockZKernel`, mtype −2, **upper** triangle CSR) to distributed MUMPS
  (`LadrunoDistBlockZKernel`, SYM=2), the stored triangle FLIPS: MUMPS symmetric
  expects entries with global **row ≥ col** (lower), matching how OpenSees's own
  `MumpsSOE`/`MumpsParallelSOE` assemble symmetric matrices (they store
  `row > vertexTag`). For the 2n block `[[aM−K,−bM],[−bM,−(aM−K)]]` the lower
  triangle is: `(i,j) j≤i` = aM−K; `(n+i, j)` for **all** j = the −bM block
  (row n+i ≥ n > j, wholly lower, supplied in full, its transpose block NOT
  supplied); `(n+i, n+j) j≤i` = −(aM−K). Verified transpose-consistent by the
  P3c-MPI adversarial gate.
