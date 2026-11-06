---
wp: PR-493
title: "493 -- 1 vanilla row(s)"
pr: "#493"
files: ["`SRC/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSolver.{h,cpp}`"]
table: "main"
legacy_seq: [228]
---
| `SRC/system_of_eqn/linearSOE/umfGEN/UmfpackGenLinSolver.{h,cpp}` | `// Ladruno` ADR-40 rank 2: default `UMFPACK_STRATEGY` **AUTO** (was hardcoded `SYMMETRIC` in both di/dl `setSize` branches — mis-orders unsymmetric tangents, e.g. `LadrunoConcrete3D`); expose `-strategy auto\|symmetric\|unsymmetric` + `-pivotTol <f>` (default keeps legacy 1.0 = maximal-threshold pivoting). Legacy path = `-strategy symmetric -pivotTol 1.0`, verified bit-identical on the Phase-0 lane-B bench. Measured on 40b lane B (66.4% linearSolve): see ADR-40 log. | [#493](https://github.com/nmorabowen/OpenSees/pull/493) |
