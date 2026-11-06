---
wp: PR-335
title: "335 -- 2 vanilla row(s)"
pr: "#335"
files: ["`SRC/system_of_eqn/linearSOE/LinearSOE.h`", "`SRC/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.{h,cpp}`"]
table: "main"
legacy_seq: [25, 26]
---
| `SRC/system_of_eqn/linearSOE/LinearSOE.h` | `// Ladruno` ADR-38 (V5, consistent parallel PCG): add four base virtuals defaulting to a serial no-op — `isDistributedDiagonal()`, `getScalingDiagonalA()`, `assembleSharedSum(Vector&)`, `globalReduceSum(double)`. Lets the shared-`OpenSeesLIB` consistent integrators drive a cross-rank PCG through the `LinearSOE*` base **without** referencing the MPI-only `MPIDiagonalSOE` (which is linked only into the MP executables, never the shared lib). Additive; every other SOE inherits the serial defaults unchanged. | [#335](https://github.com/nmorabowen/OpenSees/pull/335) |
| `SRC/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.{h,cpp}` | `// Ladruno` ADR-38 (V5): override the four hooks above — `isDistributedDiagonal()→true`, `getScalingDiagonalA()` returns the factored (1/mass) GLOBAL diagonal, `assembleSharedSum(Vector&)` sums a vector's shared-DOF entries across ranks (mirrors the solver's B-exchange), `globalReduceSum()` does `MPI_Allreduce(SUM)`. Enables the consistent (Olovsson) mass-scaling distributed PCG under `system MPIDiagonal`. MP-target-only file. | [#335](https://github.com/nmorabowen/OpenSees/pull/335) |
