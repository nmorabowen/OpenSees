---
wp: WP-260b
title: "WP-260b: packed 2-double Allreduce hook and preallocated MPI_Request/Status in MPIDiagonalSOE"
pr: "#1160"
date: 2026-11-26
files: ["SRC/system_of_eqn/linearSOE/LinearSOE.h", "SRC/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.h", "SRC/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.cpp"]
---
| `SRC/system_of_eqn/linearSOE/LinearSOE.h` | WP-260b: non-breaking virtual `globalReduceSum2(double, double, double&, double&) const` (default forwards to two `globalReduceSum` calls) so consistent-SMS parallel PCG packs paired reductions into one collective | #1160 |
| `SRC/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.h` | WP-260b: override `globalReduceSum2` and preallocate `sharedSumRequests`/`sharedSumStatuses` member vectors for `assembleSharedSum` | #1160 |
| `SRC/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.cpp` | WP-260b: reuse preallocated `MPI_Request`/`MPI_Status` scratch in `assembleSharedSum` (zero per-call heap allocs) and implement 2-double `MPI_Allreduce` in `globalReduceSum2` | #1160 |
