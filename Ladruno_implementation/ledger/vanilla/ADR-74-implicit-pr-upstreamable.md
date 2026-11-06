---
wp: ADR-74
title: "ADR-74 implicit PR -- upstreamable-table row(s)"
files: ["`SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.cpp`"]
table: "upstreamable"
legacy_seq: [351]
---
| `SRC/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.cpp` | `// Ladruno` (ADR-74 implicit lane): `dc.s.*` profiler sub-brackets in `setSize` (`nnz/rest/alloc/rowfill/colfill`) + include — attribution-only, behavior-preserving. Measured verdict: implicit setSize is LINEAR (N^1.05 over the 0.25/0.5 M Newmark+Mumps rungs); the rowfill global-eq probe loop (O(neq_global × log V_local)/rank, np-invariant) is flagged in the ADR for cluster-scale re-measurement. | ADR-74 implicit PR |
