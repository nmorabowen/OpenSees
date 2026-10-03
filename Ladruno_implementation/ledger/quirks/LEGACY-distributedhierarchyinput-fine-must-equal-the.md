---
wp: LEGACY
title: "DistributedHierarchyInput::fine must equal the MPI rank — it is not a free label"
legacy_seq: 237
---
### `DistributedHierarchyInput::fine` must equal the MPI rank — it is not a free label
- **Bites:** you try to express "rank r owns partition p" by setting `input.fine = p`, and `solveDistributedHierarchy` rejects the whole collective with `invalid distributed hierarchy input on at least one rank`. The validation is `input.fine != rank` (`LadrunoCMSHierarchy.cpp:1293`).
- **Consequence for testing:** a rank/partition permutation can only be expressed by **moving the data** (which subdomain's equations/stiffness/mass a rank carries), keeping `fine = rank`. That is also the physically honest formulation — it is what a different partitioner would hand you.
- **Workaround/status:** by design, documented here so the next agent does not read the rejection as a bug. *2026-07-26 (ADR-1000 P3d).*
