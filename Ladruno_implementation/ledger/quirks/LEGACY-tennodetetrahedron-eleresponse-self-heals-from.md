---
wp: LEGACY
title: "TenNodeTetrahedron::eleResponse self-heals from nodal trial displacement on every query"
legacy_seq: 395
---
### `TenNodeTetrahedron::eleResponse` self-heals from nodal trial displacement on every query
- **Bites:** "stresses"/"forces"/"material" responses always re-derive strain from the CURRENT nodal trial displacement and re-run the material, so a material-level Trial-state corruption (e.g. ASDP's no-op `revertToLastCommit`) is invisible through any eleResponse path after a domain-level revert; `ops.reset()` then reports a third stress value that is neither zero nor the pre-reset commit.
- **Rule:** to observe raw material Trial/Commit state after a revert, use a source-level structural pin or a recorder that reads the material directly; do not conclude "fixed" from a tet eleResponse.
