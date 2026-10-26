---
wp: WP-211
title: "modalProperties: MP-retained DOFs keep their mass under Transformation"
date: 2026-10-25
table: "main"
files: ["`SRC/domain/domain/DomainModalProperties.cpp`"]
---
| `SRC/domain/domain/DomainModalProperties.cpp` | `// Ladruno (WP-211)`: `node_map_t` marks a DOF `DOF_FIXED` only if its DOF_Group ID is -1 **and** the domain (or a load pattern) holds an SP constraint on that (node, dof); before, any -1 counted, and Transformation also gives -1 to MP-retained (slave) DOFs. A rigid diaphragm with its mass on the slave nodes therefore had `totalFreeMass` = 0, every Γ / Meff = 0, and `responseSpectrumAnalysis` returned zero. Plain/Penalty/Lagrange numbering is unaffected (their -1 DOFs are exactly the SP ones, or there are none). | — |
