---
wp: WP-149
title: "891 -- upstreamable-table row(s)"
pr: "#891"
files: ["`SRC/analysis/fe_ele/FE_Element.cpp`", "`SRC/analysis/fe_ele/transformation/TransformationFE.cpp`", "`SRC/analysis/dof_grp/TransformationDOF_Group.cpp`", "`SRC/analysis/dof_grp/DOF_Group.cpp`"]
table: "upstreamable"
legacy_seq: [667, 668, 669, 670]
---
| `SRC/analysis/fe_ele/FE_Element.cpp` | `// Ladruno WP-149`: both constructors' pool-init loops and the destructor's cleanup loop go `i<MAX_NUM_DOF` → `i<=MAX_NUM_DOF`. `theMatrices`/`theVectors` have `MAX_NUM_DOF+1` slots and the pooled branch is taken for `numDOF <= MAX_NUM_DOF`, so slot `[64]` was read UNINITIALIZED by any 64-DOF element (wild tangent/residual pointer when the heap left garbage there) and never freed. One token per loop, no behaviour change for any other DOF count. Upstream bug — requested by the owner (WP-149); do not file upstream without asking. | [#891](https://github.com/nmorabowen/OpenSees/pull/891) |
| `SRC/analysis/fe_ele/transformation/TransformationFE.cpp` | `// Ladruno WP-149`: same off-by-one in `modMatrices`/`modVectors` (`MAX_NUM_DOF` 64; ctor init loop + dtor cleanup loop → `<=`). `setID()` pools `numTransformedDOF <= MAX_NUM_DOF`, so a TransformationFE with exactly 64 transformed DOFs read slot `[64]` uninitialized. Fixed with FE_Element because under `constraints Transformation` the FE_Element fix alone would still leave a 64-DOF constrained element broken. | [#891](https://github.com/nmorabowen/OpenSees/pull/891) |
| `SRC/analysis/dof_grp/TransformationDOF_Group.cpp` | `// Ladruno WP-149`: same off-by-one (`MAX_NUM_DOF` 16; both ctors' init loops + dtor cleanup → `<=`). Reached by a transformed node with exactly 16 DOFs. Fixed by inspection (no reachable deck in the gate). | [#891](https://github.com/nmorabowen/OpenSees/pull/891) |
| `SRC/analysis/dof_grp/DOF_Group.cpp` | `// Ladruno WP-149`: same off-by-one (`MAX_NUM_DOF` 256; both ctors' init loops + dtor cleanup → `<=`). Reached only by a 256-DOF group (e.g. a Lagrange DOF_Group of that size). Fixed by inspection. | [#891](https://github.com/nmorabowen/OpenSees/pull/891) |
