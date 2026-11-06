---
wp: WP-149
title: "Analysis-object pools: slot [MAX_NUM_DOF] left uninitialized (WP-149)"
pr: "#891"
status: "fixed on the branch"
section: "table"
legacy_seq: 4
---
| **Analysis-object pools: slot `[MAX_NUM_DOF]` left uninitialized (WP-149)** — vanilla `FE_Element`, `TransformationFE`, `DOF_Group`, `TransformationDOF_Group` keep class-wide `Matrix*`/`Vector*` arrays of `MAX_NUM_DOF+1` slots indexed by DOF count and use the POOLED branch for `numDOF <= MAX_NUM_DOF`, but zeroed (and freed) only `i < MAX_NUM_DOF`. An object with exactly `MAX_NUM_DOF` DOFs (64 for the two FE classes — reachable by an ordinary coupling element) took whatever the heap left in the last slot: nonzero garbage = its tangent/residual at a wild address; zero = works, and the slot's Vector+Matrix leaked at teardown. Fix: `<=` at all 10 init/cleanup loops. Found by the material-point parallelism scoping (A3 audit of the FE_Element pool). Gate `tests/test_wp149_pool_slot_max_num_dof.py`: a 64-DOF RBE3 under `constraints Plain` and under `constraints Transformation` (an `equalDOF`-tied independent keeps the transformed count at 64), the pool array carved from a DELIBERATELY DIRTIED free list, 4 build→analyze→wipe cycles per deck in a child interpreter, RBE3 distribution checked exactly. The DOF_Group (256) / TransformationDOF_Group (16) slots are fixed by inspection only. | upstream bug fix | — (vanilla analysis objects) | `SRC/analysis/fe_ele/FE_Element.cpp`, `SRC/analysis/fe_ele/transformation/TransformationFE.cpp`, `SRC/analysis/dof_grp/{DOF_Group,TransformationDOF_Group}.cpp`, `tests/test_wp149_pool_slot_max_num_dof.py` | **fixed on the branch** | [#891](https://github.com/nmorabowen/OpenSees/pull/891) |
