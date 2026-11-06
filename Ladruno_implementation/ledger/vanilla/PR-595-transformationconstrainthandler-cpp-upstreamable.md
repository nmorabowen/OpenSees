---
wp: PR-595
title: "595 -- upstreamable-table row(s)"
pr: "#595"
files: ["`SRC/analysis/handler/TransformationConstraintHandler.cpp`", "`SRC/analysis/dof_grp/TransformationDOF_Group.cpp`"]
table: "upstreamable"
legacy_seq: [343, 344]
---
| `SRC/analysis/handler/TransformationConstraintHandler.cpp` | `// Ladruno` (ADR-74 handle fix): (a) O(1) `unordered_set`/`unordered_map` membership + first-index + per-node-SP-index structures replace every `ID::getLocation` linear scan in `handle()`'s list-build, node, element-classification (the measured N^1.94 dominant term), and FE-creation loops — searches only, every mutation/creation-order/branch quirk stock (gated byte-identical numbering + 18/18). `dc.handle` 14.40 → 0.34 s at the 2.0 M np8 rung. (b) `dc.h.*` profiler sub-brackets + includes. | [#595](https://github.com/nmorabowen/OpenSees/pull/595) |
| `SRC/analysis/dof_grp/TransformationDOF_Group.cpp` | `// Ladruno` (ADR-74 handle fix): remove the SP-only ctor's domain-wide SP sweep (O(#SP) per constrained node ⇒ O(#constrainedNodes × #SP) ~ N² per rank — the dc.h.nodes residual). Provably redundant: both in-tree callers follow the ctor with `addSP_Constraint` for every SP of the node from a superset list, setting the same `theSPs[dof]`; upstream already commented out the equivalent sweep in the MP ctor. Original kept in a comment; byte-identity + suite gated. | [#595](https://github.com/nmorabowen/OpenSees/pull/595) |
