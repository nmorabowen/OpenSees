---
wp: PR-598
title: "598 -- upstreamable-table row(s)"
pr: "#598"
files: ["`SRC/analysis/numberer/DOF_Numberer.cpp`", "`SRC/analysis/numberer/PlainNumberer.cpp`", "`SRC/analysis/handler/PlainHandler.cpp`"]
table: "upstreamable"
legacy_seq: [347, 348, 349]
---
| `SRC/analysis/numberer/DOF_Numberer.cpp` | `// Ladruno` (ADR-74 MP-index): the `-4` fixup in both `numberDOF` variants iterated ALL MP_Constraints per constrained DOF group — O(#groups × #MP), quadratic on equalDOF/tie-heavy decks. Replaced with a one-pass `constrainedNode → vector<MP*>` index (push_back preserves `getMPs()` order per node ⇒ stock multi-constraint application order). Byte-identical numbering on the tie gate (`~/ladruno_nsweep/tiegate`) + 18/18 suite. | [#598](https://github.com/nmorabowen/OpenSees/pull/598) |
| `SRC/analysis/numberer/PlainNumberer.cpp` | `// Ladruno` (ADR-74 MP-index): same `-4` fixup index in both variants. Same gates. | [#598](https://github.com/nmorabowen/OpenSees/pull/598) |
| `SRC/analysis/handler/PlainHandler.cpp` | `// Ladruno` (ADR-74 MP-index): the per-NODE full sweeps of `getMPs()` AND `getEQs()` (O(#nodes × #MP/#EQ)) replaced with `allMPs`/`allEQs` multimaps — the exact pattern the file's own `allSPs` already uses for SPs (multimap preserves insertion order for equal keys ⇒ stock visit order). Same gates. | [#598](https://github.com/nmorabowen/OpenSees/pull/598) |
