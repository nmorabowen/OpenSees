---
wp: PR-593
title: "593 -- upstreamable-table row(s)"
pr: "#593"
files: ["`SRC/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.cpp`", "`SRC/analysis/analysis/DirectIntegrationAnalysis.cpp`"]
table: "upstreamable"
legacy_seq: [330, 331]
---
| `SRC/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.cpp` | `// Ladruno` (ADR-74 setSize fix): (a) `std::sort` replaces the hand-rolled leftmost-pivot `quickSort`/`q_sort` in `setSize` — the DOF graph's `std::map` storage hands the DOF list over ALREADY ASCENDING, so the old sort was a deterministic O((N/P)²/2) per rank per domainChanged (measured N^1.99, 95-98% of `dc.setSize`; 501→3.1 s at the 2.0 M np8 rung, 161×). Output identical (unique tags ⇒ same sorted array); `q_sort` kept for provenance; byte-identity + 18/18 suite gated. (b) `dc.s.*` profiler sub-brackets (`fill/sort/bcast/shared/localize/eleid`) + `<algorithm>`/ProfilerMacros includes; one hoisted declaration (`max`) for scope balance. | [#593](https://github.com/nmorabowen/OpenSees/pull/593) |
| `SRC/analysis/analysis/DirectIntegrationAnalysis.cpp` | `// Ladruno` (ADR-74 setSize): `dc.s.graph` sub-bracket around `getDOFGraph()` inside the existing `dc.setSize` scope, so DOF-graph construction is attributed separately from the SOE's own `setSize`. Behavior-preserving. | [#593](https://github.com/nmorabowen/OpenSees/pull/593) |
