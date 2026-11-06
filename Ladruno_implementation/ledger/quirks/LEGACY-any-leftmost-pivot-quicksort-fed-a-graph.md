---
wp: LEGACY
title: "Any leftmost-pivot quicksort fed a graph-derived list is deterministically O(n²) — MapOfTaggedObjects iterates ASCENDING, so the input is already sorted"
legacy_seq: 202
---
### Any leftmost-pivot quicksort fed a graph-derived list is deterministically O(n²) — MapOfTaggedObjects iterates ASCENDING, so the input is already sorted
- **Bites:** `MPIDiagonalSOE::setSize` sorted its DOF list with a hand-rolled `q_sort` (leftmost pivot). The DOF graph is `std::map`-backed, its iterator hands tags over ascending ⇒ the sort's textbook worst case, per rank, every domainChanged: 501 s at the 2.0 M np8 rung, ~N^1.9 — masqueraded for years as "setup cost". The same trap arms ANY hand-rolled pivot sort downstream of a `MapOfTaggedObjects`/`ArrayOfTaggedObjects` iteration (both are ascending for dense tags).
- **Why:** classic quicksort worst case = sorted input + first-element pivot; map iteration guarantees sorted input.
- **Workaround/status:** FIXED in-place with `std::sort` ([#593](https://github.com/nmorabowen/OpenSees/pull/593), output provably identical). Rule: never hand-roll a sort on tag/equation lists; grep found no other `q_sort` copies. *2026-07-22.*
