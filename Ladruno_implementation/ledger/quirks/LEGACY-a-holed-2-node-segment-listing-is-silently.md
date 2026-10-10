---
wp: LEGACY
title: "A \"holed\" 2-node segment listing is silently legal -- -slave-segments 2 / -master 2 take a FLAT STRIDE-2 PAIR LIST, not a node chain"
legacy_seq: 314
---
### A "holed" 2-node segment listing is silently legal -- `-slave-segments 2` / `-master 2` take a FLAT STRIDE-2 PAIR LIST, not a node chain
- **Bites:** anyone declaring a multi-segment 2D contact surface from a bare node list (the natural-looking `contactSurface(20, "-slave-segments", 2, n0, n1, n2, n3)` declares TWO DISJOINT segments (n0,n1),(n2,n3) with a hole where (n1,n2) should be -- 3 chained segments need SIX tags: `n0,n1, n1,n2, n2,n3`). The T3 gate authoring hit this: the holed patch deck CONVERGED, BALANCED its reactions and transferred the load through the wrong distribution (master row [5/18, 8/18, 5/18]P instead of [1/4, 1/2, 1/4]P) -- the exact ADR-78 P0 shape, and the parser cannot object (an even tag count with no shared nodes is indistinguishable from intentional disjoint contact patches, which are legitimate).
- **Why:** segments are flat nps-blocks (the 3D facet convention collapsed to nps=2); the NTS chain-integrity scan only FATALs on MIS-ORDERED shared nodes, and a holed list shares none. Root-caused via the oracle: the "wrong" forces were the EXACT discrete solution of the holed surface as declared (kernel and handler exonerated by an MSVC parity driver + hand solve).
- **Workaround/status (2026-08-18, ADR-85 T3):** convention documented here and owed to the T4 user guide in bold; test decks build pair lists with an explicit `_seg_pairs()` helper. No guard is possible without forbidding legitimate disjoint patches.
