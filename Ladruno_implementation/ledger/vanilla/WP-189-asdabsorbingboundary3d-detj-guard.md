---
wp: WP-189
title: "ASDAbsorbingBoundary3D refuses det J <= 0"
date: 2026-10-14
files: ["`SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp`"]
table: "upstreamable"
---
| `SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp` | `// Ladruno WP-189`: `setDomain` — **a left-handed or degenerate node ordering is now refused at setup.** The free-field column (`addKff` / `addKffToSoil`) is integrated with the nodes in the user's order, so det J <= 0 at a gauss point gave a negative (or zero) volume weight and a silently negative stiffness (Paper 1 rev3 Tier 3: 161 left-handed skin elements per cell, unstable model, long diagnosis). New helper `H8_minDetJ` evaluates the minimum signed det J over the 8 gauss points on the unsorted `m_nodes`; if it is not > 0, `setDomain` prints the element tag, the nodes as given and the minimum det J through `opserr` and `exit(-1)` (the element's other setup refusals do the same). Checked for every boundary type containing a vertical face (L, R, F, K and corners); a pure BOTTOM element never uses det J and is left untouched. Inert for right-handed elements: byte-identical output against the base binary. | WP-189 |
