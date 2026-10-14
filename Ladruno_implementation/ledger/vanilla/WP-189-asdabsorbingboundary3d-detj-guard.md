---
wp: WP-189
title: "ASDAbsorbingBoundary3D: node-order independent volume measure + degenerate-brick guard"
date: 2026-10-14
files: ["`SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp`"]
table: "upstreamable"
---
| `SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp` | `// Ladruno WP-189`: **the free-field column no longer depends on the handedness of the written node order.** `addKff` and `addRff` integrate over the H8 built from `m_nodes` in the written order with `dV = detJ * gw`; `dNdX = dN * inv(J)` is right for either handedness and the DOF mapping (`ffMapping`) and the mass use the coordinate sort (`m_node_map`), so the signed `detJ` was the only order-dependent term: a mirrored (left-handed) element integrated a NEGATIVE free-field stiffness and residual (Paper 1 rev3 Tier 3: 161 such skin elements per cell, unstable model; base binary: singular matrix at the first transient step). Now `dV = std::abs(detJ) * gw [* scale]`; positive-det elements are bitwise unchanged. `addKffToSoil` / `addRffToSoil` never used det J; `handleDistortion` works on the coordinate-sorted nodes, so order independent. **Guard:** `setDomain` (vertical boundary types only) calls new `H8_detJRange` and refuses a degenerate brick with `exit(-1)` and a message naming the tag, the nodes as given and min/max det J: det J changing sign between the 8 gauss points (folded / inverted) or min `|det J|` < 1e-6 max `|det J|` (collapsed). Either handedness of a sound brick is accepted. | WP-189 |
