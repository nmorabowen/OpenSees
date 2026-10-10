---
wp: LEGACY
title: "H5DRMLoadPattern::CalculateBoundaryForces had a fixed 8-node buffer (BoundaryNodes/ExteriorNodes) that silently overflows for any element with more than 8 node…"
legacy_seq: 381
---
### `H5DRMLoadPattern::CalculateBoundaryForces` had a fixed 8-node buffer (`BoundaryNodes`/`ExteriorNodes`) that silently overflows for any element with more than 8 nodes on one side of the DRM boundary
- **Bites:** any element with >8 total nodes used on an H5DRM boundary (`BezierTet10`=10, `LadrunoBrick20`/`Twenty_Node_Brick`=20, `TenNodeTetrahedron`=10). Undefined behavior, not a controlled error -- discovered by code audit, not by a crash report.
- **Why:** `constexpr int MaxNodes = 8` sized `BoundaryNodes`/`ExteriorNodes` once, OUTSIDE the per-element loop (`BoundaryNodes.resize(MaxNodes)`), then the loop wrote via `BoundaryNodes(boundaryCount++) = nodeIndex` with no bounds check. Every OTHER per-element buffer in the same function (`M_be`, `K_be`, `Peff_b`, `Peff_e`, `u_b`, `a_b`, `u_e`, `a_e`) was already correctly resized per-element based on `boundaryCount`/`exteriorCount` -- only these two were missed.
- **Workaround/status (2026-08-29, ADR-88 pre-implementation audit -- FIXED same PR):** moved the resize inside the per-element loop, sized to the element's own `numElementNodes` (not a fixed constant) -- `ID::resize()` reallocates safely when growing (`ID.cpp:394-422`), confirmed by reading the implementation before relying on it. | ADR-88 (h5drm higher-order elements)
