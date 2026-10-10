---
wp: PR-564
title: "564 -- 1 vanilla row(s)"
pr: "#564"
files: ["`SRC/recorder/GmshRecorder.{h,cpp}`"]
table: "main"
legacy_seq: [297]
---
| `SRC/recorder/GmshRecorder.{h,cpp}` | `// Ladruno` hex20 MSH fix: `ELE_TAG_Twenty_Node_Brick` was mapped to MSH type 12 (the **27-node** hex per the header's own MSH table) — the 20-node serendipity hex is type **17**. Added `GMSH_HEXAHEDRON_20=17` to the `GmshType` enum, remapped `Twenty_Node_Brick`, added `ELE_TAG_LadrunoBrick20` (33018, identical shp3dv node ordering), and reorder the connectivity write for type-17 elements only: corners 0-7 pass through, mid-edges permuted `[8,11,16,9,17,10,18,19,12,15,13,14]` (OpenSees edge order (0-1)(1-2)(2-3)(3-0)(4-5)(5-6)(6-7)(7-4)(0-4)(1-5)(2-6)(3-7) → Gmsh edge order (0-1)(0-3)(0-4)(1-2)(1-5)(2-3)(2-6)(3-7)(4-5)(4-7)(5-6)(6-7), Gmsh ref. manual "Node ordering"). All other element types keep the raw `getExternalNodes()` write path. Test: `tests/test_gmsh_hex20_recorder.py`. | [#564](https://github.com/nmorabowen/OpenSees/pull/564) |
