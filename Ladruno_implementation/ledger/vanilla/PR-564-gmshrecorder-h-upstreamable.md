---
wp: PR-564
title: "564 -- upstreamable-table row(s)"
pr: "#564"
files: ["`SRC/recorder/GmshRecorder.{h,cpp}`"]
table: "upstreamable"
legacy_seq: [316]
---
| `SRC/recorder/GmshRecorder.{h,cpp}` | 20-node hex written as MSH type 12 (27-node hex) with unpermuted node order — Gmsh cannot read the mesh. Fix: type 17 + mid-edge permutation (Twenty_Node_Brick part is upstreamable as-is; the sibling type-12 mappings `TwentyNodeBrick`/`TwentyNodeBrick_u_p_U`/`TotalLagrangianFD20NodeBrick` still carry the same bug upstream). | [#564](https://github.com/nmorabowen/OpenSees/pull/564) |
