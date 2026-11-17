---
wp: WP-234
title: "region -node ... -ele ... keeps only the element nodes"
pr: "#1104"
date: 2026-11-16
---
### `region $t -node ... -ele ...` silently keeps only the elements' nodes (`both elements & nodes set, ONLY set using elements`)

- **Bites:** any deck that puts `-node` and `-ele` on one `region` line to restrict a
  recorder to "these elements plus these other nodes". `MeshRegion::setNodes` and
  `MeshRegion::setElements` (`SRC/domain/region/MeshRegion.cpp`) each delete BOTH lists,
  and with `-ele` the node list is rebuilt from the elements' connectivity. The only
  signal is the warning `WARNING region - both elements & nodes set, ONLY set using
  elements`; the explicit nodes are gone. Emitting `-ele` before `-node` does not help:
  `setNodes` then wipes the elements.
- **Seen in:** apeGmsh `ops.recorder.Ladruno(nodes_pg=..., elements_pg=...)` on an
  embedded-rebar beam (apeGmsh_Examples 02): the 432 bar and stirrup nodes were never
  recorded and the viewer drew the cage undeformed.
- **Do instead:** two regions, one per list, and `recorder ladruno ... -R $eleRegion
  -RN $nodeRegion` (WP-234). MPCO and the stock recorders have no `-RN`: give them a
  node-only region and an element-only region in two recorders.
