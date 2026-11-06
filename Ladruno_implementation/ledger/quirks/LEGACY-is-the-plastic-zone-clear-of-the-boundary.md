---
wp: LEGACY
title: "\"Is the plastic zone clear of the boundary?\" cannot be tested as a fraction of the domain extent on a GRADED mesh"
legacy_seq: 259
---
### "Is the plastic zone clear of the boundary?" cannot be tested as a fraction of the domain extent on a GRADED mesh
- **Bites:** any field post-processing that asks whether a mechanism is
  contained. The obvious test — is the outermost yielded element's centroid
  within 90 % of the half-width? — is silently wrong when the mesh grades
  outward, because the far-field element is enormous. On the `hypo_bearing`
  benchmark the x-columns run `..., 4.00, 5.83, 8.45` in a 10 m half-domain:
  the LAST element is 3.1 m wide and its centroid sits at 8.45 m, so an element
  physically TOUCHING the roller passes a "< 9 m" test. Measured: the collapse
  leg reported "contained" while 16 of the 352 elements in the
  boundary-adjacent column were at m > 0.99.
- **Tell:** the reported extent equals one of the coarse outer centroids
  exactly, and equals it for several different legs — because it is a mesh
  coordinate, not a mechanism coordinate.
- **Rule:** test **element-column / row membership**, not distance. Build
  `np.unique(np.round(np.abs(centroid[:,0]), 4))` and compare against
  `cols.max()`; report `n_yielded / n_in_column` so the reader sees how much of
  the boundary layer is engaged rather than a yes/no. Same trap applies to any
  "did the wave reach the absorbing boundary" or "is the damage localized"
  check on a graded or adaptively refined mesh.
  *2026-07-30 (ADR-79 collapse study).*
