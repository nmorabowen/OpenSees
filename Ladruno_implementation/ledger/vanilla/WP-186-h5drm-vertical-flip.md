---
wp: WP-186
title: "H5DRM vertical sign flip restored (ShakerMaker z-down rows -> z-up model)"
date: 2026-10-11
pr: ""
files: ["`SRC/domain/pattern/drm/H5DRMLoadPattern.cpp`"]
table: "main"
---
| `SRC/domain/pattern/drm/H5DRMLoadPattern.cpp` | `// Ladruno WP-186`: uncomment `d1[2] = -d1[2]; d2[2] = -d2[2]; a1[2] = -a1[2]; a2[2] = -a2[2];` in `drm_direct_read` and `drm_differentiate_displacements` (commented out by #296), and negate `d_final[2]` in the hold-past-tend branch of `drm_direct_read` (#296 / jaabell 99b7c2c11) so the held vertical does not jump sign at `tend`. Convention, now stated in the source: rows (E, N, vertical) with vertical positive DOWN (ShakerMaker); T/x0 act on the station coordinates only; the model gets u_X = +d0, u_Y = +d1, u_Z = -d2, i.e. correct for a z-up model (decks' T = (0 1 0)(1 0 0)(0 0 -1)). Author instruction (Patricio, ShakerMaker/DRM, 2026-10-11). Upstream status: OpenSees/OpenSees master has the flip active in both read paths (as of 2026-10-11, last touched 037d5d947); jaabell `merged-Explicit-ASDP-master` (60248752d) has it commented in `drm_direct_read` and active in `drm_differentiate_displacements`; this row brings the fork back to upstream behaviour plus the hold-branch fix (upstream has no hold branch). Test: `tests/test_h5drm_vertical_sign.py`. | WP-186 |
