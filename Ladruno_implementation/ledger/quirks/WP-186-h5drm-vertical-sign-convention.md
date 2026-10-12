---
wp: WP-186
title: "H5DRM: T rotates station coordinates, not motions; vertical row is negated"
date: 2026-10-11
pr: ""
---
### H5DRM vertical motion upside down / which sign does the .h5drm vertical row need

`pattern H5DRM ... T(9) x0(3)` applies `xyz_model = T (xyz_file - drmbox_x0) * crd_scale + x0` to the
station **coordinates only**. The motion rows are never rotated: the model receives
u_X = +d0, u_Y = +d1, u_Z = **-d2** (the flip, restored by WP-186 after #296 had commented it out).
ShakerMaker writes rows (E, N, vertical) with the vertical positive DOWN (2.0 writer `ve, vn, vz`), and
station xyz (N, E, depth). So with the usual `T = 0 1 0 1 0 0 0 0 -1` (model X = E, Y = N, Z up) a native
ShakerMaker file enters with the right sign and needs no edit. Two cases need row 2 negated **in the file**:
a model built in the z-down frame itself (identity T; this was the #296 drm_study setup, which is why it
"validated" the removal), and a file whose row 2 is positive UP (ShakerMaker 1.1 / the Paper 1 Legacy
`legacy_cut40_rest_drm` set, built for the no-flip binary). Horizontal rows are not touched by T either:
a T that is not the E/N swap above needs the rows rearranged in the file. Binaries without the flip:
fork #296 (29a717041) up to WP-186, incl. pins fd87e396d, c4db8b21c, `pin/paper1-single` 1a4141ee4, and
José's ASDPlastic build of 2025-12-18 (flip commented in `drm_direct_read`). Sign test:
`tests/test_h5drm_vertical_sign.py` (rigid-translation free field, exact).
