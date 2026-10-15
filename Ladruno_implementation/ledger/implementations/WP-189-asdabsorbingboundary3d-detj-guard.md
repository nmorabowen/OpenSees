---
wp: WP-189
title: "ASDAbsorbingBoundary3D node-order independence + degenerate-brick guard"
date: 2026-10-14
status: "draft"
---
| **WP-189 — ASDAbsorbingBoundary3D is independent of node-order handedness, and refuses a degenerate brick** — the free-field stiffness/residual use `|det J|`, so a mirrored (left-handed) write order gives the response of its right-handed twin (the Paper 1 rev3 Tier-3 skin defect no longer needs a deck fix); `setDomain` refuses a brick whose det J changes sign between gauss points or collapses (min `|det J|` < 1e-6 max), printing tag, nodes and det J range. Vertical boundary types only. Right-handed elements bitwise unchanged. | Element fix + guard | — (no new tag) | `SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp`, `tests/test_wp189_asd_detj_guard.py` (23 cases) | draft — PR open, not merged; mirrored = original to round-off for L, LB, R, F, K, LF, LK, RF, RK, LFB (base binary fails: singular); right-handed byte-identical to the base binary; collapsed and folded refused | — |
