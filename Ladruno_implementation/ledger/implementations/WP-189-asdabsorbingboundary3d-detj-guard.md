---
wp: WP-189
title: "ASDAbsorbingBoundary3D det J <= 0 guard"
date: 2026-10-14
status: "draft"
---
| **WP-189 — ASDAbsorbingBoundary3D refuses det J <= 0** — `element ASDAbsorbingBoundary3D` with a left-handed or degenerate node ordering (min signed det J over the 8 gauss points <= 0, on the nodes as given) now stops at `setDomain` with the element tag, the node tags and the det J, instead of integrating a negative free-field stiffness. Vertical boundary types only (a pure `B` element never uses det J). No behaviour change for valid elements. | Element guard | — (no new tag) | `SRC/element/absorbentBoundaries/ASDAbsorbingBoundary3D.cpp`, `tests/test_wp189_asd_detj_guard.py` (14 cases) | draft — PR open, not merged; valid decks byte-identical to the base binary, flipped decks refused (L, LB, R, F, LF), degenerate refused | — |
