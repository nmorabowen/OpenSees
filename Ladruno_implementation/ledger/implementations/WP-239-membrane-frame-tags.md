---
wp: WP-239
title: "LadrunoMembrane holds its supporting frames by tag; hinge-twin frame ends tie to the original nodes"
date: 2026-11-19
status: "draft"
section: "table"
---
| **WP-239 — `LadrunoMembrane` frames by tag, twin ends mapped to the original nodes** — (1) the membrane no longer caches an `Element*` per supporting frame: `addLoad` looks each frame up with `Domain::getElement(tag)` on every call, so a frame removed and re-created with the same tag (the hinge step of a nonlinear model) receives its slab share, and a frame removed for good is reported once (`supporting frame <t> is no longer in the domain`) instead of being written through a freed pointer. (2) A frame end on a node that is not an external node of the membrane but coincides with one (within `1e-6 L`) — a hinge twin tied by `equalDOF` / `zeroLength` — is matched to that external node, so the frame keeps its link line, its WP-209 end loads and its overhang-torque recipient, and no "lies on the panel but is not one of its nodes" warning is printed. A model restored from a database (which re-resolves against the twins) now equals the live one. No new state: `sendSelf`/`recvSelf` unchanged. Models that never remove frames and have no frame end coincident with a non-matching membrane node give bit-identical results. | element fix | — | `SRC/element/ladrunoMembrane/LadrunoMembrane.{h,cpp}`; `tests/_testbed/membrane_twin_model.py`, `tests/test_wp239_membrane_frame_replace.py` | draft | — |
