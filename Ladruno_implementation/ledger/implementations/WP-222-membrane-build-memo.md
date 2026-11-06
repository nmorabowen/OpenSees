---
wp: WP-222
title: "LadrunoMembrane transfer build: bit-identical memo + allocation-free nearest-feature query"
pr: "#1061"
date: 2026-11-02
status: "draft"
section: "table"
---
| **WP-222 — `LadrunoMembrane` build speed, bit-identical** (YOO Guayaquil RevB, `apeSees.run` ~88 s): (a) memo of the pure stages of the eager transfer build (validation mesh, `distribute_surface_load`, diaphragm mesh + `condense_lst`). The key is the exact bit image of each stage's inputs, and the whole key is stored, so a collision is impossible. The memo is emptied when the last membrane is destroyed (`wipe`) and has a 256 MB budget per stage. **`-noMemo`** builds a panel cold; **`transferMemo`** response = [validation, distribution, diaphragm hits, entries]. (b) Hypot-screened, allocation-free nearest-feature query and allocation-free feature sets in the `tributary` recursion; `-DLADRUNO_MEMBRANE_KERNEL_REFERENCE` keeps the old path for A/B. Captured RevB T2 (2795 membranes) replayed through the pyd: outputs and `K u`, `M a` bytewise identical to base `d56fbfc48` (38.4 M values; T1 9.6 M), 80.9 s → 5.8 s, (b) alone 44.4 s (quiet box; two more paired rounds under other sessions' load: 131.7 → 8.8 s and 205.5 → 10.9 s, so 14–19×) | element (follow-up to ADR-98 / WP-204) | `ELE_TAG_LadrunoMembrane` 33022 (unchanged) | `SRC/element/ladrunoMembrane/{LadrunoMembraneMemo.h (new), LadrunoMembraneKernel.h, LadrunoMembrane.h, LadrunoMembrane.cpp, OPS_LadrunoMembrane.cpp, CMakeLists.txt}`, `tests/test_wp222_membrane_build_memo.py` | draft | [#1061](https://github.com/nmorabowen/OpenSees/pull/1061) |
