---
wp: WP-133
title: "WP-133 (TIMs F23) — PDMY03 critical-state constants as optional flags + dilation-brake note"
pr: "#866"
status: "draft PR open"
class_tags: ["ND_TAG_PressureDependMultiYield03"]
section: "table"
legacy_seq: 185
---
| **WP-133 (TIMs F23) — PDMY03 critical-state constants as optional flags + dilation-brake note** ([[133_pdmy_notes]]) — `nDMaterial PressureDependMultiYield03 ... <-ei $e0> <-cs1 $v> <-cs2 $v> <-cs3 $v>` in Tcl and Python (defaults 0.6/0.9/0.02/0.7 = the former hard-coded values; byte-identical when omitted, proven against a baseline captured from the unmodified build). Fixes the cross-material leak in the `matCount%20` reallocation that the hard-coding had hidden. The note checks the TIMs reading of the PDMY "dilation brake": `isCriticalState()` is a CROSSING detector (returns 1 only for an increment that crosses the CSL), so retuning the CSL cannot produce a volumetric-rate plateau; with defaults a dense sand needs 10–18 % dilative volumetric strain to reach the line. | vanilla-material option + guide note | — (vanilla `ND_TAG_PressureDependMultiYield03` 112) | `SRC/material/nD/soil/PressureDependMultiYield03.{h,cpp}`, `SRC/material/nD/TclModelBuilderNDMaterialCommand.cpp`, `tests/test_wp133_pdmy03_cs_params.py`, `tests/wp133_pdmy03_deck.{py,tcl}`, `tests/wp133_pdmy03_byteid_baseline.json`, `tests/wp133_pdmy03_tcl_baseline.txt`, `Ladruno_implementation/133_pdmy_notes.md` | draft PR open | [#866](https://github.com/nmorabowen/OpenSees/pull/866) |
