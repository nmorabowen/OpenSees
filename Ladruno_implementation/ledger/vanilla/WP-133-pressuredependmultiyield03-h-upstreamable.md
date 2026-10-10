---
wp: WP-133
title: "866 -- upstreamable-table row(s)"
pr: "#866"
files: ["`SRC/material/nD/soil/PressureDependMultiYield03.{h,cpp}`", "`SRC/material/nD/TclModelBuilderNDMaterialCommand.cpp`"]
table: "upstreamable"
legacy_seq: [663, 664]
---
| `SRC/material/nD/soil/PressureDependMultiYield03.{h,cpp}` | `// Ladruno WP-133` (TIMs F23a) — the critical-state constants `ei, cs1, cs2, cs3` (hard-coded `0.6, 0.9, 0.02, 0.7` in the constructor body) become four trailing DEFAULTED constructor arguments with the same defaults (header), and `OPS_PressureDependMultiYield03` takes them as flags after all positional args (`-ei -cs1 -cs2 -cs3`; token peek, then the flag tail is parsed after the positional loops; unknown option / missing value refuses). Also fixes the `matCount%20` reallocation, which wrote the NEW material's constants into every existing slot (`einitx[i] = ei`, etc.) and leaked the old four arrays: it now copies each slot's own values and frees them. Strictly additive otherwise; the old locals are kept as comments. Byte-identical when the flags are omitted (baseline from the unmodified build, `tests/test_wp133_pdmy03_cs_params.py` G1). See [[133_pdmy_notes]]. | [#866](https://github.com/nmorabowen/OpenSees/pull/866) |
| `SRC/material/nD/TclModelBuilderNDMaterialCommand.cpp` | `// Ladruno WP-133` — the Tcl `PressureDependMultiYield03` ladder branch accepts the same `-ei -cs1 -cs2 -cs3` flags (first recognised flag cuts `argc`; pairs parsed with `Tcl_GetDouble`; bad option → `TCL_ERROR`) and passes them to the constructor. No other branch touched. | [#866](https://github.com/nmorabowen/OpenSees/pull/866) |
