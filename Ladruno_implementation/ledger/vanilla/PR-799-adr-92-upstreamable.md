---
wp: PR-799
title: "ADR-92 [#799] -- upstreamable-table row(s)"
pr: "#799"
files: ["`SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}`"]
table: "upstreamable"
legacy_seq: [516]
---
| `SRC/{tcl/tclMain.cpp,interpreter/PythonModule.cpp}` | Splash-banner feature list regen (`FEATURES-START/END`) via `patch_banner.py` — the `LadrunoSANISAND` line (last touched by the `-maxSubsteps` cap, [#792](https://github.com/nmorabowen/OpenSees/pull/792)) now names `-implex` explicitly: `LadrunoSANISAND — -implex (IMPL-EX), p_r/p_min, -maxSubsteps cap` (64 cols, ADR-87 ≤70 rule). Source of truth `Ladruno_scripts/banner_features.txt`; never hand-edited, +1 line in each SRC file. IMPL-EX itself (`ManzariDafalias`/`LadrunoSANISAND`'s explicit tangent-operator return map) shipped in ADR-92 P1 ([#798](https://github.com/nmorabowen/OpenSees/pull/798), merge commit `8e8a69aee`) with no banner row of its own — this closes that gap. | ADR-92 [#799] |
