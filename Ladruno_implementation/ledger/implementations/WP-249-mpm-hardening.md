---
wp: WP-249
title: "WP-249: MPM hardening before phase 7 -- cut-retry regression test (#1095), failDump switch, restart pack v1, small review follow-ups"
pr: "#1143"
date: 2026-11-23
status: "shipped"
section: "table"
---
| **WP-249 — MPM hardening before phase 7** (program #1075; not a §9 phase of [[232_ladruno_mpm_general_driver_adr]]; record `249_mpm_hardening.md`). #1095 (not reproduced, stays open): a 3D DruckerPrager cut-retry guard test across connectivity changes (`tests/test_ladruno_mpm_wp249_hardening.py`, child process, bit for bit against the direct sub-steps). #1106 item 2: the process-wide `mpm debug failDump` switch is cleared at mode-A teardown (`mpm reset`, `remove loadPattern`, `wipe`) and re-asserted from the live host's setting; `mpm info failDump`. #1106 item 4: `addFbarTangent`'s scratch is `thread_local`. #1106 item 5: `detN`'s 3×3 branch is `ladruno_hypo::det3` (bit-identical). #1103 item 2: the mode-A restart reader accepts pack version 1 (ladruno-v20261116) with no body force; `mpm debug packVersion 1` writes it, and falls back to the current layout for a body force or a B-spline basis. #1110: checked and narrowed (`Ladruno_scripts/mpm_249/`). | LoadPattern host + element fixes | **none — no new classes** | `SRC/domain/pattern/ladrunoMPM/{LadrunoMPMCore.h,LadrunoMPMImplicit.cpp,OPS_LadrunoMPM.cpp}`, `SRC/element/ladrunoMPMPoint/LadrunoMPMPoint.cpp`, `tests/test_ladruno_mpm_wp249_hardening.py`, `Ladruno_scripts/mpm_249/{s2_traction_probe.py,s2_traction_probe_results.json}` | shipped | [#1143](https://github.com/nmorabowen/OpenSees/pull/1143) |
