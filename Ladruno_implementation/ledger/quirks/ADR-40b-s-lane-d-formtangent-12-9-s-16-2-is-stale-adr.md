---
wp: ADR-40b
title: "ADR-40b's lane-D \"formTangent = 12.9 s (16.2%)\" is STALE — ADR-67 P-NEW-1's constant-mass tangent cache is on by default and removed it"
legacy_seq: 223
---
### ADR-40b's lane-D "formTangent = 12.9 s (16.2%)" is STALE — ADR-67 P-NEW-1's constant-mass tangent cache is on by default and removed it
- **Bites:** re-running lane D expecting the ADR-40b explicit phase mix, seeing `formTangent ≈ 0.00%`, and hunting for a broken model or a lost scope. Nothing is broken: `CentralDifferenceLadruno` now ships `massCache = true` by default (`:89`, `:154`, `formTangent` override at `:259`) — the ADR-67 P-NEW-1 constant-mass tangent cache, i.e. "`-factorOnce` behaviour with safe invalidation", which is exactly the fix ADR-40b's Finding-3 item 1 recommended. Measured 2026-07-25: lane D `formTangent` **1.3 ms of a 29.5 s step** (was 12.9 s of 79.9 s), the rest of the mix shifting to `formUnbalance` 34.0% / `newStep` 32.0% / `update` 29.9%.
- **Workaround/status:** general rule for this fork — **a phase baseline in a dated report may have been optimized away by a later ADR.** Re-measure before quoting; check `LEDGER_implementations` for a shipped fix on that path first. *2026-07-25 (ADR-75b L3-0).*
