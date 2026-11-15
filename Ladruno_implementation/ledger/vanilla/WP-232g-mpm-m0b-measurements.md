---
wp: WP-232g
title: "WP-232g: M0b measurements -- profiler scopes in the analysis set-up"
pr: "#1089"
date: 2026-11-15
table: "main"
files: ["`SRC/analysis/analysis/StaticAnalysis.cpp`", "`SRC/analysis/analysis/DirectIntegrationAnalysis.cpp`"]
---
| `SRC/analysis/analysis/StaticAnalysis.cpp` | `// Ladruno WP-232g`: profiler scopes only — `domainChanged` around the call in `analyze()`, and `dc.clearAll`, `dc.handle`, `dc.numberDOF`, `dc.setSize` ⊃ `dc.s.graph`, `dc.integrator`, `dc.algorithm` inside `domainChanged()` (mirroring `DirectIntegrationAnalysis`), so the M0b harness can split a quasi-static mode-A MPM step's rebuild. No call or return changes; `OPS_PROFILE_SCOPE` is a no-op when the profiler is off (`ProfilerMacros.h:77`). | WP-232g [#1089](https://github.com/nmorabowen/OpenSees/pull/1089) |
| `SRC/analysis/analysis/DirectIntegrationAnalysis.cpp` | `// Ladruno WP-232g`: `dc.clearAll`, `dc.integrator` and `dc.algorithm` scopes in `domainChanged()`, beside the existing `dc.handle` / `dc.numberDOF` / `dc.setSize`, so a transient mode-A step's rebuild splits with no remainder (Newmark re-reads every node in its `domainChanged`, ADR 232 W1). Scopes only. | WP-232g [#1089](https://github.com/nmorabowen/OpenSees/pull/1089) |
