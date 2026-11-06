---
wp: PR-35
title: "35, #52 -- 2 vanilla row(s)"
pr: "#35, #52"
files: ["`SRC/analysis/analysis/DirectIntegrationAnalysis.cpp`", "`SRC/analysis/analysis/StaticAnalysis.cpp`"]
table: "main"
legacy_seq: [122, 123]
---
| `SRC/analysis/analysis/DirectIntegrationAnalysis.cpp` | Profiler phase seams (P2, `OPS_PROFILE_SCOPE` step/newStep/solveCurrentStep/commit in `analyzeStep`) + `// Ladruno P0#3` per-step series hook (`OPS_PROFILE_STEP(getCommitTag, getCurrentTime, dT, getNumIterations)` after commit). +`#include <profiler/ProfilerMacros.h>`. Additive/behavior-preserving. | [#35](https://github.com/nmorabowen/OpenSees/pull/35), [#52](https://github.com/nmorabowen/OpenSees/pull/52) |
| `SRC/analysis/analysis/StaticAnalysis.cpp` | Profiler phase seams (P2, `OPS_PROFILE_SCOPE` step/newStep/solveCurrentStep/commit in `analyze`) + `// Ladruno P0#3` per-step series hook (`OPS_PROFILE_STEP`, dt=0 for static load-stepping) at the end of each loop iteration. +`#include <profiler/ProfilerMacros.h>`. Additive/behavior-preserving. | [#35](https://github.com/nmorabowen/OpenSees/pull/35), [#52](https://github.com/nmorabowen/OpenSees/pull/52) |
