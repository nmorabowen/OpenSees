---
wp: PR-295
title: "CriticalTimeStep per-element kernel"
pr: "#295"
status: "shipped — behavior-preserving (18/18 regression)"
section: "table"
legacy_seq: 41
---
| **CriticalTimeStep per-element kernel** — extracted `lumpElementMass` / `elementLambdaMax` / `elementCriticalDt` (self-report-aware) so the dt_cr query AND the SMS pass share one D8-safe lump+eigensolve | Integrator util | — | `SRC/analysis/integrator/CriticalTimeStep.{h,cpp}` | shipped — behavior-preserving (18/18 regression) | [#295](https://github.com/nmorabowen/OpenSees/pull/295) |
