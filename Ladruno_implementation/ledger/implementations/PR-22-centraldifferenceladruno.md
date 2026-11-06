---
wp: PR-22
title: "CentralDifferenceLadruno"
pr: "#22, #28, #468"
status: "shipped — full Linux link OK (#28); Zone-A CDL-1..10 battery; review-P1 #468: r…"
section: "table"
legacy_seq: 53
---
| **CentralDifferenceLadruno** — explicit leap-frog central difference done right (correct first-step starter, built-in `dt_cr`, clean full-step velocity, βK guard); coupled mode dropped → use `NewmarkExplicit 0.5` | Integrator | 33003 | `SRC/analysis/integrator/CentralDifferenceLadruno.{cpp,h}`, `tests/test_centralDifferenceLadruno_integrator.py` | shipped — full Linux link OK (#28); Zone-A CDL-1..10 battery; **review-P1 #468**: `revertToLastStep()` (re-seed from committed + re-arm starter ⇒ dt-retry == fresh restart; `tests/test_integrator_revert_to_last_step.py`) + NaN-capable breaker (`vectorIsFinite` — `pNorm(0)` is NaN-blind; `tests/test_explicit_nan_breaker.py`) | [#22](https://github.com/nmorabowen/OpenSees/pull/22), [#28](https://github.com/nmorabowen/OpenSees/pull/28), [#468](https://github.com/nmorabowen/OpenSees/pull/468) |
