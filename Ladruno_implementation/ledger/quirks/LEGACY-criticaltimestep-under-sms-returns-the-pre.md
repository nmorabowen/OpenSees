---
wp: LEGACY
title: "criticalTimeStep() under SMS returns the PRE-scaling element pencil — the post-scaling effective limit is report-only with NO getter"
legacy_seq: 6
---
### `criticalTimeStep()` under SMS returns the PRE-scaling element pencil — the post-scaling effective limit is report-only with NO getter
- **Bites:** any driver/script that queries `criticalTimeStep()` on a mass-scaled run (CentralDifferenceSMS/SMSConsistent, ExplicitBathe `-sms`) expecting the stable-Δt after scaling. It gets the un-augmented sliver limit instead — `dt = safety*criticalTimeStep()` collapses Δt to the pre-scaling value and silently defeats SMS (conservative, so it "works", just wastes the entire scaling benefit).
- **Why:** `CentralDifferenceLadruno::getCriticalTimeStep()` (`CentralDifferenceLadruno.cpp:742–746`) and `ExplicitBathe::getCriticalTimeStep()` (`ExplicitBathe.cpp:1435–1441`) return the raw element-pencil value, which cannot see nodal injected mass (scaling writes `Node::setMass`; the pencil reads `ele->getMass()`). The post-scaling limit exists (`setSMSEffectiveLimit` ← `minDtSelfReport`, P3 #475) but is consumed ONLY by the `newStep()` "[PRE-SCALING estimate]" report — protected setter, no getter, no command plumbing.
- **Workaround/status (2026-07-02, ADR-65 adversarial gate):** treat the query as a conservative lower bound, or don't query under SMS at all (SMS runs are sized to `dtTarget` at construction — just use `dtTarget`). An SMS-aware adaptive driver needs a new `getSMSEffectiveLimit()` accessor first. See [[65_ladruno_explicit_dt_strategies_adr]] §Route B.
