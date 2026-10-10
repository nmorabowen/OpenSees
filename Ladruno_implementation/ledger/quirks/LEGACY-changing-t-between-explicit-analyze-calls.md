---
wp: LEGACY
title: "Changing Δt between explicit analyze() calls without revertToLastStep() mis-centers the leap-frog by (Δt_new−Δt_old)/2·a per change"
legacy_seq: 7
---
### Changing Δt between explicit `analyze()` calls without `revertToLastStep()` mis-centers the leap-frog by `(Δt_new−Δt_old)/2·a` per change
- **Bites:** any variable-Δt driver loop over `CentralDifferenceLadruno` (and the SMS subclasses). `newStep(dt_new)` after a step committed at `dt_old` advances `Vhalf += Δt_new·a_n` (`CentralDifferenceLadruno.cpp:574`) with no previous-Δt memory — but `v_{n−1/2}` is staggered for `Δt_old`, so the correct advance is `((Δt_old+Δt_new)/2)·a_n`. Each un-reseeded Δt change injects a systematic velocity error; small per change, compounding during ramps. No warning, nothing fails.
- **Why:** the leap-frog kernel is uniform-Δt by design; only the failure path was ever exercised with Δt changes, and that path happens to be correct because `revertToLastStep()` re-arms `firstStep` and rebuilds `v_{−1/2}` from committed state at the new Δt (`:420–436`, valid standalone per the `:413–419` comment).
- **Workaround/status (2026-07-02, ADR-65 adversarial gate):** call `revertToLastStep()` after the last commit before ANY `newStep` at a different Δt — grow or shrink, not just on failure. See [[65_ladruno_explicit_dt_strategies_adr]] §Route B (stagger reseed).
