---
wp: LEGACY
title: "nDMaterial DruckerPrager cannot have pressure-dependent moduli AND plasticity, and updateMaterialStage does not reach it"
legacy_seq: 257
---
### `nDMaterial DruckerPrager` cannot have pressure-dependent moduli AND plasticity, and `updateMaterialStage` does not reach it
- **Bites:** anyone reaching for `DruckerPrager` as a cheap perfectly plastic
  soil next to PDMY. `mElastFlag` is `0 = elastic+no param update, 1 =
  elastic+param update, 2 = elastoplastic` (default), and `updateElasticParam`
  only rescales `K, G` by `sqrt(1 + p/p_atm)` when `mElastFlag == 1` — which is
  an ELASTIC state. So the PDMY-style `G ~ sqrt(p)` and plasticity are mutually
  exclusive; a plastic DruckerPrager has constant moduli, full stop. Separately,
  `setParameter` returns -1 for `"updateMaterialStage"` (it answers to
  `"materialState"` instead), so the usual `ops.updateMaterialStage(...)` staged
  gravity idiom silently does nothing here.
- **Rule:** for a collapse load this does not matter (a limit load is
  independent of the elastic constants), but the SETTLEMENT at which it arrives
  is not — so any `q` quoted at a fixed `s/B` criterion is affected, and that
  must be stated. For the gravity stage, check whether the elastic K0 state is
  admissible (previous entry) and just solve it plastically in one step rather
  than hunting for a stage flip that is not wired up. Perfect plasticity needs
  `Kinf = Ko = delta1 = delta2 = H = theta = 0`; `mHprime = (1-theta)*H` so
  `H = 0` alone suffices. The conversion from a measured cone is
  `rho = sqrt(2)*alpha` and the sqrt(J2) cohesion intercept is `sigma_y/sqrt(3)`
  (yield is `||s|| + rho*I1 - sqrt(2/3)*sigma_y`); the tension cutoff `mTo` is
  placed exactly at the cone apex, so it adds nothing.
  *2026-07-30 (ADR-79 collapse study).*
