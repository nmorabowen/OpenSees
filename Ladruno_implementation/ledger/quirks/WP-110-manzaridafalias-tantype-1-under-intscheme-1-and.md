---
wp: WP-110
title: "ManzariDafalias TanType 1 under IntScheme 1 (and 0) was a STALE matrix — ModifiedEuler never wrote mCep — FIXED (WP-110, F15c)"
legacy_seq: 466
---
### `ManzariDafalias` TanType 1 under IntScheme 1 (and 0) was a STALE matrix — `ModifiedEuler` never wrote `mCep` — FIXED (WP-110, F15c)
- **Bites:** you set `TanType 1` (continuum elastoplastic tangent) with the recommended `IntScheme 1` and get modified-Newton convergence, or a tangent that is plainly Ce at a plastic state. Nothing warns.
- **Why:** `ManzariDafalias::ModifiedEuler` computes `aCep1`/`aCep2` for its `aCep_Consistent` chain (TanType 2) but never assigned `aCep` itself. `mCep` therefore kept whatever the last writer left: Ce from the last elastic step (`elastic_integrator` / the elastic branch of `explicit_integrator`), or a `Stress_Correction` leftover on a step whose last substep needed a correction. Scheme 0 (`MaxEnergyInc`) inherits it through `nCep`. `RungeKutta4` (scheme 3) has the same hole and is NOT fixed (scheme 3 is already warned against at construction).
- **Workaround/status (WP-110, #847):** fixed — `ModifiedEuler` starts from `aCep = Ce` and, on its normal exit, writes `GetElastoPlasticTangent` at the end-of-increment state (the `BackwardEuler_CPPM` convention). Observable now via `eleResponse(ele, 'tangent')`.
