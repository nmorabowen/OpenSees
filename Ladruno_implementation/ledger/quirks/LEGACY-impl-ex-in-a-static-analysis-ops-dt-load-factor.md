---
wp: LEGACY
title: "IMPL-EX in a STATIC analysis: ops_Dt (load-factor pseudo-time) is erratic — guard the extrapolation time-factor"
legacy_seq: 78
---
### IMPL-EX in a STATIC analysis: `ops_Dt` (load-factor pseudo-time) is erratic — guard the extrapolation time-factor
- **Bites:** porting the ASDConcrete3D IMPL-EX recipe (`tf = dtime_n / dtime_n_commit * alpha`, `dtime_n = ops_Dt`) into a material and running it under a STATIC `DisplacementControl`/`LoadControl` analysis (e.g. a quasi-static cyclic wall). The IMPL-EX extrapolation `x_ext = x_n + tf·(x_n − x_{n-1})` detonates — damage jumps to garbage, the first step diverges immediately — even though the same material is fine in dynamics.
- **Why:** in a STATIC analysis the domain "time" IS the load factor λ, not a physical Δt. With `DisplacementControl`, λ's increment is whatever satisfies the controlled DOF (can be ~1e5 for a stiff structure), and `loadConst('-time',0.0)` (the standard gravity-then-pushover idiom) RESETS λ to 0 mid-run. So `ops_Dt = λ_n − λ_{n-1}` is huge / tiny / negative across steps, and `tf = ops_Dt/ops_Dt_commit` becomes a wild multiplier on the threshold extrapolation. ASDConcrete3D ships no clamp because it is typically exercised in dynamics (or with `-dtime`/user-defined Δt) where Δt is smooth.
- **Fix (proven):** clamp the time-factor in `implexTimeFactor()`: fall back to `alpha` whenever `!commitDone` or `dtime_n<=0` or `dtime_n_commit<=0` or the ratio is non-finite; clamp the result to `[_, 2·alpha]` so a single pseudo-time spike cannot blow up the extrapolation. For static IMPL-EX this effectively makes `tf≈alpha` (uniform extrapolation) — correct, since static load steps are meant to be uniform. Learned 2026-06-17 pulling Phase-4 IMPL-EX forward for `LadrunoRCConcrete`.
