---
wp: WP-138
title: "A SANISAND BVP that walls on loadingNonPosH refusals is at DM04's h-singularity — TolR, TanType, -maxSubsteps and the Krylov rung do not lift it (WP-138)"
legacy_seq: 511
---
### A SANISAND BVP that walls on `loadingNonPosH` refusals is at DM04's h-singularity — TolR, TanType, -maxSubsteps and the Krylov rung do not lift it (WP-138)
- **Bites:** anyone tuning the integrator to push a SANISAND footing past its step-floor wall. On the WP-138 strip footing (build 7936ed6e0), SAS-ME (IntScheme 129, TanType 0, TolR 1e-4) floors at s/B 0.0508 on `loadingNonPosH` (first at 0.0363). Each knob made it worse:
  - TolR 1e-3 floors EARLIER, at 0.0410;
  - TanType 1 + `-maxSubsteps 20000` collapses ds to ~1e-6 m and floors at 0.0114, at 13× the cost per unit s/B;
  - ModifiedEuler floors at 0.0292.
  The refusing points are pre-peak (ρ_α < 1). The set is {(α − α_in):n = 0, b:n ≤ 0}: h = b0/a hits its 1e10 cap and K_p → −∞ (WP-150 memo, #892). The WP-134 oracle hits the same 0/0.
- **Workaround/status:** none on the integrator side. No material switch removes it (final ladders, `138_footing_sas_me_ab.md` §11): dilatancy off (A0 = 0.001) only delays the onset (0.0363 → 0.0426) and walls at 0.0499; Presidual 0.5–20 kPa and e_init 0.65–0.85 do not remove it. (An interim snapshot had suggested A0 = 0.001 cleared it -- withdrawn.) The wall states need the concave extension meridian, c = 0.71 < 7/9: DM04 at c = 0.80 takes them from 102/320 to 0/320 failing trials (WP-151 memo §2.5). Two routes, TIMs' call: R1 (WP-151, #893: `-sasHFloor 1 -sasReseatHyst 1 [-sasSoftCap 0.5]`, two coupled opt-in flags, an h floor everywhere plus a hysteretic α_in re-seat), or a calibration with c ≥ 0.78.
  - The R1 oracle: 0/320 failures with both flags; each flag alone fails; a floor only where b:n ≤ 0 fails 102/320.
  - The wall is a Zeno accumulation of re-seats on b:n → 0⁺. NonPosH is only its b:n < 0 exit.
  - The owner approved R1 as opt-in; TIMs decide on its use. *2026-09-28.*
