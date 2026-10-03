---
wp: WP-127
title: "ManzariDafalias::ModifiedEuler ACCEPTS a substep that FAILED its error test once dT reaches dT_min — silently, with a clamp to Mc (finding C, counted WP-127)"
legacy_seq: 481
---
### `ManzariDafalias::ModifiedEuler` ACCEPTS a substep that FAILED its error test once `dT` reaches `dT_min` — silently, with a clamp to `Mc` (finding C, counted WP-127)
- **Bites:** `ManzariDafalias.cpp`, the `curStepError > TolE` branch: when `dT == dT_min` (1e-6) the substep is **accepted anyway** — elastic tangent flag, `NextStress = nStress`, then if `eta > m_Mc` a RADIAL clamp of the deviator to `m_Mc` (the critical-state ratio, not `M^b`, and compression-side only: `eta` is taken from `|dev|/tr`, no Lode dependence), and `alpha` re-derived as `CurAlpha + 3(dev(s)/tr(s) - dev(s_n)/tr(s_n))`, i.e. NOT integrated. `T` advances, so the update reports success. Nothing counted it, nothing printed. Next door, `p < p_r` at `dT == dT_min` makes ModifiedEuler RETURN with `T < 1` — the increment is partially integrated and the caller is not told either.
- **Workaround/status:** WP-127 counts both per integration point: `substepStats` columns `forcedAtDTmin` / `forcedClampMc` / `abandonedLowP` (and `last*`), and the replay trace marks each such substep (codes 2/3/6). **Read them before trusting a ring-point stress.** The behaviour itself is unchanged (byte-identical); whether to refuse instead is F18 (WP-129/130).
