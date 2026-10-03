---
wp: LEGACY
title: "Tension stiffening is monotonic-scope (live-ε1 floor re-inflates on unload); TS+interlock mixes live-p1 with frozen-crack"
legacy_seq: 84
---
### Tension stiffening is monotonic-scope (live-ε1 floor re-inflates on unload); TS+interlock mixes live-p1 with frozen-crack
- **Bites:** `-tensStiff` floors the stress to `σ_ts(ε1)` using the LIVE membrane principal strain `ε1` (no `ε1max` memory). Because `σ_ts` *decreases* with `ε1`, on UNLOADING the floor **re-inflates** (tracks `σ_ts(live ε1)` back UP) — a load→unload analysis shows the floored stress rising as strain drops. This is correct on a monotone loading branch but is NOT a hysteretic cyclic-tension model. Separately, TS uses the LIVE principal axis `p1` while the fixed-crack `-interlock` uses the FROZEN crack normal; once principal axes rotate, the TS normal-stress injection leaks a small shear `Δ·sin θ_rel cos θ_rel` onto the frozen crack plane that the interlock then bounds.
- **Fix / scope:** use `-tensStiff` for MONOTONIC / pushover analyses; combined TS+interlock is validated for PROPORTIONAL (non-rotating) loading only. The cyclic upgrade (deferred) is an `ε1max`-envelope floor + secant unload + flooring along the frozen crack plane when cracked. Documented in the kernel comment + [[LadrunoRCConcrete_guide]] §4.7. Learned 2026-06-18, [[19_ladruno_rc_shell_adr|LadrunoRCConcrete]] Phase 3a.
