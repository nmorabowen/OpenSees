---
wp: LEGACY
title: "timeSeries Path returns 0 BEYOND its last time node — float-accumulated pseudo-time overshoots and collapses prescribed strains"
legacy_seq: 76
---
### `timeSeries Path` returns 0 BEYOND its last time node — float-accumulated pseudo-time overshoots and collapses prescribed strains
- **Bites:** a multi-stage prescribed-strain cyclic driver built from `timeSeries('Path', ...)` + `LoadControl(1/nper)` over N stages. The intended end pseudo-time is `N`, but `nper×N` accumulations of `1/nper` (e.g. `1/80`, not exact in binary) land at `N + epsilon`. `Path` returns **0** outside `[t_first, t_last]`, so at the FINAL step every `sp`-prescribed DOF drops to 0 → the element snaps to ~zero strain. For the RC cyclic interlock this looked like a phantom `-3.18` crack-shear spike (crack "closed" at `en≈0` ⇒ `v_ci,max` jumps to its max `0.18√fc/0.31`) — a TEST artifact, not a kernel bug.
- **Why:** `PathSeries::getFactor(t)` returns 0 for `t > t_last` (and `< t_first`). The classic mis-diagnosis is "the cap is wrong"; the real tell is reading the GP strain at the offending step (it's ~0, not the held value).
- **Fix (robust):** pad the Path with an extra HOLD node beyond the analysis end — `times=[0..N, N+1]`, `values=[...,last, last]` — so any overshoot interpolates between two equal final values. (`-useLast` as a trailing openseespy arg did NOT take effect in this build; the pad is reliable.) Learned 2026-06-16 building the Phase-2b cyclic shear driver `_path_stages`.
