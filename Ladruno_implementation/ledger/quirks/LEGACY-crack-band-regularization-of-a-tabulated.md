---
wp: LEGACY
title: "Crack-band regularization of a TABULATED backbone must re-apply adjust() after scaling"
legacy_seq: 85
---
### Crack-band regularization of a TABULATED backbone must re-apply adjust() after scaling
- **Bites:** cloning ASDConcrete3D's `HardeningLaw::regularize` (RC stack Phase 3b, `-autoRegularization`). The post-peak strain is rescaled to hit the target fracture energy `g_reg=G_f0·(lch_ref/lch)`. ASDConcrete3D calls `adjust()` INSIDE the scaling loop and once after — easy to omit because the **fracture energy comes out correct without it** (g uses only x,y). But the stress-update reads `q`/damage, and stretching x while preserving the plastic-to-inelastic ratio can drive the tail plastic strain `x−q/E` BACKWARD (non-physical) → wrong `q`/`dt_bar`/`dc_bar` (~6% on steep-mid-softening-damage backbones). Energy-objectivity gates (and the numpy/g++ energy gates) CANNOT catch it.
- **Fix (proven):** mirror the reference — call `adjustBackbone()` (E-cap + monotone plastic strain + non-decreasing damage + `q=y/(1-d)`) after the q/d update each iteration AND once after the loop. Gate it with a steep-damage backbone asserting the post-peak plastic strain is monotone non-decreasing (the energy gate alone is blind). Also: `adjustBackbone` derives `d` from the stored `q`, so the FIRST point must have `q[0]=0` (⇒ d0=0), matching what `buildBackbone` always produces — a hand-built oracle backbone with `q[0]=1e-12` makes adjust read `d0=1` and corrupts everything. Learned 2026-06-18 (3-agent review of [[19_ladruno_rc_shell_adr|LadrunoRCConcrete]] Phase 3b).
