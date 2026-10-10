---
wp: WP-134
title: "explicit_integrator's unload test and ModifiedEuler's loading index use n, not the yield-function gradient (U10, WP-134 → WP-129)"
legacy_seq: 494
---
### `explicit_integrator`'s unload test and ModifiedEuler's loading index use `n`, not the yield-function gradient (U10, WP-134 → WP-129)
- **Bites:** the start-on-surface branch decides plastic vs unload-then-plastic with `n : Δσ_trial / |Δσ_trial| > −√TolF`, and each stage's numerator is `2G n:de_dev − K dε_v (n:r)`. The true gradient is `∂f/∂σ = n − ⅓(n:α + √(2/3) m) I`, so for a volumetric trial (`n:Δσ ≈ 0` on isotropic compression) the test decides on round-off / the wrong term; WP-134 saw isotropic-compression probes enter the err = 0 path through it.
- **Workaround/status:** SAS-ME classifies with `∂f/∂σ : Δσ_trial` (predictor) and `∂f/∂σ : C : dε` (stages), the exact gradient. ModifiedEuler unchanged.
