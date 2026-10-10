---
wp: WP-128
title: "ManzariDafalias::Stress_Correction gives up SILENTLY and hands back f > 0 as a successful update (WP-128)"
legacy_seq: 497
---
### `ManzariDafalias::Stress_Correction` gives up SILENTLY and hands back `f > 0` as a successful update (WP-128)
- **Bites:** when neither correction direction reduces `|f|`, the loop takes a bare `return;` ("Couldn't decrease the yield function", printed only under the compile-time `debugFlag`). The first direction is `λ` along `C:R` with `h b`; the second is `λ` along `∂f/∂σ`. The return leaves `NextStress/NextAlpha` uncorrected. ModifiedEuler advances `T`, and nothing downstream checks f. At low p and high η, where `∂f/∂σ` is dominated by `(n:r)/3·I` and the cone radius is `√(2/3)·m·p ≈ 1e-3 kPa`, this is routine. TIMs b8 1950/3 with `γ_xy = 1e-5` returns rc 0 with `f = +0.0139`, about 10× the cone. From the floor with `dε_yy = 3e-4` it returns `f = 11.2 kPa`. 4 of 640 ring replays return `f > 1e-6` with rc 0. It breaks the material checklist's "a non-converged return map must FAIL".
- **Workaround/status:** read `f` yourself: the replay returns `f after`. Fix direction for WP-129: the give-up and the `i == maxIter` "still outside" branches raise a refusal that ModifiedEuler turns into a failed update.
