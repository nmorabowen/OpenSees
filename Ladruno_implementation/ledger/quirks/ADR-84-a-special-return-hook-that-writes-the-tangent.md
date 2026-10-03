---
wp: ADR-84
title: "A special_return hook that writes the tangent itself SILENTLY OVERRIDES tangent_type (FIXED, ADR-84 P3)"
legacy_seq: 303
---
### A `special_return` hook that writes the tangent itself SILENTLY OVERRIDES `tangent_type` (FIXED, ADR-84 P3)

- **Bites:** any ASDP yield function opting into `yf_has_special_return`
  (today only `MohrCoulombTensionCutoff_YF`). `Backward_Euler` normally ends at
  `ComputeTangentStiffness()`, which is the ONLY place the material's configured
  `INT_OPT_tangent_operator_type` (`Elastic` / `Continuum` / `Secant` /
  `Numerical_Algorithmic_*`) is consulted. The `special_return` hook returns
  BEFORE that call and assigned `Stiffness` directly, so every Gauss point the
  hook resolved got whatever the YF happened to write — regardless of the deck.
  ADR-84 P0 wrote a Secant blend `(E+D)/2` there, so `tangent_type Continuum`
  was a **no-op on exactly the Gauss points that needed it**.
- **The symptom is diagnostic and easy to misread:** changing `tangent_type`
  changes NOTHING. The Cerro Lindo M5 report tried shipped defaults,
  `strict_convergence` with a scaled tolerance, and `tangent_type = Continuum`,
  and all three stalled at the same place (λ = 0.50 / 0.45 / 0.475) — which
  reads as "the material is broken in a way no setting can reach", when in fact
  one of those settings was never applied. **If a knob provably does nothing,
  suspect an early `return` upstream of where the knob is read.**
- **Why it mattered here (not just tidiness):** the blend is not a conservative
  choice at a multi-surface corner. Measured against numerical differentiation
  of the hook's own return map, the true tangent there is a rigid attractor
  (`min|eig| = 9.4e-4`, `cond = 6.8e9` — the stress does not move for a 20x
  range of strain increment), the raw Koiter tangent reproduces it to 12%, and
  the shipped blend is **87% wrong with `min|eig| = 2.0e6` and `cond = 3.2`**:
  it is *well-conditioned where the truth is rank-deficient*. A global Newton
  told "push here and stress rises by 2e6*deps", whose stress then does not move
  at all, does not fail — it **stalls**, at a vanishing displacement increment
  against a residual that never falls. Plain MC on the same path hands out
  `cond ~1e16` and converges. A wrong-but-invertible tangent is worse than a
  singular one, because it never announces itself.
- **Workaround/status (2026-08-13, ADR-84 P3, PR):** the YF now returns the RAW
  active-set tangent (`SPECIAL_RETURN`'s documented contract) and the integrator
  applies the configured operator at the call site, as the generic path does.
  Default `Secant` keeps it byte-identical; `Continuum` finally delivers the raw
  tangent (3.5x fewer Newton iterations on `test_adr84_p3_confined_corner.py`).
  **Rule for any future hook that returns early from `Backward_Euler`: return
  the raw operator and let the integrator apply policy — never bake a tangent
  choice into a yield function.**
