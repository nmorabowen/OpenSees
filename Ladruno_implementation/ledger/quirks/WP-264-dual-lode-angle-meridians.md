---
wp: WP-264
title: "A dual pushed through the Lode angle (acos of cos 3 theta) is NaN on every uniaxial and equal-biaxial state"
pr: "#1175"
date: 2026-11-28
status: "open"
---
### Dual-number gradient of a Menetrey-Willam / Willam-Warnke function is NaN or inf at uniaxial, equal-biaxial and confined states (repeated principal stresses)

- **Bites:** anyone templating a function of the Lode angle on `ladruno_dual::Dual<N>`:
  `invariants()` → `theta = acos(c3)/3` → `lodeR(theta, e)` (Concrete3D's surface, potential and damage drivers,
  LadrunoRCConcrete, the HS / SANISAND Lode factors). It is also the eigen-solver degeneracy hazard in other clothes:
  a repeated principal pair puts the state on a meridian, c3 = cos 3θ = ±1.
- **Why:** dθ/dc3 = −1/(3 sin 3θ) is infinite on the meridians, and dr/dθ = 0 there. The dual computes 0 · inf =
  NaN. Near the meridian it computes a large times a small number, losing the digits that `sqrt(1 − c3²)` cancels.
  Measured (`tests/_testbed/concrete3d_dual_check.cpp`, mutant M4): NaN on 56 of 80 exactly repeated-pair states.
  A second trap sits beside it: clamping the whole Dual at |c3| > 1 drops a real derivative. c3 is smooth and has a
  zero gradient at ±1, but round-off pushes it past ±1 near a meridian. Measured: up to 5e-6 of the ε̃ gradient
  against an eigen-free `Dual<6>` pass.
- **Fix (WP-264, `LadrunoConcrete3DKernel.h`):** do not differentiate θ. Keep θ a double, `acos(value(c3))/3`.
  Carry the derivative through c3 with the closed-form dr/dc3 (`lodeRT` → `drDcos3thetaRegularized`, finite on both
  meridians). Clamp only the VALUE of c3 (`clampValueKeepDerivative`). Do not push a dual through `eig3sym` either.
  The drivers are symmetric in the principal values, so differentiate in w (`Dual<3>`) and lift with Σ f_a v_a⊗v_a,
  which is basis-independent inside a repeated pair (gate G3: 1.1e-12).
- **Related:** the cone vertex of ε̃ (`rad == 0`) is a third 0 · inf (`sqrt`'). It is reached off the hydrostatic
  axis by ε̃_c, whose argument wc vanishes whenever every tensile principal is cracked and the rest are zero, as in
  plane-stress tension. Guarded in `equivStrainGeneralT`: the linear part's derivative, a finite subgradient. There
  α_c = 0, so `dkdc = ac·detc + …` never reads it, but a NaN times 0 would still have poisoned the tangent.
