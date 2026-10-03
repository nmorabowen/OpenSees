---
wp: LEGACY
title: "‖D_current‖/‖D_initial‖ is NOT a usable material-degradation proxy — the Frobenius norm is bulk-dominated (so deviatoric plasticity barely moves it) and sign-b…"
legacy_seq: 371
---
### `‖D_current‖/‖D_initial‖` is NOT a usable material-degradation proxy — the Frobenius norm is bulk-dominated (so deviatoric plasticity barely moves it) and sign-blind (so a softening tangent makes it GROW)
- **Bites:** the natural-looking generalization of a damage-scaled stabilization
  (`LadrunoBrick`/`LadrunoQuad` Tier-A `Kstab`) from "scale on the damage scalar"
  to "scale on the tangent norm, which covers plasticity models too". It sounds
  strictly more general. It is strictly worse, in both directions at once.
- **Why (measured on the consistent J2 tangent at a pure-shear flow state,
  E=1, flow direction n):** `‖C_ep‖/‖C_e‖` = **0.943** at ν=0.2 with *zero*
  hardening, 0.967 at ν=0.3, 0.998 at ν=0.45, **1.000** at ν=0.499 — while the
  shear entry `D(3,3)` that actually carries the hourglass modes has gone to
  **0.000**. `‖C‖²= 9K² + 20G²` elastic vs `9K² + 16G²` at full plastic flow:
  the bulk block dominates and plasticity is deviatoric, so the norm is nearly
  blind to it, and blindest exactly in the near-incompressible regime where
  plastic flow lives.
- **Why (the other direction):** a norm has no sign. For a softening tangent with
  flow-direction slope `−(1+h)·2G`, the ratio goes 0.957 / 1.000 / 1.155 / 1.915 /
  6.733 at h = 0.5 / 1 / 2 / 5 / 20 — it **exceeds 1** past h≈1, gets clipped, and
  yields no degradation at all *exactly at localization*. A secant-returning
  material hides this; a true consistent damaged tangent (e.g. `LadrunoConcrete3D`
  P3b, with the `−σ⊗dω` rank update) does not.
- **Also:** a tangent-based scale is non-monotone in load history — it snaps back
  to 1 on elastic unloading, so any floored stabilization would toggle every load
  reversal, with `∂s/∂u` missing from the tangent. Damage ratchets; tangents do not.
- **Rule:** if you need a degradation proxy for a *mode-specific* stabilization,
  use the mode-specific entry (shear → `D(3,3)`, which `formUri` already does),
  monotonize it over history, and keep the floor. And degrade where the material
  **softens**, not merely where it yields — a hardening element has no localization
  to enable and still needs its hourglass modes controlled.
- **Status (2026-07-30):** challenge investigated and closed with NO code change;
  full study in [[11_brick_asdconcrete_integration]] §3.1 (incl. the measurement
  that the frozen-`Kstab` bias converges away at ≈O(h): `Ehg/W_ext` = 29 / 15.5 /
  5.7% and load bias 1.77 / 1.45 / 1.15 at 2 / 4 / 8 elements through the bending
  depth). Regression test `tests/test_ladrunoBrick_kstab_plasticity.py`.
- *2026-07-30 (Tier-A `Kstab` scope challenge).*
