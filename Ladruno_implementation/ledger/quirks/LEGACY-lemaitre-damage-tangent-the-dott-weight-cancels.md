---
wp: LEGACY
title: "Lemaitre damage tangent: the dotT weight cancels the tensor↔engineering shear factor"
legacy_seq: 58
---
## Lemaitre damage tangent: the `dotT` weight cancels the tensor↔engineering shear factor

When assembling `∂D/∂ε` for the Lemaitre coupled-damage consistent tangent
(`LadrunoJ2Kernel.h::returnMapDamaged`), the `∂(p̄)/∂ε` term goes through
`∂‖M‖/∂M = wt·M/‖M‖` where `wt = (1,1,1,2,2,2)` are the `dotT` factor-2 weights on
the off-diagonal (shear) pairs. The strain variable being differentiated is then
converted tensor→engineering by `w = (1,1,1,½,½,½)`. **`wt_j · w_j = 1` for every
component** (normal: 1·1; shear: 2·½), so the net engineering derivative is just
`(2G/h)·n_j` with **NO shear half-factor**. The natural-but-wrong instinct is to
carry the `w_j=½` on shear (mirroring the strain-mapping elsewhere in the file);
that silently halves the shear columns of `∂p̄/∂ε` and breaks the tangent ONLY on
shear-coupled (3D-mixed) states — uniaxial paths pass clean, so a uniaxial-only FD
check misses it. Caught by the 3D-mixed FD case in `tests/ladruno_damage_check.cpp`
(error was a fixed 1.5e-5 independent of the FD step ⇒ a real bug, not truncation).

- **Rule:** for a `dotT`-norm gradient differentiated w.r.t. engineering strain, the
  `dotT` weight and the tensor→engineering factor cancel — use the bare component.
  Always FD-check the consistent tangent on a **shear-inclusive** state, not just
  uniaxial. Learned 2026-06-02 (Lemaitre damage, [[15_lemaitre_ductile_damage_adr]]).
