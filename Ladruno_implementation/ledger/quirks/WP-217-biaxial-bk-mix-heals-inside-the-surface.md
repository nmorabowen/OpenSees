---
wp: WP-217
title: "LadrunoCohesiveHingeBiaxial with -bk eta != 1: D followed the current mix while unloading (non-conservative, healing, wrong tangent)"
date: 2026-10-29
---
### `LadrunoCohesiveHingeBiaxial -bk` unloading tangent fails an FD check, or a closed loop inside the surface does work

- **Bites:** only `eta != 1`. frame-ladruno-validation 00b (B3): at (1.5, 3) after a
  square path, the diagonal unloading tangent missed the central FD of M(alpha) by 3.5e-4
  (4.3e-2 at penalty ratio 10). A closed loop strictly inside r < r_max did +1.15e-4 of
  sum|M·dalpha|, and the reversed loop created that energy. D also decreased (healed) on
  non-radial unloading.
- **Why:** the law evaluated `D = D(r_max, w(alpha))` at the CURRENT mode mix, also while
  `r < r_max`. With eta = 1, D does not depend on the mix (A = S/Esoft = 2/(ratio-1)), so
  this was harmless. With a B-K mix, the secant `(1-D)Kpen` changed with the jump direction
  inside the elastic domain. M(alpha) was then not a gradient field, the frozen diagonal
  dropped `dD/dw·dw/dalpha`, and D could fall.
- **Fix (WP-217, ADR 34 2026-10-29 amendment):** D is the internal variable. `(r_max, w_z, w_y)`
  is the state at D's last growth, and `D_c = D(r_max, w)` at THAT mix. Loading is
  `D_law(alpha) >= D_c`; otherwise D stays frozen. The stress is continuous (a max of
  continuous functions), conservative and never heals, and the diagonal tangent is exact.
  Loading (and so every radial path, and the B-K dissipation) is unchanged. IMPL-EX
  extrapolates r_max at the frozen mix, so its step is linear in alpha.
  eta = 1 keeps the old `r >= r_max` branch literally. The two are the same law there.
- **Freeze-mix alone is not enough:** freezing the mix while keeping `r >= r_max` as the
  loading test makes the stress JUMP on re-entry at a different mix. Adding the missing
  tangent term instead fixes the FD check but keeps the non-conservative, healing law.
- **Meaning of `rMax` changed:** for eta != 1 the `rMax` response is r at the last damage
  growth, which can be below the largest r seen if the mix moved to a tougher mode.
- **Gate:** `tests/test_ladruno_cohesive_hinge_implex_bk.py` (FD of the unloading
  tangent, closed-loop work, monotone D, eta = 1 against the ADR 34 closed form).
