---
wp: WP-129
title: "ModifiedEuler's step factor has no upper cap: an err = 0 substep sets q = +inf and the next substep swallows the rest of the increment (WP-129)"
legacy_seq: 488
---
### `ModifiedEuler`'s step factor has no upper cap: an `err = 0` substep sets `q = +inf` and the next substep swallows the rest of the increment (WP-129)
- **Bites:** on acceptance `q = max(0.8√(TolE/err), 0.5)`, then `dT = min(q·dT, 1 − T)`. When both stages are elastic in stress (`err = 0`, see the stress-only row), `q = inf`, and whatever the stages did to α is extrapolated over the whole remaining increment in one step. There is also no "no growth right after a rejection" rule, so accept/reject can oscillate. WP-128 found it inert on the escape chains (they escape at `dT = 1` already), so it is hygiene, not a mechanism.
- **Workaround/status:** SAS-ME: `q = clamp(0.9√(TolR/err), 0.1, 1.1)` with `err ≥ DBL_EPSILON`, and `q ≤ 1` on the substep after a rejection.
