---
wp: LEGACY
title: "equationConstraint refuses zero coefficients (\"WARNING invalid rcoef inputs\") — skip the zero-arm terms when emitting plane-section rows"
legacy_seq: 157
---
### `equationConstraint` refuses zero coefficients ("WARNING invalid rcoef inputs") — skip the zero-arm terms when emitting plane-section rows

- **Bites:** emitting `u_i − u0 − θ·z_i = 0` rows for a node ON the reference
  plane (z_i = 0) fails parse: the EQ parser hard-rejects `coef == 0.0`.
- **Workaround (2026-07-07, ADR-66 G7):** drop the θ term when `|z_i| < tol`
  (the row degenerates to `u_i = u0`, which is exactly right) — the same filter
  the LadrunoTie shell-solid tests use for near-zero shape weights.
