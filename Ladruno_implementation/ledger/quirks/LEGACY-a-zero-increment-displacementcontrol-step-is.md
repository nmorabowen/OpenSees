---
wp: LEGACY
title: "A zero-increment DisplacementControl step is DEGENERATE — the load-factor correction dLambda = −dUabar/dUahat is unbounded, so \"hold everything\" must be LoadCo…"
legacy_seq: 441
---
### A zero-increment `DisplacementControl` step is DEGENERATE — the load-factor correction `dLambda = −dUabar/dUahat` is unbounded, so "hold everything" must be `LoadControl 0.0`
- **Bites:** the natural way to hold a state while something else converges (an augmentation
  sweep, a staged activation, a settling pass) is "same integrator, zero increment". Under
  `DisplacementControl` that is a trap: `newStep` sets `dlambda = theIncrement/dUahat = 0`, but
  the corrector `update()` still computes `dLambda = -dUabar/dUahat`
  (`SRC/analysis/integrator/DisplacementControl.cpp:332`) from whatever residual displacement
  is left, and with no increment to normalise it the load factor runs away. Measured on a
  2×2×2 elastic gate: a single zero-increment step moved the load factor from ~1.0 to
  **377.19**; a reviewer's fixture reached **−1.6e38**. The step often still returns `ok = 0`,
  so nothing announces it.
- **Why:** `DisplacementControl` is formulated to solve for `lambda` such that the control DOF
  moves by `theIncrement`. At `theIncrement = 0` the constraint "move the control DOF by zero"
  is satisfied by the trivial solution only if the residual is already zero; otherwise the
  method is free to buy that zero displacement with an arbitrarily large load factor.
- **Workaround/status (2026-09-14):** to hold a converged state, switch to **`LoadControl 0.0`**
  — it freezes the load factor and leaves the control DOF free, which is what "hold" actually
  means here. Note the consequence: a `DisplacementControl` target is **released** while you
  hold (measured 1.08 % drift on the WP-101 gate while the AL multipliers tightened), so read
  the control displacement back afterwards rather than assuming it. See
  [[LadrunoKinematicCoupling_guide]] §4.3 / [PR #839](https://github.com/nmorabowen/OpenSees/pull/839).
