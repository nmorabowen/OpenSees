---
wp: LEGACY
title: "The Newton must run in the surface's OWN variable, not in y1"
legacy_seq: 412
---
### The Newton must run in the surface's OWN variable, not in `y1`

The Hoek-Brown surface exists only for `arg = s - mb*y1/sigma_ci >= 0`. With `y1` as
the Newton unknown the first step from the elastic predictor OVERSHOOTS (measured
`arg = -2.2456e-03` at iteration 1 on the oracle's own near-apex trial) and the next
Jacobian is singular. Substituting `arg = (w^2)^(1/a)` — so
`y1 = T - (sigma_ci/mb)(w^2)^(1/a)` and `f = y1 - y3 - sigma_ci*w^2` — is
polynomial-smooth and feasible for ANY real `w`: no clipping, no line search, no
feasibility guard. Write `(w*w)^(1/a)` rather than `w^(2/a)`: the latter is NaN for a
negative iterate, and the two agree for `w > 0`.

Related: NORMALIZE the flow direction in the residual (multiplier rescaled by `|m|`,
the returned stress is invariant). `|m| ~ arg^(a-1)` blows up exactly where the
near-apex returns land — 460.6 there against 3.9 on an ordinary face point — and the
worst-case Newton count over the oracle's 400-trial scan is **6** un-normalized and
**5** normalized.
