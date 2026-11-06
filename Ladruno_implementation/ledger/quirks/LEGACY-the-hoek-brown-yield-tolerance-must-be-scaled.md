---
wp: LEGACY
title: "The Hoek-Brown yield tolerance must be scaled by the GRADIENT, not by sigma_ci"
legacy_seq: 411
---
### The Hoek-Brown yield tolerance must be scaled by the GRADIENT, not by sigma_ci

`|df/dy1| = 1 + a*mb*arg^(a-1)` **diverges** at the apex (`arg -> 0`), so a last-ulp
error in the returned `y1` carries `eps*|y1|*|df/dy1|` into `f`: an ABSOLUTE 1e-10
admissibility gate is unattainable within ~1e-2 kPa of the vertex. The P0 oracle
measures `max|f| = 1.08e-10` against its own round-off floor of `1.16e-08` there.
This is the Hoek-Brown instance of ADR-94 M5's `f_relative_tol` lesson, and it is
sharper: scaling by `sigma_ci` (or by `strength_scale`) is not enough on its own,
because the offending factor is the gradient's own conditioning.
