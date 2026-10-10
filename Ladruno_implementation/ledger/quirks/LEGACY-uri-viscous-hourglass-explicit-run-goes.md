---
wp: LEGACY
title: "uri viscous hourglass: explicit run goes silently unstable at large eps + dt (CDL doesn't trap it)"
legacy_seq: 50
---
### uri `viscous` hourglass: explicit run goes silently unstable at large eps + dt (CDL doesn't trap it)
- **Bites:** `LadrunoBrick -hourglass viscous` adds a velocity-proportional damping
  force but NO hourglass stiffness. Under CentralDifferenceLadruno the viscous term
  has its own explicit stability bound; at `eps≈0.5` with `dt = 0.1·le/c` the
  hourglass mode blew up to `nodeDisp ~ 1e+99` — yet `analyze()` still returned 0
  (CDL does not check for NaN/Inf), so the run *looks* like it completed. Smaller
  eps tolerates the larger dt; large eps needs a smaller dt.
- **Fix / rule:** for viscous-hourglass explicit runs use a conservative step
  (`dt ≈ 0.02–0.03·le/c`) and modest `eps` (≤0.1). Don't trust a clean `analyze()`
  return alone — assert `isfinite(nodeDisp)`. (Element-level: the viscous tangent
  is rank-deficient ⇒ statics is singular; it is explicit-only by construction.)
  Learned 2026-06-01.
