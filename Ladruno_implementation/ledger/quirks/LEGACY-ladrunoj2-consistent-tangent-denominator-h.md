---
wp: LEGACY
title: "LadrunoJ2 consistent-tangent denominator h assumes NON-softening hardening"
legacy_seq: 54
---
### LadrunoJ2 consistent-tangent denominator `h` assumes NON-softening hardening
- **Bites:** the analytic consistent tangent in `LadrunoJ2Kernel.h` divides by
  `h = dtheta + (2/3)·sig_y'(pbar) − n:Mp` (commented `= −df/ddG > 0`) to form
  `betaNN` and `betaMpN`. For standard hardening (`Hiso,Qinf,Cₖ,γₖ ≥ 0`) `h>0`
  always. But the material accepts arbitrary user params: a **negative `Hiso`
  (linear softening) or `Qinf<0`** makes `sig_y'` negative and can drive `h→0` or
  `h<0` ⇒ an `inf`/`NaN` or a sign-flipped (non-physical) tangent that poisons the
  global Newton with no diagnostic. The local scalar-Newton residual `dR` has the
  same exposure.
- **Why it's not "fixed":** this is **pre-existing** behaviour inherited verbatim
  from the original `integrate()` (the 2026-06-02 kernel extraction was deliberately
  bit-identical, so it was preserved, not introduced). The model is designed for
  hardening/perfect plasticity; softening is out of its intended envelope.
- **Rule:** do not feed `LadrunoJ2` softening hardening params. If softening is ever
  wanted, the right fix is **parameter validation at construction** (reject/warn on
  `Hiso<0`/`Qinf<0` that can violate `h>0`) plus a kernel guard
  (`if (fabs(h) < eps·stressScale)` → fall back to the elastic-predictor tangent,
  mirroring the existing `‖M‖→0` `normFloor` treatment) — a separate PR, not a
  bit-identical-extraction change. Surfaced by the PR #97 adversarial review,
  2026-06-02.
