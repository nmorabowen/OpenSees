---
wp: LEGACY
title: "cot(phi) is Inf at phi == 0, and Inf * 2*sin(phi) is the IEEE-754 indeterminate NaN -- phi == 0 is a NORMAL, legitimate input (undrained clay), not an edge case"
legacy_seq: 417
---
### `cot(phi)` is `Inf` at `phi == 0`, and `Inf * 2*sin(phi)` is the IEEE-754 indeterminate NaN -- `phi == 0` is a NORMAL, legitimate input (undrained clay), not an edge case
- **Bites:** `StiffSoilShear_YF.h`'s `qf = (c*cot(phi) + sigma3) * 2*sin(phi) / (1 - sin(phi))` NaNs on the very FIRST yield-function call for any deck with `MC_phi = 0.0` -- a completely standard cohesive-only / undrained (`phi_u = 0`) clay model, not a degenerate corner case anyone would think to guard against by inspection. `cot(0) = cos(0)/sin(0) = 1/0 = Inf`; the SAME expression then multiplies that `Inf` by `2*sin(phi) == 0`, and `Inf * 0` is the one indeterminate form IEEE-754 refuses to resolve to anything but NaN.
- **Why:** the formula was transcribed from a geotechnical reference in `cot`-form without noticing that `cot` is exactly the term that blows up at the model's most common "no friction" configuration.
- **Fix (ADR-97 P5, `StiffSoilShear_YF.h`):** multiply the ORIGINAL expression through by `sin(phi)` algebraically -- `qf = 2*(c*cos(phi) + sigma3*sin(phi)) / (1 - sin(phi))` -- which removes `cot` entirely, is numerically identical to the old formula for any `phi != 0` (matched to ~1e-14 relative in a standalone probe), and gives the physically correct Tresca limit `qf -> 2c` as `phi -> 0` instead of NaN. **General lesson:** any `cot(x)`/`tan(x)`/`1/sin(x)` term that later gets multiplied by `sin(x)` (or a factor containing it) in the SAME expression is a candidate for exactly this bug -- always check whether the singularity cancels algebraically before accepting it at `x == 0`.
