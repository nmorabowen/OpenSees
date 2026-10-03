---
wp: WP-134
title: "The ModifiedEuler \"err = 0 path\": a loading stage read as elastic + the uncapped q makes the substep error EXACTLY zero (WP-134 → WP-129)"
legacy_seq: 495
---
### The ModifiedEuler "err = 0 path": a loading stage read as elastic + the uncapped `q` makes the substep error EXACTLY zero (WP-134 → WP-129)
- **Bites:** when both Heun stages take the `dγ < 0` "elastic" branch (finding F above -- a LOADING stage with a negative denominator, or genuine unloading), both stress increments are the same elastic increment, the stress-only error is exactly 0, `q = +inf`, and the next substep covers the rest of the increment while α is dragged along with the stress ratio. No TolE can reject it (WP-134: 20-65 % stress errors on benign 20-100 kPa states even at TolE 1e-8; all 25 campaign-ME ring escapes from admissible starts).
- **Workaround/status:** SAS-ME: stage classification from the elastic-trial numerator (a loading stage with H ≤ 0 is refused / cut, never elastic), α and z in the error, `err ≥ DBL_EPSILON`, `q ≤ 1.1`. ModifiedEuler unchanged.
