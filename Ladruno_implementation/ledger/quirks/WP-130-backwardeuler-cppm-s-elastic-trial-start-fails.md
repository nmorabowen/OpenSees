---
wp: WP-130
title: "BackwardEuler_CPPM's elastic-trial start fails exactly where SANISAND is stiffest -- the first plastic increment after a reversal or the stage flip (WP-130)"
legacy_seq: 484
---
### `BackwardEuler_CPPM`'s elastic-trial start fails exactly where SANISAND is stiffest -- the first plastic increment after a reversal or the stage flip (WP-130)
- **Bites:** with `alpha_in = alpha` (the stage flip under `-flipAlphaIn init`, or any loading
  reversal) `(alpha - alpha_in):n = 0`, `GetStateDependent` returns its `h = 1e10` sentinel and
  `NewtonSol` zeroes `dGamma` in the Jacobian; the condensed system is near singular from the
  elastic trial, the local Newton returns -1 and vanilla goes straight to the explicit fallback, or
  returns 0 and halves up to 2^9 times. On the one-quad free-DOF decks of WP-130 the failing step
  is the first plastic push step after the flip.
- **Status (WP-130, #868):** `-cppmStart explicit` retries the local Newton ONCE per level from a
  50-substep ForwardEuler guess (vanilla carries this rung as `SchemeControl == 1`, dead code that
  would also pass uninitialised K, G). On the 100 kPa quad it returns 3 of 4 failing increments
  (census `cppmGuessOk/cppmGuessTries`) but the global step still fails; it is opt-in and
  unqualified beyond the WP-130 decks. `-cppmLineSearch on` (halving on ||R||) cut the local-Newton
  failures on the 10 kPa quad from 36 to 2 but did not make that global step converge either.
