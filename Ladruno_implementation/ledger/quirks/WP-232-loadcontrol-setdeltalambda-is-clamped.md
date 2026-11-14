---
wp: WP-232
title: "LoadControl::setDeltaLambda does not set the next step: newStep clamps the increment to [dLambdaMin, dLambdaMax], which default to the user's dLambda"
date: 2026-11-14
---
### `LoadControl::setDeltaLambda` does not set the next step — `newStep` clamps Δλ to `[dLambdaMin, dLambdaMax]`, and both default to the user's Δλ

- **Bites:** a driver that cuts or sets a quasi-static step through `setDeltaLambda` (WP-232f's MPM quasi-static `mpm analyze`, the retry's cut): the call returns 0, the next `newStep` multiplies by `specNumIncrStep / numIncrLastStep` and then clamps to `[dLambdaMin, dLambdaMax]`. `integrator LoadControl $dLambda` without the optional arguments sets `min = max = dLambda` (`OPS_LoadControlIntegrator`, `LoadControl.cpp:72-86`; classic Tcl `commands.cpp`, the `else` branch), so the step silently stays the user's Δλ.
- **Why:** `newStep` is `deltaLambda *= factor; if (deltaLambda < dLambdaMin) deltaLambda = dLambdaMin; else if (deltaLambda > dLambdaMax) deltaLambda = dLambdaMax;` (`LoadControl.cpp:119-125`); the members are private with no getter.
- **Workaround / status (2026-11-14, WP-232f):** `LoadControl.h` gains `ladrunoGetIncrement` / `ladrunoSetIncrement` (Δλ, min, max, vanilla-ledgered); the MPM driver reads the user's three values, steps with `min = max = dt`, and restores them when `mpm analyze` returns, printing the override once.
