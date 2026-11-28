---
wp: WP-260
title: "TRBDF2 uses different effective-stiffness weights on its two legs, and LadrunoFactorReuse keys only on dt"
date: 2026-11-25
---
### `TRBDF2` uses different effective-stiffness weights on its two legs — `LadrunoFactorReuse` keys only on `dt`, so it reuses a TR factor for a BDF2 step

- **Bites:** `TRBDF2` alternates full trapezoidal (TR) and BDF2 steps across separate `analyze()` calls, and the two legs have different effective stiffness (`c2 = 2/dt, c3 = 4/dt²` vs `c2 = 1.5/dt, c3 = 2.25/dt²`, `TRBDF2.cpp:128-146`). At constant `dt`, `LadrunoFactorReuse` keeps the factor across the switch (its key is `(SOE, integrator, neq, dt)`, `LadrunoFactorReuse.cpp:494-498`), so every other step iterates on a factor of the wrong `K_eff`: slower convergence or more failed steps, not a wrong converged answer (the residual is exact). Separately, any `dt` change (`:109`) or a `revertToLastStep()` (`:179` sets `step = 1`) restarts the pair with a TR leg, so an adaptive-`dt` run that changes `dt` every step never takes a BDF2 leg and gets no numerical damping.
- **Why:** the two-leg scheme is spread over two `analyze()` steps instead of one composite step, and the factor-reuse key has no way to see the integrator's coefficients.
- **Workaround / status (2026-11-25, WP-260 — disclosed):** WP-260 Contract #2 (`IncrementalIntegrator::getEffectiveTangentKey`) lets `LadrunoFactorReuse` key on `(c1, c2, c3)`; `TRBDF2`/`TRBDF3` implement it in WP-260d. Until then, do not pair `TRBDF2` with `LadrunoFactorReuse`. An implicit Bathe replacement is closed by ADR-77 G3 unless the owner re-opens it (WP-260 decision D6).
