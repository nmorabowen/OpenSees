---
wp: WP-260
title: "TRBDF2 breaks factor reuse across substeps and silently degrades to TR on dt changes"
date: 2026-11-25
---
### `TRBDF2` breaks factor reuse across substeps and silently degrades to TR on `dt` changes

- **Bites:** `SRC/analysis/integrator/TRBDF2.cpp` cannot share LU factors across its two substeps, silently drops to undamped Trapezoidal Rule whenever `dt` changes, and corrupts its substep state if `revertToLastStep()` is called.
- **Why:** `TRBDF2` alternates between TR and BDF2 across separate `analyze()` calls using different effective stiffness weights (`c2=2/dt, c3=4/dt^2` vs `c2=1.5/dt, c3=2.25/dt^2` at `:128-146`), resets `step = 0` whenever `deltaT != dt || step == 1` (`:109`), and hardcodes `step = 1` in `revertToLastStep()` (`:179`).
- **Workaround / status (2026-11-25, WP-260 — disclosed):** use `LadrunoImplicitBathe` ($\rho_\infty$-Bathe with optimal $\gamma_0(\rho_\infty)$, identical $c_2,c_3$ across both substeps, composite step execution, and `getEffectiveTangentKey`) for composite implicit wave propagation and factor reuse.
