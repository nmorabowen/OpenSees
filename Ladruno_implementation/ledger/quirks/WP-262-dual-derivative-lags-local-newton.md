---
wp: WP-262
title: "A dual number carried through a local Newton has the derivative of the PREVIOUS iterate"
date: 2026-11-27
status: "open"
---
### Dual-number tangent of a return map is slightly off even though the stress matches: the derivative lags the local Newton by one iterate

- **Bites:** anyone templating a return map with a local Newton (J2 Δγ, HS regimes,
  Concrete3D 3×3 / 4×4, SANISAND CPPM) on `ladruno_dual::Dual<N>` and reading the tangent off
  the derivative part after the loop. The value is exactly the double run's, but the derivative
  is not quite the algorithmic tangent.
- **Why:** one Newton update x ← x − R/R_x, differentiated, gives
  x' ← −R_ε/R_x + R·R_x'/R_x², all evaluated at the iterate it started from. So after the loop
  stops (the test reads only the value of R), the derivative is the implicit-function-theorem
  derivative of the previous iterate, not of the converged one. With quadratic convergence and
  a tight local tolerance the gap is small but real. Measured on the LadrunoJ2 strong-Voce path
  (b = 300, `tests/ladrunoj2_dual_check.cpp`): 7.5e-13 relative tangent error without the fix,
  2.8e-14 with it. That is close to a 1e-12 gate, and it grows with the local tolerance and with
  how nonlinear R is.
- **Fix (WP-262):** after the loop converges, take one more Newton step at the converged point
  with its value part zeroed: `T corr = R/dR; ladruno_dual::zero_value(corr); dG -= corr;`
  inside `if constexpr (ladruno_dual::is_dual<T>::value)`. The derivative becomes −R_ε/R_x at
  the converged point and the value stays bit-identical to the double run
  (`LadrunoJ2Kernel.h`, `returnMapT`). For a vector local Newton the same step needs the
  converged local Jacobian solve (Albany/LCM `computeFADInfo` does exactly this).
- **Related traps** (do not rediscover):
  - **Mixed scalars.** Template argument deduction never converts `double → Dual`, so every
    mixed operator must be written out (`LadrunoDual.h` does it).
  - **Branches and tolerances** must read `value(x)`. The header has no implicit conversion to
    `double`, on purpose.
  - **Bit-identity of the double path** is the first gate when templating a kernel. Compare
    the existing C++ checks' `%.17g` output against the pre-change headers.
