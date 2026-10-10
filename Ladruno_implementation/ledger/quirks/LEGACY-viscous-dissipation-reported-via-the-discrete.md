---
wp: LEGACY
title: "Viscous dissipation reported via the discrete work integral ∫f·du is a DIAGNOSTIC, not an exact energy balance"
legacy_seq: 51
---
### Viscous dissipation reported via the discrete work integral `∫f·du` is a DIAGNOSTIC, not an exact energy balance
- **Bites:** `LadrunoBrick::hourglassEnergy()` for uri-viscous returns a committed
  accumulator `hgDissipated += c_visc·Σ q̇·Δq` (work against the FB rate damper).
  For LIGHT damping this tracks the true dissipated energy well (≈82% of the
  hourglass KE recovered over a long run); for STRONG damping the per-step velocity
  collapse in the leapfrog stagger makes `f·Δu` UNDER-count, so "more damping ⇒
  more reported dissipation" is FALSE as measured (it is non-monotone in eps).
- **Rule:** treat it as a monotone, energy-bounded spurious-mode diagnostic
  (GLSTAT-style), and write tests around the robust properties — non-decreasing,
  positive under hourglass excitation, `0 < E ≤ KE_imparted` across eps, exactly 0
  for a rigid/constant-strain velocity (γ⟂linear) — NOT exact energy convergence.
  Learned 2026-06-01.
