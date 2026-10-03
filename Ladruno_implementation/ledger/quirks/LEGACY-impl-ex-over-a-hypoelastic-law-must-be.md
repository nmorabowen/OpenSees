---
wp: LEGACY
title: "IMPL-EX over a HYPOELASTIC law must be INCREMENTAL, and the extrapolated stress needs the p_min clamp"
legacy_seq: 368
---
## IMPL-EX over a HYPOELASTIC law must be INCREMENTAL, and the extrapolated stress needs the `p_min` clamp

**Found 2026-09-05 — the first by the ADR-92 Fable review before the oracle ran, the second by
the oracle's corner gate (`_adr92_p0_oracle_results` §5, §8).**

- **Incremental, not total.** ASD-style IMPL-EX is usually written `sigma~ = C:(eps - eps_p~)`.
  `ManzariDafalias` integrates `dsigma = Ce(p):deps_e` with moduli at the **committed** stress
  (`elastic_integrator` `:1008-1011`, `BackwardEuler_CPPM` `:2223-2226`); rebuilding the stress
  from a total elastic strain discards the pressure-dependent history. Measured **at a ZERO
  increment**: the total form returns `p = 99.50` on a committed `100.00` (0.5 %) and
  **`1.11` on a committed `5.00` (78 %)**. Its error does not vanish as `dt -> 0`, so a
  convergence gate fails for a reason that has nothing to do with IMPL-EX. Write
  `sigma~ = sigma_n + Ce(p_n):((eps_{n+1} - eps_n) - f*d_eps_p(n))`. Same family as the
  ADR-90 P0b entry above (a wrapper's proof silently assuming a constant elastic operator).
- **Clamp `sigma~`.** Nothing in the extrapolation knows about the floor: on a path driven onto
  `p_min = 0.0101 kPa` the extrapolated mean stress reached **-1.37 / -0.16 / -0.09 kPa** at
  40 / 80 / 160 steps (first order, O(1)-O(10) relative) while the committed state sat at
  `+0.0101`. Apply the code's own device — `sigma~ = dev(sigma~) + p_min*I1` when
  `tr(sigma~)/3 < p_min` — or the free-surface element receives tensile mean stress every
  iteration. Related: the static-`ops_Dt` IMPL-EX entry (§ "IMPL-EX in a STATIC analysis") —
  ADR-92 D2 rediscovered it from the `Domain.cpp:2054` side.
