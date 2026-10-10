---
wp: LEGACY
title: "An explicit SOFT/Courant penalty k ∝ 1/dt² keeps the contact event at a CONSTANT step count — energy error converges in SOFSCL, not dt"
legacy_seq: 26
---
### An explicit SOFT/Courant penalty `k ∝ 1/dt²` keeps the contact event at a CONSTANT step count — energy error converges in SOFSCL, not dt
- **Bites:** validating the ADR-39 B1 SOFT=1 penalty's energy balance. The instinct is to
  refine `dt` and watch the impact restitution `e → 1` (energy conservation). It does NOT
  improve — `1−e` is flat across `dt`. The naive read is "the penalty leaks energy."
- **Why:** SOFT sizes `k_soft = SOFSCL·4·m_eff/dt²`, so the contact period
  `T_contact = 2π√(m_eff/k_soft) ∝ dt`. The steps spanning a contact event,
  `T_contact/dt = π/√SOFSCL`, are **independent of dt** — refining dt shrinks the period
  and the step in lockstep, so the contact is always resolved by the same ~`π/√SOFSCL`
  steps. The discrete one-sided-contact engagement error (EITHER sign at coarse
  resolution — the chatter D2 viscous damps) is therefore a function of **SOFSCL**, not dt.
- **Workaround/status (2026-06-24):** test energy convergence by refining **SOFSCL** (not
  dt): `proto_b1_soft_penalty.py` T3/T4 sweep SOFSCL∈{0.1, 0.025, 0.00625} and show
  `|1−e|`, `|ΔKE/KE₀|` → 0. Frame it as "bounded & SOFSCL-convergent (not a formulation
  leak)" — NOT "conservative at the shipped SOFSCL=0.1" (there the bounded error is ~1–2%,
  sign-indefinite; that's the chatter, by design left for `-visc`). SOFSCL is the
  accuracy/stability knob: smaller = stiffer + better-resolved + less penetration.
