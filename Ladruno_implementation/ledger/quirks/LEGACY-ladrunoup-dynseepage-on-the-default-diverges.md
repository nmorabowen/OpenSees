---
wp: LEGACY
title: "LadrunoUP -dynSeepage on (the default) DIVERGES under Δt-refinement in quasi-static consolidation runs — the ü-term feeds integrator noise into the seepage sou…"
legacy_seq: 173
---
### LadrunoUP `-dynSeepage on` (the default) DIVERGES under Δt-refinement in quasi-static consolidation runs — the ü-term feeds integrator noise into the seepage source
- **Bites:** ZS84-class consolidation column, Newmark γ=0.6: with `-dynSeepage off` the error converges 1.3e-3 → 7e-4 as Δt shrinks 0.08 → 0.005; with the default `on` it GROWS 1.8e-2 → 8.7e-1. Smaller Δt is WORSE: trial accelerations of numerically-damped compressible-wave modes are noise, and f_seep integrates them.
- **Why:** the dynamic-seepage drive (b − ü) is physically right for genuine dynamics (B5-class, P4-gated) but quasi-static consolidation has no meaningful ü — the term is pure noise amplification there.
- **Rule (AMENDED at P4):** the default is now **`off`** — the B5 Simon gate measured the failure in genuine dynamics too (wandering post-front p ≈ 1.7–2.0 vs β = 0.973; unbounded shallow-station growth). `-dynSeepage on` is an explicit research opt-in (ADR §12 log 2026-07-13). Companion: `-stab auto` adds ~10% spurious ringing on wave problems — wave runs use `-stab off`. Measured in `tests/test_ladruno_up_element_analytic.py` sweep (ADR-71 P1, 2026-07-11).
