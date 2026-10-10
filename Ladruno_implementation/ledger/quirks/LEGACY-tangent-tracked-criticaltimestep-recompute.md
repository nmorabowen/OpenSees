---
wp: LEGACY
title: "Tangent-tracked criticalTimeStep() (-recompute/-tangent) OVER-reports the safe explicit step by √(E0/E_tan) for any reloadable material — unsafe as an adaptive…"
legacy_seq: 8
---
### Tangent-tracked `criticalTimeStep()` (`-recompute`/`-tangent`) OVER-reports the safe explicit step by √(E0/E_tan) for any reloadable material — unsafe as an adaptive-growth target
- **Bites:** any adaptive-Δt driver that grows the step toward `safety·criticalTimeStep()` on a model with plasticity or damage-with-reload. As the material softens/yields the tangent stiffness drops, so `criticalTimeStep()` returns a *larger* Δt_cr (Δt_cr = 2/ω_max ∝ √(m/k_tan)). Grow into it and the run blows up the instant an elastic-reload wave hits the softened element — because reload happens at the UNDAMAGED modulus E0, not the current tangent.
- **Why:** the pencil sees only the current (tangent) stiffness; a reloadable material can stiffen back to E0 in one step. Measured over-report = √(E0/E_tan): exactly 7.07× on Steel01 b=0.02 (=√(1/b)) and 12.5× on a Concrete01 softening bar (ADR-65 P0 oracle, `adr65_headroom_oracle/`, 2026-07-03). The tangent tracking itself is CORRECT and stays positive on the softening branch (0 bad samples / 2200 steps) — it's the *interpretation as a safe step* that's wrong.
- **Workaround/status (2026-07-03, ADR-65 Route B oracle):** an adaptive driver must clamp Δt growth to the UNDAMAGED/elastic Δt_cr, not the instantaneous tangent — which makes growth headroom ~1.00× (nil) for reloadable nonlinearity, i.e. adaptive *growth* buys nothing there. The safe use of the tangent query is the SHRINK direction (Δt_cr dropping mid-run) + permanent non-reloadable stiffness loss (erosion). See [[65_ladruno_explicit_dt_strategies_adr]] §Route B P0 oracle RUN.
