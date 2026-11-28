---
wp: WP-266
title: "LadrunoUP corot/hypo: Newmark scales the p-row rate block by gamma/(beta dt), the incremental residual by 1/dt"
pr: "#1177"
date: 2026-11-28
status: "open"
---
### LadrunoUP `-geom corot|hypo` under Newmark: Newton on the p-rows converges linearly at 1 − β/γ, whatever the element tangent

- **Bites:** anyone counting Newton iterations on a transient `LadrunoUP -geom corot` or
  `-geom hypo` deck, or judging how good the u-p tangent is from them. Also anyone who expects
  a better element Jacobian (the wave-2 P4 full-Jacobian option) to give quadratic convergence
  under the battery's `integrator Newmark 0.6 0.3025`.
- **Why:** under corot / hypo, `getResistingForceIncInertia` replaces the whole p-row rate term
  by the INCREMENTAL form `[α Σ N_p tr(Δε) dv + (S + αH̃)(p − p_n)] / Δt` (ADR 78 §3.3: the
  chord-defect and mixed-scheme-pump cure). Its derivative carries 1/Δt. The integrator still
  builds the effective tangent as `K + c2·C` with `getDamp`'s p-rows `[Q_curᵀ, S + αH̃]` and
  Newmark's c2 = γ/(βΔt). With γ = 0.6, β = 0.3025 that is 1.98/Δt: the p-row rate block of
  the tangent is twice the residual's derivative. The code comment calls it "a uniform γ/β
  factor ... absorbed as an inexact-tangent term". In a consolidation step where the rate
  block dominates the p-rows, Newton contracts the residual by about 1 − β/γ = 0.50 per iteration.
- **Measured (WP-266):** the `Dual<32>` oracle (`tests/ladruno_up_hypo_dual_check.cpp`, table row
  "Newmark c2") puts this at ‖K_N − K‖/‖[J_pu J_pp]‖ = 0.94–1.00 on every state. That is larger
  than every dropped Jacobian block: 3× the fast-draining state's, 50–5000× the others. On the large-strain consolidation deck
  (`tests/_testbed/up_hypo_consolidation_newton.py`, pre-WP engine f0b9a1df), Newmark(0.6, 0.3025)
  needs 1045 Newton iterations over 40 steps and Newmark(0.5, 0.5) needs 385. The per-iteration
  residual ratio is 0.49 with the first and 0.22 with the second. The 0.22 is what the dropped
  blocks cost.
- **Workaround:** for quasi-static consolidation use an integrator whose c2 is 1/Δt, such as
  `integrator Newmark 0.5 0.5` (γ = β; unconditionally stable since 2β ≥ γ ≥ ½, with no numerical
  damping and a large period error, so a quasi-static choice).
  Expect the inertia residual to change slightly with it (settlement −0.3839 vs −0.3839 on the
  deck, 0.4 % after 6 steps).
- **Open:** a fix belongs to the element (scale the p-row rate block of `getDamp` by
  (1/Δt)/c2 when the integrator's factors are known) or to the wave-2 P4 full-Jacobian option.
  Neither is in WP-266. Tracked in issue #1180.
