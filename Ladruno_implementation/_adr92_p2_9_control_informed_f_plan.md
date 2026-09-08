---
title: "ADR 92 / P2-9 — control-informed extrapolation factor: plan"
project: Ladruno
type: plan
status: "MEASURED 2026-09-08 — `control` REFUTED (fork R3), `controlIter` PASSES R3, Esmeralda dense-refuse arm owed"
priority: high
owner: nmora
related:
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[93_ladruno_sanisand_zero_confinement_adr]]"
  - "[[_adr93_seat_replay]]"
  - "[[LadrunoSANISAND_implex_guide]]"
  - "[[_adr92_p2_9_oracle_results]]"
  - "[[_adr92_p2_9_r3_results]]"
tags: [adr, sanisand, implex, p2, plan, measurement]
updated: 2026-09-08
---

# ADR 92 / P2-9 — control-informed extrapolation factor

> [!abstract] **Question.** Is the P2-2 guard (`f = 0` after a committed `Kp ≤ 0` or an `α_in`
> reset) general enough? **No**, in three ways measured on Esmeralda (#807 acceptance): it is
> SANISAND-specific in its trigger, it is a threshold at zero that under-fires (the engine that
> fired at ~2× the points by accident reached 10 % further on the twin's curve at no accuracy
> cost: honest wall 0.0177 vs 0.0169), and it is reactive (reads the committed predecessor; a
> trial that first reaches a softening point is extrapolated at full strength, and the P2-6
> retry rescues only when the elastic predictor already passes tol — 322 rescues against
> thousands of refusals). This plan replaces the guard's *degree* and *timing* with a
> closed-form, material-agnostic factor chosen from the companion the control already computes.

## 1. The operator

Under `-implexControl` the companion stress `σ_impl` is computed at every trial. The extrapolated
stress is linear in the factor:

```
σ~(f)  = σ_n + Ce:(Δε − f·Δε_p(n))
σ~(f) − σ_impl = A − f·B,   A = σ_n + Ce:Δε − σ_impl,   B = Ce:Δε_p(n)
```

The factor minimising `‖A − f·B‖` is closed form:

```
f* = clamp( (A·B) / (B·B), 0, f_max ),   f_max = alpha·dt_{n+1}/dt_n (today's f)
```

- `f* → 0` exactly where the history points the wrong way (today's P2-2 case) or is stale;
- `f* → f_max` where the history is right (today's default);
- in between it is a graded version of the P2-8 threshold, with no model-specific trigger.

**Frozen per step.** `f*` is computed from the FIRST iterate's `Δε` and held for the step, so
the global step stays linear (frozen `Ce`) — the property that removed the ladder. Later
iterates reuse it. (Alternative, to be priced: recompute per iterate; costs linearity.)

**Reported.** `implexDetail[5]` = the `f` actually used; a new census slot counts steps where
`f* < 0.5·f_max` (the operator "backed off").

**Material-agnostic.** Uses only `σ_n`, `Ce`, `Δε`, `Δε_p(n)`, `σ_impl` — every IMPL-EX material
holds them. First target `LadrunoSANISAND`; the Lemaitre wrapper can adopt it unchanged.

**Limits.** Without control there is no companion at the trial: the committed-predecessor
guard (P2-2) remains the fallback there. The loose wall at s/B 0.039 is the material's
(ADR 93), untouched by any factor.

## 2. Pre-registered predictions (write the numbers before the run)

| test | prediction | refutes P2-9 if | measured |
|---|---|---|---|
| G0 oracle rows (`adr92_p0_oracle`) | byte-identical when `B·B = 0` (elastic) and where `A ∥ B` | any G0 row changes without a plastic history | **PASS.** `--gate G0`/`--gate G2` byte-identical before/after the variant-D commit (3 seeds, 4.2e-14…3.9e-13); GD.2 (a real elastic-unloading path) `max\|σ_A − σ_D\| = 0.0`. See `_adr92_p2_9_oracle_results.md` §3/GD.2 |
| Seat replay (`_adr93_seat_replay`, step 331) | error 0.46 → ≤ 0.05 (the `f = 0` value was 0.029; `f*` lands at or below it) | error > 0.1 | **PASS, by ~60x.** `f* = 0.058515` → error 0.00082 (`cos(A,B) = 0.9996`). New finding not in this table: `f*` **frozen** on a bad first iterate is 100–450x worse than today's `f` on the G2 sweep; recomputed on the converged `d_eps` it is 1.0–2.1x better. See `_adr92_p2_9_oracle_results.md` §5/GD.4 |
| Fork R3 registered arm (`adr92_bvp_fix/ctl`, tol 0.1) | depth ≥ P2's 0.076, refusals/step ↓, overlay ≤ 2 % | depth < 0.076 or overlay > 5 % | **`control` (Leg 1, `a6a53948e`): REFUTED.** depth 0.0521 < 0.076, overlay 11.12 % > 5 %, `n_guard_ctlf` 196,549 over 366 steps; refusals/step DID fall (17.97 → 2.31) but the curve drifted double digits before seizing. **`controlIter` (Leg 2, `9e73060a9`): PASS.** depth 0.1149 ≥ 0.076, overlay 1.70 % ≤ 2 % (best of any arm this campaign); cost: refusals/step rose to 7.39, wall 1858.4 s (~13x Leg 1's 141.6 s), Newton max-iter stalls 1 → 89. See `_adr92_p2_9_r3_results.md` |
| Esmeralda dense refuse arm (TIMs, q10, fork push, control 0.1/0.01, cap 20000) | honest wall **≥ 0.0177** (the accidental engine's) on the twin within 0.3 % | wall < 0.0169 (P2-7c's) | **owed** — not yet run; per §4 below, run on HEAD (`aeb1da492`) with `-implexFactor controlIter` (the only mode that cleared R3) |
| Esmeralda loose arm | wall unchanged at 0.039 (material) | — (a change there would be a *finding*, not a pass) | **owed** — not yet run, same handoff as the dense arm |

Decision rule: P2-9 ships if the dense refuse wall reaches ≥ 0.0177 on the twin and no
oracle row regresses; otherwise the ADR records the factor as a graded guard with its measured
gain and P2-8's fixed threshold is the fallback. **As measured so far:** `control` (frozen f*) is
REFUTED by the fork R3 arm and is not a shipping candidate; `controlIter` (per-trial f*) is the
only mode to clear R3, at a measured wall-time/Newton-churn cost, and the Esmeralda dense-refuse
arm on `controlIter` is now the decision rule's outstanding input.

## 3. Lanes

- **Oracle (numpy, first):** implement `f*` in `sanisand_implex_oracle.py` as variant D; run
  G0, the G2 rows at p0 = 5/100, and the seat replay; report `f*` along the seat path.
- **C++:** `-implexFactor fixed|control` (default `fixed` until the gate passes) in
  `LadrunoSANISAND`; `f*` at the first iterate under control; census slot; wire; echo.
- **Tests:** a test where the history points the wrong way (a reversal) must give `f* ≈ 0`
  and one where it is right must give `f* = f_max`; byte-identity under `fixed`.
- **Fork regression:** the R3 registered arm; then TIMs' pair (they run it the hour the hash
  lands; their harness reads `implexDetail[5]` per point already).
- **Docs:** ADR 92 row; guide; quirks if any.

## 4. Handoff

State lives in memory (`ladruno-adr92-sanisand-implex`) and in ADR 92's P2 table; the TIMs
agent (session `nonlinear-response-curve-planning-2b60bc-95`, ledger ESMERALDA.md §41–99)
has the acceptance harness and the reference numbers. Branch `wp/92f-implex-control-f`.

**Handoff to TIMs:** run the dense refuse arm on HEAD (`aeb1da492`) with `-implexFactor
controlIter` — the only mode that cleared the fork R3 bars (§2 above); their harness already
reads `implexDetail[5]` per point, no change needed there; `implexGuards` is now 7 slots
(`[6]` is the new P2-9 back-off census, fires in both control modes).

## Log

- 2026-09-07 — opened after #807/#801 merged; no code.
- 2026-09-07 — C++ landed: `-implexFactor fixed|control` (`a6a53948e`), the closed-form
  `f* = clamp((A:B)/(B:B), 0, f_max)` frozen at the first trial of the step, gated on
  `-implexControl`.
- 2026-09-07 — oracle lane run (`_adr92_p2_9_oracle_results.md`): G0/G2 byte-identical,
  seat prediction PASS at 0.00082 (bar was ≤ 0.05); found and priced the freeze-on-a-bad-
  predictor collapse (frozen f* up to 450x worse than `fixed` on the G2 sweep; recomputed
  on the converged `d_eps` it is 1.0–2.1x better) — a failure mode the plan did not name.
- 2026-09-08 — fork R3 registered arm, Leg 1 (`a6a53948e`): `-implexFactor control`
  **REFUTED** (depth 0.0521 < 0.076, overlay 11.12 % > 5 %) — the oracle's freeze bias
  reproduced on a real BVP. `_adr92_p2_9_r3_results.md`.
- 2026-09-08 — `-implexFactor controlIter` added (`9e73060a9`): f* recomputed at every
  Newton trial from that trial's own `d_eps` instead of frozen at the first; `f_max` still
  stored once per step, back-off census still counted once per step.
- 2026-09-08 — fork R3 registered arm, Leg 2 (`9e73060a9`, restamped to `aeb1da492`):
  `-implexFactor controlIter` **PASSES** (depth 0.1149 ≥ 0.076, overlay 1.70 % ≤ 2 %) at
  ~13x Leg 1's wall time and a jump from 1 to 89 Newton max-iteration stall markers.
  `_adr92_p2_9_r3_results.md`.
