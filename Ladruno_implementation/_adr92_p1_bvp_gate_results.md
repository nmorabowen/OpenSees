---
title: "ADR 92 / P1 gate 1 — the BVP gate: the ladder is gone, and the leg reached its target"
project: Ladruno
type: results
status: "SUPERSEDED by the red/blue review — headline claims withdrawn pending re-run; see _adr92_p1_redblue_review"
priority: high
owner: nmora
related:
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[_adr92_p1_execution_plan]]"
  - "[[_adr92_cp1_surcharge_results]]"
  - "[[_adr92_p0_oracle_results]]"
tags: [adr, sanisand, implex, bvp-gate, measurement]
updated: 2026-09-06
---

# ADR 92 / P1 gate 1 — the BVP gate

> [!danger] **Read `_adr92_p1_redblue_review.md` first (2026-09-06).** The adversarial
> review found that (1) the extrapolation factor ran at `f ≡ 1` on this deck (ratio gated
> `> 0.0` on a negative clock), (2) the `-implexControl` floor was dead for the same reason,
> so the ctl seizure is a bug not a tolerance, (3) the committed gate returns PARTIAL on the
> registered arm, (4) "not one step left rung 1" is false as written (25/25/25 rung entries
> on the ctl arm, all material refusals), (5) `tail = 95.9 %` is a 4-point-fit artefact
> (21.8–29.5 % on a matched window), and (6) over the checkable overlap the control-OFF
> curve is 11.4 % off on average and +38 % at s/B 0.02. The table below is retained as the
> record of what was claimed; the review's §1 and §4 state what the data licenses.

> [!success] **The ladder claim is CONFIRMED. Nothing else is.**
> Three legs, identical deck (`h1.0_e0.6944`, `Q = 10 kPa`, `-maxSubsteps 20000`,
> R3 domain), one flag apart:
>
> | | control | `-implex` + control 0.05 | `-implex`, control OFF |
> |---|---|---|---|
> | steps | 64 | 110 | **142** |
> | rung 1 / rung 2 / rung 3 | 33 / 15 / 16 | **110 / 0 / 0** | **142 / 0 / 0** |
> | past rung 1 | **48.4 %** | **0.0 %** | **0.0 %** |
> | subdivisions | 0/80 | 25/80 | **0/80** |
> | termination | WALL @ `s/B` 0.0553 | **FLOOR @ 0.00504** | **TARGET @ 0.25000** |
> | tail % of initial | 45 | — | **95.9** |
>
> **Not one step in either IMPL-EX leg left rung 1.** The control needed a fallback
> rung on 48.4 % of its steps. A frozen SPD operator makes the step linear and a
> linear step has no ladder to fail — IMPL-EX's structural claim, measured on the
> boundary-value problem for the first time. P0 could not test this; CP1 could only
> show that the ladder was where the wall clock went.
>
> **And the arm with control OFF reached `s/B = 0.25000` on mode TARGET** — the first
> leg in this entire campaign (GATE U, T8, CP1, the deep legs) to finish its push
> rather than stop on a wall, a budget or a floor. Every prior leg died below 0.10.

> [!warning] **Three things this does NOT say, each of which must be quoted with it**
> 1. **It still did not plateau.** `tail = 95.9 %` of the initial slope at `s/B = 0.25`.
>    Reaching the target is not reaching a capacity; `q_u` remains unmeasured, exactly
>    as at GATE U. What changed is the *reason* — no longer "the solver stopped", now
>    "the load is still climbing at a quarter of the footing width".
> 2. **That leg ran with `-implexControl` OFF**, so nothing policed the extrapolation
>    error. P0 measured IMPL-EX-A unusable from `d eps = 5e-4` at `p0 = 5 kPa`, and
>    this deck's free-surface corner is that regime. **A fast curve that nothing is
>    checking is the failure mode this campaign already met once**, on the dynamic
>    relaxation legs that passed every gate on a wrong path.
> 3. **With the control ON at `tol = 0.05` the leg seizes at `s/B = 0.00504`** — a
>    tenth of the control's depth — on 25 refusals down to the DS_MIN floor. So the
>    two arms bracket the truth: **the mechanism works, the tolerance as set does not.**
>    The operating point is somewhere between "refuses everything" and "polices
>    nothing", and it has not been found.

**Provenance, stated exactly.** Engine built from `wp/92c-implex-p1` after the D2
sign-guard fix; the binary **contains** that fix but **reports** hash `8fe4f5630`,
because the build was configured before the fix was committed (`3c788778f`). The fix
was verified **behaviourally, not by hash**: the negative-`dt` refusal fired 4896
times before and **0** times after. The baseline leg reports `80e65a4de`, one commit
earlier, differing by the same fix plus tests. Driver `sanisand_tau0_band.py --implex`,
gate `adr92_bvp_gate.py`, outputs in `adr92_bvp/`.

---

## 1. What the gate was, and that it could have failed

Written before the run, as a prediction with a refutation clause: *past rung 1 falls
from ~48 % to near zero and the failed-rung share to single digits; if the ladder
still fires at that rate the cost case is refuted and P1 closes as a correctness-only
feature.* Self-tested by feeding it CP1's non-IMPL-EX legs, where it returned
**REFUTED** — the instrument could produce the negative result.

## 2. The one number in the table that is NOT meaningful

Arm A's `failed-rung iterations = 85.0 %`. It cannot be a ladder cost: **no step in
that leg left rung 1.** It is an artefact of the gate charging every subdivided step
the full three-rung cap, which is right for a step that *ground* through the ladder
and wrong for one the material *declined* on its first evaluation.

The cost model is now split (`ladder_fails` vs `refused`), but it needs a per-leg
refusal count the driver does not yet record, and the material's refusal message is
throttled to 10 per process so the engine log cannot supply it either. **Rather than
manufacture a number, arm A's failed-rung column is reported as not meaningful.**
Owed: `n_material_refused` in the leg payload. Arm B is unaffected — it has no
refusals and no subdivisions, so its `0.0 %` is exact.

## 3. What this licenses

**MAY say:**
- On this deck, `-implex` removes the solver ladder entirely: 0 of 252 steps across
  both arms needed a fallback rung, against 48.4 % of the control's.
- With control off it completes a `0.25 s/B` push in 142 steps with zero
  subdivisions — the campaign's first completed push.
- The ADR §2 item-1 cost claim, "warranted by CP1 but unproven at BVP level", is now
  **proven at BVP level for the ladder term**.

**MUST NOT say:**
- That a capacity, plateau or `q_u` has been measured. It has not; `tail = 95.9 %`.
- That `-implex` is usable as configured. The only leg that ran far had no error
  control, and the controlled one seized ten times shallower.
- Anything about wall-clock speedup: the arms had different terminations and
  different work, and no timing comparison at matched settlement was made.
- That the deeper reach is physically right. It is an equilibrium of the
  **extrapolated** stress and must be confirmed on the implicit material — which
  currently cannot get past `s/B = 0.055` to do the confirming. That asymmetry is
  itself a finding and is the ADR §8 reading hazard, live.

## 4. Owed next, in order

1. **Find the `-implexControl` operating point.** A tolerance sweep (0.05 is
   unusable; try 0.2, 0.5, 1.0) against depth reached and `implexError` at the corner.
   Until then `-implex` ships either unpoliced or unusable, and neither is shippable.
2. **`n_material_refused` in the leg payload**, so §2's column can be computed rather
   than disclaimed.
3. **A confirmation strategy for a curve the implicit material cannot reach.** The
   ADR requires implicit confirmation up to the last settlement the implicit solver
   reaches; that is `s/B ≈ 0.055` against IMPL-EX's 0.25. What licenses the other 78 %
   of the curve is an open question and the ADR does not currently answer it.
4. **Oracle parity (gate 2)** — still not run.

## Log

- 2026-09-06 — Gate 1 run in two arms after the D2 sign-guard fix (itself found by
  this gate's first attempt, which died at `s/B = 0.00000`). Ladder claim confirmed;
  no plateau; control unusable at 0.05.
