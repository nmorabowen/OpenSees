---
title: "ADR 92 / P1 gate 1 — the BVP gate: the ladder is gone, and the leg reached its target"
project: Ladruno
type: results
status: "WITHDRAWN 2026-09-06 — see _adr92_p1_bvp_gate_rerun"
priority: high
owner: nmora
related:
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[_adr92_p1_execution_plan]]"
  - "[[_adr92_p1_redblue_review]]"
  - "[[_adr92_p1_bvp_gate_rerun]]"
  - "[[_adr92_cp1_surcharge_results]]"
  - "[[_adr92_p0_oracle_results]]"
tags: [adr, sanisand, implex, bvp-gate, measurement, withdrawn]
updated: 2026-09-06
---

# ADR 92 / P1 gate 1 — the BVP gate (WITHDRAWN)

> [!danger] **Every headline claim below is withdrawn.** The adversarial review
> (`_adr92_p1_redblue_review.md`, §1 B1-B3) found the C++ this data was measured on was
> defective on the campaign's own decks — the extrapolation factor ran at `f ≡ 1`, the
> `-implexControl` floor was dead, and the companion refusal was silently swallowed by
> `Domain::commit()` — and separately (§1 B4) that this memo's prose overstated what
> even that broken run's own data showed. Current numbers, on the fixed binary, are in
> `_adr92_p1_bvp_gate_rerun.md`. This file is kept only as the record of what was
> claimed and why each claim fell.

## What was claimed, and what the review found (B4)

| claimed | what the same run's data actually shows |
|---|---|
| "PREDICTION MET" | Feeding the gate the **registered** arm (`-implex -implexControl`) alone returns **PARTIAL** (0.0 % / 85.0 %; both clauses must be <= 10 %). "MET" is reachable only by scoring the **unregistered** control-OFF arm instead. |
| "Not one step in either IMPL-EX leg left rung 1" | False as written: the ctl engine log records **25 / 25 / 25** `solveCurrentStep` failures on NR / NLS / AcceleratedNewton (`nfail = 75 = 3 x 25`). What survives: zero of those 75 were `CTestNormUnbalance` failures — all were material refusals, so the ladder-as-*iteration-cost* claim holds even though the sentence as written does not. |
| (undisclosed) | **Three pre-registration deviations**: the headline arm run was unregistered (control OFF, not the pre-registered `-implexControl` arm); the subdivision cap was raised 1000 -> 20000 without flagging it as a variable (a sibling commit calls it "controlled"; `adr86b_t8` shows a 38 % swing in `q_u` across caps on this mesh); the baseline "past rung 1" figure was restated 61 % -> 48.4 % inside the very sentence claiming to quote a number written before the run. |
| "warranted by CP1 but unproven at BVP level" (quoted as an ADR line) | **Occurs nowhere in the ADR.** Section 2 item 1 (constitutive work) was never measured by this gate; the gate measured item 2 (the linear step), which the ADR already states as a derivation, not an open question. |
| `implexError` reported | **Never recorded.** ADR §8 mandates it "beside every WP1 verdict"; the driver on this run has zero references to it and no artefact carries a value. |
| `tail = 95.9 %` of initial slope | **Inflated 3.3-4.4x.** The 4-point fit spans 60 um and returns **-6709** on the ctl arm (an undisclosed fit failure). A matched-window refit (`s <= 1 mm`) agrees across arms to 12.6 % and gives tails of **21.8-29.5 %** — still no plateau, at the honest number. |
| overlay never compared | Over the checkable overlap (`s/B <= 0.0553`) the control-OFF curve runs **11.4 % off on average, +38 % at s/B 0.02**, non-monotone (back to +0.7 % by 0.04), carrying 6.6-9.4 % *less* plastic strain while carrying up to 38 % *more* load. Step 1 is 2.05x the control (pure elastic predictor); step 2 *drops* q by 18.6 % under increasing settlement. |

## Why this memo was withdrawn rather than corrected in place

The C++ it measured had three confirmed defects (B1-B3) that change what the numbers
*mean*, not just how they are described — a fixed binary was required before any
re-run could be trusted, and the fixed-binary run mixes a different build per arm for
a legitimate reason (see the rerun memo's build-mismatch callout). See
`_adr92_p1_redblue_review.md` §5 for the required order of fixes, and
`_adr92_p1_bvp_gate_rerun.md` for the re-run itself, including the registered-arm
verdict (`PARTIAL`) this memo never reached.

## Log

- 2026-09-06 — Withdrawn following the red/blue adversarial review (5 blockers,
  B1-B5). Body rewritten as a record of the original claims and their refutation;
  superseded by `_adr92_p1_bvp_gate_rerun.md` for current numbers.
- 2026-09-06 (earlier) — Gate 1 run in two arms; the ladder-removal claim itself
  stands (confirmed again on the fixed binary, both arms — see the rerun memo and
  the ADR §2 item 2 update); the plateau, tolerance-operating-point and cost claims
  did not survive review.
