---
title: "ADR 92 / P2-9 — Esmeralda dense/loose-refuse arm results (CLOSED — otherwise branch)"
project: Ladruno
type: results
status: "CLOSED 2026-09-08 — dense refuse wall 0.01755 clears refutation (< 0.0169) but misses ship (>= 0.0177) by 0.85%; decision rule's 'otherwise' branch applies: controlIter recorded as a graded guard with measured gain, NOT a default; fixed stays default, P2-8 stays the fallback"
priority: high
owner: nmora
related:
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[_adr92_p2_9_control_informed_f_plan]]"
  - "[[_adr92_p2_9_r3_results]]"
  - "[[_adr92_p2_9_oracle_results]]"
  - "[[LadrunoSANISAND_implex_guide]]"
tags: [adr, sanisand, implex, p2-9, measurement, esmeralda, closed]
updated: 2026-09-08
---

# ADR 92 / P2-9 — Esmeralda dense/loose-refuse arm results

## Provenance

TIMs act response-curve-matrix session, reported 2026-09-08, ledger `ESMERALDA.md` §§106, 109,
115. Engine `179da6ffb` (the P2-9 merge, PR #822), flags `-implexFactor controlIter -implexFloor
refuse`, cap 20000, q10, fork push, vanilla twin `--presidual 1.01`, ds 5e-5 -> 1e-4 floor
1.95e-7, budget 200. Outputs under
`labs/ape/response-curve-matrix/level3/D-L-dl-vt-{dense,gorini}-q10-sp-p29-iter*/coarse/out/`
(`curve.csv`, `run.json`, `argmax_events.csv`, slurm log). Reference legs: 146599 (P2-7c fixed
f, build `4e07ef014`) honest wall 0.01689; 146570 (build `87b9cf846`) 0.01774.

## Results table

| leg | control | end s/B | how | overlay vs implicit twin past s/B 0.005 (mean / max) | refusals | failed attempts | f=0 guards (push) | implexGuards[6] back-offs (push) | it/step | wall s |
|---|---|---|---|---|---|---|---|---|---|---|
| 146607 dense | 0.1/0.01 | 0.01755 | honest wall, refuse floor | +0.28% / 0.42% | 102 | 11 | 8 747 | 2 697 987 (= 28% of 34 560 pts x 276 steps) | 17.2 | 5 620 |
| 146608 dense | 0.01/0.01 | 0.01932 | subdivision budget spent (201) | +1.44% / 1.77% | 8 124 | 201 | 43 303 | -- | 11.3 | 9 484 |
| 146609 loose | 0.1/0.01 | 0.03921 | honest wall, refuse floor | +0.08% / 0.25% | 380 | 127 | 20 375 | -- | 17.0 | 16 998 |
| 146599 dense ref (fixed f) | 0.1/0.01 | 0.01689 | honest wall | +0.31% | 42 545 | 248 | ~38 000 | n/a | -- | 1 933 |

## Extra facts

- The worst `|implexError|` on both tol-0.1 legs (0.069 / 0.066) is the **exempt first plastic
  step after the flip** (P2-7c's un-primed first step, ADR-92 row P2-7c), not a row past tol.
- On 146607's last rows the `f=0` fires climb 6 -> 448 -> 715 per row alongside the back-offs,
  then the control refuses to the floor -- the seat is still the softening-reversal one (the
  P2-2/P2-9 mechanism, unchanged).
- The tol-0.01 variant (146608) reaches further (0.01932) but sits **1.44% HIGH** against the
  twin (stiffer, not closer) -- see "The tol-0.01 finding" below.
- The loose arm (146609) ends at 0.03921, unchanged at 0.039 to the third figure exactly as
  pre-registered: **PASSES, not a finding.**
- TIMs' harness now carries `--implex-factor {fixed,control,controlIter}`, a seven-slot ledger
  with `guard_backoff`, and refuses the control modes without `--implex-control` before the mesh.

## The verdict against the pre-registered rule

Decision rule (`_adr92_p2_9_control_informed_f_plan.md` §2): "P2-9 ships if the dense refuse
wall reaches >= 0.0177 on the twin and no oracle row regresses; otherwise the ADR records the
factor as a graded guard with its measured gain and P2-8's fixed threshold is the fallback."

Measured dense refuse wall (146607) = **0.01755**.

- Refutation bar: `< 0.0169`. `0.01755 > 0.0169` -- **clears refutation.**
- Ship bar: `>= 0.0177`. `(0.0177 - 0.01755) / 0.0177 = 0.85%` short -- **misses ship by 0.85%.**

`0.01755` sits strictly between the two bars, so neither `if` branch of the rule fires and the
**"otherwise" branch is decisive**:

- **P2-9 does NOT ship as a default.** `-implexFactor fixed` remains the default, which is
  already the code state, so **no code change is owed** by this wave -- it is documentation
  only.
- **`controlIter` is recorded as a graded guard with its measured gain** and stays an explicit,
  documented opt-in for the TIMs campaign, not a recommendation for ordinary decks.
- `control` (frozen f*) remains **REFUTED** (fork R3, `_adr92_p2_9_r3_results.md`) and is not a
  candidate in any form.
- **P2-8's fixed threshold (`-implexGuardKp`, listed, not built) is the documented fallback**
  the ADR records if a graded factor is wanted without `controlIter`'s cost.

### Gain, measured

- **Reach:** 0.01689 -> 0.01755, `(0.01755 - 0.01689) / 0.01689 = 3.9%` over P2-7c's fixed-f
  honest wall.
- **Overlay:** comparable to slightly better -- +0.28% mean (146607) vs +0.31% (146599 ref).
- **Refusal churn collapse:** refusals 42 545 -> 102 (a ~417x drop); failed attempts 248 -> 11
  (a ~23x drop).

### Cost, measured

- **Wall time:** ~2.9x on this deck -- 1 933 s -> 5 620 s (`5620 / 1933 = 2.907`).
- **Iterations/step:** 17.2 for 146607 (146599's per-step iteration count is not in this
  table's columns).
- **Explicit correction of the R3 figure:** the fork R3 leg (`_adr92_p2_9_r3_results.md`)
  measured `controlIter` at ~13x the wall time of `fixed` on its own deck. That figure **does
  NOT generalise** to Esmeralda's dense-refuse deck: the measured multiplier here is **2.9x**
  (1 933 s -> 5 620 s), not 13x. Any future citation of "~13x wall" must be scoped to the R3
  deck specifically.

## The tol-0.01 finding

146608 (`control` 0.01/0.01) reaches further than the tol-0.1 dense arm (0.01932 vs 0.01755) by
spending its full subdivision budget (201 failed attempts, the budget cap) rather than refusing
honestly. Its overlay against the twin past s/B 0.005 is **+1.44% mean / +1.77% max** -- HIGH,
i.e. the curve runs **stiffer**, not closer to the twin. TIMs flags this as a finding, not a
number to keep: a tighter control tolerance does not buy a cleaner match, it buys a
budget-exhaustion crawl that happens to land further out while drifting off the twin in the
wrong direction. This variant is not a candidate for any use and is recorded here only as a
measured data point, not as an alternative configuration.

## What this means for use

`-implexFactor fixed` stays the only default and needs no change: it was already the shipped
state, and this wave confirms rather than revises that. `controlIter` is a real, measured,
graded guard -- it collapses refusal churn by two orders of magnitude and buys a few percent of
extra reach at comparable-or-better overlay accuracy -- but its cost (roughly 3x wall time on a
production-scale deck, not the smaller R3 deck's 13x) and its miss of the pre-registered ship
bar by well under a percentage point mean it stays an explicit opt-in for campaigns that are
already paying for dense, deep pushes near a softening seat and that value fewer refused steps
over wall time. It is not a fix for the p = 0 ring (ADR 93) and is not a substitute for P2-8's
simpler fixed-threshold guard, which remains the documented fallback if a future measurement
wants a graded guard without `controlIter`'s Newton-churn and wall-time cost.
