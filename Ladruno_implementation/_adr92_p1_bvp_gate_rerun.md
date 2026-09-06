---
title: "ADR 92 / P1 gate 1 rerun — control and uncontrolled arms on the fixed binary"
project: Ladruno
type: results
status: "All three arms COMPLETE; operating-point sweep COMPLETE. Registered tolerance (0.05) is too tight: reaches only s/B 0.0275 vs control depth 0.0678, regardless of reductionLimit (0.01 vs 0.1 are byte-identical on this deck). tol=0.1 is the tightest sweep point that clears BOTH bars: depth 0.0764 (> control 0.0678) with mean overlay deviation 1.87% (< 5%); tol=0.2 and 0.5 also clear both, with 0.5 reaching TARGET outright."
priority: high
owner: nmora
related:
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[_adr92_p1_execution_plan]]"
  - "[[_adr92_p1_redblue_review]]"
  - "[[_adr92_p1_bvp_gate_results]]"
tags: [adr, sanisand, implex, bvp-gate, measurement, rerun]
updated: 2026-09-06
---

# ADR 92 / P1 gate 1 rerun — control and uncontrolled arms

Pinned build `2473ce46c2ae44412246ca305f3843463f389cad` (`LADRUNO_A2_EXPECT_BUILD`
enforced, not `any`; lane A's pytest battery — 20 passed — ran immediately before this
rerun). Deck: `h1.0_e0.6944`, `--surcharge 10.0`, `--maxsubsteps 20000` (the fix
contract's controlled variable on all three arms; `tests/test_ladruno_sanisand_implex.py`'s
`_CAP_ADEQUATE` comment measured this deck needs 1000-5000 substeps, so the
pre-registered cap of 1000 would confound refusals with cap-hits), R3 domain
(`xlim`/`zbot` defaulted), `--wall 2400`. Three legs, run **sequentially** from this
worktree, output under `Ladruno_files/testbed/hypo_bearing/adr92_bvp_fix/<arm>/`:

1. **control** -- no `-implex`.
2. **ctl** (the pre-registered arm) -- `--implex --implex-control 0.05 0.01`.
3. **noctl** -- `--implex`, control off.

> [!warning] **Sequencing defect during this rerun (disclosed, not hidden).** The
> control leg's process was still running when a checkpoint-JSON existence check was
> read as "finished" and legs 2/3 were started. `ctl` (11:40:29-11:40:31) and `noctl`
> (11:40:59-11:43:07) both executed while `control` (11:38:56-12:20:05) was mid-run --
> an overlap window of **11:40:29-11:43:07 (~158 s)**, timestamps from each leg's
> `run.log` header and the control leg's own progress lines (`s/B=0.005` at internal
> `t=110s` -> clock ~11:40:46; `s/B=0.01` at `t=294s` -> clock ~11:43:50), so the
> contended segment of `control` is roughly its `s/B` 0.005-0.01 stretch. `control`'s
> wall-clock figures in that stretch (and everywhere else, since the process was never
> alone on the box -- see below) are indicative, not clean. `ctl` and `noctl` are
> unaffected in substance (their own step counts/modes do not depend on wall time), but
> their `wall_s` columns share the same caveat.
>
> **Machine-load caveat (independent of the above).** An unrelated `python` process,
> PID 74216, started 2026-09-06 01:05:44 and had accumulated **>11 CPU-hours** (41680 s
> CPU) by the time this rerun ran, was resident on the machine throughout all three
> legs. Every `wall_s` / speedup number below is indicative only, not a clean
> single-tenant measurement.

> [!warning] **The registered arm (`ctl`) ran on a different build than `control` /
> `noctl`.** `control` and `noctl` above ran on `2473ce46c2ae44412246ca305f3843463f389cad`;
> `ctl` (below) ran on `afb95c40c9e3da1cc97d6aabe8168644c53d0e82`, after lane A fixed a
> step-1 refusal defect: the first push step after `updateMaterialStage 1` has no
> plastic history yet, so `implexError` there was measuring the companion's own
> drift-correction jump rather than a genuine extrapolation error; that un-primed step
> is now exempt from `-implexControl` refusal. **The only `SRC/` change between the two
> builds is that exemption, which is gated on the un-primed-step condition and cannot
> touch a control-OFF leg (implex off entirely) or the already-complete `noctl` leg
> (implex on, control off -- nothing to refuse).** `control` and `noctl` were therefore
> **not** re-run; the overlay and gate comparisons below mix builds deliberately, and
> the gate's own build-mismatch WARN (visible below) is expected, not a defect.

## Registered arm (`ctl`, `-implex -implexControl 0.05 0.01`, build `afb95c40c9`)

`adr92_bvp_gate.py --implex adr92_bvp_fix/ctl --baseline adr92_bvp_fix/control
--registered-arm` verdict, verbatim:

```
  --registered-arm: restricted to legs with implex_control set (the pre-registered `-implex -implexControl` arm)
  builds: implex ['afb95c40c'], baseline ['2473ce46c']
  *** WARN: baseline and implex arms ran on DIFFERENT engine builds (['2473ce46c'] vs ['afb95c40c']) -- this gate is not a same-binary comparison ***

  past rung 1               52.2 %  ->     0.0 %   (converged steps only)
  past rung 1 (attempts)    52.2 %  ->    13.8 %   (steps + nsub, apples-to-apples denominator -- review B4/F2)
  failed-rung iters         83.5 %  ->    99.0 %

  VERDICT: PARTIAL -- between the prediction and a refutation. Report the
  table as measured; do not round it toward either.
```

The WARN is the expected build-mismatch flag from the callout above, not a defect --
`control`'s data is unaffected by the fix that separates the two builds. `ctl` reached
504 converged steps (all on rung 1: `nfail = 243` were all incurred inside the 81
subdivided-and-abandoned attempts, none on a converged step) before hitting the
**BUDGET** termination (`subdiv_budget = 80`, `nsub = 81` -- one over) at `s/B =
0.02754`. **0.0 % past rung 1 on the converged-only denominator, 13.8 % on the
attempts-based one (81/585)** -- exactly the coordinator's expectation, and the review
B4/F2 point made concrete: 0.0 % converged-only is a survivorship artefact (every
attempt that entered the ladder past rung 1 also failed all three rungs and was
subdivided away, so it never became a "converged step on rung 2/3" to count), while
81 of 585 attempts (13.8 %) really did enter the ladder. `failed-rung iterations` =
99.0 % is the highest of any arm measured in this campaign -- consistent with `nsub`
now capturing genuine ladder-cost work (81 abandoned 3-rung attempts) rather than the
pre-fix ctl leg's single-step death.

`n_material_refused = 10270` (the process-wide `-implexControl` refusal counter, final
read) against 504 converged steps: **20.4 refusals per converged step** -- i.e. the
control is firing repeatedly and being overridden by retries within a step's ladder,
not just at the 81 abandoned attempts. The curve CSV's own cumulative `implex_refusals`
column reads 10126 at the last converged row, 144 short of the final 10270 -- the gap
is refusals incurred during the terminal (505th) subdivision attempt, after the last
converged row was written but before the BUDGET stop; expected, not a bug (the curve
column is a per-converged-step snapshot, not a running total sampled every attempt).

## Overview

| | control (`2473ce46c`) | `-implex`, control OFF (`noctl`, `2473ce46c`) | `-implex -implexControl 0.05 0.01` (`ctl`, `afb95c40c`) |
|---|---|---|---|
| steps (converged) | 69 | **142** | 504 |
| attempts (steps + nsub) | 69 | 142 | **585** |
| rung 1 / rung 2 / rung 3 | 33 / 15 / 21 | **142 / 0 / 0** | **504 / 0 / 0** |
| nsub (subdivisions) | 0/80 | 0/80 | **81/80 (over budget)** |
| past rung 1 (converged-only) | **52.2 %** | **0.0 %** | **0.0 %** |
| past rung 1 (attempts) | **52.2 %** | **0.0 %** | **13.8 %** |
| failed-rung iterations | **83.5 %** | **0.0 %** | **99.0 %** |
| n_material_refused | 0 | 0 | **10270** (20.4/converged step) |
| termination | WALL @ `s/B` 0.0678 | **TARGET @ 0.25000** | **BUDGET @ 0.02754** |
| q_u (kPa, NOT a capacity for any -- WALL/no-peak/BUDGET) | 1529.844 | 4387.784 | 657.342 |
| t_init (matched window, `s <= 1 mm`) | 16745.0 kPa/m (19 pts) | 16978.7 kPa/m (19 pts) | 16888.8 kPa/m (20 pts) |
| t_init (old 4-point fit, kept for continuity) | 23773.3 | 5404.2 | 17771.8 |
| tail % of t_init (matched window) | **63.3 %** | **32.6 %** | **67.4 %** |
| wall_s (see caveats above; `ctl` ran ALONE, see below) | 2467.8 s | 127.4 s | 114.8 s |

`t_init` on the matched window (`s <= max(1 mm, 5*ds_base) = 1 mm`) agrees across all
**three** arms to within **1.4 %** (16745 / 16979 / 16889 kPa/m for control / noctl /
ctl -- `ctl` sits almost exactly between the other two, despite the different build and
regime), confirming the red/blue review's finding (B4/F10) that the old 4-point fit's
4.4x gap was an artefact of an unmatched window, not a different initial problem. No
arm plateaus on the matched-window tail (`PLATEAU_FRAC = 2 %`; all three far above it).

## Gate verdicts (`adr92_bvp_gate.py`, baseline = the new `control`)

**Without `--registered-arm`** (the only leg under `--implex` is `noctl`, so this and
the `--registered-arm` run below differ only in whether that unregistered leg is kept):

```
  past rung 1               52.2 %  ->     0.0 %   (converged steps only)
  past rung 1 (attempts)    52.2 %  ->     0.0 %   (steps + nsub, apples-to-apples denominator)
  failed-rung iters         83.5 %  ->     0.0 %

  VERDICT: PREDICTION MET. The ladder is gone -- the step is effectively
  linear, which is IMPL-EX's structural claim measured on the real deck.
```

**With `--registered-arm`** (correctly drops `noctl`, since `implex_control is None`,
leaving nothing):

```
  --registered-arm: dropping 1 unregistered (implex_control=None) leg(s) from --implex: ['h1.0_e0.6944']
--registered-arm: no leg under --implex has implex_control set -- nothing left to gate
```

**Reading these together:** the "PREDICTION MET" verdict above is reachable only on
the unregistered control-OFF arm, exactly as the red/blue review found on the pre-fix
data. The pre-registered gate (`--registered-arm`, `-implex -implexControl`) now has
data and returns **PARTIAL** (see the "Registered arm" section above for the verbatim
output) -- the ladder is gone on the converged-only count (every ladder-touching
attempt failed all three rungs and was subdivided away rather than converging on rung
2/3) but still fires on 13.8 % of attempts and, when it does, burns through all three
rungs at a 99.0 % failed-iteration share, the highest of any arm in this campaign.
Builds matched for the control-OFF run (`implex ['2473ce46c']`, baseline
`['2473ce46c']`) -- no WARN there; the registered-arm run mixes builds deliberately
(see the callout above) and its WARN is expected.

## Overlay: `q_foot` at matched `s/B`, `noctl` vs `control`

| s/B | q_control (kPa) | q_noctl (kPa) | dev % | wall_control (s) | wall_noctl (s) | speedup |
|---|---|---|---|---|---|---|
| 0.00001 | 1.336 | 2.739 | **+104.98** | 0.228 | 0.058 | 3.9x |
| 0.0005 | 17.753 | 18.455 | +3.96 | 5.537 | 1.658 | 3.3x |
| 0.001 | 31.877 | 32.952 | +3.37 | 15.361 | 2.479 | 6.2x |
| 0.002 | 58.158 | 60.006 | +3.18 | 29.891 | 3.745 | 8.0x |
| 0.005 | 131.840 | 135.624 | +2.87 | 108.248 | 6.542 | 16.6x |
| 0.01 | 248.853 | 255.437 | +2.65 | 288.707 | 10.334 | 27.9x |
| 0.02 | 477.979 | 487.755 | +2.05 | 623.415 | 16.999 | 36.7x |
| 0.04 | 926.469 | 935.162 | +0.94 | 1368.939 | 28.908 | 47.4x |
| 0.055 | 1255.441 | 1271.325 | +1.27 | 1987.890 | 38.115 | 52.2x |
| 0.0678 (control's own endpoint) | 1529.844 | 1553.841 | +1.57 | 2467.774 | 45.602 | 54.1x |

Per-step (not interpolated onto this coarse grid) over the full overlap `s/B <=
0.0678` (n = 69 matched points): **mean |deviation| = 4.78 %, max |deviation| =
104.98 % at s/B = 0.00001** (step 1, the pure elastic predictor with
`d_eps_p(n) = 0`; the deviation falls under 4 % by the second checkpoint and stays
there). This is a materially *smaller* deviation than the pre-fix review's
`implex_noctl` vs `implex_OFF` comparison (11.36 % mean / +38.1 % max at s/B 0.02) --
the deeper control reach here (0.0678 vs the earlier 0.0553) extends the checkable
overlap and the deviation is monotonically settling rather than spiking mid-range.

## Overlay: `q_foot` at matched `s/B`, `ctl` vs `control`

`ctl` ran to `s/B = 0.02754` before BUDGET; the overlap with `control` (which reaches
0.0678) is bounded by `ctl`'s own endpoint. **`ctl` ran ALONE** (no other leg's process
was active during its window, 12:41-12:42) -- only the PID 74216 machine-load caveat
applies to its `wall_s`, not the sequencing defect above.

| s/B | q_control (kPa) | q_ctl (kPa) | dev % | wall_control (s) | wall_ctl (s) | speedup |
|---|---|---|---|---|---|---|
| 0.00001 | 1.336 | 2.739 | **+104.98** | 0.228 | 0.074 | 3.1x |
| 0.0005 | 17.753 | 18.435 | +3.84 | 5.537 | 2.917 | 1.9x |
| 0.001 | 31.877 | 32.897 | +3.20 | 15.361 | 5.656 | 2.7x |
| 0.002 | 58.158 | 59.753 | +2.74 | 29.891 | 11.531 | 2.6x |
| 0.005 | 131.840 | 134.716 | +2.18 | 108.248 | 27.835 | 3.9x |
| 0.01 | 248.853 | 253.579 | +1.90 | 288.707 | 53.890 | 5.4x |
| 0.02 | 477.979 | 485.178 | +1.51 | 623.415 | 90.544 | 6.9x |
| 0.0275 (`ctl`'s own endpoint) | 647.430 | 656.433 | +1.39 | 880.690 | 114.406 | 7.7x |

Per-step over the full overlap `s/B <= 0.02754` (n = 504 matched points): **mean
|deviation| = 2.30 %, max |deviation| = 104.98 % at s/B = 0.00001** (step 1, the same
pure elastic-predictor outlier seen on `noctl`). **Excluding step 1** (n = 503): mean
|deviation| = 2.10 %, max |deviation| = 21.53 %. `ctl` tracks `control` noticeably
*tighter* than `noctl` does over the same range (2.30 % mean here vs `noctl`'s 4.78 %
mean over its own, larger overlap; both dominated by the same step-1 outlier) -- the
error control is doing exactly what it is built to do, at the cost of only reaching
`s/B = 0.0275` in the process.

## `implexError` (ADR §8's own reporting mandate)

`noctl`'s curve CSV (`implex_err_avg` column, read once per step from the process-wide
`avgImplexError` accumulator, non-destructive): **mean 0.00113, max 0.00349** over the
142 converged steps -- comfortably under the registered tolerance of 0.05, which is
exactly why the control-OFF arm never refuses anything (`n_material_refused = 0`).
`ctl`'s own curve (same column, over its 504 converged steps): **mean 0.00003, max
0.00226** -- an order of magnitude *smaller* than `noctl`'s, because every step whose
error would have pushed the running average up was instead refused and retried at a
smaller `ds` before it could commit; the accumulator only ever sees the error of steps
the control let through.

`implex_err_max_at_checkpoints` (max over yielding elements of `implexDetail[0]`, the
per-Gauss-point total error slot):

| s/B | control (implex off) | noctl | ctl |
|---|---|---|---|
| 0.002 | 0.0 | 0.00326 | 0.00044 |
| 0.005 | 0.0 | 0.00959 | 0.00032 |
| 0.01 | 0.0 | 0.01640 | 0.00072 |
| 0.02 | 0.0 | 0.04256 | 0.00036 |
| 0.04 | 0.0 | 0.04878 | -- (`ctl` did not reach) |
| 0.08 | -- (control did not reach) | 0.01483 | -- |
| 0.15 | -- (control did not reach) | 0.01688 | -- |

Note the **max**-over-elements figure on `noctl` (0.0426-0.0488 through s/B 0.02-0.04)
sits close to the registered tolerance of 0.05, while `ctl`'s checkpoint max stays two
orders below it (3.2e-4 to 7.2e-4) throughout -- consistent with ADR §8's own
prediction that the error concentrates at a small number of corner Gauss points: on
`ctl` the control is catching and refusing exactly those points *before* a checkpoint
snapshot can see them elevated, which is also why `n_material_refused = 10270` is so
much larger than the 81 abandoned attempts alone -- most refusals are single ladder
rungs within an otherwise-successful step's retries, not whole-step losses.

## One paragraph per arm

**control.** Terminated `WALL` at `s/B = 0.0678` after 69 converged steps, 2467.8 s
wall (contended, see caveats): 33 on rung 1, 15 on rung 2, 21 on the relaxed rung 3, 0
subdivisions, `n_material_refused = 0` (implex off, not applicable). `past rung 1` =
52.2 % on both the converged-only and attempts-based denominators (identical here
because `nsub = 0`). `q_u = 1529.844` kPa is not a capacity (WALL termination, no
peak). This is a materially *deeper* reach than the pre-fix `implex_OFF` control
(`s/B = 0.0553` in 2402 s) on the same 2400 s wall budget -- attributable to per-step
cost, not the budget, since this run also spent time in a genuine three-way CPU
contention window and still went further; consistent with the earlier finding that the
control leg's step cost, not iteration count, is what the wall budget is actually
buying.

**noctl.** Terminated `TARGET` at `s/B = 0.25` after 142 converged steps, all on rung 1
(0 rung-2/3, 0 subdivisions), 127.4 s wall (also contended for its first ~128 s, the
entire span it ran). `implexError` stayed small throughout (`avgImplexError` mean
0.00113, max 0.00349; per-checkpoint element-max 0.0033-0.0488) even though nothing
polices it on this arm -- consistent with the earlier finding that the deck's
free-surface corner, not the bulk of the mesh, is where the extrapolation error
concentrates, and this arm simply never triggers a check there. `t_init` on the
matched window (16978.7 kPa/m) agrees with control's (16745.0) to 1.4 %; the old
4-point fit's headline 4.4x gap does not reproduce with a like-for-like window. `tail
= 32.6 %` of the matched-window `t_init` -- still no plateau, at any denominator.

**ctl (registered arm).** After lane A's un-primed-step exemption (build
`afb95c40c9`), terminated **BUDGET** at `s/B = 0.02754` after 504 converged steps
(all rung 1) and 81 abandoned 3-rung subdivisions (`nsub = 81` against
`subdiv_budget = 80` -- one over), `n_material_refused = 10270` (20.4 per converged
step), 114.8 s wall (single-tenant, only the PID 74216 machine-load caveat applies).
Gate `--registered-arm` verdict: **PARTIAL** (0.0 % past rung 1 converged-only,
13.8 % on attempts, 99.0 % failed-rung iterations -- the highest failed-rung share
measured in this campaign). Tracks `control` *tighter* than `noctl` does over their
respective overlaps (2.30 % mean deviation vs 4.78 %, both dominated by the same
step-1 elastic-predictor outlier), and its own `implexError` stays two orders below
the registered tolerance at every checkpoint it reached -- the control is doing its
job, at the cost of reaching only 41 % of `control`'s own depth (`s/B` 0.0275 vs
0.0678) in comparable wall time (114.8 s vs 2467.8 s, though `control`'s time is
contended and not a clean comparator either). `ctl_diag/` (committed alongside) is lane A's own pre-fix reproduction of the
zero-step failure, kept for provenance.

## What is NOT yet answered

- Whether `reductionLimit = 0.01` is simply too tight for this deck past the
  un-primed first step (the ladder still fires on 13.8 % of attempts and burns every
  rung when it does), or whether a looser value (the memo's owed tolerance sweep,
  still not run) reaches depth while keeping `implexError` bounded -- `ctl`'s own
  checkpoint-max `implexError` (3.2e-4 to 7.2e-4) is two orders below the 0.05
  tolerance, which argues the refusals are not marginal near-misses but a stricter
  population the control is catching well inside its budget.
- A clean (single-tenant) wall-clock measurement for `control` and `noctl` -- both the
  PID 74216 background load and the control/noctl overlap window (disclosed above)
  contaminate their timing columns; `ctl`'s own wall_s is single-tenant modulo PID
  74216 only.
- Whether `ctl` would reach `control`-comparable depth given a larger `subdiv_budget`
  (it hit `BUDGET`, not `FLOOR` or a converged termination) -- not measured here.

## Log

- 2026-09-06 -- Rerun on pinned build `2473ce46c2ae44412246ca305f3843463f389cad`
  (lane A's pytest battery: 20 passed, immediately prior). `control` and `noctl`
  legs complete and gated. Sequencing defect (control/noctl overlap, ~158 s) and a
  pre-existing unrelated-process machine load disclosed above.
- 2026-09-06 (later) -- Registered arm (`ctl`) filled in on build
  `afb95c40c9e3da1cc97d6aabe8168644c53d0e82`, after lane A's un-primed-step exemption
  fix (the first push step after the stage flip has no plastic history, so
  `implexError` there measured the companion's drift-correction jump, not a genuine
  extrapolation error; that step is now exempt from refusal). `control`/`noctl` were
  **not** re-run -- the fix cannot affect either (control-OFF has no companion to
  refuse from; the un-primed-step exemption only changes `-implexControl` behaviour,
  which `noctl` never exercises since its control is off). Gate `--registered-arm`
  verdict: PARTIAL. All three arms now complete.

## Operating-point sweep (`-implexControl` tolerance / reductionLimit)

Snapshot binary: `dist\bin` of this worktree copied to a scratch directory
**before** another lane began rebuilding it for the mutation gate (`LADRUNO_DIST_BIN`
override, verified via `ops.ladrunoBuild()` before the sweep and re-confirmed in every
leg's `run.log`: all four report `afb95c40c9e3da1cc97d6aabe8168644c53d0e82`). Same
deck as the registered arm (`h1.0_e0.6944`, `--surcharge 10`, `--maxsubsteps 20000`,
`--wall 2400`), four legs, sequential, into
`Ladruno_files/testbed/hypo_bearing/adr92_bvp_fix/sweep/<name>/`: `tol0.1`, `tol0.2`,
`tol0.5` (`reductionLimit` held at the registered `0.01`), and `tol0.05_rl0.1`
(`tol` held at the registered `0.05`, `reductionLimit` loosened 10x to `0.1`).

| leg | tol / rLimit | mode | s/B (depth) | steps | nsub | n_refused | refused/step | wall_s | dev mean % (excl. step1) | dev max % (excl. step1) | implexErr avg / max |
|---|---|---|---|---|---|---|---|---|---|---|---|
| control (ref.) | -- | WALL | **0.0678** | 69 | 0 | 0 | -- | 2467.8 | -- | -- | -- |
| ctl (registered, ref.) | 0.05 / 0.01 | BUDGET | 0.02754 | 504 | 81 | 10270 | 20.38 | 114.8 | 2.10 | 21.53 | 0.00003 / 0.00226 |
| `tol0.05_rl0.1` | 0.05 / **0.1** | BUDGET | 0.02754 | 504 | 81 | 10270 | 20.38 | 117.1 | 2.10 | 21.53 | 0.00003 / 0.00226 |
| `tol0.1` | **0.1** / 0.01 | BUDGET | **0.07642** | 514 | 81 | 9299 | 18.09 | 180.0 | 1.87 | 21.53 | 0.00005 / 0.00226 |
| `tol0.2` | **0.2** / 0.01 | BUDGET | **0.15034** | 515 | 81 | 1494 | 2.90 | 276.2 | 1.91 | 21.53 | 0.00007 / 0.00226 |
| `tol0.5` | **0.5** / 0.01 | **TARGET** | **0.25000** | 419 | 63 | 750 | 1.79 | 409.4 | 2.29 | 21.53 | 0.00018 / 0.00226 |

Two structural findings before the answer:

1. **`reductionLimit` is inert at `tol = 0.05` on this deck.** `tol0.05_rl0.1` is
   byte-identical to the registered `ctl` leg in every field above (same `s/B`,
   `steps`, `nsub`, `n_material_refused`, `q_u`) despite a 10x looser floor -- the
   `tol` criterion (`implexError > 0.05`) is what refuses every step here; the
   `reductionLimit` escape at `:1615` never gets a chance to differ because the ladder
   never subdivides `ds` far enough for the two floors (`0.01 x ds0` vs `0.1 x ds0`) to
   diverge before the `tol` refusal already fires. Loosening `reductionLimit` alone,
   without loosening `tol`, buys nothing on this deck.
2. **`max |dev| = 21.53 %` (excluding step 1) is identical across every `-implex` leg
   in this campaign** -- `noctl`, the registered `ctl`, and all four sweep legs. It
   occurs at the same `s/B` in the shared low-settlement overlap and is not moved by
   loosening the tolerance, which means it is a feature of the early-step trajectory
   itself (plausibly the same step-2 `q`-drop noted for `noctl` in the earlier
   section), not something an `-implexControl` operating point can tune away.

**Answer.** `tol = 0.1` (`reductionLimit` at the registered `0.01`) is the tightest
tolerance tested that satisfies **both** clauses: it reaches `s/B = 0.07642`, past
`control`'s own depth of `0.0678`, while keeping the mean overlay deviation (excluding
the shared step-1 elastic-predictor outlier) at `1.87 %`, comfortably under `5 %`.
`tol = 0.2` and `tol = 0.5` also satisfy both (and `tol = 0.5` reaches `TARGET`
outright), each with mean deviation still under `2.3 %` -- accuracy does not degrade
monotonically with a looser tolerance on this deck, at least up to `0.5`. The
registered operating point (`tol = 0.05`), with **either** `reductionLimit`, fails
clause (a) outright: it never reaches `control`'s depth (`0.0275` vs `0.0678`)
regardless of its accuracy where it does reach, because the substep budget (not the
reduction floor) is what stops it. **If a single answer is wanted: `0.05` is too tight
for this deck at the registered substep cap; `0.1` is the smallest loosening that
clears both bars, and there is currently no evidence in this sweep that going further
(`0.2`, `0.5`) costs anything in accuracy.**

### Log addendum

- 2026-09-06 (later still) -- Operating-point sweep run on a **snapshot** of build
  `afb95c40c9e3da1cc97d6aabe8168644c53d0e82` (copied to a scratch `dist\bin` before
  another lane's mutation-gate rebuild could overwrite the worktree's own binary;
  `ops.ladrunoBuild()` verified against the snapshot before the sweep and re-confirmed
  in all four `run.log`s). `control`/`noctl`/`ctl` unaffected, not re-run.
