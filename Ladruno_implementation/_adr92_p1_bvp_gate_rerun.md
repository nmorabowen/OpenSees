---
title: "ADR 92 / P1 gate 1 rerun — control and uncontrolled arms on the fixed binary"
project: Ladruno
type: results
status: "Registered arm (-implex -implexControl 0.05 0.01): NO VERDICT -- gate refuses verbatim `no leg records under adr92_bvp_fix/ctl` (0 converged steps, 0 JSON, before lane A's re-fix; PENDING RE-RUN). Control and -implex/control-OFF arms: COMPLETE, gated, VERDICT (control-OFF arm only) = PREDICTION MET."
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

## Registered arm (`ctl`, `-implex -implexControl 0.05 0.01`) -- PLACEHOLDER, PENDING

**This section is intentionally incomplete.** The registered arm produced **zero
converged steps and zero JSON leg records** (the driver's `assert len(rows) > 8`
fired before any checkpoint or final write). `adr92_bvp_gate.py --implex
adr92_bvp_fix/ctl --baseline adr92_bvp_fix/control` refuses verbatim:

```
no leg records under adr92_bvp_fix/ctl
```

This is a genuine regression in *reach* relative to the pre-fix `implex_ctl` leg (which
got 110 converged steps to `s/B = 0.005` before seizing): the B2 floor fix landed by
lane A now arms `mImplexDt0` correctly and enforces the reduction limit as designed,
and at `reductionLimit = 0.01` that enforcement refuses **every** ladder rung on the
very first attempted step (`implexError` 0.129-0.295 vs `tol 0.05`, `|dt| = 2e-05`
already at the `floor 2e-07`) -- the material's own throttled log records 10 explicit
`-implexControl REFUSES` lines (budget-suppressed after that; `implexRefusals` was
never read since the process exited via the driver's own control-assert before any
`eleResponse` call) and `7 NewtonRaphson + 7 NewtonLineSearch + 7 AcceleratedNewton`
ladder-rung failures (7 subdivisions x 3 rungs = 21, consistent with the FLOOR
termination). Lane A is re-fixing this and will re-run into `adr92_bvp_fix/ctl/`
separately; **do not read the control/noctl comparison below as a "with vs without
control" result** -- it is control vs the unregistered control-OFF arm only, exactly
as the red/blue review found for the pre-fix data.

## Overview

| | control | `-implex`, control OFF (`noctl`) |
|---|---|---|
| steps (converged) | 69 | **142** |
| attempts (steps + nsub) | 69 | 142 |
| rung 1 / rung 2 / rung 3 | 33 / 15 / 21 | **142 / 0 / 0** |
| nsub (subdivisions) | 0/80 | 0/80 |
| past rung 1 (converged-only) | **52.2 %** | **0.0 %** |
| past rung 1 (attempts) | **52.2 %** | **0.0 %** |
| failed-rung iterations | **83.5 %** | **0.0 %** |
| n_material_refused | 0 | 0 |
| termination | WALL @ `s/B` 0.0678 | **TARGET @ 0.25000** |
| q_u (kPa, NOT a capacity for either -- WALL/no-peak) | 1529.844 | 4387.784 |
| t_init (matched window, `s <= 1 mm`, 19/19 pts) | 16745.0 kPa/m | 16978.7 kPa/m |
| t_init (old 4-point fit, kept for continuity) | 23773.3 | 5404.2 |
| tail % of t_init (matched window) | **63.3 %** | **32.6 %** |
| wall_s (see caveats above) | 2467.8 s | 127.4 s |

`t_init` on the matched window (`s <= max(1 mm, 5*ds_base) = 1 mm`) agrees between the
two arms to **1.4 %** (16745 vs 16979 kPa/m), confirming the red/blue review's finding
(B4/F10) that the old 4-point fit's 4.4x gap was an artefact of an unmatched window, not
a different initial problem. Neither arm plateaus on the matched-window tail
(`PLATEAU_FRAC = 2 %`; both far above it).

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

**Reading these together:** the "PREDICTION MET" verdict above is, exactly as the
red/blue review found on the pre-fix data, reachable **only** on the unregistered
control-OFF arm. The pre-registered gate (`--registered-arm`, `-implex
-implexControl`) currently has no data to report at all -- not PARTIAL, not REFUTED,
literally empty, because the registered leg cannot complete a single step on this
binary. Builds matched on both runs (`implex ['2473ce46c']`, `baseline
['2473ce46c']`) -- no WARN.

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

## `implexError` (ADR §8's own reporting mandate)

`noctl`'s curve CSV (`implex_err_avg` column, read once per step from the process-wide
`avgImplexError` accumulator, non-destructive): **mean 0.00113, max 0.00349** over the
142 converged steps -- comfortably under the registered tolerance of 0.05, which is
exactly why the control-OFF arm never refuses anything (`n_material_refused = 0`) while
the registered arm above refuses on step 1.

`implex_err_max_at_checkpoints` (max over yielding elements of `implexDetail[0]`, the
per-Gauss-point total error slot):

| s/B | control (implex off) | noctl |
|---|---|---|
| 0.002 | 0.0 | 0.00326 |
| 0.005 | 0.0 | 0.00959 |
| 0.01 | 0.0 | 0.01640 |
| 0.02 | 0.0 | 0.04256 |
| 0.04 | 0.0 | 0.04878 |
| 0.08 | -- (control did not reach) | 0.01483 |
| 0.15 | -- (control did not reach) | 0.01688 |

Note the **max**-over-elements figure (0.0426-0.0488 through s/B 0.02-0.04) sits close
to the registered tolerance of 0.05, while the **average** above is two orders smaller
-- confirming ADR §8's own prediction that the error concentrates at a small number of
corner Gauss points rather than being uniform, which is exactly the population the
registered `-implexControl 0.05` arm is built to catch, and exactly why it refuses
immediately once armed correctly.

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

**ctl (registered arm).** See the placeholder section above. Zero converged steps,
zero leg records, gate refuses with no data to average. PENDING a lane-A re-run once
the immediate first-step refusal (an `implexError` of 0.13-0.30 against `tol 0.05` at
the deck's very first, smallest step) is understood -- whether that is a correctly
strict floor doing its job on a genuinely too-tight `reductionLimit`, or a residual
defect, is not decided by this data and must not be guessed at here.

## What is NOT yet answered

- The registered gate verdict (PASS/PARTIAL/REFUTE on `-implex -implexControl`) --
  no data.
- Whether `reductionLimit = 0.01` is simply too tight for this deck's first step, or
  whether a looser value (the memo's owed tolerance sweep, still not run) reaches
  depth while keeping `implexError` bounded.
- A clean (single-tenant) wall-clock measurement -- both the PID 74216 background load
  and the control/noctl overlap window contaminate every number in this memo's timing
  columns.

## Log

- 2026-09-06 -- Rerun on pinned build `2473ce46c2ae44412246ca305f3843463f389cad`
  (lane A's pytest battery: 20 passed, immediately prior). `control` and `noctl`
  legs complete and gated; `ctl` (the registered arm) produced zero data and is
  PENDING a lane-A re-fix and re-run. Sequencing defect (control/noctl overlap,
  ~158 s) and a pre-existing unrelated-process machine load disclosed above.
