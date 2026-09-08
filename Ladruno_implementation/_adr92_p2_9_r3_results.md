---
title: "ADR 92 / P2-9 — control-informed factor: Fork R3 registered-arm results (REFUTED)"
project: Ladruno
type: results
status: "Leg 1 (`control`, first-iterate f*) REFUTED 2026-09-08; Leg 2 (`controlIter`, per-trial f* recompute, build 9e73060a9) PASSES the depth/overlay bars 2026-09-08"
priority: high
owner: nmora
related:
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[_adr92_p2_9_control_informed_f_plan]]"
  - "[[_adr92_p1_bvp_gate_rerun]]"
tags: [adr, sanisand, implex, p2-9, measurement, refuted]
updated: 2026-09-08
---

# ADR 92 / P2-9 — Fork R3 registered arm: results (REFUTED)

## Provenance

- Build: `a6a53948ef316f7eb6662ea09f2e1fb52d2b607f` (HEAD of `wp/92f-implex-control-f`;
  `python3.12 -c "import opensees as o; print(o.ladrunoBuild())"` verified before the run).
- Driver: `Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py`, uncommitted edit
  that adds `--implex-factor {fixed,control}` (appends `-implexFactor control` to the
  `nDMaterial LadrunoSANISAND` call, and records `implexGuards[6]` as `guard_ctlf` /
  `n_guard_ctlf` per-step in the curve CSV and the leg JSON). `git diff --stat` on the
  driver: 32 insertions / 10 deletions, one file.
- Invocation (same deck as the P1/P2 `tol0.1` legs, plus the new flag), run from the
  worktree root with `PYTHONPATH`/`PATH` pointed at `dist\bin`,
  `LADRUNO_OPENSEES_QUIET=1`, and `LADRUNO_A2_EXPECT_BUILD` set to the HEAD hash above
  (the driver's `assert_engine()` pins a hardcoded `EXPECTED_BUILD` that predates this
  worktree's rebuild, so the override was necessary and deliberate):

  ```
  python3.12 -u Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py \
    --out Ladruno_files/testbed/hypo_bearing/adr92_bvp_fix/p2_9_ctl \
    --legs h1.0_e0.6944 \
    --implex --implex-control 0.1 0.01 --implex-factor control \
    --surcharge 10 --maxsubsteps 20000 --wall 2400
  ```

  Output: `Ladruno_files/testbed/hypo_bearing/adr92_bvp_fix/p2_9_ctl/`. Actual wall time
  141.6 s (well inside the 2400 s budget — the leg ended by seizing, not by timing out).

## Comparison table

| | P1 registered `tol0.1` (`afb95c40c`) | P2 `tol0.1` (`87b9cf846`) | **P2-9 `control` (`a6a53948e`, this run)** |
|---|---|---|---|
| mode | BUDGET | BUDGET | **FLOOR (SEIZED)** |
| depth reached, `s/B` | 0.07642 | 0.07594 | **0.05210** |
| converged steps | 514 | 514 | **366** |
| attempts (steps+nsub) | -- | 595 | 429 |
| nsub | 81 | 81 | 63 |
| `n_material_refused` | 9299 | 9235 | **845** |
| refused / converged step | 18.09 | 17.97 | **2.31** |
| `n_guard_floor` (P2-1) | -- | 0 | 0 |
| `n_guard_f0` (P2-2, committed-predecessor `f=0`) | -- | 2913 (5.67/step) | 3722 (10.17/step) |
| `n_guard_hold` (P2-3) | -- | 0 | 0 |
| `n_guard_res` (P2-5c, hold-skip commits, `+N`/hold) | -- | 0 | **1,196,916** |
| `n_guard_ctlf` (P2-9, `f*<0.5·f_max` backoff) | n/a (predates P2-9) | n/a | **196,549** |
| overlay vs `control/`, mean \|dev\| (excl. step 1) | 1.87 % | 1.92 % | **11.12 %** |
| overlay vs `control/`, max \|dev\| (excl. step 1) | 21.53 % | 23.99 % (at `s/B=0.00002`) | **24.89 % (at `s/B=0.00002`)** |
| wall_s | 180.0 | 192 | 141.6 |

Overlay methodology (unchanged from the P1/P2 rerun memo): `control/`'s 69-row curve
(ending `s/B=0.0678`) is linearly interpolated onto the finer leg's own `s/B` grid over
the overlap window (`control`'s endpoint bounds it; all 366 of this leg's rows fall
inside it since the leg itself only reached `s/B=0.0521`), `dev % = (q_arm -
q_control_interp)/q_control_interp*100` at each matched point, and step 1
(`s/B=1e-5`, the shared pure-elastic-predictor outlier at +104.98 % on every `-implex`
arm in this campaign) is excluded from the mean/max headline.

Coarse checkpoint-grid overlay (this run vs `control/`):

| s/B | q_control (kPa) | q_p2_9_ctl (kPa) | dev % |
|---|---|---|---|
| 0.0005 | 17.753 | 20.254 | +14.09 |
| 0.001 | 31.336 | 34.655 | +10.59 |
| 0.002 | 58.669 | 63.017 | +7.41 |
| 0.005 | 131.361 | 138.117 | +5.14 |
| 0.01 | 249.315 | 273.137 | +9.55 |
| 0.02 | 476.637 | 535.747 | +12.40 |
| 0.04 | 925.149 | 1051.431 | +13.65 |

The R3 gate (`adr92_bvp_gate.py --implex adr92_bvp_fix/p2_9_ctl --baseline
adr92_bvp_fix/control --registered-arm`, build-mismatch WARN expected since baseline
predates HEAD): past-rung-1 (converged-only) 52.2 % -> 1.1 %, past-rung-1 (attempts)
52.2 % -> 15.6 %, failed-rung iters 83.5 % -> 92.1 %. Own verdict string: PARTIAL — not
independently decisive here; the depth/overlay bars below are.

## Verdict: REFUTED

Per the plan's rule ("Fork R3 registered arm ... refutes P2-9 if depth < 0.076 or
overlay > 5 %"), **both** clauses trip:

- **Depth**: 0.05210 < 0.076 — the leg seized about 31 % short of P2's own depth, and
  in a different, worse failure mode (`FLOOR`: every ladder rung failed at the
  `DS_MIN` floor) rather than the `BUDGET` wall-clock timeout the `tol0.1` legs hit
  after reaching their target.
- **Overlay**: 11.12 % mean \|dev\| >> 5 % (and >> the plan's 2 % target), roughly 6x
  the `tol0.1` legs' ~1.9 %.

The one part of the pre-registered prediction that *did* land is "refusals/step down":
`n_material_refused`/converged-step fell from 17.97 (P2) to 2.31 — an ~8x drop,
consistent with the closed-form factor extrapolating through cases the old `f=0`
guard used to refuse outright. But that reduction in explicit refusals bought a worse
outcome on every load-path metric that matters: fewer refusals, more of the P2-2 guard
firing (`n_guard_f0` per-step nearly doubled, 5.67 -> 10.17), a new and very large
`n_guard_res` (hold-skip-commit) count (0 -> 1.196M) that appears wherever the graded
factor lets a near-zero-`dt` commit through repeatedly, a P2-9 backoff census
(`n_guard_ctlf`) firing 196,549 times over only 366 steps, and a curve that visibly
diverges from the non-implex baseline by double digits well before the seizure (e.g.
+12-14 % by `s/B=0.02-0.04`, vs ~2 % for `tol0.1`). Read together, the control-informed
`f*` factor is letting through extrapolation that the fixed-threshold `f=0` guard used
to block, and on this deck that extra extrapolation both drifts the curve away from
the reference and precipitates an earlier, harder numerical seizure — the opposite of
what P2-9 predicted (`depth >= 0.076, overlay <= 2 %`). Per the plan's own decision
rule, P2-8's fixed threshold remains the guard to keep; the graded factor should be
recorded as measured, not shipped as the new default.

## Caveats

- A pytest battery may have shared this machine's cores during the run (a
  `tasklist`-visible pool of other `python3.12.exe` processes was present partway
  through); wall-clock numbers above (141.6 s / 192 s / 180 s) are indicative only,
  not a clean single-tenant measurement — same caveat class as the earlier P1/P2
  rerun memos in this campaign. The depth/overlay/refusal/guard-census numbers are
  step-indexed and unaffected by contention.
- Baseline (`control/`) and this leg ran on different engine builds
  (`2473ce46c` vs `a6a53948e`), same as every prior arm in this campaign; the gate
  script flags this itself.
- Only the `control` factor arm was run here (the WP's ask); `fixed` byte-identity
  under the new flag was not re-verified in this session (it is covered by the
  C++-level unit tests referenced in the P2-9 plan, lane "Tests").

## Leg 2 -- `controlIter` (build 9e73060a9)

### Provenance

- Build: `9e73060a90e81a35ebd842bcd10a315a980ffac4` (restamp of `wp/92f-implex-control-f`
  adding the third `-implexFactor controlIter` mode: f* is recomputed at EVERY trial of
  the step from that iterate's `d_eps`, not frozen from the first iterate; `f_max` is
  still stored per step and the back-off census (`n_guard_ctlf`) still fires once per
  step). Verified via `ladrunoBuild()` before the run.
- Driver: `Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py` extended so
  `--implex-factor` now takes `{fixed, control, controlIter}`; the new token is passed
  through verbatim as the `-implexFactor` value (the `fixed`/`control` code paths are
  unchanged -- only the tuple of accepted tokens and the flag-emission condition grew a
  branch).
- Invocation (identical to Leg 1's, `--implex-factor controlIter` substituted, plus the
  new build's expected-hash override):

  ```
  python3.12 -u Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py \
    --out Ladruno_files/testbed/hypo_bearing/adr92_bvp_fix/p2_9_iter \
    --legs h1.0_e0.6944 \
    --implex --implex-control 0.1 0.01 --implex-factor controlIter \
    --surcharge 10 --maxsubsteps 20000 --wall 2400
  ```
  (`LADRUNO_A2_EXPECT_BUILD=9e73060a90e81a35ebd842bcd10a315a980ffac4` set to satisfy the
  driver's hardcoded `assert_engine()` pin, same reason as Leg 1.)

  Output: `Ladruno_files/testbed/hypo_bearing/adr92_bvp_fix/p2_9_iter/`. Actual wall time
  1858.4 s (well inside the 2400 s budget -- again ended by seizing, not by timing out).

### Comparison table (all four arms)

| | P1 `tol0.1` (`afb95c40c`) | P2 `tol0.1` (`87b9cf846`) | Leg 1 `control` (`a6a53948e`) | **Leg 2 `controlIter` (`9e73060a9`, this run)** |
|---|---|---|---|---|
| mode | BUDGET | BUDGET | FLOOR (SEIZED) | **FLOOR (SEIZED)** |
| depth reached, `s/B` | 0.07642 | 0.07594 | 0.05210 | **0.11486** |
| converged steps | 514 | 514 | 366 | **327** |
| attempts (steps+nsub) | -- | 595 | 429 | **384** |
| nsub | 81 | 81 | 63 | **57** |
| `n_material_refused` | 9299 | 9235 | 845 | **2418** |
| refused / converged step | 18.09 | 17.97 | 2.31 | **7.39** |
| `n_guard_floor` (P2-1) | -- | 0 | 0 | **0** |
| `n_guard_f0` (P2-2) | -- | 2913 (5.67/step) | 3722 (10.17/step) | **5146 (15.74/step)** |
| `n_guard_hold` (P2-3) | -- | 0 | 0 | **0** |
| `n_guard_res` (P2-5c) | -- | 0 | 1,196,916 | **1,626,519 (4974/step)** |
| `n_guard_ctlf` (P2-9 backoff) | n/a | n/a | 196,549 | **594,583 (1818/step)** |
| overlay vs `control/`, mean \|dev\| (excl. step 1) | 1.87 % | 1.92 % | 11.12 % | **1.70 %** |
| overlay vs `control/`, max \|dev\| (excl. step 1) | 21.53 % | 23.99 % | 24.89 % (at `s/B=0.00002`) | **24.89 % (at `s/B=0.00002`, same shared elastic-predictor artifact)** |
| wall_s | 180.0 | 192 | 141.6 | **1858.4** |
| engine-log max-iter stalls (`after: 25 iterations`) | -- | -- | 1 (over 366 steps) | **89 (over 327 steps)** |

Overlay recomputed with the same methodology as Leg 1 (linear interpolation of
`control/`'s 69-row curve onto the leg's own `s/B` grid over the overlap window,
`s/B=1e-5` step-1 excluded, `dev % = (q_arm - q_control_interp)/q_control_interp*100`);
the Leg-1 numbers above were independently reproduced with the same script (11.12 % /
24.89 % at `s/B=0.00002`) before trusting the Leg-2 numbers.

Coarse checkpoint-grid overlay, Leg 2 vs `control/`:

| s/B | q_control (kPa) | q_controlIter (kPa) | dev % |
|---|---|---|---|
| 0.0005 | 17.753 | 18.448 | +3.92 |
| 0.001 | 31.877 | 32.559 | +2.14 |
| 0.002 | 58.158 | 59.294 | +1.95 |
| 0.005 | 131.840 | 133.670 | +1.39 |
| 0.01 | 248.853 | 251.584 | +1.10 |
| 0.02 | 477.979 | 483.932 | +1.25 |
| 0.04 | 926.469 | 936.801 | +1.12 |

(Small numeric differences from Leg 1's own checkpoint table above are the two
`control/`-side interpolations picking slightly different bracketing rows off the same
69-row baseline curve at each arm's own `s/B` grid -- not a baseline change.)

Gate (`adr92_bvp_gate.py --implex adr92_bvp_fix/p2_9_iter --baseline
adr92_bvp_fix/control --registered-arm`, same build-mismatch WARN as Leg 1): past-rung-1
(converged-only) 52.2 % -> 27.2 %, past-rung-1 (attempts) 52.2 % -> 38.0 %, failed-rung
iters 83.5 % -> 97.6 %. Gate's own verdict string: PARTIAL -- again not independently
decisive; the depth/overlay bars below are.

### Verdict: PASS

Per the plan's rule ("depth >= 0.076 and overlay <= 2 % PASS; depth < 0.076 or overlay >
5 % REFUTED; between = partial"), Leg 2 clears **both** clauses:

- **Depth**: 0.11486 >= 0.076 -- the leg reached 51 % further than the depth bar, and
  in fact 2.2x deeper than Leg 1's `control` arm (0.05210) before its own FLOOR seizure.
- **Overlay**: 1.70 % mean \|dev\| <= 2 % -- inside the plan's target band, and better
  than every other arm measured in this campaign including the non-implex-informed
  `tol0.1` legs (1.87 %/1.92 %).

### Leg 1 vs Leg 2 -- did per-iterate recompute remove the bias?

Yes, on the metric that mattered: recomputing f* at every Newton trial (instead of
freezing it from the step's first, elastic-predictor-biased iterate) removed almost all
of the curve drift Leg 1 showed -- mean overlay deviation fell from 11.12 % to 1.70 %,
and the +12-14 % drift Leg 1 exhibited by `s/B=0.02-0.04` collapses to +1.1-1.3 % at the
same depths in Leg 2's checkpoint table. The diagnosis behind Leg 2 (elastic-predictor
bias in the first-iterate f*) is confirmed by this outcome: fixing the bias source fixed
the curve.

That said, per-trial recompute is not free. Three costs show up:

1. **Newton-iteration cost.** The engine log's max-iteration stall marker
   (`WARNING: CTestNormUnbalance::test() - failed to converge ... after: 25
   iterations`) fires 89 times in Leg 2's 327-step run vs only once in Leg 1's 366-step
   run -- roughly two orders of magnitude more non-convergent retries per step. Neither
   leg's JSON/curve-CSV records a clean "Newton iterations per converged step" field, so
   this stall count is the best available proxy; it points at controlIter buying its
   overlay accuracy with materially more solver churn per step, not fewer refusals.
2. **More material refusals, not fewer.** `n_material_refused`/converged-step rose from
   2.31 (Leg 1) to 7.39 (Leg 2) -- the opposite direction from Leg 1's own headline
   result (18 -> 2.3 refusals/step going from `tol0.1`/`fixed` to `control`). Recomputing
   f* every trial re-exposes more trials to the `-implexControl` tolerance check that a
   frozen-for-the-step f* was implicitly shielding.
3. **Guard census still large.** `n_guard_res` (hold-skip commits) and `n_guard_ctlf`
   (P2-9 backoff) both grew per-step relative to Leg 1 (4974 vs 3270/step, 1818 vs
   537/step respectively) -- the graded factor is still triggering its own back-off
   machinery constantly, just from a less-biased anchor point.

Net: `controlIter` is the first `-implexFactor` mode in this campaign to independently
clear both the depth and overlay bars, at the cost of running ~13x longer in wall time
(1858 s vs 141.6 s) than Leg 1 for a shallower budget headroom, and with visibly more
per-step solver retries. Whether that trade is worth shipping as a default depends on
whether the campaign is optimizing for curve fidelity (controlIter wins outright) or for
wall-clock/refusal-count cost (Leg 1's frozen `control` factor was cheaper per step, just
wrong).
