---
title: "WP-193 — cost-aware refactorization trigger for factorization reuse (KrylovNewton -reuseCost)"
project: Ladruno
type: performance work package (solve lane)
status: "implemented, draft PR #980 — not merged; default off"
priority: high
owner: nmora
related:
  - "[[188_ladruno_factor_reuse_adr]] (the reuse mechanism this trigger steers)"
  - "[[40_ladruno_performance_adr]] (lever 1: factorization reuse)"
  - "paper repo: review/engstruct-2026/30_campaigns/timehistory/perf/validation_v3/table_v3_1D.md (the 0.66-0.80x 1D result)"
tags: [performance, solver, algorithm, newton, mpi, wp-193]
updated: 2026-10-14
---

# WP-193 — cost-aware refactorization trigger

> [!summary] The short version
> WP-188's triggers refactor on age, slow convergence, rate or failure. They cannot see
> what an iteration and a factorization *cost*, so the same `-reuseFactor m` that wins on
> 4D (a factorization ≈ 20 iterations) loses on the Paper 1 1D cell (0.66–0.80× of full
> Newton: 5.4–6.9 iterations per step on a kept factor instead of 2, and there an
> iteration costs more than the factorization it saves). `-reuseCost on` measures both
> costs while the analysis runs and keeps a factor only while that is cheaper:
>
> ```tcl
> algorithm KrylovNewton -maxDim 20 -reuseFactor 1000 -reuseCost on \
>     ?-reuseCostSafety s? ?-reuseCostWarmup w? ?-reuseCostProbe p?
> ```
>
> Off by default; off is byte-identical to WP-188 (and to vanilla without `-reuseFactor`).
> The WP-188 triggers stay as hard limits. On the test frame it picks "refactor every step"
> when an iteration costs 17× a factorization and "keep" when a factorization costs 200×
> an iteration, each time within 1.03–1.24× of the better fixed policy; on a 52 k-equation
> MUMPS block it finds WP-188's best `m = 20` result (2.84× full Newton) with no `m` to tune.

## 1. Syntax

| Option | Meaning | Default |
|---|---|---|
| `-reuseCost on\|off` | Enable the cost model. Needs `-reuseFactor m` (fail closed, like every `-reuse*` option); `m` stays the hard age cap, so `-reuseFactor 1000 -reuseCost on` lets the cost model alone decide. | off |
| `-reuseCostSafety s` | `s > 0`. Multiplies the projected cost of continuing on the kept factor, in both rules (§3). `s > 1` refactors earlier (closer to full Newton), `s < 1` keeps longer. | 1.0 |
| `-reuseCostWarmup w` | The first `w` steps follow the WP-188 age rule only, while the costs are learned. | 4 |
| `-reuseCostProbe p` | After `p` consecutive step starts at which the cost rule refused the kept factor, keep it once, so a stale "reuse does not pay" estimate is re-measured (the tangent and the costs change during a record). `0` = never. | 20 |

The tuning options without `-reuseCost on` refuse the command. Statistics: the WP-188 line,
unchanged, plus one line (rank 0 at the end / `-reuseReport`; every rank with `print -algorithm`).
From the brick block of §6.4:

```
KrylovNewton reuse cost model (WP-193) so far [rank 0]: t_factor 1.5502 s, t_iter 0.0892805 s,
ratio 17.3633, n_fresh 3, n_kept1 7.80851, slope 0.0529379, starts kept by cost 35,
starts formed by cost 1, in-step cost refreshes 0, probes 0
```

## 2. What is measured

All wall clock (`std::chrono::steady_clock`), inside the step driver, solver-independent
(it does not use the WP-190 MUMPS timers, so it works for PARDISO, UmfPack, BandGeneral and
MUMPS alike, and the PR does not depend on #978):

| Symbol | Measured as |
|---|---|
| `t_I` | one iteration on a kept factor: triangular solve + Krylov acceleration + `update` (element state determination) + `formUnbalance` + the test. Every iteration whose solve is not the first after a formation. |
| `t_F` | one factorization: `formTangent` (assembly) + the first iteration after it − `t_I` (the first `solve()` after `formTangent` is where every OpenSees SOE factors). In-step formations (`-reuseIter`, `-reuseRate`, the `-maxDim` refresh, the retry, the cost refresh) are sampled the same way. |
| `n_F` | effective iterations of a step that started on a fresh tangent. |
| `n_K(a)` | effective iterations of a step that started on a kept factor of age `a`, modelled as a line `c0 + c1·a`. |

*Effective* iterations count each factorization formed inside the step (beyond the one at its
start) as `t_F / t_I` iterations, so a kept-factor step that had to refresh mid-way is booked
with its true cost. `t_F`, `t_I`, `n_F` are exponential running means (weight 0.25 per sample;
the first sample initialises). `n_K(a)` is a least-squares line through the `(age, effective
iterations)` pairs of the kept-factor steps, with exponential forgetting (0.98 per sample, a
memory of about 50 kept steps), so its slope follows the current phase of the response;
`c1` is clamped at ≥ 0 (a factor does not improve with age); with one distinct age the line is
flat; with no kept step yet, `n_K = n_F` (optimistic, so a kept factor gets tried and measured).

A first version used one running mean per age and the cycle's *measured* iterations on the
right-hand side of §3. On the frame it kept the factor in exactly the phases where a fresh
tangent also needed many iterations (a hard phase inflates the running cycle mean), and
over-refactored on a long-lived factor (each age bucket is visited once per cycle, so the
marginal at age `a` was one noisy sample from another phase): 12 factorizations against the
optimal 2 at `t_F/t_I ≈ 300`. The fitted line with forgetting fixed both (§6.2).

## 3. The decision rule

A factor's life is a renewal cycle: it is formed (cost `t_F`), then serves steps of ages
0, 1, 2, … with costs `t_I · n(a)` (`n(0) = n_F`, `n(a) = n_K(a)`). Over a cycle of `S`
steps the mean cost per step, in units of `t_I`, is

    A(S) = ( t_F/t_I + n_F + Σ_{a=1}^{S-1} n_K(a) ) / S

and the long-run cost per step is minimised by the cycle length that minimises `A`. Since
`A(S+1) < A(S)` ⇔ `n_K(S) < A(S)`, and `n_K` is non-decreasing in age, `A` is unimodal and
the optimal policy is the classical age-replacement rule: **keep while the marginal cost of
the next step is below the cycle's mean cost per step**. With the linear `n_K`, the optimal
cycle is `S* ≈ sqrt(2 (t_F/t_I + n_F − c0) / c1)` (the EOQ form).

**At each step start** (only when the WP-188 hard rules would keep the factor):

    keep  iff  s · n_K(S)  ≤  A(S)

where `S` is the number of steps the current factor has served (including the one that formed
it), so `n_K(S)` is the projected cost of the next step on it.

Limits it reproduces:
- **Factorization cheap** (`t_F/t_I` below `n_K(1) − n_F`, the 1D cell): at `S = 1`,
  `A(1) = t_F/t_I + n_F < n_K(1)`, so every step refactors — full Newton, with one probe step
  every `p` steps.
- **Factorization dominant** (`t_F/t_I ≫ n_K`, 4D): `A(S)` stays above `n_K(S)` for many
  steps; the factor is kept to `S*` or `m`.
- **Linear model**: `n_K = n_F`, `c1 ≈ 0`; the factor is kept to `m` (as WP-188), up to
  timing noise in the fitted slope (§6.1).

**After each non-converged iteration on a kept factor** (ski-rental rule): refactor once the
iterations on it in this step overrun the projection `n_K(a)` by `t_F / t_I`:

    refresh  iff  s · t_I · ( i − n_K(a) )  ≥  t_F

i.e. once the extra iterations have cost as much as the factorization that would have avoided
them (deterministic ski rental: at most 2× the cost of the clairvoyant choice). On 4D
(`t_F/t_I ≈ 20`) it almost never fires before WP-188's `-reuseIter 10`; on a cheap
factorization it fires after the first overrun. `-reuseIter`, `-reuseRate` and the failure
retry still apply on top; a cost refresh is booked under `iter` in the WP-188 line and counted
separately in the WP-193 line.

**Warm-up and probes.** The first `w` steps run on the age rule only (WP-188 forms at steps 1
and 2, so `t_F`, `t_I`, `n_F` are measured; steps 3… start on a kept factor and measure
`n_K(1), n_K(2)`). The rule is active once `t_F`, `t_I` and `n_F` exist. A probe keeps the
factor once after `p` consecutive cost refusals, so a stale `n_K` (measured in a strongly
nonlinear phase, say) is re-measured later in the record.

## 4. MPI

Every rank must take the same keep/refactor decision (WP-188 §5: a split decision puts MUMPS
into job 5 on some ranks and job 3 on others). The inputs here are **wall times**, which
differ per rank. So the per-step timing samples are `MPI_Allreduce(MAX)`-reduced at the step
end (two doubles per step: the slowest rank sets the wall time) before they enter the running
means; the iteration counts are rank-identical already. Every rank therefore holds
bit-identical `t_F, t_I, n_F, n_K(·)` and computes the same decision, which is then
max-reduced anyway, exactly as WP-188 does. Verified: on 2 ranks the WP-193 line, timings
included, is identical on both ranks (`test_mp_ranks_agree_on_the_cost_model`). The double
reduction lives in the per-target TU `LadrunoFactorReuseMPI.cpp` beside the integer one (the
OBJECT-library `#ifdef` trap). With `-reuseCost off` no reduction is added.

The reduction runs at the two common exits of a step (converged; failed after the retry
decision). The early error returns of WP-188 (a `formTangent` / `formUnbalance` failure before
the iteration loop) skip it, as they already skip WP-188's own votes.

## 5. Determinism

With `-reuseCost on` the iteration path depends on measured time, so **two runs of the same
deck can take different keep/refactor decisions** and stop at different points inside the
test tolerance (ADR-188 §4: differences of the order of the tolerance). Results stay within
the test tolerance of full Newton in the same sense as WP-188; they are not bit-reproducible
run to run. For bit-reproducible serial runs use `-reuseCost off`. Where the two costs are
far apart (the regimes that matter, §6.2) the decisions are stable run to run; where they are
close (the frame without hooks, `t_F ≈ t_I`) they flip, and then either choice costs about the
same.

## 6. Tests and measurements

`tests/test_wp193_reuse_cost.py` (Zone-A, 18 cases, about 20 s) reuses the WP-188 frame deck
(3-story 2-bay fiber RC frame, Newmark, `-maxDim 20`, NormDispIncr 1e-6, first 300 steps of the
record). Two **test-only hooks**, read once in the constructor and inert when unset, make the
regimes reproducible on a model that runs in a fraction of a second:
`LADRUNO_WP193_TEST_FACTOR_DELAY_MS` (every factorization the reuse driver forms burns that many
milliseconds of wall time, a busy wait: a large model's factorization) and
`LADRUNO_WP193_TEST_ITER_DELAY_MS` (every iteration burns it after `update`: an expensive state
determination). All runs on the Windows desktop, P-cores only (affinity 0x5555, one thread per
process), OMP/MKL threads 1, 2026-10-14.

### 6.1 Off and the WP-188 tests

- `-reuseCost off` gives the same bytes as no `-reuseCost` (EDP + roof history md5 and the
  statistics line), and no cost line. Against a dist built from the base `origin/ladruno`
  `0c4c6a1d8` in the same worktree: option absent, `-reuseFactor 30` and the elastic model with
  `-reuseFactor 30` are md5-identical, statistics lines identical (`LADRUNO_WP193_BASE_DIST`).
- `tests/test_wp188_factor_reuse.py`: 18/18 on this build (base-dist bit-identity included,
  MP 2/4 ranks included).
- Linear model (`-reuseFactor 50 -reuseCost on`): equals full Newton to round-off; > 90 % of the
  steps on a kept factor (a few cost refactorizations from the noise in the fitted slope).

### 6.2 The trigger follows the cost ratio (frame, 300 steps)

| Regime (hook) | measured `t_F/t_I` | `n_F` / `n_K(1)` | always refactor | always keep (`-reuseFactor 1000`) | **`-reuseCost on`** | steps on kept factor | factorizations |
|---|---|---|---|---|---|---|---|
| iterations expensive (1 ms / iteration) | 0.06 | 2.97 / 3.93 | **1.04 s** | 1.75 s | **1.08 s** (1.03× best) | 26 % | 241 (219 start by cost + 20 in-step cost + 2) |
| factorization expensive (20 ms / factorization) | 197 | 2.25 / 5.27 | 6.15 s | **0.17 s** | **0.21 s** (1.24× best; 0.17–0.21 over reruns) | 99 % | 3 |
| no hook, 600 steps, 3 runs each (`t_F ≈ t_I ≈ 0.06 ms`) | 0.8–1.0 | 2.9 / 3.6–3.9 | 0.18–0.22 s (598 fact., 1881 it.) | 0.29–0.35 s (2 fact., 3213 it.) | 0.19–0.25 s (≈ 265 fact., ≈ 2060 it.) | 56–58 % | ≈ 265 |

"Always refactor" is the same driver with `-reuseCostSafety 1e6 -reuseCostProbe 0`. The keeps in
the iteration-expensive arm are the quiet ramp of the record (first 1 s), where a kept factor
converges as fast as a fresh one, plus 9 probes. Every arm is within 1.2e-5 of full Newton
(peak/residual roof, peak IDR, roof history; tolerance-scale, as WP-188). Without hooks the
frame's `t_F/t_I` is about 1 and `n_K(1) − n_F ≈ 0.7`: the model says keeping one step and
refactoring the next is marginally cheapest, and the trigger alternates. In wall time it is
level with always-refactor and full Newton (0.17–0.24 s) inside the ±20 % run-to-run noise of a
0.2 s run, and 1.5× faster than always-keep.

### 6.3 The test assertions

`test_trigger_follows_the_cost_ratio[iterations_expensive]`: `t_F/t_I < 0.5`, < 40 % of the steps
on a kept factor, cost-driven factorizations > half the steps, wall ≤ 1.2 × always-refactor +
50 ms. `[factorization_expensive]`: `t_F/t_I > 50`, > 90 % kept, ≤ 10 factorizations, wall ≤
1.2 × always-keep + 50 ms. Both: within 1e-4 of full Newton.

### 6.4 A larger model, no hooks (OpenSeesMP + MUMPS, 4 ranks)

`tests/wp188_factor_reuse/perf_block.tcl` (ADR-188 §6.6): 32 × 32 × 16 bricks (≈ 52 k equations),
J2 top 3 layers, 40 steps at 2 g, tol 1e-6, 4 ranks pinned to P-cores 0,2,4,6.

| Arm | wall | factorizations | iterations | corner history vs full Newton |
|---|---|---|---|---|
| full Newton | 92.7 s | 40 | — | — |
| `-reuseFactor 20` (WP-188) | 32.5 s (2.85×) | 5 (3 start + 2 iter) | 242 | 3.7e-7 |
| `-reuseFactor 1000 -reuseCost on` | 32.7 s (2.84×), rerun 35.3 s | 5 (3 start + 2 iter) | 246 | 3.3e-7 |

The model measured `t_F` = 1.55 s, `t_I` = 0.089 s (ratio 17.4, the 4D ratio 19.4 in miniature),
`n_F` = 3, `n_K(1)` = 7.8 with slope 0.05/step: it kept the factor at 35 step starts, refactored
once by cost, and WP-188's `-reuseIter 10` fired twice in the plastic phase. With no `m` to tune
it lands on the WP-188 `m = 20` result.

### 6.5 Model check against the 1D cell (inf.)

From `table_v3_1D.md`: full Newton (metis_blr) 1.08 s/step at 2 iterations and one
factorization; reuse20 1.72 s/step at 6.91 iterations. Solved for the two costs: `t_I ≈ 0.24 s`,
`t_F ≈ 0.6 s` (`t_F/t_I ≈ 2.5`), `n_K − n_F ≈ 4.9 > 2.5`, so at `S = 1` the rule refactors: the
trigger would run the 1D cell as full Newton (≈ 1.08 s/step plus a probe step every 20 steps,
≈ +2 %) instead of 0.66–0.80×. Not run (no esmeralda work in this WP).

## 7. What the 4D validation arms should be compared against

The w7 restage (paper repo, `validation_v3/README.md`) runs, in submission order 904–909:
base_rerun, newton_tol1e-5, metis, metis_blr, metis_blr_reuse20_tol1e-5, metis_blr_reuse20_tol1e-4
— so 908/909 are the two reuse arms (inf., from the order; check `w7/submitted.txt`). Compare:

1. **Speed**: 908/909 against **907 (metis_blr, no reuse)** — the only arm with the same
   ordering, BLR and rank layout; the reuse gain is that ratio, not the ratio to 870.
2. **EDPs / energy**: 908 against **905 (newton_tol1e-5)** at the same tolerance and against the
   870/904 base-rerun scatter; 909 against 907 at the campaign tolerance (ADR-188 §4).
3. **For this trigger**: read from 908's WP-188 line the iterations per kept step and from 907
   `t_F/t_I` (MUMPS factor vs one iteration; WP-190 `-stats` gives the factor phase). If 908's
   kept steps need fewer than about `n_F + t_F/t_I ≈ 2 + 19` iterations on average, a
   `-reuseCost on -reuseFactor 1000` arm should match 908's wall within a few percent (as on the
   block, §6.4) and could keep factors longer than 20 steps where `n_K` stays flat; if 908 is
   slower than 907, the trigger would have refactored (the 1D outcome). A 4D `-reuseCost` arm
   is the follow-up once 907/908 are in.

## 8. Risks and open items

- Not bit-reproducible run to run with the trigger on (§5). Off stays deterministic.
- A wall-clock model sees machine load: a busy node inflates `t_I` and `t_F` together (the ratio
  is what matters), but a burst during one factorization skews `t_F` until the running mean
  forgets it (weight 0.25).
- The probe costs at most one kept step per `p` refactored steps (≈ +2 % on the 1D estimate).
- `t_F` is inferred from the first solve after a formation; an SOE that factored inside
  `formTangent` would book it there instead (same sum).
- Not measured on 4D (no esmeralda in this WP); §7 says what to compare when the arms land.

## 9. Implementation log

- 2026-10-14 — WP number: the brief named 191, but 192 appeared on the remote while this was
  being set up and another agent was reported to be opening a WP for Windows MUMPS + METIS; 193
  was taken to avoid a collision.
- First version: per-age running means and the measured cycle iterations (§2, last paragraph);
  replaced by the fitted line with forgetting after the frame runs.
- The frame without hooks is not an "always refactor" model: its factorization (~60 µs) costs
  about one iteration, and the 1D cell's regime needed the iteration hook.
- Files: `LadrunoFactorReuse.{h,cpp}` (the model, the two rules, the hooks, the statistics line),
  `LadrunoFactorReuseMPI.{h,cpp}` (`ladrunoFactorReuseAgreeMaxDouble`); no vanilla file touched
  (WP-188's `AcceleratedNewton` hooks and parsers already forward every `-reuse*` option).
