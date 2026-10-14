---
title: "WP-193 — cost-aware refactorization trigger for factorization reuse (KrylovNewton -reuseCost)"
project: Ladruno
type: performance work package (solve lane)
status: "in progress — draft PR, not merged; default off"
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
> The WP-188 triggers stay as hard limits.

## 1. Syntax

| Option | Meaning | Default |
|---|---|---|
| `-reuseCost on\|off` | Enable the cost model. Needs `-reuseFactor m` (fail closed, like every `-reuse*` option); `m` stays the hard age cap, so `-reuseFactor 1000 -reuseCost on` lets the cost model alone decide. | off |
| `-reuseCostSafety s` | `s > 0`. Multiplies the projected cost of continuing on the kept factor, in both rules (§3). `s > 1` refactors earlier (closer to full Newton), `s < 1` keeps longer. | 1.0 |
| `-reuseCostWarmup w` | The first `w` steps follow the WP-188 age rule only, while the costs are learned. | 4 |
| `-reuseCostProbe p` | After `p` consecutive step starts at which the cost rule refused the kept factor, keep it once, so a stale "reuse does not pay" estimate is re-measured (the tangent and the costs change during a record). `0` = never. | 20 |

The tuning options without `-reuseCost on` refuse the command. Statistics: the WP-188 line,
unchanged, plus one line

```
KrylovNewton reuse cost model (WP-193) final: t_factor 0.0123 s, t_iter 0.00061 s, ratio 20.1,
n_fresh 2.0, n_kept1 3.1, starts kept by cost 960, starts formed by cost 3,
in-step cost refreshes 1, probes 0
```

## 2. What is measured

All wall clock (`std::chrono::steady_clock`), inside the step driver, solver-independent
(it does not use the WP-190 MUMPS timers, so it works for PARDISO, UmfPack, BandGeneral):

| Symbol | Measured as |
|---|---|
| `t_I` | one iteration on a kept factor: triangular solve + Krylov acceleration + `update` (element state determination) + `formUnbalance` + the test. Every iteration whose solve is not the first after a formation. |
| `t_F` | one factorization: `formTangent` (assembly) + the first iteration after it − `t_I` (the first solve after `formTangent` is where every OpenSees SOE factors). In-step formations (`-reuseIter`, `-reuseRate`, the `-maxDim` refresh, the retry) are sampled the same way. |
| `n_F` | effective iterations of a step that started on a fresh tangent. |
| `n_K(a)` | effective iterations of a step that started on a kept factor of age `a` (one table entry per age, `a` = 1…256, ages above 256 share the last). |

*Effective* iterations count each factorization formed inside the step (beyond the one at its
start) as `t_F / t_I` iterations, so a kept-factor step that had to refresh mid-way is booked
with its true cost. All five are exponential running means with weight 0.25 (the first
sample initialises). Ages not yet measured are extrapolated linearly from the two highest
measured ages below (with `n_F` as age 0), never decreasing; with nothing measured,
`n_K = n_F` (optimistic, so a kept factor gets tried and measured).

## 3. The decision rule

A factor's life is a renewal cycle: it is formed (cost `t_F`), then serves steps of ages
0, 1, 2, … with costs `t_I · n(a)` (`n(0) = n_F`, `n(a) = n_K(a)`). Over a cycle of `L`
steps the mean cost per step is

    A(L) = ( t_F + t_I · Σ_{a<L} n(a) ) / L

and the long-run cost per step is minimised by the `L` that minimises `A`. Since
`A(L+1) < A(L)` ⇔ `t_I · n(L) < A(L)`, and `n` is non-decreasing in age (the tangent drifts
from the factor), `A` is unimodal and the optimal policy is the classical age-replacement
rule: **keep while the marginal cost of the next step is below the running mean of the
current cycle**.

**At each step start** (only when the WP-188 hard rules would keep the factor):

    keep  iff  s · t_I · n_K(a+1)  ≤  ( t_F + t_I · I_cyc ) / S_cyc

where `a` is the factor's current age, `S_cyc` the steps it has served (including the one that
formed it) and `I_cyc` the iterations actually run on it. The right-hand side uses the
measured iterations of this cycle, not the projection, so a cycle that went badly ends
sooner.

Limits it reproduces:
- **Factorization cheap** (`t_F ≪ t_I`, the 1D cell): the right-hand side at age 1 is
  ≈ `t_I · n_F`, and `n_K(1) > n_F`, so every step refactors — full Newton, with one probe
  step every `p` steps.
- **Factorization dominant** (`t_F ≫ t_I · n`, 4D): the right-hand side stays large for many
  steps; the factor is kept until `n_K(a)` has grown to about the cycle mean, or `m`.
- Linear model: `n_K = n_F = 1` at every age, the factor is kept to `m` (as WP-188).

**After each non-converged iteration on a kept factor** (ski-rental rule): refactor once the
iterations on it in this step overrun the projection `n_K(a)` by `t_F / t_I`:

    refresh  iff  s · t_I · ( i − n_K(a) )  ≥  t_F

i.e. once the extra iterations have cost as much as the factorization that would have
avoided them (the deterministic ski-rental rule; at most 2× the cost of the clairvoyant
choice). On 4D (`t_F/t_I ≈ 20`) it almost never fires before WP-188's `-reuseIter 10`; on a
cheap factorization it fires after the first overrun. `-reuseIter`, `-reuseRate` and the
failure retry still apply on top; a cost refresh is booked under `iter` in the WP-188 line
and counted separately in the WP-193 line.

**Warm-up and probes.** The first `w` steps run on the age rule only (step 1 and step 2
always form under WP-188, so `t_F`, `n_F` and `t_I` are measured; steps 3… start on a kept
factor and measure `n_K(1), n_K(2)`). The rule is active once all three of `t_F`, `t_I`, `n_F`
exist. A probe keeps the factor once after `p` consecutive cost refusals, so a stale
`n_K(1)` (measured in a strongly nonlinear phase, say) is re-measured later in the record.

## 4. MPI

Every rank must take the same keep/refactor decision (WP-188 §5: a split decision puts MUMPS
into job 5 on some ranks and job 3 on others). The inputs here are **wall times**, which
differ per rank. So the per-step timing samples are `MPI_Allreduce(MAX)`-reduced at the step
end (two doubles per step: the slowest rank sets the wall time), before they enter the
running means; the iteration counts are rank-identical already. Every rank therefore holds
bit-identical `t_F, t_I, n_F, n_K(·)` and computes the same decision, which is then
max-reduced anyway, exactly as WP-188 does. The double reduction lives in the per-target TU
`LadrunoFactorReuseMPI.cpp` beside the integer one (the OBJECT-library `#ifdef` trap). With
`-reuseCost off` no reduction is added.

The reduction runs at the two common exits of a step (converged; failed after the retry
decision). The early error returns of WP-188 (a `formTangent` / `formUnbalance` failure
before the iteration loop) skip it, as they already skip WP-188's own votes.

## 5. Determinism

With `-reuseCost on` the iteration path depends on measured time, so **two runs of the same
deck can take different keep/refactor decisions** and stop at different points inside the
test tolerance (ADR-188 §4: differences of the order of the tolerance). Results stay within
the test tolerance of full Newton in the same sense as WP-188; they are not bit-reproducible
run to run. For bit-reproducible runs use `-reuseCost off` (WP-188 is deterministic in
serial). The `-reuseCostSafety` knob shifts the decisions but cannot make them deterministic.

## 6. Tests and measurements

`tests/test_wp193_reuse_cost.py` (Zone-A, needs the dist binaries). A test-only hook,
`LADRUNO_WP193_TEST_FACTOR_DELAY_MS`, makes every factorization the reuse driver forms burn
that many milliseconds of wall time (a busy wait), standing in for a large model's
factorization on a small one; it is read once in the constructor and is inert when unset.

(Results: §6.1–6.4 below, filled in from the runs on this branch.)

## 7. What the 4D validation arms should be compared against

(see §6.5)

## 8. Implementation log

- 2026-10-14 — WP number: 191 was named by the brief, but 192 appeared on the remote while
  this was being set up and another agent was reported to be opening a WP for Windows
  MUMPS + METIS; 193 was taken to avoid a collision.
