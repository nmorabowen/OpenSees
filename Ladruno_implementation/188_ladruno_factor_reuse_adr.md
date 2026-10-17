---
title: "WP-188 — factorization reuse across time steps (KrylovNewton -reuseFactor)"
project: Ladruno
type: performance work package (solve lane)
status: "implemented, PR open (draft) — not merged; default off"
priority: high
owner: nmora
related:
  - "[[40_ladruno_performance_adr]] (lever 1: factorization reuse)"
  - "[[76_ladruno_tangent_reuse_adr]] (why an engine-side 'tangent unchanged' shortcut was withdrawn)"
  - "paper repo: review/engstruct-2026/30_campaigns/timehistory/perf/PROFILE_4D.md (§4.4, the prototype)"
tags: [performance, solver, algorithm, newton, mumps, mpi, wp-188]
updated: 2026-10-11
---

# WP-188 — factorization reuse across time steps

> [!summary] The short version
> The Paper 1 Tier 3–4 models (linear soil box + nonlinear building, OpenSeesMP + MUMPS,
> about 208 k equations) spend 97.6 % of each step in MUMPS, and that is one LU per step
> (12.8 s) plus 0.66 s per iteration (PROFILE_4D §2). `algorithm KrylovNewton` already
> reuses the LU *inside* a step; WP-188 lets it reuse the LU *across* steps:
>
> ```tcl
> algorithm KrylovNewton -maxDim 20 -reuseFactor 20 ?-reuseIter 10? ?-reuseRate 1.0? ?-reuseReport n?
> ```
>
> A step may start on the factor of an earlier step (modified Newton, Krylov-accelerated),
> with explicit refactorization triggers and a refactor-and-retry before any failure is
> reported. Off by default and byte-identical when off. It is not a change of the
> equations: every step still passes the user's convergence test; it stops at a different
> point inside the tolerance than full Newton does. On the test frame, reuse at a 10×
> tighter test differs from full Newton by no more than full Newton at the original test
> differs from the exact solution, at 4.6–5.5× less solver work (4D cost ratio).

## 1. What

- A fork option on the `KrylovNewton` algorithm command, in both parsers (classic Tcl
  `SRC/tcl/commands.cpp`, which `OpenSees.exe`, `OpenSeesSP.exe` and **`OpenSeesMP.exe`**
  use, and the interpreter `OPS_KrylovNewton` in `SRC/interpreter/OpenSeesCommands.cpp`,
  which openseespy and openseesmp use). Both build `AcceleratedNewton` + `KrylovAccelerator`
  (`KrylovNewton.cpp` is never built by `algorithm KrylovNewton`; PROFILE_4D §4.4 lost a
  probe to that).
- A fork class `LadrunoFactorReuse` (`SRC/analysis/algorithm/equiSolnAlgo/LadrunoFactorReuse.{h,cpp}`)
  that holds the policy and the statistics and drives the step when enabled.
- A per-target translation unit `LadrunoFactorReuseMPI.{h,cpp}` (in
  `OPS_MPI_PER_TARGET_SOURCES`) with the `MPI_Allreduce` that makes every rank take the
  same refactorization decision.
- Vanilla footprint: `AcceleratedNewton.{h,cpp}` (a null-by-default pointer, an early
  branch at the top of `solveCurrentStep`, a `domainChanged()` override, `Print`), the two
  parsers, two CMake lists. Every edit is marked `// Ladruno WP-188`.

Not in scope: the other `AcceleratedNewton` spellings (`RaphsonNewton`, `SecantNewton`,
`PeriodicNewton`, `MillerNewton`) — the class supports them, the parsers do not expose the
option (a one-line follow-up each if wanted); `Newton` / `ModifiedNewton` (different
classes; `ModifiedNewton` has no Krylov acceleration and the profile says the accelerated
path is what the decks run).

## 2. Why the knob is an algorithm option

The decision "may this step start on the old factor, and when must it be replaced" needs
three things: the iteration count and the convergence-test norms of the current step, the
ability to re-form the tangent at the current iterate, and the ability to restart the
Krylov subspace (whose vectors belong to the old factor). Only the solution algorithm has
all three. The alternatives were considered and rejected:

| Candidate | Why not |
|---|---|
| SOE / MUMPS option (`system Mumps -keepFactor`) | The SOE sees `zeroA`/`addA`/`solve` and its `factored` flag, never the convergence history; it cannot tell a slow step from a fast one, cannot re-form the tangent, and would have to be taught to *ignore* a `formTangent` the algorithm asked for. It would also be MUMPS-only, while the policy is solver-independent (PARDISO, UmfPack, BandGeneral gate on `factored` the same way). |
| Integrator option (`integrator Newmark ... -reuse`) | The integrator forms the tangent on request; it does not own the iteration loop or the test. Every implicit integrator would need the option. |
| A new algorithm class (`LadrunoKrylovNewton`) | Clean in principle, but `algorithm KrylovNewton` is what the decks already say; a new name means every deck edit plus duplicating `AcceleratedNewton::solveCurrentStep`. The option keeps decks one token away from on/off and lets the off path stay literally the vanilla function. |

This also fits ADR-76's conclusion: the engine should not guess whether a tangent "has
changed" (it cannot, §3 there); the reuse is an explicit, user-chosen modified-Newton
strategy whose correctness never depends on the tangent being current.

## 3. How

### 3.1 Syntax

```tcl
algorithm KrylovNewton ?-iterate t? ?-increment t? ?-maxDim n? \
    -reuseFactor m ?-reuseIter k? ?-reuseRate r? ?-reuseReport n?
```

| Option | Meaning | Default |
|---|---|---|
| `-reuseFactor m` | Enable. A factor may serve at most `m` steps after the step that formed it. `0` = off (same as absent). | off |
| `-reuseIter k` | A step running on a kept factor that has done `k` iterations on it without converging forms the current tangent, restarts the Krylov subspace and continues. | 10 |
| `-reuseRate r` | Same, as soon as the test-norm ratio `n_i / n_{i-1}` on a kept factor exceeds `r` (i ≥ 2). `0` disables. Uses whatever norm the test records (displacement increment, unbalance, energy), so with `NormUnbalance` it is the residual criterion. | 1.0 (divergence only) |
| `-reuseReport n` | Print the running statistics every `n` steps (rank 0). | 0 (summary at the end only) |

Fail-closed parsing: a `-reuse*` tuning option without `-reuseFactor`, a negative or
non-numeric value, or a missing value refuses the algorithm command (Tcl error / Python
exception).

### 3.2 The step

`LadrunoFactorReuse::solveCurrentStep` is `AcceleratedNewton::solveCurrentStep` with
three decisions added; the residual, update, test, profiler scopes and their order are
the vanilla ones.

1. **Form at step start, or start on the kept factor.** Form if any of: no valid factor;
   it has already started `m` steps; the step size `dt` (recovered as the current domain
   time minus the time of the last converged step) differs from the one it was formed at
   by more than 1e-8 relative; the SOE, the integrator or the number of equations differs;
   `domainChanged()` was called since (the override invalidates). The first step of an
   analysis always forms; so does the second (the algorithm cannot know step 1's `dt`).
   Otherwise the step skips `formTangent`, the SOE stays `factored`, and MUMPS runs job 3.
2. **Replace a kept factor that converges too slowly** (`-reuseIter`, `-reuseRate`), at
   the point where `KrylovAccelerator::updateTangent` would replace it; the vanilla
   `-maxDim` refresh still applies on top.
3. **Refactor and retry.** If a step that *started* on a kept factor fails (test returns
   failure, or a solve / update / accelerate error), the step's accumulated increments are
   undone (`update(-ΣΔU)`; exact for every integrator whose update is linear in ΔU — all
   transient integrators and `LoadControl`; path-following controls retry from the last
   iterate), the tangent is formed at the predictor, and the step is run once more as plain
   KrylovNewton with the full iteration budget. Only if that fails too is the failure
   returned to the analysis (and the deck's own step-cutting takes over, which changes `dt`
   and therefore forces a factorization anyway).

Statistics (rank 0, one line; `print -algorithm` prints it on every rank):

```
KrylovNewton factor reuse (WP-188) final: steps 1000, started on a kept factor 966 (96.6 %),
factorizations 34 = start 34 + iter 0 + rate 0 + krylov 0 + retry 0, retries 0 (failed 0),
iterations 3201 (on kept factors 3105)
```

The line is printed by the destructor (`wipe`, `wipeAnalysis`, a new `algorithm`), with
`-reuseReport n` every `n` steps, and by `print -algorithm` ("so far"). A one-line
announcement at the first step states the settings and the tolerance recommendation (§4).

## 4. Correctness: why the answers differ from full Newton, and what controls it

Nothing in the equations changes. The residual of every iteration is the exact residual of
the current trial state, and a step ends only when the user's test passes, so the converged
state satisfies **the same convergence test** as without reuse. What changes is the
iteration path, hence the point inside the test's tolerance ball where iteration stops.

- **Full Newton** converges quadratically: when the increment norm first drops below `tol`,
  the remaining error is of order `|Δu|²` — usually far below `tol`.
- **Modified Newton** (an older factor, Krylov-accelerated) converges linearly with some
  contraction `ρ`: the step stops at an increment just below `tol`, with a remaining error
  of order `ρ/(1−ρ) · tol`.

So two different, equally valid solutions of the same test differ by an amount of the order
of the tolerance per step, and in a path-dependent model (cracking, yielding) that
difference propagates into the response history. It is **controlled by the tolerance**, not
by `m`: the measurements below show the reuse-vs-full-Newton difference scaling linearly
with `tol` (each decade of `tol` takes off about a decade), saturating in `m` once most
steps start on a kept factor.

**Recommendation (printed at the first step):** for result-neutral runs, tighten the test
about 10× when reuse is on. Reuse at `tol/10` then differs from full Newton at `tol` by no
more than full Newton at `tol` differs from the converged-to-round-off solution (frame:
3.3e-4 vs 3.7e-4 of the peak roof displacement), and costs 4.6–5.5× less solver work than
full Newton at `tol`. The extra iterations of the tighter test are cheap where reuse pays
(one iteration is 1/19 of a factorization on 4D). A `-reuseTighten f` option that tightens
only the reused steps was considered and not built: the global tighter test costs at most
one extra iteration on the full-Newton steps, and a single test is easier to report.

This also explains the prototype's 4D F-window result (PROFILE_4D §4.4: element forces and
energy 3–5× outside the rerun scatter at NormDispIncr 1e-4): the prototype ran at the
campaign tolerance, where the "difference" is the tolerance, not an error of the method.

## 5. MPI

Under OpenSeesMP every rank runs the algorithm in lock-step on its own `Domain`, and a
rank that formed the tangent while another did not would put MUMPS into job 5 on some
ranks and job 3 on others — a hang or garbage. The inputs of every decision are
rank-identical in practice (the test norms are computed from the broadcast solution, `dt`
and the counters are the same), but they are **reduced, not assumed**: each of the three
decisions is an `MPI_Allreduce(MAX)` over `MPI_COMM_WORLD`, in
`LadrunoFactorReuseMPI.cpp`, which is compiled per target (the OBJECT-library `#ifdef`
trap, ADR-78 P1 / `LadrunoContactAbort.h`). The reductions run only when the option is on;
one `int` per step plus one per non-converged iteration on a kept factor.

Verified live, not by reading: the `build.ninja` DEFINES give that TU
`_PARALLEL_INTERPRETERS` in OpenSeesMP / OpenSeesPyMP and `_PARALLEL_PROCESSING` in
OpenSeesSP; and a test hook (`LADRUNO_WP188_TEST_VOTE_RANK=1`: rank 1 votes "form" at every
step start) makes **rank 0** form at all 1000 steps — rank 0 can only learn that through
the reduction.

Limits: `system Mumps -commSplit` (several independent models on sub-communicators in one
job) is not supported with reuse on — the reduction spans `MPI_COMM_WORLD`. OpenSeesSP was
built (the TU compiles with its define) but not run (inf.).

## 6. Measurements

All on the 24-core Windows desktop, this branch's build, 2026-10-11.

### 6.1 Off = vanilla (bit-identity)

The base binary is `origin/ladruno` `c50a281a0` built in the same worktree before the edits.
The frame deck without the option, serial (`OpenSees.exe`): md5 of the EDP line + the roof
history **identical** to the base build in 7/7 configurations (fiber at tol 1e-8 / 1e-4 /
1e-5 with `-maxDim 20`, elastic, UmfPack, a dt change, and `-reuseFactor 0`). An existing
KrylovNewton deck (`tests/wp133_pdmy03_deck.py`, PDMY03) gives md5-identical JSON on both
builds. OpenSeesMP (MUMPS, 2 and 4 ranks) is **not byte-reproducible run to run on either
build** (three reruns of the base binary give three md5s, differences 1e-16 relative —
M39), so the MP off-path is checked to that level only.

### 6.2 Linear model (elastic frame, tol 1e-10)

Reuse equals full Newton to 1.7e-15 of the peak roof displacement, and forms exactly
`1 + ceil((N−1)/(m+1))` factors: 2 at `m = 1000`, 201 at `m = 4` (N = 1000; full Newton
forms 1000).

### 6.3 Nonlinear frame: EDP difference vs `m` and the test tolerance

3-story, 2-bay RC frame (fiber forceBeamColumn, Concrete01/Steel02, P-Delta), synthetic
ground motion, 0.6 g, Newmark, 1000 steps of 0.01 s, peak roof drift 0.21 m (3 % IDR,
plastic), `KrylovNewton -maxDim 20`, NormDispIncr `tol` (m) / 100. EDPs: peak roof
displacement, peak IDR, residual roof displacement, peak base shear, peak column-base
moment. Differences are relative to the full-Newton value (residual: to the peak roof
displacement); the history column is the maximum roof-displacement difference over the
record. Cost units weight a factorization 19.4 iterations (the 4D ratio 12.8 s / 0.66 s),
per step. Default triggers (`-reuseIter 10 -reuseRate 1.0`); none fired — the frame's kept
factors always converged within 10 iterations, every factorization is an `m` expiry.

| tol | m | factorizations | iterations | steps on kept factor | cost units | speed-up vs full Newton at 1e-4 | max EDP diff vs full Newton, same tol | roof history diff | full Newton at this tol vs tol 1e-12 |
|---|---|---|---|---|---|---|---|---|---|
| 1e-04 | off | 1000 | 2006 | - | 21.41 | 1.0x | - | - | 3.7e-04 |
| 1e-04 | 1 | 501 | 2026 | 49.9 % | 11.75 | 1.8x | 2.4e-04 | 3.0e-04 | 3.7e-04 |
| 1e-04 | 3 | 251 | 2056 | 74.9 % | 6.93 | 3.1x | 5.7e-04 | 5.8e-04 | 3.7e-04 |
| 1e-04 | 10 | 92 | 2113 | 90.8 % | 3.90 | 5.5x | 1.9e-03 | 1.3e-03 | 3.7e-04 |
| 1e-04 | 30 | 34 | 2210 | 96.6 % | 2.87 | 7.5x | 1.6e-03 | 3.0e-03 | 3.7e-04 |
| 1e-04 | 100 | 11 | 2306 | 98.9 % | 2.52 | 8.5x | 2.1e-03 | 4.1e-03 | 3.7e-04 |
| 1e-05 | off | 1000 | 2285 | - | 21.68 | 1.0x | - | - | 1.5e-05 |
| 1e-05 | 1 | 501 | 2439 | 49.9 % | 12.16 | 1.8x | 1.5e-05 | 1.9e-05 | 1.5e-05 |
| 1e-05 | 3 | 251 | 2620 | 74.9 % | 7.49 | 2.9x | 3.1e-05 | 3.4e-05 | 1.5e-05 |
| 1e-05 | 10 | 92 | 2904 | 90.8 % | 4.69 | 4.6x | 6.8e-05 | 8.7e-05 | 1.5e-05 |
| 1e-05 | 30 | 34 | 3201 | 96.6 % | 3.86 | 5.5x | 2.6e-04 | 3.1e-04 | 1.5e-05 |
| 1e-05 | 100 | 11 | 3394 | 98.9 % | 3.61 | 5.9x | 3.0e-04 | 3.3e-04 | 1.5e-05 |
| 1e-06 | off | 1000 | 2945 | - | 22.34 | 1.0x | - | - | 2.0e-06 |
| 1e-06 | 1 | 501 | 3149 | 49.9 % | 12.87 | 1.7x | 1.8e-06 | 7.8e-07 | 2.0e-06 |
| 1e-06 | 3 | 251 | 3394 | 74.9 % | 8.26 | 2.6x | 3.0e-06 | 1.4e-06 | 2.0e-06 |
| 1e-06 | 10 | 92 | 3857 | 90.8 % | 5.64 | 3.8x | 6.6e-06 | 3.0e-06 | 2.0e-06 |
| 1e-06 | 30 | 34 | 4380 | 96.7 % | 5.04 | 4.2x | 4.0e-06 | 5.1e-06 | 2.0e-06 |
| 1e-06 | 100 | 12 | 4687 | 98.9 % | 4.92 | 4.4x | 9.3e-06 | 8.0e-06 | 2.0e-06 |
| 1e-08 | off | 1000 | 4222 | - | 23.62 | 0.9x | - | - | 6.4e-09 |
| 1e-08 | 1 | 501 | 4495 | 49.9 % | 14.21 | 1.5x | 1.1e-08 | 5.2e-09 | 6.4e-09 |
| 1e-08 | 3 | 251 | 4899 | 74.9 % | 9.77 | 2.2x | 2.2e-08 | 5.9e-09 | 6.4e-09 |
| 1e-08 | 10 | 92 | 5662 | 90.8 % | 7.45 | 2.9x | 1.2e-08 | 1.7e-08 | 6.4e-09 |
| 1e-08 | 30 | 34 | 6631 | 96.7 % | 7.29 | 2.9x | 3.9e-08 | 5.8e-08 | 6.4e-09 |
| 1e-08 | 100 | 15 | 7080 | 99.1 % | 7.37 | 2.9x | 8.6e-08 | 8.7e-08 | 6.4e-09 |

Reading: (1) the difference scales with `tol` (about one decade per decade) and is
0.5–10× `tol / peak` (`tol/peak` = 4.8e-4 at 1e-4); (2) it grows with `m` up to about
`m = 30` and then saturates; (3) reuse at `tol/10` stays within the full-Newton-vs-exact
error at `tol` (row 1e-5, any m ≤ 100: ≤ 3.3e-4 against 3.7e-4); (4) the speed-up saturates
because a kept factor needs more iterations (2.0 → 3.4 per step at 1e-5); on 4D that is
cheap. The same sweep with the vanilla `-maxDim 3` (the KrylovNewton default) gains less at
tight tolerances because the `-maxDim` refresh fires after 4 iterations (tol 1e-8: 483
factorizations at any m ≥ 10, 3.3e-9 difference) — use `-maxDim 20` (the decks already do).

### 6.4 Triggers

- `dt` change: halving `dt` for 20 steps at step 300 with `m = 1000` forms 4 factors
  (step 1, step 2, the halving, the restore) instead of 2.
- `-reuseIter 3` (tol 1e-6): 316 iteration-triggered factorizations, run completes.
  `-reuseRate 0.3`: 46 rate-triggered ones.
- **Forced non-convergence and retry** (tol 1e-8, test max 7 iterations, `-reuseIter 100
  -reuseRate 0` so nothing refreshes in-step; full Newton needs at most 6): 52 steps fail on
  a kept factor, all 52 are undone and retried and converge; the run completes all 1000
  steps, roof history within 1e-5 of full Newton. With max 6 iterations (full Newton itself
  needing exactly 6 on 48 steps) one retry fails at step 118, where the retried trajectory
  needs a 7th iteration: the retry is plain KrylovNewton with the user's budget, so a deck
  on the edge of its iteration limit stays on the edge.

### 6.5 OpenSeesMP + MUMPS (2 and 4 ranks, `mpiexec -localonly`)

The frame split by element over the ranks (Transformation, `LadrunoParallelPlain`,
`system Mumps -ICNTL14 200`): full Newton MP = serial to 1.5e-15; reuse (`m = 30`,
tol 1e-6) completes on 2 and 4 ranks, every rank reports identical statistics (96.7 % of
the steps on a kept factor, 34 factorizations), roof history within 5.1e-6 of full Newton
MP. The forced-retry case on 2 ranks (tol 1e-8, max 7) retries and converges with identical
per-rank counts; the vote hook (§5) makes rank 0 form at every step.

### 6.6 Wall time, OpenSeesMP + MUMPS, 4 ranks (`tests/wp188_factor_reuse/perf_block.tcl`)

A 3D brick block (stdBrick, elastic below, J2Plasticity in the top layers) under a base
acceleration — the 4D shape in miniature (one large factor per step, a plastic zone that
keeps the tangent changing). `KrylovNewton -maxDim 20`, NormDispIncr, Newmark, 4 ranks.

32 × 32 × 16 hexes (≈ 52 k equations after Transformation), top 3 layers plastic, 40 steps
of 0.01 s, 2 g base acceleration (peak corner drift 0.21 m), `-reuseFactor 20`, default
triggers. One run per arm, dedicated desktop, no other load (timings ±5 %, (inf.) — single
samples).

| Arm | tol | wall (40 steps) | s/step | factorizations | iterations | corner history vs full Newton 1e-6 |
|---|---|---|---|---|---|---|
| full Newton | 1e-6 | 89.6 s | 2.24 | 40 | (not counted) | — |
| full Newton | 1e-7 | 89.9 s | 2.25 | 40 | (not counted) | 3.5e-7 |
| reuse m = 20 | 1e-6 | 36.0 s | 0.90 (**2.5×**) | 5 (3 start + 2 iter) | 242 | 3.7e-7 |
| reuse m = 20 | 1e-7 | 44.1 s | 1.10 (**2.0×** vs full Newton 1e-6) | 7 (2 start + 5 iter) | 267 | 3.6e-7 (3.4e-8 vs full Newton 1e-7) |

On this block a factorization is worth roughly 13 iterations (inf., from the arms' s/step
and iteration counts; the 4D factor is 208 k equations and 19 iterations' worth), and a kept
factor needs about 6 iterations a step, so the frame's cost-unit table (§6.3) is the better
predictor for 4D; the block shows the mechanism end to end under
MUMPS: 92–95 % of the steps on a kept factor, the iteration trigger firing in the plastic
phase, and the tightened test (1e-7) staying 2× faster than full Newton at 1e-6 while
landing 10× closer to full Newton at 1e-7.

## 7. Recommended settings for the 4D validation matrix

The frame numbers say what the knob does; they do not say how long a kept factor survives
in 4D's cracked strong phase, which is the open number (PROFILE_4D §5.3). Recommended
matrix, each arm a full record on the rev3 deck with gravity (inf. until run):

| Arm | algorithm | test | Purpose |
|---|---|---|---|
| A0 | `KrylovNewton -maxDim 20` | NormDispIncr 1e-4 (as today) | reference |
| A0' | same, rerun | same | the M39 rerun scatter |
| A1 | `... -reuseFactor 20` | NormDispIncr **1e-5** | **the candidate**: result-neutral setting (§4) |
| A2 | `... -reuseFactor 20` | NormDispIncr 1e-4 | the difference at the campaign tolerance |
| A3 | `... -reuseFactor 5` | NormDispIncr 1e-5 | sensitivity to m |
| A4 | A1 + `-reuseIter 6` | NormDispIncr 1e-5 | cheaper refresh in the strong phase |

Defaults `-reuseIter 10 -reuseRate 1.0` are the recommendation for A1: on 4D one iteration
costs 1/19 of a factorization, so a kept factor that needs up to ~10 extra iterations still
pays, and the retry covers the steps where it does not. Accept A1 if its EDP and energy
differences to A0 are within the A0/A0' scatter, or within the A0-at-1e-4 vs A0-at-1e-5
difference; read the reuse summary line for the fraction of steps on a kept factor. Use
`-reuseReport 500` to follow it in the log. BLR (`-BLR`) composes with reuse (cheaper
refreshes); combine only after A1 is accepted on its own.

## 8. Risks and open items

- **The 4D strong-phase survival of a kept factor is unmeasured** (inf.); it decides
  whether the record gain is 1.8× or 4×+ (PROFILE_4D §5.1). The validation matrix measures it.
- Anything that rewrites the SOE matrix between steps without a domain change, an `dt`
  change or a new analysis object (none known in the decks: recorders, EnergyBalance, H5DRM
  do not touch `A`) would leave a stale *matrix*, not a wrong answer: the residual is still
  exact and the test still decides; only convergence speed would suffer.
- Integrators whose tangent coefficients change at constant `dt` (TRBDF2/3 sub-steps) still
  converge but reuse a factor of the other sub-step's coefficients (slower, not wrong).
- `-commSplit` unsupported with reuse on (§5); OpenSeesSP not run (inf.).
- `ladrunoSchema` v0 lists material options only; the algorithm options are not in it.
  Feature detection is `ladrunoBuild` ≥ the merge commit.
- The statistics line is printed by the destructor; a process killed by the scheduler
  prints only the `-reuseReport` lines.

## 9. Upstream status

Not filed upstream (by instruction). The change is self-contained, default off, and touches
`AcceleratedNewton` only through a null pointer and an early branch, so it is a candidate
for an upstream PR once the 4D validation is in; the per-target MPI TU is fork
infrastructure that upstream would replace with its own `_PARALLEL_*` handling.

## 10. Implementation log

- 2026-10-11 — prototype read (`scratch/perf-4d-scope`, environment switches, static
  globals, refresh after `m` iterations, no retry, no MPI agreement). Kept: the place
  (AcceleratedNewton), the "refresh once and restart the subspace" move. Replaced: the
  environment switches (option), the static guard (members + `domainChanged` + `dt` +
  integrator identity), the iteration-count semantics of `m` (now the step age; the
  iteration limit is `-reuseIter`), and added the rate trigger, the retry, the reduction
  and the statistics.
- Surprise 1: `OpenSeesMP.exe` parses `algorithm` in `SRC/tcl/commands.cpp` (classic), not
  in `OpenSeesCommands.cpp`; both needed the option.
- Surprise 2: the prototype's "not round-off identical" was the tolerance (§4), not a bug.
- Surprise 3: OpenSeesMP + MUMPS is not byte-reproducible run to run even on a 300-DOF
  frame at 2 ranks; the MP bit-identity gate is therefore serial-only.
- Pre-existing, unrelated: `tests/test_wp135_pdmy_substep_cap.py::test_h1_two_element_pair_past_crossing_returns`
  fails on the base build too (`n2 = 83`, expected 15–25).
