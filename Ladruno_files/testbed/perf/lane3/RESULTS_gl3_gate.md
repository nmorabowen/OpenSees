---
title: ADR-75b G-L3 — the Lane 3 re-entry gate, measured at production scale
project: Ladruno
date: 2026-07-26
host: esmeralda, 1 node (32 cores / 60 GB), np=16
binary: OpenSeesMP @ 45a5b55fe
tags: [adr, performance, lane3, openmp, gate]
---

# G-L3 — measured. **Lane 3 FAILS the gate by ~40×, and should be CLOSED for the production path.**

## What G-L3 asked

ADR-75b §12 parked Lane 3 pending exactly one measurement: on a **≥500k DOF**
production-scale deck, what fraction of step does each **individual** assembly loop
hold?

- **Gate: largest SINGLE loop ≥ 40%** ⇒ re-authorize Lane 3 at that loop.
- **No loop reaches 40%** ⇒ close Lane 3 as Amdahl-irrelevant.

The gate is on a **single loop**, not the aggregate element-kernel fraction — Lane D
reads kernel 85.30% but gate **FAIL** at loop A 38.95%.

## Method

`gl3_deck.tcl` (this directory): parametric 3D solid cantilever, `stdBrick` +
**`J2Plasticity`**, `numberer LadrunoParallelRCM`, `system Mumps -matrixType 2`,
`LoadControl 0.35` × 2 steps under full Newton, **coarse** `profiler start -perStep`
(not `-deep`: `enabled()` gates coarse phase timing while `deep()` additionally gates
per-element buckets, `Profiler.h:327-334`, so `-perStep` gives the per-loop fractions
untaxed). Loops map to rollup nodes: **A** = `update`, **B** = `formTangent`,
**C** = `formUnbalance`. Parser: `parse_gl3.py`. Drivers: `gl3.sbatch`,
`gl3_trend.sbatch`.

**Nonlinearity verified, not assumed:** the 540 675-DOF run took **11 `formTangent`
calls for 2 steps = 5.5 Newton iterations/step**. An elastic run would be 1–2. The
material yields and `update` does genuine constitutive work.

## Result

| N (DOF) | loop A `update` | loop B `formTangent` | loop C `formUnbalance` | `linearSolve` | **max loop** | gate |
|---|---|---|---|---|---|---|
| 7 623 | 1.40% | 6.28% | 0.70% | 87.44% | 6.28% | **FAIL** |
| 28 611 | 1.83% | 8.23% | 1.06% | 87.77% | 8.23% | **FAIL** |
| 143 811 | 0.58% | 2.07% | 0.31% | 96.70% | 2.07% | **FAIL** |
| **540 675** | **0.26%** | **0.95%** | **0.13%** | **98.54%** | **0.95%** | **FAIL** |

Rank spread at 540 675 DOF is negligible (ranks 0/1/8/15 give max-loop
0.95/1.13/0.93/1.11%), so this is not a one-rank artifact.

**Per call at 540 675 DOF: one MUMPS solve is 40.2 s; one element tangent assembly
is 0.39 s — a 104× ratio.** The element loops are not a minority of the step, they
are a rounding error on it.

## Verdict — CLOSE Lane 3 for the cluster/production path

The gate fails at the production-scale point by a factor of **~42×** (0.95% vs 40%),
and the trend is **monotonic in N**: from 28.6k → 540.7k the largest loop falls
8.23% → 0.95% while the solve rises 87.8% → 98.5%. Growing the model makes Lane 3
*less* attractive, not more — exactly what ADR-75 §1 predicted ("the solve fraction
is a floor that grows with N").

### Why the deck's caveats do not rescue it — quantified

The Tcl `nDMaterial` command is a hand-written `strcmp` ladder that does **not**
carry `LadrunoJ2` (Python-only), so this deck uses vanilla `stdBrick` +
`J2Plasticity`. Both are **cheaper per element** than `LadrunoBrick` +
`LadrunoJ2` (which runs 16 function-scope statics per `update()`), so these element
fractions are a **lower bound** on the fork's own decks.

**That caveat cannot change the answer.** To lift the largest loop from 0.95% to the
40% gate, the element kernel would have to be **~66× more expensive** at equal solve
cost. `LadrunoBrick` is not 66× `stdBrick`. The margin is wide enough that the proxy
is immaterial.

### What this does NOT close

- **The desktop/PARDISO path.** L3-0's Lane B measured 74.85% element fraction at
  **11 520 DOF** under `Pardiso -matrixType 2` — a far cheaper solve per DOF than
  distributed MUMPS. That number is not contradicted here; it is a different solver
  at a different scale. But the desktop cannot *reach* this regime (P1c: UmfPack
  OOM'd at 86 490 DOF, PARDISO reached 136 080), so any Lane-3 win there is confined
  to models ≲136k DOF.
- ⇒ **Lane 3 is now, at best, a desktop-only optimization on mid-size models** —
  against ADR-75 §1's standing statement that the production regime is huge solid
  nonlinear models, which live on the cluster. Its addressable scope has collapsed.

## Caveats

- **MUMPS at np=16.** The np-sweep (`phase2/RESULTS_p2h_npsweep.md`) shows np=16 is
  near-optimal for symmetric at 143.8k (17.9 s at np=8 vs 18.1 s at np=16), so this
  is not a pathologically bad rank count. At tiny N the small-model rows above are
  dominated by MPI overhead rather than real solve work, which is why only the two
  largest rows should carry weight.
- **n=1 per size.** The effect is ~40× the gate; the phase2 noise floor is ~2.7%.
- Single element type, single material, single box.
