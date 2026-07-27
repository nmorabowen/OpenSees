---
title: Which solver flags for which deck — the ADR-75 results, as instructions
project: Ladruno
owner: nmora
relates: 75_ladruno_sparse_direct_strategy_adr
tags: [guide, performance, solver, mumps, pardiso]
---

# Which solver flags for which deck

ADR-75 measured a large amount of free performance. **None of it is on by default.**
This page exists so the measurements turn into speed instead of sitting in a results
file — the failure mode this program has already paid for once, when
`-commitSolveState` shipped, no deck set it, and it cost **−29.6%** until a review
tripped over it a year later.

Everything below is measured, with the source named. Nothing here is a projection.

---

## The one-line answer

| you are on | use |
|---|---|
| desktop / workstation, ≲136k DOF | `system Pardiso -matrixType 2` (add `-krylov 6` under full Newton) |
| cluster, any size | `system Mumps -ICNTL14 200 -matrixType 2` |
| either, and the tangent is **not** symmetric | drop `-matrixType 2` — see the trap below |

---

## `-matrixType 2` — the single biggest free win

**What it buys** (`phase2/RESULTS_p2h_cluster_blr.md`, `phase2/RESULTS_p2h_npsweep.md`):

| where | speed | memory |
|---|---|---|
| cluster MUMPS, 143.8k DOF, np=16 | **2.32× faster** | −42.5% per-proc peak |
| cluster MUMPS, 408.5k DOF, np=16 | **1.73× faster** | −43.1% per-proc peak |
| cluster MUMPS, 884.8k DOF, np=32 | **1.69× faster** | −43.2% total factor |
| desktop PARDISO, 51k DOF (P1d) | **~1.25×** vs unsymmetric | −41.8% peak |

It is **exact**, not approximate: at 884.8k DOF the symmetric and full-rank answers are
**bit-identical**. So unlike `-BLR` it is legal on oracle/QA paths.

### ⚠ The trap — when it is WRONG

Half-storage reads only one triangle. On a **genuinely unsymmetric tangent it silently
solves the wrong system** — no error, no warning, a plausible wrong answer. Do **not**
use it with:

- **contact** (`LadrunoContact`) — frictional tangents are unsymmetric
- **non-associated plasticity** (Drucker-Prager with ψ≠φ, cap models)
- **follower / pressure loads** (load-stiffness terms)
- **`LadrunoUP`** — the Biot u–p tangent is unsymmetric by construction
- any element whose `getTangentStiff()` you have not confirmed symmetric

Associated-flow J2, elasticity, and standard solid/shell decks are fine.
**This is why the default stays unsymmetric and this is opt-in.**

### Quoting the memory number honestly

The two memory measures move in **opposite** directions with rank count
(`RESULTS_p2h_npsweep.md`), so a bare "−43%" is not meaningful:

- **"Will it fit on my nodes?"** → per-proc peak `INFOG(21)`. Saving **improves** with
  np (−36.6% at np=4 → −52.8% at np=32).
- **"What does the whole job consume?"** → total `INFOG(22)`. Saving **erodes** with
  np (−38.0% → −22.6%).

Always state which measure, and at what np.

⚠ **`INFOG(21)` can badly understate the saving when the factorization is imbalanced.**
At 884 835 DOF it reported only −8.8% while `INFOG(22)` reported a normal −43.2% —
because one rank held **4.16×** the average (`RESULTS_p2h_885k_anomaly.md`). If a
large symmetric job OOMs on a single rank while total memory looks comfortable, that
is the reason, and **adding nodes will not fix it**.

---

## `-BLR <eps>` — almost never

Measured on the cluster from 19.6k to 408.5k DOF:

- **On the unsymmetric path it never once saved peak memory** (+57.1 / +10.9 / +54.0 /
  +16.4%) and was **always slower** (1.21×–3.43×).
- It did **not** rescue an out-of-memory run — at 884.8k DOF on one node it drove the
  node to **323 MB free** without completing, where full-rank had merely failed fast.
- The **one** case where it helped: composed with `-matrixType 2` at **408.5k DOF** it
  bought a further **−18.9%** peak memory for ~1.28× wall. The effect flips sign with
  size (worse at 143.8k, better at 408.5k), consistent with its front-size mechanism.

**So:** reach for `-matrixType 2` first. Only if the model still does not fit, and only
at **≥ ~400k DOF**, and only where an approximate factorization is legal, add
`-BLR 1e-8` on top. Never use it on the unsymmetric path to save memory.

---

## Rank count — more is NOT better

`RESULTS_p2h_npsweep.md`, 143.8k DOF on one 32-core node:

| np | full-rank wall | symmetric wall |
|---|---|---|
| 4 | 39.5 s | 23.3 s |
| **8** | **32.7 s** | **17.9 s** |
| 16 | 45.5 s | 18.1 s |
| 32 | 65.9 s | 21.4 s |

**np=32 is 2.0× SLOWER than np=8** on the same node and model. MUMPS stops scaling
around 8–16 ranks per node here; past that you pay memory bandwidth and parallel
overhead for nothing. Size the rank count to the model, and **measure** rather than
maxing out cores. (The optimum is model- and box-specific; the *shape* is not.)

Multi-node has its own cost: at 408.5k DOF a 2-node np=32 run was several times slower
than 1-node np=16. Use extra nodes when you need the **memory**, not for speed.

---

## `-stats` — how to check any of this yourself

`system Mumps ... -stats` dumps `INFOG(9)` (factor entries), **`INFOG(21)` (peak
MB/proc — the number that decides whether a model fits)**, `INFOG(22)` (total MB), and
`RINFOG(3)` (elimination flops), once per factorization.

⚠ This only works from **Tcl** since ADR-75 **P2h**. Before that the Tcl `system Mumps`
ladder silently discarded `-BLR` and `-stats`, so an older binary gives you a full-rank
solve and no stats while *looking* like it honoured the flags. If you see no `INFOG`
lines, your binary predates P2h.

---

## apeGmsh

apeGmsh #864's typed emitter already supports the full surface —
`ops.system.Mumps(icntl14=..., matrix_type="symmetric", blr=..., stats=True)` and
`ops.system.Pardiso(...)` — and emits them into Tcl decks.

**Its default is `matrix_type="unsymmetric"`**, which is the safe default and should
stay that way. But it means **every generated deck leaves the 1.7–2.3× on the table
unless the author opts in.** If your model is in the "safe" category above, set
`matrix_type="symmetric"` explicitly when you build the job.
