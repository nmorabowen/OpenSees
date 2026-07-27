---
title: ADR-75 P2h follow-up — the np-sweep at fixed N (closes a confound, refutes a hypothesis)
project: Ladruno
date: 2026-07-26
host: esmeralda, ONE node (32 cores / 60 GB), OpenMPI 4.0.5
binary: OpenSeesMP @ 45a5b55fe
tags: [adr, performance, mumps, symmetric, scaling]
---

# The np-sweep — what `-matrixType 2` actually saves, as a function of rank count

## Why this run exists

`RESULTS_p2h_cluster_blr.md` §2d-bis reported that symmetric's **per-process peak**
memory saving shrank as ranks grew (−43.1% at np=16 → −8.8% at np=32), and proposed
"more ranks ⇒ smaller per-proc saving" as the mechanism. But **N and np moved
together** across every row supporting it, so it was recorded as an observed pattern
with the settling experiment named. This is that experiment.

## Design — one variable

- **N fixed** at 143 811 DOF (28×28×56 `stdBrick`, `ElasticIsotropic`).
- **One node** for every point, so np ∈ {4, 8, 16, 32} all fit in 32 cores and the
  **interconnect never enters** — the 2-node runs elsewhere confound node count with
  rank count, and this deliberately does not.
- `numberer LadrunoParallelRCM` throughout (the numberer the first sweep should have
  used).
- 143.8k chosen for memory headroom: full-rank at np=16 was ~16 GB total, so every
  np fits 60 GB. (408.5k would OOM at np=32 on one node: ~2 GB/proc × 32 = 64 GB.)

Driver: `p2h_npsweep.sbatch`. Deck: `p2h_smoke2.tcl`.

## Results

| np | `full` INFOG21 | `sym` INFOG21 | **Δ per-proc** | `full` INFOG22 | `sym` INFOG22 | **Δ total** | `full` wall ms | `sym` wall ms | **sym speedup** |
|---|---|---|---|---|---|---|---|---|---|
| 4 | 1 002 | 635 | **−36.6%** | 3 755 | 2 327 | **−38.0%** | 39 492 | 23 275 | 1.70× |
| 8 | 738 | 480 | **−35.0%** | 4 461 | 2 857 | **−36.0%** | 32 653 | 17 920 | 1.82× |
| 16 | 1 046 | 571 | **−45.4%** | 4 822 | 3 417 | **−29.1%** | 45 482 | 18 102 | 2.51× |
| 32 | 709 | 335 | **−52.8%** | 6 390 | 4 947 | **−22.6%** | 65 874 | 21 387 | **3.08×** |

## Finding 1 — the hypothesis is REFUTED, and backwards

**Per-proc peak saving does not shrink with rank count; it grows** — −36.6% at np=4
to **−52.8%** at np=32. The proposed mechanism ("the peak is set by the largest
distributed front, so more ranks dilute the benefit") predicts the opposite of what
happens and should not be repeated.

⇒ **The −8.8% at 884.8k DOF is not an np effect.** That run differed from these rows
in **two** ways — 6.2× the DOF **and** 2 nodes rather than 1 — and node count has
never been isolated. It stays unexplained. The only thing now established is what it
is *not*.

## Finding 2 — the two memory measures move in OPPOSITE directions with np

| | np=4 | np=32 |
|---|---|---|
| Δ per-proc peak (`INFOG(21)`) | −36.6% | **−52.8%** (better) |
| Δ total factor (`INFOG(22)`) | −38.0% | **−22.6%** (worse) |

Total factor memory grows with rank count for **both** modes (full 3 755 → 6 390 MB;
sym 2 327 → 4 947 MB) — the usual duplication/overhead of more ranks — and symmetric's
total grows proportionally faster, eroding its advantage. So:

- Asking **"will this model fit on my nodes?"** → quote **per-proc** `INFOG(21)`, and
  more ranks help symmetric.
- Asking **"how much memory does the job consume?"** → quote **total** `INFOG(22)`,
  and more ranks hurt symmetric.

**Never quote a symmetric memory saving without saying which measure and at what np.**

## Finding 3 (unlooked-for, and the most actionable) — MUMPS does NOT scale past ~8 ranks here

Full-rank wall by np: 39.5 s (4) → **32.7 s (8)** → 45.5 s (16) → **65.9 s (32)**.
**np=8 is the optimum and np=32 is 2.0× SLOWER than np=8** on the same node and the
same model. Symmetric is flatter but shows the same shape: 23.3 → **17.9** → 18.1 →
21.4 s.

Throwing ranks at a MUMPS solve on one node is actively counterproductive past ~8–16.
This is one node, so it is memory-bandwidth and MUMPS parallel efficiency, not
network. It also explains why symmetric's *speedup ratio* looks best at np=32
(3.08×): the denominator is degrading faster than the numerator, so **that ratio
flatters symmetric and should not be quoted as a headline** — the honest headline
remains the 1.7–2.5× measured at sane rank counts.

## Caveats

- Single N (143 811 DOF) and a single node type. Finding 3's optimum (np≈8) is
  specific to this model and box; the *shape* is likely general, the number is not.
- Elastic material (`ElasticIsotropic`): this measures the **solver**, which is the
  intent — element cost would only add a constant to both modes.
- n=1 per cell. The phase2 noise floor was measured at ~2.7% on a repeated
  configuration, so the ≥20% effects above are safe; the np=4-vs-np=8 per-proc
  difference (−36.6% vs −35.0%) is **not** separable from noise.
