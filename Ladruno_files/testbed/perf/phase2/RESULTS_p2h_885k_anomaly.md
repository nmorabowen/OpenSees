---
title: ADR-75 P2h — closing the 884.8k per-proc anomaly (it is LOAD IMBALANCE, not np, nodes, or N)
project: Ladruno
date: 2026-07-26
host: esmeralda, np=32 throughout, 1 node (32 cores) and 2 nodes (16+16)
binary: OpenSeesMP @ 32bc9c115
tags: [adr, performance, mumps, symmetric, load-balance]
---

# The 884.8k anomaly — closed

## The loose end

`-matrixType 2` cut the per-process peak (`INFOG(21)`) by 43–53% everywhere **except**
at 884 835 DOF / np=32 / 2 nodes, where it managed only **−8.8%** — while cutting
*total* factor memory (`INFOG(22)`) by a perfectly normal −43.2%.

Two hypotheses were floated and **both were wrong**:

1. *"Per-proc saving shrinks with rank count."* **Refuted** by the fixed-N np-sweep
   (`RESULTS_p2h_npsweep.md`): the saving **grows** with np, −36.6% (np=4) → −52.8%
   (np=32).
2. *"Then it must be node count or N"* — the two variables still confounded in the
   original observation (884.8k differed from the np-sweep rows by **6.2× the DOF**
   *and* **2 nodes instead of 1**). Measured below. **Also not the cause.**

## The isolation runs

np=32 fixed throughout, `LadrunoParallelRCM`, same binary and deck.
`max/avg` = `INFOG(21) / (INFOG(22)/np)` — a **direct measure of load imbalance**,
since MUMPS reports `INFOG(21)` as the max over procs and `INFOG(22)` as the sum.

| N (DOF) | nodes | full max | full max/avg | sym max | **sym max/avg** | Δ per-proc | Δ total |
|---|---|---|---|---|---|---|---|
| 143 811 | 1 | 709 | 3.55 | 335 | 2.17 | −52.8% | −22.6% |
| 143 811 | 2 | 779 | 3.83 | 392 | 2.51 | −49.7% | −23.1% |
| 408 483 | 1 | 1 857 | 2.73 | 1 017 | 2.69 | −45.2% | −44.5% |
| 408 483 | 2 | 1 496 | 2.20 | 902 | 2.30 | −39.7% | −42.2% |
| **884 835** | 2 | 4 382 | 2.59 | **3 996** | **4.16** | **−8.8%** | **−43.2%** |

### Node count is not the cause

At fixed N and fixed np, going 1 → 2 nodes barely moves the saving:
**−52.8% → −49.7%** at 143 811 DOF, **−45.2% → −39.7%** at 408 483 DOF. Nothing like
a collapse to −8.8%.

### N is not the cause either

The per-proc saving does decline gently with N (−52.8% → −45.2% at fixed 1 node), but
**408 483 DOF still returns −45.2%** — squarely in the normal band. There is no
gradual slide toward −8.8%; the 884 835 row is a discontinuity, not the end of a trend.

## What it actually is: **symmetric's load imbalance inverts at 884.8k**

Read the `max/avg` columns:

- At 143 811 DOF, symmetric is **markedly better balanced** than full-rank
  (2.17 vs 3.55; 2.51 vs 3.83).
- At 408 483 DOF the two are **level** (2.69 vs 2.73; 2.30 vs 2.20).
- At 884 835 DOF the ordering **inverts and the gap is large**: symmetric **4.16** vs
  full-rank **2.59**.

So symmetric is still halving the factor **globally** — `INFOG(22)` is −43.2%, exactly
in line with the −42.2%/−44.5% measured at 408 483 DOF. The saving has not gone
anywhere. It has **concentrated**: one rank now holds ~4.2× the average, so
`INFOG(21)`, which is a **max over procs**, stops tracking the average and reports
−8.8%.

**The anomaly is a property of the metric under imbalance, not a property of
symmetric storage.**

### Mechanism — measured vs hypothesised, stated separately

- **Measured:** the `max/avg` imbalance ratio for symmetric jumps from ~2.2–2.7 to
  **4.16** between 408 483 and 884 835 DOF, while full-rank's stays flat at ~2.2–2.7.
  That is a direct arithmetic consequence of the reported `INFOG(21)`/`INFOG(22)`.
- **Hypothesised, NOT proven:** the concentration is most likely the top separator's
  dense frontal block landing on one rank, whose working space does not halve the way
  distributed factor storage does. Confirming that needs per-rank `INFOG(16)`/`INFOG(18)`
  or MUMPS's own tree diagnostics, which the fork's `-stats` does not currently expose.
  **Do not repeat the mechanism as established** — only the imbalance is.

## Consequences

1. **The `-matrixType 2` recommendation is unaffected.** It halves the factor at every
   size tested (`INFOG(22)` −42% to −45% for N ≥ 408k), is exact, and is 1.7–2.5×
   faster. Nothing here weakens it.
2. **Quote `INFOG(22)` when you mean "how much memory does symmetric save".**
   `INFOG(21)` is a max and is hostage to balance; it answers a *different* question
   ("will the worst rank fit?") and can understate the saving badly, as it does here.
3. **A memory-limited run at ~1M DOF may be limited by ONE rank, not by the average.**
   That is the practically useful finding: raising node count does not help if a single
   rank holds 4× the mean. If a large symmetric job OOMs on one rank while total memory
   looks comfortable, this is why.
4. **Open, and now well-posed:** does symmetric's imbalance keep growing past 884.8k,
   and is it fixable by MUMPS ordering (`-ICNTL7`) or a different partition? Worth one
   run *if* a production model actually hits the wall — not before.

## Caveats

- n=1 per cell; the ~2.7% noise floor is far below the effects discussed, but the
  408 483-vs-884 835 imbalance jump rests on **one** run at 884 835 DOF.
- `max/avg` is a coarse balance proxy. It is exactly the ratio of two numbers MUMPS
  reports, which is its virtue; it says *that* work concentrated, not *where*.
- 884 835 DOF was only ever measured on 2 nodes — it does not fit one (60 GB), so the
  bottom row cannot be node-count-controlled the way the others are.
