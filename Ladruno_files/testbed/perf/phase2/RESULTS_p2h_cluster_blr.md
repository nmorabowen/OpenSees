---
title: ADR-75 P2h — MUMPS BLR vs full-rank vs symmetric half-storage, ON THE CLUSTER
project: Ladruno
date: 2026-07-26
host: esmeralda (SLURM, node[1-18], 32 cores / 60 GB per node, OpenMPI 4.0.5 /opt/openmpi)
binary: OpenSeesMP @ claude/adr75-continue (P2h), /mnt/nfshare MUMPS + ScaLAPACK
tags: [adr, performance, mumps, blr, cluster]
---

# ADR-75 items 1 + 2 — measured on esmeralda

Closes the two cluster-gated items the ADR-75 handoff had been carrying:
**item 1** (is BLR the memory lever at cluster scale?) and **item 2** (does P1d's
symmetric half-storage win transfer from desktop PARDISO to MUMPS?).

## 0. The measurement was impossible until P2h — and would have LIED

`OpenSeesMP` is Tcl-driven. The Tcl `system Mumps` ladder understood only
`-ICNTL14`/`-ICNTL7`/`-matrixType`; `-BLR` and `-stats` were wired into the
**Python** ladder only, and the Tcl ladder's terminal `else currentArg++` arm
discarded unknown tokens **silently**. The handoff's prescribed command would
have run full-rank twice with no stats and reported "`INFOG(21)` did not move".

This was **not** hypothetical: `apeGmsh` #864 (`feat(system): PARDISO solver +
full MUMPS option surface`) emits exactly `-ICNTL14 / -ICNTL7 / -matrixType /
-BLR / -ICNTL35 / -CNTL7 / -commSplit / -stats` **into Tcl decks**
(`src/apeGmsh/opensees/analysis/system.py`). The deck generator was already
producing options the solver silently ignored.

### Parse-path proof (execution, on the cluster binary)

| deck | result |
|---|---|
| `system Mumps` (bare) | runs, no stats — pre-existing behaviour preserved |
| `system Mumps -ICNTL14 200 -stats` | **prints INFOG(9)/(21)/(22), RINFOG(3)** — never worked in Tcl before |
| `system Mumps ... -BLRR 1e-8 -stats` (typo) | `Mumps Warning: unrecognized option '-BLRR' -- ignored` **+** `'1e-8' -- ignored` |
| `system Mumps -ICNTL14 200 -stats -ICNTL7` (trailing) | `Mumps Warning: -ICNTL7 needs a value -- ignored` (was a read past the end of `argv`) |

**BLR provably engages** (this is the check that a bit-identical answer at loose
eps would have failed): at `-BLR 1e-4` the tip displacement moves from
`1.17964412133551431e+01` to `1.17964412133573440e+01` (~1.9e-13 relative) and
`INFOG(9)` drops 12%.

## 1. Model

Parametric 3D solid cantilever, `stdBrick` + `ElasticIsotropic`, k-slab
partitioned one slab per rank, `numberer ParallelRCM`, `constraints Plain`,
2 × `LoadControl` steps under full Newton (`NormDispIncr 1e-7`).
Deck: `p2h_smoke.tcl` (in this directory). Sweep drivers: `p2h_sweep.sbatch`
(1 node), `p2h_big2.sbatch` (2 nodes).

`INFOG(21)` = **peak factor MB per process, max over procs** — the number that
decides whether a model fits. Reported per solve; the max across the run is quoted.

## 2. Results

### 2a. ~19.6k DOF, np=2 (the size class the desktop verdict #626 covered)

| mode | INFOG(21) MB/proc | wall ms | vs full |
|---|---|---|---|
| `full` (unsym) | 126 | 1 972 | — |
| `-BLR 1e-8` | **198** | 3 324 | **+57.1% mem, 1.69× slower** |
| `-BLR 1e-4` | 187 | 6 762 | +48.4% mem, **3.43× slower** |
| **`-matrixType 2`** | **79** | **1 538** | **−37.3% mem, 1.28× faster** |
| `-matrixType 2 -BLR 1e-8` | 122 | 1 804 | worse than symmetric alone |

**Cross-validates the desktop P2b run almost exactly** on the time axis — #626
measured 1.70× (eps 1e-9) and 3.17× (1e-4) slower; here 1.69× and 3.43×. Two
different machines, two different MUMPS builds, same verdict.

### 2b. ~19.6k DOF, np=16 — same model, more ranks

| mode | INFOG(21) MB/proc | wall ms | vs full |
|---|---|---|---|
| `full` | 110 | 1 248 | — |
| `-BLR 1e-8` | 122 | 2 139 | +10.9% mem, 1.71× slower |
| **`-matrixType 2`** | **99** | **720** | −10.0% mem, **1.73× faster** |
| sym + BLR | 99 | 1 157 | — |

Note the **memory advantage of symmetric shrinks with rank count at fixed N**
(−37.3% at np=2 → −10.0% at np=16): with 16 ranks the per-proc peak is no
longer dominated by this model's factor. That effect reverses once N grows (2c).

### 2c. ~143.8k DOF, np=16 — the decision-grade point

| mode | INFOG(21) MB/proc | wall ms | vs full |
|---|---|---|---|
| `full` | 993 | 47 805 | — |
| `-BLR 1e-8` | **1 529** | 58 064 | **+54.0% mem, 1.21× slower** |
| **`-matrixType 2`** | **571** | **20 581** | **−42.5% mem, 2.32× faster** |
| sym + BLR | 657 | 29 924 | −33.8% mem, 1.60× faster — still worse than sym alone |

### 2d. ~408.5k DOF, np=16 — past every "small end" caveat

| mode | INFOG(21) MB/proc | wall ms | vs full |
|---|---|---|---|
| `full` | 3 020 | 179 661 | — |
| `-BLR 1e-8` | **3 515** | 231 890 | **+16.4% mem, 1.29× slower** |
| **`-matrixType 2`** | **1 719** | **103 740** | **−43.1% mem, 1.73× faster** |
| `-matrixType 2 -BLR 1e-8` | **1 395** | 132 918 | **−53.8% mem**, 1.35× faster |

At 48 GB of the node's 60 GB, full-rank *fit*; BLR's +16.4% is the direction that
puts a model **over** a node, not under it.

> **⚠ The one result that cuts the other way — do not skip it.** At this size, and
> only at this size, **BLR composed with symmetric storage finally SAVED memory**:
> 1 395 vs symmetric-alone's 1 719 MB/proc, **−18.9%**, for a 1.28× time cost.
> The sym+BLR-vs-sym comparison changes sign with N:
>
> | N (DOF) | sym | sym+BLR | Δ |
> |---|---|---|---|
> | 19.6k | 99 | 99 | 0.0% |
> | 143.8k | 571 | 657 | **+15.1%** (worse) |
> | 408.5k | 1 719 | 1 395 | **−18.9%** (better) |
>
> This is precisely the front-size-driven sign flip the BLR hypothesis predicted —
> it just needs **half-storage plus a large front** before it appears, and it does
> **not** appear on the unsymmetric path at any size tested. Any claim that "BLR
> never pays" is therefore **too strong**; the defensible claim is the narrower one
> in §3.

### 2e. CONFOUND CONTROL — the numberer (a miss, caught in review, then measured)

The sweep above used vanilla **`numberer ParallelRCM`**. It should have used
**`LadrunoParallelRCM`**: the fork's own `phase1/mp_blr_smoke.py` prefers it with a
fallback, ADR-74 exists precisely because the vanilla parallel numberer degrades at
scale, and it **is** wired into the Tcl ladder (`commands.cpp:4387`) — so it was
available and simply not used.

Rather than argue about whether that matters, it was **run**. Same model
(143.8k DOF), same binary, np=16, `step1` = first `analyze` (carries numbering +
`setSize`), `step2` = a pure second step:

| numberer | mode | wall ms | step1 ms | step2 ms | setup ≈ step1−step2 |
|---|---|---|---|---|---|
| `ParallelRCM` | full | 46 518 | 24 139 | 22 379 | **1 760** |
| `ParallelRCM` | sym | 18 892 | 10 289 | 8 603 | **1 686** |
| `LadrunoParallelRCM` | full | 45 487 | 23 136 | 22 351 | **785** |
| `LadrunoParallelRCM` | sym | 17 845 | 9 170 | 8 675 | **495** |

- **`LadrunoParallelRCM` cuts setup 2.2×–3.4×** (1 760→785 ms full; 1 686→495 ms
  sym) — ADR-74's claim reproduces.
- **It changes no verdict here.** At 143.8k the solve dominates, so setup is only
  1.7–3.8% of wall and the whole-run gain is 2.2% (full) / 5.5% (sym).
- **The sym-vs-full ratio moves 2.46× → 2.55×**, i.e. the vanilla numberer was
  *diluting* the symmetric win. Every wall ratio in §2 is therefore **conservative**,
  never inflated — which is the direction that cannot manufacture a false positive.
- Reproducibility bonus: this re-run's `full` (46 518 ms) vs the original sweep's
  (47 805 ms) differ by **2.7%**, a useful read on the noise floor.

⚠ **Not extrapolated.** ADR-74's degradation concern is about *much* larger decks
(its kill-run is 18.6M DOF). Nothing here bounds the numbering share at 885k+; the
sweep's absolute walls at the top end may carry more setup than these percentages
suggest. Use `LadrunoParallelRCM` in any follow-up.

## 3. Verdicts

### Item 1 — BLR is not the memory lever; but the claim has to be narrower than it first looked

**The precise finding, in three parts:**

1. **On the UNSYMMETRIC path, `-BLR` never saved peak memory at any size or rank
   count tested, and was always slower.** Four configurations, sign invariant.
2. **On the SYMMETRIC path, BLR's memory effect flips sign with N** — worse at
   143.8k (+15.1%), **better at 408.5k (−18.9% vs symmetric alone)**. This is the
   front-size mechanism finally showing up, and it means "BLR never pays" would be
   an **overclaim**.
3. **BLR was never a time win** in any configuration: 1.21×–3.43× slower, always.

So BLR is not a drop-in memory lever and must not be the first thing reached for —
but it is not dead either, and the ADR's front-size reasoning is vindicated in the
one corner where fronts are large *and* storage is already halved.

`-BLR` lost on **both** axes at **every** size and rank count measured:

`-BLR 1e-8` lost on **both** axes at **every** size and rank count measured — the
sign never flips:

| N (DOF) | np | BLR Δ peak mem | BLR wall ratio |
|---|---|---|---|
| 19.6k | 2 | **+57.1%** | 1.69× slower |
| 19.6k | 16 | **+10.9%** | 1.71× slower |
| 143.8k | 16 | **+54.0%** | 1.21× slower |
| **408.5k** | 16 | **+16.4%** | **1.29× slower** |

The handoff's escape hatch was that #626 "bounds the small end only", since BLR's
payoff is governed by **front size** (~N^(2/3) DOF, front memory ~N^(4/3)). That
hypothesis is now tested across a **21× DOF range up to 408.5k** — past the 136k
that ADR-75 §1 sets as its non-representative threshold — and it is **not
supported**. The magnitude of the memory penalty is noisy (+10.9 / +54.0 / +16.4%)
but its **sign is invariant: BLR never once saved peak memory.**

Front memory itself behaved as predicted (110 → 993 → 3 020 MB/proc), so the
*premise* is right and the *conclusion* still fails: **BLR compresses stored
factors (`INFOG(9)`) but adds working space, and `INFOG(21)` is dominated by the
working space.** At 408.5k the +16.4% is the difference between fitting a node
and not.

⚠ **Scope, stated honestly.** This bounds BLR through ~408k DOF. It does **not**
prove BLR is useless at 1M–10M DOF, where fronts reach ~1.7 GB–36 GB and the
compressible fraction of a front grows. An 885k-DOF point was attempted; see §5.

**Recommendation — an ORDER, not a yes/no:**

1. **Reach for `-matrixType 2` first.** It is exact, it is faster, and it is the
   only lever that saved memory at every size (§Item 2).
2. **If the model still does not fit, then add `-BLR` on top of it** — but only at
   **≥ ~400k DOF**, where it bought a further **−18.9%**, and only where an
   approximate factorization is legal (never on an oracle path). Expect to pay
   ~1.3× wall for it.
3. **Do not use `-BLR` on the unsymmetric path** to save memory. It did not, once,
   at any size here.

The next candidate for a genuine memory wall remains **MUMPS out-of-core
`ICNTL(22)`** (trades memory for I/O, not accuracy), still unexposed by the fork.

### Item 2 — `-matrixType 2` IS the lever, and it transfers from PARDISO to MUMPS

`-matrixType 2` won on **both** axes at **every** size and rank count:

| N (DOF) | np | sym Δ peak mem | sym wall ratio |
|---|---|---|---|
| 19.6k | 2 | −37.3% | 1.28× faster |
| 19.6k | 16 | −10.0% | 1.73× faster |
| 143.8k | 16 | **−42.5%** | **2.32× faster** |
| **408.5k** | 16 | **−43.1%** | **1.73× faster** |

The memory saving **converges on ≈ −43%** once the model is large enough to
dominate the per-proc peak, landing almost exactly on P1d's desktop PARDISO
measurement of **−41.8%** — an independent confirmation across a different
solver, machine, and parallel model. Unlike BLR it is **exact** (answers agree to
~1e-16), so it is legal on oracle paths.

**Recommendation:** `-matrixType 2` should be the **default for symmetric-tangent
cluster decks**, exactly as on the desktop — with the same standing caveat that
half-storage silently solves the wrong system on a genuinely unsymmetric tangent
(contact / non-associated / follower loads / `LadrunoUP`), so it stays opt-in.

## 4. Gotchas banked here

- **`mpirun` works inside a 1-node SLURM allocation and dies instantly across 2
  nodes** — `FORCE-TERMINATE ... plm_slurm_module.c(471)`, internal ORTE error.
  Use `srun --cpu-bind=cores --mpi=pmix_v3` for multi-node, per
  `02_esmeralda_linux_build_guide.md` §7. The sweep script still printed its
  final `..._DONE` banner and the job exited 0, so **rc proved nothing** — only
  the absent `P2H_RESULT` artifacts revealed that all 6 runs had died. This is
  the banked "assert on artifacts, not rc" rule earning its keep again.
- **Do not rebuild the binary while a sweep is running.** Each mode re-execs
  `openseesmp.sh`, so replacing the file mid-job silently mixes binaries across
  the comparison.
- **SLURM is at `/opt/slurm/bin`, not on `PATH`.** `which sbatch` returning
  nothing is *not* evidence the cluster is down — which is how it came to be
  recorded as down for a whole session while it had 33 days of uptime.
