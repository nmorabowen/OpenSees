# ADR-75 P2 — MUMPS `-BLR`: wired, verified engaging, and **slower at desktop scale**

Run 2026-07-24, `nmb` desktop, `dist/openseesmp/openseesmp.pyd`, `mpiexec -n 2`.
Partitioned 3D LadrunoBrick + LadrunoJ2 (z-slab per rank), Newton, 2 LoadControl steps,
`numberer LadrunoParallelRCM`, `system Mumps -ICNTL14 200 [-BLR eps]`.

## Small model — 8×8×16 (1,024 ele, ~4k DOF)

| config | wall | rel diff vs full-rank |
|---|---|---|
| full-rank | 0.618 s | — |
| `-BLR 1e-9` | 0.546 s | 1.202e-16 |
| `-BLR 1e-4` | 0.527 s | **1.202e-16** |

**A loose tolerance producing the SAME machine-epsilon difference as a tight one means BLR never
engaged** — the fronts are too small to compress. This model cannot test BLR at all; recording it
because "it ran and the answer was right" would otherwise look like a passing test.

## Larger model — 20×20×24 (9,600 ele, ~30k DOF) — BLR genuinely engages

| config | wall | rel diff vs full-rank |
|---|---|---|
| full-rank | **6.419 s** | — |
| `-BLR 1e-9` | 10.936 s (**1.70× slower**) | 0.000e+00 |
| `-BLR 1e-4` | 20.362 s (**3.17× slower**) | **1.864e-12** |

The **non-zero** deviation at `eps=1e-4` is the proof of engagement: BLR is an approximate
factorization and is now actually dropping rank. At `1e-9` the compression is tight enough to leave
the answer bit-identical.

## Finding — BLR is a MEMORY lever, not a speed lever, and it can COST time

At this scale BLR is **slower**, and — counterintuitively — **the looser tolerance is slower still**.
Two compounding causes:

1. **Compression overhead exceeds flop savings on small fronts.** BLR must compute low-rank
   approximations of each frontal matrix; that only pays once fronts are large enough that the
   avoided flops dominate. ~30k DOF on 2 ranks is far below that crossover.
2. **In a NONLINEAR loop, an inexact solve degrades Newton.** A looser BLR tolerance returns a less
   accurate correction, so `NormDispIncr 1e-7` needs more iterations — i.e. *more solves*. This is
   why 1e-4 is worse than 1e-9, the opposite of the "looser = more compression = faster" intuition
   that holds for a single linear solve. **This interaction is specific to nonlinear analysis and is
   the main trap in using BLR here.**

**Consequence for ADR-75:** BLR's justification is **factor memory at large scale** (the ADR-75 P1c
capability wall), *not* time. Do not enable it expecting a speedup, and expect to re-tune the Newton
test if you do. It stays opt-in and off by default.

## P2b — MEASURED memory, via the new `-stats` (INFOG/RINFOG dump)

`system Mumps ... -stats` now prints MUMPS's own factorization telemetry, so BLR's effect is
finally observable. Same 20×20×24 model (n = 31,752), `mpiexec -n 2`:

| config | factor entries `INFOG(9)` | **MB/proc `INFOG(21)`** | MB total `INFOG(22)` | elim flops `RINFOG(3)` | BLR `RINFOG(14)` |
|---|---|---|---|---|---|
| full-rank | 4.2369e7 | **263** | 498 | 5.0251e10 | — |
| `-BLR 1e-9` | 4.2273e7 (**−0.2%**) | **285 (+8.4%)** | 553 (+11%) | 5.0251e10 | 5.2778e10 |
| `-BLR 1e-4` | 3.3129e7 (**−21.8%**) | **251 (−4.6%)** | 478 (−4.0%) | 5.0251e10 | 2.7786e10 (**−45%**) |

### The non-obvious result: BLR shrinks the FACTORS but barely moves PEAK MEMORY

- **At `eps=1e-9` BLR is strictly worse on every axis**: factor entries essentially unchanged
  (−0.2% — the tolerance is too tight to drop meaningful rank), peak memory **+8.4% per proc** from
  BLR's own bookkeeping, no flop saving, and 1.70× slower.
- **At `eps=1e-4` BLR really compresses**: factor entries **−21.8%** and BLR flops `RINFOG(14)`
  **45% below** the full-rank elimination flops. **Yet peak memory per proc falls only 4.6%** — and
  the run was 3.17× *slower*.

**Why:** peak memory during factorization is dominated by the **active frontal/working space**, not
by the stored factors. BLR compresses what it *stores*; it does not shrink the working set by the
same proportion. So "BLR saves memory" is true of factor storage and largely false of the peak
allocation that actually decides whether a model fits.

**Consequence:** at ~32k DOF / np2 BLR is **not a win on any axis** except factor-entry count. It did
not move the constraint that matters (peak memory) and it cost 1.7–3.2× wall time. The regime where
BLR is designed to pay — fronts large enough that compression dominates its overhead — is **above
what this desktop can reach**, so this does *not* refute BLR at cluster scale; it does refute
enabling it by default or expecting relief at these sizes.

## Status: what is and is not verified

✅ **Verified here:** the `-BLR <eps>` / `-ICNTL35` / `-CNTL7` options parse; `ICNTL(35)`/`CNTL(7)`
reach MUMPS at both the analysis and factorization phases; the controls propagate to subordinate
ranks via `sendSelf`/`recvSelf` (without that, rank 0 would factor BLR while the others factored
full-rank — an inconsistent distributed factorization); results are correct; BLR demonstrably
engages above a size threshold; default-off is byte-identical to the pre-BLR solver.

✅ **Now also verified (P2b):** `-stats` surfaces `INFOG(9)/(21)/(22)`, `RINFOG(3)` and the BLR
`RINFOG(14)/(15)`, so compression is directly observable; and the memory question is *answered at
this scale* — BLR shrinks factor entries by up to 21.8% but peak memory per proc by only 4.6%.

❌ **NOT verified — needs the cluster:** whether BLR crosses over into a genuine win on
**production-size fronts**. Everything here is ~32k DOF on 2 ranks, far below BLR's design regime;
the measured losses bound the *small* end only. A production deck should be run with `-stats` and
BLR on/off, comparing `INFOG(21)` (the number that decides whether a model fits) and wall time.

## Reproduce
```bash
BLR_NX=20 BLR_NY=20 BLR_NZ=24 dist/openseesmp/mpiexec.exe -n 2 \
  python3.12 -S Ladruno_files/testbed/perf/phase1/mp_blr_smoke.py dist/openseesmp
```
