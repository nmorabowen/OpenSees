# ADR-75 P1c — the desktop-PARDISO win COMPOUNDS with model size (+ a capability wall)

Run 2026-07-24, `nmb` desktop, worktree `dist/bin/opensees.pyd`, **MKL_NUM_THREADS=4**.
Cubic LadrunoBrick + LadrunoJ2, 3 LoadControl steps, Newton, `numberer RCM`.
Motivated by the **production-regime correction**: this fork's real workload is *huge solid
nonlinear models*, so every earlier ADR-75 gate (set at 11.5k DOF) was measured two-plus orders of
magnitude below production.

## Results

| n | elements | nDOF | UmfPack | PARDISO | speedup | rel err |
|---|---|---|---|---|---|---|
| 15 | 3,375 | 11,520 | 2.535 s | 1.578 s | 1.61× | **0.0** |
| 20 | 8,000 | 26,460 | 9.267 s | 4.304 s | 2.15× | **0.0** |
| 25 | 15,625 | 50,700 | 38.588 s | 11.363 s | **3.40×** | **0.0** |
| 30 | 27,000 | 86,490 | ❌ **out of memory** | **30.365 s** | **— (∞)** | — |
| 35 | 42,875 | 136,080 | ❌ (not attempted) | **68.586 s** | — | — |

(n = 15–25 are median-of-3; n = 30/35 are single-sample capability probes.)

## Finding 1 — the speedup COMPOUNDS: 1.61× → 2.15× → 3.40×

This is structural, not a constant factor. Fitting the measured points:

| solver | empirical scaling |
|---|---|
| UmfPack | ~**O(N^2.0)** (time ×3.66 then ×4.16 for nDOF ×2.30 then ×1.92) |
| PARDISO | ~**O(N^1.45)** over 15→25, drifting to ~**N^1.8** at 30→35 (normal 3D-direct asymptotics) |

PARDISO uses a better fill-reducing ordering (METIS nested dissection) and a more efficient
multifrontal kernel, so **the gap widens as the model grows** — exactly the direction that matters
for production-size solid models. The P1 headline (1.71× at 11.5k DOF) **understated the win for
the real regime by roughly 2×**.

## Finding 2 — the capability wall (the more important result)

**UmfPack ran OUT OF MEMORY at 86,490 DOF** (`numeric analysis returns -1` =
`UMFPACK_ERROR_out_of_memory`) — a hard wall, not slowness. **PARDISO solved the same model in
30.4 s**, and went on to solve **136,080 DOF in 68.6 s**.

So on this desktop, wiring PARDISO did not merely make existing models faster — it **raised the
largest solvable single-machine model by at least 1.6× in DOF** (and the true ceiling is untested
above 136k). For huge solid nonlinear work that is worth more than any speed ratio: *models that
previously required the cluster may now fit on a workstation.*

## Finding 3 — correctness holds everywhere it can be checked

Bit-identical tip displacement to UmfPack (rel err **exactly 0.0**) at every size where UmfPack
survived, at 4 threads. Consistent with the P1 thread sweep: MKL PARDISO's threaded factorization is
deterministic here, so ADR-75 §7's determinism concern remains a Lane-3-only problem.

## Consequences for ADR-75

1. **The solver lane is worth more than the ADR assumed** for this fork's actual workload. The 66%
   `linearSolve` fraction measured at 11.5k DOF is a *floor*; it grows with N.
2. **Memory is the binding constraint at scale — empirically, now.** UmfPack's wall is the first
   hard evidence. This promotes **MUMPS BLR** (`ICNTL35`) and makes the **LP64 nnz < 2³¹** ceiling a
   real, measurable risk rather than a theoretical footnote.
3. **P1b (symmetric `mtype ±2`) gets more attractive**, since it targets memory as well as time —
   but still must be *measured* (`SparseSYM` is 2.10× slower than UmfPack; symmetric ≠ automatically
   better).
4. **Extrapolation is NOT a substitute for measuring at production size.** The exponents suggest
   ~5–7× at a few hundred thousand DOF, but that is a projection. A representative production deck
   (DOF, element/material mix, steps × iterations) should be benched directly.

## Caveats
- n = 30/35 are single-sample; treat their absolute times as indicative, the *capability* result as solid.
- 3 load steps (not 15) — this measures per-solve cost trend, not a full pushover.
- UmfPack's failure is machine-RAM-dependent; the *ordering* of the wall (UmfPack first) is the point.

## Reproduce
```bash
cd Ladruno_files/testbed/perf/phase1 && OPS_SCALE_SIZES=15,20,25,30 MKL_NUM_THREADS=4 python3.12 -S laneB_scaling.py
```
Knobs: `OPS_SCALE_SIZES`, `OPS_SCALE_SOLVERS`, `OPS_SCALE_STEPS`, `OPS_SCALE_REPEATS`, `OPS_PYD`.
Raw: `laneB_scaling_both.json` (n=15–30 both solvers), `laneB_scaling_pardiso_big.json` (n=30/35).
