# ADR-75 P1 — desktop PARDISO: MEASURED, gate PASSED

Run 2026-07-24 on the `nmb` desktop against the **newly built**
`.claude/worktrees/mumps-opensees-study-f833bf/dist/bin/opensees.pyd` (the first serial OpenSees
binary in this fork with a threaded sparse-direct solver: 477 pardiso symbols, was 0). Lane-B `15³`
LadrunoBrick + LadrunoJ2, ~11.5k DOF, 15 steps, Newton, `numberer RCM`. Median of 5 interleaved
rounds per thread count; UmfPack re-measured **inside the same binary** so the comparison is
apples-to-apples.

## Results

| MKL threads | UmfPack (s) | **PARDISO (s)** | PARDISO vs UmfPack | tip `ux` rel err |
|---|---|---|---|---|
| 1 | 18.546 | **15.611** | **1.19×** | **0.0** |
| 2 | 18.291 | **12.230** | **1.50×** | **0.0** |
| 4 | 17.775 | **10.396** | **1.71×** | **0.0** |
| 8 | 17.330 | **9.870** | **1.76×** | **0.0** |

PARDISO thread scaling against its own 1-thread time: **1.00× / 1.28× / 1.50× / 1.58×**.

## Verdict

**The P1 gate is PASSED.** The gate was "beat UmfPack's locked 22.711 s on Lane B at 4 threads":
PARDISO at 4 threads is **10.396 s**. Within the same binary the honest figure is **1.71× faster
than UmfPack at 4 threads** (1.76× at 8).

Three findings worth keeping:

1. **PARDISO wins even single-threaded (1.19×).** Part of the gain is a better solver, not just
   threading — so the desktop win does not depend entirely on core count.

2. **Correctness is exact at every thread count — rel err `0.0`, bit-identical to UmfPack.**
   ~~Threading introduced **no** floating-point drift here: MKL PARDISO's threaded factorization is
   deterministic for this problem. That matters for the fork's byte-identical/1e-12 discipline —
   the ADR-75 §7 determinism concern is a *Lane-3* (assembly-reduction) problem, and does **not**
   apply to this solver. No ordered-reduction flag is needed for PARDISO.~~

   > **⚠ CORRECTED 2026-07-25 (ADR-75 P1f).** The accuracy statement stands — every configuration
   > agrees with UmfPack to `0.0`. **The determinism conclusion does not.** It rested on ONE run per
   > thread count, which cannot distinguish "deterministic" from "the same value came up twice".
   > Re-measured with **10 runs of the same binary at a fixed thread count**: at 4 threads a 14³
   > Lane-B model returns **two distinct tip displacements in a 5/5 split** (`…834687` vs `…834732`,
   > ~1 ULP); at 1 thread it is 10/10 identical. Verified on the *pre-P1f* binary too, so it is a
   > property of threaded MKL PARDISO, not of any fork change. It is size-dependent — an 8³ model is
   > reproducible even at 4 threads — which is why a single-sample check at 15³ passed.
   >
   > **Consequence:** threaded `system Pardiso` is **not** byte-reproducible run-to-run, so the
   > ADR-75 §7 determinism concern is *not* Lane-3-only. Any byte-identical CI gate or oracle
   > comparison must pin `MKL_NUM_THREADS=1` (reproducible, and P1e/P1f measurements confirm the
   > 1-thread path is stable) — or use PARDISO's own CNR control, `iparm[33]`, which is available to
   > this fork because we set `iparm[1]=2` (serial METIS); Intel only forbids CNR when `iparm[1]=3`.
   > `iparm[33]` is **not** currently exposed as an option. Drift is last-bit, so this is a
   > *reproducibility* problem, not an accuracy one.

3. **Scaling flattens past 4 threads (1.50× → 1.58×), exactly as predicted.** Sparse factorization
   is memory-bandwidth-bound, not compute-bound; 4 threads captures ~95% of the available win.
   **Recommend 4 threads as the practical desktop default.**

### Amdahl cross-check (the numbers are self-consistent)
ADR-75 predicted a ~2.2× ceiling on *total* runtime because PARDISO threads only the factor/solve
(~66% of Lane B per ADR-40b), leaving assembly and element state determination single-threaded. The
measured 1.76× sits just under that ceiling, and PARDISO's own 1.58× internal scaling accounts for
the gap. **The remaining ~34% is untouchable by any solver — it is ADR-75 Lane 3 (threaded element
assembly) territory.**

### Side effect worth noting
UmfPack also drifts faster with more threads (18.55 → 17.33 s). UMFPACK is single-threaded at the
algorithm level but calls BLAS for its dense frontal kernels, and this build now links **threaded**
MKL — so those kernels pick up some threading for free. Small (~7%) but real, and a bonus of the
threaded-MKL link.

## Reproduce

```bash
cd Ladruno_files/testbed/perf/phase1 && ./sweep_p1.sh
```

Raw data: `laneB_p1_threads{1,2,4,8}.json` (all samples retained). Harness knobs:
`OPS_BENCH_SOLVERS`, `OPS_BENCH_REPEATS`, `OPS_PYD`, `MKL_NUM_THREADS`.
