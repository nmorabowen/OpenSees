---
title: "WP-232g — M0b measurements (ADR 232 phase 2)"
project: Ladruno
tags: [mpm, wp-232, m0b, measurements]
---

# WP-232g — M0b measurements

ADR 232 §9 phase 2 ([[232_ladruno_mpm_general_driver_adr]]): the per-step rebuild of mode A against one Newton iteration, explicit (mode F) throughput across materials, memory per particle-element, and the 232f data owed on the cantilever. The definitions and the pre-registered decision table are in `.claude/skills/ladruno-mpm/references/wp-briefs.md`, section "232g". This document records the numbers. The section "Applying the pre-registered table" is the worker's **proposal**; the orchestrator's application is the decision.

Data: `Ladruno_scripts/mpm_m0b/results/full_20261115/` (the aggregates `*/m*_results.json` and `index.json` are committed; per-run files, profiler HDF5 and child logs are not). Harness: `Ladruno_scripts/mpm_m0b/` (its README states every model choice and its reason). Every number below is derived from those aggregates. No run was repeated after the fact, except the M2 driver legs, which the orchestrator re-ran (see M2).

## Method

- **Instrumentation.** Profiler scopes only, in build `9db5803f6` (PR #1089). `OPS_PROFILE_SCOPE` is a no-op when the profiler is off (`ProfilerMacros.h:77`).
  - `StaticAnalysis` gets a `domainChanged` scope and `dc.*` sub-scopes.
  - `DirectIntegrationAnalysis` gets `dc.clearAll`, `dc.integrator` and `dc.algorithm`.
  - Mode A gets `mpm.a.step`, then `mpm.a.prepare` with six children, and `mpm.a.g2p`.
  - Mode F gets the deep `mpm.f.material`: one `ElemScope` per material call.
- **One fresh process per run.**
  - M1, M3 mode A and M4: `MKL_NUM_THREADS` and `OMP_NUM_THREADS` = 8.
  - M2 and M3 mode F: 1, plus one elastic row with both unset ("default").
  - Every child pins its process affinity to the P-cores (`0xFFFF`) before `opensees` loads. Each run records the mask the OS reports back: `0xffff` of system mask `0xffffff`.
- **M1 (per step).** Each `mpm analyze 1 dt` is one profiler run. From its tree, following the brief:
  - **t_rb** = `mpm.a.prepare` + `domainChanged` + every `soe.symbolic` outside `domainChanged`. PARDISO's symbolic runs in the first solve; UmfPack's runs inside `setSize`.
  - **t_it** = (formTangent + linearSolve − symbolic + formUnbalance + update + convTest) / k, where k = the step's `linearSolve` calls.
  - **R** = t_rb / t_it. **f** = t_rb / wall(`mpm.a.step`).
  - Per run, over the timed steps: t_it = Σ iteration time / Σ k, and f = Σ t_rb / Σ wall.
  - Across repeats: the median, with the spread (max − min) / median in parentheses. n = repeats that produced data.
- **Normalisation of cut steps** (the orchestrator's instruction, 2026-11-16). One harness step (`mpm analyze 1 dt`) that fails and is cut runs several rebuild cycles: one `Domain::domainChange()` per attempt, plus the accepted sub-steps.
  - The brief does not define how a cut step is normalised. It defines k per accepted step, and t_rb as "everything a mode-A step spends because the Domain changed".
  - **The per-rebuild-cycle normalisation was the orchestrator's choice, made after the run (2026-11-16).** The headline R is (Σ t_rb / Σ `domainChanges`) / t_it.
  - Two other readings are given beside it, labelled: the aggregate's **R per harness step** (t_rb summed over every attempt of `mpm analyze 1 dt`), and **R per accepted sub-step** (Σ t_rb / Σ accepted sub-steps, over t_it).
  - On a step without a cut, all three are identical. No row depends on the choice: R ≤ 1 already fails on cut-free cells (10⁴, all four cases).
  - f is the share of all step time spent rebuilding, failed attempts included. That is the cost a user pays, and it is reported as is.
- **Models** (README). M1 and M3 mode A use a laterally confined 3D soft block (1 m cube, ppc 2, uGIMP, n = 5 / 11 / 23 cells: N = 1,000 / 10,648 / 97,336), in two loadings:
  - QS: a `LoadControl` ramp to 0.75 g in 8 steps.
  - Transient: Newmark 0.5 / 0.25, Δt = T/16.

  Materials: ElasticIsotropic, and DruckerPrager (σY 400 Pa, ρ 0.21, non-associated). Analysis: `test NormDispIncr 1e-6 50`, Newton, RCM, cutoff 5e-2. Steps: 2 warm-up + 6 timed (at least 4 timed at 10⁵).

  M2 and M3 mode F use a 3D quarter-column collapse (a × a × 2a, frictional floor, N = 1,024 / 11,664 / 93,312) in three twins: elastic, DP φ 30°, LadrunoJ2 2 kPa. Steps: 5 warm-up + 20 timed. The SANISAND projection comes from a single-point stdBrick driver (the fork's SANISAND deck).

  M4 is the 232f cantilever QS r1–r3, built by the 232f harness.

## Pinned conditions

| | |
|---|---|
| Build | `ladrunoBuild` = `9db5803f6d20672e7f73e89615af33a42a449dac` in all 233 per-run files that exist; the orchestrator verified the same hash on `OpenSees.exe` and `dist\bin\opensees.pyd` (all five targets, 19:02) |
| Machine | i7-13700KF (8 P + 8 E), 64 GB (seen: Intel64 Family 6 Model 183, 24 logical, 68.6 GB, Windows 11 26100) |
| Python | 3.12.10 |
| Threads | implicit 8 (MKL and OMP), explicit 1 (+ the default row); `ladrunoThreads` = 1 |
| Affinity | `0xFFFF` (the eight P-cores with HT), set before `opensees` loads |
| Repeats | per the orchestrator's plan (`run_all.py --plan full`): M1 default 3; UmfPack 3 / 3 / 1 at 10³ / 10⁴ / 10⁵ (2700 s timeout); ghost 3 / 3 / 1; exclude 1 (10³, 10⁴); M2 3; M3 2; M4 1 |
| Warm-up | M1 2 steps (excluded), M2 5 steps (excluded), M4 none (every step is data) |
| Run | started 2026-11-15 19:17; 7.08 h; 235 runs, 230 ok, 5 failed (below); no other heavy job (orchestrator) |

The failures, as the orchestrator decided to treat them:

- The **10⁵ QS DP** block fails to converge at step 7 (zero-based step index; load factor ≈ 0.68). Newton stalls just above the tolerance at 10⁵ only: the attempts at the smaller cut increments run 50 iterations and end at 1.1e-6 to 6.1e-6 against 1e-6 (renormalize; ghost 1.05e-6 to 1.0e-5). The first two attempts of the failing step, at the full and half increment, end at 1.2e-3 to 7.6e-3. N = 10³ and 10⁴ pass the same load, in 10 and 21 iterations at their last step. This is not a limit load: a confined block has none (`Ladruno_scripts/mpm_m0b/README.md:107`).
  - PARDISO renormalize r0 and r1, and PARDISO ghost r0, fail at step 7 after 4 cuts. Their accepted timed steps 2–6 are the data, labelled **"partial: Newton stalls at step 7 (10⁵ only)"**.
  - **PARDISO renormalize r2 crashed** during the cut 3 → 4 retry, with `Matrix::operator*(Vector): incompatible sizes` (end of its log), and wrote no result. The orchestrator is filing it as its own issue. The default 10⁵ QS DP cell therefore has **2 partial repeats, not 3**.
  - **UmfPack renormalize QS DP 10⁵**: a timeout at 2700 s, recorded as a timeout. UmfPack decides nothing alone.

## M1 — the rebuild against one Newton iteration

Columns: k = Newton iterations per harness step (all attempts); t_rb per rebuild cycle; R per cycle (headline) and per harness step (the aggregate's); f; the symbolic share of t_rb; formTangent's share of t_it; retries per repeat. Both shares are ratios of sums over the timed steps. The aggregate's own `formTangent_share_of_t_it` field is a per-step mean and differs where k varies (10⁵ QS DP: 0.760 against 0.68 here); both are above 0.30.

### M1 — PARDISO, renormalize (default)

| N | case | n | k / step | t_rb / cycle [ms] | t_it [ms] | **R** (per cycle) | R per harness step | f | symbolic / t_rb | formTangent / t_it | retries | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1,000 | qs elastic | 3 | 5.5 (0 %) | 15.9 (5 %) | 14.2 (2 %) | **1.13 (5 %)** | 1.13 (5 %) | 0.17 (5 %) | 0.18 (3 %) | 0.85 (0 %) | 0/0/0 |  |
| 1,000 | qs dp | 3 | 6.5 (0 %) | 23.3 (2 %) | 23.2 (3 %) | **1.01 (5 %)** | 1.01 (5 %) | 0.134 (4 %) | 0.14 (2 %) | 0.84 (1 %) | 0/0/0 |  |
| 1,000 | transient elastic | 3 | 4.33 (0 %) | 15.6 (1 %) | 15.6 (2 %) | **0.998 (2 %)** | 0.998 (2 %) | 0.185 (2 %) | 0.17 (3 %) | 0.84 (0 %) | 0/0/0 |  |
| 1,000 | transient dp | 3 | 4.67 (0 %) | 32.2 (2 %) | 32.2 (1 %) | **0.997 (1 %)** | 0.997 (1 %) | 0.174 (1 %) | 0.11 (4 %) | 0.85 (1 %) | 0/0/0 |  |
| 10,000 | qs elastic | 3 | 6 (0 %) | 268.5 (1 %) | 211.9 (2 %) | **1.26 (1 %)** | 1.26 (1 %) | 0.173 (1 %) | 0.1 (1 %) | 0.87 (0 %) | 0/0/0 |  |
| 10,000 | qs dp | 3 | 9 (4 %) | 381.3 (0 %) | 371.1 (1 %) | **1.03 (1 %)** | 1.03 (1 %) | 0.102 (4 %) | 0.088 (2 %) | 0.86 (1 %) | 0/0/0 |  |
| 10,000 | transient elastic | 3 | 6 (0 %) | 260.6 (0 %) | 224.9 (2 %) | **1.16 (2 %)** | 1.16 (2 %) | 0.161 (2 %) | 0.1 (2 %) | 0.87 (1 %) | 0/0/0 |  |
| 10,000 | transient dp | 3 | 7.67 (2 %) | 486.4 (1 %) | 453.1 (1 %) | **1.07 (0 %)** | 1.07 (0 %) | 0.122 (2 %) | 0.088 (1 %) | 0.86 (0 %) | 0/0/0 |  |
| 100,000 | qs elastic | 3 | 6.33 (0 %) | 4,852.6 (0 %) | 4,034.7 (5 %) | **1.2 (5 %)** | 1.2 (5 %) | 0.159 (5 %) | 0.081 (1 %) | 0.77 (4 %) | 0/0/0 |  |
| 100,000 | qs dp | 2 | 26.6 (0 %) | 5,572.6 (2 %) | 5,920.0 (4 %) | **0.941 (1 %)** | 1.32 (1 %) | 0.0472 (1 %) | 0.079 (4 %) | 0.68 (2 %) | 1/1 | **partial**: r0, r1 Newton stalls at step 7, just above tolerance (steps 2–6 used); r2 crashed (`Matrix::operator*(Vector): incompatible sizes`), no result |
| 100,000 | transient elastic | 3 | 6.83 (0 %) | 4,002.2 (0 %) | 3,494.6 (1 %) | **1.15 (1 %)** | 1.15 (1 %) | 0.143 (0 %) | 0.083 (2 %) | 0.82 (0 %) | 0/0/0 |  |
| 100,000 | transient dp | 3 | 37.7 (0 %) | 6,208.8 (3 %) | 5,989.9 (5 %) | **1.02 (5 %)** | 2.37 (5 %) | 0.0589 (5 %) | 0.078 (3 %) | 0.71 (3 %) | 4/4/4 |  |

### M1 — UmfPack, renormalize

| N | case | n | k / step | t_rb / cycle [ms] | t_it [ms] | **R** (per cycle) | R per harness step | f | symbolic / t_rb | formTangent / t_it | retries | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1,000 | qs elastic | 3 | 5.5 (0 %) | 13.5 (3 %) | 17.7 (3 %) | **0.758 (5 %)** | 0.758 (5 %) | 0.121 (4 %) | 0.03 (4 %) | 0.7 (1 %) | 0/0/0 |  |
| 1,000 | qs dp | 3 | 6.17 (0 %) | 20.9 (3 %) | 26.9 (1 %) | **0.777 (4 %)** | 0.777 (4 %) | 0.112 (4 %) | 0.023 (4 %) | 0.73 (1 %) | 0/0/0 |  |
| 1,000 | transient elastic | 3 | 4.33 (0 %) | 13.2 (6 %) | 18.8 (4 %) | **0.705 (4 %)** | 0.705 (4 %) | 0.138 (3 %) | 0.028 (15 %) | 0.68 (2 %) | 0/0/0 |  |
| 1,000 | transient dp | 3 | 4.67 (0 %) | 28.4 (3 %) | 37.4 (2 %) | **0.759 (1 %)** | 0.759 (1 %) | 0.138 (1 %) | 0.019 (7 %) | 0.76 (0 %) | 0/0/0 |  |
| 10,000 | qs elastic | 3 | 6 (0 %) | 266.8 (0 %) | 285.2 (1 %) | **0.937 (0 %)** | 0.937 (0 %) | 0.135 (0 %) | 0.078 (2 %) | 0.7 (0 %) | 0/0/0 |  |
| 10,000 | qs dp | 3 | 9.17 (0 %) | 364.1 (0 %) | 587.5 (0 %) | **0.62 (0 %)** | 0.62 (0 %) | 0.0633 (0 %) | 0.054 (2 %) | 0.56 (0 %) | 0/0/0 |  |
| 10,000 | transient elastic | 3 | 6 (0 %) | 257.3 (1 %) | 292.6 (1 %) | **0.885 (1 %)** | 0.885 (1 %) | 0.128 (1 %) | 0.072 (1 %) | 0.71 (0 %) | 0/0/0 |  |
| 10,000 | transient dp | 3 | 7.67 (0 %) | 455.1 (1 %) | 738.0 (3 %) | **0.617 (2 %)** | 0.617 (2 %) | 0.0741 (1 %) | 0.037 (3 %) | 0.56 (1 %) | 0/0/0 |  |
| 100,000 | qs elastic | 1 | 6.33 | 6,340.5 | 8,244.3 | **0.769** | 0.769 | 0.108 | 0.31 | 0.39 | 0 |  |
| 100,000 | qs dp | 0 | — | — | — | — | — | — | — | — | — | r0: **timeout** at 2700 s |
| 100,000 | transient elastic | 1 | 6.83 | 4,208.5 | 5,889.1 | **0.715** | 0.715 | 0.0944 | 0.14 | 0.48 | 0 |  |
| 100,000 | transient dp | 1 | 37.7 | 8,213.0 | 10,391.6 | **0.79** | 1.84 | 0.0466 | 0.34 | 0.42 | 4 |  |

### M1 — PARDISO, ghost 1e-2

| N | case | n | k / step | t_rb / cycle [ms] | t_it [ms] | **R** (per cycle) | R per harness step | f | symbolic / t_rb | formTangent / t_it | retries | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1,000 | qs elastic | 3 | 5.5 (0 %) | 15.5 (2 %) | 13.5 (3 %) | **1.16 (2 %)** | 1.16 (2 %) | 0.174 (2 %) | 0.17 (3 %) | 0.86 (0 %) | 0/0/0 |  |
| 1,000 | qs dp | 3 | 6.67 (0 %) | 27.8 (3 %) | 26.0 (0 %) | **1.07 (4 %)** | 1.07 (4 %) | 0.138 (3 %) | 0.15 (3 %) | 0.85 (0 %) | 0/0/0 |  |
| 1,000 | transient elastic | 3 | 4.33 (0 %) | 15.6 (1 %) | 15.0 (2 %) | **1.05 (3 %)** | 1.05 (3 %) | 0.192 (2 %) | 0.17 (3 %) | 0.85 (0 %) | 0/0/0 |  |
| 1,000 | transient dp | 3 | 4.67 (0 %) | 38.5 (0 %) | 36.7 (3 %) | **1.05 (3 %)** | 1.05 (3 %) | 0.182 (2 %) | 0.13 (2 %) | 0.85 (0 %) | 0/0/0 |  |
| 10,000 | qs elastic | 3 | 6 (0 %) | 315.7 (1 %) | 244.8 (2 %) | **1.29 (1 %)** | 1.29 (1 %) | 0.177 (0 %) | 0.1 (3 %) | 0.87 (0 %) | 0/0/0 |  |
| 10,000 | qs dp | 3 | 8.83 (4 %) | 435.1 (1 %) | 407.1 (1 %) | **1.07 (1 %)** | 1.07 (1 %) | 0.108 (4 %) | 0.094 (2 %) | 0.86 (1 %) | 0/0/0 |  |
| 10,000 | transient elastic | 3 | 6 (0 %) | 412.7 (0 %) | 361.7 (0 %) | **1.14 (0 %)** | 1.14 (0 %) | 0.159 (0 %) | 0.086 (3 %) | 0.88 (0 %) | 0/0/0 |  |
| 10,000 | transient dp | 3 | 7.67 (0 %) | 546.8 (0 %) | 502.7 (1 %) | **1.09 (0 %)** | 1.09 (0 %) | 0.124 (0 %) | 0.085 (1 %) | 0.86 (1 %) | 0/0/0 |  |
| 100,000 | qs elastic | 1 | 6.33 | 5,110.8 | 4,227.2 | **1.21** | 1.21 | 0.16 | 0.083 | 0.77 | 0 |  |
| 100,000 | qs dp | 1 | 15 | 4,802.5 | 5,094.3 | **0.943** | 0.943 | 0.0591 | 0.085 | 0.72 | 0 | **partial**: Newton stalls at step 7, just above tolerance (steps 2–6 used) |
| 100,000 | transient elastic | 1 | 6.83 | 4,906.8 | 4,465.2 | **1.1** | 1.1 | 0.138 | 0.079 | 0.78 | 0 |  |
| 100,000 | transient dp | 1 | 37.8 | 6,374.1 | 6,055.9 | **1.05** | 2.46 | 0.0608 | 0.08 | 0.73 | 4 |  |

### M1 — PARDISO, exclude

| N | case | n | k / step | t_rb / cycle [ms] | t_it [ms] | **R** (per cycle) | R per harness step | f | symbolic / t_rb | formTangent / t_it | retries | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1,000 | qs elastic | 1 | 5.5 | 14.8 | 12.8 | **1.15** | 1.15 | 0.173 | 0.18 | 0.86 | 0 |  |
| 1,000 | qs dp | 1 | 6.5 | 22.6 | 22.1 | **1.02** | 1.02 | 0.136 | 0.14 | 0.84 | 0 |  |
| 1,000 | transient elastic | 1 | 4.33 | 15.1 | 14.1 | **1.07** | 1.07 | 0.195 | 0.17 | 0.84 | 0 |  |
| 1,000 | transient dp | 1 | 4.67 | 30.3 | 31.0 | **0.978** | 0.978 | 0.171 | 0.11 | 0.84 | 0 |  |
| 10,000 | qs elastic | 1 | 6 | 268.4 | 212.2 | **1.26** | 1.26 | 0.174 | 0.099 | 0.87 | 0 |  |
| 10,000 | qs dp | 1 | 8.5 | 379.0 | 362.0 | **1.05** | 1.05 | 0.11 | 0.091 | 0.86 | 0 |  |
| 10,000 | transient elastic | 1 | 6 | 259.4 | 224.8 | **1.15** | 1.15 | 0.16 | 0.1 | 0.87 | 0 |  |
| 10,000 | transient dp | 1 | 7.67 | 484.4 | 454.9 | **1.06** | 1.06 | 0.121 | 0.088 | 0.86 | 0 |  |

### Cut steps, normalised per rebuild cycle

Every cell with retries > 0, in any item. "test" means a failed attempt that ran into the convergence test's limit; the log shows 50 iterations. "refusal" means an attempt ended by a particle's trial refusal (an inverted particle, or a centre leaving the grid) after 1 to 50 iterations; the log does not print the count, hence the ranges. Iterations per accepted sub-step = (all iterations − failed-attempt iterations) / accepted sub-steps.

| item | N | case | run | timed steps | rebuild cycles (domainChange) | accepted sub-steps | failed attempts (test / refusal) | t_rb per harness step [ms] | t_rb per cycle [ms] | t_rb per accepted sub-step [ms] | t_it [ms] | R per cycle | R per harness step | k per harness step | k per accepted sub-step | iterations in failed attempts |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| pardiso_renormalize | 100,000 | qs dp | r0 | 5 | 7 | 6 | 1 / 0 | 7,891.5 | 5,636.8 | 6,576.3 | 6,026.4 | 0.935 | 1.31 | 26.6 | 13.8 | 50 |
| pardiso_renormalize | 100,000 | qs dp | r1 | 5 | 7 | 6 | 1 / 0 | 7,711.7 | 5,508.3 | 6,426.4 | 5,813.6 | 0.947 | 1.33 | 26.6 | 13.8 | 50 |
| pardiso_renormalize | 100,000 | transient dp | r0 | 6 | 14 | 10 | 2 / 2 | 14,595.4 | 6,255.2 | 8,757.2 | 5,981.0 | 1.05 | 2.44 | 37.7 | 2.6 – 12.4 | 100 + 2 refusal(s) of 1–50 each |
| pardiso_renormalize | 100,000 | transient dp | r1 | 6 | 14 | 10 | 2 / 2 | 14,217.3 | 6,093.1 | 8,530.4 | 5,989.9 | 1.02 | 2.37 | 37.8 | 2.7 – 12.5 | 100 + 2 refusal(s) of 1–50 each |
| pardiso_renormalize | 100,000 | transient dp | r2 | 6 | 14 | 10 | 2 / 2 | 14,487.3 | 6,208.8 | 8,692.4 | 6,251.0 | 0.993 | 2.32 | 37.7 | 2.6 – 12.4 | 100 + 2 refusal(s) of 1–50 each |
| umfpack_renormalize | 100,000 | transient dp | r0 | 6 | 14 | 10 | 2 / 2 | 19,163.6 | 8,213.0 | 11,498.1 | 10,391.6 | 0.79 | 1.84 | 37.7 | 2.6 – 12.4 | 100 + 2 refusal(s) of 1–50 each |
| pardiso_ghost | 100,000 | transient dp | r0 | 6 | 14 | 10 | 2 / 2 | 14,873.0 | 6,374.1 | 8,923.8 | 6,055.9 | 1.05 | 2.46 | 37.8 | 2.7 – 12.5 | 100 + 2 refusal(s) of 1–50 each |

So at 10⁵ the large aggregate R values (2.3–2.5 transient DP; 1.31–1.33 QS DP) come from the cuts:

- One harness step runs 2–3 rebuild cycles: a failed attempt, then two accepted halves. Normalised per rebuild cycle, **R is 0.99–1.05** (transient DP) and **0.94–0.95** (QS DP).
- The rebuild cycle on these late steps is **30–60 % dearer**: 5.5–6.4 s against 4.0–4.3 s. t_it rises too, from 3.4–3.7 s to 6.0 s per iteration. The cut is not the cause. The uncut QS elastic block shows the same rise at the same stage (10⁵, step 6: `dc.s.graph` 1.93 → 2.97 s, formTangent 2.73 → 3.76 s per call, private memory 2,302 → 3,606 MB). Support growth causes both (see "Steps with and without a connectivity change").
- **R per accepted sub-step** (the third reading): transient DP 10⁵ **1.42** (1.39–1.46), QS DP 10⁵ 1.10.
- Transient DP k: 37.7 per harness step, but 2.6–12.4 per accepted sub-step. The failed attempts take 100 iterations plus two refusals.

### QS DP at 10⁵: with and without step 6 (the last accepted step before the stall)

Step 6 is kept in the headline numbers, as instructed. Its iterations (37 in the accepted attempt, 103 in the harness step, 46 in ghost) pull f down. They come from Newton approaching the stall of step 7, not from a limit load (see Pinned conditions).

| item | run | timed steps used | R per cycle | f | k per harness step | R per cycle without step 6 | f without step 6 | k without step 6 | step 6: iterations (testIter), retries | failed step 7: cuts, iterations |
|---|---|---|---|---|---|---|---|---|---|---|
| pardiso_renormalize | r0 | 2, 3, 4, 5, 6 | 0.935 | 0.0469 | 26.6 | 1.18 | 0.136 | 7.5 | 103 (37), 1 | 5 cuts, 326 iterations (all attempts) |
| pardiso_renormalize | r1 | 2, 3, 4, 5, 6 | 0.947 | 0.0475 | 26.6 | 1.11 | 0.129 | 7.5 | 103 (37), 1 | 5 cuts, 362 iterations (all attempts) |
| pardiso_renormalize | r2 | — | — | — | — | — | — | — | — | error: child exited 3221225477 without a result (see m1_N100000_qs_dp_Pardiso_renormalize_r2.log) |
| umfpack_renormalize | r0 | — | — | — | — | — | — | — | — | timeout: child exited timeout without a result (see m1_N100000_qs_dp_UmfPack_renormalize_r0.log) |
| pardiso_ghost | r0 | 2, 3, 4, 5, 6 | 0.943 | 0.0591 | 15.0 | 1.22 | 0.144 | 7.2 | 46 (46), 0 | 5 cuts, 359 iterations (all attempts) |

### The rebuild's sub-phases (default, per rebuild cycle)

### t_rb per rebuild cycle by sub-phase, default, N = 1,000 (ms, and % of t_rb; median over repeats)

| sub-phase | qs elastic | qs dp | transient elastic | transient dp |
|---|---|---|---|---|
| `mpm.a.support` | 0.72 (5 %) | 0.93 (4 %) | 0.72 (5 %) | 1.14 (4 %) |
| `mpm.a.sp` | 0.03 (0 %) | 0.03 (0 %) | 0.03 (0 %) | 0.03 (0 %) |
| `mpm.a.churn` | 0.01 (0 %) | 0.01 (0 %) | 0.01 (0 %) | 0.01 (0 %) |
| `mpm.a.p2g` | 0.19 (1 %) | 0.23 (1 %) | 0.20 (1 %) | 0.29 (1 %) |
| `mpm.a.repoint` | 0.73 (5 %) | 1.07 (5 %) | 0.71 (5 %) | 1.20 (4 %) |
| `mpm.a.domainChange` | 0.00 (0 %) | 0.00 (0 %) | 0.00 (0 %) | 0.00 (0 %) |
| `mpm.a.prepare.other` | 0.03 (0 %) | 0.04 (0 %) | 0.03 (0 %) | 0.03 (0 %) |
| `dc.clearAll` | 0.08 (0 %) | 0.15 (1 %) | 0.07 (0 %) | 0.30 (1 %) |
| `dc.handle` | 0.34 (2 %) | 0.63 (3 %) | 0.32 (2 %) | 1.04 (3 %) |
| `dc.numberDOF` | 4.35 (27 %) | 6.74 (29 %) | 4.38 (28 %) | 9.37 (29 %) |
| `dc.s.graph` | 5.84 (37 %) | 9.08 (39 %) | 5.70 (37 %) | 13.4 (42 %) |
| `dc.setSize.soe` | 0.62 (4 %) | 1.15 (5 %) | 0.61 (4 %) | 1.62 (5 %) |
| `dc.integrator` | 0.00 (0 %) | 0.00 (0 %) | 0.01 (0 %) | 0.02 (0 %) |
| `dc.algorithm` | 0.00 (0 %) | 0.00 (0 %) | 0.00 (0 %) | 0.00 (0 %) |
| `dc.other` | 0.07 (0 %) | 0.08 (0 %) | 0.07 (0 %) | 0.09 (0 %) |
| `soe.symbolic` | 2.80 (18 %) | 3.19 (14 %) | 2.68 (17 %) | 3.41 (11 %) |
| **t_rb / cycle** | 15.9 | 23.3 | 15.6 | 32.2 |

### t_rb per rebuild cycle by sub-phase, default, N = 10,000 (ms, and % of t_rb; median over repeats)

| sub-phase | qs elastic | qs dp | transient elastic | transient dp |
|---|---|---|---|---|
| `mpm.a.support` | 10.8 (4 %) | 12.8 (3 %) | 11.0 (4 %) | 14.7 (3 %) |
| `mpm.a.sp` | 0.10 (0 %) | 0.11 (0 %) | 0.16 (0 %) | 0.17 (0 %) |
| `mpm.a.churn` | 0.04 (0 %) | 0.04 (0 %) | 0.06 (0 %) | 0.07 (0 %) |
| `mpm.a.p2g` | 2.90 (1 %) | 3.33 (1 %) | 2.96 (1 %) | 3.97 (1 %) |
| `mpm.a.repoint` | 10.3 (4 %) | 13.6 (4 %) | 11.0 (4 %) | 15.6 (3 %) |
| `mpm.a.domainChange` | 0.00 (0 %) | 0.00 (0 %) | 0.00 (0 %) | 0.00 (0 %) |
| `mpm.a.prepare.other` | 0.33 (0 %) | 0.33 (0 %) | 0.33 (0 %) | 0.36 (0 %) |
| `dc.clearAll` | 1.13 (0 %) | 1.61 (0 %) | 1.23 (0 %) | 3.12 (1 %) |
| `dc.handle` | 5.06 (2 %) | 9.87 (3 %) | 5.06 (2 %) | 15.5 (3 %) |
| `dc.numberDOF` | 75.0 (28 %) | 104.4 (27 %) | 72.1 (28 %) | 132.3 (27 %) |
| `dc.s.graph` | 118.5 (44 %) | 174.2 (46 %) | 114.6 (44 %) | 223.7 (46 %) |
| `dc.setSize.soe` | 16.4 (6 %) | 26.2 (7 %) | 14.4 (6 %) | 32.7 (7 %) |
| `dc.integrator` | 0.00 (0 %) | 0.00 (0 %) | 0.11 (0 %) | 0.11 (0 %) |
| `dc.algorithm` | 0.00 (0 %) | 0.00 (0 %) | 0.00 (0 %) | 0.00 (0 %) |
| `dc.other` | 0.93 (0 %) | 1.05 (0 %) | 0.90 (0 %) | 1.06 (0 %) |
| `soe.symbolic` | 27.0 (10 %) | 33.4 (9 %) | 27.0 (10 %) | 42.8 (9 %) |
| **t_rb / cycle** | 268.5 | 381.3 | 260.6 | 486.4 |

### t_rb per rebuild cycle by sub-phase, default, N = 100,000 (ms, and % of t_rb; median over repeats)

| sub-phase | qs elastic | qs dp | transient elastic | transient dp |
|---|---|---|---|---|
| `mpm.a.support` | 142.6 (3 %) | 153.6 (3 %) | 128.2 (3 %) | 161.7 (3 %) |
| `mpm.a.sp` | 1.06 (0 %) | 1.10 (0 %) | 1.28 (0 %) | 0.97 (0 %) |
| `mpm.a.churn` | 0.73 (0 %) | 0.71 (0 %) | 0.98 (0 %) | 0.64 (0 %) |
| `mpm.a.p2g` | 38.1 (1 %) | 40.0 (1 %) | 35.1 (1 %) | 44.2 (1 %) |
| `mpm.a.repoint` | 145.1 (3 %) | 161.9 (3 %) | 132.5 (3 %) | 176.0 (3 %) |
| `mpm.a.domainChange` | 0.00 (0 %) | 0.00 (0 %) | 0.00 (0 %) | 0.00 (0 %) |
| `mpm.a.prepare.other` | 2.62 (0 %) | 2.91 (0 %) | 2.47 (0 %) | 3.23 (0 %) |
| `dc.clearAll` | 34.7 (1 %) | 41.2 (1 %) | 28.6 (1 %) | 61.4 (1 %) |
| `dc.handle` | 146.0 (3 %) | 169.1 (3 %) | 106.5 (3 %) | 220.6 (4 %) |
| `dc.numberDOF` | 1,333.6 (27 %) | 1,511.7 (27 %) | 1,130.2 (28 %) | 1,635.8 (26 %) |
| `dc.s.graph` | 2,266.7 (47 %) | 2,647.4 (48 %) | 1,869.8 (47 %) | 2,947.1 (47 %) |
| `dc.setSize.soe` | 340.1 (7 %) | 390.6 (7 %) | 222.5 (6 %) | 453.0 (7 %) |
| `dc.integrator` | 0.00 (0 %) | 0.00 (0 %) | 0.93 (0 %) | 0.91 (0 %) |
| `dc.algorithm` | 0.00 (0 %) | 0.00 (0 %) | 0.00 (0 %) | 0.00 (0 %) |
| `dc.other` | 11.7 (0 %) | 12.3 (0 %) | 11.0 (0 %) | 12.0 (0 %) |
| `soe.symbolic` | 392.2 (8 %) | 439.9 (8 %) | 331.9 (8 %) | 476.8 (8 %) |
| **t_rb / cycle** | 4,852.6 | 5,572.6 | 4,002.2 | 6,208.8 |

**Largest sub-phase, everywhere: `dc.s.graph`**, the analysis model's DOF-graph construction: 37–48 % of t_rb. With `dc.numberDOF` (RCM over that graph, 26–29 %), the two make up **65–75 % of every rebuild**. `soe.symbolic` (PARDISO's METIS reorder) takes 8–18 %. The host's own preparation (`mpm.a.*`) is under 12 % together: support 3–5 %, repoint 3–5 %, P2G about 1 %, SP rebuild and grid-node churn under 0.1 %. The raise itself, the integrator's and the algorithm's `domainChanged` are negligible; Newmark's node re-read costs 0.9 ms at 10⁵.

### Steps with and without a connectivity change

| item | N | case | n changed / unchanged (per run) | t_rb / cycle, changed [ms] | t_rb / cycle, unchanged [ms] | ratio | symbolic share of t_rb, changed steps |
|---|---|---|---|---|---|---|---|
| pardiso_renormalize | 1,000 | qs elastic | 3/3, 3/3, 3/3 | 16.0 | 15.7 | 1.02 | 0.18 (2 %) |
| pardiso_renormalize | 1,000 | qs dp | 3/3, 3/3, 3/3 | 31.7 | 14.7 | 2.16 | 0.12 (4 %) |
| pardiso_renormalize | 1,000 | transient elastic | 2/4, 2/4, 2/4 | 16.0 | 15.3 | 1.05 | 0.17 (6 %) |
| pardiso_renormalize | 1,000 | transient dp | 4/2, 4/2, 4/2 | 40.4 | 15.3 | 2.63 | 0.096 (3 %) |
| pardiso_renormalize | 10,000 | qs elastic | 6/0, 5/1, 5/1 | 268.5 | 261.3 | 1.03 | 0.1 (1 %) |
| pardiso_renormalize | 10,000 | qs dp | 5/1, 6/0, 5/1 | 402.6 | 262.1 | 1.54 | 0.087 (3 %) |
| pardiso_renormalize | 10,000 | transient elastic | 4/2, 4/2, 4/2 | 259.8 | 262.1 | 0.991 | 0.11 (2 %) |
| pardiso_renormalize | 10,000 | transient dp | 4/2, 4/2, 4/2 | 598.5 | 262.3 | 2.28 | 0.085 (1 %) |
| pardiso_renormalize | 100,000 | qs elastic | 6/0, 6/0, 6/0 | 4,852.6 | — | — | 0.081 (1 %) |
| pardiso_renormalize | 100,000 | qs dp | 5/0, 5/0 | 5,572.6 | — | — | 0.079 (4 %) |
| pardiso_renormalize | 100,000 | transient elastic | 5/1, 5/1, 5/1 | 3,998.8 | 4,000.8 | 1 | 0.082 (3 %) |
| pardiso_renormalize | 100,000 | transient dp | 5/1, 5/1, 5/1 | 6,378.2 | 4,022.4 | 1.59 | 0.078 (3 %) |
| umfpack_renormalize | 1,000 | qs elastic | 4/2, 4/2, 4/2 | 13.6 | 13.4 | 1.02 | 0.03 (6 %) |
| umfpack_renormalize | 1,000 | qs dp | 5/1, 5/1, 5/1 | 22.4 | 12.9 | 1.73 | 0.023 (5 %) |
| umfpack_renormalize | 1,000 | transient elastic | 2/4, 2/4, 2/4 | 13.2 | 13.3 | 0.989 | 0.031 (12 %) |
| umfpack_renormalize | 1,000 | transient dp | 4/2, 4/2, 4/2 | 35.9 | 13.3 | 2.71 | 0.017 (8 %) |
| umfpack_renormalize | 10,000 | qs elastic | 5/1, 5/1, 5/1 | 268.9 | 255.4 | 1.05 | 0.078 (0 %) |
| umfpack_renormalize | 10,000 | qs dp | 5/1, 5/1, 5/1 | 385.3 | 254.5 | 1.51 | 0.05 (3 %) |
| umfpack_renormalize | 10,000 | transient elastic | 4/2, 4/2, 4/2 | 257.4 | 259.8 | 0.991 | 0.068 (4 %) |
| umfpack_renormalize | 10,000 | transient dp | 4/2, 4/2, 4/2 | 553.4 | 262.0 | 2.11 | 0.026 (8 %) |
| umfpack_renormalize | 100,000 | qs elastic | 6/0 | 6,340.5 | — | — | 0.31 |
| umfpack_renormalize | 100,000 | transient elastic | 5/1 | 4,194.2 | 4,279.9 | 0.98 | 0.14 |
| umfpack_renormalize | 100,000 | transient dp | 5/1 | 8,515.6 | 4,279.1 | 1.99 | 0.35 |
| pardiso_ghost | 1,000 | qs elastic | 2/4, 2/4, 2/4 | 15.7 | 15.5 | 1.01 | 0.16 (3 %) |
| pardiso_ghost | 1,000 | qs dp | 5/1, 5/1, 5/1 | 30.3 | 15.0 | 2.02 | 0.14 (3 %) |
| pardiso_ghost | 1,000 | transient elastic | 2/4, 2/4, 2/4 | 15.4 | 15.6 | 0.986 | 0.16 (3 %) |
| pardiso_ghost | 1,000 | transient dp | 4/2, 4/2, 4/2 | 50.3 | 15.1 | 3.34 | 0.12 (3 %) |
| pardiso_ghost | 10,000 | qs elastic | 5/1, 6/0, 5/1 | 316.7 | 291.0 | 1.09 | 0.1 (2 %) |
| pardiso_ghost | 10,000 | qs dp | 6/0, 5/1, 5/1 | 461.2 | 288.9 | 1.6 | 0.094 (2 %) |
| pardiso_ghost | 10,000 | transient elastic | 4/2, 4/2, 4/2 | 473.0 | 291.9 | 1.62 | 0.083 (2 %) |
| pardiso_ghost | 10,000 | transient dp | 4/2, 4/2, 4/2 | 674.5 | 291.2 | 2.32 | 0.083 (2 %) |
| pardiso_ghost | 100,000 | qs elastic | 6/0 | 5,110.8 | — | — | 0.083 |
| pardiso_ghost | 100,000 | qs dp | 5/0 | 4,802.5 | — | — | 0.085 |
| pardiso_ghost | 100,000 | transient elastic | 5/1 | 5,087.5 | 4,002.9 | 1.27 | 0.078 |
| pardiso_ghost | 100,000 | transient dp | 5/1 | 6,556.5 | 4,003.6 | 1.64 | 0.08 |
| pardiso_exclude | 1,000 | qs elastic | 2/4 | 14.9 | 14.8 | 1.01 | 0.17 |
| pardiso_exclude | 1,000 | qs dp | 5/1 | 24.2 | 14.9 | 1.62 | 0.13 |
| pardiso_exclude | 1,000 | transient elastic | 2/4 | 15.6 | 14.8 | 1.06 | 0.16 |
| pardiso_exclude | 1,000 | transient dp | 4/2 | 38.2 | 14.5 | 2.64 | 0.099 |
| pardiso_exclude | 10,000 | qs elastic | 5/1 | 271.1 | 254.5 | 1.07 | 0.1 |
| pardiso_exclude | 10,000 | qs dp | 6/0 | 379.0 | — | — | 0.091 |
| pardiso_exclude | 10,000 | transient elastic | 4/2 | 261.1 | 256.0 | 1.02 | 0.1 |
| pardiso_exclude | 10,000 | transient dp | 4/2 | 597.7 | 257.9 | 2.32 | 0.086 |

- **The full rebuild runs whether or not the connectivity changed.** On the default treatment, for elastic, t_rb per cycle on unchanged steps equals t_rb on changed steps within noise (ratio 0.98–1.05) at 10³, 10⁴ and 10⁵. Not in every item: ghost transient elastic gives 1.62 (10⁴) and 1.27 (10⁵), ghost 10⁴ QS elastic 1.09, exclude 1.06–1.07. In those items, too, the larger values sit on the late steps. The rebuild does not look at whether anything moved. This is an observation; it triggers no row, because f ≤ 0.25.
- **Why DP's t_rb on changed steps is 1.5–2.7× its unchanged t_rb** (QS DP 10³: 31.7 ms against 14.7), with no retries. The phase breakdown shows the cause, but the cause is the step index, not the connectivity change.
  - DP's changed steps are its late, plastic steps. In QS DP 10³ r0, t_rb is 14.6–16.4 ms at steps 2–5 and 39.5–40.1 ms at steps 6–7.
  - The jump is in the phases whose cost grows with the size of the DOF graph: `dc.s.graph` 5.8 → 16.7 ms (×2.9), `dc.numberDOF` 4.2 → 11.8 (×2.8), `dc.setSize.soe` 0.6 → 2.4 (×4), `dc.handle` 0.3 → 1.4 (×4).
  - The phases that are linear in the particle–node pairs grow less: `mpm.a.support` 0.8 → 1.4 (×1.75), `mpm.a.repoint` 0.9 → 2.0.
  - That pattern fits particles whose uGIMP support grows from 8 or 12 nodes to 18 or 27 as the plastic block deforms. The graph's edge count grows with the square of the element's DOFs (a ×1.75 support gives about ×3 in edges). An element with more than 21 support nodes (numDOF > 64) also gets its own tangent in `FE_Element` instead of the pool (`dc.handle`).
  - **Elastic runs jump too, later.** At 10⁵ QS elastic step 6, all three repeats show it: `dc.s.graph` 1.93 → 2.97 s, `dc.numberDOF` 1.16 → 1.72 s, formTangent 2.73 → 3.76 s per call, private 2,302 → 3,606 MB. Ghost transient elastic jumps as well (its changed/unchanged ratio of 1.62 and 1.27). So the late cost jumps come from support growth as the block deforms, in elastic and DP runs alike, with or without a cut. DP reaches that state earlier.
  - The support histogram per step was not recorded in M1, so the mechanism is an inference from the phase pattern and the memory jump. M3's profiler counts support it (see M3).
- The **symbolic share of t_rb on connectivity-changed steps**, default: 0.10–0.18 at 10³, 0.085–0.10 at 10⁴, **0.078–0.082 at 10⁵**.

### Scaling fits (default)

| case | N points (particles) | R per cycle | slope R | R per harness step | slope | f | slope f | t_rb / cycle [ms] | slope t_rb | t_it [ms] | slope t_it |
|---|---|---|---|---|---|---|---|---|---|---|---|
| qs elastic | 1,000, 10,000, 100,000 | 1.13, 1.26, 1.2 | 0.0133 | 1.13, 1.26, 1.2 | 0.0133 | 0.17, 0.173, 0.159 | -0.0141 | 15.9, 268.5, 4,852.6 | 1.24 | 14.2, 211.9, 4,034.7 | 1.23 |
| qs dp | 1,000, 10,000, 100,000 | 1.01, 1.03, 0.941 | -0.0144 | 1.01, 1.03, 1.32 | 0.0587 | 0.134, 0.102, 0.0472 | -0.226 | 23.3, 381.3, 5,572.6 | 1.19 | 23.2, 371.1, 5,920.0 | 1.2 |
| transient elastic | 1,000, 10,000, 100,000 | 0.998, 1.16, 1.15 | 0.0299 | 0.998, 1.16, 1.15 | 0.0299 | 0.185, 0.161, 0.143 | -0.0561 | 15.6, 260.6, 4,002.2 | 1.2 | 15.6, 224.9, 3,494.6 | 1.18 |
| transient dp | 1,000, 10,000, 100,000 | 0.997, 1.07, 1.02 | 0.00433 | 0.997, 1.07, 2.37 | 0.188 | 0.174, 0.122, 0.0589 | -0.235 | 32.2, 486.4, 6,208.8 | 1.14 | 32.2, 453.1, 5,989.9 | 1.13 |

| N | worst R per cycle (case) | worst R per harness step (case) | worst f (case) |
|---|---|---|---|
| 1,000 | 1.13 (qs elastic) | 1.13 (qs elastic) | 0.185 (transient elastic) |
| 10,000 | 1.26 (qs elastic) | 1.26 (qs elastic) | 0.173 (qs elastic) |
| 100,000 | 1.2 (qs elastic) | 2.37 (transient dp) | 0.159 (qs elastic) |

Slope of the worst-case envelope (the worst R at each N): 0.0133 per cycle (points 1.13, 1.26, 1.2) and 0.161 per harness step (points 1.13, 1.26, 2.37). **The worst slope of any single case per cycle is 0.030** (transient elastic). Per accepted sub-step the slopes are 0.013 / 0.019 / 0.030 / **0.077** (QS elastic / QS DP / transient elastic / transient DP; points for transient DP 0.997, 1.073, 1.424). All are ≤ 0.1.

## M1-t — empty-node treatments

Observable: the top-layer settlement, after the last step that both runs accepted. Cost: the median per-step `mpm.a.step` wall. Ghost has 3 repeats at 10³ and 10⁴ and 1 at 10⁵; exclude has 1 repeat.

| N | case | settlement, renormalize [m] | settlement, ghost [m] | difference | step wall, renormalize [s] | step wall, ghost [s] | ghost cheaper by | note |
|---|---|---|---|---|---|---|---|---|
| 1,000 | qs elastic | -0.054053 | -0.053917 | -0.251 % | 0.094 | 0.089 | +4.9 % | settlement after step 7; wall = median per-step `mpm.a.step` |
| 1,000 | qs dp | -0.054737 | -0.054591 | -0.266 % | 0.104 | 0.103 | +0.8 % | settlement after step 7; wall = median per-step `mpm.a.step` |
| 1,000 | transient elastic | -0.141211 | -0.140932 | -0.197 % | 0.081 | 0.076 | +5.7 % | settlement after step 7; wall = median per-step `mpm.a.step` |
| 1,000 | transient dp | -0.143974 | -0.143704 | -0.188 % | 0.222 | 0.275 | -23.4 % | settlement after step 7; wall = median per-step `mpm.a.step` |
| 10,000 | qs elastic | -0.053997 | -0.053929 | -0.125 % | 1.533 | 1.792 | -16.9 % | settlement after step 7; wall = median per-step `mpm.a.step` |
| 10,000 | qs dp | -0.054654 | -0.054581 | -0.133 % | 1.748 | 2.005 | -14.7 % | settlement after step 7; wall = median per-step `mpm.a.step` |
| 10,000 | transient elastic | -0.140855 | -0.141332 | +0.338 % | 1.541 | 1.756 | -13.9 % | settlement after step 7; wall = median per-step `mpm.a.step` |
| 10,000 | transient dp | -0.144073 | -0.144424 | +0.244 % | 4.456 | 5.047 | -13.3 % | settlement after step 7; wall = median per-step `mpm.a.step` |
| 100,000 | qs elastic | -0.053880 | -0.054091 | +0.391 % | 25.632 | 26.657 | -4.0 % | settlement after step 7; wall = median per-step `mpm.a.step` |
| 100,000 | qs dp | -0.047721 | -0.047460 | -0.547 % | 27.522 | 25.535 | +7.2 % | settlement after step 6; wall = median per-step `mpm.a.step` |
| 100,000 | transient elastic | -0.140025 | -0.140022 | -0.002 % | 26.878 | 27.167 | -1.1 % | settlement after step 7; wall = median per-step `mpm.a.step` |
| 100,000 | transient dp | -0.142709 | -0.142761 | +0.036 % | 251.568 | 253.616 | -0.8 % | settlement after step 7; wall = median per-step `mpm.a.step` |

| N | case | settlement, exclude [m] | vs renormalize | step wall, exclude [s] |
|---|---|---|---|---|
| 1,000 | qs elastic | -0.053917 | -0.251 % | 0.085 |
| 1,000 | qs dp | -0.054591 | -0.266 % | 0.103 |
| 1,000 | transient elastic | -0.140932 | -0.197 % | 0.072 |
| 1,000 | transient dp | -0.143704 | -0.188 % | 0.218 |
| 10,000 | qs elastic | -0.053929 | -0.125 % | 1.563 |
| 10,000 | qs dp | -0.054581 | -0.133 % | 1.725 |
| 10,000 | transient elastic | -0.141489 | +0.450 % | 1.530 |
| 10,000 | transient dp | -0.144495 | +0.293 % | 4.510 |

- Exclude gives the same answer as ghost at 10³: the same particles, the same support. At 10⁴ transient it differs by 0.3–0.45 %.
- Ghost's step costs more at 10⁴ (13–17 %). It keeps every support node, so the graph is larger: its t_rb per cycle on changed steps is 1.82× renormalize's at 10⁴ transient elastic (473.0 against 259.8 ms).

## M2 — explicit throughput, material share, the SANISAND projection

**Throughput.** One thread unless stated. "Window" = N × 20 / Σ profiler `step` wall (coarse run). "Material share" comes from the deep run.

| N (particles) | material | threads | n | throughput [particle-steps/s] (window) | from the per-step median | ns / particle-step | material share (deep) | non-material share | material ns / particle-step | yield fraction (timed mean) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1,000 | elastic | 1 | 3+3 | 1,273,838 (4 %) | 1,273,315 (3 %) | 785 (4 %) | 0.0627 (1 %) | 0.937 (0 %) | 53.6 (3 %) | — |
| 1,000 | dp | 1 | 3+3 | 560,376 (0 %) | 572,179 (1 %) | 1,785 (0 %) | 0.553 (1 %) | 0.447 (2 %) | 1,036 (2 %) | 0.389 (0 %) |
| 1,000 | j2 | 1 | 3+3 | 1,108,213 (1 %) | 1,115,165 (1 %) | 902 (1 %) | 0.133 (3 %) | 0.867 (0 %) | 129 (3 %) | 0.312 (0 %) |
| 1,000 | elastic | default | 3+3 | 1,263,768 (3 %) | 1,261,550 (2 %) | 791 (3 %) | 0.0619 (6 %) | 0.938 (0 %) | 53.7 (6 %) | — |
| 10,000 | elastic | 1 | 3+3 | 1,368,474 (0 %) | 1,361,122 (1 %) | 731 (0 %) | 0.0841 (2 %) | 0.916 (0 %) | 68.1 (2 %) | — |
| 10,000 | dp | 1 | 3+3 | 662,387 (1 %) | 669,345 (1 %) | 1,510 (1 %) | 0.508 (0 %) | 0.492 (0 %) | 811 (0 %) | 0.0711 (0 %) |
| 10,000 | j2 | 1 | 3+3 | 1,038,203 (1 %) | 1,034,299 (1 %) | 963 (1 %) | 0.234 (1 %) | 0.766 (0 %) | 254 (3 %) | 0.111 (0 %) |
| 10,000 | elastic | default | 3+3 | 1,368,471 (3 %) | 1,344,546 (3 %) | 731 (3 %) | 0.0846 (4 %) | 0.915 (0 %) | 68.3 (5 %) | — |
| 100,000 | elastic | 1 | 3+3 | 1,384,270 (0 %) | 1,387,291 (0 %) | 722 (0 %) | 0.0897 (3 %) | 0.91 (0 %) | 71.4 (3 %) | — |
| 100,000 | dp | 1 | 3+3 | 698,955 (0 %) | 704,678 (1 %) | 1,431 (0 %) | 0.494 (1 %) | 0.506 (1 %) | 755 (0 %) | 0.0171 (0 %) |
| 100,000 | j2 | 1 | 3+3 | 1,026,191 (2 %) | 1,026,922 (2 %) | 974 (2 %) | 0.242 (1 %) | 0.758 (0 %) | 266 (1 %) | 0.0417 (0 %) |
| 100,000 | elastic | default | 3+3 | 1,370,607 (0 %) | 1,377,515 (1 %) | 730 (0 %) | 0.0906 (2 %) | 0.909 (0 %) | 72.9 (3 %) | — |

- **The default-threads elastic row equals the 1-thread row** (1.37 against 1.38 × 10⁶ at 10⁵). Mode F's particle loop is serial in v1, so threads change nothing.
- **What the DP and J2 columns measure at large N.** The yield fraction over the timed window (the DP column's `yield_fraction_timed_mean`) is **0.389 at 10³, 0.071 at 10⁴ and 0.017 at 10⁵**. In 25 steps of Δt = 0.25 Δt_c, a column 36 cells high hardly starts to collapse. So **the 10⁵ DP and J2 throughputs are mostly elastic-branch throughputs**: 96.6 % of the DP particles never yield in the window (end yield fraction 0.0343; the window mean is 0.017), and J2 is the same (0.042). **The N representative of a collapsing column is 10³** (39 % of the particles yielding).
  - The DP material cost per particle-step reflects this: 1,036 ns at 10³, 811 at 10⁴, 755 at 10⁵. The driver's DP cost is 661 ns per elastic call and 1,792 ns per plastic call.
- **Non-material share, elastic, 10⁵, 1 thread: 0.910** (spread 0 %). It is 0.937 at 10³ and 0.916 at 10⁴.

**The strain increment the driver matches.** The DP column's per-step particle strain-increment norm, one step after the timed window:

| N | DP: median per-step particle strain-increment norm, all particles | over yielded particles | yielded particles |
|---|---|---|---|
| 1,000 | 0.000859 (0 %) | 0.000944 (0 %) | 803 (0 %) |
| 10,000 | 1.25e-05 (0 %) | 0.000505 (0 %) | 1,726 (0 %) |
| 100,000 | 8.82e-12 (0 %) | 0.000256 (0 %) | 3,396 (0 %) |

**The driver (single-point stdBrick), per call (trial + commit).** Re-run by the orchestrator at the yielded increments, same build, 3 repeats each.

| material | increment | × 0.5 | × 1 | × 2 | elastic per call |
|---|---|---|---|---|---|
| LadrunoSANISAND, plastic | 2.561e-4 (yielded median, 10⁵) | 397,352 ns (2 %) | **611,392 ns** (3 %) | 850,578 ns (1 %) | 1,051 ns |
| DruckerPrager, plastic | 2.561e-4 | 1,793 ns | **1,792 ns** | 1,784 ns | 661 ns |
| LadrunoJ2, plastic | 2.561e-4 | 202 ns | 205 ns | 203 ns | 127 ns |
| ElasticIsotropic | — | — | — | — | 84 ns |
| LadrunoSANISAND, plastic | 9.436e-4 (yielded median, 10³) | 822,292 ns (1 %) | **1,159,174 ns** (4 %) | 1,649,080 ns (4 %) | 1,027 ns |
| DruckerPrager, plastic | 9.436e-4 | 1,781 ns | **1,790 ns** | 1,793 ns | 657 ns |

**SANISAND / DP plastic cost ratio: 341 (band 222–477) at 2.56e-4, and 647 (band 462–920) at 9.44e-4.** The yielded median at 10⁴ is 5.05e-4, between the two. SANISAND's cost grows with the increment, about as increment^0.5–0.6 over the band; DP and J2 are flat. A SANISAND figure is therefore meaningful only together with its increment.

**The defective first measurement** (kept for the record; labelled defective).

- In the full run, `matched_increment()` matched the driver to the DP column's median over **all** particles at the largest N. At 10⁵ that median is 8.8e-12, because most of the column is at rest in the timed window (yield fraction 0.017 on average, 0.034 at its end). The driver legs were measured at 8.8e-12 (band 4.4e-12 / 1.8e-11), which gave:
  - SANISAND plastic 7,270 / 7,140 / 6,996 ns, so SANISAND/DP = 4.03 / 3.96 / 3.83;
  - DP 1,802 ns, flat;
  - the SANISAND projection of `m2/m2_results.json` (`projection`): 7,946 ns per particle-step, about 126,000 particle-steps/s.
- **These SANISAND numbers are not the plastic cost** of a loading particle. DP's validation passed anyway, because DP's cost does not depend on the increment.
- `matched_increment()` now uses the median over the yielded particles and records the all-particle median beside it; a quirk entry records the trap.
- The pre-yield leg of the same full run (increment 2.45e-5) gave SANISAND 97,000–100,000 ns per call, consistent with the reruns' trend.

**The projection, recomposed with the yielded reruns.** c(X, φ) = (1 − φ) c_el + φ c_pl.

- **literal**: t = non-material(elastic column, deep) + c(X, φ).
- **difference**: t = t(elastic column, coarse) + c(X, φ) − c_el(elastic).

The DP validation uses φ(DP) = the DP column's yield fraction at that N.

| column N | phi(DP) | driver increment (source) | DP projected, literal (miss) | DP projected, difference (miss) | SANISAND, phi = 1, literal [ns] (×0.5 – ×2) | particle-steps/s | SANISAND, phi = phi(DP), literal [ns] | SANISAND / elastic step time (phi = 1) |
|---|---|---|---|---|---|---|---|---|
| 1,000 | 0.3887 | 9.436e-04 (yielded, N = 1000) | 1904 (9.0 %) | 1799 (3.0 %) | 1,159,981 (823,099 – 1,649,887) | 862 | 451,973 | 1,477 |
| 1,000 | 0.3887 | 2.561e-04 (yielded, N = 100000) | 1907 (9.1 %) | 1802 (3.1 %) | 612,198 (398,159 – 851,384) | 1,633 | 239,080 | 780 |
| 10,000 | 0.0711 | 9.436e-04 (yielded, N = 1000) | 1485 (0.6 %) | 1389 (7.0 %) | 1,159,922 (823,040 – 1,649,828) | 862 | 84,088 | 1,579 |
| 10,000 | 0.0711 | 2.561e-04 (yielded, N = 100000) | 1489 (0.3 %) | 1393 (6.8 %) | 612,139 (398,100 – 851,325) | 1,634 | 45,178 | 833 |
| 100,000 | 0.0171 | 9.436e-04 (yielded, N = 1000) | 1396 (1.6 %) | 1314 (7.4 %) | 1,159,894 (823,012 – 1,649,800) | 862 | 21,592 | 1,609 |
| 100,000 | 0.0171 | 2.561e-04 (yielded, N = 100000) | 1400 (1.3 %) | 1318 (7.1 %) | 612,111 (398,072 – 851,297) | 1,634 | 12,229 | 849 |

- **DP projection miss: ≤ 9.1 % (literal) and ≤ 7.4 % (difference)**, at every N and with either increment. The composition holds on DP.
- **SANISAND (φ = 1, all particles loading plastically): 612,000 ns per particle-step at 2.56e-4 (band 398,000–851,000), about 1,630 particle-steps/s; and 1,160,000 ns at 9.44e-4 (band 823,000–1,650,000), about 860 particle-steps/s.** That is 780–1,600 times the elastic column's step.
- With the measured φ(DP), the cost falls to 12,000–452,000 ns, because most particles are elastic.

## M3 — memory per particle-element

The slope of private bytes against N between the two largest N (10⁴ → 10⁵), after the first step and after the second. Two repeats each (they agree within 2 %).

| mode | material | N | private bytes / particle after step 1 (**the rule's quantity**) | after step 2 | build (s1 − s0) | step 1 (s2 − s1) |
|---|---|---|---|---|---|---|
| A | elastic | 10,648 → 97,336 | **17,754** | 22,023 | 893 | 16,862 |
| A | DP | 10,648 → 97,336 | **20,313** | 24,607 | 3,431 | 16,883 |
| F | elastic | 11,664 → 93,312 | 2,907 | 2,907 | 939 | 1,972 |
| F | DP | 11,664 → 93,312 | 5,442 | 5,442 | 3,468 | 1,973 |

**What holds it.** The analytical count next to the measured slope, per particle, at the 10⁵ run's measured support sizes.

- Mode A mean support: 12.9 nodes per particle. 3.8 % of the particles have 27 nodes, so numDOF = 81 > 64.
- Mode F mean support: 11.5 nodes.

| item | A elastic | A DP | F elastic | F DP |
|---|---|---|---|---|
| material copy (ElasticIsotropicThreeDimensional / DruckerPrager3D) | 272 | 2,320 | 272 | 2,320 |
| particle state + `mpm reset` snapshot | 1,038 | 1,038 | 1,038 | 1,038 |
| `LadrunoMPMPoint` element | 3,938 | 3,938 | — | — |
| `FE_Element` | 2,479 | 2,476 | — | — |
| … of which its OWN tangent (numDOF > 64 only, `FE_Element.cpp:188-219`) | 2,040 | 2,038 | — | — |
| grid nodes (0.14 per particle × Node + DOF_Group + host maps) | 198 | 198 | — | — |
| support arrays / trial arrays / dense grid (2.47 nodes per particle) | — | — | 466 / 344 / 381 | 466 / 344 / 381 |
| **analytic total** | **7,925** | **9,971** | **2,501** | **4,549** |
| **measured slope** | **17,754** | **20,313** | **2,907** | **5,442** |
| measured − analytic | 9,829 (55 %) | 10,343 (51 %) | 406 (14 %) | 893 (16 %) |

Where they disagree:

1. **Mode A: about 10 kB per particle (51–55 %) is not in the count.** The count leaves out the SOE matrix and PARDISO's factor (MKL-internal).
   - The SOE has 0.39 equations per particle, so its CSR matrix is about 2 kB per particle at 3D GIMP's stencil. The factor's fill would be most of the rest. That is an inference: MKL's memory is not itemised.
   - The **growth from step 1 to step 2 (+4.3 kB per particle)** is in the same unattributed share: the second rebuild's factor, or allocator retention.
2. **The FE_Element's own tangent is the largest counted item in mode A: 2.0 kB per particle.** That is 52 kB × the 3.8 % of particles with 27 support nodes.
   - This reconciles with the profiler. An FE_Element owns an 81×81 tangent (52,488 B) once numDOF > 64 (`FE_Element.cpp:188-219`). The `-memory` window holds 115.75 MB of `Matrix` storage after step 1 (about 2,205 own tangents, 1.33 kB per particle by the slope). Step 2 adds 80.4 MB (about 1,532 more): 3,737, against the 3,729 elements with 27 nodes in the histogram. The histogram, and with it the 2.0 kB, is taken after step 2; the slope, after step 1.
   - The share of 27-node particles grows as the block deforms, and this item grows with it (below).
3. **The DP material copy:** the measured build slope difference DP − elastic is 2,538 B, against 2,048 B counted. DruckerPrager3D holds about 490 B more than its header count.
4. **Mode F:** the count is 14–16 % short. That is allocator slack and the id registry (`unordered_set`), not modelled.

**The same quantity later in the run.** The M1 runs of the same block record private bytes after every step. Their slope (10,648 → 97,336, default treatment, median of the repeats), in kB per particle:

| case | step 0 (= M3's "after the first step") | steps 1–5 | late steps |
|---|---|---|---|
| QS elastic | 17.9 | 22.3–22.8 | **37.5** (steps 6–7) |
| QS DP | 20.5 | 24.9–25.2 | **39.0** (step 6) |
| transient elastic | 18.0 | 21.3–22.6 | 20.8–21.6 (steps 6–7: no jump) |
| transient DP | 20.5 | 24.7–25.2 (steps 1–3) | **37.5–39.3** (steps 4–7) |

- **Cause.** The jump of about +15 kB per particle happens at the step where `dc.s.graph`, `dc.numberDOF` and formTangent jump (M1). It is the same +1.3–1.4 GB at 10⁵ under UmfPack, so it is not PARDISO's factor. It fits elements whose support grows to 27 nodes and which then own a 52,488 B tangent each (numDOF 81 > 64, `FE_Element.cpp:188-219`). Roughly 29 % more particles would cross into that state.
- **Bound.** If the whole jump were own tangents, they would hold (0.29 + 0.038) × 52.5 kB ≈ 17 kB of the 37.5 kB: **at most about 46 %**.
- **The pre-registered instant is the model's minimum.** "After the first step" is the state with the smallest supports: particles start on a GIMP support boundary and the block deforms in one direction.

## M4 — the 232f data owed (cantilever QS, UmfPack, renormalize)

| case | particles | cutoff | completed steps / 20 | stops at step | Newton iterations per step (testIter) | tip at 20 [m] |
|---|---|---|---|---|---|---|
| r1 | 640 | 5e-2 | 20 | — | 4 5 5 5 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 | **−2.374756** |
| r1 | 640 | 1e-6 | 6 | 7 (5 cuts, 150 iterations) | 4 5 6 6 6 6 | — |
| r1 | 640 | 1e-10 | 6 | 7 (5 cuts, 125 iterations) | 4 5 8 6 7 6 | — |
| r2 | 2,560 | 5e-2 | 20 | — | 5 5 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 7 8 8 | **−2.387597** |
| r2 | 2,560 | 1e-6 | 8 | 9 (5 cuts) | 5 6 29 6 6 6 6 6 | — |
| r2 | 2,560 | 1e-10 | 8 | 9 (5 cuts; 2 retries already at steps 2 and 4) | 5 5 15 6 6 6 6 6 | — |
| r3 | 10,240 | 5e-2 | 20 | — | 5 6 6 6 6 6 6 6 6 6 7 7 7 7 7 7 8 8 8 8 | **−2.396258** |
| r3 | 10,240 | 1e-6 | 4 | 5 (5 cuts) | 5 6 6 6 | — |
| r3 | 10,240 | 1e-10 | 3 | 4 (5 cuts) | 5 15 6 | — |

- **The tip at the small cutoffs, against 5e-2, at the last step both reached.** For every overlapping step, see "Applying the pre-registered table", M4: the largest difference is at step 2.

  | case | step | 5e-2 | 1e-6 | 1e-10 |
  |---|---|---|---|---|
  | r1 | 6 | −0.7189 | −0.7444 (+3.5 %) | −0.7444 (+3.5 %) |
  | r2 | 8 | −0.9799 | −0.9969 (+1.7 %) | −0.9970 (+1.7 %) |
  | r3 | 4 | −0.4933 | −0.5020 (+1.8 %) | — |
  | r3 | 3 | −0.3668 | — | −0.3768 (+2.7 %) |

  At no resolution does renormalize converge through the ramp at 1e-6 or 1e-10.
- **The 5e-2 tips reproduce 232f** (r1 −2.37476). The series −2.3748, −2.3876, −2.3963 passes the FE hypo value −2.3917 (`tests/mpm_implicit_cases.txt:130`; EICR −2.4046) between r2 and r3.
- **formTangent's share of t_it.**
  - Cantilever (2D, UmfPack): 0.43–0.49 at 5e-2.
  - **At 10⁵ (M1, 3D, the default PARDISO): 0.68 (QS DP) – 0.82 (transient elastic).** With UmfPack at 10⁵: 0.39–0.48.
  - At 10³ and 10⁴ (PARDISO): 0.84–0.87.

## Applying the pre-registered table — the worker's proposal

This is a proposal; the orchestrator applies the table. Each item quotes the rule's inputs, the deciding number and the row it falls in.

**M1 — whether mode C is ever needed.**

- **Read on:** the default (PARDISO, `renormalize` 5e-2), at every N, taking the worst of QS / transient and elastic / DP.
- **Worst run at each N:**

  | N | worst R per cycle | worst f |
  |---|---|---|
  | 10³ | **1.13**, QS elastic (spread 5 %) | **0.185**, transient elastic |
  | 10⁴ | **1.26**, QS elastic (spread 1 %) | **0.173**, QS elastic |
  | 10⁵ | **1.20**, QS elastic (spread 5 %) | **0.159**, QS elastic |

  The aggregate's per-harness-step R, labelled, is worst at 10⁵ for transient DP: 2.37, from cut steps.
- **Slope** of R against N, per rebuild cycle: **worst single case 0.030** (transient elastic). By case: 0.013 (QS elastic), −0.014 (QS DP), 0.030 (transient elastic), 0.004 (transient DP); the worst-case envelope gives 0.013. Per accepted sub-step the worst is 0.077 (transient DP). Per harness step it is 0.188 (transient DP, envelope 0.16), driven by the cut steps at 10⁵.
- **Symbolic share of t_rb on connectivity-changed steps:** 0.078–0.082 at 10⁵ (≤ 0.18 at any N). **Largest rebuild sub-phase:** `dc.s.graph` (47–48 % at 10⁵), then `dc.numberDOF` (26–28 %).
- **C-closed** needs f ≤ 0.25 ✓ (≤ 0.185 everywhere), R ≤ 1 ✗ (1.13 / 1.26 / 1.20) and slope ≤ 0.1 ✓ (≤ 0.030 per cycle, ≤ 0.077 per accepted sub-step; 0.188 per harness step ✗). **Not met:** R exceeds 1 at every N, by 13–26 % in the worst case. At 10⁴ it does so in all four cases (1.03–1.26, 12 cut-free runs), under every normalisation.
- **C-cheaper-first** needs f > 0.25 ✗. Not met (f ≤ 0.185). Its avoidable-without-C argument covers only steps whose connectivity did not change, where the full rebuild (graph and numbering first) runs anyway. At 10⁵ those are 0–3 of 18 timed steps per case (QS elastic 0 of 18, transient 3 of 18; QS DP 0 of 10 kept).
- **C-needed** needs f > 0.5 at 10⁵ ✗ (0.159) and a symbolic share ≥ 50 % ✗ (8 %). Not met.
- → **Proposed row: C-deferred** ("anything else"). f is small: rebuilding is 5–19 % of a step at every N. For elastic it holds at 0.17 / 0.17 / 0.16; the DP fall at 10⁵ comes from the iteration count (26.6, 37.7). R stays flat in N (slope ≤ 0.030 per cycle). But R is just above 1: one rebuild costs 1.0–1.3 Newton iterations. Two near misses:
  - C-closed fails only on R ≤ 1;
  - C-cheaper-first fails only on its f > 0.25 entry condition.

  Per the table, C stays deferred; the numbers are recorded and revisited at phase 11.
- **Conditions of this reading.**
  - f ≤ 0.25 was measured with **k ≥ 4 Newton iterations on every timed step**. f ≈ R / (R + k): the excluded first step at 10⁵ QS elastic has k = 4 and f 0.253 (all three repeats); at 10⁴ (k = 3–4, with one-time element creation) it is 0.35–0.39. A workload converging in 3 iterations or fewer per step, with the same t_rb and t_it, would sit above 0.25.
  - R and f are measured for **today's formTangent**, which is 68–87 % of t_it (ratio of sums; 73–87 % as the per-step mean) and which M4's rule sends to an issue (item 15). A what-if from the measured cells, not a measurement: halving formTangent puts f at about **0.27 (10⁴ QS elastic) and 0.28 (10³ transient elastic)**, above 0.25.
  - So "revisit" means at phase 11 **and after item 15**.
- *UmfPack* (decides nothing alone): R per cycle is 0.62–0.94 at 10³–10⁴ and 0.72–0.79 at 10⁵ (QS DP 10⁵ timed out); f ≤ 0.14; flat in N. It is not in a worse row than PARDISO: its iteration is dearer, so R is lower. → No guide change from the UmfPack clause.

**M1-t — whether the default treatment changes.**

- **Inputs:** ghost against renormalize, the top-layer settlement and the per-step cost, at every N.
- **Answer clause** (within 0.5 %): met at 10³ and 10⁴ (≤ 0.34 %) and for three of four cases at 10⁵ (≤ 0.39 %). **Not met for 10⁵ QS DP: −0.547 %**, compared at step 6 (the last accepted step before the stall), single repeat.
- **Cost clause** (≥ 25 % cheaper at every N): **not met at 10³ or 10⁴**. By the median per-step wall, ghost ranges from 23 % dearer to 7 % cheaper. By the mean, from 60 % dearer (10⁴ transient elastic, 2.60 against 1.62 s) to 51 % cheaper (10⁵ QS DP, partial: one ghost run, 81 against 165 s). The clause fails either way at 10³ and 10⁴.
- → **Proposed row: the default stands** (`renormalize` 5e-2). `exclude` is cost-only and anchors free surfaces: never a candidate.

**M2 — baselines (cannot reopen C).**

- **Non-material share, elastic, 10⁵, 1 thread: 0.910 > 0.70** → **a performance issue on F's transfer loops is filed** (P2G, grid update, the two projections; the material is 9 % of the step).
- **DP projection miss: ≤ 9.1 % (literal) / ≤ 7.4 % (difference) < 25 %**, with the yielded-increment drivers → **the SANISAND row is reported as reliable** by the rule. SANISAND at φ = 1: 612,000 ns per particle-step (2.56e-4) or 1,160,000 ns (9.44e-4).
- **SANISAND / DP cost ratio, for phase 5:** 341 (222–477) at 2.56e-4; 647 (462–920) at 9.44e-4.
- The rule's DP validation cannot see the increment dependence (DP is flat): see below.

**M3 — the memory side of C.**

- **Input:** mode A, private bytes per particle-element after the first step, 10⁴ → 10⁵.
- **Deciding numbers:** **17.8 kB (elastic) and 20.3 kB (DP)**; 22.0 and 24.6 kB after step 2.
- → **Proposed row, by the letter: between 8 and 32 kB → recorded, revisit at phase 11.**
- **For the owner:** the pre-registered instant is the model's minimum. The same slope, read later in the same runs, is **37.5–39.3 kB** in three of the four default cases (QS elastic, QS DP, transient DP); transient elastic stays at 21–23 kB.
  - Read late, the row would be **"> 32 kB"**. Its escape clause (one storage item ≥ 50 % removable without C) is then not met by the measured state: FE_Element's own tangents are at most about 46 % of 37.5 kB.
  - A block with 27-node supports everywhere was not measured.
  - The orchestrator applies "between" by the letter and puts the late reading to the owner.

**M4 — data.**

- **Bias clause:** renormalize does **not** converge at 1e-6 or 1e-10 at any of r1–r3; it stops at steps 7/7, 9/9 and 5/4. The clause's precondition ("converges … at all three resolutions") is not met → **the default 5e-2 stands.** **The default stands because the small-cutoff runs could not finish, not because the tips agree.** Where the runs overlap, the small-cutoff tips are larger in magnitude by up to **8.7 % (r1), 4.9 % (r2) and 3.2 % (r3)**, each at step 2. They are 7.3 % and 7.9 % at r1 steps 3–4, and 1.7–3.6 % at the last common step. Step 1 is identical. r1 is outside the 5 % tolerance on three steps.
- **Item 15:** formTangent / t_it at 10⁵ = **0.68–0.82 > 30 %** (default PARDISO; 0.39–0.48 with UmfPack) → **item 15 (the tangent B-build cleanup) is filed as an issue.**

**What the table cannot decide, or decides only by a catch-all.**

1. **M1's outcome** (f ≤ 0.19, R 1.13–1.26 flat in N) fits C-deferred ("anything else"). It is a near miss of C-closed on R ≤ 1, and of C-cheaper-first on f > 0.25. The latter also depends on the iteration count and on formTangent's cost (Conditions, above).
2. **R per rebuild cycle against R per harness step:**
   - The worst single-case slope is 0.030 per cycle, 0.077 per accepted sub-step and 0.188 per harness step. The headline follows the orchestrator's per-cycle choice, made after the run.
   - The per-harness-step slope would fail C-closed's slope clause too, at the 10⁵ transient DP point, where 4 of the 6 timed steps are cut. No row changes: R ≤ 1 fails on cut-free cells.
3. **SANISAND's cost depends on the strain increment per call:**
   - The driver gives 7 µs at 8.8e-12, 0.61 ms at 2.6e-4 and 1.16 ms at 9.4e-4. Which increment a SANISAND MPM particle would see cannot be measured while SANISAND is not admitted at particles. The DP composition check, the rule's only safeguard, is blind to it.
   - The proposed phase-5 ratio is therefore a range (222–920), not a number.
   - Also: at 10⁵ the DP and J2 columns barely yield (1.7 %, 4.2 %), so their 10⁵ throughputs are elastic-branch figures; 10³ is the collapsing-column N.
4. **The 10⁵ QS DP cell:** Newton stalls just above tolerance at step 7, at 10⁵ only (10³ and 10⁴ pass the same load). The cell holds 2 partial repeats, a crash (filed by the orchestrator) and a UmfPack timeout. Its numbers carry step 6's 37–46 iterations; without step 6, R per cycle is 1.11–1.18 and f is 0.13–0.14.
5. **The 10⁵ M1-t QS DP comparison** (−0.547 %) is a single step just before the stall. It does not change the M1-t outcome, which the cost clause decides.
6. **M3's instant.** Read by the letter (after the first step), M3 is "between". Read a few steps later in the same runs, it is "> 32 kB" in three of four cases, with the escape clause unmet (own tangents ≤ ~46 %). To the owner.
