---
title: WP-232i — B-spline bases and the manufactured-solution orders (results)
project: Ladruno
tags:
  - mpm
  - verification
---

# WP-232i — B-spline bases and the manufactured-solution orders

ADR 232 §8 (P3) and §9 phase 3. Frozen interface: `tests/mpm_bspline_interface.md` (B1–B12, the
orchestrator's freeze with F1–F4, the round-2 and critic-round notes §5b). Pre-registered gates:
`.claude/skills/ladruno-mpm/references/wp-briefs.md` ("232i", with the rows the freeze and the
critic round added before any number of theirs), copied into `tests/mpm_bspline_cases.txt`. Every
row below is evaluated against that table; no band was widened, and nothing reported here is
re-gated.

## Gates applied

| Row | Gate (pre-registered) | Result | Verdict |
|---|---|---|---|
| G1 | basis properties, cardinal values, C¹ | g++ core 92/92; binary `mpm debug shape` ≥ 10⁴ positions per basis and ndm | **pass** |
| G1-K4 | folded S′ near each `fix` face | g++ core and binary, rollers and clamped faces; min S′ = 0 exactly | **pass** |
| G2 | bspline2, small amplitude: p ≥ 1.8, e non-increasing | p = 2.015 (1/h 8–128) | **pass** |
| G3 | bspline2, large amplitude: p ≥ 1.5, e non-increasing | p = 1.828 (8–128; 679,048 crossings at 1/h = 128) | **pass** |
| G4 | mode A quasi-static, bspline2: p ≥ 1.8, e non-increasing | p = 2.005 (8–128) | **pass** |
| G5 | G2–G4 for bspline3 | p = 2.018 / 2.017 / 2.475 | **pass** |
| G6 | ZERO and SCALE each turn G2 (or G3) red | ZERO: G2 and G3 red; SCALE: G3 red (G2 survives, as predicted in Round 1) | **pass** |
| G7 | `info basis`, `.ladruno`, mode-A restart, K2 refusal, inadmissible particle | `tests/test_ladruno_mpm_bspline_basis.py`, Python and classic Tcl | **pass** |
| G8 | images against the full domain, mode F, rel 1e-10, ≥ 100 steps, two amplitudes | 120 steps: core ≤ 4.5e-15; binary within 1e-10 | **pass** |
| G8-A | G8 in mode A (Newmark, default cutoff `renormalize` 5e-2) — added by the critic round | on `113ada8a0`: stress 7.1 % / 2.8 % off (small / large); on `4673949fe`: within 1e-10 | **pass** (after the A1 fix) |
| G8-F-cut | G8 in mode F at `massCutoff 0.05` — added by the critic round | on `113ada8a0`: stress 53 % / 8.5 % off; on `4673949fe`: within 1e-10 (core 4.9e-13) | **pass** (after the A1 fix) |

**Suite on `4673949fe`.** Every `tests/test_*mpm*.py` (Python and classic Tcl): 553 passed,
8 skipped (the 232f fine gates without `LADRUNO_MPM_GATES`), 5 xfailed. The 232f fine-resolution
gates (`LADRUNO_MPM_GATES=1`): 5 passed, 3 xfailed — the same outcomes as on `ladruno` (the
footing stays a strict xfail).

**Mutation gate (ADR-87 D2), family MPM, on `4673949fe`** (`Ladruno_scripts/mutation_gate.py`,
Python half, `LADRUNO_TCL_EXE` pointing nowhere in every run; baseline 87 passing):

| Mutant | Passing | Score | Floor | 232i tests killed |
|---|---|---|---|---|
| ZERO | 43 | **0.506** | 0.46 | all 7 (the six MMS rows on 1/h 8–32 and the e(T/4) oracle) |
| SCALE | 53 | **0.391** | 0.38 | G4, G5-G4 and the e(T/4) oracle; G2, G3, G5-G2, G5-G3 survive |

SCALE's margin rests on the e(T/4) oracle (critic A2): without it SCALE would score
33/87 = 0.379, under its floor. The one-period rows survive SCALE on the coarse rungs because at
one period e is a drift norm; SCALE turns the warrant run's G3 ladder red (G6 below).

## Method and pinned conditions

- **Builds.** Everything on `ladrunoBuild` `4673949fefe6c63757501535221014bf8bef678a`: the warrant
  run (exe and pyd), the mutants (`dist\mut_MPM_ZERO`, `dist\mut_MPM_SCALE`,
  `-DLADRUNO_MUTATE_FAMILY=MPM`), the suite. MSVC release, the nmb workstation (Windows 11),
  `opensees.pyd` through `tests/_testbed`. Results: `tests/data/mpm/bspline/{warrant,g6_zero,g6_scale}.json`,
  written by `tests/mpm_bspline_warrant.py`, checked by `tests/test_mpm_bspline_warrant_results.py`
  (which recomputes p and monotonicity from e and h and requires one build). The earlier run on
  `113ada8a0` (before the A1 fix) gave the same numbers for every MMS row: those bodies fill the
  box, so no face-layer node falls near the cutoff.
- **Manufactured solution** (Wallstedt & Guilkey 2008), `tests/mpm_mms.py`: unit square (cube in
  R1-3d), rollers on every face (K4 images on B-spline grids), `u_i = A sin(2πX_i) sin(cπt)`,
  c = 2 (period T = 1), ElasticIsotropic E = 1, ν = 0.3, ρ₀ = 1, A = max|F_ii − 1| / (2π). The body
  force uses the stress the kernel integrates: for a coaxial stretch the Hughes–Winget increments
  sum to ln F (to O(Δt²)) and the material's σ = λ tr(ln F) I + 2μ ln F is the Cauchy stress, so
  b = a − (1/ρ₀) ∇_X·(J σ F⁻ᵀ). The derivation is checked by centred differences through an
  independent route (`tests/test_mpm_mms_reference.py`, ≤ 1e-6 relative to the balance's terms).
- **Driver.** `mpm particleBodyForce` every step at t_n (mode F) or t_{n+1} / λ_{n+1} (mode A);
  `mpm particleVelocity` for v(X, 0) = ωA sin(2πX_i). Error e = RMS_p |u_p − u(X_p, t_end)| / A at
  t_end = T; p = least-squares slope of log e against log h over the three finest rungs. Also
  recorded: e(T/4) (the step nearest T/4), max_t e(t), cell crossings (particle-steps whose cell
  index changed), N, wall time, and the local (pairwise) orders.
- **Δt rule.** Explicit rows: Δt fixed across the ladder at 0.1 × the finest rung's critical Δt,
  h / (c_p + ωA√ndm). Validity: at the finest rung, |e(Δt/2) − e(Δt)| / e(Δt) < 5 %, else Δt is
  halved (at most 3 times; a gate row that never meets it is void). Every gated row met it at the
  rule Δt. Mode-A transient (R1): Newmark average acceleration, Δt = 0.005.
- **1/h = 128.** Added to a gated ladder when 8× (quasi-static 4×) the 64-rung wall time is under
  10 min — on every gated row; the ladder then runs at the 128 rung's Δt.

## Results per row

### G2 — gate

MMS, mode F, 2D plane strain, small amplitude (max |F_ii - 1| = 1e-3): bspline2 p >= 1.8 and e non-increasing along the whole ladder  
basis `bspline2`, mode explicit, ndm 2, ppc 2, max |F_ii − 1| = 0.001, α_F = 1.0; Δt = 6.7253e-04; validity change at the finest rung 0.13%; 1/h = 128 added (64-rung probe 7.2 s).

| 1/h | N | e | e(T/4) | max_t e | crossings | wall (s) |
|---|---|---|---|---|---|---|
| 8 | 256 | 8.818e-01 | 6.786e-03 | 8.818e-01 | 0 | 0.3 |
| 16 | 1024 | 2.136e-01 | 4.444e-03 | 2.187e-01 | 0 | 1.1 |
| 32 | 4096 | 5.088e-02 | 1.515e-03 | 5.383e-02 | 0 | 4.2 |
| 64 | 16384 | 1.252e-02 | 4.062e-04 | 1.340e-02 | 0 | 16.9 |
| 128 | 65536 | 3.114e-03 | 1.045e-04 | 3.341e-03 | 0 | 83.9 |

p (three finest) = **2.015**; pairwise 2.046, 2.070, 2.022, 2.008; p of e(T/4) = 1.929; p of max_t e = 2.005. Verdict: **pass** (p=2.015 (>= 1.8), non-increasing=True).

### G3 — gate

MMS, mode F, 2D plane strain, large amplitude (max |F_ii - 1| = 0.1; crossings per rung reported): bspline2 p >= 1.5 and e non-increasing  
basis `bspline2`, mode explicit, ndm 2, ppc 2, max |F_ii − 1| = 0.1, α_F = 1.0; Δt = 6.0020e-04; validity change at the finest rung 1.78%; 1/h = 128 added (64-rung probe 9.5 s).

| 1/h | N | e | e(T/4) | max_t e | crossings | wall (s) |
|---|---|---|---|---|---|---|
| 8 | 256 | 8.817e-01 | 9.282e-03 | 8.817e-01 | 0 | 0.4 |
| 16 | 1024 | 2.098e-01 | 6.030e-03 | 2.149e-01 | 744 | 1.4 |
| 32 | 4096 | 5.067e-02 | 4.112e-03 | 5.290e-02 | 11088 | 5.6 |
| 64 | 16384 | 1.297e-02 | 1.567e-03 | 1.355e-02 | 85240 | 22.5 |
| 128 | 65536 | 4.018e-03 | 7.330e-04 | 4.018e-03 | 679048 | 100.9 |

p (three finest) = **1.828**; pairwise 2.071, 2.050, 1.967, 1.690; p of e(T/4) = 1.244; p of max_t e = 1.859. Verdict: **pass** (p=1.828 (>= 1.5), non-increasing=True).

### G4 — gate

MMS, mode A, quasi-static, 2D, small amplitude (u = lam A sin(2 pi X_i), lam ramped in >= 4 load steps, body force updated per step): bspline2 p >= 1.8 and e non-increasing  
basis `bspline2`, mode quasiStatic, ndm 2, ppc 2, max |F_ii − 1| = 0.001, α_F = —; quasi-static, 4 load steps; 1/h = 128 added (64-rung probe 3.6 s).

| 1/h | N | e | e(T/4) | max_t e | crossings | wall (s) |
|---|---|---|---|---|---|---|
| 8 | 256 | 1.392e-02 | — | 1.392e-02 | 0 | 0.1 |
| 16 | 1024 | 3.278e-03 | — | 3.278e-03 | 0 | 0.2 |
| 32 | 4096 | 8.071e-04 | — | 8.071e-04 | 0 | 0.8 |
| 64 | 16384 | 2.010e-04 | — | 2.010e-04 | 0 | 3.6 |
| 128 | 65536 | 5.012e-05 | — | 5.012e-05 | 0 | 16.3 |

p (three finest) = **2.005**; pairwise 2.087, 2.022, 2.006, 2.003; p of e(T/4) = —; p of max_t e = 2.005. Verdict: **pass** (p=2.005 (>= 1.8), non-increasing=True).

### G5-G2 — gate

cubic (K1, delivered): G2 at the same thresholds for bspline3  
basis `bspline3`, mode explicit, ndm 2, ppc 2, max |F_ii − 1| = 0.001, α_F = 1.0; Δt = 6.7253e-04; validity change at the finest rung 0.10%; 1/h = 128 added (64-rung probe 8.9 s).

| 1/h | N | e | e(T/4) | max_t e | crossings | wall (s) |
|---|---|---|---|---|---|---|
| 8 | 256 | 1.067e+00 | 2.569e-02 | 1.067e+00 | 0 | 0.4 |
| 16 | 1024 | 2.747e-01 | 2.484e-03 | 2.785e-01 | 0 | 1.3 |
| 32 | 4096 | 6.505e-02 | 1.285e-03 | 6.840e-02 | 0 | 5.2 |
| 64 | 16384 | 1.596e-02 | 3.659e-04 | 1.701e-02 | 0 | 21.4 |
| 128 | 65536 | 3.966e-03 | 9.552e-05 | 4.242e-03 | 0 | 90.0 |

p (three finest) = **2.018**; pairwise 1.958, 2.078, 2.027, 2.009; p of e(T/4) = 1.875; p of max_t e = 2.006. Verdict: **pass** (p=2.018 (>= 1.8), non-increasing=True).

### G5-G3 — gate

cubic (K1, delivered): G3 at the same thresholds for bspline3  
basis `bspline3`, mode explicit, ndm 2, ppc 2, max |F_ii − 1| = 0.1, α_F = 1.0; Δt = 6.0020e-04; validity change at the finest rung 0.06%; 1/h = 128 added (64-rung probe 10.5 s).

| 1/h | N | e | e(T/4) | max_t e | crossings | wall (s) |
|---|---|---|---|---|---|---|
| 8 | 256 | 1.075e+00 | 2.721e-02 | 1.075e+00 | 0 | 0.4 |
| 16 | 1024 | 2.765e-01 | 3.532e-03 | 2.805e-01 | 744 | 1.5 |
| 32 | 4096 | 6.529e-02 | 1.657e-03 | 6.854e-02 | 11088 | 6.0 |
| 64 | 16384 | 1.602e-02 | 4.755e-04 | 1.704e-02 | 85208 | 24.1 |
| 128 | 65536 | 3.988e-03 | 1.387e-04 | 4.254e-03 | 679064 | 108.3 |

p (three finest) = **2.017**; pairwise 1.959, 2.082, 2.027, 2.006; p of e(T/4) = 1.789; p of max_t e = 2.005. Verdict: **pass** (p=2.017 (>= 1.5), non-increasing=True).

### G5-G4 — gate

cubic (K1, delivered): G4 at the same thresholds for bspline3  
basis `bspline3`, mode quasiStatic, ndm 2, ppc 2, max |F_ii − 1| = 0.001, α_F = —; quasi-static, 4 load steps; 1/h = 128 added (64-rung probe 9.2 s).

| 1/h | N | e | e(T/4) | max_t e | crossings | wall (s) |
|---|---|---|---|---|---|---|
| 8 | 256 | 3.653e-04 | — | 3.653e-04 | 0 | 0.1 |
| 16 | 1024 | 1.682e-05 | — | 1.682e-05 | 0 | 0.5 |
| 32 | 4096 | 1.094e-06 | — | 1.094e-06 | 0 | 1.9 |
| 64 | 16384 | 1.513e-07 | — | 1.513e-07 | 0 | 9.2 |
| 128 | 65536 | 3.536e-08 | — | 3.536e-08 | 0 | 39.7 |

p (three finest) = **2.475**; pairwise 4.441, 3.943, 2.854, 2.097; p of e(T/4) = —; p of max_t e = 2.475. Verdict: **pass** (p=2.475 (>= 1.8), non-increasing=True).

### R1-gimp-G2 — report

R1: G2 for gimp (literature: near 2 at small amplitude; an order below 1.5 is filed as an issue against shipped code)  
basis `gimp`, mode explicit, ndm 2, ppc 2, max |F_ii − 1| = 0.001, α_F = 1.0; Δt = 1.3451e-03; validity change at the finest rung 0.11%.

| 1/h | N | e | e(T/4) | max_t e | crossings | wall (s) |
|---|---|---|---|---|---|---|
| 8 | 256 | 5.390e-01 | 5.675e-02 | 5.390e-01 | 0 | 0.2 |
| 16 | 1024 | 1.177e-01 | 1.298e-02 | 1.196e-01 | 0 | 0.5 |
| 32 | 4096 | 2.683e-02 | 2.348e-03 | 2.777e-02 | 0 | 1.9 |
| 64 | 16384 | 6.473e-03 | 4.960e-04 | 6.741e-03 | 0 | 7.8 |

p (three finest) = **2.092**; pairwise 2.195, 2.133, 2.051; p of e(T/4) = 2.355; p of max_t e = 2.075.

### R1-gimp-G3 — report

R1: G3 for gimp  
basis `gimp`, mode explicit, ndm 2, ppc 2, max |F_ii − 1| = 0.1, α_F = 1.0; Δt = 1.2004e-03; validity change at the finest rung 0.11%.

| 1/h | N | e | e(T/4) | max_t e | crossings | wall (s) |
|---|---|---|---|---|---|---|
| 8 | 256 | 6.613e-01 | 4.734e-02 | 6.613e-01 | 0 | 0.2 |
| 16 | 1024 | 1.719e-01 | 7.981e-03 | 1.750e-01 | 744 | 0.6 |
| 32 | 4096 | 4.537e-02 | 2.887e-03 | 4.752e-02 | 11088 | 2.4 |
| 64 | 16384 | 1.289e-02 | 1.024e-03 | 1.348e-02 | 85016 | 9.6 |

p (three finest) = **1.868**; pairwise 1.944, 1.922, 1.815; p of e(T/4) = 1.481; p of max_t e = 1.849.

### R1-linear-G2 — report

R1: G2 for linear (literature: erratic once particles cross cells)  
basis `linear`, mode explicit, ndm 2, ppc 2, max |F_ii − 1| = 0.001, α_F = 1.0; Δt = 1.3451e-03; validity change at the finest rung 0.12%.

| 1/h | N | e | e(T/4) | max_t e | crossings | wall (s) |
|---|---|---|---|---|---|---|
| 8 | 256 | 5.378e-01 | 5.687e-02 | 5.378e-01 | 0 | 0.1 |
| 16 | 1024 | 1.171e-01 | 1.302e-02 | 1.189e-01 | 0 | 0.5 |
| 32 | 4096 | 2.653e-02 | 2.334e-03 | 2.744e-02 | 0 | 1.8 |
| 64 | 16384 | 6.374e-03 | 6.381e-04 | 6.598e-03 | 0 | 7.1 |

p (three finest) = **2.099**; pairwise 2.200, 2.142, 2.057; p of e(T/4) = 2.176; p of max_t e = 2.086.

### R1-linear-G3 — report

R1: G3 for linear  
basis `linear`, mode explicit, ndm 2, ppc 2, max |F_ii − 1| = 0.1, α_F = 1.0; Δt = 1.2004e-03; validity change at the finest rung 3.08%.

| 1/h | N | e | e(T/4) | max_t e | crossings | wall (s) |
|---|---|---|---|---|---|---|
| 8 | 256 | 5.533e-01 | 5.975e-02 | 5.533e-01 | 0 | 0.2 |
| 16 | 1024 | 2.054e-01 | 3.960e-02 | 2.054e-01 | 672 | 0.6 |
| 32 | 4096 | 9.768e-02 | 6.817e-02 | 2.171e-01 | 10008 | 2.2 |
| 64 | 16384 | 1.153e-01 | 4.856e-02 | 3.107e-01 | 81016 | 8.8 |

p (three finest) = **0.417**; pairwise 1.430, 1.072, -0.239; p of e(T/4) = -0.147; p of max_t e = -0.299.

### R1-pic-G2 — report

R1: the PIC leg of G2 (alpha_F = 0)  
basis `bspline2`, mode explicit, ndm 2, ppc 2, max |F_ii − 1| = 0.001, α_F = 0.0; Δt = 1.3451e-03; validity change at the finest rung 79.73% (NOT reached in 3 halvings).

| 1/h | N | e | e(T/4) | max_t e | crossings | wall (s) |
|---|---|---|---|---|---|---|
| 8 | 256 | 3.397e-02 | 9.375e-01 | 1.052e+00 | 0 | 0.2 |
| 16 | 1024 | 1.970e-02 | 7.917e-01 | 1.081e+00 | 0 | 0.5 |
| 32 | 4096 | 2.680e-01 | 4.818e-01 | 8.737e-01 | 0 | 2.0 |
| 64 | 16384 | 2.408e-01 | 1.758e-01 | 4.216e-01 | 0 | 8.4 |

p (three finest) = **-1.806**; pairwise 0.786, -3.765, 0.154; p of e(T/4) = 1.085; p of max_t e = 0.679.

Validity search at the finest rung: Δt 1.345e-03: e 2.408e-01 → e(Δt/2) 3.046e-01 (26.5%); Δt 6.725e-04: e 3.046e-01 → e(Δt/2) 2.694e-01 (11.6%); Δt 3.363e-04: e 2.694e-01 → e(Δt/2) 1.444e-01 (46.4%); Δt 1.681e-04: e 1.444e-01 → e(Δt/2) 2.926e-02 (79.7%).

### R1-pic-G3 — report

R1: the PIC leg of G3 (alpha_F = 0)  
basis `bspline2`, mode explicit, ndm 2, ppc 2, max |F_ii − 1| = 0.1, α_F = 0.0; Δt = 1.2004e-03; validity change at the finest rung 81.56% (NOT reached in 3 halvings).

| 1/h | N | e | e(T/4) | max_t e | crossings | wall (s) |
|---|---|---|---|---|---|---|
| 8 | 256 | 2.949e-02 | 9.437e-01 | 1.048e+00 | 0 | 0.2 |
| 16 | 1024 | 1.920e-02 | 8.089e-01 | 1.083e+00 | 0 | 0.6 |
| 32 | 4096 | 2.551e-01 | 5.123e-01 | 9.033e-01 | 976 | 2.4 |
| 64 | 16384 | 2.567e-01 | 1.940e-01 | 4.570e-01 | 56072 | 10.1 |

p (three finest) = **-1.870**; pairwise 0.619, -3.732, -0.009; p of e(T/4) = 1.030; p of max_t e = 0.623.

Validity search at the finest rung: Δt 1.200e-03: e 2.567e-01 → e(Δt/2) 3.074e-01 (19.7%); Δt 6.002e-04: e 3.074e-01 → e(Δt/2) 2.535e-01 (17.5%); Δt 3.001e-04: e 2.535e-01 → e(Δt/2) 1.235e-01 (51.3%); Δt 1.500e-04: e 1.235e-01 → e(Δt/2) 2.276e-02 (81.6%).

### R1-3x3-G2 — report

R1: the 3 x 3 particles-per-cell leg of G2  
basis `bspline2`, mode explicit, ndm 2, ppc 3, max |F_ii − 1| = 0.001, α_F = 1.0; Δt = 1.3451e-03; validity change at the finest rung 0.04%.

| 1/h | N | e | e(T/4) | max_t e | crossings | wall (s) |
|---|---|---|---|---|---|---|
| 8 | 576 | 7.981e-01 | 1.750e-02 | 7.981e-01 | 0 | 0.3 |
| 16 | 2304 | 1.939e-01 | 7.433e-04 | 1.982e-01 | 0 | 1.2 |
| 32 | 9216 | 4.649e-02 | 5.329e-04 | 4.893e-02 | 0 | 4.7 |
| 64 | 36864 | 1.147e-02 | 1.642e-04 | 1.219e-02 | 0 | 19.4 |

p (three finest) = **2.040**; pairwise 2.041, 2.060, 2.020; p of e(T/4) = 1.089; p of max_t e = 2.012.

### R1-3d — report

R1: the 3D small-amplitude MMS at 1/h = 4, 8, 16  
basis `bspline2`, mode explicit, ndm 3, ppc 2, max |F_ii − 1| = 0.001, α_F = 1.0; Δt = 5.3788e-03; validity change at the finest rung 0.14%.

| 1/h | N | e | e(T/4) | max_t e | crossings | wall (s) |
|---|---|---|---|---|---|---|
| 4 | 512 | 9.027e-01 | 2.939e-01 | 1.992e+00 | 0 | 0.1 |
| 8 | 4096 | 1.081e+00 | 1.559e-02 | 1.081e+00 | 0 | 1.1 |
| 16 | 32768 | 2.609e-01 | 3.015e-03 | 2.675e-01 | 0 | 8.4 |

p (three finest) = **0.895**; pairwise -0.260, 2.050; p of e(T/4) = 3.303; p of max_t e = 1.449.

### R1-modeA-transient — report

R1: mode A transient (Newmark average acceleration) MMS  
basis `bspline2`, mode implicit, ndm 2, ppc 2, max |F_ii − 1| = 0.001, α_F = 1.0; Δt = 5.0000e-03.

| 1/h | N | e | e(T/4) | max_t e | crossings | wall (s) |
|---|---|---|---|---|---|---|
| 8 | 256 | 8.821e-01 | 6.674e-03 | 8.821e-01 | 0 | 1.5 |
| 16 | 1024 | 2.141e-01 | 4.508e-03 | 2.191e-01 | 0 | 6.9 |
| 32 | 4096 | 5.132e-02 | 1.539e-03 | 5.427e-02 | 0 | 30.8 |

p (three finest) = **2.052**; pairwise 2.043, 2.061; p of e(T/4) = 1.058; p of max_t e = 2.011.

### Summary of orders

| row | status | p (e, gate) | p of e(T/4) | p of max_t e | verdict |
|---|---|---|---|---|---|
| G2 | gate | 2.015 | 1.929 | 2.005 | pass |
| G3 | gate | 1.828 | 1.244 | 1.859 | pass |
| G4 | gate | 2.005 | — | 2.005 | pass |
| G5-G2 | gate | 2.018 | 1.875 | 2.006 | pass |
| G5-G3 | gate | 2.017 | 1.789 | 2.005 | pass |
| G5-G4 | gate | 2.475 | — | 2.475 | pass |
| R1-gimp-G2 | report | 2.092 | 2.355 | 2.075 | reported |
| R1-gimp-G3 | report | 1.868 | 1.481 | 1.849 | reported |
| R1-linear-G2 | report | 2.099 | 2.176 | 2.086 | reported |
| R1-linear-G3 | report | 0.417 | -0.147 | -0.299 | reported |
| R1-pic-G2 | report | -1.806 | 1.085 | 0.679 | reported |
| R1-pic-G3 | report | -1.870 | 1.030 | 0.623 | reported |
| R1-3x3-G2 | report | 2.040 | 1.089 | 2.012 | reported |
| R1-3d | report | 0.895 | 3.303 | 1.449 | reported |
| R1-modeA-transient | report | 2.052 | 1.058 | 2.011 | reported |

### R1-column — the slice-1 self-weight column per basis and resolution

`column_selfweight_small_2d_gimp` with the basis replaced (width 1.0 so a B-spline grid has ≥ 2 cells; h = 0.5 / r, Δt / r, steps × r):

| basis | r = 1: selfweight_error / meter | r = 2 | r = 4 |
|---|---|---|---|
| gimp | 0.0128 / 0.0122 | 0.0065 / 0.0059 | 0.0034 / 0.0028 |
| linear | 0.0144 / 0.0125 | 0.0082 / 0.0063 | 0.0051 / 0.0033 |
| bspline2 | 0.0070 / 0.0018 | 0.0036 / 0.0007 | 0.0021 / 0.0003 |
| bspline3 | 0.0059 / 0.0013 | 0.0030 / 0.0005 | 0.0015 / 0.0002 |

`column_selfweight_2d_gimp` with the basis replaced (width 1.0 so a B-spline grid has ≥ 2 cells; h = 0.5 / r, Δt / r, steps × r):

| basis | r = 1: selfweight_error / meter | r = 2 | r = 4 |
|---|---|---|---|
| gimp | 0.0717 / 0.0205 | 0.0758 / 0.0224 | 0.0975 / 0.0321 |
| linear | 0.4721 / 0.1431 | 0.7666 / 0.2364 | 0.8311 / 0.3330 |
| bspline2 | 0.0407 / 0.0147 | 0.0789 / 0.0247 | 0.1349 / 0.0497 |
| bspline3 | 0.0136 / 0.0045 | 0.0158 / 0.0054 | 0.0129 / 0.0034 |

### G6 — mutation (ladders 1/h = 8–64, no validity search)

ZERO drops the P2G internal force; SCALE multiplies the G2P position increment by 1.5 (ADR §8).
The gate needs each mutant to turn G2 or G3 red.

| Mutant | G2 (p; verdict) | G3 (p; verdict) | G6 |
|---|---|---|---|
| ZERO | −0.033; red | −0.054; red | killed |
| SCALE | 1.980; survives (predicted in Round 1: at one period e is a drift norm) | 1.121; red | killed |

G6 `MPM=ZERO` (build `4673949fe`):
- G2: e = 6.978e+00, 8.059e+00, 8.356e+00, 8.433e+00; p = -0.033; non-increasing False → RED
- G3: e = 6.685e+00, 7.847e+00, 8.290e+00, 8.452e+00; p = -0.054; non-increasing False → RED

G6 `MPM=SCALE` (build `4673949fe`):
- G2: e = 1.323e+00, 3.226e-01, 7.836e-02, 2.074e-02; p = 1.980; non-increasing True → green (survives)
- G3: e = 1.312e+00, 3.004e-01, 8.809e-02, 6.348e-02; p = 1.121; non-increasing True → RED

## Findings (reported, not re-gated)

1. **The quadratic B-spline's local order falls at the fine rungs of the large-amplitude row.**
   G3 (bspline2) pairwise orders 2.07 / 2.05 / 1.97 / 1.69 (8→16 … 64→128): the least-squares
   1.828 passes its 1.5 threshold, but the finest pair is below 1.7 as cell crossings grow
   (85,240 → 679,048). The cubic holds 2.01 at the same pair (G5-G3: 1.96 / 2.08 / 2.03 / 2.01).
2. **The quarter-period error at large amplitude converges below second order for the quadratic
   B-spline.** Over the three finest rungs, p of e(T/4) is 1.244 for G3 (bspline2) against 1.789
   for G5-G3 (bspline3); at small amplitude both are near 2 (G2 1.929, G5-G2 1.875). The gate
   measures e at one period, where the exact displacement is zero and e is a drift norm; at the
   peak displacement, with particles crossing cells, the quadratic B-spline's error is limited by
   the particle-quadrature error under cell crossing (Steffen, Kirby & Berzins 2008), and the cubic,
   smoother across knots, does better. The gate stands as pre-registered; this is disclosed, not
   re-gated. Other rows: GIMP G3 1.481, mode-A transient 1.058 and the 3×3 leg 1.089 — their
   small-amplitude e(T/4) sits near 1e-4 to 1e-3 of A with irregular pairwise ratios, so those
   slopes are fits to a noisy small quantity.
3. **PIC dissipates; its error grows with the number of steps.** R1-pic-G2 / G3: p = −1.81 / −1.87
   and the Δt validity condition is never reached in 3 halvings (e at the finest rung 0.241 → 0.305
   → 0.269 → 0.144 → 0.029 as Δt halves, G2). Pure PIC removes a little energy every step, so at a
   fixed time a smaller Δt or a finer grid (more steps, more transfers) damps the solution more;
   the one-period drift norm then mixes amplitude loss with phase. A finding about PIC, not a
   defect: the gated rows are FLIP (α_F = 1).
4. **The 3D R1 ladder was too coarse for a three-point fit.** R1-3d: e = 0.903, 1.081, 0.261 at
   1/h = 4, 8, 16, so the least-squares p = 0.895 over all three; the 4 → 8 pair (−0.26) is
   pre-asymptotic (two cells per half wavelength), the 8 → 16 pair is 2.05 (a factor 4.14). At
   1/h = 8 the 3D error is 1.226 × the 2D one ≈ √(3/2), the third component. A ladder-design
   weakness of the pre-registered R1 row (4, 8, 16), not a 3D defect.
5. **Linear hat at large amplitude fails as the literature says** (p = 0.42; e grows from 1/h = 32
   to 64 once particles cross cells); at small amplitude, with no crossing, it converges at 2.10.
   **GIMP** gives 2.09 / 1.87 — above the 1.5 threshold that would have opened an issue against
   shipped code, so none is filed.
6. **The 20 %-compression self-weight column** (interface A4: "the acceptance test for phase 3's …
   B-spline basis"): **bspline3 meets the frozen 0.03 at every resolution** (0.0136 / 0.0158 /
   0.0129); **bspline2 does not** (0.041 / 0.079 / 0.135, growing with refinement, like GIMP's
   0.072 / 0.076 / 0.098 and the linear hat's 0.47–0.83). At small compression every basis
   converges, and the B-splines' cell-crossing meter is 7–17× below GIMP's at r = 4 (bspline2
   0.0003, bspline3 0.0002, GIMP 0.0028). Reported; the A4 strict xfail on GIMP is unchanged.
7. **Mode-A memory (F2).** Bytes per particle-element (M0b M3 method, slope 10³ → 10⁴, elastic):
   GIMP mode A 18,287 (base `6ee968060`) → 18,503 (232i), +1.2 %, inside the ±3 % run-to-run noise;
   GIMP mode F 2,918 → 2,933 (+0.5 %); bspline2 mode A 81,726, mode F 3,402. bspline2 mode A is
   past M0b M3's 32 kB "memory case for C" line because 78 % of its supports have 27 nodes and each
   such `FE_Element` owns an 81 × 81 tangent — the item #1098 (the `FE_Element`'s own tangent)
   already tracks. M3 was a decision about the GIMP default and is not applied here.
8. **The images and the cutoff (critic A1).** A face-layer node of a `fix` face is its own image,
   so its one-sided mass is 1/2^k of the body plus its reflections; the active-node cutoff dropped
   in the half model nodes the full model keeps (G8-A and G8-F-cut, 2.8–53 % in stress on
   `113ada8a0`). The cutoff now tests m × 2^k on B-spline grids (mode F step 4 and `activeNodes`,
   mode A's mass map); GIMP and linear are untouched. Quirk
   `WP-232i-mpm-images-face-node-flip-increment`.

## Follow-ups

- **#1103** — the code-review items left open at the ready-flip: mode A builds the per-component
  temporaries for every particle on a B-spline grid (only odd-folded particles need them); the
  restart pack v2 refuses a version-1 database (loud, documented in the guide); the particle body
  force is keyed by its length only (a re-ordered particle set of the same size is not detected);
  one helper for the face rules (`imageOfNode`, `zeroFaceOdd`, `imageCutoffMass` each re-derive the
  face test).
- **#1098** — the `FE_Element`'s own tangent (the bspline2 mode-A memory above).
- **nmorabowen/apeGmsh#3013** — the apeGmsh readers accept `GRID/BASIS = bspline2 | bspline3`
  (Contract 2 amended in 232a §2); Contract 1 v2 for B-spline particle sets is T's.

## Core-level pre-check (g++, before the build; not the gate)

The OpenSees-free core with the interface's hypoelastic law, mode F, 1/h = 8, 16, 32 at Δt = 0.1 ×
the 1/h = 32 critical Δt (rates between consecutive rungs): bspline2 2.04, 2.07 (small) / 2.07,
2.05 (large); bspline3 1.96, 2.08 / 1.96, 2.08; gimp 2.19, 2.13 / 1.94, 1.92; linear 2.20, 2.14 /
1.42, 0.99 — the binary ladders above reproduce them.
