---
title: WP-232i — B-spline bases and the manufactured-solution orders (results)
project: Ladruno
tags:
  - mpm
  - verification
---

# WP-232i — B-spline bases and the manufactured-solution orders

ADR 232 §8 (P3) and §9 phase 3. Frozen interface: `tests/mpm_bspline_interface.md` (B1–B12, the
orchestrator's freeze with F1–F4). Pre-registered gates: `.claude/skills/ladruno-mpm/references/wp-briefs.md`
("232i"), copied into `tests/mpm_bspline_cases.txt`. A band is never widened; a miss is reported as
a miss, with its diagnosis.

**Status: the warrant run is pending** (the orchestrator launches `tests/mpm_bspline_warrant.py` on
the 232i build). The tables below are filled from `tests/data/mpm/bspline/warrant.json` and the
mutant runs `g6_zero.json`, `g6_scale.json`.

## Method

- **Manufactured solution** (Wallstedt & Guilkey 2008), `tests/mpm_mms.py`: unit square / cube,
  rollers on every face, `u_i = A sin(2πX_i) sin(cπt)`, c = 2 (period 1), ElasticIsotropic E = 1,
  ν = 0.3, ρ₀ = 1. The body force uses the stress the kernel integrates: for a coaxial stretch the
  Hughes–Winget increments sum to ln F (to O(Δt²)), the material returns σ = λ tr(ln F) I + 2μ ln F
  and the engine treats it as Cauchy, so b = a − (1/ρ₀) ∇_X·(J σ F⁻ᵀ). The derivation is checked by
  centred differences through an independent route (`tests/test_mpm_mms_reference.py`).
- **Driver** (`tests/ladruno_mpm_bspline_harness.py`): the body force is set every step with
  `mpm particleBodyForce` — at t_n in mode F (the symplectic-Euler / leapfrog pairing of MUSL), at
  t_{n+1} or λ_{n+1} in mode A; the initial velocity with `mpm particleVelocity`. Error
  e = RMS_p |u_p − u(X_p, t_end)| / A at t_end = one period; p is the least-squares slope over the
  three finest rungs. Δt is fixed across the ladder at 0.1 × the finest rung's critical Δt; the
  validity check halves Δt at the finest rung and requires a change in e below 5 %.
- **Boundary treatment (K4):** images on `fix` faces (interface B5, B6, B8).

## Results

_Pending the warrant run._ Per row: e per rung, p, e(T/4), max_t e(t), crossings per rung, Δt,
wall time, the build hash.

## Core-level pre-check (g++, not the gate)

The OpenSees-free core with the interface's hypoelastic law (the ElasticIsotropic law), mode F,
1/h = 8, 16, 32 at Δt = 0.1 × the 1/h = 32 critical Δt (rates between consecutive rungs):

| basis | max \|F_ii − 1\| = 1e-3 | 0.1 |
|---|---|---|
| bspline2 | 2.04, 2.07 | 2.07, 2.05 |
| bspline3 | 1.96, 2.08 | 1.96, 2.08 |
| gimp | 2.19, 2.13 | 1.94, 1.92 |
| linear | 2.20, 2.14 | 1.42, 0.99 |

G8 (images against the full domain) on the core: rel 0 / 3e-16 (positions / stresses) at small
amplitude, 3e-16 / 4.5e-15 at large amplitude (max |F_ii − 1| = 0.108), after 120 steps.
