---
title: "WP-105 / F12 — qualifying `IntScheme 2` (BackwardEuler_CPPM) for LadrunoSANISAND"
project: Ladruno
type: results
status: "PHASES A AND B COMPLETE — verdict PARTIAL: QUALIFIED as a prescribed-increment integrator (more accurate AND cheaper than IntScheme 1 at the campaign increment); REFUTED as the primary integrator of a load-controlled BVP (475x shallower than the baseline for the same wall clock)"
priority: high
owner: nmora
related:
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[_adr92_p0_oracle_results]]"
  - "[[93_ladruno_sanisand_zero_confinement_adr]]"
  - "[[LadrunoSANISAND_implex_guide]]"
tags: [adr, sanisand, integrator, cppm, measurement, wp-105]
updated: 2026-09-16
---

# F12 — is `IntScheme 2` usable on `LadrunoSANISAND`?

> [!success] **VERDICT — PARTIAL**
>
> **QUALIFIED where the strain increment is GIVEN.** On a replayed strain path
> (the ADR-92 G0 discipline, zero free DOF, so nothing but the constitutive
> integrator differs) scheme 2 converges to **the same answer** as scheme 1 —
> `1.3e-3` / `2.9e-3` maximum relative stress deviation over the whole path at
> `Δε_z = 1e-5`, terminal `η` within `4.2e-4` / `2.0e-4` at `p0 = 100` / `20 kPa` —
> and at the campaign's own increment it is **3.7–4.3× MORE accurate than the
> baseline and 4.2–7.6× cheaper per point**; at `4.6e-4` it is **7–30× more
> accurate and 10–13× cheaper**.
>
> **REFUTED as the primary integrator under a global Newton at the campaign
> increment.** On the free-standing drained-triaxial deck (5 free DOF) scheme 2
> stalls in **8 of 8** arms at `Δε_z ≥ 1e-4` where scheme 1 stalls in 1 of 8, and
> each failing step burns **12–134 s** against a 30 ms normal step (up to **4400×**)
> grinding `BackwardEuler_CPPM`'s recursive halving ladder. Loosening the global
> test from `1e-9` to `1e-7` does not rescue it, so this is not a solver-tolerance
> artefact. **Phase (b) settles it on the real deck:** on the CP1/ADR-95 bearing
> leg scheme 2 committed **11 steps to `s/B = 4e-5` in 1347 s** against the
> baseline's **51 steps to `s/B = 0.019` in 1267 s** — **475x shallower for the
> same wall clock**, with `ds` pinned at 25x the subdivision floor and **100 % of
> its steps on the relaxed rung 3**. The 1 % load-settlement bar could not be
> evaluated: the arms do not overlap in `s/B` at all (§7).
>
> **ADR-92 D3's stated rationale does not survive measurement.** D3 says scheme 2
> is "not an implicit return where the campaign's problem lives" because 58–74 %
> of its calls take the low-`p` branch and integrate explicitly. Measured on this
> build: **0 of 1820 steps** on every replayed triaxial path at either confinement,
> and **0 of 80** on the descent of the `p → p_min` path. The 58–74 % figure
> reproduces **only once the point is already pinned at `p_min` with a zero
> deviator** (85 of 160 steps there) — i.e. on steps where nothing is being
> integrated. **D3's conclusion (scheme 1 is the default) survives; its reason does
> not, and the correct reason is the ladder cost above, which D3 never measured.**

**Build.** Every number below is from
`634824e1fbcf802bf27c7fbd29c113e5a2d63cb6` (`ops.ladrunoBuild()`, asserted at the
top of every run) read out of
`.claude/worktrees/release-build-634824e1f/dist/bin/opensees.pyd`, interpreter
`python3.12`. **No source file was edited.** Work files:
`scratchpad/f12/` (`f12_matpoint.py`, `f12_bvp.py`, `run_*.sh`, `analyse.py`,
`data/`, `logs/`, `bvp/`).

---

## 0. Assumptions recorded (decisions taken without asking)

1. **The certificate.** ADR-92 publishes no scheme-vs-scheme tolerance. G0's
   `1e-8` is an *oracle-vs-binary* bar on the *same* algorithm and cannot be
   asked of a different integrator. I therefore stated a two-part bar and
   measured against it (§2):
   - **C1 CONSISTENCY** — at the finest increment both schemes agree to
     `≤ 1e-2` relative stress over the path and `≤ 1e-3` in terminal `η`
     (they integrate the same model to the same limit).
   - **C2 NO-DEGRADATION** — at the campaign increment scheme 2's own
     discretisation error, measured against its own fine reference, is **not
     larger** than scheme 1's against its own.
2. **The instrument.** The free-standing probe
   (`adr92_p0_oracle/probe_binary_triaxial.py`) has 5 free DOF and a global
   Newton, so a stall there confounds integrator with solver. I kept it (§4,
   it is the robustness measurement) but added a **strain-path replay** arm —
   the same cube with every face normal prescribed from the *baseline's own
   recorded path*, zero free DOF — which is where §2's accuracy and cost numbers
   come from. This is G0's own discipline with scheme 2 in the oracle's seat.
3. **`p0 = 20 kPa`** is not an ADR-92 case (P0 used 100, 5 and 1). It was run as
   specified; `p0 = 5 kPa` is the seed of the floor path in §5.
4. **The `p → p_min` case** is a prescribed volumetric **extension** ramped to
   `ε_v = +2.2e-4` by `t = 0.5` and then held, with deviatoric shear `6e-3·t` on
   throughout, from a committed `p0 = 5 kPa` — the ADR-92 G3 / ADR-93 ring
   construction, sized so the point *descends through* the low-`p` regime rather
   than arriving at the floor in one step (`2.2e-4` is the elastic strain that
   takes `p` from 5 kPa to 0 under `K ∝ √p`; P0's `3e-3` is 13× that and pins the
   point at the floor from step 1).
5. **The fallback detector.** `ManzariDafalias::debugFlag` is a compile-time
   `const bool = false` (`:57`), so the CPPM prints nothing when it falls back.
   No source edit was made to change that. Instead the run reads the shipped
   `substeps` response (`LadrunoSANISAND::setResponse`, `:4044-4058`), whose
   first component is `mSubstepsTakenInME` — non-zero **only** if `ModifiedEuler`
   ran. On scheme 2 that is an exact detector of the CPPM's explicit fallback.
6. Constants, `-Pmin 0.0101`, `-Presidual 0.0`, `-honorTolR 0`, `TolF = TolR =
   1e-10`, `TanType 2`, `LadrunoBrick -formulation bbar`, and the
   consolidate-then-push sequence are the probe's, unchanged. `-maxSubsteps` is
   0 (uncapped, vanilla) except in the §6 control.

---

## 1. What was run

| arm | deck | what it measures |
|---|---|---|
| **replay** (`rp_*`) | cube, **all six face normals prescribed** from the scheme-1 nominal path, interpolated to `N` steps | accuracy and constitutive cost, integrator-only |
| **free-standing** (`tx_*`, `ct_*`) | the ADR-92 probe: drained triaxial, 5 free DOF, `KrylovNewton` + `NormDispIncr` | robustness under a global Newton |
| **floor** (`fl_*`) | prescribed extension onto `p_min` from `p0 = 5 kPa`, zero free DOF | the ADR-93 ring regime |
| **cap control** (`ms_*`) | floor path with `-maxSubsteps 100` | does the `-maxSubsteps` seam reach scheme 2? |

---

## 2. The certificate — replayed strain path (integrator-only)

Deviation metric is G0's own `σ rel` = `‖σ_a − σ_b‖ / ‖σ_b‖`, taken as the
**maximum over the whole path**, not just the terminal step. `dEz` is the axial
strain increment; the path runs to `ε_z = 0.0182`.

| p0 | N | dEz | scheme | steps | ms/step | it/step | subME>0 | eta_end | p_end | q_end | rel.dev vs own-scheme fine | rel.dev s2-vs-s1 (same N) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 1820 | 1.0e-05 | 1 | 1820/1820 | 0.7 | 1.00 | 1820 | 1.4842 | 196.69 | 291.92 |  |  |
| 100 | 1820 | 1.0e-05 | 2 | 1820/1820 | 1.1 | 1.00 | 0 | 1.4848 | 196.78 | 292.18 |  | 1.331e-03 |
| 100 | 182 | 1.0e-04 | 1 | 182/182 | 5.5 | 1.00 | 182 | 1.4852 | 198.04 | 294.12 | 3.386e-02 |  |
| 100 | 182 | 1.0e-04 | 2 | 182/182 | 1.3 | 1.00 | 0 | 1.4867 | 197.11 | 293.04 | 7.789e-03 | 3.690e-02 |
| 100 | 40 | 4.6e-04 | 1 | 40/40 | 22.6 | 1.00 | 40 | 1.5277 | 221.74 | 338.77 | 1.915e-01 |  |
| 100 | 40 | 4.6e-04 | 2 | 40/40 | 1.7 | 1.00 | 0 | 1.4962 | 201.56 | 301.58 | 2.714e-02 | 1.890e-01 |
| 20 | 1820 | 1.0e-05 | 1 | 1820/1820 | 1.3 | 1.00 | 1820 | 1.8593 | 51.87 | 96.45 |  |  |
| 20 | 1820 | 1.0e-05 | 2 | 1820/1820 | 1.1 | 1.00 | 0 | 1.8597 | 51.73 | 96.21 |  | 2.911e-03 |
| 20 | 182 | 1.0e-04 | 1 | 182/182 | 11.4 | 1.00 | 182 | 1.8582 | 52.86 | 98.22 | 7.476e-02 |  |
| 20 | 182 | 1.0e-04 | 2 | 182/182 | 1.5 | 1.00 | 0 | 1.8583 | 50.86 | 94.52 | 2.037e-02 | 9.232e-02 |
| 20 | 40 | 4.6e-04 | 1 | 40/40 | 37.8 | 1.00 | 40 | 2.0373 | 130.12 | 265.11 | 1.617e+00 |  |
| 20 | 40 | 4.6e-04 | 2 | 40/40 | 3.8 | 1.00 | 1 | 1.8508 | 49.06 | 90.80 | 5.370e-02 | 6.398e-01 |

**C1 CONSISTENCY — PASS at both confinements.** At `Δε_z = 1e-5` the two schemes
differ by `1.3e-3` (`p0 = 100`) and `2.9e-3` (`p0 = 20`) at the worst point of the
path, and by `6.6e-4` / `2.6e-3` at the terminal step; terminal `η` differs by
`4.2e-4` / `2.0e-4`, i.e. the fourth significant figure. **They are integrating
the same model to the same limit.** (Both bars are met with an order of margin.)

**C2 NO-DEGRADATION — PASS, and by a wide margin in scheme 2's favour.** Each
scheme's own discretisation error against its own fine reference:

| `p0` | `Δε_z` | scheme 1 | scheme 2 | ratio s2/s1 |
|---|---|---|---|---|
| 100 | `1.0e-4` | `3.39e-2` | `7.79e-3` | **0.23** |
| 100 | `4.6e-4` | `1.92e-1` | `2.71e-2` | **0.14** |
| 20 | `1.0e-4` | `7.48e-2` | `2.04e-2` | **0.27** |
| 20 | `4.6e-4` | `1.62e+0` | `5.37e-2` | **0.033** |

At `p0 = 20 kPa` and `Δε_z = 4.6e-4` scheme 1's answer is **160 % wrong** in
stress-norm and its terminal `η = 2.037` exceeds `M^b = 1.929` — it is outside
its own bounding surface — while scheme 2 is `5.4e-2` off with `η = 1.851`
against `M^b = 1.929`. That is the single most consequential row in this memo:
**at the coarse end of the campaign's increment range the baseline, not the
candidate, is the one producing an inadmissible stress state.**

### 2.1 `η` against `M^b`

| p0 | N | eta_end s1 | eta_end s2 | d(eta)/eta | M^b (s1) | eta/M^b s1 | eta/M^b s2 | rel.dev terminal |
|---|---|---|---|---|---|---|---|---|
| 100 | 1820 | 1.48417 | 1.48479 | 4.242e-04 | 1.8762 | 0.7910 | 0.7914 | 6.558e-04 |
| 100 | 182 | 1.48516 | 1.48665 | 1.005e-03 | 1.8755 | 0.7919 | 0.7925 | 4.388e-03 |
| 100 | 40 | 1.52775 | 1.49620 | 2.065e-02 | 1.8630 | 0.8200 | 0.7986 | 9.781e-02 |
| 20 | 1820 | 1.85933 | 1.85970 | 2.000e-04 | 1.9999 | 0.9297 | 0.9298 | 2.608e-03 |
| 20 | 182 | 1.85825 | 1.85835 | 5.245e-05 | 1.9987 | 0.9297 | 0.9286 | 3.767e-02 |
| 20 | 40 | 2.03733 | 1.85083 | 9.154e-02 | 1.9294 | 1.0559 | 0.9239 | 6.398e-01 |

At the fine increment `η/M^b` agrees to four figures (`0.7910` vs `0.7914` at
`p0 = 100`; `0.9297` vs `0.9298` at `p0 = 20`). The path stops at `ε_z = 0.0182`,
short of peak `η`, so neither arm reaches the bounding surface: `η/M^b < 1`
everywhere except scheme 1's `Δε_z = 4.6e-4`, `p0 = 20` row (`1.0559`).

### 2.2 Cost per point (constitutive only, zero free DOF)

| `p0` | `Δε_z` | s1 ms/step | s2 ms/step | speed-up |
|---|---|---|---|---|
| 100 | `1.0e-5` | 0.7 | 1.1 | **0.64× (scheme 2 slower)** |
| 100 | `1.0e-4` | 5.5 | 1.3 | **4.2×** |
| 100 | `4.6e-4` | 22.6 | 1.7 | **13.3×** |
| 20 | `1.0e-5` | 1.3 | 1.1 | 1.2× |
| 20 | `1.0e-4` | 11.4 | 1.5 | **7.6×** |
| 20 | `4.6e-4` | 37.8 | 3.8 | **9.9×** |

Scheme 1's cost per point grows ~linearly with the increment (its substep count
does); scheme 2's is nearly flat, because a 19-unknown Newton costs the same
whatever the increment — until it fails, which is §4. **Scheme 2 is more
expensive only at increments so fine that `ModifiedEuler` needs almost no
substeps.**

### 2.3 The explicit fallback did not fire

`subME>0` counts steps whose last update ran `ModifiedEuler`. Scheme 2: **0 of
1820, 0 of 182, 0 of 40** at `p0 = 100`; **0, 0, 1 of 40** at `p0 = 20`. Scheme 1:
every step, by construction. **ADR-92 D3's 58–74 % does not reproduce on a
drained triaxial path at either confinement.**

---

## 3. The ADR-93 ring regime (`p → p_min`)

| tag | scheme | N | done | ms/step | subME max | steps with subME>0 | eta peak | p min | p end |
|---|---|---|---|---|---|---|---|---|---|
| fl_p5_s1_n160 | 1 | 160 | 160 | 14.03 | 1320 | 160 | 1.9325 | 0.0101 | 0.0101 |
| fl_p5_s1_n40 | 1 | 40 | 40 | 53.15 | 4045 | 40 | 1.8551 | 0.0101 | 0.0101 |
| fl_p5_s2_n160 | 2 | 160 | 160 | 2.63 | 7 | 85 | 1.9649 | 0.0101 | 0.0101 |
| fl_p5_s2_n40 | 2 | 40 | 40 | 6.46 | 1282 | 23 | 1.8893 | 0.0101 | 0.0101 |

| step | t | p s1 | p s2 | eta s1 | eta s2 | rel.dev | subME s1 | subME s2 |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.0000 | 5 | 5 | 0.0000 | 0.0000 | 0.000e+00 | 0.0 | 0.0 |
| 8 | 0.0500 | 3.0239 | 3.0207 | 0.5873 | 0.5561 | 1.447e-02 | 101.0 | 0.0 |
| 16 | 0.1000 | 2.3639 | 2.3579 | 0.8469 | 0.8228 | 1.167e-02 | 120.0 | 0.0 |
| 24 | 0.1500 | 1.8443 | 1.8363 | 1.0416 | 1.0225 | 1.069e-02 | 137.0 | 0.0 |
| 32 | 0.2000 | 1.4298 | 1.4203 | 1.2026 | 1.1872 | 1.116e-02 | 157.0 | 0.0 |
| 40 | 0.2500 | 1.0732 | 1.0584 | 1.3430 | 1.3310 | 1.675e-02 | 185.0 | 0.0 |
| 48 | 0.3000 | 0.71519 | 0.6964 | 1.4736 | 1.4654 | 2.814e-02 | 233.0 | 0.0 |
| 56 | 0.3500 | 0.40468 | 0.38959 | 1.6020 | 1.5979 | 3.820e-02 | 313.0 | 0.0 |
| 64 | 0.4000 | 0.17842 | 0.16827 | 1.7344 | 1.7348 | 5.679e-02 | 476.0 | 0.0 |
| 72 | 0.4500 | 0.041689 | 0.036682 | 1.8811 | 1.8907 | 1.182e-01 | 955.0 | 0.0 |
| 80 | 0.5000 | 0.0101 | 0.0101 | 0.0000 | 0.0000 | 0.000e+00 | 7.0 | 7.0 |
| 88 | 0.5500 | 0.0101 | 0.0101 | 0.0000 | 0.0000 | 0.000e+00 | 7.0 | 7.0 |
| 96 | 0.6000 | 0.0101 | 0.0101 | 0.0000 | 0.0000 | 0.000e+00 | 7.0 | 7.0 |
| 104 | 0.6500 | 0.0101 | 0.0101 | 0.0000 | 0.0000 | 0.000e+00 | 7.0 | 7.0 |
| 112 | 0.7000 | 0.0101 | 0.0101 | 0.0000 | 0.0000 | 0.000e+00 | 7.0 | 7.0 |
| 120 | 0.7500 | 0.0101 | 0.0101 | 0.0000 | 0.0000 | 0.000e+00 | 7.0 | 7.0 |
| 128 | 0.8000 | 0.0101 | 0.0101 | 0.0000 | 0.0000 | 0.000e+00 | 7.0 | 7.0 |
| 136 | 0.8500 | 0.0101 | 0.0101 | 0.0000 | 0.0000 | 0.000e+00 | 7.0 | 7.0 |
| 144 | 0.9000 | 0.0101 | 0.0101 | 0.0000 | 0.0000 | 0.000e+00 | 7.0 | 7.0 |
| 152 | 0.9500 | 0.0101 | 0.0101 | 0.0000 | 0.0000 | 0.000e+00 | 7.0 | 7.0 |
| 160 | 1.0000 | 0.0101 | 0.0101 | 0.0000 | 0.0000 | 0.000e+00 | 7.0 | 7.0 |

Reading:

- **On the descent** (`p` 5 → 0.04 kPa, steps 8–72) the two schemes track each
  other to `1.1–1.7e-2`, degrading to `1.2e-1` at the last step before the floor
  — the expected loss of significance as `p` approaches `p_min` and `‖σ‖ → 0`.
  Peak `η` `1.9325` (s1) vs `1.9649` (s2), 1.7 % apart at `N = 160`.
- **Scheme 2 needs no explicit fallback on the descent**: `subME = 0` on every
  one of the first 80 steps, while scheme 1 spends 101 → 955 substeps per update
  (max 1320 at `N = 160`, 4045 at `N = 40`). Scheme 2 is **5.3× / 8.2×** cheaper.
- **At the floor both schemes are identical and both are `ModifiedEuler`.** Once
  `p` is clamped at `p_min = 0.0101` the deviator collapses to zero, and
  `BackwardEuler_CPPM`'s `p < m_Pmin` branch takes over: its Newton is disabled
  by a literal `errFlag = 0` (`ManzariDafalias.cpp:2418` on this checkout; the
  ADR cites `:2264`, which has moved) and it calls `explicit_integrator`
  unconditionally. Scheme 2's `subME > 0` on 85 of 160 steps — **53 %, squarely
  inside P0's 58–74 % band** — and every one of those is a floor-pinned step.
  **So P0's measurement was right about the place it measured, and wrong as a
  general statement about scheme 2.**

---

## 4. Robustness under a global Newton — the refutation

The free-standing probe (5 free DOF, `KrylovNewton`, `NormDispIncr` 200 iters).
`t9`/`t7` are the test tolerances `1e-9` (the probe's own) and `1e-7` (a control,
to separate a solver artefact from a constitutive one).

| run | scheme | N | dEz | global tol | steps done | stalled@ | ms/step (passed) | wall of the FAILING step (s) | it/step | it max | subME max | eta_end |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `ct_p100_s1_n200_t7` | 1 | 200 | 1.0e-04 | 1e-07 | 200 | None | 23.9 | 0.00 | 7.61 | 13 | 62 | 1.5196 |
| `ct_p100_s1_n200_t9` | 1 | 200 | 1.0e-04 | 1e-09 | 182 | 183 | 61.4 | 1.00 | 15.90 | 34 | 60 | 1.4852 |
| `ct_p100_s1_n40_t7` | 1 | 40 | 5.0e-04 | 1e-07 | 40 | None | 200.4 | 0.00 | 11.32 | 21 | 282 | 1.5668 |
| `ct_p100_s1_n40_t9` | 1 | 40 | 5.0e-04 | 1e-09 | 40 | None | 394.3 | 0.00 | 21.43 | 84 | 281 | 1.5668 |
| `ct_p100_s2_n200_t9` | 2 | 200 | 1.0e-04 | 1e-09 | 176 | 177 | 30.2 | 133.75 | 8.94 | 23 | 0 | 1.4730 |
| `ct_p100_s2_n40_t7` | 2 | 40 | 5.0e-04 | 1e-07 | 3 | 4 | 141.0 | 32.92 | 7.33 | 8 | 402 | 0.5658 |
| `ct_p100_s2_n40_t9` | 2 | 40 | 5.0e-04 | 1e-09 | 2 | 3 | 13.5 | 20.59 | 9.50 | 11 | 0 | 0.3836 |
| `ct_p20_s1_n200_t7` | 1 | 200 | 1.0e-04 | 1e-07 | 200 | None | 81.0 | 0.00 | 9.87 | 17 | 135 | 1.8585 |
| `ct_p20_s1_n200_t9` | 1 | 200 | 1.0e-04 | 1e-09 | 200 | None | 132.9 | 0.00 | 16.65 | 37 | 128 | 1.8586 |
| `ct_p20_s1_n40_t9` | 1 | 40 | 5.0e-04 | 1e-09 | 40 | None | 770.7 | 0.00 | 22.90 | 33 | 534 | 1.9163 |
| `ct_p20_s2_n200_t7` | 2 | 200 | 1.0e-04 | 1e-07 | 108 | 109 | 140.6 | 129.85 | 8.79 | 57 | 0 | 1.6731 |
| `ct_p20_s2_n200_t9` | 2 | 200 | 1.0e-04 | 1e-09 | 107 | 108 | 492.3 | 12.39 | 9.74 | 23 | 0 | 1.6676 |
| `ct_p20_s2_n40_t7` | 2 | 40 | 5.0e-04 | 1e-07 | 9 | 10 | 19033.0 | 22.91 | 49.44 | 84 | 613 | 1.3019 |
| `ct_p20_s2_n40_t9` | 2 | 40 | 5.0e-04 | 1e-09 | 1 | 2 | 36.9 | 16.51 | 11.00 | 11 | 0 | 0.3857 |

- **Scheme 1 completes 7 of 8** arms (counting the `ct_p20_s1_n40_t7` re-run; the
  one exception is `p0 = 100`, `Δε_z = 1e-4`, `tol 1e-9`, at 182 of 200).
- **Scheme 2 stalls in 8 of 8** (counting the `ct_p100_s2_n200_t7` re-run at 194
  of 200). At `Δε_z = 4.6e-4` it reaches 1–9 steps of 40; at `Δε_z = 1e-4`,
  107–194 of 200. **Loosening the global tolerance
  does not rescue it** (`p0 = 100`, `N = 40`: 2 → 3 steps; `p0 = 20`, `N = 40`:
  1 → 9), so the cause is not the solver's test.
- **The cost of a failing step is the finding.** `ct_p100_s2_n200_t9` spends
  **133.75 s** in the single step it fails, against 30 ms for a step it passes —
  4400×. `ct_p20_s2_n200_t9`: 12.4 s. `ct_p100_s2_n40_t9`: 20.6 s.
  Scheme 1's failing step costs 1.0 s.
- The mechanism is §5.2's ladder: a global-Newton trial iterate at a coarse step
  proposes a strain increment the CPPM cannot return, and `BackwardEuler_CPPM`
  then recurses (`:2538`, `implicitLevel++`, cap 10 at `:2352-2357`) through up
  to `2^9` halvings, each a 19-unknown Newton of up to 30 iterations, before
  falling through to `explicit_integrator`. Nothing reports that this happened.
- **The replay arm is the control that makes this unambiguous.** On the *same*
  path at the *same* `Δε_z = 4.6e-4`, with the increment given rather than
  proposed by a global Newton, scheme 2 completes 40 of 40 in 1.7 ms/step. The
  failure is on **off-path trial iterates**, not on the solution path.

Two runs (`ct_p100_s2_n200_t7`, `ct_p20_s1_n40_t7`) died with no traceback, no
CSV and no summary during the first sweep. **Both re-ran clean in isolation
(exit 0, 194/200 and 40/40).** Recorded as non-reproducible, attributed to
machine contention, **not** claimed as a finding.

---

## 5. Read-only code study (for phases b/d)

All line numbers verified on this checkout (`634824e1f`).

### 5.1 Does `TanType 2` give an algorithmic tangent under scheme 2? — **YES**, with two caveats

`LadrunoSANISAND3D::getTangent()` (`LadrunoSANISAND3D.cpp:170-178`, shadowing the
identical `ManzariDafalias3D.cpp:134-141`) returns `mCep_Consistent` for
`mTangType == 2`. Under scheme 2 `mCep_Consistent` is written at
`ManzariDafalias.cpp:2609` (`Cep_Consistent = aCepConsistent;`) from
`NewtonIter2(Delta0, InVariants, Delta, aCepConsistent)` (`:2457`), which fills it
in `NewtonSol` as the condensation of the 19×19 CPPM Jacobian:

```cpp
// ManzariDafalias.cpp:3448-3472
DSigma        = aC * DSigma;
if (DSigma.Invert(CSigma) != 0) { ... return -1; }
else  CSigma = CSigma * aC;
...
    Cep            = -1.0 * CSigma;
```

That is a genuine **algorithmic (consistent) tangent** — `∂σ_{n+1}/∂ε_{n+1}` of the
converged return map, not a continuum tangent.

**Caveat 1 — it is one iterate stale.** `NewtonIter2`'s loop (`:3052-3115`) tests
convergence *before* calling `NewtonSol`, so on the iterate that converges
`NewtonSol` is never called again: the tangent handed back was evaluated at the
**second-to-last** iterate. Harmless at tight `TolR`, not exact.

**Caveat 2 — it silently degrades.** Every path out of the CPPM that reaches
`explicit_integrator` (`:2431-2436` low-`p`; `:2584-2588` ladder exhaustion)
**overwrites** `aCepConsistent` with `ModifiedEuler`'s substep-chained product
(`:1835`, `aCep_Consistent = aCep_thisStep * (aD * aCep_Consistent + T * mIImix)`).
So `TanType 2` returns an algorithmic tangent on steps where the implicit return
succeeded and a chained explicit tangent on steps where it did not, **with no
diagnostic distinguishing them**.

**`TanType 2` is not scheme-2-only.** `ModifiedEuler` maintains its own
`aCep_Consistent` (`:1490`, `:1835`), so `TanType 2` is meaningful on scheme 1
too — the deck (`sanisand_tau0_band.py:349`) already uses it. The two objects are
different: scheme 1's is a product of continuum elastoplastic tangents over the
substep chain; scheme 2's is the return map's own Jacobian. **Under `-implex`
both are inert** — the material hands out `Ce(p_n)` and says so
(`LadrunoSANISAND.cpp:2097-2099`: *"TanType … is INERT under -implex"*).

### 5.2 What happens on CPPM non-convergence? — **it always succeeds, by falling back**

```cpp
// ManzariDafalias.cpp:2345
int errFlag = 1, SchemeControl = 2, mMaxSubStep = 10;
//   0 : newton did not converge in MaxIter iterations
//  -1 : the jacobian is singular
//  -2 : converged stress has p < 0
//  -3 : max number of sub-stepping reached
```

`NewtonIter2` (`:2457`) returns 0 / −1; `Check()` may return −2 / −4. Then
(`:2472-2586`):

```cpp
while(errFlag != 1)
{
    if (errFlag == -1) SchemeControl = 3; // do an explicit integration
    if (errFlag == -2) SchemeControl = 2; // do sub-stepping
    ...
    } else if (SchemeControl == 2) {
        implicitLevel++;                                    // :2538
        nStrain = cStrain + StrainInc / 2;
        errFlag = BackwardEuler_CPPM(..., implicitLevel);   // recursion, first half
        if (errFlag == -3) { SchemeControl += 1; continue; }
        ...                                                 // then second half
    } else {
        explicit_integrator(...);                           // :2584
        errFlag = 1;                                        // :2588
    }
}
```

The recursion gives up at `implicitLevel > mMaxSubStep` (`:2352-2357`, `return
-3`); the caller then bumps `SchemeControl` to 3 and integrates the **whole**
increment explicitly. The header default is `implicitLevel = 1`
(`ManzariDafalias.h:331`) and `integrate()` passes no argument, so the tree is up
to 9 levels deep — **up to 512 half-increments, each a 19-unknown Newton of up to
30 iterations with a 19×19 solve** — before the fallback. That is §4's 12–134 s.

**The loop cannot exit with `errFlag != 1`.** Its terminal branch always assigns
1. The same is true of the low-`p` branch (`:2418` `errFlag = 0;` → `:2431`
explicit → `:2436` `errFlag = 1;`).

### 5.3 Does any of that reach the element as a refusal? — **NO**

```cpp
// ManzariDafalias.cpp:1022-1027
        if ((mScheme == INT_BackwardEuler))
            BackwardEuler_CPPM(mSigma_n, mEpsilon_n, mEpsilonE_n, mAlpha_n, ...);
```

**`integrate()` discards the return value** — there is no assignment. Combined
with §5.2 (it is always 1 anyway) and `debugFlag = false` (`:57`), a CPPM
non-convergence is invisible in every channel the fork has: no return code, no
`opserr` line, no response.

The fork's only trial-time refusal for this material is

```cpp
// LadrunoSANISAND.cpp:3993-3996
int
LadrunoSANISAND::ladrunoUpdateStatus(void) const
{
    return mSubstepCapHitInME ? LADRUNO_MATERIAL_REFUSED : 0;
}
```

and `mSubstepCapHitInME` is set at exactly one site, inside `ModifiedEuler()`
(`ManzariDafalias.cpp:1578-1600`), under `mMaxSubstepsInME`. So on scheme 2 a
refusal can only arise **after** the CPPM has already given up and fallen back to
`ModifiedEuler`, and only if `-maxSubsteps > 0`. The F7 refusal roster
(`LEDGER_quirks.md`, "element refusal roster", PR #838) is therefore irrelevant
to a CPPM failure: there is nothing for the element to forward.

### 5.4 A defect found on the way: the `-maxSubsteps` inertness warning is WRONG for scheme 2

`LadrunoSANISAND::schemeReachesModifiedEuler()` (`:1035-1048`) returns **false**
for `mScheme == 2`, so the constructor prints

> `WARNING LadrunoSANISAND tag 1: -maxSubsteps N has NO EFFECT with IntScheme 2.`
> `The seam it sets (ManzariDafalias mMaxSubstepsInME) is read at exactly one site,`
> `inside ModifiedEuler(), and this scheme does not route there.` — `:1206-1213`

(and the identical claim for `-honorTolR`, `:1187-1199`). **Scheme 2 does route
there**: `explicit_integrator`'s `switch` (`:1070-1101`) does not enumerate
`INT_BackwardEuler`, so it falls to `default:` → `ModifiedEuler`, which reads both
seams. **Measured (§6): a `-maxSubsteps 100` cap on scheme 2 changed a run that
completed 40 of 40 uncapped into one that refuses at step 18.** The warning is
false, and it contradicts the guide's own §3, which *requires* `-maxSubsteps > 0`
on scheme 2 under `-implex` and refuses the deck without it
(`LadrunoSANISAND.cpp:2047-2058`). One of the two has to move; the measurement
says it is the warning.

---

## 6. Control — does `-maxSubsteps` reach scheme 2?

| run | scheme | cap | uncapped result | capped result |
|---|---|---|---|---|
| `ms_fl_s2_n40_cap100` | 2 | 100 | 40/40 (max 1282 substeps) | **refuses at step 18** |
| `ms_fl_s1_n40_cap100` | 1 | 100 | 40/40 (max 4045 substeps) | **refuses at step 1** |

Both arms refuse, through `LADRUNO_MATERIAL_REFUSED` → `LadrunoBrick` →
`Domain::update` → step cut. The seam is live on scheme 2. §5.4 stands.

---

## 7. Phase (b)/(c) — the bearing BVP

> [!failure] **PHASE (b)/(c) — REFUTED, and not narrowly.**
> On the CP1/ADR-95 bearing deck scheme 2 committed **11 steps to `s/B = 4e-5` in
> 1347 s**, against the scheme-1 baseline's **51 steps to `s/B = 0.019` in 1267 s**
> — **475× shallower for the same wall clock**. Its step size collapsed to 25× the
> subdivision floor (`ds = 0.005 mm` against the baseline's `2.5 mm`, `free =
> False`), **every one of its 11 committed steps needed the relaxed rung 3**
> (100 % past rung 1, against CP1's 61–83 %), and it spent 4 of 80 subdivisions
> doing it. **The 1 % load–settlement bar could not be evaluated: there is no
> overlap in `s/B` between the arms at all.**

**Deck.** `sanisand_tau0_band.py`, leg `h1.0_e0.6944`, purpose-sized domain
`--xlim 10 --zbot 8` (624 DOF, 84 hexes — the `adr92_deck/x10z8` configuration),
`--maxsubsteps 20000`, `NormUnbalance @ 1e-5·γV`, `LoadControl` (stock, no
`-tangentPredictor`), 1200 s budget per arm, arms run **sequentially** so the wall
column means something. `INT_SCHEME` is a module constant with no CLI or env seam
(`:348`, read at call time at `:696` and `:1090`); `f12_bvp.py` loads the driver as
a module and sets it, so **no repo file was edited** and every one of the driver's
own controls, assertions and JSON provenance still ran. All three arms cleared
the geostatic controls identically (`resultant 4.44e-16`, `patch 1.89e-14`,
`η/M_c 0.6425`, `OutsideBounding 0`, `CLAMPING 0`).

### 7.1 Leg outcome (x10z8, `h1.0_e0.6944`, 624 DOF, 84 hexes, 1200 s budget each)

| arm | scheme | TanType | steps | mode | s/B reached | q at end (kPa) | wall s | s per step | nfail | nsub | nrelax | CLAMPING | OutsideBounding |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| IntScheme 1, TanType 0 | 1 | 0 | 52 | WALL | 0.02155 | 668.78 | 1209 | 23.24 | 50 | 1 | 15 | 0 | 0 |
| IntScheme 1, TanType 2 | 1 | 2 | 51 | WALL | 0.01905 | 593.02 | 1267 | 24.84 | 61 | 1 | 22 | 0 | 0 |
| IntScheme 2, TanType 2 | 2 | 2 | 11 | WALL | 0.00004 | 2.29 | 1347 | 122.50 | 34 | 4 | 11 | 0 | 0 |

### 7.2 Committed load-settlement at matched s/B

| s/B | IntScheme 1, TanType 0 | IntScheme 1, TanType 2 | IntScheme 2, TanType 2 |
|---|---|---|---|
| 0.001 | 34.19 | 34.76 | -- |
| 0.002 | 64.53 | 65.44 | -- |
| 0.005 | 153.55 | 155.17 | -- |
| 0.010 | 303.38 | 305.86 | -- |
| 0.015 | 458.34 | 462.52 | -- |
| 0.020 | 618.30 | -- | -- |

**scheme 2 vs scheme 1 (both TanType 2) at matched s/B — the 1 % bar:**

| s/B | q s1_T2 | q s2_T2 | rel. diff | inside 1 %? |
|---|---|---|---|---|
| 0.001 | 34.76 | -- | -- | (not reached by both) |
| 0.002 | 65.44 | -- | -- | (not reached by both) |
| 0.005 | 155.17 | -- | -- | (not reached by both) |
| 0.010 | 305.86 | -- | -- | (not reached by both) |
| 0.015 | 462.52 | -- | -- | (not reached by both) |
| 0.020 | -- | -- | -- | (not reached by both) |

### 7.3 Wall clock to matched s/B

| s/B | IntScheme 1, TanType 0 | IntScheme 1, TanType 2 | IntScheme 2, TanType 2 |
|---|---|---|---|
| 0.001 | 27 s | 44 s | not reached |
| 0.002 | 68 s | 102 s | not reached |
| 0.005 | 258 s | 278 s | not reached |
| 0.010 | 549 s | 633 s | not reached |
| 0.015 | 842 s | 883 s | not reached |
| 0.020 | 1150 s | not reached | not reached |

### 7.4 Ladder decomposition (`adr92_bvp_gate.py`'s own arithmetic)

| arm | steps | rung1 | rung2 | rung3 | past rung 1 % | failed-rung iteration share % |
|---|---|---|---|---|---|---|
| IntScheme 1, TanType 0 | 52 | 20 | 17 | 15 | 61.5 | 85.4 |
| IntScheme 1, TanType 2 | 51 | 15 | 14 | 22 | 70.6 | 88.2 |
| IntScheme 2, TanType 2 | 11 | 0 | 0 | 11 | 100.0 | 95.7 |

**Reading.**

1. **The two scheme-1 arms are a working control.** `TanType 0` and `TanType 2`
   agree to **0.8–1.7 %** at every matched `s/B` and reproduce CP1's ladder rates
   (61.5 % / 70.6 % past rung 1, 85.4 % / 88.2 % failed-rung iteration share
   against CP1's 61–83 % / 89–93 %). So the reduced deck is representative and the
   instrument is sound. *(Aside: on this 84-element deck `TanType 0` is not the
   trap the driver's comment describes — it is marginally **faster** to every
   checkpoint. The driver's measurement was at `h0 = 0.25`, where the linear
   solve is not negligible; it is not contradicted, only bounded.)*
2. **Scheme 2 never leaves the subdivision floor.** `relaxed = 1` on every
   committed row from the first; `ds` fell from `0.02 mm` to `0.005 mm` and stayed
   there. This is §4's mechanism on a real BVP: the global Newton's trial iterates
   are off-path, the CPPM cannot return them, it grinds its recursive-halving
   ladder, the rung fails, the controller halves `ds`, and the smaller step does
   not help because the *iterate*, not the step, is what the CPPM chokes on.
   122.5 s per committed step against the baseline's 23–25 s.
3. **`base_foot_mismatch` is 7.9e-2 on the scheme-2 arm against 2–3e-4 on the
   others.** Do **not** read this as an equilibrium defect: the arm is at
   `q = 2.29 kPa`, a factor ~260 below the others, so the same absolute residual
   is a far larger relative one. It is reported because the driver reports it, not
   as a finding.
4. **`-implex` was OFF in all three arms.** Whether scheme 2 is the better
   *companion* — the question §9.2 leaves open — is **not** answered here and
   cannot be: the companion runs at `commitState` on an increment that is already
   given, which is §2's replay regime (where scheme 2 wins), not this one.

---

## 8. What this licenses — and what it does not

**MAY say:**

- Scheme 2 integrates the same model to the same limit as scheme 1 (§2, C1).
- At `Δε_z ≥ 1e-4` on a given strain path scheme 2 is 3.7–4.3× more accurate and
  4.2–7.6× cheaper than scheme 1; at `4.6e-4`, 7–30× more accurate and 10–13×
  cheaper (§2).
- At `p0 = 20 kPa`, `Δε_z = 4.6e-4`, scheme 1 leaves its own bounding surface
  (`η/M^b = 1.056`) and is 160 % wrong in stress norm; scheme 2 is not (§2).
- ADR-92 D3's premise ("58–74 % of scheme 2's calls integrate explicitly") holds
  **only** on floor-pinned steps and is **0 %** everywhere else on these paths
  (§2.3, §3).
- `TanType 2` under scheme 2 is a genuine algorithmic tangent, one iterate stale,
  silently replaced by the explicit chained tangent whenever the CPPM falls back
  (§5.1).
- A CPPM failure is **invisible**: return value discarded, `debugFlag` off, no
  refusal path (§5.3).
- `-maxSubsteps` **does** work on scheme 2 and the class's warning to the contrary
  is false (§5.4, §6).

**MAY say (phase b):**

- On the CP1/ADR-95 bearing deck, `IntScheme 2` is **not usable as the primary
  integrator**: 475x shallower than the baseline for the same wall clock, 122.5 s
  per committed step against 23-25 s, `ds` pinned at the subdivision floor, 100 %
  of steps past rung 1 (§7).
- The two scheme-1 arms (`TanType 0` and `TanType 2`) agree to 0.8-1.7 % at every
  matched `s/B` and reproduce CP1's ladder rates, so the reduced deck and the
  instrument are sound (§7).

**MUST NOT say:**

- That scheme 2 is faster *in general*. It is slower at `Δε_z = 1e-5`
  (`p0 = 100`, 0.64×) and its failing steps cost 12–134 s (§4).
- That scheme 2 is more robust. Under a global Newton at the campaign increment
  it is measurably **less** robust (§4).
- Anything about plane strain, `LadrunoUP`, cyclic/reversal paths, `-implex` +
  scheme 2, or parallel. None were run. In particular **the companion question is
  open**: `-implex` was OFF in every arm here, and the companion runs at
  `commitState` on an increment that is already given — §2's regime, where
  scheme 2 wins — not §4/§7's.
- That the phase-(b) refutation would survive `-tangentPredictor`, a
  displacement-controlled push, or a per-iterate cap on the CPPM ladder. None
  were tried.
- Anything about a capacity, plateau or limit point. Single Gauss point in §2–§6.

---

## 9. The exact ledger and guide text this would add

### 9.1 `Ladruno_implementation/LEDGER_quirks.md` — two new rows

> **`ManzariDafalias::integrate()` discards `BackwardEuler_CPPM`'s return value,
> and the CPPM can never fail anyway.** `ManzariDafalias.cpp:1023-1027` calls the
> implicit return without an assignment. It does not matter, because the
> `while(errFlag != 1)` ladder at `:2472-2586` always terminates by calling
> `explicit_integrator` and setting `errFlag = 1` (`:2584-2588`), as does the
> low-`p` branch (`:2418` `errFlag = 0` → `:2431` explicit → `:2436` `= 1`). With
> `debugFlag` a compile-time `false` (`:57`), an `IntScheme 2` step whose Newton
> diverged, whose Jacobian was singular, or which recursed through up to 512
> half-increments (`:2352-2357`, `:2538`) and then gave up, is **indistinguishable
> from a clean implicit return** in every channel: no code, no log line, no
> response. The only observable is `mSubstepsTakenInME` via the `substeps`
> response — non-zero iff `ModifiedEuler` ran — and that is silent when the ladder
> succeeds by substepping rather than by falling through. Measured cost of an
> invisible failure on a single-element drained triaxial at `Δε_z = 1e-4`:
> **133.75 s for one step against 30 ms for its neighbours** (WP-105 / F12).

> **`-maxSubsteps` / `-honorTolR` are NOT inert on `IntScheme 2`, but
> `LadrunoSANISAND` says they are.** `schemeReachesModifiedEuler()`
> (`LadrunoSANISAND.cpp:1035-1048`) returns false for scheme 2, so the
> constructor prints "`-maxSubsteps N has NO EFFECT with IntScheme 2`"
> (`:1206-1213`) — while `LadrunoSANISAND_implex_guide` §3 and the `-implex`
> parser (`:2047-2058`) *require* `-maxSubsteps > 0` on that very scheme. The
> warning is wrong: `explicit_integrator`'s switch (`:1070-1101`) does not
> enumerate `INT_BackwardEuler`, so scheme 2 falls to `default: ModifiedEuler`
> whenever the CPPM falls back, and both seams are read there. **Measured
> (WP-105 / F12): on the `p → p_min` path a `-maxSubsteps 100` cap turned a
> scheme-2 run that completed 40 of 40 uncapped into one that refuses at step 18.**
> Fix is one line in `schemeReachesModifiedEuler()` (`s == 2` must return true);
> until then, read the warning as false and the guide as authoritative.

### 9.2 `Ladruno_implementation/92_ladruno_sanisand_implex_adr.md` — amend D3

> **D3, amended by WP-105 / F12 (2026-09-15).** The *decision* stands — scheme 1
> with `-maxSubsteps` is the companion default. **The stated reason does not.**
> D3 says scheme 2 "is not an implicit return where the campaign's problem lives"
> because 58–74 % of its calls take the low-`p` branch. Measured on
> `634824e1f`: **0 of 1820 steps** on a replayed drained-triaxial path at
> `p0 = 100` and `20 kPa`, and **0 of 80** on the descent of a `p → p_min` path;
> the 58–74 % reproduces (53 %) **only once the point is pinned at `p_min` with a
> zero deviator**, i.e. on steps where nothing is integrated. On a *given* strain
> increment — which is exactly what the commit-time companion receives — scheme 2
> is **3.7–4.3× more accurate and 4.2–7.6× cheaper** than scheme 1 at the campaign
> increment.
>
> **The real reason to keep scheme 1 as the primary integrator is one D3 never
> measured:** under a global Newton at `Δε_z ≥ 1e-4` scheme 2 stalls where scheme 1
> does not (8 of 8 arms vs 1 of 8) and each failing step burns 12–134 s in the
> recursive-halving ladder (`ManzariDafalias.cpp:2538`, up to 512 half-increments)
> before falling back to `ModifiedEuler` **silently** (`integrate()` discards the
> return value at `:1023-1027`; `debugFlag` is compiled off at `:57`). Loosening the
> global tolerance does not rescue it. **On the CP1/ADR-95 bearing leg
> (`x10z8`, `h1.0_e0.6944`, 1200 s) scheme 2 reached `s/B = 4e-5` against the
> baseline's `0.019` — 475x shallower for the same wall clock — with `ds` pinned at
> 25x the subdivision floor and 100 % of its steps on the relaxed rung 3.** That is
> the measurement D3 should have carried, and it supports D3's conclusion far more
> strongly than D3's own argument did.
>
> **Whether scheme 2 is the better COMPANION remains open and is worth a
> measurement**, because the companion runs at `commitState` on an increment that
> is already given, with no global Newton to feed it off-path iterates — the regime
> in which it is 3.7–4.3x more accurate and 4.2–7.6x cheaper. `-implex` was OFF in
> every arm of WP-105 / F12.

### 9.3 `Ladruno_implementation/LadrunoSANISAND_implex_guide.md` §3 — replace the scheme-2 bullet

> - **Scheme 2 (`BackwardEuler_CPPM`).** *Permitted* but not the default. Same
>   `-maxSubsteps > 0` requirement, refused the same way — and note that the
>   constructor's "`-maxSubsteps` has NO EFFECT with IntScheme 2" warning is
>   **false** (WP-105 / F12; `schemeReachesModifiedEuler()` is wrong for `s == 2`,
>   see `LEDGER_quirks`), because the CPPM's fallback routes through
>   `explicit_integrator`'s `default:` to `ModifiedEuler`, where both seams live.
>   What scheme 2 buys, measured at a material point on a given strain increment:
>   **3.7–4.3× less discretisation error and 4.2–7.6× less wall time** than
>   scheme 1 at `Δε_z = 1e-4`, rising to 7–30× / 10–13× at `4.6e-4`; the low-`p`
>   explicit fallback fires on **0 %** of steps until the point is pinned at
>   `p_min`, where it fires on ~53 % and both schemes become the same operator.
>   What it costs: a CPPM step that cannot return recurses through up to 512
>   half-increments before falling back, **silently** — `integrate()` discards the
>   return value, `debugFlag` is compiled off, and nothing refuses. Under a global
>   Newton at `Δε_z ≥ 1e-4` that is measured at 12–134 s for a single failing step
>   and at outright stalls where scheme 1 completes; on the CP1/ADR-95 bearing leg
>   it reached `s/B = 4e-5` in 1347 s against the baseline's `0.019` in 1267 s, with
>   `ds` pinned at 25x the subdivision floor. **Use it where the increment is given
>   (a companion return, a prescribed-strain probe); do not make it the primary
>   integrator of a load- or displacement-controlled BVP without a cap on the
>   ladder.**

---

## Log

- 2026-09-16 — WP-105 / F12 phases (a) and (b)/(c) run and written. Build `634824e1f`, no
  source edited. 16 free-standing + 12 replay + 4 floor + 18 control
  material-point runs, plus the phase-(b) bearing arms in §7.
