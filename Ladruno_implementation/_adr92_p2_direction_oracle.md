---
title: "ADR 92 / P2 — does correcting the FLOW DIRECTION fix the corner error? (oracle, LANE E)"
project: Ladruno
type: results
status: "NO — the crossing is the QUIETEST step on the path (A's implexError there is 500x below its own path max), the trial direction does not even fix the SIGN, and it forfeits the frozen tangent: KEEP operator A."
priority: high
owner: nmora
related:
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[_adr92_p0_oracle_results]]"
  - "[[_adr92_p1_redblue_review]]"
tags: [adr, implex, sanisand, oracle, measurement, non-associated-flow]
updated: 2026-09-06
---

# ADR 92 / P2 — the direction question

> [!danger] **VERDICT — keep A.** The premise is right; the conclusion does not follow. A's
> extrapolated volumetric plastic increment **does** carry the wrong sign across phase transformation
> at every resolved crossing — and it costs nothing, because `D -> 0` there so the wrong-signed
> quantity is itself vanishing: at `p0 = 100 kPa`, `d eps_z = 1e-5` A's `implexError` at the crossing
> is **2.80e-6** against a path max of **1.39e-3**, the *quietest* step by 500x. **B fixes neither**:
> its volumetric increment has the **same wrong sign** as A's at every resolved crossing (`p100`,
> `1e-5`: A `+3.73e-10`, B `+9.24e-11`, true `-2.79e-10`), because `D`'s sign is set by the back
> stress `alpha` against `sqrt(2/3)*aD` and `alpha` is lagged in *every* variant whatever stress `R`
> is read at. B is 1.8x-22x worse on path-max error, misses the terminal `q/p` by 9-54 % where A
> matches the implicit to four figures, and forfeits the frozen tangent: **0.44-2.29** vs A's 1e-10.

## 1. What was run

`--gate GE` of `adr92_p0_oracle/sanisand_implex_oracle.py` (numpy 2.4.6, `python3.12`): the G2
protocol — drained TX compression from the P0 seeds (`tx_p100_e0.6944_s1_n40`, `tx_p5_…`), companion
scheme 1, to `eps_z = 0.02`, increments `{1e-5, 5e-5, 1e-4, 5e-4, 1e-3}` (**P0's G2 set was
`1e-4 … 1e-2`;** the window moves down a decade because the crossing must be *resolved* before its
error can be attributed). `kX` = first step whose committed `D` changes sign. **B** (oracle
`implex_T`): `g~ = f*|d_eps_p(n)|` (tensor norm of the committed plastic strain increment, *not*
`mDGamma`), `d_eps_p~ = g~ R_tr/|R_tr|`, `R_tr = B n - C (n^2 - I/3) + D I/3` at the elastic trial
stress `sigma_n + Ce(p_n):d_eps` with committed `alpha_n, z_n, e_n, alpha_in_n`. **C** (`implex_H`):
deviator from history as A, volumetric part `ONE3 * g~ * D_tr/|R_tr|`. Stress form and (absent) clamp
as A; the oracle's own form `"B"` is P0's `dGamma` variant, untouched. **A is byte-identical** —
`--gate G0`/`G4` before and after the edit are `diff`-clean.

## 2. The crossing, per increment — `implexError` max over path / at `kX`

| `p0` | `d eps_z` | `kX` | `ez(kX)` | A | **B** | C | `q/p` A / B / C (implicit) |
|---|---|---|---|---|---|---|---|
| 100 | 1.0e-5 | 215 | 2.16e-3 | **1.39e-3 / 2.80e-6** | 3.09e-2 / 4.00e-5 | 1.40e-3 / 2.67e-6 | **1.5500** / 1.4041 / 1.5500 (1.5500) |
| 100 | 5.0e-5 | 43 | 2.20e-3 | **1.37e-2 / 7.45e-5** | 6.86e-2 / 1.03e-4 | 1.37e-2 / 7.14e-5 | **1.5500** / 1.1070 / 1.5500 (1.5500) |
| 100 | 1.0e-4 | 20 | 2.10e-3 | **5.68e-2 / 2.83e-4** | 1.34e-1 / 2.95e-4 | 5.67e-2 / 2.89e-4 | **1.5514** / 1.2594 / 1.5514 (1.5516) |
| 100 | 5.0e-4 | 1 † | 1.00e-3 | **1.87e-1 / 1.86e-1** | 4.26e-1 / 1.37e-1 | 1.90e-1 / 1.90e-1 | **1.6013** / 0.7329 / 1.6012 (1.5968) |
| 5 | 1.0e-5 | 34 | 3.50e-4 | **1.18e-2 / 9.22e-5** | 7.06e-2 / 1.04e-4 | 1.25e-1 / 8.91e-5 | **2.0155** / 1.1968 / 1.3983 (2.0154) |
| 5 | 5.0e-5 | 4 | 2.50e-4 | **1.29e-1 / 1.24e-2** | 2.37e-1 / 1.03e-2 | 3.58e-1 / 1.24e-2 | **2.0159** / 1.1949 / 1.8704 (2.0159) |
| 5 | 1.0e-4 | 1 † | 2.00e-4 | **1.82e-1 / 1.82e-1** | 3.87e-1 / 1.38e-1 | 5.25e-1 / 1.86e-1 | **2.0175** / 1.7373 / 1.3278 (2.0173) |

† crossing inside step 1-2, so `5e-4`/`1e-3` cannot resolve it; at `p100`, `1e-3` the committed `D` never
changes sign at all on A's or C's path (max 2.02e-1 / 7.06e-1 / 5.69e-1). **The two `p0 = 5` breakdown
rows are omitted:** there B's and C's max error is 50x and 230x *below* A's (6.6e-1, 1.6 vs 3.3e+1,
3.6e+2) — the only place B wins — but the *implicit companion* has itself collapsed there
(`q/p = 0.0000` at `1e-3`), so they measure breakdown, not accuracy.

## 3. The table that decides it

| `p0` | `d eps_z` | A at `kX` | B at `kX` | **A/B** | A max | B max | **A/B (max)** |
|---|---|---|---|---|---|---|---|
| 100 | 1.0e-5 | 2.80e-6 | 4.00e-5 | **0.07** | 1.39e-3 | 3.09e-2 | **0.045** |
| 100 | 5.0e-5 | 7.45e-5 | 1.03e-4 | **0.73** | 1.37e-2 | 6.86e-2 | **0.20** |
| 100 | 1.0e-4 | 2.83e-4 | 2.95e-4 | **0.96** | 5.68e-2 | 1.34e-1 | **0.43** |
| 100 | 5.0e-4 | 1.86e-1 | 1.37e-1 | 1.35 | 1.87e-1 | 4.26e-1 | **0.44** |
| 5 | 1.0e-5 | 9.22e-5 | 1.04e-4 | **0.89** | 1.18e-2 | 7.06e-2 | **0.17** |
| 5 | 5.0e-5 | 1.24e-2 | 1.03e-2 | 1.21 | 1.29e-1 | 2.37e-1 | **0.55** |
| 5 | 1.0e-4 | 1.82e-1 | 1.38e-1 | 1.32 | 1.82e-1 | 3.87e-1 | **0.47** |

`A/B < 1` = A better. **B never wins on a non-breakdown path** (`A/B(max)` 0.045-0.55); where it wins at the crossing *step* (1.2-1.35x) the crossing is unresolved and `q/p` is 14-54 % off.

## 4. The sign — `volp~ = tr(d_eps_p~)` against the companion's `tr(d_eps_p)`

| path | A | **B** | C | true | verdict |
|---|---|---|---|---|---|
| `p100` 1e-5 `kX` 215 | +3.73e-10 | +9.24e-11 | +8.54e-11 | −2.79e-10 | **all WRONG** |
| `p100` 5e-5 `kX` 43 | +8.78e-9 | +1.07e-9 | +1.14e-9 | −7.54e-9 | **all WRONG** |
| `p5` 1e-5 `kX` 34 | +2.48e-9 | +5.91e-10 | +5.91e-10 | −1.37e-9 | **all WRONG** |
| `p100` 5e-4 `kX` 1 † | +5.55e-6 | −7.90e-7 | −7.72e-7 | −2.34e-6 | A wrong, B/C right |
| `p5` 1e-4 `kX` 1 † | +8.95e-7 | −1.49e-7 | −1.48e-7 | −5.25e-7 | A wrong, B/C right |

**On every resolved crossing all three get the sign wrong**; B recovers it only on the two unresolved rows, where the step is coarse enough that the elastic trial lands past the crossing. The magnitudes are `1e-10`-`1e-9`, four to five orders below the same increment's deviator — hence invisible in `implexError`.

## 5. How non-linear a B step is

`J` = numerical secant `d sigma~/d eps` (`h = 1e-9`); `Ce` = the frozen `Ce(p_n)` the C++ would hand the element; chord = how far the *whole* step is from `J`-linear. Entries A / **B** / C.

| `p0` | `d eps_z` | warm | `max\|J-Ce\|/max\|Ce\|` | asym(`J`) | chord non-lin |
|---|---|---|---|---|---|
| 100 | 1e-4 | 20 | 9.7e-11 / **6.98e-1** / 2.7e-4 | 8.3e-11 / 7.6e-4 / 6.1e-6 | 2.4e-10 / **4.10e-2** / 3.5e-6 |
| 5 | 1e-4 | 2 | 2.8e-11 / **4.41e-1** / 1.4e-4 | 4.0e-11 / 1.1e-4 / 4.4e-7 | 7.1e-11 / **2.00e-2** / 5.7e-7 |
| 5 | 1e-4 | 20 | 1.3e-10 / **2.29e+0** / 3.0e-4 | 1.9e-11 / 1.6e-3 / 4.9e-6 | 1.2e-10 / **7.79e-1** / 7.4e-6 |

A is exact by construction. **B's true tangent differs from the operator the element would be given by
44-229 % of `max|Ce|`, and its step is 2-78 % non-linear** (`p100`, `5e-4`, warm 4, omitted: 7.40e-1 /
2.75e-2). The asymmetry is small (`1e-4`-`1.6e-3`): B does not cost the *symmetric* solver, it costs
the tangent being right at all — visible as robustness, since B's step broke the oracle's lateral
secant 8 times in 2000 (`p100`, `1e-5`) where A's never broke.

## 6. Interpretation

The premise is exactly right and its consequence is not. Non-associated flow does make A carry a
stale direction, and at phase transformation that direction's volumetric part is stale in *sign*, not
merely in size (§4) — but phase transformation is *defined* by `D = 0`, so the wrong-signed quantity
is `O(1e-10)` while the deviator of the same increment is `O(1e-5)`: the crossing is the one place on
the path where the extrapolation is least wrong (§2). Nor does evaluating `R` at the elastic trial
stress repair the sign, for a structural reason — `D`'s sign is `sqrt(2/3)*aD - alpha:n`, the back
stress against the dilatancy surface, and `alpha` is committed in every IMPL-EX variant, so moving the
stress at which `R` is read cannot move a lag living in `alpha`. What it *does* do is evaluate the flow
direction at a state the material never visits: the elastic predictor overshoots the real increment by
roughly the ratio of elastic to elastoplastic modulus, so B trades a one-step-old direction for a
never-visited one — hence B's terminal `q/p` of 1.11-1.26 against the implicit 1.55. A's whole-tensor
lag keeps magnitude and direction consistent with each other; B's split does not.

## 7. MUST NOT say

- **Not** "the direction correction is within noise" — 1.8-22x worse, 9-54 % off in `q/p`. Equally
  **not** "B never wins": on the two `p0 = 5` breakdown rows its max error is 50x and 230x below
  A's, in a regime where the implicit companion has itself collapsed.
- **Not** "A gets the volumetric sign right at phase transformation". It does not (§4); what survives
  is that the sign error is *inconsequential*, warranted by §2's 500x. And **not** "the crossing is
  where the corner error lives" — it is the quietest step, and this says nothing about the `p_min`
  floor (P0's G3 found the O(1) error there) nor about a BVP corner whose path is not monotone TX.
- **Not** "B is unsymmetric": asymmetry is `1e-4`-`1.6e-3`; the objection is that `Ce(p_n)` is not B's
  tangent. **Not** any cost or wall-clock claim — nothing here is timed (`~100x` stays #792 T8 / P3).
  **Not** anything about reversal: `alpha_in` reversal detection was not exercised, and every row is
  monotone drained compression at one Gauss point.

## 8. Recommendation for P2

**Adopt nothing; keep operator A as ADR §2 has it.** The deciding number is **`A/B(max) = 0.045`** at
`p0 = 100`, `d eps_z = 1e-5` — B is 22x worse on exactly the path the question was asked about —
backed by `max|J-Ce|/max|Ce| = 0.70`: B forfeits ADR §2's benefit #2 (a global step linear on a frozen
operator) to buy an error *increase*. C is the same verdict for a smaller reason — indistinguishable
from A at `p0 = 100` (1.40e-3 vs 1.39e-3), 11x worse at `p0 = 5` — so it costs a trial-state
`state_dependent()` per global iteration and buys nothing. **What would reopen it is BVP evidence:** a
Gauss point in the `hypo_bearing` corner whose *committed* `D` oscillates in sign step to step, not the
single monotone crossing a triaxial ramp gives — if P1's per-element `implexError` shows the worst
points straddling `D = 0` rather than on the `p_min` floor, this study does not cover that case.

## Log

- 2026-09-06 — LANE E. `python3.12 …/sanisand_implex_oracle.py --gate GE` (also in `--gate all`); it
  needs the gitignored `data/` probe CSVs, the unchanged P0 seeds. The oracle gained `Implex` forms
  `"T"`/`"H"` and `gate_GE`; A left byte-identical. §2-§4 are one `--gate GE` run; §5's `asym(J)` column
  came from a second `_tangent_probe` call, same seeds, its other columns reproducing exactly.
