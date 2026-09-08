---
title: "ADR 93 / ADR 92 — the seat replayed: step 331, element 4095 GP 8"
project: Ladruno
type: results
status: "REPRODUCTION GATE MET at 1e-12 (this lane's first): the oracle reproduces the binary's step 332 at the seat exactly, and the IMPL-EX committed error is EXACTLY `||Ce:(dep(n+1) - f dep(n))||/den`. Hypothesis (e) is half-confirmed (the state IS softening: `Kp = -0.757 G`, `t4 = -0.0996 G`) and half-refuted (the plastic increment vanishes superlinearly as ds -> 0, so the term is NOT O(sigma) at any ds); the 0.46 needs the seat's own strain increment, 164x the step's kinematic strain."
priority: high
owner: nmora
related: ["[[93_ladruno_sanisand_zero_confinement_adr]]", "[[_adr93_p0_ring_replay]]", "[[92_ladruno_sanisand_implex_adr]]", "[[_adr92_p0_oracle_results]]"]
tags: [adr, sanisand, implex, oracle, seat, measurement]
updated: 2026-09-07
---

**Build:** numpy under `python3.12`, no binary, no C++. Oracle
`adr92_p0_oracle/sanisand_implex_oracle.py` (variant A), driver `adr93_p0/seat_replay.py`
(`--only gate|diag|probe|repro|prevent`), input `adr93_p0/data/census/` from the TIMs census leg
`D-L-dl-vt-dense-q10-sp-implex-census` (engine `c162833ed`, `LoadControl -ds`); material from its
`run.json` (`e_init 0.6271`, `-Presidual 1.01`, `-Pmin 0.0101`, `IntScheme 1`, tol 0.1, cap 20000).

**Recovered / inferred / absent.** **READ, exact** — `sig`, `epsE`, `alpha`, `z`, `alpha_in`, `e`, `dGamma` at rows **331 and 332**
(`argmax_events.csv` carries all 26 `getState` slots); the sign pairing (`sig` negated, state slots
not) is CHECKED by `F(sig,alpha)` = 4.6e-11 / 3.6e-11 < TolF, any flip giving 1.2e+02, and `tr(eps)`
follows from `e = e_init - (1+e_init) tr(eps)`. **INVERTED, gated** — `d_eps` of step 332 (§1).
**INFERRED** — the deviatoric `eps`, from the wall hold census (s/B 0.017842, 76 rows later); it
feeds only `den`, whose `P_atm|eps|` leg is **1.69 %**. **ABSENT** — the state at step **330**,
hence `d_eps_331` and `d_eps_p(331)`: bounded and calibrated in §3. **This memo starts from the
row-331 committed state**, read not fitted; step 331's own 0.5309 is **not** reproduced.

## 1. Reproduction gate — MET

Solved for: `ModifiedEuler(state_331, d_eps) = sig_332`, 6 unknowns / 6 residuals. Held out:

| held out | oracle vs binary |
|---|---|
| `sig` residual (the inversion's target) / `dGamma` 6.7732e-07 | 1.6e-14 / 1.1e-9 rel |
| `tr(d_eps)` vs the void-ratio value — independent of the inversion | **2.2e-8** rel |
| `alpha` (‖·‖ 1.1325) / `z` (‖·‖ 4.2488) | **1.3e-12** / **2.3e-12** rel |
| `alpha_in` / substeps | **0.0** (the reset fired, `alpha_in_332 == alpha_331`) / **43 vs 43**, 0 force-accepts, 0 clamps |

The lane's **first met gate** (the ring replay reached a fitted 3e-5). **The increment is the
story:** `|d_eps| = 1.281e-4` on `ds = 3.906e-7 m`, so `|d_eps|/ds` = **327.8 /m** at the seat
against **8.44** at the ring on that step and **0.77** at the ring's steady rows 327-328 —
**164x** the step's own kinematic strain (`ds`/element = 7.8e-7).

## 2. Diagnosis at the row-331 state — no integration

`psi -0.18219` (dense), `e 0.626824`, `p 57.710 kPa` (with `p_r`), `eta = q/p 1.3003`,
`M^b 2.5181`, `M^d 0.4669`, **`eta/M^b 0.516`**, **`eta/M^d 2.785`**, `‖alpha_b-alpha‖ 1.3863`,
`‖alpha_d-alpha‖ 0.7975`, `(alpha_b-alpha):n +1.2213`, `(alpha_d-alpha):n -0.4536`,
**`D -0.02268`**, `A 0.0500`, `z:n -3.498`, `G 67480`, `K 157838 kPa`. With the **committed**
`alpha_in`: `aain = (alpha-alpha_in):n` **-0.16449**, `h = b0/aain` **-1086.7** (`b0` 178.76),
**`Kp = 2/3 p h (alpha_b-alpha):n` = -5.106e+04 kPa = -0.757 G**, denominator
**`t4 = Kp + 2G(B - C tr n^3) - K D (n:r)` = -6.718e+03 = -0.0996 G**; after the reset step 332
performed (`alpha_in := alpha_n`, bit-exact), `aain` **0.0** exactly, `h` the 1e10 guard,
`Kp = t4 = +4.699e+11`. The seat is **not at the bounding surface** (`eta/M^b 0.52`) but far past
the dilatancy one (`eta/M^d 2.79`, `D < 0`); its softening is `aain < 0` — a stale `alpha_in`, a
*pending reversal* — and `t4` is a near-cancellation **20x smaller than the elastic `2G`**.
`dGamma_331 = 0` exactly: ModifiedEuler's `dg < -SMALL -> dg := 0` branch fired. The return to a
small `d_eps` from row 331, along `d_eps_332`:

| `\|d_eps\|` | 1e-7 | 1e-6 | 1e-5 | 1e-4 | **1.281e-4** (actual) | 1e-3 |
|---|---|---|---|---|---|---|
| `\|d_eps_p\|/\|d_eps\|` | 0.0006 | 0.0087 | 0.0675 | 0.4616 | **0.4972** | 0.7067 |
| `\|Ce:d_eps_p\|/\|sig\|` (substeps) | 1e-7 (1) | 1e-5 (1) | 0.0009 (1) | 0.0540 (64) | **0.0745** (81) | 0.8271 (606) |

**Hypothesis (e): half-confirmed, half-refuted.** Confirmed — the seat *is* a dense point in a
softening state with a near-vanishing algorithmic modulus. Refuted — the plastic increment is
**not** of order `|sig|/G` for a tiny `d_eps`, it collapses superlinearly. It is O(sigma) **only
because the seat's strain increment is 1.3e-4, not 1e-6**: the loop supplies the size, the softening state the superlinearity.

## 3. The error, reproduced

`ManzariDafalias` holds `G, K` committed for the whole substep loop, so
`d_sigma = Ce:(d_eps - d_eps_p)` exactly and **`sigma~ - sigma_impl = Ce:(d_eps_p(n+1) -
f d_eps_p(n))`** — residual **1.8e-12 kPa**: the committed error measures the **step-to-step
change of the plastic increment**, nothing else. `den = 119.693 + 1.990 = 121.683 kPa`:

| quantity | value |
|---|---|
| binary numerator `err*den` | **56.278 kPa = 47 % of `\|sig\|`** |
| `\|Ce:d_eps\|` / `\|Ce:d_eps_p(332)\|` (reproduced) | 8.847 / **3.499** kPa (2.9 % of `\|sig\|`) |
| implied `\|Ce:d_eps_p(331)\|` | **[52.8, 59.8] kPa = 44-50 % of `\|sig\|` = 15-17x step 332's** |
| calibrated `\|d_eps_331\|` reaching it; reproduced `err` | 1.279e-3 (**10.0x** step 332's), 411 substeps (binary 240) -> `err` **0.4625** vs binary **0.4625**, both steps at the same `dt` |

## 4. What would have prevented it — one number each (P2 candidates, not decisions)

Same state, same `d_eps_332`, §3's calibrated `d_eps_p(n)`; only the operator varies, so (0) is exact.

| arm | `err` | x binary |
|---|---|---|
| (0) as shipped: `f = 1`, history direction (form A) | 0.4625 | 1.000 |
| (i) `alpha = 0.5` (extrapolation halved) | **0.2169** | 0.469 |
| (ii) lane-E variant B — trial-flow direction, form T | **0.4623** | 1.000 |
| (iii) `f` capped at 1 after a refusal chain | **0.4625** | 1.000 (`f = 1` exactly: INERT) |
| (iv) history + `dt` from step 329 (last full `ds` 2.5e-5 under tol) | **0.0287** | 0.062 |

(ii) **refutes the direction hypothesis at this seat** — magnitude carries the error. (iv) works
only because `f = ds/ds_full = 1.6e-2` makes the history term negligible: it lands on the `f = 0`
floor (0.0288), i.e. IMPL-EX switched off there; its strain-scales-with-`ds` variant is
degenerate, `|d_eps| = 8.2e-2` going elastic.

## MUST NOT say

- That (e) was confirmed: `Kp < 0` and `t4 = -0.0996 G` **are** measured, but "a tiny `d_eps`
gives a plastic increment of order `|sig|/G`" is **refuted** (ratio 6e-4 at `|d_eps|` 1e-7).
- That step **331**'s 0.5309 was reproduced (state 330 is absent; only **332**'s 0.4625 is), or
that `|d_eps_331| = 1.279e-3` is measured — it is the value making the *binary's own* error come
out at this state and the substep count wants ~5e-4. Quote the range, not the digit.
- That the seat is "at peak" in the bounding-surface sense (`eta/M^b = 0.52`) — the softening is a
stale `alpha_in`; that the loop (d) is out at the seat (it carries **164x** the step's kinematic
strain, the ring 4.2x); or that `den` is exact (its 1.69 % `P_atm|eps|` leg is inferred).
- Any absolute `Q`, curve or capacity from this leg — it is the census leg, whose 0.005 hold
already moved the curve 29 % (ADR 93 log, 2026-09-07).

## Log

- 2026-09-07 — Seat replayed from the row-331 committed state; gate MET at 1e-12, error identity
`Ce:(dep(n+1) - f dep(n))` established to 1.8e-12 kPa. (e) half-confirmed (`Kp = -0.757 G`,
`t4 = -0.0996 G`, `psi = -0.182`, `D = -0.0227`), half-refuted (the plastic increment vanishes
superlinearly). P2 candidates: `alpha 0.5` 0.47x, variant B 1.00x, `f`-cap inert, full-step
history 0.06x but only by switching the extrapolation off.
