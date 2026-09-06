---
title: "ADR 92 / CP1 — the surcharge legs (owner decision B): results"
project: Ladruno
type: results
status: "MEASURED — no plateau, no peak; the binding resource is now WALL CLOCK; decision B's own effect on the clamp is UNMEASURABLE with the shipped diagnostic"
priority: high
owner: nmora
related:
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[_adr92_p0_oracle_results]]"
  - "[[_adr90_tau0_qu_band]]"
  - "[[86_ladruno_sanisand_handoff]]"
tags: [adr, sanisand, surcharge, gate-u, cp1, measurement]
updated: 2026-09-05
---

# ADR 92 / CP1 — decision B measured

> [!warning] **VERDICT — still no capacity, and the reason has moved again**
> Four legs (2 densities × `Q` = 2 and 10 kPa, cap 1000, 40-min budget each) on the
> merged engine. **Not one plateaus or peaks**, exactly as at GATE U and at T8. But the
> three things that changed are worth more than the repeat verdict:
>
> 1. **The binding resource is now the WALL CLOCK.** Every leg ends `WALL` with **6–23
>    of 80** subdivisions spent and the step still 100–3125× above its floor. The
>    integrator stopped being the blocker at T8; the controller never was; now it is
>    simply time. The deepest leg needed **2411 s to reach `s/B = 0.0616` of a `0.25`
>    target**.
> 2. **Decision B's own effect on the clamp cannot be measured with the shipped
>    diagnostic.** The low-`p` warning is throttled to a **process-wide budget of 10**
>    (`ManzariDafalias.cpp:1453-1454`), and every dense log — T8's and these — holds
>    exactly 10 events plus the suppression notice. `CLAMPING 10` means **"≥ 10, count
>    unknown"**. The handoff §5b's own prescription (*"capture the log and count
>    CLAMPING"*) is defeated by that throttle. **The quantity the owner's decision was
>    chosen to fix is the one quantity this campaign cannot report.**
> 3. **The mechanism has not formed, and the surcharge itself proves it.** A surcharge
>    adds `Q·N_q` *at collapse*. Between `Q` = 2 and 10 that is a predicted **+209 kPa**
>    (`φ = 33.0°`, `N_q = 26.1`). Measured at the deepest common checkpoint: **+14.0 kPa
>    on Gorini's density (index 0.067) and −12.1 kPa on the dense one (index −0.058)**.
>    At most **7 %** of the surcharge's collapse contribution is realised — an
>    independent confirmation of "no capacity" that uses the surcharge as a *probe*
>    rather than as a cure.

**Provenance.** Engine `c4d08f1874060812ef25cfddf49724e7175f2978` (the merged 86b branch),
`ops.ladrunoBuild()`-asserted by the driver. T8's rows were measured on `5c01bb5f7`, which
differs by the sentinel narrowing and tests — §7b.3 argues the two are indistinguishable on
a SANISAND-only deck, and that argument is a caveat on every delta below, not a defect.
Driver `sanisand_tau0_band.py --surcharge`, rollup `adr92_cp1_rollup.py`, outputs in
`Ladruno_files/testbed/hypo_bearing/adr92_cp1/`. Four concurrent processes,
`OMP_NUM_THREADS=4`, so wall columns are upper bounds under contention.

---

## 1. What decision B is, and where it goes in the deck

`--surcharge Q` applies `Q × tributary area` downward on every free-surface node **not**
under the footing, as a third `Plain` pattern — the idiom of
`tests/test_r3_prandtl_collapse_gate.py:405-408`, **whose own `Q0` is 10 kPa**. So the
`Q = 10` leg puts this deck in the confinement configuration the testbed's one reliably
plateauing gate already uses.

Three placement choices, each load-bearing:

- **After** the geostatic controls — the 1-D patch check asserts `σ_zz = γz` and
  `σ_xx = K0γz`, which a surcharge invalidates. (Patch error stays `1.2–1.8e-13`.)
- **Before** the stage flip and still at stage 0, so the surcharge is carried *elastically*
  into the state SANISAND starts plastic life from: a confinement condition, not a load
  history.
- `loadConst("-time", 0.0)` first, or the gravity pattern's `Linear` series keeps scaling
  `γ` past 1.0 while the surcharge ramps.

Footing nodes are **excluded**, so there is no `r_corr` to add back the way the R3 gate
needs; the measured `q` is the footing's own. Base-reaction identity re-asserted against
`γV + Q·A_surch`: **4.44e-16** on every leg.

## 2. Termination and depth

| leg | `Q` | mode | `s/B` end | plateau | peaked | tail % | subdiv | clamp | wall s |
|---|---|---|---|---|---|---|---|---|---|
| `e0.6944` | 0 *(T8)* | WALL | 0.05875 | NO | NO | 107 | 18/80 | 0 | 2705 |
| `e0.6944` | **2** | WALL | 0.04968 | NO | NO | **49** | 8/80 | 0 | 2400 |
| `e0.6944` | **10** | WALL | **0.06155** | NO | NO | **45** | 6/80 | 0 | 2411 |
| `e0.60` | 0 *(T8)* | BUDGET | 0.03201 | NO | NO | 276 | 81/80 | ≥10 | 1774 |
| `e0.60` | 0 *(T8)* | WALL | 0.03374 | NO | NO | 111 | 16/80 | ≥10 | 2718 |
| `e0.60` | **2** | WALL | **0.04326** | NO | NO | 162 | 23/80 | ≥10 | 2401 |
| `e0.60` | **10** | WALL | **0.04327** | NO | NO | 110 | 19/80 | ≥10 | 2409 |

**Depth improves** — `0.0588 → 0.0616` on Gorini's density and `0.0337 → 0.0433` (+28 %)
on the dense one. **The Gorini tail flattens sharply**: 107 % of its initial slope at T8,
**45–49 %** with a surcharge. The dense tail moves 276 → 110–162 %, still steepening.
Against the WP1 criterion (**tail ≤ 2 %** over `Δ(s/B) ≥ 0.01`) every leg is an order of
magnitude away. **No leg is admissible as a capacity**; the `q_u` column of the rollup is
where each run stopped, not where the soil failed.

## 3. Matched settlement — the surcharge's price, pre-peak

| leg | `Q` | 0.002 | 0.005 | 0.010 | 0.020 | 0.040 |
|---|---|---|---|---|---|---|
| `e0.6944` | 0 | 55.76 | 120.94 | 233.17 | 454.60 | 875.27 |
| `e0.6944` | 2 | 60.68 | 130.15 | 249.74 | 479.50 | 920.21 |
| `e0.6944` | 10 | 62.75 | 133.28 | 253.52 | 484.07 | 934.17 |
| `e0.60` | 0 | 84.12 | 194.39 | 402.88 | 882.84 | — |
| `e0.60` | 2 | 89.70 | 205.09 | 420.60 | 888.11 | 1906.74 |
| `e0.60` | 10 | 91.85 | 207.37 | 420.25 | 875.45 | 1894.69 |

A surcharge buys **+5 to +9 %** on Gorini's density and **+4 to +7 %** on the dense one at
matched settlement — nothing like `Q·N_q`, because `N_q` is a collapse quantity and these
rows are pre-peak. **These legs are not comparable to §6.2's values without that
understanding, and the capacity column is not comparable at all until a leg is admissible.**

## 4. The mechanism-formation index — the finding worth keeping

Between `Q` = 2 and `Q` = 10 the *collapse* loads must differ by `8 × N_q = 209 kPa`. What
is actually measured:

| `s/B` | Gorini Δq | index | dense Δq | index |
|---|---|---|---|---|
| 0.002 | +2.08 | 0.010 | +2.15 | 0.010 |
| 0.005 | +3.12 | 0.015 | +2.28 | 0.011 |
| 0.010 | +3.78 | 0.018 | −0.35 | −0.002 |
| 0.020 | +4.57 | 0.022 | −12.66 | −0.061 |
| 0.040 | **+13.96** | **0.067** | −12.06 | −0.058 |

Gorini's index rises monotonically and reaches **6.7 %**; the dense one shows no systematic
effect at all and goes negative. **Read with its own error bar:** as a fraction of `q` these
differences are `+0.95 … +3.4 %` (Gorini) and `−1.4 … +2.4 %` (dense), against §5.3's
**0.8–1.4 %** solver-configuration floor — so the dense column is *inside* the scatter and
is not resolved, and Gorini's is only just outside it. The honest statement is therefore
one-sided: **the index is nowhere near 1, on either density, at every settlement reached.**
The bearing mechanism has not formed, and no amount of reading the load curve will make a
capacity out of these legs.

## 5. What this licenses, and what it does not

**MAY say:**
- Decision B is implemented, its equilibrium identity is exact, and it costs nothing in
  the controls (patch `1.8e-13`, `η/M_c` 0.64, `OutsideBounding` 0).
- It buys **depth** (+28 % on the dense density) and, on Gorini's, a markedly **flatter
  tail** (107 → 45 %).
- It does **not** produce a plateau or a peak within a 40-minute leg.
- The binding resource is now **wall clock**, not the constitutive integrator (T8) and not
  the stepping controller (6–23 of 80).
- The failure mechanism is **≤ 7 % formed** by the surcharge's own probe at `s/B ≈ 0.05`.

**MUST NOT say:**
- That the surcharge did or did not reduce clamping. **The counter saturates at 10.**
- Any capacity number. No leg is admissible; `q_u raw` is a stopping point.
- That `Q = 10` is better than `Q = 2`. It is deeper on Gorini and identical on the dense
  one, and the differences sit at the scatter floor.
- That the dense density behaves like Gorini's. Its tail is still steepening.

## 6. Owed next, in order

1. **Give the legs the wall clock.** Every leg died on time, none on a solver failure, with
   the step still 100–3125× above its floor. `s/B = 0.0616` in 2411 s against a `0.25`
   target means a plateau — if there is one — is hours away, not minutes. **This is the
   cheapest decisive experiment and it is launched** (§7).
2. **A cheaper deck**, T8 §7b.4's undone half. The domain is R3's `60 × 20 m` for a 2 m
   footing — sized for a *weightless* half-space. A weighted mechanism is far more local,
   and a purpose-sized domain is the difference between hours and minutes per leg. It
   changes the BVP, so it needs its own extent study; it is the real fix.
3. **A clamp COUNTER, not a throttled message.** `n_clamping` saturates at 10 and cannot
   answer the question decision B was taken to settle. A `getResponse`-style integer per
   material (or a per-process counter that keeps counting after it stops printing) is a
   small, contained change. Owed a `LEDGER_quirks` row either way — the throttle silently
   caps a *measurement*, not just a diagnostic.

## 7. The deep run

Launched at the close of this session on the same engine and legs, `--wall 9000`
(2.5 h each, four concurrent), outputs in `adr92_cp1/deep_*`. The driver writes its leg
record at **every checkpoint**, atomically, marked `partial: true`, so an external kill —
which took both GATE U launches at ~1 h — loses nothing but the tail. Results append here
as §8.

## Log

- 2026-09-05 — Owner took decision **B (small surcharge)** at CP1 and asked for the run.
  Driver gained `--surcharge` (`af2d1324d`); four legs measured; no plateau; the clamp
  counter found saturated; the mechanism-formation index introduced and reported with its
  scatter floor. Deep run launched.
