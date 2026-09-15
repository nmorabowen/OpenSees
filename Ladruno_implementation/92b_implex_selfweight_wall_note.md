---
title: "ADR-92 F10 — the IMPL-EX self-weight wall: the refusal storm is the STEPPING CONTROLLER, not the soil"
project: Ladruno
type: measurement note
status: "MEASURED 2026-09-14 — verdict delivered, no C++ owed"
priority: high
owner: nmora
amends: 92_ladruno_sanisand_implex_adr
requested_by: "TIMs Workbench, F10 (self-weight strip footing act)"
engine: "ladruno tip 9c2f964eae3bbd1a055c3ede81381a6c601a982b (pinned release build, nothing built for this WP)"
related:
  - "[[92_ladruno_sanisand_implex_adr]]"
  - "[[LadrunoSANISAND_implex_guide]]"
  - "[[93_ladruno_sanisand_zero_confinement_adr]]"
  - "[[95_prandtl_reissner_campaign_report]]"
  - "[[_adr90_tau0_qu_band]]"
  - "[[LEDGER_quirks]]"
tags: [adr, material, sanisand, implex, soil, stepping, tims, measurement]
aliases: [ADR-92 F10, "implex self-weight wall"]
updated: 2026-09-14
---

# ADR-92 F10 — the IMPL-EX self-weight wall

> [!info] Numbering
> **`92b`, not a new ADR number.** This amends ADR 92; the folder's ADR high-water
> mark is **97** (`97_ladruno_asdp_closest_point_adr.md`) and is untouched. The
> `<n><letter>` suffix follows `40a`/`40b`, `75a`/`75b`/`75c`, `80a`–`80d`, `49a`.
> `92b` rather than `92a` because `_adr92_p1_*` / `_adr92_p2_9_*` already occupy
> the underscore-prefixed working-note space for P1/P2.

**A DIAGNOSIS. No C++, no banner, no flag defaults changed.** Test bed:
`Ladruno_files/testbed/hypo_bearing/adr92_f10/` (driver, reducer, artefacts,
README with the full deck table and the declared differences from the TIMs deck).

---

## 0. The report, and the one-line answer

TIMs measured, on a plane-strain strip footing on SELF-WEIGHT SANISAND under
`-implex -implexControl`, that the control's refusal counter climbs from the
first push step (204 120 at step 3 on 9 720 Gauss points), that every refusal
cuts the step, and that the harness step falls below its floor at `s/B = 0.0125`
— while the same material class and the same build carried the weightless ADR-95
slab campaign to `s/B = 0.15`. Minimum `p'` in the block is 3.8 kPa, so this is
not the ADR-93 apex wall.

**Reproduced** on the fork's own deck (leg B below: `FLOOR` at `s/B = 0.0085`,
724 refusals, every one of them in the `-implexControl` bucket).

**The answer is none of the three candidates the brief offered, and it is not
the soil.** `implexError` is *first order in the step* — measured, on a
controlled step-size refinement at the seizure settlement — and `-implexControl`
is an **absolute bound** on it, so every deck has a maximum admissible step at a
given tolerance. The ADR-63 D16 adaptive controller **can only discover that step
by overshooting it**: it doubles `ds` after six clean steps, the clock ratio
`f = dt_{n+1}/dt_n` becomes 2 on that step (doubling the extrapolated plastic
increment on top of the doubled strain increment), the error crosses `tol`, the
step is refused, halved, and six steps later doubled again. Pinning the growth
factor at 1.0 and changing **nothing else** — same tolerance, same guards, same
soil, same mesh — takes the refusal count from **724 to 6** and the reach from
**0.0085 to 0.0265** (leg K).

Self-weight enters only by *setting where that maximum step sits*: it is what
puts a low-`p'`, strongly dilatant ring immediately under the footing edge
(5.7–11 kPa) while the rest of the block sits at 30–107 kPa, and an absolute
tolerance is set by the worst point in the block.

---

## 1. The deck

`B = 1.5` m strip, `15B x 12B` box, plane strain, 285 `LadrunoBrick
-formulation bbar`, **2 280 Gauss points**; `gamma' = 9.81` kN/m3 as a body
force, 7.65 kPa surcharge outside the footprint, 18.4 kN/m of footing dead load;
`LadrunoSANISAND` on the ADR-92 CP1 / ADR-86 §5 Gorini set
(`tests/test_ladruno_sanisand.py::_PARAMS`, `e_init = 0.6944`), `IntScheme 1
TanType 2`, `-Presidual 0 -Pmin 0.0101 -maxSubsteps 1000`; rough rigid footing;
`confine -> updateMaterialStage 1 -> push`; `LoadControl(-ds)` on an `sp`
pattern under a `Linear` series with the `Transformation` handler (the guide's
push-idiom precondition — **never `DisplacementControl`**); R3's adaptive
controller with `DS_BASE = 2e-5`, `DS_MIN = 2e-7`, `DS_MAX = 1e-3`,
`GROW_AFTER = 6`, `SUBDIV_BUDGET = 80`.

Controls, per leg: gravity resultant identity, the **1-D geostatic patch**
(`3.1e-13` measured — the check the resultant is structurally blind to), and
`eta_max/M_c = 0.7275 < 1` at the stage flip.

**The K0 device.** `nu* = K0/(1+K0)` is passed as the material's Poisson ratio,
so the stage-0 elastic gravity state is exactly `K0` at every depth (exact
because `nu` is constant, even though SANISAND's moduli go as `sqrt(p)`).
**Declared difference:** the fork's `LadrunoSANISAND` takes `nu` positionally and
offers no way to put it back after the K0 stage, so `nu*` is held for the whole
leg rather than being temporary. The rest of the declared differences from the
TIMs deck (calibration, mesh density, footing kinematics) are in the test bed's
README; **measured `psi` on this deck is `-0.128 … -0.109`** against TIMs' quoted
`~ -0.13`.

**The at-rest census matches the report.** At `K0 = 0.455`, **2 280/2 280
= 100.00 %** of Gauss points sit on the DILATANT side of `M^d` (median
`eta = 0.857`, median `M^d = 0.666`) — TIMs' 99.96 %. Minimum `p'` is 6.4 kPa
(TIMs: 3.8).

---

## 2. Which criterion fires — the answer the responses give directly

`-implexControl` has exactly two refusal branches
(`SRC/material/nD/LadrunoSANISAND.cpp`):

| line | test | counted in |
|---|---|---|
| `:2907` | `implexPrimed && mImplexError > errorTol` | — the gate on both |
| `:2919` → `:3019` | `... && dtAbs >= reductionLimit * mImplexDt0` — **above** the reduction floor | `implexRefusals[2]` (`control`) via `noteRefusalControl()` |
| `:2919` else → `:3154` | **at** the floor, `-implexFloor implicit` (the default) delivers the companion and does **not** refuse | `implexGuards[0]` via `noteFloorFallback()` |
| `:2716` / `:3235` | the companion hit `-maxSubsteps` | `implexRefusals[3]` (`companion`) |

Measured across all fifteen legs:

* **every refusal is the 0.05 error tolerance above the floor.** Leg B: 724
  total, 724 `control`, **0 `companion`, 0 `signChange`**;
* **the reduction limit never binds at the shipped `0.01`.** `implexGuards[0] = 0`
  on every leg except **F2**, where `reductionLimit` was deliberately raised to
  `0.5`: it then fired **151** floor fallbacks and bought **+71 % reach**
  (0.0085 → 0.0145). This is the first fork measurement in which
  `reductionLimit` is not inert (the guide §7 records it inert on the CP1 deck)
  — it is inert **because it is set two decades below the working step**, not
  because the branch does not work;
* **the companion never failed** at `-maxSubsteps 1000` — 0 on B, C, D, E, K, and
  ≤ 42 anywhere (M).

The throttled warning line (`:3040`–`:3046`) names both quantities, and on this
deck reads e.g.

```
-implexControl REFUSES this step -- implexError 0.0609871 > tol 0.05 at
|dt| = 4e-05 (floor 2e-07, f = 0, |d_eps_p(n)| = 2.55013e-05)
```

`floor 2e-07` against `|dt| = 4e-05` — two decades of headroom. **The floor is
not in play.**

---

## 3. The leg table

`h0 = 0.5` m, 2 280 Gauss points, wall budget 400 s/leg, target `s/B = 0.05`.
Modes: `TARGET`/`BUDGET` are admissible, `FLOOR`/`WALL` are seizure/stop.

| leg | what it changes against B | s/B | mode | steps | subdiv | wall s | refusals tot / ctl / comp | floor fallbacks | f=0 guard |
|---|---|---|---|---|---|---|---|---|---|
| **A** | weightless + uniform 10 kPa (the ADR-95 campaign condition) | 0.0189 | BUDGET | 656 | 81 | 381 | 1545 / 1545 / 0 | 0 | 123 055 |
| **B** | **the reported configuration** (`tol 0.05`, `K0 = 0.455`) | **0.0085** | **FLOOR** | 364 | 53 | 142 | **724 / 724 / 0** | **0** | 31 491 |
| **C** | `K0 = 0.818` via `nu* = 0.45` (the slab's own `nu`) | 0.0209 | BUDGET | 511 | 81 | 188 | 1727 / 1727 / 0 | 0 | 13 484 |
| **D** | `-implexFactor controlIter` | 0.0089 | WALL | 102 | 13 | 408 | 725 / 725 / 0 | 0 | 5 289 |
| **E** | a hold at the flip + 2e-6 first step, growth 1.5 | 0.0139 | FLOOR | 431 | 40 | 122 | 506 / 506 / 0 | 0 | 17 697 |
| **F1** | `tol = 0.1` (the shipped default) | 0.0117 | FLOOR | 236 | 37 | 110 | 342 / 313 / 29 | 0 | 24 393 |
| **F2** | `reductionLimit = 0.5` (floor raised so it can bind) | 0.0145 | BUDGET | 615 | 81 | 201 | 1246 / 1246 / 0 | **151** | 45 016 |
| **F3** | `tol = 0.5` | **0.0500** | **TARGET** | 114 | **2** | 169 | **18 / 15 / 3** | 0 | 2 991 |
| **G** | **implicit** (no `-implex`) — the reference | 0.0023 | WALL | 29 | 0 | 413 | **0 / 0 / 0** | 0 | 0 |
| **H** | + 100 kPa uniform surcharge (min `p'` up two decades) | 0.0135 | BUDGET | 496 | 81 | 129 | 2754 / 2748 / 6 | 0 | 24 254 |
| **I** | `-implexGuard off` | 0.0157 | FLOOR | 522 | 78 | 166 | 1032 / 1020 / 12 | 0 | 0 |
| **J** | `-implexGuard off` **and** `-implexTrialGuard off` | 0.0185 | BUDGET | 499 | 81 | 183 | 12 758 / 12 756 / 2 | 0 | 0 |
| **K** | **growth factor pinned at 1.0** (`ds` never doubles) | **0.0265** | WALL | **1 984** | **0** | 400 | **6 / 6 / 0** | 0 | 57 109 |
| **L** | constant `ds0 = 8e-5` with growth 1.0 | 0.0185 | WALL | 2 765 | 3 | 400 | 144 / 144 / 0 | 0 | 51 484 |
| **M** | growth factor 1.25 instead of 2 | 0.0113 | FLOOR | 371 | 25 | 109 | 253 / 211 / 42 | 0 | 21 281 |

Reading, in order of size of effect **at the registered `tol = 0.05`**:

* **K (growth 1.0): 724 → 6 refusals, 0.0085 → 0.0265, zero subdivisions in
  1 984 converged steps.** Nothing about the soil, the guards or the tolerance
  changed. This is the finding.
* **M (growth 1.25): 0.0113** — a *gentler* growth is not a fix. Any growth
  re-enters the loop; it just takes longer to overshoot.
* **L (constant 8e-5, no growth): 144 refusals, all in the first `s/B ~ 0.001`,
  then 2 765 steps with 3 subdivisions.** The early push is where a
  too-large constant step is punished; after three halvings it runs clean —
  and with growth pinned at 1.0 it can never recover the step it lost, which is
  why K (which starts at the base step) out-reaches it.
* **F3 (`tol = 0.5`) reaches the target** `s/B = 0.05` with **18 refusals and 2
  subdivisions in 169 s** — the cheapest reach in the table, and the only leg
  that terminated on `TARGET`.
* **I / J (guards off): +85 % / +117 %.** The P2-2 `f = 0` guard is a *refusal
  source* on this deck (see §5).
* **D (`controlIter`) is not a fix here**: 0.0089 against 0.0085 for 2.9x the
  wall time — the same "does not generalise" caveat the guide §12 already
  carries.
* **G (implicit) is refusal-free and slow**: `s/B = 0.0023` in 413 s, 29 steps,
  25 failed attempts. Leg K delivers **11x more reach per wall second** at the
  registered tolerance.

---

## 4. The decisive measurement — a controlled step-size refinement

Walk to `s/B = 0.00851` on a **refusal-free constant-`ds` path** (638 steps of
2e-5 m, `q = 164.686` kPa, reproduced to the digit by all seven runs), then take
**one** step of size `ds` and census every Gauss point's `implexDetail`. Each
`ds` is a separate process on the same deterministic path, so all seven probe
steps start from the same committed state.

| `ds` (m) | median err | max err | GPs > 0.05 | max err / previous |
|---|---|---|---|---|
| 8e-5 | 2.14e-05 | **1.126e-02** | 0 | — |
| 4e-5 | 5.74e-06 | 5.778e-03 | 0 | 1.95 |
| 2e-5 | 1.35e-06 | 3.039e-03 | 0 | 1.90 |
| 1e-5 | 1.90e-06 | 1.672e-03 | 0 | 1.82 |
| 5e-6 | 2.54e-06 | 9.902e-04 | 0 | 1.69 |
| 2e-6 | 2.95e-06 | 5.851e-04 | 0 | 1.69 |
| 5e-7 | 1.66e-05 | 4.309e-04 | 0 | 1.36 |

Two things follow, and they are the spine of the verdict:

1. **The error is first order in the step** (ratios 1.95 / 1.90 / 1.82 for a
   factor-2 refinement) down to about `ds = 5e-6`, then flattens onto a
   `dt`-independent floor of `~4e-4` — **three orders below the 0.05
   tolerance**. That is textbook IMPL-EX (`O(dt)` plus a round-off/substep
   residual), and it means *there is no irreducible companion jump at this
   state*: a smaller step is always a way out, until the floor, and the floor is
   nowhere near the tolerance.
2. **At the very settlement where leg B seized, the maximum error over all 2 280
   Gauss points at a step FOUR TIMES the base is 0.011 — 4.4x under `tol`.** The
   state is not the problem. Leg B's seizure is therefore **path-induced**: it
   arrived at `s/B = 0.0085` having repeatedly overshot the tolerance, refused,
   halved and re-grown, and equilibrium on the extrapolated stress at each of
   those steps put it on a *different strain path* from the clean walk.

---

## 5. Where the refusing Gauss points are, at the step size that actually refuses

The per-Gauss-point census at a fixed `ds` (control tolerance set so high that
nothing can refuse; step 3, i.e. primed, so the un-primed exemption at `:2907`
does not apply):

| deck | `p'` min/med/max kPa | dilatant at rest | `ds = 2e-5` | `4e-5` | `8e-5` | `2e-4` | max err at `2e-4` |
|---|---|---|---|---|---|---|---|
| **B** self-weight `K0 = 0.455` | 6.4 / 29.2 / 106 | **100 %** | 0 | **4** (0.2 %) | 4 (0.2 %) | 81 (3.6 %) | 0.223 |
| **A** weightless, uniform 10 kPa | 6.3 / 6.5 / 7.0 | 100 % | 0 | 4 (0.2 %) | 8 (0.4 %) | 199 (8.7 %) | 0.157 |
| **C** self-weight `K0 = 0.818` | 8.9 / 40.4 / 146 | **0 %** | 0 | **0** | **0** | 20 (0.9 %) | 0.086 |
| **H** self-weight + 100 kPa | 70 / 93 / 170 | 100 % | 0 | **0** | **0** | **0** | 0.036 |

(counts are Gauss points with `implexError > 0.05`, out of 2 280.)

**"Wrong at every point at once" is not what happens at the operating step.** At
`|dt| = 4e-5` — the step size at which leg B's warnings actually fire — the
over-tolerance population is **four Gauss points out of 2 280**, all at
`|z| = 0.25` m directly under the footing edge (`dx = 0`), at `p' ~ 9.9` kPa. One
refusing Gauss point refuses the whole step, which is how four points buy 724
refusals. The block-wide picture (31 % of the >0.05 population deeper than 2B,
17 % more than 2B from the edge, `p'` up to 68 kPa) only appears at
`ds >= 2e-4` — a step **ten times** the deck's own base increment. A harness
whose first push increments are that large *will* see the report's "refusals from
the first step on thousands of points"; that is the same measurement, read at a
step size the tolerance cannot carry.

**Confinement is a modifier, not the mechanism.** `implexError` is normalised
(`ladrunoImplexMeasureError`, `:2513`-`:2530`; the denominator is `:2520`):

```
implexError = ||sigma~ - sigma_impl|| / ( ||sigma_impl|| + P_atm*||eps|| )
```

so a low-`p'` point's *relative* error is large for the same absolute
discrepancy. Leg H (min `p'` 70 kPa) has **zero** Gauss points over tolerance at
every step size up to 2e-4 — and still refuses **2 754** times in its own leg,
because its adaptive controller grows the step past 2e-4.

**The dilatant-at-rest state costs accuracy — measurably, and second-order.**
The contractant-at-rest twin (leg C, `K0 = 0.818`, 0 % dilatant) has a strictly
smaller error at every matched step and refuses nothing below `ds = 2e-4`, where
leg B already refuses 4 points at 4e-5. It also reaches 2.5x further. So
candidate (1) is *real as an aggravator*.

**The P2-2 guard is itself a refusal source here.** Of leg B's throttled refusal
lines, **30 of 49 report `f = 0`** — the extrapolation term was already switched
off by the P2-2 guard (`:2243`), so `sigma~` is a *pure elastic predictor* and
the "extrapolation error" the control refuses on is the plastic correction
itself, which no choice of `f` can improve. The guide's statement that the guard
trades "the prediction's accuracy, not the step" is true with `-implexControl`
**off**; with it **on**, lost prediction accuracy *is* a refused step. Turning
the guard off (leg I) buys +85 % reach. `implexGuards[1]` runs to 31 491 on leg
B — about 3.8 % of all Gauss-point-steps.

---

## 6. Verdict against the three candidates the brief named

| candidate | verdict | evidence |
|---|---|---|
| **(1)** the extrapolation error of the DILATANT-AT-REST state — "wrong at every point at once" | **Aggravator, not the mechanism.** The dilatant fraction is confirmed (100 % vs TIMs' 99.96 %) and it does cost accuracy — the contractant twin (C) refuses nothing where B refuses. But at the step size that actually refuses, the over-tolerance population is **4 of 2 280 points**, all under the footing edge. "Every point at once" is a `ds >= 2e-4` phenomenon, not a `4e-5` one. | §5 census table; legs B vs C |
| **(2)** the substepper's error control vs the control's tolerance at low `p'` | **REFUTED as the mechanism.** The companion never failed at `-maxSubsteps 1000` (`implexRefusals[3] = 0` on B/C/D/E/K); the reduction floor never binds at the shipped `0.01` (`implexGuards[0] = 0`); and the controlled refinement shows the error is first order in `ds` with a floor three orders under `tol` — no `dt`-independent jump. Low `p'` **is** a modifier through the error's own normalisation (leg H: 0 over-tolerance points at every `ds` tested), but raising `p'` does not remove the wall (H still refuses 2 754 times). | §2, §4, §5 |
| **(3)** the `nu*` device leaving an elastic operator / `alpha_in` / `alpha` history the first plastic extrapolation cannot use | **REFUTED.** At `K0 = 0.455` the device's `nu* = 0.31271` **is** the material's own calibrated `nu = 0.3129` to three decimals, so the K0 state it produces is the state stage-0 gravity reaches natively — leg B is simultaneously the "K0 reached by the material's own elasticity" control, and there is nothing anomalous left behind. And the device pushed *far* from the material's own value (leg C, `nu* = 0.45`) makes the run **better**, not worse. | §1, legs B and C |
| **(4) — the measured cause** | **The stepping controller against an absolute error bound.** `implexError` is `O(ds)`; `-implexControl` bounds it absolutely; the ADR-63 D16 controller discovers the bound only by overshooting it, and its `x2` recovery latch makes the clock ratio `f = 2` on exactly the step that overshoots. Pin the growth factor at 1.0 and the same deck, tolerance, guards and mesh go from 724 refusals / `FLOOR` at 0.0085 to **6 refusals / 0.0265 with zero subdivisions**. | leg K vs B; M and L as the dose-response |

**A caveat we cannot close.** TIMs report `K0 = 0.82` walling *earlier* (0.0082)
than `K0 = 0.455` (0.0125). This deck measures the opposite ordering by a factor
of 2.5 (C 0.0209 vs B 0.0085), and the at-rest census explains why on the fork's
deck (the `K0 = 0.818` state is 0 % dilatant, so it is further from `M^d` and its
extrapolation is better). Their `K0 = 0.82` leg must differ from ours in more
than `K0`; the comparison should not be quoted as agreement or disagreement until
that leg's deck is on the table.

---

## 7. Accuracy — what is NOT established

Per ADR 92 §8's reporting condition, an `-implex` curve must be confirmed against
the implicit twin **over the overlap only**. The implicit leg (G) reached
`s/B = 0.00227` in 413 s, so the overlap is short and lies on the steeply rising
part of the curve, where a small `s/B` mismatch between two adaptive step
sequences is worth several percent (the `_adr90_tau0_qu_band` caveat).

| arm | mean \|dev\| % vs implicit, `0.00045 <= s/B <= 0.00227` | max \|dev\| % |
|---|---|---|
| B (`tol 0.05`) | 6.00 | 8.93 |
| F1 (`tol 0.1`) | 5.03 | 6.06 |
| F3 (`tol 0.5`) | 5.38 | 6.67 |
| K (growth 1.0) | 5.31 | 6.30 |
| I (`-implexGuard off`) | 7.09 | 9.94 |
| J (both guards off) | 5.70 | 6.71 |
| D (`controlIter`) | 4.24 | 4.96 |

Every `-implex` arm agrees with every other to 1–3 % and sits 4–7 % above the
implicit twin over this window. **That window is too short and too steep to be an
accuracy verdict, and this note does not issue one** — in particular, `tol = 0.5`
reaching `s/B = 0.05` (leg F3) is a REACH result, not a capacity, and no peak or
plateau is claimed on any leg here. The honest statement is: the tolerance
changes reach by a factor of 6 and moves the curve by less than the
implicit-vs-IMPL-EX offset already present at `tol = 0.05`.

Also not established: mesh convergence of any of this (one mesh, `h0 = 0.5`);
behaviour past `s/B = 0.05`; parallel behaviour (ADR-92's `sendSelf` limit
stands); and leg E's load column, whose reference reaction was captured before
its hold and is therefore offset — **leg E's reach is usable, its `q` is not.**

---

## 8. Recommendation

**For a self-weight deck, in order.**

1. **Stop growing the step, or grow it gently and shrink it pre-emptively.**
   This is the single highest-value change and it is in the *harness*, not in
   `SRC/`. Pin the growth factor at 1.0 (measured: 724 → 6 refusals,
   +212 % reach, 0 subdivisions) or, better, drive the controller off
   `avgImplexError` — read it once per step and shrink when it approaches `tol`
   — so the bound is approached from below instead of discovered by refusal.
   A growth factor of 1.25 is **not** enough (leg M).
2. **Do not run a self-weight SANISAND push at `tol = 0.05` with a doubling
   controller.** If the controller cannot change, raise the tolerance: `0.1`
   (the shipped default) buys +37 %, `0.5` reaches the target with 18 refusals
   — but then §7's caveat binds and the curve is not a capacity until an
   implicit twin covers it.
3. **Raise `reductionLimit`.** It is inert at the shipped `0.01` on this deck
   because the floor sits two decades below the working step. At `0.5` it fires
   151 times, converts refusals into floor fallbacks, and buys +71 % with the
   tolerance untouched — the cheapest single-flag change that keeps `tol = 0.05`.
4. **`-implexGuard off` is a real lever (+85 %) and should not be reached for
   blind:** the guard exists for ADR-93's softening seat, and this deck does not
   test that. Record it; do not recommend it.
5. **`-implexFactor controlIter` is not the tool for this** (0.0089 vs 0.0085 at
   2.9x the wall time), exactly as the guide §12's "does not generalise" caveat
   predicts.
6. **The implicit lane is NOT the only lane for a self-weight deck.** Measured on
   the same box, same deck: implicit 413 s for `s/B = 0.0023` (29 steps, 25
   failed attempts, 0 refusals); IMPL-EX with a pinned growth factor 400 s for
   `s/B = 0.0265`. That is **11x more reach per wall second** at the registered
   tolerance, with 6 refusals. A cluster implicit run remains the right answer
   when the curve must be a capacity; for reach on a desktop it is the expensive
   option, not the safe one.

---

## 9. Two defects found in passing (recorded, not fixed)

* **`implexPrimed` is a bare sign test.** `LadrunoSANISAND.cpp:2905` reads
  `const bool implexPrimed = (this->GetNorm_Cov(mImplexDEpsP) > 0.0);`. A leg-B
  refusal was measured at `|d_eps_p(n)| = 6.66e-12` with `f = 1` and
  `implexError = 0.211`: a Gauss point whose committed plastic history is
  numerically zero counts as *primed*, loses the un-primed exemption — which
  exists precisely because "a pure elastic predictor's error must scale with
  `d_eps`; this one does not" (`:2896`) — and is refused on an error
  that is the elastic predictor's own drift. This is the same defect shape
  P2-5b already fixed for the reversal reset (an absolute threshold that cannot
  separate a real increment from noise, replaced by a *relative* one), left
  unfixed here. A relative test (`||d_eps_p(n)|| > rel * ||d_eps||`) would be the
  matching fix. **Not fixed in this WP** — F10 is a diagnosis, and the fix needs
  its own gate and mutation score.
* **The P2-2 guard and `-implexControl` work against each other.** The guard
  trades prediction accuracy for safety by forcing `f = 0`; the control then
  refuses the step *because* the prediction is inaccurate. 30 of leg B's 49
  throttled refusal lines are at `f = 0`. Either the guard should exempt the
  point it just guarded from the control's tolerance for that step (the same
  shape as the un-primed exemption), or the control should measure a guarded
  point against the elastic predictor's own expected error. **Design question
  for ADR 92, not a fix here.**

Both are in `LEDGER_quirks.md`.
