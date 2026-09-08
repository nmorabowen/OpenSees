---
title: "ADR 93 / P0 — Gauss-point replay of the Esmeralda ring: results"
project: Ladruno
type: results
status: "P0 INCONCLUSIVE for the ADR's question — the instrumented point is CONFINED throughout (p 6 → 51 kPa, zero clamps, D_factor never fires), so it cannot price I.1, II.1 or II.2 and the §4 decision rule cannot be applied. Reproduction gate NOT MET (3e-5, not 1e-6: the dump carries no state slots). The tables stand as a CONTROL ARM — what these candidates cost at a confined ring point."
priority: high
owner: nmora
related: ["[[93_ladruno_sanisand_zero_confinement_adr]]", "[[_adr92_p0_oracle_results]]", "[[92_ladruno_sanisand_implex_adr]]", "[[86_ladruno_sanisand_pr3_tripwire_memo]]"]
tags: [adr, sanisand, oracle, low-confinement, measurement]
updated: 2026-09-06
---

# ADR 93 / P0 — what the ring dump can and cannot decide

> [!warning] **VERDICT — INCONCLUSIVE for the ADR's question**
> **The instrumented point is confined throughout.** It was chosen by GEOMETRY — first element
> outboard of the footprint, top row — and it loads: `p` rises monotonically from `≈ 6 kPa` to 51
> (dense) / 35 (gorini), `‖σ‖` from 11.6 to 114 / 75, `clamp_fired` is **0 on every one of the 362
> dumped steps**, and the `D_factor` sigmoid never fires. A 1.00× substep cut for I.1 on a point
> that never seizes says nothing about I.1; the II.1 / II.2 numbers are properties of a *confined*
> point. **The floor-bound ring, if it exists on this deck at all, is elsewhere or later.**
>
> **What P0 needs instead:** the Gauss point with **minimum committed `p` at the wall — chosen by
> STATE, not geometry** — from a leg that actually walls (the DisplacementControl legs at
> `s/B 0.001–0.002`, or these LoadControl legs once they stop), **with the internal state dumped
> every step: `alpha`(6), `z`(6), `alpha_in`(6), void ratio `e`** — so the reproduction gate can be
> met without fitting.
>
> **What this run does establish: a source correction to ADR-93 §1.** `m_Presidual` does **not**
> enter the elastic moduli — all three `GetElasticModuli` overloads (`ManzariDafalias.cpp:4834`,
> `:4876`, `:4896`) take `pn = tr(σ)/3` and floor it at **`m_Pmin`**. The §1 row reading
> `G(p + p_r)` is wrong: `p_r` is **purely plastic-side**, so today there is **no stiffness floor at
> all** — `G, K → 0` with `p` whatever `p_r` is. **`p_r,e` is a NEW elastic-only parameter, not a
> split of an existing coupling.**

**Build:** no binary, no C++; numpy under `python3.12`. Oracle `adr92_p0_oracle/…_oracle.py`,
driver `adr93_p0/replay_ring.py`, input TIMs `D-L-dl-vt-{dense,gorini}-q10-sp-{implex,implicit}/…/
ring_point.csv` (engine `c162833ed`, element 4047 GP 5), copied to `adr93_p0/data/`. The IMPL-EX
legs were still running: last steps used **154** (dense, `s/B 7.933e-3`) and **170** (gorini,
`8.200e-3`); the implicit twins are complete at **19** (`1.667e-3`).

## 0. Parameter match

No `-sp-` leg had written `run.json`; constants come from the sibling
`D-L-dl-vt-dense-q10-grow-implex` (`declared.material_common` + `declared.column.flags`). The
`dense`↔0.6271 / `gorini`↔0.6944 assignment is **measured**: the dilatancy the first committed step
requires matches each leg's own `e_init` to 8–14 %, the cross assignment not at all (dense at
0.6944 needs `D = −0.0068` against the model's `+0.0327`).

| quantity | oracle built-in | as used | action |
|---|---|---|---|
| `G0 nu Mc c lambda_c e0 ksi P_atm m h0 ch nb A0 nd z_max cz` | 264.32 0.3129 1.3309 0.71 0.027 0.83 0.45 101.0 0.005 1.3 0.968 3.5 0.05 5.75 12.5 1100.0 | identical | none |
| `e_init` | **0.6944** | **0.6271** (dense) / **0.6944** (gorini) | set per leg |
| `-Presidual` | **0.0** | **1.01** = `1e-2 P_atm` (vanilla) | set to 1.01 |
| `-Pmin` / `IntScheme` / `TanType` / `Den` | 0.0101 / 1 / — / — | 0.0101 / 1 = `ModifiedEuler` / 2 / 2.0 | none; last two inert here |
| `-maxSubsteps` | none | 20000 | inert (max seen 92 substeps/step) |

## 1. Reproduction gate — NOT MET (best 3e-5 against a 1e-6 target)

The dump's header states what it carries — *"available at S4: detail, refusals, stress, strain"* —
so `alpha`, `z`, `alpha_in` at the first step are **unknown (18 dof)**. Reconstruction: `e` by the
material's own rule on the dumped `eps`; `epsE` arbitrary (it never feeds back into the stress
path); `alpha` pinned to the yield surface `alpha = r − sqrt(2/3) m n`, leaving only the unit
normal `n`; `alpha_in`, `z` free. Stage 1 fits `n` (5 dof) plus the scalar
`aain = (alpha − alpha_in):n` that sets `h = b0/aain` — 6 dof against 6 residuals; stage 2 opens 18.

| leg / kind | seed | steps | step 1 | median | max | last |
|---|---|---|---|---|---|---|
| dense/implicit | naive | 18 | 6.01e-02 | 1.26e-01 | 1.58e-01 | 1.58e-01 |
| dense/implicit | **fitted** | 18 | **1.75e-04** | 6.14e-03 | 1.65e-02 | 1.65e-02 |
| gorini/implicit | **fitted** | 18 | **2.97e-05** | 7.76e-04 | 1.28e-03 | 1.28e-03 |
| dense/implex | **fitted** | 153 | 1.35e-02 | 2.03e-01 | 2.68e-01 | 2.67e-01 |
| gorini/implex | **fitted** | 169 | 1.57e-02 | 6.33e-02 | 8.69e-02 | 8.69e-02 |

Stage 1 alone reaches **2.78e-5** (dense/implicit, `aain = 0.431`) and **2.53e-5**
(gorini/implicit, `aain = 0.345`) on the first step. **(a) Model, units, component order, sign
convention and strain path are all correct**: with `n` from the step's own plastic strain
increment, the model's `D` matches what the dumped step requires to 8–14 % at every step tried —
the expected error of a frozen-coefficient estimate on a step that moves `p` by 8 %. **(b) The gate
cannot be met from this input**: 1e-6 presumes an exact seed, and an exact seed would not survive
either — ADR-92 P0's Control A measured the oracle moving *from itself* by 2.6e-5 under a `1e-15`
seed perturbation over 400 steps, because ModifiedEuler's accept/reject is discontinuous in the
state. **1e-6 is below this algorithm's own conditioning on a 150-step path.** The 3e-5 → 1e-2
drift is that same amplification, seeded by residual state error. So **everything below is
DIFFERENTIAL** — every arm runs from the same seed on the same strain path, so seed error is
common-mode and no absolute stress from it may be quoted.

> **OPEN ITEM for ADR 92 (not a finding).** On fitted seeds the IMPL-EX legs reproduce ~100× worse
> than their implicit twins — best 2.4e-3 (dense) / 6.8e-3 (gorini) against 2.8e-5 / 2.5e-5 — and
> opening all 18 dof did not close it, although the committed stress under `-implex` is meant to be
> exactly `ModifiedEuler(state_n, Δε)` (`LadrunoSANISAND.cpp:1949`). Best hypothesis: the strain
> increment handed to the commit-time companion is found on the *frozen elastic* operator, so it is
> larger and differently directed than Newton's and the unknown seed bites harder; a real
> difference in the committed state is not excluded. Decides on a state dump, not on a fit.

### 1b. IMPL-EX slot 0 split into numerator and denominator — exact, seed-free

`implex_err = ‖σ~ − σ_impl‖ / den`, `den = ‖σ_impl‖ + P_atm‖ε‖` in the doubled-shear contravariant
norm (`LadrunoSANISAND.cpp:1533-1555`); `σ_impl` and `ε` are both dumped, so the numerator comes
out offline — no oracle, no seed.

| leg | step | d(s/B) | slot 0 | den [kPa] | ‖σ~ − σ_impl‖ [kPa] | num / d(s/B) |
|---|---|---|---|---|---|---|
| dense | 1 | 3.33e-05 | 6.33e-02 | 11.57 | 7.32e-01 | 2.2e+04 |
| dense | 20 | 3.33e-05 | 6.74e-04 | 21.72 | 1.46e-02 | 4.4e+02 |
| dense | 80 | 6.67e-05 | 5.57e-04 | 52.36 | 2.92e-02 | 4.4e+02 |
| dense | 154 | 6.67e-05 | 1.31e-04 | 115.7 | 1.51e-02 | 2.3e+02 |
| gorini | 20 | 3.33e-05 | 5.71e-04 | 19.40 | 1.11e-02 | 3.3e+02 |
| gorini | 170 | 6.67e-05 | 9.41e-05 | 75.05 | 7.06e-03 | 1.1e+02 |

From step ~20 on, `num / d(s/B)` is flat to a factor 2–4 across a 4× change in `d(s/B)` and a 5×
growth in `den`: **the extrapolation error is O(ds) in the NUMERATOR, not only in the ratio** — the
ADR's 2026-09-07 03:00 entry confirmed independently of its denominator. The first ~10 steps are a
startup transient (`num/ds` 100× larger, `den` smallest, the growth ladder's first attempts), not
an O(1) event.

## 2. P0 item 2 — I.1, the substep norm's stress reference *(control arm: a confined point)*

`ManzariDafalias.cpp:1746-1753` reads `err = e` if `‖σ‖ < 0.5` else `e/(2‖σ‖)`: the two branches
are the one expression **`e / (2 max(‖σ‖, σ_ref))`** at `σ_ref = 0.5`, since below the switch the
divisor is exactly `1.0`. The oracle now carries `σ_ref` and is **bit-identical to the hardcoded
branch at 0.5** on both ring paths and on two synthetic paths that exercise the absolute branch
(`min ‖σ‖ = 0.189`). *(ADR-92's G0 gate could not be re-run — its `probe_binary_triaxial.py` CSVs
are not in the repo and this P0 runs no binary — so that bit-identity is its substitute.)*

| leg | σ_ref [kPa] | /P_atm | substeps | mean/step | max/step | substep cut | median Δσ | max Δσ |
|---|---|---|---|---|---|---|---|---|
| dense | 0.5 (today) | 4.95e-3 | 9450 | 61.8 | 92 | — | — | — |
| dense | 5 | 4.95e-2 | 9450 | 61.8 | 92 | **1.00x** | 0.0 | 0.0 |
| dense | 50 | 4.95e-1 | 9433 | 61.7 | 92 | **1.00x** | 4.2e-04 | 5.6e-04 |
| gorini | 0.5 (today) | 4.95e-3 | 10722 | 63.4 | 89 | — | — | — |
| gorini | 5 | 4.95e-2 | 10722 | 63.4 | 89 | **1.00x** | 0.0 | 0.0 |
| gorini | 50 | 4.95e-1 | 10592 | 62.7 | 88 | **1.01x** | 2.4e-02 | 7.6e-02 |

`σ_ref = 5 kPa` changes **nothing at all** — not one substep, not one bit of committed stress —
since `max(‖σ‖, 5) = ‖σ‖` throughout; `σ_ref = 50 kPa` (**50× the ADR's own defensibility bound**)
buys 0–1 % and costs up to 7.6e-2 in committed stress. **On a confined point I.1's saving is zero
by construction — which is what a control arm is for.** Nothing here prices I.1.

## 3. P0 item 3 — II.1 decoupled floors and II.2 / D5a *(control arm: a confined point)*

`p_r,e` is *added* inside `GetElasticModuli`, where `m_Presidual` is absent today; `p_r,p`
replaces `m_Presidual` everywhere else. The coupled control `p_r = 1.01` is today's binary.

| leg | arm | substeps | mean/step | median Δσ | max Δσ | min λ(Ce) [kPa] | replayed p |
|---|---|---|---|---|---|---|---|
| dense | `p_r = 1.01` coupled (control) | 9450 | 61.8 | — | — | 44125 | 6.53 – 64.8 |
| dense | `p_r,p = 0`, `p_r,e = 0` | 9812 | 64.1 | 1.3e-02 | 5.2e-02 | 44125 | 6.54 – 65.9 |
| dense | `p_r,p = 0`, `p_r,e = 0.1` | 9831 | 64.3 | 1.6e-02 | 5.2e-02 | 44488 | 6.54 – 66.1 |
| dense | `p_r,p = 0`, `p_r,e = 1.0` | 9875 | 64.5 | 4.2e-02 | 5.2e-02 | 47626 | 6.57 – 67.8 |
| gorini | `p_r = 1.01` coupled (control) | 10722 | 63.4 | — | — | 39695 | 6.43 – 45.5 |
| gorini | `p_r,p = 0`, `p_r,e = 0` | 11163 | 66.1 | 1.8e-02 | 5.7e-02 | 39695 | 6.43 – 46.6 |
| gorini | `p_r,p = 0`, `p_r,e = 0.1` | 11173 | 66.1 | 2.1e-02 | 5.7e-02 | 40026 | 6.44 – 46.8 |
| gorini | `p_r,p = 0`, `p_r,e = 1.0` | 11223 | 66.4 | 5.3e-02 | 5.7e-02 | 42887 | 6.47 – 48.3 |

`p` is 6.4–68 kPa, so `p_r,e ∈ {0, 0.1, 1.0}` moves `G ∝ sqrt(p + p_r,e)` by at most +8 % and
`λ_min(Ce) = min(3K, 2G)` never falls below **3.97e4 kPa**: **II.1's own question — does the floor
keep the tangent regular — cannot be asked here.** Dropping `p_r,p` from 1.01 to 0 is the larger
move (17 % on `p` in `psi`, `M^b`, `M^d`, `D` at the first step, 2 % at the last) and it costs **4 %
MORE** substeps, not fewer. **II.2 / D5a is likewise unmeasurable here**: `steps below
0.05 P_atm` is **0 on all eight arms**, so `min/max D` are identical to every printed digit with
the sigmoid on and off (dense control `−0.0569 / +0.0150`, gorini `−0.0371 / +0.0408`). Dropping
`p_r,p` to 0 does move gorini's `max D` by 44 % (0.0408 → 0.0230) — the `psi(p)` shift, **not** the
sigmoid.

## 4. The §4 decision rule — why it cannot be applied here

> *"if I.1 at a defensible `σ_ref` (≤ 1e-2 P_atm) cuts **the ring's** substeps by ≥ 10× … P1 = I.1
> … If it does not … P1 = II.1 + II.2 together, with I.1's declared `σ_ref` kept as the cost
> bound."*

The rule's subject is a Gauss point **at the floor**. Mechanically, `1e-2 P_atm = 1.01 kPa` gives a
cut of **exactly 1.00×** here — but only because `‖σ‖ ≥ 11.5 kPa`, i.e. because the antecedent's
*subject* is absent, not because its *test* failed. **Reading that as "P1 = II.1 + II.2" would be
inferring a verdict about the floor from a point that never reaches it, so this memo does not draw
it.** The rule stays unapplied until P0 is re-run on a floor-bound point.

**What to ask TIMs for, so it can be applied.** `ring_point.csv` picks by geometry (top row, first
element outboard of the footprint) — the *band peak*, which loads and confines. Instead: **(a)** a
per-step census of committed `p` and `clamp_fired` over the whole top element row, so the point can
be **chosen by state** — minimum committed `p` at the wall; **(b)** that point taken from a leg
that **walls** (the DisplacementControl legs at `s/B 0.001–0.002`, or these LoadControl legs once
they stop); **(c)** the internal state at every step — `alpha`(6), `z`(6), `alpha_in`(6), `e`, the
26 `getState` slots — without which no replay is gated better than ~1e-4, fitted.

## MUST NOT say

- That the §4 rule was applied, or that **P1 = II.1 + II.2** follows from this run. The rule needs
a floor-bound point; this one is confined. The run is **INCONCLUSIVE**.
- That I.1 was refuted or priced (**inert** here), that II.1 does or does not keep the tangent
regular (`λ_min(Ce) ≥ 3.97e4 kPa`), or that D5a was settled (the sigmoid never fired).
- That the oracle reproduces this BVP path. **NOT MET**: 3e-5 on one step, 1e-2 over 18.
- Any **absolute** stress, `q`, substep count or error from §2/§3 as a property of the binary —
only the arm-to-arm *differences* are seed-independent.
- That the ring "carries no load" or "sits at the floor" on these legs: `p` **rises** 6 → 51 kPa
and the clamp never fires.
- That the IMPL-EX / implicit reproducibility gap is understood, or is a finding of this memo — it
is an **open item for ADR 92**.

## Reproducing

```bash
cd Ladruno_implementation/adr93_p0 && python3.12 replay_ring.py --gate path
#  … --gate identity | calibrate (~5 min) | repro | implex | I1 | II1 (~8 min)
```

## Log

- 2026-09-06 — P0 run on the four `-sp-` legs as dumped (IMPL-EX legs still in flight).
**INCONCLUSIVE**: the instrumented point is confined, so I.1 / II.1 / II.2 could not be priced and
the §4 rule stays unapplied; the tables are kept as the confined-point control arm. Reproduction
gate not met, cause identified. ADR §1's `G(p + p_r)` row corrected against the source (`p_r` is
plastic-side only; there is no stiffness floor today). The IMPL-EX-vs-implicit reproducibility gap
is logged as an open item for ADR 92.
