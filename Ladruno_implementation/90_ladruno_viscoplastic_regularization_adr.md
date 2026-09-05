---
title: "ADR 90 — Rate (overstress) regularization for non-associated softening: scope, evidence, and the reopened D2"
project: Ladruno
type: ADR
status: "P0 + P0b complete; GATE U run (inconclusive on q_u, width non-convergent); D2 REOPENED; close-out PROPOSED at CP2; WP-C NOT opened"
priority: high
owner: nmora
related:
  - "[[_adr90_regularization_planning_brief]]"
  - "[[_adr90_orchestration_plan]]"
  - "[[_adr90_a0_results]]"
  - "[[_adr90_p0b_results]]"
  - "[[86_ladruno_sanisand_adr]]"
  - "[[31_ladruno_concrete3d_adr]]"
  - "[[59_ladruno_gradient_concrete_adr]]"
  - "[[32_ladruno_dispbeamcolumn_regularization_adr]]"
  - "[[73_ladruno_porous_overlay_adr]]"
  - "[[87_ladruno_depth_with_width_adr]]"
  - "[[testbed/00_canonical_testbed]]"
tags:
  - adr
  - material
  - wrapper
  - regularization
  - viscoplastic
  - overstress
  - localization
  - sanisand
  - tims
aliases:
  - ADR-90
  - LadrunoOverstress
  - duvaut-lions wrapper
updated: 2026-09-05
---

# ADR 90 — Rate (overstress) regularization for non-associated softening

> [!warning] Status — **P0 + P0b complete; GATE U run; D2 REOPENED; close-out PROPOSED at CP2; WP-C NOT opened**
> **D2 was decided at CP1 on 2026-09-04 on evidence since shown incomplete, and is REOPENED.**
> A three-lens adversarial pass (strategy / numerics / constitutive) found that P0's theorem
> carries an unwritten hypothesis and that three exposures had not been measured. Phase **P0b**
> measured them (`[[_adr90_p0b_results]]`, commits `892c22770` / `1052ecd36`). The four verdicts
> are in §4.5; they do not support the generic wrapper, and the P0b recommendation is
> **disclosure-only by default, WP-F conditionally, generic wrapper not at all** (§9 D2).
>
> **No C++ opens on this ADR.** WP-C is not opened; ND class tag **33022 remains RESERVED and
> absent from `SRC/classTags.h`**.
>
> **GATE U ran on 2026-09-05** (WP-A2, `Ladruno_implementation/_adr90_tau0_qu_band.md`) and its
> answer reframes the ADR: `q_u` is **not measurable** on the campaign's deck (every leg seizes
> inside the constitutive integrator), the localization **width does not converge** at τ = 0
> (0.675–0.824), and the **matched-settlement load — the quantity vault 65 D6 actually asks for —
> already contracts under refinement WITHOUT regularization** (7.80 → 5.28 → 2.86 %). §8.1
> therefore carries a **PROPOSED close-out**: P0 closes, width regularization is disclosure-only,
> and the actionable work is an ADR-86 integrator follow-up. The owner decides at **CP2**.
>
> **Number 90 on purpose.** 88 is cited throughout `SRC/` by PR #778; 89 is proposed for Track T.
> 33022 is free in the **ND** registry (highest is 33021, `SRC/classTags.h:594`) but is **live in
> code** as `PATTERN_TAG_LadrunoPorousOverlay` — used at `SRC/domain/domain/Domain.cpp:2274` — and
> as `EigenSOE_TAGS_FeastEigenSOE` (`:57`). Per-registry namespaces; deliberately not collisions.

---

## 1. Driver and problem — and the gate that decides whether there is a problem at all

The TIMs / APE macroelement campaign fits ultimate-surface loci from radial pushovers on
`LadrunoSANISAND`. Those loci live on the post-peak branch of a **non-associated softening**
material, where the rate-independent BVP loses ellipticity and the answer starts depending on the
mesh. The fork has no gate for this: every objectivity gate it ships regularizes **dissipated
energy** through a characteristic length `lch`, and one of them says so in its own docstring —
*"lch regularizes the dissipated ENERGY, not the localization WIDTH"*
(`tests/test_lemaitre_notched_bar.py:21-23`). Greps for `rudnicki|bifurcat|shear band|band width`
in `tests/` return nothing.

### 1.1 The size of the problem in 1-D — and its two caveats

A0 (`[[_adr90_a0_results]]` §5.1), 1-D softening bar, `τ = 0`, **one-element imperfection**,
`nsteps = 2000` at every mesh:

| N | h [mm] | w1 | w2 | w3 | w2/h | W₅₀ | W₅₀ ratio |
|---|---|---|---|---|---|---|---|
| 20 | 5.000 | 5.000 | 5.000 | 5.000 | 1.00 | 12.150 | — |
| 40 | 2.500 | 2.500 | 2.500 | 2.500 | 1.00 | 6.075 | 2.000 |
| 80 | 1.250 | 1.250 | 1.250 | 1.250 | 1.00 | 3.038 | 2.000 |
| 160 | 0.625 | 0.625 | 0.625 | 0.625 | 1.00 | 1.519 | 2.000 |

> **Two labels the earlier draft of this section did not carry.**
> (i) **The `W₅₀` column and the step-count column belong to different runs.** `W₅₀` above is the
> *one-element* convention at a fixed 2000 steps; the 4000 → 64000 step counts quoted below are the
> *graded fixed-length* convention. They are different decks and must not be read as one table.
> (ii) **Steps-to-converge is solver-dependent.** The graded-imperfection `τ = 0` problem needed
> 4000 / 16000 / 32000 / 64000 uniform steps at N = 20/40/80/160 while every regularized run
> finished at 250 on every mesh. That is a real and useful signature, but it is a property of *this*
> Newton implementation with *this* predictor and line search — a different algorithm would give
> different numbers. It is reported as a **symptom of ill-posedness**, not as a measurement of it.

### 1.2 The un-descope gate — is the wrapper needed at all?

**The need has been asserted, not measured.** ADR-59 died partly on that. Before P1 may open:

> **GATE U (un-descope).** Measure the **τ = 0 three-mesh band on the SANISAND /
> `LadrunoBrick -formulation bbar` deck** — the deck the campaign actually uses (S3, D7).
>
> **RUN AND REPORTED 2026-09-05** — `Ladruno_implementation/_adr90_tau0_qu_band.md` (WP-A2,
> branch `wp/90a2-tau0-qu-band`), engine `548fe911427e90a2edfead05cb3672a738d25b6d`, driver
> `Ladruno_files/testbed/hypo_bearing/sanisand_tau0_band.py`. Plane-strain strip footing,
> `h0 = 1.0 / 0.5 / 0.25` × `e_init = 0.6944` (Gorini) / `0.60` (dense), weighted soil, rough
> footing, patterned `sp` under the 1-argument `LoadControl` — i.e. **lane (b) of §3.1, already**.
>
> **VERDICT: INCONCLUSIVE on `q_u` (seizure); MEASURED on width (non-convergent) and on the
> matched-settlement load (convergent, contracting).**

**(i) `q_u` could not be measured at all.** All six legs seize far short of peak — the deepest
reached `s/B = 0.0228` of a `0.25` target (9 %), with its last four increments still rising
**54.3 / 53.8 / 53.5 kPa per 5 mm** — i.e. essentially **linear in settlement at 10.7 kPa/mm**,
with no sign of turning over. **No leg is a capacity**, so nothing from that WP may be quoted as a
collapse load, and this gate returns no `q_u` value to compare against a tolerance.

**And the seizure is not the stepping controller.** Every leg used **0 of its 80** pinned
subdivisions and ended with its step still **6400× / 12800× / 25000×** above the floor
(`h0 = 0.25 / 0.5 / 1.0`). The binding resource is the **constitutive integrator**: SANISAND's
substepped `ModifiedEuler` return collapses to `dT_min = 1e-6` (`ManzariDafalias.cpp:1380`) with
**no substep-count cap**, and the longest gap between two consecutive *converged* steps is
**2056 s = 34.3 min** — on `h0.5_e0.6944`, the one leg that terminated cleanly, **59 % of its
entire 3300 s budget went into that single step**. That is the boundary-value-problem form of §1.1's
steps-to-converge blow-up — ill-posedness presenting as **non-termination** rather than as a wrong
number.

**(ii) The localization width does NOT converge — C4's negative half, measured on the frozen
element.** `w2(h/2)/w2(h) = **0.675–0.824**` at every settlement and both densities (1.0 would be
converged, 0.5 a one-element band). Meanwhile the **yielding volume** is near-objective
(**3.95 %** at `s/B = 0.005`). *How much soil yields* is already nearly mesh-objective; *how thick
the mechanism is* is not — which is exactly why the fork's shipped `lch` energy gates cannot serve
this consumer (brief F6, §3.3).

*(All terminal numbers above are re-derived from the per-step curve CSVs, not from the leg JSONs:
a reaped leg's JSON stops at its last checkpoint and understated the deepest settlement by 12 %.
See `LEDGER_quirks`; corrected in WP-A2 `cc90f3d5e`.)*

**(iii) The matched-settlement load band CONTRACTS under refinement, without any regularization.**

| `s/B` | `e_init = 0.6944`: h1.0 / h0.5 / h0.25 | band | `e_init = 0.60` | band |
|---|---|---|---|---|
| 0.002 | 55.76 / 53.04 / 51.57 kPa | **7.80 %** | 84.12 / 79.68 / 81.78 | 5.42 % *(unresolved)* |
| 0.005 | 120.60 / 115.61 / 114.40 | **5.28 %** | 194.39 / 187.61 / — | **3.55 %** |
| 0.010 | 232.13 / 225.59 / — | **2.86 %** | — | — |

Monotone from above for Gorini's density. **Resolution caveats, from the WP's own controls:** the
`s/B = 0.002` row is *not* resolved (run-to-run scatter reaches **5.93 %** on relaxed ladder rungs,
against **0.000 %** at six shared clean checkpoints); the `0.005` and `0.010` rows **are**. The
solver-configuration floor is **0.8–1.4 %**, so no band below 1.4 % is a physical difference. The
**direction** is evidence; the **exponent** is not.

### 1.2.1 The reframing this forces — and the thing TIMs most needs to hear

Vault 65 **D6** already sets the campaign's ultimate criterion as **settlement-based**: when no
peak forms, take `q` at a fixed `s/B`. **Then the P1 deliverable is a matched-settlement load, not
a collapse load — and GATE U has just measured that quantity's three-mesh band CONTRACTING under
refinement with τ = 0: 7.80 % → 5.28 % → 2.86 %, monotone from above, on the campaign's own
material, element, mesh sequence and lane.**

Said plainly: **on the quantity the campaign will actually report, the unregularized problem is
already converging.** The regularizer is not needed to make that number stable. What is *not*
converging is the band **width** (0.675–0.824), and no deliverable in vault 65 is a width.

Two things this does **not** say. It does not say the band is *inside* the campaign's tolerance —
**OQ2 was never supplied**, so "inside tolerance" is unavailable by construction and this ADR does
not substitute a number of its own. And it does not say anything about `q_u`: a capacity was never
reached, so the collapse load remains unmeasured — **unreachable, not merely uncertain**.

### 1.2.2 Two deck defects GATE U had to fix first — both fork-wide

Neither is about this ADR, and both change how every existing fork SANISAND deck should be read
(details and reproducers in the WP-A2 note §5 and in that branch's `LEDGER_quirks` entries):

1. **`TanType` defaults to the ELASTIC tangent.** The parser default is `0`
   (`LadrunoSANISAND.cpp:117`, copied from `OPS_ManzariDafaliasMaterial`), so
   `ManzariDafalias3D::getTangent()` returns `mCe` — **every existing fork SANISAND deck runs
   `algorithm Newton` as modified Newton.** It is invisible on the zero-free-DOF material-point
   decks the fork's SANISAND tests all use, and costs ~7× on a BVP.
2. **A tight `NormDispIncr` is unreachable on SANISAND** — the substepped return makes the discrete
   stress–strain map only piecewise smooth, and the residual displacement norm stalls at a median
   **6.6e-7 m**. It is also **not mesh-neutral** (the norm runs over free DOFs, 3.6× more at
   `h0 = 0.25`), which disqualifies it inside a mesh-convergence study. Use `NormUnbalance`
   relative to the model's own weight.

Also from that WP, and a prerequisite for *any* deeper SANISAND push: at the ADR-86 default
**`-Presidual 0.0`** the dense coarse leg drives a free-surface Gauss point onto the low-`p` floor
at `s/B ≈ 0.0153` and starts logging *"the result at this integration point is set by the clamp,
not by the model"*. It is past every number reported, but a faster deck hits it next.

## 2. Named consumer, and what is settled

**Named consumer:** vault 65 **P1** — the ultimate-surface loci fitted from radial pushovers on
SANISAND. This ADR serves that consumer and no other; if P1 is withdrawn, this ADR is withdrawn.

| # | Settled (brief §1 — not re-opened) | Source |
|---|---|---|
| S1 | Method = a rate-regularizing **wrapper at the `NDMaterial` level**, not a modified `ManzariDafalias`. | vault 64 §6; vault 65 D5 |
| S2 | **τ is a declared numerical parameter, not a soil property**; the deliverable is *"mesh-independent given a declared τ"*, characterised on a Deborah number; uniqueness is **not** restored. | vault 10 line 75; vault 65 D3 + D5 |
| S3 | Element frozen: `LadrunoBrick -formulation bbar`; tetrahedra prohibited on failure legs. | vault 65 D4; vault 71; `tests/test_r3_prandtl_collapse_gate.py` (#722) |
| S4 | **Provenance is output, not documentation.** | vault 65 D3; `ops.ladrunoBuild()` (#718) |
| S5 | Vary the relaxation parameter and the ramp duration **independently**; matched pairs collapse onto the ratio. | vault 14 |
| S6 | Primary instrument = the displacement-controlled radial monotonic curve; radial probes, **sequential swipe** as cross-check. | vault 65 D2, D7 |
| S7 | The wrapper is the smaller half; the validation campaign is the larger half. | vault 64 §6 |
| S8 | Papers in hand; no acquisition phase. | skill `tim-macroelement/references/library_map.md` |
| S9 | The staggered u–p overlay's **amended clause 2** governs how a rate-dependent constitutive law composes with the fixed-stress split — any transient lane here inherits it. | vault 73, amended clause 2; [[73_ladruno_porous_overlay_adr]] |

> **S6 is in tension with this ADR's own model.** A *sequential swipe* is by construction a
> reversal path, and §4.3/§4.5 show the wrapper's declared domain excludes reversals. The swipe
> cross-check therefore lands exactly where the formulation is weakest. Recorded here rather than
> buried in a risk row.

## 3. The quasi-static contradiction, and the lane decision

Vault 65 makes the displacement-controlled radial probe the primary instrument (D2, D7) and rate
regularization the mechanism (D5). Four measured tensions:

1. **There is no intrinsic length in quasi-statics** (de Borst et al. 1993). A0 quantifies it: the
   converged band width moves a **factor 12 for a factor 3 in De**. P0b-(c) adds that it also moves
   **×4.3 with the imperfection amplitude and ×3.5 with the zone length at fixed De** — so the
   width is not a τ property at all (§4.5 V3).
2. **`ops_Dt` is not uniform under the campaign's own integrator.** It is `Domain::dT` = current −
   committed pseudo-time (`Domain.cpp:2080-2082`, `:2125`, `:2392`); under `DisplacementControl`
   the λ-increment swings by ~1e5 and goes negative after `loadConst`.
3. **β changes per Newton ITERATION under `DisplacementControl` and `ArcLength`.**
   `DisplacementControl::update()` (`SRC/analysis/integrator/DisplacementControl.cpp:298`) calls
   `applyLoadDomain(currentLambda)` at **`:346`**, and `ArcLength::update()` (`:226`) calls it at
   **`:302`** — i.e. the pseudo-time (and therefore `ops_Dt`, and therefore
   `β = Δt/(τ+Δt)`) is re-set *inside* the iteration, not once per step. A residual that depends on
   the iteration number is not a function of `u`, and Newton's convergence argument is void.
4. **The transient is Δt-dependent even at fixed De** (A0 §5.5: `w2` moves 2.06 % over a 32× step
   refinement), and **ellipticity is a β criterion, not a De criterion** (§4.5 V4).

### 3.1 Alternatives table

| Option | What it means | Cost | Verdict |
|---|---|---|---|
| (a) Transient lane, τ replaces β_K | Slow damped transients with Rayleigh **β_K = 0** and τ as the one declared relaxation time. Δt physical, De a deck quantity. | Wall clock (vault 14: ~17 min/run). **And an unpriced one:** removing β_K removes the only damping on the high-frequency content the softening branch injects — vault 14's β_K was doing stabilization work that τ does not replace, so undamped ringing at β_K = 0 is a cost nobody has costed. | **SECONDARY** (demoted 2026-09-05) |
| **(b) Uniform-pseudo-time static lane** | Prescribe the pushover displacement as an **`sp` inside a load pattern** and drive it with the **1-argument `LoadControl(dλ)`**. `SP_Constraint::applyConstraint` sets `valueC = loadFactor * valueR` when the constraint is not constant (`SRC/domain/constraints/SP_Constraint.cpp:331-337`), so a patterned `sp` under a uniform λ-increment **is displacement control with exactly uniform pseudo-time**: `Δt = dλ` every step, `β` constant, `De` a deck quantity. It keeps limit-point capability, keeps vault 60's D16 stepping guards, and keeps comparability with the vault's existing runs. **The 4-argument `LoadControl(dλ, numIncr, min, max)` adapts `dλ` and destroys the uniformity — it must not be used.** | Re-verify D16's guards; the D16 adaptive subdivision still varies Δt (see m11 below). | **PRIMARY** (promoted 2026-09-05) |
| (c) Strain-driven internal clock | Replace Δt by an accumulated strain measure. | A different constitutive model (Perzyna-in-arclength); no literature anchor, no closed-form gate. | **REJECTED** |
| (d) Do nothing to the lane; report the three-mesh band | Run as-is and disclose the spread (vault 10's stance). | Zero engine work. | **The NEGATIVE CONTROL — and, after P0b, the recommended default (§9 D2).** |

**Consequence for the wrapper (numerics B3):** the wrapper must **latch β and Δt at `newStep()` /
the first trial evaluation of a step and hold them for the whole iteration**, rewinding on
`revertToLastCommit`; or it must **hard-refuse** `DisplacementControl` and `ArcLength`. Only the
1-argument `LoadControl` is uniform without a latch.

### 3.2 The headline claim, worded to survive review

> *The regularized quantity `q(De, Δt; h)` converges in `h` at a **declared (De, Δt, imperfection
> field)**, and its dependence on each of those three is measured and disclosed.*

Not "mesh-independent collapse load"; not even "converges at fixed De". A0 §5.2 shows curve
convergence and width convergence are different gates (at De = 1e-4 the curve converges and the
width does not), and P0b-(c) shows the converged width is a function of the imperfection field as
well as of De.

### 3.3 Reconciliation with the fork's shipped alternatives, and the honest-framing test on τ

| Shipped alternative | What it regularizes | Why it does not serve this consumer |
|---|---|---|
| **Crack-band `lch`** | Dissipated **energy** per unit band area. | Leaves the band one element wide; its own gate says so (`tests/test_lemaitre_notched_bar.py:21-23`) and §1.1 is that statement in numbers. |
| **`LadrunoConcrete3D -eta`** (the fork's shipped rate regularization, `SRC/material/nD/LadrunoConcrete3DKernel.h:1492-1503`) | The same blend, but **inside** the CDPM2 kernel on the effective stress and `κ_p`. | Needs the committed internal effective stress and `κ_p` — nothing a generic seam exposes. Tier-1 only (`!mp.implex`). And it is a *damage* model: ADR-31 §4.4 refused a nominal-stress blend because it relaxes damage as well as plasticity (inherited here as **D1**). |
| **ADR-32 embedded discontinuity** | Localization in a **frame element**, by a kinematic enhancement. | 1-D element kinematics; no 3-D continuum version exists. |
| **ADR-59 gradient / non-local** | Width, correctly, by a true internal length. | **Descoped.** Its kill list is the standard this ADR is written against. |

**Conclusion:** nothing shipped regularizes localization *width* for a 3-D plasticity-type
continuum material. This ADR is not "mostly re-wiring".

> **The honest-framing test, applied to τ itself (ADR-59's own instrument, turned inward).**
> ADR-59 was descoped partly because its internal length `ℓ` was a free parameter that would be
> tuned to the answer. **τ is in exactly that position, and worse:** it is unbounded above, it moves
> the band width ×12 for ×3, and §7.5 asks TIMs to *supply* the De they will run at. A collapse load
> obtained that way is **calibrated, not predicted**. Therefore:
> - **The guide must forbid tuning τ (or De) against a target `q_u` or a target band width.** τ is
>   declared *a priori* from the loading rate and the intended Deborah number, and then reported.
> - **The honest deliverable is `q_u(De)` as a curve with a declared floor** (§7.4's De rule), not
>   a single regularized number.
> - Failing that, this ADR repeats ADR-59's failure with a different Greek letter.

## 4. Formulation

### 4.0 What this model actually is — and what it must therefore be called

The update in §4.1 relaxes the stress toward `σ̄(ε)`, a **function of the strain history evaluated
by a separate, rate-independent material**. That is an **overstress model on a rate-independent
backbone** — the Maxwell / Krempl viscoplasticity-based-on-overstress (VBO) family,

```
    sigma_dot = C_e : eps_dot - (sigma - sigma_bar(eps)) / tau
```

— and it is **not** Duvaut–Lions. In Duvaut–Lions the relaxation target is the **projection of the
current state** `(σ, q)` onto the elastic domain, which is why the model retains an elastic domain
and a normal-cone structure. The two-track wrapper's target is a strain-history function, so:

- it has **no elastic domain after first yield** — every step relaxes, including "elastic" ones;
- its `σ − σ̄` is **not** a normal-cone direction, which is exactly why the dissipation gate fails
  (§4.5 V2);
- the two coincide **only** under the hypotheses in §4.3, all of which SANISAND violates.

**Therefore the class must not be called `LadrunoDuvautLions`.** If anything is built, it is
`LadrunoOverstress` (D10). Reserving the Duvaut–Lions name for an implementation that actually
performs the projection keeps WP-F's name free and stops the ADR's own title from making a
correctness claim.

### 4.1 The two-track update, with the off-switch as an early-out

```
    if (tau <= 0.0 || dt <= 0.0) {                 # the fork convention -- an EARLY-OUT, not beta=1
        inner->setTrialStrain(eps_{n+1});          # forward, then COPY the inner's own answer
        sigma_vp = copy(inner->getStress());       # byte-identical by CONSTRUCTION, not by algebra
        D_vp     = copy(inner->getTangent());
        return;
    }
    C_e       = <the elastic operator, see D3>     # deep-copied at the WRAPPER's commitState
    sigma_bar = copy(inner->getStress());          # inner runs INVISCIDLY on the total strain path
    C_bar     = copy(inner->getTangent());
    sigma_tr  = sigma_vp,n + C_e : delta_eps
    beta      = dt / (tau + dt)                    # LATCHED at newStep -- see 3.(3)
    sigma_vp  = (1 - beta) * sigma_tr + beta * sigma_bar
    D_vp      = (1 - beta) * C_e      + beta * C_bar
```

**The early-out is required, not stylistic** (numerics B2). The shipped `-eta` obtains its byte
identity by *skipping the block entirely* — `if (!mp.implex && mp.eta > 0.0 && dt > 0.0) { … }`,
`SRC/material/nD/LadrunoConcrete3DKernel.h:1492`. The earlier draft of this ADR wrote the blend
unconditionally and relied on `(1−1.0)·σ_tr + 1.0·σ̄` being exact. It is exact for finite values —
but **`0.0 * NaN = NaN` and `0.0 * Inf = NaN`**, so an inner that returns a non-finite trial (a
diverging Newton iterate, an unconverged substep) would have that poison multiplied by zero and
*survive* as NaN in the wrapper where the inner alone would have been recoverable. The gate is
therefore restated (C1/D9) as **bit-identical stress, tangent AND committed-state trajectory over a
multi-step path**, with an explicit NaN-propagation case, not as an "instruction path" claim that
no test can express.

### 4.2 `τ ≤ 0` or `Δt ≤ 0` ⇒ inviscid — and why "inert" was the wrong word

The fork gates to the **inviscid** branch, never to elastic, whenever `τ ≤ 0` or `Δt ≤ 0`
(`LEDGER_quirks.md:1612`; `LadrunoConcrete3DKernel.h:1486-1488`).

**This is not "inert" (constitutive 9).** Taking the inviscid branch for one committed step **dumps
the entire accumulated overstress in that step**: `σ_vp` jumps from `σ̄ + Δσ_over` to `σ̄`, a finite
stress drop with no strain increment. Three routes hit it in normal use:

- `loadConst` makes the pseudo-time increment negative;
- a **held-load** step (`analyze` with no load increment) has `Δt` positive but no strain rate, so
  the model relaxes fully toward the inviscid backbone with time constant τ — **staged geostatic
  steps do exactly this**;
- **`Domain::revertToLastCommit()` sets `dT = 0.0` and re-applies the load**
  (`SRC/domain/domain/Domain.cpp:2334-2339`), so **the first evaluation of every retried step after
  a step cutback runs with `β = 1`.**

The wrapper must therefore **count the steps that committed with `τ > 0` but `β = 1`, and report
that count in the provenance block. A non-zero count fails acceptance** — it means the run silently
mixed regularized and unregularized steps.

### 4.3 The identity, and the three hypotheses it needs

> **Theorem (1-D, from rest, monotonic loading, CONSTANT elastic operator).** For the 1-D
> associated model with any hardening function `K(α)`, the two-track update and true Duvaut–Lions
> produce the identical stress path.
>
> *Proof.* The TDL update gives `σ_{n+1} + Eα_{n+1} = σ_tr + Eα_n`, so `ψ = σ + Eα` advances by
> exactly `EΔε` every step. From rest `ψ_0 = 0`, hence **`α_n = ε_n − σ_n/E`**. Substituting into
> the projection equation yields `σ̄ = K(ᾱ)` with `ᾱ = ε_{n+1} − σ̄/E` — the definition of the
> inviscid stress at `ε_{n+1}`. So TDL's projection target *is* `inner.getStress()` on the total
> strain path. ∎

**The hypotheses, all three of which must be stated wherever the theorem is:**

1. **A constant elastic operator.** `ψ = σ + Eα` cannot advance by `EΔε` unless E is the same number
   on both tracks. Added 2026-09-05 after P0b-(a). **SANISAND's `G(p)`, `K(p)` violate it by
   construction**, as does Drucker–Prager with `G(p)` and every hypoelastic inner.
2. **Holonomicity of the inner** — no memory beyond `(ε, σ)`. This is the sharp form of what the
   earlier draft called "proportional and monotonic": what the proof actually needs is that the
   inner's internal state be *recoverable* from the current strain and stress. Linear **kinematic**
   hardening is holonomic even on non-proportional paths; **isotropic** hardening is holonomic only
   on proportional monotonic ones. **SANISAND is holonomic nowhere**: `mAlpha` follows a rate law,
   `mFabric` accumulates, `mAlpha_in` is defined by the last reversal; only the void ratio is
   recoverable.
3. **No unloading.** The relation `α = ε − σ/E` survives elastic steps in 1-D, but the *inner*'s
   response after a reversal is not a function of the current `(ε, σ)`.

**Hypothesis 3 is fatal for the deliverable, not merely limiting.** The boundary of a localization
band is an **elastic unloading zone** (Rice) — that is the mechanism that selects the band's width.
So the wrapper's declared validity domain **never holds where the answer is set**.

### 4.4 The closed-form 1-D anchor — and its scope

The steady overstress of the discrete backward-Euler fixed point is `σ* − σ_Y = E·ε̇·τ`, exactly and
independently of Δt, measured to `4.26e-15` over `Δt ∈ {1, 0.25, 0.05}`. **This closed form is
derived for PERFECT plasticity** (the backbone is a constant `σ_Y`); for a hardening or softening
backbone the same algebra gives the overstress *above the current backbone*, which is not a
closed-form function of time. Cite it as the perfect-plasticity anchor only. It is the ADR's
non-self-referential oracle for the ADR-87 verification manifest.

### 4.5 P0b — the four measurements, and what they did to §4.3

`[[_adr90_p0b_results]]`, git `892c22770`, 2026-09-05.

| # | Measurement | Result |
|---|---|---|
| **V1** | Generic two-track over a **pressure-dependent inner**, `E(σ) = E₀(1 + β_E|σ|/σ_Y)`, monotonic proportional isotropic hardening | rel `|σ_TT − σ_TDL|` = **9.2e-4 … 7.3e-3** at the working De (β_E = 0.3…1.2), rising to **0.28…0.94** at De = 0.3. Constant E: **3.0e-14**. The obvious repair — the predictor modulus at the wrapper's own stress — buys 1–2 orders and **is not implementable at the seam** (`getInitialTangent()` returns one sampled matrix, not the function `E(·)`). **⇒ the identity is exact only for a constant elastic operator; on the named consumer's material class the wrapper is inexact on EVERY path.** |
| **V2** | **Dissipation** `D_w = σ_vp · C_e⁻¹(σ_vp − σ̄)Δt/τ` on non-proportional J2, 1-D cyclic, and a swipe surrogate | **VIOLATED on load/unload/reload**: worst step **−2.26e-5**, worst cumulative-negative fraction **−7.9e-3** at De = 0.1, growing monotonically with De (round-off at 1e-3, −1.4e-5 at 1e-2). Non-proportional J2 and the swipe surrogate are clean. Every run is net-positive overall and the inner's own dissipation is positive, so **an energy-balance check cannot see this**. |
| **V3** | **Is the width τ-set or imperfection-set?** At De = 3e-4, amplitude and zone length varied | The converged `w2` moves **×4.28 with amplitude** and **×3.47 with zone length** at fixed De; `w2` is mesh-converged in every case. **⇒ imperfection-set.** Also: `w3` (FWHM) does **not** converge outside A0's original parameter point (2.500 → 1.875, tracking `h`, in 3 of 9 configurations) — A0's "both w2 and w3 converged at De = 3e-4" was parameter-lucky. And **`w2 ∈ [h, L]` by construction**, so a `w2` near L means *no band*, not *a wide band*. |
| **V4** | **Blended acoustic tensor** `det[n·((1−β)C_e + βC^ep)·n]` on a plane-strain **non-associated softening Drucker–Prager** point model | The inviscid tangent **does** lose ellipticity (min normalised det Q = **−0.0145** at **51.5°**, from step 75/800). The blend regains it for **β < 0.9857**, i.e. **De > 1.45e-2 / nsteps**. Every De in the working window clears it with margin. **But the criterion is on β = Δt/(τ+Δt) — the same De is elliptic or not depending on the step count.** Also: TT ≡ TDL to 2e-14 on a **proportional** path over a non-associated model (non-associativity alone does not break the theorem); **path rotation** does (1.1e-4 at De 1e-4 → 7.5e-2 at De 0.1). |

**V4 is the one positive result**: the mechanism does restore well-posedness at the material point
for the consumer's material class, which is the thing A0's 1-D bar could not see. **V1, V2 and V3
are why D2 is reopened.**

## 5. Architecture — *if* anything is built

Retained because it is the reviewed record of what a correct implementation must do; **it is not an
instruction to build.** ND tag **33022**, modelled on `StagedStrainNDMaterial` (33014).
Command (flags after positionals, ADR-86 rule): `nDMaterial LadrunoOverstress $tag $innerTag -tau $tau`.

| Item | Rule | The source that forces it |
|---|---|---|
| Adopting ctor | Take the inner by tag, `getCopy("ThreeDimensional")` once, own it, delete in the dtor. | `StagedStrainNDMaterial` pattern. |
| `setTrialStrainIncr` | Keep an **absolute** committed strain; reconstruct `ε_{n+1} = ε_committed + Δε`. Never forward an increment to the inner. | The blend needs the total strain for the inner and the increment for the predictor; the rate overloads discard the rate (`ManzariDafalias3D.cpp:88-91`). |
| `getCopy(const char*)` | Resolve **every** supported type explicitly (`"ThreeDimensional"`, `"3D"`, `"PlaneStrain"`, `"PlaneStrain2D"`), build the inner's view with `inner->getCopy(type)`, wrap it, return `0` gracefully otherwise. **Never** route through the inner's `getCopy(void)`. | `StagedStrainNDMaterial::getCopy(const char*)` (`SRC/material/nD/StagedStrainNDMaterial.cpp:293-309`) has **no** `"ThreeDimensional"` special case; only `InitStrainNDMaterial.cpp:290-292` has one. **Corrected rationale (numerics m12):** `getCopy(void)` *is* overridden in both concrete views (`ManzariDafalias3D.h:49`, `ManzariDafaliasPlaneStrain.h:52`) — the base-class `exit(-1)` is not the live hazard. The real reason is that `getCopy(const char*)` **constructs from parameters and carries no state**, so routing a *stateful* copy through it silently drops the committed state; and the plane-strain view's null ctor writes the wrong classTag (§10 OQ7). |
| *(do not lean on this)* | `FluidSolidPorousMaterial::getCopy(const char *code)` (`SRC/material/nD/soil/FluidSolidPorousMaterial.cpp:415-425`) **ignores** `code` and routes through the copy ctor, whose `:156` calls the soil's `getCopy(void)`. The only plane-strain path in the tree that reaches a wrapped material's void `getCopy` — a coincidence, recorded as a quirk, never a dependency of this wrapper or its tests. | |
| **`C_e` for the predictor (D3 — CORRECTED)** | **Deep-copy `C_e` inside the WRAPPER's `commitState()`, immediately after `inner->commitState()`.** Not "read before `setTrialStrain`". | **The earlier rule was wrong.** `mCe` is rewritten *inside* `ManzariDafalias::integrate()` — the integrator calls at `SRC/material/nD/UWmaterials/ManzariDafalias.cpp:979 / :987 / :992` all pass `mCe` by reference — so it carries the **previous Newton iterate's** value, not the committed one. A predictor built on it makes the residual a function of the iteration history rather than of `u`, which breaks Newton. **Gate: `Newton` and `ModifiedNewton` must produce bit-identical committed states on a fixed path.** |
| **Copy, never alias — corrected citation (numerics M5)** | Deep-copy every `Vector`/`Matrix` the inner returns, **before any other call on any instance of that class**. | The earlier citation was wrong in the specific. `ManzariDafalias3D.h:74-79` declares only `mSigma_M` and `mEpsilon_M` static — its **tangents are per-instance**. It is the **plane-strain** view that is dangerous: `ManzariDafaliasPlaneStrain.h:78-82` declares `static Matrix mTangent;` and `static Matrix mTangent_init;`, and `getInitialTangent()` **writes into the shared `mTangent_init`** (`ManzariDafaliasPlaneStrain.cpp:157-171`). That is precisely what D3 reads, on precisely leg A2. **Gate: a two-instance interleave test on plane strain** — construct two wrapped plane-strain materials, interleave their `getInitialTangent()` calls, and require both to see their own moduli. |
| **`setParameter` (M6)** | Match `"tau"` **only** when `argc > 1 && atoi(argv[1]) == this->getTag()`. On a miss, forward to the inner **with `argv[1]` rewritten to the inner's tag** (or accept both tags). | `ManzariDafalias::setParameter` (`:791-799`) reads `theMaterialTag = atoi(argv[1])` and does nothing unless it equals its own tag. Forwarding the wrapper's argv verbatim therefore **cannot reach the inner** and the static `mElastFlag` never flips — a staged geostatic run stays elastic, silently. The design's reliance on that **static, class-wide** flag must be documented: it is shared by every `ManzariDafalias` instance in the process. |
| τ as a `Parameter` | `setParameter("tau")` / `updateParameter`, so De sweeps need no reconstruction. | D6. |
| `revertToLastCommit` (M10) | The **inner's** stub is benign — `integrate()` restarts from the `_n` state each call. The real exposures are the **wrapper's own**: it must rewind the cached `C_e`, the committed `σ_vp`, and the latched Δt/β. | `ManzariDafalias.cpp:513-517`; `Domain.cpp:2334-2339`. |
| Response tokens | `overstress`, `beta`, `dt` — with **`dt` exposing the whole Δt series**, not a single number. | S4; the Δt diagnostic below. |
| **Δt diagnostic (M8)** | A **file-scope latch keyed on `ops_Dt` transitions**, reset on `revertToStart`. Report **min/max Δt ratio and the D16 subdivision count** in the provenance block. Do **not** try to name the integrator. | A per-instance warn-once fires once per material and floods; the integrator is not identifiable from inside a material. **And m11: D16's adaptive subdivision makes Δt non-uniform on lane (b) too**, so this is a reporting obligation, not an exception. |
| `sendSelf` / `recvSelf` | `ID(4)` = `{tag, inner classTag, inner dbTag, nstate}` → `Vector` of wrapper state (τ, committed `σ_vp`, committed `ε`, latched Δt, the `β = 1`-with-`τ > 0` counter) → delegate to `inner->sendSelf`. | `StagedStrainNDMaterial.cpp:316-…`. |
| Registration | ×3 (Tcl, functionMap, `FEM_ObjectBroker`) + CMake. | The Tcl/Python double-dispatch trap is the first entry in `LEDGER_quirks.md`. |
| **Inner-tangent type (constitutive 11)** | Refuse, or loudly warn, when the inner returns a **continuum** rather than an algorithmic tangent — `ManzariDafalias` with `mTangType = 0` returns `mCe`, so the "blend" degenerates to `C_e` and the wrapper is elastic in the tangent while inelastic in the stress. | C3's claim assumes an algorithmic tangent. |
| **Scope guard (D1, constitutive 11)** | Damage / plastic-damage inners are out of scope. They cannot be *detected*, but the fork's own can be **refused by a classTag denylist**: `ND_TAG_LadrunoConcrete3D` (33017), `ND_TAG_LadrunoRCConcrete` (33015), and `ASDConcrete3D`. | Cheap, and it removes the three most likely misuses. |
| Unsymmetric solver | For a non-associated inner the blended tangent is **unsymmetric**; the guide and the decks must require a general solver. | ADR-31 R13 precedent; V4's model is non-associated by construction. |

## 6. Claims, each with its gate — and what is explicitly not claimed

| # | Claim | Gate | Status after P0b |
|---|---|---|---|
| **C1** | `τ = 0` reproduces the inner **bit-identically in stress, tangent and committed-state trajectory** over a multi-step path, including a **NaN/Inf** case | early-out + forward + copy (§4.1); trajectory gate at the material point and on `LadrunoBrick`; g++ byte check | restated (was "instruction path") |
| **C2** | 1-D steady overstress `= E·ε̇·τ`, Δt-independent, **for PERFECT plasticity** | PV3 clone; `4.26e-15` measured | scoped |
| **C3** | Tangent `= (1−β)C_e + β C_alg`, matches a central FD — **provided the inner returns an algorithmic tangent** | PV4/PV5 clones + an inner-tangent-type guard | scoped |
| **C4** | Band width converges under h-refinement at a declared **(De, Δt, imperfection field)**, and does not at τ = 0 | §7 width gate + the §7.4 De rule | **negative half MEASURED** by GATE U (`w2(h/2)/w2(h) = 0.675–0.824` on the frozen element, both densities, every settlement) — but **weakened by V3**: the positive half is imperfection-set. GATE U's ratios are also **pre-peak**; a post-peak width is owed |
| **C5** | ~~`q_u` at fixed De converges in h~~ → **the matched-settlement load `q(s/B)` at fixed De converges in h — and `q(De)/q(τ=0)` is reported per leg** | slow-tier gate; capacity rule | **RE-BASED and largely pre-empted.** No capacity is reachable on this deck (GATE U (i)), so under vault 65 D6 the deliverable is `q` at fixed `s/B` — and GATE U measured **that** band already contracting at τ = 0 (7.80 → 5.28 → 2.86 %). C5 therefore has little left to add, and the `q_u` upward bias below applies to whatever replaces it |
| **C6** | The De-dependence of width and `q_u` is measured and disclosed, **at different step counts** | vault-14 sweep; matched pairs with different `nsteps` | unchanged |
| **C7** | Transparent to `updateMaterialStage`, `getCopy(type)`, database round-trip, MP wire | Zone-A + roundtrip + `test_fspm_over_manzari_family` twin + the **two-instance plane-strain interleave** + **Newton vs ModifiedNewton bit-identity** | strengthened |
| **C8** *(new)* | The **blended acoustic tensor is elliptic** at the declared (De, Δt) at every committed state of the A1/A2/A3 legs, and the minimum normalised det Q is reported | V4's metric, run on the deck | new gate |

**NOT claimed.**

1. No restoration of uniqueness or the bound theorems (S2).
2. No intrinsic length in statics — corroborated and quantified by A0 and V3.
3. No mesh independence without a declared **(De, Δt, imperfection field)**.
4. No applicability to damage or plastic-damage inners (D1).
5. No change to the tetrahedron prohibition (S3).
6. **No correctness claim of any size is established for the named consumer's material.** The
   identity of §4.3 needs a constant elastic operator and a holonomic inner; SANISAND satisfies
   neither, and V1 measures the resulting error on a caricature only.
7. **The regularized `q_u` is UPWARD-BIASED, and the bias is of the same order as the artefact the
   campaign already removed.** From A0's own peak loads against the rate-independent 18.0 N:
   **+9.7 % at De = 3e-4, +12.5 % at De = 1e-3, +17.0 % at De = 3e-3.** Vault 14's β_K artefact,
   which the campaign deliberately eliminated, was **11.75 %**. Substituting a declared overstress
   for an accidental one does not make the strength gain physical. Every regularized `q_u` must be
   reported **beside** `q_u(τ = 0)` and as the ratio.
8. **No claim that one global τ preserves the shape or tilt of the ultimate-surface ellipsoid.**
   Under a direction-dependent strain-rate field a single τ produces a direction-dependent
   overstress, which distorts the fitted surface — vault 14's "what it does not establish", vault
   65 D9. A WP-D leg on **≥ 2 probe directions** is required before any fitted surface is cited.
9. **No claim of non-negative incremental dissipation** (V2).
10. **No claim that regularization is what stands between the campaign and a mechanism.** GATE U
    measured the binding obstacle, and it is **the SANISAND constitutive integrator, not
    ill-posedness**: a substepped `ModifiedEuler` that collapses to `dT_min = 1e-6` with **no
    substep-count cap** (`ManzariDafalias.cpp:1380`), an elastic-tangent parser default that makes
    `Newton` behave as modified Newton (`LadrunoSANISAND.cpp:117`), and a convergence test that
    cannot be met on this material. Every leg used **0 of 80** subdivisions and stopped
    **6400–25000×above** the step floor. Adding a relaxation time to a material that spent
    **2056 s (34.3 min) inside one converged step — 59 % of that leg's whole budget** — does not
    make the mechanism reachable; the integrator work does.
    Any statement that the wrapper "lets the problem finish" must be measured against a deck with
    that work already done, not against this one.
11. **No claim on `q_u` in value, from any source.** GATE U reached no capacity on any leg; the
    unregularized collapse load is **unreachable, not merely uncertain**, and nothing in this ADR
    quotes one.

## 7. Acceptance case

| Leg | Deck | Material | Element | Status |
|---|---|---|---|---|
| **A0** | 1-D softening bar, N = 20…160 | numpy oracle | — | **DONE** (`cc7c7f7a5`, `8863a468c`) |
| **P0b** | point models + DP + acoustic tensor | numpy oracle | — | **DONE** (`892c22770`, `1052ecd36`) |
| **GATE U** | τ = 0 three-mesh `q_u` band | `LadrunoSANISAND` | B-bar hex | **PENDING — §1.2; blocks P1** |
| **A1** | Plane-strain biaxial, symmetric half, smooth ends, graded imperfection field; 3 meshes + 2 orientations | `DruckerPragerPlaneStrain` (non-associated) | `quad PlaneStrain` | P2 |
| **A2** | Same deck | `LadrunoSANISANDPlaneStrain` (33021) | `quad PlaneStrain` | P2, with the OQ7 and M9 constraints |
| **A3** | Same specimen, one element thick | `LadrunoSANISAND` | **`LadrunoBrick -formulation bbar`** | P2 — the only admissible leg |
| **A4** | R3 Prandtl–Reissner strip | SANISAND / DP | B-bar hex | P2 regression — **as a two-sided gate, see 7.1(v)** |

### 7.1 Corrections to the acceptance protocol, from A0 and P0b

**(i) The imperfection must be a mesh-convergent FIELD — and now also a *declared* one.** A
one-element defect never converges in width; a *flat* fixed-length zone makes every mesh exact
(a gate that cannot fail); a **graded** fixed-length notch converges. **P0b-(c) adds the harder
requirement:** because the converged width moves ×4.3 with amplitude and ×3.5 with zone length, the
imperfection field is a **first-class declared input of every reported number**, on the same footing
as τ. It goes in the provenance block.

**(ii) The De window.** The brief's `De ∈ {0.01, 0.1, 1}` is unusable on the bar deck (at De = 1 the
load never peaks; at 0.01 and 0.1 the "band" is the whole specimen). A0 found `3e-4 … 1e-3`
on that deck — and P0b narrows even that: **only De = 3e-4 converged in both `w2` and `w3`, and
P0b-(c) shows that reading was parameter-lucky.** The window is deck-dependent; **the SANISAND-deck
window must be MEASURED in P2** by the rule in §7.4.
> **And the recipe presupposes material softening (constitutive 11).** "Sweep De down until a
> post-peak branch appears" works for the softening bar. **Drucker–Prager localizes on a *hardening*
> branch** (Rudnicki–Rice), so A1 and A2 need a different window criterion: use **C8's acoustic
> tensor** — the smallest De whose blended tangent stays elliptic along the whole path — not the
> appearance of a load peak.

**(iii) The negative control — primary and secondary (numerics M7).** The width-∝-h control needed
4000 → 64000 uniform steps at τ = 0 and is **probably unrunnable at the slow tier on a 3-D deck**.
Therefore:
- **PRIMARY negative control: steps-to-converge.** The τ = 0 deck must fail to complete, or complete
  only at a step count that grows with refinement, while the regularized deck completes at a
  mesh-independent step count. Cheap, and it fails loudly.
- **SECONDARY: width ∝ h**, run at whatever mesh range the budget allows.

**(iv) The De collapse must use DIFFERENT step counts.** `β = 1/(1 + nsteps·De)` and the strain
increment is `u_max/nsteps`, so equal `(De, nsteps)` runs are **bit-identical** whatever `(τ, T)`
produced that De — measured to 12 digits. With different step counts the collapse holds to 0.47 %
in width, 0.009 % in peak load.

**(v) A4 must be a two-sided gate (strategy M6).** "Unchanged within the De family" is not a gate.
A4 passes only if **(a)** at `τ = 0` the result is bit-identical to the pre-wrapper R3 number, **and
(b)** the shift at the declared De **matches §4.4's closed-form overstress prediction to a stated
tolerance**. A shift of the right sign but the wrong size is a defect, and the current wording
cannot see it.

**(vi) Staging (strategy m10).** The wrapper must be **`τ = 0` during the staged geostatic phase**,
and a **post-gravity byte-identity gate** must show that the staged state with the wrapper present
equals the staged state without it. Otherwise §4.2's relaxation-under-held-load silently changes
the initial stress field the whole campaign is fitted to.

**(vii) The imperfection field on the real material (numerics M9).** On SANISAND the physically
meaningful imperfection knob is the **initial void ratio `e_init`** — a *construction* parameter.
The `voidRatio` `Parameter` is **tag-addressed**, so a spatially varying field needs **one material
tag (and one wrapper tag) per element**. Name the perturbed variable and the field form in the
deck, and either **price the tag explosion** (a 3-mesh × 2-orientation A2/A3 study is thousands of
tags) or add an **`-ele`-routed field** as a **P2 prerequisite**. This is not a detail: without it
requirement (i) cannot be met on the consumer's own material.

### 7.2 Measurements per leg

Band width `w2` (threshold-free), `w3`, **`w2/h`**, peak load, `q_u`, **`q_u(De)/q_u(τ=0)`**,
post-peak curve, dissipated work, **min normalised det Q (C8)**, `De`, **`Δt` min/max ratio and
subdivision count**, **the `β = 1`-with-`τ > 0` step count**, the imperfection field, and
`ladrunoBuild()`.

### 7.3 The width metric

`w2 = √(12·Var)`, `Var = Σ p_e[(x_e − x̄)² + h²/12] / Σ p_e`, over the post-peak plastic-strain
increment. A one-element band reads exactly `h`; a k-element top hat reads exactly `k·h`. Unit-pinned.
**Two limits must always be quoted with it:** `w2 ∈ [h, L]` **by construction** — a `w2` near `L`
means *no band*, not *a wide band* — and threshold metrics (`w1`, `w3`) disagree by up to 40× on the
same run, so they are never the headline.

### 7.4 Gates — including the De RULE

- **τ = 0**: primary and secondary negative controls per 7.1(iii) must fail objectivity.
- **The De rule (strategy M4 / constitutive 6).** The limits `h → 0` and `De → 0` **do not commute**:
  refining the mesh at fixed De eventually resolves a band that De alone would not have set, and
  lowering De at fixed mesh eventually returns the one-element band. So De is not a free choice —
  it is a function of the mesh:
  > **Declare `De_min(h)` = the smallest De whose converged band spans at least 3–4 elements on the
  > finest affordable mesh — i.e. `w2/h ≥ 3`.** Run at that De. Report `q_u` at **{½De, De, 2De}**
  > *beside* the τ = 0 three-mesh band, so the reader sees the regularization's price next to the
  > disease.
- **C8**: the blended acoustic tensor must be elliptic at the declared (De, Δt) at every committed
  state; report the minimum normalised det Q and the step count it was computed at.
- **Imperfection independence (constitutive 5 / P0b-(c))**: report the converged width for at least
  two imperfection amplitudes. A width that moves more than the declared band with the imperfection
  must be disclosed as imperfection-set.
- **De collapse** at different step counts (7.1(iv)).

### 7.5 Parameters TIMs must supply (OQ2)

Target band width relative to B; ramp duration / strain rate; tolerance bands; the ultimate
criterion (OQ1). **The fork supplies `De_min(h)` from §7.4 — TIMs no longer "sets the De".** That
change is deliberate: §3.3's honest-framing test forbids choosing De to hit a target.

## 8. Phases and exit gates

| Phase | Content | Exit gate |
|---|---|---|
| **P0** ✅ **COMPLETE 2026-09-04** | ADR + numpy oracle + A0 bar | PV1–PV6 green; A0 verdicts H1–H4. `cc7c7f7a5`, `8863a468c` |
| **P0b** ✅ **COMPLETE 2026-09-05** | State-dependent `C_e` leg; dissipation gate; imperfection study; non-associated DP + blended acoustic tensor | V1–V4 (§4.5). `892c22770`, `1052ecd36`. 13 `zone_a` cases, 30.3 s |
| **GATE U** ✅ **RUN 2026-09-05 (WP-A2)** | τ = 0 three-mesh band on the SANISAND / B-bar deck; `_adr90_tau0_qu_band.md` | **INCONCLUSIVE on `q_u`** (all six legs seize at ≤ 9 % of target, no capacity); **MEASURED on width** (`w2(h/2)/w2(h) = 0.675–0.824`, non-convergent); **MEASURED on the matched-settlement load** (band contracts 7.80 → 5.28 → 2.86 %, monotone from above, **without regularization**). No "inside tolerance" verdict exists — OQ2 unsupplied. §1.2 |
| **P1 — the C++ wrapper** | `LadrunoOverstress.{h,cpp}` per §5 | **Entry conditions, all required before a line is written:** (a) GATE U fired; (b) **written consumer sign-off** on the lane (§3.1), on the `q_u` upward bias (§6 NOT-claimed 7), and on whether results may be stamped *provisional* pending OQ1; (c) D2 re-decided on the P0b evidence. **Exit:** trajectory byte gate, PV3 overstress, non-tautology guard, stage-forwarding, round-trip, MP wire, Newton-vs-ModifiedNewton bit-identity, the plane-strain interleave test, **mutation gate red**, and the **out-of-family verdict — at P1, not P3** (it is the family PR; ADR-87 D9) |
| **P2 — acceptance case** | A1 → A3 (slow tier) + A4; the De window measured by §7.4; the **non-proportionality error leg** on the SANISAND deck (OQ9); the ≥ 2-probe-direction leg (§6 NOT-claimed 8) | C4, C5, C6, C8 green; A4 two-sided. **Failure branch (strategy M9): if C5 is red — the collapse load does not converge at the declared De — the ADR does NOT iterate on τ. It falls to lane (d) and closes.** |
| **P3 — TIMs integration** | Provenance fields; one radial probe on the chosen lane | Joint sign-off; vault note records the case as *verified* |
| **P4 — optional** | Explicit/transient tier; `-implex` interplay; mesh-aware τ | Only on a named consumer |
| **WP-F — conditional** | **True Duvaut–Lions inside `LadrunoSANISAND`** (projection of the current state; relaxes σ, `α`, `z` using the `protected` base members) | Judged only under the close-out's condition below. **P0b makes this the preferred build**, not the fallback: it fixes V1 (real `C_e(p)`) and restores the convexity argument behind V2 |

### 8.1 PROPOSED close-out — the owner decides at CP2

> **P0 closes.** Regularization for localization **width** is **DISCLOSURE-ONLY**: the width does
> not converge with τ = 0 (GATE U: 0.675–0.824) and P0b shows what does converge is set by the
> imperfection field (V3), not by τ — so no wrapper delivers a τ-declared width, and the ADR
> discloses the width band rather than regularizing it.
>
> **The deliverable is already converging.** Under vault 65 D6 the campaign reports `q` at fixed
> `s/B`; that band contracts 7.80 → 5.28 → 2.86 % at τ = 0, monotone from above. A regularizer
> would bias it upward by +9.7…+17 % (R9) for no convergence gain.
>
> **The actionable engine work is an ADR-86 follow-up WP** — small, fork-owned, no new ADR: a
> **substep-count cap** on `ManzariDafalias`'s `ModifiedEuler` (`dT_min = 1e-6`, uncapped), the
> **`TanType` parser default** (elastic today — every fork SANISAND deck runs modified Newton),
> and **tolerance guidance** (`NormUnbalance` relative to weight; `NormDispIncr` is unmeetable and
> not mesh-neutral). Settle `-Presidual` while there.
>
> **Then re-run GATE U.** **WP-F is judged only if a peak becomes reachable AND the
> matched-settlement band is outside OQ2's tolerance** — which requires TIMs to supply OQ2 first.

**PROPOSED, not decided.** The owner decides at **CP2**. Nothing in §5 is authorised by it, and
33022 stays reserved either way.

## 9. Decisions

| # | Decision | Status |
|---|---|---|
| **D1** | Plasticity-type inners only; damage documented out of scope **and refused by classTag denylist** (33017, 33015, ASDConcrete3D) | DECIDED, strengthened |
| **D2** | **REOPENED 2026-09-05** — "decided at CP1 on evidence since shown incomplete". The CP1 decision (generic two-track, declared domain = proportional-and-monotonic) rested on §4.3's theorem, which P0b showed needs a **constant elastic operator** and a **holonomic** inner. SANISAND has neither. **P0b recommendation: disclosure-only (lane d) by default; WP-F conditionally; the generic wrapper not at all** (`[[_adr90_p0b_results]]` §6) | **REOPENED** |
| **D3** | **CORRECTED.** `C_e` is **deep-copied in the wrapper's `commitState()` after `inner->commitState()`** — *not* read before `setTrialStrain`. `mCe` is rewritten inside `integrate()` (`ManzariDafalias.cpp:979/987/992`), so the old rule read the previous Newton iterate's tangent and made the residual a non-function of `u`. Gated by a Newton-vs-ModifiedNewton bit-identity test | **CORRECTED** |
| **D4** | Δt from `ops_Dt`; `Δt ≤ 0 ⇒` inviscid **as an early-out**; **β and Δt LATCHED at `newStep`** and rewound on `revertToLastCommit`; a **file-scope** Δt diagnostic; the `β = 1`-with-`τ > 0` counter reported and **required to be zero** | REVISED |
| **D5** | **Lane (b) — a patterned `sp` under the 1-argument `LoadControl` — is PRIMARY** (uniform pseudo-time by `SP_Constraint.cpp:331-337`; keeps limit points, D16 guards, vault comparability). Lane (a) demoted to secondary (unpriced undamped-ringing cost at β_K = 0). The **4-argument `LoadControl` is forbidden** | **REVISED** |
| **D6** | τ as a `Parameter` | DECIDED |
| **D7** | B-bar hex for the final acceptance leg | DECIDED |
| **D8** | Flags after positionals | DECIDED |
| **D9** | Off-switch: `τ = 0` **bit-identical in stress, tangent and committed-state trajectory**, implemented as an **early-out**, with a NaN/Inf case | RESTATED |
| **D10** | **Name: `LadrunoOverstress`, not `LadrunoDuvautLions`.** The model relaxes toward a *strain-history function*, not toward a projection of the current state — the Maxwell/Krempl VBO family. The Duvaut–Lions name is reserved for WP-F, which would actually perform the projection | **CHANGED** |
| **D11** *(new)* | For a non-associated inner, an **unsymmetric solver is mandatory** in the guide and every deck | DECIDED |

## 10. Open questions

- **OQ1** [Prof. Gorini] — the ultimate criterion. Results may be stamped *provisional* pending it,
  with consumer sign-off (P1 entry condition (b)).
- **OQ2** [TIMs] — the acceptance-case numbers (§7.5), **minus the De**, which the fork now derives.
- **OQ3** [fork] — **RE-OPENED by P0b-(a).** Closed on 2026-09-04 as "two-track is adequate"; that
  answer was conditional on a hypothesis SANISAND violates.
- **OQ4** [both, P2] — lane (b) vs (a) on the real probe; whether vault 14's β_K family and the τ
  family collapse onto one De.
- **OQ5** [fork] — **CLOSED**: the metric is §7.3, with its `[h, L]` bounds now stated.
- **OQ6** [fork, P2] — interaction with SANISAND's substepping / `-Pmin` resets; verify on A2.
- **OQ7** [fork] — **REWRITTEN after WP-B (PR #785).** The `ManzariDafaliasPlaneStrain` null-ctor
  classTag anomaly (`LEDGER_quirks.md:4826`) is a **recorded owner decision NOT to fix**. Therefore
  **leg A2 must not depend on a broker / database restore of the plane-strain material**: build A2's
  plane-strain models directly in the deck and keep round-trip coverage (C7) on the 3-D view.
- **OQ8** — **CLOSED**, see D10 (and the answer changed).
- **OQ9** [fork] — **PROMOTED TO A P1 BLOCKER.** How large is the two-track approximation on
  SANISAND itself? P0b measured surrogates only: 9.2e-4…7.3e-3 (pressure-dependent 1-D, working
  De), 1.1e-4…7.5e-2 (rotating non-associated DP), 4.4e-2 (non-proportional J2). SANISAND is
  holonomic nowhere, so **no correctness statement of any size exists for the consumer's own
  material**. It is no longer acceptable to measure this after P1.
- **OQ10** [fork, P2 — new] — does the composite (wrapper + inner) satisfy a Clausius–Duhem
  inequality under any free-energy split, given V2? If not, what is the admissible envelope in De?
- **OQ11** [fork, P2 — new] — the imperfection field on SANISAND (§7.1(vii)): `-ele`-routed field,
  or per-element tags, and at what cost?

## 11. Risks

| # | Sev | Risk | Mitigation / status |
|---|---|---|---|
| **R1** | **BLOCKING (re-rated 2026-09-05)** | Two-track ≠ true DL **for any inner with a state-dependent elastic operator** — i.e. SANISAND, DP with `G(p)`, every hypoelastic model — on **every** path, monotonic and proportional included. Retired on 2026-09-04 by the theorem; **reinstated** by V1, which shows the theorem's hypothesis fails for the whole consumer class. | D2 reopened; OQ9 is a P1 blocker; §6 NOT-claimed 6. |
| **R2** | BLOCKING | Regularization inert or erratic on the campaign's lane. | D5 lane (b) primary; D4 latching; the `β = 1`-with-`τ > 0` counter. |
| **R3** | MAJOR | Band width is De-dependent and a reviewer reads "converges in h" as "mesh-independent". | §3.2 wording; C6; A0's ×12-per-×3. |
| **R4** | MAJOR | Width metric threshold-dependent / bounded. | §7.3, with the `[h, L]` bound and the `w2/h ≥ 3` floor. |
| **R5** | MAJOR | A 2-D `quad` result does not transfer to the B-bar hex. | A3 is the admissible leg. |
| **R6** | MAJOR | Static-buffer aliasing, `getCopy` routing, stage-flag forwarding, stale `C_e` — each a silent wrong answer. | §5, each with the verified source and a named gate. |
| **R7** | MINOR | Slow-tier cost. | Wall times in docstrings; nightly Zone-B. |
| **R8** | **MAJOR** | **The error on rotating stress paths is unmeasured on the real material.** A footing pushover rotates principal directions continuously — exactly what §4.3 excludes. | OQ9 is now a P1 blocker; V4's rotating-path numbers (1.1e-4 … 7.5e-2) are the only bracket. |
| **R9** | **MAJOR — new** | **The regularized `q_u` is upward-biased by +9.7 / +12.5 / +17.0 % at De = 3e-4 / 1e-3 / 3e-3** — the same order as the **11.75 %** β_K artefact the campaign deliberately removed. Trading an accidental artefact for a declared one is not progress unless it is reported. | §6 NOT-claimed 7; C5/C6 report `q_u(De)/q_u(τ=0)` per leg; §7.4 reports `q_u` at {½De, De, 2De}. |
| **R10** | **MAJOR — new** | **Incremental dissipation goes negative on unloading** (V2: −7.9e-3 of cumulative at De = 0.1), and the violated region — elastic unloading — **is the band boundary**. A total-energy check cannot see it. | §6 NOT-claimed 9; OQ10; WP-F removes it. |
| **R11** | **MAJOR — new** | **β changes per Newton iteration** under `DisplacementControl` / `ArcLength` (`DisplacementControl.cpp:346`, `ArcLength.cpp:302`), so the residual is not a function of `u`. | D4 latching, or a hard refusal of those integrators; D5 makes the uniform lane primary. |
| **R12** | **MINOR — new** | **One global τ distorts the fitted ellipsoid** under a direction-dependent rate field. | §6 NOT-claimed 8; the ≥ 2-probe-direction WP-D leg. |
| **R13** | **BLOCKING — new (GATE U)** | **The campaign is blocked by the SANISAND integrator, not by ill-posedness — and this ADR could absorb effort that belongs there.** GATE U's six legs seized at 9 % of target settlement with **0/80** subdivisions used and steps **6400–25000×** above the floor, one of them burning **2056 s (34.3 min) — 59 % of its budget — on a single converged step**; the cost is `dT_min = 1e-6` with no substep cap (`ManzariDafalias.cpp:1380`), compounded by an elastic-tangent parser default (`LadrunoSANISAND.cpp:117`) and an unmeetable convergence test. A wrapper cannot fix any of those, and P2's A2/A3 legs would seize the same way. | §6 NOT-claimed 10; the §8 **PROPOSED close-out** routes the actionable work to an ADR-86 follow-up WP and re-runs GATE U after it. |
| **R14** | **MAJOR — new (GATE U)** | **The deliverable may already be converged without any regularization.** Under vault 65 D6 the reported quantity is `q` at fixed `s/B`, and its τ = 0 three-mesh band contracts 7.80 → 5.28 → 2.86 % monotonically from above. Shipping a regularizer for a quantity that is already converging — and one that biases it upward by +9.7…+17 % (R9) — would make the campaign's numbers worse, not better. | §1.2.1; C5 re-based; the close-out makes width-regularization **disclosure-only**. |
| **R15** | **MAJOR — new (GATE U)** | **Every existing fork SANISAND result was produced with the ELASTIC tangent.** The parser default `TanType 0` makes `algorithm Newton` behave as modified Newton; it is invisible on the zero-free-DOF material-point decks that constitute the fork's whole SANISAND test surface, and costs ~7× on a BVP. Convergence claims on any prior SANISAND BVP deck are suspect. | §1.2.2; owned by the ADR-86 follow-up WP, not by this ADR. |
| **R16** | **MINOR — new (GATE U)** | **`-Presidual 0.0` (the ADR-86 default) clamps a free-surface Gauss point** on the dense coarse leg at `s/B ≈ 0.0153`, past which the answer is the clamp's, not the model's. A faster deck hits it next. | A declared `-Presidual`, a surcharge, or an accepted-and-disclosed clamp — a **choice**, made before any deeper push. |

## 12. Ledger obligations

**Discharged in PR #783:** the ADR; the `LEDGER_implementations` row (33022 **RESERVED, not yet
built**); the `README` index line; dated CORRECTION callouts in the planning brief (§2 F2, §3 lane,
§6 P0b, §7 D2, §8 OQ7); and the `LEDGER_quirks` entries for the two tautological A0 gates, the
`H_res` default, the post-peak Newton predictor, **the `E(σ)` hypothesis, `w2` saturation at L,
β-per-Newton-iteration under `DisplacementControl`, `revertToLastCommit` setting `dT = 0` ⇒ β = 1,
and the plane-strain static `mTangent_init`**.

**Owed by P1 (not before):** `SRC/classTags.h` 33022 with the Ladruno comment and the per-registry
non-collision note (**PATTERN 33022 is live at `Domain.cpp:2274`**); the verification manifest
(`ledger/ladruno_overstress.yml`, ≥ 1 non-self-referential oracle — §4.4 scoped to perfect
plasticity); `-DLADRUNO_MUTATE_OVERSTRESS`; the guide stub, **including the prohibition on tuning τ
against a target**; both interpreters; the banner line **on ship**; and two further quirks:
`FluidSolidPorousMaterial::getCopy(const char*)` ignoring its `code` argument, and
`ManzariDafalias::setParameter` gating on `atoi(argv[1]) == getTag()`.

## 13. Implementation log

| Date | Event |
|---|---|
| 2026-09-04 | **Planning brief** from three read-only source sweeps. Number 90 allocated; 33022 reserved. `3feec0fc9`. |
| 2026-09-04 | **WP-A / A0.** Oracle + bar; PV1–PV6; the TT ≡ TDL theorem found and proved; H1–H4. `cc7c7f7a5`, `8863a468c`. |
| 2026-09-04 | **CP1 — owner decision.** D2 approved as recommended (generic two-track; declared domain proportional-and-monotonic; OQ9 added to WP-D; WP-F parked with a trigger). ADR written against it. `819b69022`. |
| 2026-09-04 | WP-B (PR #785) landed the SANISAND plane-strain prerequisites and established that the plane-strain classTag anomaly is a recorded owner decision **not** to fix — rewriting OQ7. |
| **2026-09-05** | **Three-lens adversarial pass** (strategy / numerics / constitutive). Findings: the need is asserted not measured; the theorem needs a constant elastic operator and a holonomic inner; `q_u` is upward-biased; D3 reads a stale tangent; β changes per Newton iteration; the model is an overstress model, not Duvaut–Lions. |
| **2026-09-05** | **Owner decision: RE-SCOPE (option A). D2 REOPENED**, "decided at CP1 on evidence since shown incomplete". No C++ opens. |
| **2026-09-05** | **P0b executed** — V1 state-dependent `C_e`; V2 dissipation; V3 imperfection; V4 blended acoustic tensor. `892c22770`, `1052ecd36`. Recommendation: disclosure-only default, WP-F conditional, generic wrapper not recommended. |
| **2026-09-05** | **GATE U executed (WP-A2, branch `wp/90a2-tau0-qu-band`, `_adr90_tau0_qu_band.md`)** on engine `548fe911427e90a2edfead05cb3672a738d25b6d`, 3 meshes x 2 densities on `LadrunoBrick -formulation bbar` under lane (b). All six legs **seized** at <= 9 % of target settlement with **0/80** subdivisions used and steps 6400x/12800x/25000x above the floor (h0.25/h0.5/h1.0) — the binding resource is SANISAND's uncapped `ModifiedEuler` substepping, whose longest gap between consecutive converged steps is 2056 s = 34.3 min (59 % of that leg's entire budget in one step). So **`q_u` is unreachable, not merely uncertain**. Width **does not converge** (0.675-0.824); yielding **volume** is near-objective (3.95 %); the **matched-settlement load band contracts** 7.80 -> 5.28 -> 2.86 % at tau = 0. Two fork-wide deck defects fixed first: the `TanType` elastic parser default and an unmeetable, mesh-non-neutral `NormDispIncr`. |
| **2026-09-05** | **Close-out PROPOSED (§8.1)** on the GATE U evidence — P0 closes; width regularization disclosure-only; the actionable work is an ADR-86 integrator follow-up WP (substep cap, `TanType` default, tolerance guidance, `-Presidual`); GATE U re-run after it; WP-F judged only if a peak becomes reachable *and* the matched-settlement band is outside OQ2. **Owner decides at CP2.** |
| **2026-09-05** | **This ADR revised against all three critics** — status, §1.2 GATE U, §3 lane inversion, §4.0 renaming, §4.3 hypotheses, §4.5 P0b, §5 D3/M5/M6 corrections, §6 the bias and the dissipation, §7 the De rule and the two-sided A4, §8 P1 entry conditions and the P2 failure branch, §9 D2 reopened, §11 R1 reinstated + R9–R12. |

## References

- **de Borst, R., Sluys, L. J., Mühlhaus, H.-B. & Pamin, J. (1993).** "Fundamental issues in finite
  element analyses of localization of deformation." *Engineering Computations* **10**, 99–121.
  The comparison of regularization strategies, the canonical plane-strain biaxial acceptance deck,
  and the internal length `ℓ = 2mc_e/E`. **The caveat this ADR is built around**, paraphrased: rate
  dependence works for both failure mechanisms, but its applicability is limited to transient
  loading and the regularizing effect falls away rapidly for slow loading or near the
  rate-independent limit. A0 §5.2 and P0b-(c) are that sentence in numbers.
- **Simo, J. C. & Hughes, T. J. R. (1998).** *Computational Inelasticity.* Springer. §2.7 —
  viscoplastic regularization; the Duvaut–Lions backward-Euler closed form and the relaxation of the
  internal variables. Box 3.2 — the J2 algorithmic tangent used by the oracle.
- **Simo, J. C., Kennedy, J. G. & Govindjee, S. (1988).** "Non-smooth multisurface plasticity and
  viscoplasticity." *IJNME* **26**, 2161–2185. The closest-point-projection formulation of
  Duvaut–Lions — the model this ADR is **not** implementing (§4.0).
- **Krempl, E. (1987, and the VBO literature).** Viscoplasticity based on overstress — the family
  the two-track update actually belongs to, and the reason for the rename (D10).
- **Perzyna, P. (1966).** "Fundamental problems in viscoplasticity." *Advances in Applied Mechanics*
  **9**, 243–377. The alternative overstress family; a different class, not a flag.
- **Rudnicki, J. W. & Rice, J. R. (1975).** "Conditions for the localization of deformation in
  pressure-sensitive dilatant materials." *JMPS* **23**, 371–394. The acoustic-tensor criterion
  P0b-(d) evaluates, and the reason a non-associated material localizes while still hardening.
- **Rice, J. R. (1976).** "The localization of plastic deformation." *Proc. 14th IUTAM*, 207–220.
  The band-boundary-as-elastic-unloading picture that makes §4.3 hypothesis 3 fatal rather than
  limiting.
- **Needleman, A. (1988).** "Material rate dependence and mesh sensitivity in localization
  problems." *CMAME* **67**, 69–85. Rate dependence restores well-posedness and sets a length
  through the **wave speed** — and, read carefully, the origin of the quasi-static caveat that §3
  turns into a lane decision.
