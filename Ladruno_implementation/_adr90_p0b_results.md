---
title: "ADR 90 leg P0b — results: the four measurements the 3-lens adversarial pass showed were missing"
project: Ladruno
type: measurement note (phase P0b, WP-A)
status: "COMPLETE 2026-09-05 — D2 is REOPENED and these numbers decide it"
owner: nmora
related:
  - "[[_adr90_a0_results]]"
  - "[[_adr90_regularization_planning_brief]]"
  - "[[90_ladruno_viscoplastic_regularization_adr]]"
tags: [adr, regularization, viscoplastic, duvaut-lions, localization, oracle, measurement, p0b]
updated: 2026-09-05
---

# ADR 90 leg P0b — results

> [!abstract] One line
> A0's theorem had an **unwritten hypothesis** — a constant elastic operator — and the three
> exposures the critics named are all real: the generic wrapper is inexact on **every** path over
> a pressure-dependent inner, its dissipation increment goes **negative** on unloading, and the
> "regularized" band width is set **jointly by De and the imperfection field**, not by τ. The one
> thing that does work is the one A0 could not see: the **blended acoustic tensor restores
> ellipticity** for the consumer's material class — but on a criterion in β = Δt/(τ+Δt), not De.

## 0. Provenance

| | |
|---|---|
| Oracle | `tests/_testbed/duvaut_lions_ref.py` §7–§10 (numpy only; imports no OpenSees) |
| Driver | `python3.12 tests/_testbed/run_a0_sweep.py --p0b` — regenerates every table here (~60 s) |
| Gate | `tests/test_duvaut_lions_oracle.py` — 13 `zone_a` cases, **30.3 s measured** (budget 40 s) |
| Git | `892c22770` on `wp/90a-adr-oracle` |
| Date | 2026-09-05 |
| Engine | numpy 2.4.6, python3.12; **no C++ build exists for this work package** |

## 1. Leg (a) — a STATE-DEPENDENT elastic operator breaks the identity on every path

### 1.1 The mechanism, stated before the numbers

A0's theorem rests on `ψ = σ + Eα` advancing by exactly `E·Δε` every step. That requires **E to be
the same number on both tracks**. A generic `NDMaterial` wrapper has exactly one elastic operator
available to it — `inner.getInitialTangent()`, i.e. `E` evaluated at the **inner's** committed
stress `σ̄_n`. But the stress it must predict from is its **own** `σ_vp,n`, and relaxation is
precisely the statement that `σ_vp ≠ σ̄`. So the wrapper integrates

```
    sigma_tr = sigma_vp,n + E(sigma_bar_n) * de          # what the seam can deliver
```

while true Duvaut–Lions integrates `σ_vp,n + E(σ_vp,n)·Δε`. The conserved quantity is destroyed at
first yield and the error accumulates **on every path — monotonic, proportional, single-surface**.
SANISAND (`G(p)`, `K(p)`), Drucker–Prager with `G(p)`, and every hypoelastic inner are in this
class. For the named consumer this is not a corner case; it is the normal case.

### 1.2 Measured — 1-D, isotropic hardening `H/E = +0.1`, monotonic proportional ramp, `E(σ) = E₀(1 + β_E|σ|/σ_Y)`

Max relative `|σ_TT − σ_TDL|` along the path (`nsteps = 2000`):

| De | constant E | β_E = 0.3 | β_E = 0.6 | β_E = 1.2 |
|---|---|---|---|---|
| 1e-4 | 4.3e-15 | 8.85e-08 | 2.39e-07 | — |
| 3e-4 | 1.4e-14 | 8.38e-07 | 2.24e-06 | — |
| 1e-3 | 1.9e-14 | 9.45e-06 | 2.52e-05 | — |
| 3e-3 | 3.0e-14 | 8.49e-05 | 2.26e-04 | — |
| **1e-2** | 1.5e-14 | **9.19e-04** | **2.44e-03** | **7.32e-03** |
| 3e-2 | 2.6e-15 | 7.69e-03 | 2.04e-02 | — |
| 1e-1 | 2.0e-14 | 6.80e-02 | 1.79e-01 | — |
| 3e-1 | 8.4e-15 | 2.79e-01 | **5.99e-01** | 9.44e-01 |

Summary (max over the whole ladder / over the **working window** De ∈ [1e-4, 1e-2]):

| β_E | all De | working window | working window, the **non-buildable** repair |
|---|---|---|---|
| 0.3 | 2.79e-01 | 9.19e-04 | 2.99e-05 |
| 0.6 | 5.99e-01 | 2.44e-03 | 6.27e-05 |
| 1.2 | 9.44e-01 | 7.32e-03 | 1.33e-04 |

This **independently reproduces the constitutive critic**, who measured 3.0e-3 / 2.5e-2 / 2.1e-1 /
6.3e-1 at De = 0.01 / 0.03 / 0.1 / 0.3; this oracle gets 2.44e-3 / 2.04e-2 / 1.79e-1 / 5.99e-1 on
its own deck at β_E = 0.6. Same order, same monotone growth, different code.

### 1.3 The repair, and why it is not available

Variant `tt_vp` is the obvious fix: evaluate the predictor modulus at the **wrapper's own**
committed stress. It buys **1–2 orders** (2.44e-3 → 6.27e-5 in the window; 5.99e-1 → 1.43e-2 at
De = 0.3) but does **not** close the gap, because TDL also relaxes the internal variable and the
wrapper cannot.

**And it is not implementable at the `NDMaterial` seam.** `getInitialTangent()` returns *one
sampled matrix*, not the function `E(·)`. To evaluate the modulus at a stress the inner has never
been at, the wrapper would have to either (i) call `inner->setTrialStrain(C_e⁻¹σ_vp)` and read the
tangent — which corrupts the inner's trial state and is a different material anyway, or (ii) know
the inner's elastic law, which is exactly what "generic" forbids. **No wrapper-accessible choice
of a single scalar/matrix modulus closes the gap.**

> **V1 — generic TT over a pressure-dependent inner: error 9.2e-4 … 7.3e-3 at the working De
> (β_E = 0.3…1.2), rising to 0.28…0.94 at De = 0.3. The identity is exact (3.0e-14) *only* for a
> constant elastic operator. ⇒ NOT a correctness claim — a declared approximation whose size on
> SANISAND's real `G(p)` and pressure range is UNMEASURED.**

## 2. Leg (b) — the dissipation gate

`D_w = σ_vp · Δε^vp_w`, `Δε^vp_w = C_e⁻¹(σ_vp − σ̄)Δt/τ = β C_e⁻¹(σ_tr − σ̄)` (the second form is
used, being finite at τ = 0 and algebraically identical). In **true** Duvaut–Lions `σ̄` is the
closest-point projection **of `σ_vp`**, so `σ_vp − σ̄` lies in the normal cone and `D_w ≥ 0`
follows from convexity. The two-track wrapper projects a **different** trial, so nothing enforces
that geometry.

Negatives below round-off (|worst| < 1e-9 × cumulative) are reported but not counted.

| path | H | De | `D_w` cumulative | worst step | n_neg | sum(neg)/cum | significant | inner's own |
|---|---|---|---|---|---|---|---|---|
| J2 non-proportional | +2000 | 1e-3 | 1.887e-01 | −1.53e-18 | 11 | −3.0e-17 | no | 1.812e-01 |
| J2 non-proportional | +2000 | 1e-1 | 3.070e-01 | −2.48e-22 | 8 | −2.1e-21 | no | 1.812e-01 |
| J2 non-proportional | −1000 | 1e-3 | 1.384e-01 | −1.53e-18 | 11 | −4.0e-17 | no | 1.284e-01 |
| J2 swipe (rotate at fixed \|ε\|) | ±  | all | 2.1e-01 … 3.6e-01 | **> 0** | 0 | 0 | no | 1.3e-01 … 2.0e-01 |
| **1-D cyclic load/unload/reload** | +2000 | 1e-2 | 1.537e-01 | −7.34e-11 | 189 | −1.62e-08 | no | 1.446e-01 |
| **1-D cyclic** | **+2000** | **1e-1** | 1.917e-01 | **−1.95e-05** | 156 | **−7.90e-03** | **YES** | 1.446e-01 |
| **1-D cyclic** | **−1000** | **1e-2** | 1.521e-01 | **−3.91e-08** | 143 | **−1.43e-05** | **YES** | 1.347e-01 |
| **1-D cyclic** | **−1000** | **1e-1** | 2.085e-01 | **−2.26e-05** | 119 | **−6.06e-03** | **YES** | 1.347e-01 |

Every run still dissipates net-positive overall, and the inner's own dissipation is positive
throughout — so the violation cannot be seen in a total-energy check; it is per-increment.

> **V2 — dissipation: VIOLATED on load/unload/reload. Worst single step −2.26e-5 (−1.1e-4 of that
> run's cumulative); worst cumulative-negative fraction **−7.9e-3** at De = 0.1, growing
> monotonically with De (round-off at De = 1e-3, −1.4e-5 at 1e-2, −7.9e-3 at 1e-1). Non-proportional
> J2 and the swipe surrogate are clean to round-off.**

The exposure is not academic: **the boundary of a localization band is an elastic-unloading zone**
(Rice), which is the region that sets the band's width — so the violation lives exactly where the
answer is decided.

## 3. Leg (c) — the converged width is an IMPERFECTION property

Run at De = 3e-4 — the one De at which A0 reported **both** `w2` and `w3` converging — varying the
imperfection amplitude and the zone length. `w2/h` is the band-resolution floor; `w2/L` matters
because **`w2 ∈ [h, L]` by construction** (a uniform profile returns exactly L), so a `w2` near L
means *no band*, not *a wide band*.

| amp | zone/L | N | h | w2 | w3 | w2/h | w2/L | P_peak |
|---|---|---|---|---|---|---|---|---|
| 0.05 | 0.05 | 80 / 160 | 1.250 / 0.625 | 25.098 / 25.479 | 2.500 / **1.875** | 20.1 / 40.8 | 0.251 / 0.255 | 20.082 / 20.083 |
| 0.05 | 0.10 | 80 / 160 | | 12.949 / 13.027 | 2.500 / 2.500 | 10.4 / 20.8 | 0.129 / 0.130 | 20.049 / 20.049 |
| 0.05 | 0.20 | 80 / 160 | | 6.500 / 6.426 | 2.500 / 2.500 | 5.2 / 10.3 | 0.065 / 0.064 | 19.896 / 19.895 |
| 0.10 | 0.05 | 80 / 160 | | 11.826 / 12.226 | 2.500 / **1.875** | 9.5 / 19.6 | 0.118 / 0.122 | 20.046 / 20.048 |
| **0.10** | **0.10** | 80 / 160 | | **3.629 / 3.527** | 2.500 / 2.500 | 2.9 / 5.6 | 0.036 / 0.035 | 19.752 / 19.754 |
| 0.10 | 0.20 | 80 / 160 | | 5.269 / 5.189 | 2.500 / 2.500 | 4.2 / 8.3 | 0.053 / 0.052 | 19.119 / 19.121 |
| 0.20 | 0.05 | 80 / 160 | | 2.607 / 2.176 | 2.500 / **1.875** | **2.1** / 3.5 | 0.026 / 0.022 | 19.343 / 19.392 |
| 0.20 | 0.10 | 80 / 160 | | 3.140 / 3.042 | 2.500 / 2.500 | 2.5 / 4.9 | 0.031 / 0.030 | 18.169 / 18.183 |
| 0.20 | 0.20 | 80 / 160 | | 4.333 / 4.260 | 2.500 / 2.500 | 3.5 / 6.8 | 0.043 / 0.043 | 17.403 / 17.402 |

The row in bold is A0's original parameter point.

Three findings:

1. **`w2` is mesh-converged everywhere** (N = 80 vs 160 agree to 1–3 %) — so this is a property of
   the converged solution, not a discretization artefact.
2. **The converged `w2` moves ×4.28 with imperfection amplitude** (zone 10 %, finest mesh:
   13.03 → 3.04) **and ×3.47 with zone length** (amp 10 %: 12.23 → 5.19). At fixed De. The
   imperfection is a **modelling choice with no declared status**; τ is a declared numerical
   parameter. The width is a function of both, and the imperfection is not the smaller term.
3. **`w3` (FWHM) does not converge** outside A0's original parameter point: it drops 2.500 → 1.875
   (tracking `h`) in three of the nine configurations. A0's report that "both `w2` and `w3`
   converged at De = 3e-4" was **parameter-lucky**. `w2` also hits the resolution floor
   (`w2/h = 2.1` at amp 0.20 / zone 0.05, N = 80) where the "band" is two elements wide.

> **V3 — the width is IMPERFECTION-SET, not τ-set. At fixed De the converged band width varies by
> ×4.3 with imperfection amplitude and ×3.5 with zone length. Viscosity stops the band from
> collapsing onto one element; *what it collapses to* is chosen by the imperfection field.**

## 4. Leg (d) — non-associated Drucker–Prager and the BLENDED ACOUSTIC TENSOR

New plane-strain point model (§10 of the oracle): `f = √J₂ + ρp − c(α)`, `g = √J₂ + ρ̄p` with
**ρ̄ = 0.15 < ρ = 0.40** (non-associated), cohesion softening `c = max(20 − 200α, 5)`,
E = 20000, ν = 0.25. Algorithmic tangent by central differencing the return map (exact to ~1e-8,
and it cannot get the unsymmetric non-associated structure wrong by hand). Path: isotropic
consolidation to p = −100, then a plane-strain deviatoric push, 800 steps, 725 of them plastic,
no apex returns.

### 4.1 TT vs TDL over a non-associated model

| path | De | β | max rel `|σ_TT − σ_TDL|` |
|---|---|---|---|
| proportional | 0 … 0.1 | 1.000 … 0.012 | **2.0e-14** (identity) |
| rotating (principal directions turn 63°) | 1e-4 | 0.926 | 1.13e-04 |
| rotating | 3e-4 | 0.807 | 3.38e-04 |
| rotating | 1e-3 | 0.556 | 1.12e-03 |
| rotating | 1e-2 | 0.111 | 1.08e-02 |
| rotating | 1e-1 | 0.012 | 7.49e-02 |

**Non-associativity alone does not break the theorem** — a proportional path over a non-associated
model keeps the identity to 2e-14, because the flow direction is fixed. **Path rotation does**, and
the error grows monotonically with De.

### 4.2 Ellipticity — the only measurement that speaks to well-posedness for this material class

Minimum over in-plane band normals of `det[n·C·n]` (2×2 in-plane block), normalized by the elastic
value. `≤ 0` ⇒ loss of ellipticity (Rudnicki & Rice 1975).

**Inviscid**: min normalized `det Q = −0.0145`, first lost at step 75 of 800, at a band normal of
**51.5°**. So the non-associated softening tangent does what the theory says it does.

**Blended** `(1−β)C_e + β C^ep` at the worst inviscid state:

| β | min det Q / det Q_e | θ [deg] |
|---|---|---|
| 1.00 | **−0.01451** | 51.5 |
| 0.99 | **−0.00442** | 128.5 |
| 0.95 | +0.03594 | 128.5 |
| 0.90 | +0.08642 | 51.5 |
| 0.80 | +0.18746 | 51.5 |
| 0.50 | +0.49129 | 51.5 |
| 0.10 | +0.89803 | 51.5 |
| 0.00 | +1.00000 | — |

**β_crit = 0.9857** — the blend is elliptic for β below that. Since `β = 1/(1 + nsteps·De)`, the
equivalent De is **`De > 1.45e-2 / nsteps`**:

| nsteps | 250 | 1000 | 2000 | 8000 |
|---|---|---|---|---|
| De required | 5.81e-05 | 1.45e-05 | 7.26e-06 | 1.81e-06 |

Every De in A0's working window clears it with large margin (at nsteps = 800: De = 1e-4 → β = 0.926,
det = +0.060; De = 3e-4 → β = 0.807, det = +0.181; De = 1e-3 → β = 0.556, det = +0.435).

> **V4 — the blended acoustic tensor IS elliptic above De ≈ 1.45e-2/nsteps (β < 0.9857). The
> mechanism does restore well-posedness at the material point for the consumer's material class —
> but the criterion is on β = Δt/(τ+Δt), i.e. on the STEP COUNT as much as on τ. The same De is
> elliptic or not depending on how finely the analyst steps.**

This is the sharpest form of the "the (h → 0, De → 0) limits do not commute" objection: ellipticity
is bought by making Δt small relative to τ, which is a property of the time discretization, not of
the material.

## 5. Verdict lines

- **V1** — *generic TT over a pressure-dependent inner:* error **9.2e-4 … 7.3e-3** at the working
  De (β_E = 0.3…1.2), **0.28…0.94** at De = 0.3; exact (3e-14) **only** for a constant elastic
  operator; no wrapper-accessible repair closes it. **NOT acceptable as a correctness claim.**
  Acceptable only as a declared approximation — and its size on SANISAND is unmeasured.
- **V2** — *dissipation:* **violated** by up to **−7.9e-3 of the cumulative dissipated work**
  (worst step −2.26e-5) on **load/unload/reload** paths, growing with De; clean on non-proportional
  J2 and on the swipe surrogate. The violated region — elastic unloading — is the band boundary.
- **V3** — *width:* **imperfection-set**, not τ-set (×4.3 with amplitude, ×3.5 with zone length at
  fixed De); and `w3` does not converge outside A0's original parameter point.
- **V4** — *blended acoustic tensor:* **elliptic above De ≈ 1.45e-2/nsteps** (β < 0.9857); the
  inviscid tangent loses ellipticity at −0.0145 / 51.5°.

## 6. D2 recommendation

The rule set at re-scope was: *generic TT wrapper only if V1 and V2 pass.* **They do not.** V1
removes the exactness claim on the consumer's entire material class, and V2 is a genuine
per-increment second-law violation located precisely at the band boundary.

That leaves *true DL inside `LadrunoSANISAND` (WP-F)* and *disclosure-only (lane d)*. **V3 decides
between them, and it decides against building anything yet:**

> **RECOMMENDATION — disclosure-only (lane d) is the default; WP-F is the conditional upgrade;
> the generic TT wrapper is not recommended at all.**
>
> 1. **Do not build the generic wrapper.** V1 and V2 fail on the named consumer's material class,
>    and the ADR's justification for going generic — A0's theorem — has been shown to hold only
>    under a hypothesis (constant `C_e`) that SANISAND violates by construction.
> 2. **Do not build anything until the need is demonstrated.** V3 is variant-independent: true DL
>    would fix V1 and V2, but it would **not** make the width a τ property. Whatever is built
>    delivers *a converged solution at a declared (De, Δt, imperfection field)*, not *a
>    material-length-regularized band*. The deliverable the campaign was promised is not
>    achievable by rate regularization in quasi-statics — which is what de Borst et al. 1993 say,
>    now with numbers. So the strategy critic's **un-descope gate** must fire first: measure the
>    τ = 0 three-mesh `q_u` band on the SANISAND / B-bar deck; if it is inside the campaign
>    tolerance, this is width *research*, not a P1 deliverable.
> 3. **If the gate fires, build WP-F, not the wrapper.** True DL inside `LadrunoSANISAND` can use
>    the real `C_e(p)` and re-seed `α`, `z`, which fixes V1 and restores the convexity argument
>    behind V2. The generic wrapper's advantage was that it was free; V1 and V2 are the price, and
>    it is higher than the material-specific implementation.
> 4. **Whatever is built, V4 is the gate that means something at the material point**: the blended
>    acoustic tensor must be elliptic at the declared (De, Δt) on A1/A2/A3, and the number must be
>    reported, because it depends on the step count.

## 7. What P0b does *not* settle

- The `E(σ)` leg is a 1-D caricature with a linear pressure dependence; SANISAND's `G(p) ∝ p^0.5`
  over a footing's pressure range is a different and probably larger perturbation. **Measuring it
  on the real material is still owed** (ADR 90 OQ9).
- The dissipation gate measures the wrapper's own increment. A full Clausius–Duhem statement for
  the composite (wrapper + inner) with a free-energy split is not attempted.
- The DP leg is a point model. It says the *tangent* is elliptic; it does not run a BVP, so it does
  not show that the *band* on a 2-D mesh converges for this material.
- Nothing here has been run through the C++ engine. Every number is numpy.
