---
title: "LadrunoJ2Finite — Finite-Strain-Native Combined-Hardening J2 (co-rotating backstress)"
project: Ladruno
type: reference / user guide
status: shipped (v1 native material #127 · IMPL-EX #134 · analytic channel B #139)
classTag: 33012
material: nDMaterial LadrunoJ2Finite
related:
  - "[[LadrunoJ2_guide]]"
  - "[[finite_strain_trifecta_guide]]"
  - "[[16_finite_native_j2_adr]]"
  - "[[LadrunoBrick_reference]]"
  - "[[09_finite_strain_material_wrapper]]"
  - "[[Ladruno_materials_guide]]"
  - "[[LadrunoUniaxialJ2_guide]]"
tags:
  - material
  - plasticity
  - j2
  - vonMises
  - finite-strain
  - kinematic-hardening
  - objectivity
  - reference
---

# LadrunoJ2Finite — Finite-Strain-Native Combined-Hardening J2

> [!abstract] One-line summary
> `LadrunoJ2Finite` is a **finite-strain-native** combined-hardening von Mises
> (J2) `nDMaterial` for OpenSees. It does what the wrapper path
> `LogStrain → LadrunoJ2` **cannot**: it makes **kinematic (Chaboche) hardening
> frame-indifferent under large rotation**, by owning the backstress as a
> **spatial** tensor and **co-rotating it each step** with the incremental
> material rotation. It crosses the dSNPO **§14.11** objectivity boundary that v1
> pins as a strict `xfail`. It is a [[finite_strain_trifecta_guide|FiniteStrainNDMaterial]]
> subclass driven by `setTrialF(F)`, and its **only consumer** is
> `LadrunoBrick -geom finite` (3D-only).

This document is the **descriptive reference**: why the material exists (the
§14.11 boundary the wrapper cannot cross), the one genuinely new algorithmic step
(the co-rotating backstress), the two-channel consistent tangent, the IMPL-EX
mode, the OpenSees wiring, and when to reach for this material instead of the
simpler wrapper path. For the chronological design/decision log, the numpy-oracle
evidence, and the full follow-up backlog see [[16_finite_native_j2_adr]]. For the
small-strain constitutive law it inherits (yield surface, Voce isotropic,
Chaboche kinematic, the scalar return map) see [[LadrunoJ2_guide]] — this guide
does **not** repeat it.

---

## Contents

- [[#1 Intended use cases|1. Intended use cases]]
- [[#2 Why it exists — the §14.11 boundary|2. Why it exists — the §14.11 boundary]]
- [[#3 The fix — a co-rotating spatial backstress|3. The fix — a co-rotating spatial backstress]]
- [[#4 The verified algorithm|4. The verified algorithm]]
- [[#5 The consistent tangent — two channels|5. The consistent tangent — two channels]]
- [[#6 IMPL-EX (`-implex`)|6. IMPL-EX]]
- [[#7 OpenSees implementation|7. OpenSees implementation]]
- [[#8 Command syntax and usage|8. Command syntax]]
- [[#9 LadrunoJ2Finite vs the wrapper path|9. LadrunoJ2Finite vs the wrapper path]]
- [[#10 Verification and validation|10. Verification]]
- [[#11 Limitations and scope|11. Limitations and scope]]
- [[#12 References|12. References]]

---

## 1. Intended use cases

`LadrunoJ2Finite` exists for **one specific corner of the modeling space**:
combined (kinematic) hardening metal plasticity where **both** the strain **and**
the rotation are large and the cyclic realism of the backstress must remain
physically correct as the material rotates.

| Use case | Why LadrunoJ2Finite |
|---|---|
| **Finite cyclic plasticity with large rotation** (buckling-brace-like loops, large-amplitude low-cycle fatigue with member rotation) | The backstress (centre of the yield surface) co-rotates with the material, so the **Bauschinger effect and ratcheting stay frame-indifferent** as the element spins — the property the wrapper path loses. |
| **Buckling-restrained braces / fuses / shear links taken into the large-rotation regime** | Chaboche kinematic hardening drives the cyclic loop; objectivity under large rotation keeps that loop physically meaningful through the buckled / rotated geometry. |
| **Metal forming / upsetting with cyclic or reversed straining** | Combined hardening at genuinely large strain *and* large rotation, in one native return map. |
| **A reference for the §14.11 boundary** | It is the fork's demonstrator that the wrapper's non-objectivity is a *framework* limit, not a bug — and the constructive fix. |

> [!tip] When **not** to reach for LadrunoJ2Finite
> - **Isotropic hardening at finite strain** (`-kin 0`) — the wrapper path
>   `LogStrain → LadrunoJ2` is **already exact and objective** there, and simpler
>   (one fewer concept). Use it. See [[#9 LadrunoJ2Finite vs the wrapper path|§9]].
> - **Small-strain or only-large-rotation** combined hardening — use plain
>   `LadrunoJ2` (`-geom linear` or `-geom corot`); `corot` is objective for
>   kinematic hardening because the de-rotated frame is a reference frame. See
>   [[finite_strain_trifecta_guide]] §7.3.
> - **2D finite plasticity** — there is no 2D finite element consumer, so the
>   dimensional finite views are orphaned and not built ([[#11 Limitations and scope|§11]]).
> - **Pressure-sensitive / non-metal** materials — the whole J2 family is
>   deviatoric-only (use Drucker–Prager / Mohr–Coulomb / ASDConcrete / ASDPlastic).

---

## 2. Why it exists — the §14.11 boundary

### 2.1 The wrapper path is exact for the spine but not the backstress

The shipped finite-strain route (PR #97, the trifecta) is
`LadrunoBrick -geom finite → nDMaterial LogStrain → nDMaterial LadrunoJ2`: it lifts
the **unchanged** small-strain material to finite strain by the Hencky log-strain
technique. The log-strain (MATISU) algorithm co-rotates the **elastic** state
automatically,

$$\mathbf b^e_{\text{tr}} = \mathbf f_\Delta\,\mathbf b^e_n\,\mathbf f_\Delta^{\mathsf T},$$

so the deviatoric stress $\mathbf s$ emerges in the **current** frame. That is
enough to make the **isotropic spine exact**:

- **Isotropic** yield uses only $\lVert\mathbf s\rVert$ and the scalar
  $\bar\varepsilon^p$ — both rotation-invariant ($\lVert\mathbf Q\mathbf s\mathbf Q^{\mathsf T}\rVert = \lVert\mathbf s\rVert$) — so a superposed rotation rotates
  the Cauchy stress rigidly (verified to $\sim10^{-9}$).

### 2.2 Why kinematic hardening breaks

The kinematic part uses the **relative (shifted) stress** $\boldsymbol\xi = \mathbf s - \boldsymbol\alpha$, and the backstress $\boldsymbol\alpha$ is a
**second-order tensor** — the centre of the yield surface, with its own
orientation. In the wrapper path $\boldsymbol\alpha$ lives **inside the inner
material in a fixed frame**, and the wrapper has **no channel** to rotate it:
`LadrunoJ2` exposes $\boldsymbol\alpha$ read-only via `getResponse`, with no
setter. So $\lVert\mathbf s - \boldsymbol\alpha\rVert$ subtracts two tensors in
**different frames** ⇒ the yield check is wrong, the flow direction is wrong, and
the response is **not** frame-indifferent under large rotation.

> [!important] This is a framework limit, not a wiring bug
> The identical failure appears in the direct dSNPO Box-14.4 chain (verified), and
> a **33-agent review of PR #97** independently confirmed it is **fundamental to
> the "unchanged inner" v1 architecture** — `LogStrain`'s contract is *"wrap an
> UNCHANGED small-strain `NDMaterial` as a black box,"* which structurally cannot
> reach the inner's $\boldsymbol\alpha$. This is exactly the
> kinematic-hardening-at-finite-strain case de Souza Neto, Perić & Owen defer to
> **§14.11**. v1 pins it as a strict `xfail`
> (`test_ladrunoJ2_finite.py::test_finite_combined_hardening_large_rotation_objectivity_is_v2`).

> [!note] The 1D twin is immune
> For the 1D `LadrunoUniaxialJ2` the backstress is a **scalar** $X$ along the bar
> axis — there is no frame to rotate, so the uniaxial twin is objective by
> construction. This is the precise meaning of *"kernel-sharing is
> nD↔finite-strain only, not 1D"* — there is **no** `LadrunoJ2Finite` analogue for
> the uniaxial material.

---

## 3. The fix — a co-rotating spatial backstress

`LadrunoJ2Finite` is a `FiniteStrainNDMaterial` subclass that **owns all the
finite-strain state itself** — the elastic left Cauchy–Green tensor
$\mathbf b^e_n$, the equivalent plastic strain $\bar\varepsilon^p_n$, **and the
backstresses $\boldsymbol\alpha_{n,k}$ as SPATIAL tensors**. Because it owns
$\boldsymbol\alpha$, it can do the one thing the wrapper cannot: **co-rotate the
backstress** by the incremental material rotation, *before* running the return
map:

$$\boxed{\;\tilde{\boldsymbol\alpha}_{n,k} = \mathbf R_\Delta\,\boldsymbol\alpha_{n,k}\,\mathbf R_\Delta^{\mathsf T},
\qquad \mathbf R_\Delta = \operatorname{polar}(\mathbf f_\Delta) = \operatorname{polar}(\mathbf F\,\mathbf F_n^{-1}).\;}$$

This is the **§14.11 push-forward**. $\mathbf R_\Delta = \operatorname{polar}(\mathbf f_\Delta)$ is the incremental,
**incrementally-objective** material rotation (Hughes–Winget style):

- For a **pure rigid increment** $\mathbf f_\Delta = \mathbf Q$ it gives
  $\tilde{\boldsymbol\alpha} = \mathbf Q\boldsymbol\alpha\mathbf Q^{\mathsf T}$,
  **exactly cancelling** the rotation of $\mathbf s$ — so $\boldsymbol\xi = \mathbf s - \tilde{\boldsymbol\alpha}$ stays frame-consistent.
- As $\mathbf f_\Delta \to \mathbf I$, $\mathbf R_\Delta \to \mathbf I$, so the
  step reduces **exactly** to the small-strain / non-rotating case.

> [!note] The kernel runs verbatim — that is the whole payoff
> After the push-forward, `LadrunoJ2Finite` **reuses the verified return map
> `LadrunoJ2Kernel.h` VERBATIM**. The kernel is frame-agnostic — it never knew
> which frame it ran in — so the *only* genuinely new code is the
> $\boldsymbol\alpha$ push-forward plus carrying and serializing
> $\boldsymbol\alpha$ as spatial state. This is exactly why the kernel was
> extracted in PR #97: the small-strain material, every dimensional view, the
> wrapper path, and now the finite-native material all call **one** verified core,
> and there is never a second return map to drift out of sync (see
> [[Ladruno_materials_guide]] §2.1).

**Objectivity, measured (numpy oracle vs the v1 fixed-frame control):**

| Property | co-rotating $\boldsymbol\alpha$ (this material) | fixed-frame $\boldsymbol\alpha$ (v1 wrapper) |
|---|---|---|
| superposed rigid $\mathbf Q$ on a plastic state → $\boldsymbol\sigma_{\text{rot}} = \mathbf Q\boldsymbol\sigma\mathbf Q^{\mathsf T}$ | **3.2e-13** ✓ | 6.17 ✗ |
| continuously rotating + stretching plastic path → $\boldsymbol\sigma_{\text{rot}} = \mathbf R\,\boldsymbol\sigma_{\text{body}}\mathbf R^{\mathsf T}$ | **7.9e-13** ✓ | 2.37 ✗ |
| no-rotation (symmetric $\mathbf F$) path → reduces to v1 | 2.4e-12 ✓ | — |

So the polar push gives **full path objectivity** (not merely
superposed-rotation objectivity) and the correct small-strain reduction, to
$\sim10^{-12}$ — exactly the property v1 fails.

---

## 4. The verified algorithm

The per-Gauss-point update, as prototyped and pinned by the numpy oracle
(`tests/ladrunoj2_finite_native_reference.py`):

```text
state (per GP):  F_n,  Be_n,  ebarP_n,  alpha_n[N]      # α carried as SPATIAL tensors

setTrialF(F):
    f_Δ    = F · inv(F_n)
    Be_tr  = f_Δ · Be_n · f_Δᵀ                          # elastic predictor   (reuse LogStrainKernel)
    eps_tr = ½ ln(Be_tr)                                #   spectral; degeneracy branches handled
    R_Δ    = polar(f_Δ)                                 # incremental material rotation
    α̃_k    = R_Δ · alpha_n[k] · R_Δᵀ                    # ← THE NEW STEP (§14.11 push-forward)
    (τ, Δεᵖ, alpha_{n+1}, ebarP_{n+1}, Δγ, D)
           = LadrunoJ2Kernel::returnMap(eps_tr as elastic-trial, α̃, ebarP_n)   # VERBATIM kernel
    σ      = τ / J ;  spatial tangent via LogStrainKernel::spatial_tangent_full
    stage  Be_{n+1} = exp[2(eps_tr − Δεᵖ)],  alpha_{n+1} (CURRENT frame),  ebarP_{n+1}

commit:  promote the staged state.
```

Step-by-step: form the incremental gradient $\mathbf f_\Delta$, push the stored
$\mathbf b^e_n$ forward to the elastic predictor and take its log to get the trial
Hencky strain (reusing the `LogStrainKernel.h` spectral math, including the
repeated-eigenvalue degeneracy branches), compute the incremental material
rotation $\mathbf R_\Delta$, **co-rotate every backstress term** by it, then call
the **unchanged** `LadrunoJ2Kernel.h` return map with $\tilde{\boldsymbol\alpha}$
as the committed backstress and $\bar\varepsilon^p_n$ as the committed equivalent
plastic strain. The return map yields the Kirchhoff stress, the plastic-strain
increment, the updated backstresses, and the multiplier $\Delta\gamma$; the
material recovers Cauchy stress $\boldsymbol\sigma = \boldsymbol\tau/J$ and the
spatial tangent, and stages the updated $\mathbf b^e_{n+1}$,
$\boldsymbol\alpha_{n+1}$ (in the current frame), and $\bar\varepsilon^p_{n+1}$
for commit.

> [!check] Why the polar transport is the right (and sufficient) transport
> A **step-refinement study** (numpy, no build) drove the same return map along a
> severe **non-coaxial** path $\mathbf F(s) = \exp(s\,(\mathbf D + \mathbf W))$ —
> a **114° rotation + ~19% isochoric stretch** with the stretch-rate $\mathbf D$
> and spin $\mathbf W$ deliberately **not** coaxial (the regime where the
> transport choice actually matters; a fixed-axis stretch is degenerate). Result:
> the incremental $\operatorname{polar}(\mathbf f_\Delta)$ transport and the
> **exact §14.11 exponential-map** transport **share the same limit**, and their
> difference is **second order** in the step ($\approx O(\Delta s^2)$) — **3–4
> orders of magnitude below** the unavoidable **first-order** return-map
> discretization error. **Verdict: keep the incremental polar transport; the exact
> exp-map transport would not measurably improve accuracy.** (ADR P1, resolved
> 2026-06-02, no code change.)

---

## 5. The consistent tangent — two channels

For implicit Newton solvers the material must return the **algorithmic**
(consistent) tangent. Here the Cauchy stress depends on the trial $\mathbf F$
through **two** channels:

- **(A) the standard log-strain dependence.** The trial elastic strain
  $\boldsymbol\varepsilon^{e\,\text{tr}} = \tfrac12\ln(\mathbf f_\Delta\,\mathbf b^e_n\,\mathbf f_\Delta^{\mathsf T})$ — captured by the
  **already-shipped** spatial tangent (`LogStrainKernel::spatial_tangent_full`,
  the same one the wrapper uses).
- **(B) the co-rotated backstress.**
  $\tilde{\boldsymbol\alpha} = \mathbf R(\mathbf f_\Delta)\,\boldsymbol\alpha\,\mathbf R(\mathbf f_\Delta)^{\mathsf T}$ —
  $\tilde{\boldsymbol\alpha}$ depends on $\mathbf F$ through
  $\mathbf R = \operatorname{polar}(\mathbf f_\Delta)$. This channel is **new** to
  the native material and is exactly what a naive *"reuse the log-strain tangent"*
  would omit.

### 5.1 Channel B is small but must be included

| Property | Finding |
|---|---|
| channel-B size $\lVert\partial\boldsymbol\sigma/\partial\mathbf F_{\text{full}} - \partial\boldsymbol\sigma/\partial\mathbf F_{(\text{A only})}\rVert / \lVert\partial\boldsymbol\sigma/\partial\mathbf F\rVert$ | **$\sim10^{-4}$** (stiff metals, $G \gg \sigma_y$) up to **$\sim2\times10^{-3}$** (soft elasticity / strong kinematic) |
| vs. the tangent-test tolerance (rtol **2e-4**) | **straddles it** ⇒ the exact tangent **needs** channel B |
| effect on Newton robustness | none ($\ll1\%$ off ⇒ still $\sim$quadratic) |
| frame-objectivity | the channel-B fraction is **rotation-invariant** (objective) |
| non-zero only when | the step is **plastic** (the elastic predictor is $\boldsymbol\alpha$-independent) |

> [!important] Verdict: include channel B
> Channel B crosses the tangent-test tolerance, so it must be included for an
> exact consistent tangent. It does **not** hurt robustness (the omission was
> $\ll1\%$), but a tight FD-tangent gate sees it.

### 5.2 How channel B is computed — numeric then analytic

- **PR #127 (first ship): numeric R-perturbation.** Keep channel A **analytic**
  (the shipped spatial tangent); add channel B by perturbing **only**
  $\mathbf R = \operatorname{polar}(\mathbf f_\Delta)$ — hold
  $\boldsymbol\varepsilon^{e\,\text{tr}}$ and $J$ at the base point, perturb the
  6 (sym) / 9 directions, and finite-difference the extra Cauchy-stress response
  through the **unchanged** return map ($\approx$9–18 cheap scalar-Newton solves
  per GP-tangent, only when plastic). The two channels are proven cleanly additive
  to machine precision:
  $\partial\boldsymbol\sigma/\partial\mathbf F_{\text{full}} = \partial\boldsymbol\sigma/\partial\mathbf F_{(\text{A, R-frozen})} + \partial\boldsymbol\sigma/\partial\mathbf F_{(\text{B, R-only})}$.
- **PR #139 (default): the analytic $\partial\mathbf R/\partial\mathbf F$ chain.**
  The numeric R-perturbation is replaced by the analytic chain — the
  polar-rotation derivative ($\boldsymbol\Omega$ skew solving
  $(\operatorname{tr}\mathbf U\cdot\mathbf I - \mathbf U)\,\boldsymbol\omega = \operatorname{axial}(\mathbf A - \mathbf A^{\mathsf T})$,
  $\mathbf A = \mathbf R^{\mathsf T}\mathrm d\mathbf f$) composed with
  $\partial\tilde{\boldsymbol\alpha}/\partial\mathbf R$ and the return-map
  backstress sensitivity $\partial\boldsymbol\tau/\partial\tilde{\boldsymbol\alpha}$.
  No return-map calls, no FD step. The numeric path is **retained behind a compile
  flag** (`-DLADRUNO_J2FINITE_CHANNELB_NUMERIC`) for validation.

> [!check] Channel B is purely constitutive
> Channel B carries **no geometric term** — it is purely the constitutive
> sensitivity of the stress to the co-rotated backstress, mapped through the
> $\dot{\mathbf F} = \mathbf L\mathbf F$ relation
> ($c^{\text{B}}_{ijkl} = \sum_m (\partial\sigma_{ij}/\partial F_{km})\,F_{lm}$).
> This is what lets the C++ channel-B tensor be gated cleanly against the numpy
> oracle (commit-semantics-free), catching any $\mathbf F\to\mathbf L$
> index/factor/sign error.

---

## 6. IMPL-EX (`-implex`)

`-implex` (PR #134) adds the classic **Oliver–Huespe–Cante implicit/explicit
split**, here on the **plastic MULTIPLIER** $\Delta\gamma$ (**not** on a damage
variable — that is the small-strain `LadrunoJ2 -damage lemaitre -implex` Lemaitre
path).

### 6.1 The scheme

Every step the **implicit co-rotated return map still runs and is committed** — so
the committed finite history ($\mathbf b^e_n$, $\bar\varepsilon^p_n$,
$\boldsymbol\alpha_{n,k}$) is **byte-identical** to a fully-implicit run. Only the
stress/tangent **reported to the solver** are replaced, freezing two history
quantities: the extrapolated multiplier and the committed flow direction,
**co-rotated by the SAME $\mathbf R_\Delta = \operatorname{polar}(\mathbf f_\Delta)$**
as the backstress:

$$\tilde{\Delta\gamma} = \Delta\gamma_n\cdot\frac{\Delta t_{n+1}}{\Delta t_n}
\quad(\text{uniform step} \Rightarrow \Delta\gamma_n),
\qquad \tilde{\mathbf N}_n = \mathbf R_\Delta\,\mathbf N_n\,\mathbf R_\Delta^{\mathsf T},$$

$$\boldsymbol\varepsilon^{\tilde e} = \boldsymbol\varepsilon^{e\,\text{tr}} - \tilde{\Delta\gamma}\,\tilde{\mathbf N}_n,
\qquad \tilde{\boldsymbol\tau} = \mathbb C^e : \boldsymbol\varepsilon^{\tilde e},
\qquad \tilde{\boldsymbol\sigma} = \tilde{\boldsymbol\tau}/J.$$

### 6.2 The two clean consequences

1. **$\tilde{\mathbf N}_n$ co-rotates ⇒ $\tilde{\boldsymbol\sigma}$ is
   frame-indifferent** (objective to $\sim7\times10^{-13}$).
2. **The backstress drops out of $\tilde{\boldsymbol\sigma}$** (it appears only
   through the implicit $\mathbf N_n$, already frozen) ⇒ the **reported material
   tangent is the constant SPD elastic operator** — no plastic $h$, and the
   co-rotation **channel B vanishes**.

> [!note] Price and payoff
> **Price:** an $O(\Delta t)$ consistency error (the explicit stress → the
> implicit stress at **observed order $\approx2$** under step refinement for smooth
> flow) plus a **one-step lag at first yield** ($\Delta\gamma_n = 0 \Rightarrow \tilde{\boldsymbol\sigma}$ = elastic trial).
> **Payoff:** a **factor-once constant SPD tangent** for explicit / quasi-static
> use, and the SPD-tangent **enabler when finite J2 is later paired with a
> softening law** (finite Lemaitre). Plain hardening J2 already has an SPD implicit
> tangent, so IMPL-EX is not a robustness cure *on its own* for the hardening case.

> [!warning] Uniform-Δt assumption (documented)
> The extrapolation bakes in $\Delta t_{n+1}/\Delta t_n = 1$ (i.e.
> $\tilde{\Delta\gamma} = \Delta\gamma_n$). Accuracy degrades under variable
> $\Delta t$; the factory usage string and `Print` flag this. Wiring the dt-ratio
> is future work, blocked on materials not cleanly receiving $\Delta t$.

> [!note] Default off ⇒ bit-identical
> Without `-implex`, `LadrunoJ2Finite` is **bit-identical** to the fully-implicit
> material — IMPL-EX changes only the *reported* stress/tangent, never the
> committed history.

---

## 7. OpenSees implementation

### 7.1 Class and files

| Item | Value |
|---|---|
| Class | `LadrunoJ2Finite : public FiniteStrainNDMaterial` (the F-interface base, driven by `setTrialF(F)`) |
| classTag | **`ND_TAG_LadrunoJ2Finite = 33012`** (`SRC/classTags.h`; next free ND tag after `LadrunoJ2 = 33011`) |
| Material source | `SRC/material/nD/LadrunoJ2Finite.{h,cpp}` |
| Return-map kernel (reused verbatim) | `SRC/material/nD/LadrunoJ2Kernel.h` |
| Hencky kinematics + spatial tangent + degeneracy | `SRC/material/nD/LogStrainKernel.h` |
| Shared hardening law | `SRC/material/nD/LadrunoHardening.h` (`yieldStressVoceLinear`) |
| Polar rotation primitive | from `SRC/element/solidTransformation/SolidTransformationCorot` |
| Base + element seam | `FiniteStrainNDMaterial` (`setTrialF`, spatial-tangent accessors) + `LadrunoBrick -geom finite` |

> [!note] Reuse is the architecture
> Almost everything `LadrunoJ2Finite` needs already existed: the **return map**
> (`LadrunoJ2Kernel.h`), the **Hencky kinematics + spatial tangent + degeneracy
> branches** (`LogStrainKernel.h`), the **yield law** (`LadrunoHardening.h`), the
> **polar primitive** (from `SolidTransformationCorot`), and the **base + element
> seam** (`FiniteStrainNDMaterial` + `LadrunoBrick -geom finite`). The only
> genuinely new code is the $\boldsymbol\alpha$ push-forward, carrying
> $\boldsymbol\alpha$ as spatial state, and the channel-B tangent term.

### 7.2 The only consumer — and the element-level gating

`LadrunoJ2Finite` is **3D-only** (`getType() == "ThreeDimensional"`) and is driven
exclusively by `LadrunoBrick -geom finite` through a single `dynamic_cast` to
`FiniteStrainNDMaterial` (`setTrialF(F)`). Standard 2D/3D continuum elements drive
materials via `setTrialStrain` — disabled here — so there is no other consumer.

At the element level **only `-formulation std` is gated**; `bbar` (F-bar), `uri`,
and `eas` with the native material **work** (the element drives `setTrialF`
regardless of formulation) but are not all gated by tests.

### 7.3 Serialization

The backstresses are serialized **as spatial state** (`sendSelf`/`recvSelf`),
along with $\mathbf b^e_n$, $\bar\varepsilon^p_n$, and (for IMPL-EX) the committed
$\Delta\gamma_n$ and flow direction $\mathbf N_n$. The round-trip is tested
(`test_native_finite_database_roundtrip` — FE_Datastore round-trip of a committed
plastic finite state; skips without a DB).

### 7.4 Hygiene (from the 34-agent review)

The 34-agent adversarial review (7 dimensions) found **no confirmed correctness
bug**; the hardening it applied:

- **Per-instance `K0init`** — the elastic tangent matrix was a process-shared
  function-static `Matrix` in `getInitialTangent`; made per-instance.
- **Polar eigenvalue clamp** — `polarRotation` clamps eigenvalues so a
  near-singular channel-B perturbation cannot produce a rank-deficient / NaN
  $\mathbf R$.
- **`override` on all contract methods** — signature drift becomes a compile
  error, not a silent base-default fallback.
- **Serialization of $\boldsymbol\alpha$ as spatial state** round-trip tested
  (above), and a standalone g++ channel-B tensor check
  (`tests/ladrunoj2_finite_native_tangentB_check.cpp`) gates the constitutive
  channel-B against the numpy oracle to rel-Frobenius $<10^{-3}$.

> [!note] Why element-level FD-of-K is intentionally not gated
> Channel B scales with the committed backstress: `eleResponse("stiff")` re-forms
> post-commit ($\boldsymbol\alpha \ne 0$) while a from-virgin force-FD carries
> $\boldsymbol\alpha = 0$ — different committed bases, and an FE_Datastore replay
> fights OpenSees commit/sp semantics. Channel B is therefore gated **off-element**
> (material-level exact-to-FD + the g++ tensor check); the element side is covered
> by the v1 assembly gate and Newton convergence on the rotated-plastic
> objectivity path.

---

## 8. Command syntax and usage

### 8.1 Grammar

```tcl
nDMaterial LadrunoJ2Finite $tag  $K $G  -iso voce $sig0 $Qinf $b $Hiso  -kin $N $C1 $g1 …  <-implex>
element     LadrunoBrick    … $tag -geom finite       ;# the ONLY consumer (3D, F-interface)
```

The parameter set is the **same combined-hardening law as `LadrunoJ2`** (see
[[LadrunoJ2_guide]] §3 for the full meaning of each):

| Option | Arguments | Meaning |
|---|---|---|
| (positional) | `$K $G` | bulk and shear moduli |
| `-iso voce` | `$sig0 $Qinf $b $Hiso` | Voce + linear isotropic law $\sigma_y = \sigma_0 + Q_\infty(1 - e^{-b\bar\varepsilon^p}) + H_{\text{iso}}\bar\varepsilon^p$ |
| `-kin` | `$N  $C1 $g1  $C2 $g2 …` | $N$ Chaboche AF terms, each a $(C_k,\gamma_k)$ pair |
| `-implex` | (flag) | IMPL-EX Δγ-extrapolation → constant SPD elastic tangent (§6); default off ⇒ bit-identical to the implicit material |

> [!warning] Softening-parameter guard (PR #139)
> The factory **warns** when the minimum isotropic hardening slope
> $$\sigma_y'^{\,\min} = H_{\text{iso}} + \min(0,\,Q_\infty)\cdot b < 0$$
> because the consistent tangent may then become indefinite (below $-3G$ it loses
> positive-definiteness), recommending `-implex` (constant SPD tangent) or
> crack-band regularization. It is a **warning, not a rejection** — `-implex`
> makes softening well-posed.

### 8.2 Examples

```tcl
# --- Combined hardening (Voce + 3-term Chaboche), finite strain, objective ---
nDMaterial LadrunoJ2Finite 1  $K $G  -iso voce 350.0 60.0 10.0 0.0 \
                                     -kin 3  60000.0 500.0  12000.0 120.0  2000.0 10.0
element LadrunoBrick 101 $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8  1 -geom finite

# --- Same, with IMPL-EX (constant SPD tangent for explicit / quasi-static) ---
nDMaterial LadrunoJ2Finite 2  $K $G  -iso voce 350.0 60.0 10.0 0.0 \
                                     -kin 3  60000.0 500.0  12000.0 120.0  2000.0 10.0  -implex
element LadrunoBrick 102 $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8  2 -geom finite
```

```python
# OpenSeesPy — combined-hardening finite cyclic, consumed by LadrunoBrick -geom finite
import openseespy.opensees as ops
E, nu = 200000.0, 0.3
K, G = E/(3*(1-2*nu)), E/(2*(1+nu))

ops.nDMaterial("LadrunoJ2Finite", 1, K, G,
               "-iso", "voce", 350.0, 60.0, 10.0, 0.0,
               "-kin", 3, 60000.0, 500.0, 12000.0, 120.0, 2000.0, 10.0)
ops.element("LadrunoBrick", 101, *node_tags, 1, "-geom", "finite")
```

---

## 9. LadrunoJ2Finite vs the wrapper path

This is the single most important modeling decision around this material. Both
deliver finite-strain J2; they differ on **objectivity of kinematic hardening
under large rotation**.

| You need… | Use | Why |
|---|---|---|
| **Combined (kinematic) hardening AND large rotation** (finite cyclic, buckling-brace-like loops) | **`LadrunoJ2Finite`** + `-geom finite` | the backstress co-rotates ⇒ the Bauschinger / ratcheting response stays frame-indifferent (§3) |
| **Isotropic hardening at finite strain** (`-kin 0`) | **`LogStrain → LadrunoJ2`** (`-kin 0`) + `-geom finite` | the wrapper path is **already exact and objective** for the isotropic spine, and **simpler** |
| **Combined hardening, small rotation only** (or `-geom corot`) | plain **`LadrunoJ2`** | `corot` is objective for kinematic hardening (the de-rotated frame is a reference frame); the wrapper path is fine for small rotation |

> [!tip] The rule of thumb
> **Isotropic ⇒ wrapper (`LogStrain → LadrunoJ2 -kin 0`); combined hardening +
> large rotation ⇒ native (`LadrunoJ2Finite`).** The wrapper path is the simpler,
> exact choice everywhere *except* the combined-hardening-under-large-rotation
> corner this material was built to fill. See [[finite_strain_trifecta_guide]] §7.3
> and §9 for the full decision guide.

---

## 10. Verification and validation

**Status: MERGED.** v1 native material (classTag 33012) via **PR #127**, IMPL-EX
via **PR #134**, the analytic channel B + softening guard via **PR #139**. A
**34-agent adversarial review** (7 dimensions, each finding refuted-or-confirmed
by an independent skeptic) found **no confirmed correctness bug**.

| Test / oracle | What it pins |
|---|---|
| `test_native_objective_under_superposed_rotation` | $\boldsymbol\sigma_{\text{rot}} = \mathbf Q\boldsymbol\sigma\mathbf Q^{\mathsf T}$ on a plastic combined-hardening state — **the property v1 fails** (the passing native twin of the v1 strict `xfail`) |
| `ladrunoj2_finite_native_reference.py` | the numpy oracle: matches step-for-step along the rotating + stretching plastic path, reduces to v1 with no rotation |
| `test_ladrunoJ2_finite_native_tangent.py` (3) | the A+B tangent recipe exact-to-FD; channel B small-but-present; channel B objective |
| `ladrunoj2_finite_native_tangentB_check.cpp` + `…_tangentB_cpp.py` | the **C++ channel-B tensor** vs the oracle (rel-Frobenius $<10^{-3}$) — catches any $\mathbf F\to\mathbf L$ index/factor/sign error |
| `test_native_finite_database_roundtrip` | FE_Datastore round-trip of a committed plastic finite state (serialization of $\boldsymbol\alpha$ as spatial state) |
| step-refinement (`…_steprefine_reference.py` + test, 5) | polar transport converges at order $\approx1$; polar-vs-exact transport gap is $O(\Delta s^2)$, 3–4 orders below the discretization error (§4 callout) |
| IMPL-EX oracles (`…_finite_implex_reference.py` + test, 5) | committed history == implicit (byte); explicit→implicit order $\approx2$; objective $7\times10^{-13}$; constant SPD elastic tangent; finite→small-strain reduction; documented yield-onset lag |
| `ladrunoj2_finite_implex_check.cpp` + `…_implex_cpp.py` | standalone g++ mirror of the C++ `useImplex` math vs the oracle to $10^{-9}$ |
| analytic channel B (`ladrunoj2_finite_channelB_reference.py` + test, 5) | polar-deriv vs FD $\sim10^{-9}$; analytic vs numeric $\sim3\times10^{-8}$ (stiff & soft); A+B = full tangent; elastic ⇒ 0 |
| `test_ladrunoJ2Finite_element.py` | `LadrunoBrick -geom finite` over the native material (incl. `-implex`) — built-binary acceptance gate |

> [!check] The gate that defines done
> The v1 strict `xfail` in `test_ladrunoJ2_finite.py` (the LogStrain-wrapper
> non-objectivity under large rotation) **correctly stays `xfail`**, and gets a
> **passing native twin** in `test_native_objective_under_superposed_rotation`.
> The full J2 battery is green (34 passed / 1 xfailed).

---

## 11. Limitations and scope

> [!caution] Scope
> - **3D-only.** `getType() == "ThreeDimensional"`; the **only** consumer is
>   `LadrunoBrick -geom finite`. A plane-strain study is done with a 3D brick +
>   fixed out-of-plane DOFs.
> - **In scope:** isotropic + Chaboche AF kinematic (the full `LadrunoJ2` law),
>   the co-rotated backstress, the channel-A+B consistent tangent, IMPL-EX, and
>   serialization of $\boldsymbol\alpha$ as spatial state.

> [!caution] Deferred (out of scope) — each its own future PR
> - **Plane-stress / dimensional finite views** (§14.7 nested route) — **orphaned:
>   there is no 2D finite element consumer** (no element drives a
>   `FiniteStrainNDMaterial` in 2D), so these are deferred until one exists.
> - **Tabulated / Bézier isotropic curve** — gated on the small-strain `LadrunoJ2`
>   tabulated mode landing first (shared `LadrunoHardening.h`).
> - **Thermomechanical coupling.**
> - **`-formulation bbar/uri/eas` coverage** — the element drives `setTrialF`
>   regardless of formulation, so these should work with the native material, but
>   only `std` is gated by tests (`bbar+finite = F-bar` is shipped on the element
>   side).

> [!caution] Carried-in constitutive boundaries
> The deviatoric (pressure-independent) and rate-independent boundaries of the J2
> family carry in unchanged — see [[LadrunoJ2_guide]] §12. The **softening tangent**
> ($\sigma_y'^{\,\min} < 0$) caveat is now **guarded** with a factory warning
> (§8.1), recommending `-implex` or crack-band regularization.

> [!note] No 1D analogue exists
> Kernel-sharing is **nD↔finite-strain only, not 1D**. There is no
> `LadrunoJ2Finite` analogue for the uniaxial twin
> ([[LadrunoUniaxialJ2_guide]]) — a 1D backstress is a scalar with no frame to
> co-rotate.

---

## 12. References

**Theory**
- de Souza Neto, Perić & Owen (2008), *Computational Methods for Plasticity* —
  **Ch. 14** (finite-strain multiplicative / logarithmic-strain plasticity: §14.3
  log strain, §14.4 exponential map, Box 14.3 `MATISU`, §14.5 spatial tangent,
  §14.6 principal-stress form) and **§14.11** (the finite-strain
  kinematic-hardening boundary this material crosses). **The primary framework
  reference.**
- Oliver, Huespe & Cante — the IMPL-EX implicit/explicit integration split (§6).
- Hughes & Winget — incrementally-objective incremental rotation (the
  $\operatorname{polar}(\mathbf f_\Delta)$ transport).
- Armstrong & Frederick (1966); Chaboche (1986); Kobayashi & Ohno (2002) — the
  nonlinear-kinematic hardening law and its scalar return-map reduction (inherited
  from `LadrunoJ2`).
- Lee (1969); Miehe (1998) — the multiplicative split; closed-form repeated-
  eigenvalue treatment for the log-map tangent.

**Within this repo**
- [[16_finite_native_j2_adr]] — the full ADR: design, §14.11 boundary, the
  co-rotation algorithm, channels A/B, IMPL-EX, the PR history (#127/#134/#139) and
  follow-ups. **The decision log and oracle evidence.**
- [[LadrunoJ2_guide]] — the small-strain combined-hardening law this material
  inherits (yield surface, Voce/Chaboche, the scalar return map, the tangent). The
  partner reference (classTag 33011).
- [[finite_strain_trifecta_guide]] — where this material sits relative to the
  `LogStrain` wrapper path (§3) and the §14.11 objectivity boundary (§7.3).
- [[09_finite_strain_material_wrapper]] — the Hencky adaptor (`LogStrain`) design
  and the GREEN/YELLOW/RED material matrix.
- [[LadrunoBrick_reference]] — the only consumer (`-geom finite`, the F-interface).
- [[Ladruno_materials_guide]] — the material catalog (§3.3).
- [[LadrunoUniaxialJ2_guide]] — the uniaxial twin (immune to the §14.11 boundary).
- Source: `SRC/material/nD/LadrunoJ2Finite.{h,cpp}`, reusing `LadrunoJ2Kernel.h`,
  `LogStrainKernel.h`, `LadrunoHardening.h`, and the polar primitive from
  `SolidTransformationCorot`. classTag **33012**.
