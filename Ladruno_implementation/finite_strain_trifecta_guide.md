---
title: "Finite-Strain Trifecta — Element Transformation × Material Wrapper × Material"
project: Ladruno
type: reference / user guide
status: shipped (all three legs on `ladruno`; V&V P1–P5 complete)
related:
  - "[[LadrunoBrick_reference]]"
  - "[[LadrunoJ2_guide]]"
  - "[[solid_transformation_wrapper]]"
  - "[[09_finite_strain_material_wrapper]]"
  - "[[16_finite_native_j2_adr]]"
  - "[[17_finite_strain_validation_plan]]"
  - "[[18_finite_strain_validation_report]]"
tags:
  - finite-strain
  - large-deformation
  - element
  - material
  - wrapper
  - corotational
  - plasticity
  - reference
  - guide
---

# Finite-Strain Trifecta — the large-deformation stack

> [!abstract] One-paragraph summary
> Large-deformation simulation in the Ladruno fork is delivered not by one
> monolithic finite-strain element but by **three orthogonal, independently
> verifiable pieces that compose**: a **geometry-method layer** on the element
> (`SolidTransformation` → `LadrunoBrick -geom linear|corot|finite`), a
> **material-side finite-strain adaptor** (`nDMaterial LogStrain`, the Hencky /
> logarithmic-strain wrapper built on the `FiniteStrainNDMaterial` base), and a
> **flagship small-strain constitutive law** that feeds it (`nDMaterial
> LadrunoJ2`, combined-hardening von Mises). Stack them and you get genuine
> large-rotation **and** large-strain elastoplasticity —
> `element LadrunoBrick … -geom finite` over `nDMaterial LogStrain $t $j2` over
> `nDMaterial LadrunoJ2 …` — with **zero new return-map code** and every
> small-strain material in the library reusable verbatim. This guide is the
> single entry point: the theory of each leg, the architectural decisions that
> let them snap together at one clean seam, the OpenSees wiring, the intended use
> cases, and modeling/solver guidance.

This document **ties the three together**. Each leg also has its own deep-dive —
the element [[LadrunoBrick_reference]], the material [[LadrunoJ2_guide]], the
geometry layer [[solid_transformation_wrapper]], the wrapper
[[09_finite_strain_material_wrapper]] — and the end-to-end V&V record is
[[18_finite_strain_validation_report]]. Read this one first; descend into the
others for full per-leg detail.

---

## Contents

1. [[#1 · Why a trifecta and not a finite-strain element|Why a trifecta]]
2. [[#2 · Intended use cases|Intended use cases]]
3. [[#3 · The stack at a glance|The stack at a glance]]
4. [[#4 · The one idea — three seams, one contract|Three seams, one contract]]
5. [[#5 · Leg 1 — the element geometry layer (SolidTransformation)|Leg 1: SolidTransformation]]
6. [[#6 · Leg 2 — the material wrapper (LogStrain / Hencky)|Leg 2: LogStrain wrapper]]
7. [[#7 · Leg 3 — the constitutive law (LadrunoJ2)|Leg 3: LadrunoJ2]]
8. [[#8 · The seam-3 contract (what binds them)|The seam-3 contract]]
9. [[#9 · Putting it together — modeling recipes|Modeling recipes]]
10. [[#10 · Solver & analysis guidance|Solver guidance]]
11. [[#11 · Validation summary|Validation]]
12. [[#12 · Limitations & boundaries|Limitations]]
13. [[#13 · References|References]]

---

## 1 · Why a trifecta and not a finite-strain element

The textbook way to add large-strain plasticity is to write a dedicated
finite-strain element with a finite-strain material baked in. The fork
deliberately did **not** do that. The driving observation:

> [!important] The three concerns are mathematically independent
> 1. **How the element extracts a strain measure from nodal motion** (kinematics
>    + geometric nonlinearity) is an *element* concern.
> 2. **How a strain measure becomes stress** (the constitutive return map) is a
>    *material* concern, and the *small-strain* return map already exists,
>    verified, for dozens of materials.
> 3. **How a small-strain return map is lifted to finite strain** (the
>    logarithmic-strain pre/post-processing) is a *material-independent adaptor*
>    concern.

De Souza Neto, Perić & Owen (dSNPO 2008) Ch. 14, **Remark 14.5**, is the formal
license for the split: with the logarithmic (Hencky) elastic strain and the
exponential-map plastic update, *every* small-strain stress-update and tangent
subroutine is reused **without modification** — only kinematic pre/post-processing
is new, and that processing is material-independent. Abaqus uses exactly this for
its finite-strain plasticity.

So instead of one element × one material, the fork ships **three reusable layers**:

| Concern | Layer | Command surface |
|---|---|---|
| element kinematics / geometric nonlinearity | `SolidTransformation` (linear / corot / finite) | `element LadrunoBrick … -geom <…>` |
| small-strain → finite-strain lift | `LogStrainNDMaterial` (Hencky adaptor) | `nDMaterial LogStrain $t $inner` |
| the stress law itself | any GREEN 3D `NDMaterial`; flagship `LadrunoJ2` | `nDMaterial LadrunoJ2 …` |

**The payoff is combinatorial.** Any of the ~20 GREEN 3D materials (J2, soils,
Drucker–Prager, ASDPlastic, …) gets large-strain capability the moment
`LogStrain` wraps it — no per-material work. Any future Ladruno solid element
that consumes `SolidTransformation` gets all three geometry regimes for free.
And each leg is **testable in isolation**: the element against bit-identical
upstream bricks, the material against an independent return-map oracle, the
wrapper against a numpy Hencky/MATISU mirror.

---

## 2 · Intended use cases

The trifecta targets problems where **geometry changes enough that
small-displacement theory is wrong**, spanning two distinct regimes that the
`-geom` axis separates cleanly:

| Regime | `-geom` | Canonical problems | Material strain |
|---|---|---|---|
| **Large rotation, small strain** | `corot` | member/column buckling, P-Δ and soft-story sway, snap-through of shallow arches/shells, pounding rotation, slender cantilever roll-up, deployable/rotating mechanisms | < a few % |
| **Large strain (and rotation)** | `finite` | metal forming & upsetting, ductile **necking** and Taylor-bar impact, base-isolation elastomer pads, soft-soil large-deformation & post-failure, localization/shear bands, rubber-block compression | unrestricted (J2: exact incompressibility) |
| **Small everything** | `linear` (default) | ordinary stress analysis — the trifecta degrades to the plain verified small-strain hex | infinitesimal |

> [!tip] When you need each leg
> - **Only large rotation?** `-geom corot` + the *unchanged* small-strain
>   material. You do **not** need `LogStrain`. The full small-strain library
>   works as-is in the de-rotated frame.
> - **Large strain?** `-geom finite` + `LogStrain` wrapping your material. This is
>   the leg that *requires* the material wrapper.
> - **Plasticity with cyclic/Bauschinger realism?** Use `LadrunoJ2` as the inner
>   law (small or finite strain). For monotonic ductile metals, isotropic-only
>   `LadrunoJ2` is the clean, fully-objective finite-strain choice.

> [!caution] What it is **not** for
> - **Pressure-sensitive yield read on true Cauchy pressure at large volumetric
>   strain** (some soils) needs the deferred `-cauchyYield` J-seam flag — the
>   wrapper currently feeds the inner law Kirchhoff stress (J≈1 ⇒ non-issue for
>   metals; see §12).
> - **Plane-stress finite plasticity**, beam/shell finite kinematics, anisotropic
>   / single-crystal finite plasticity — out of scope (the wrapper needs a genuine
>   3D order-6 `σ(ε)` inner law fed a full **F**).
> - **Combined (kinematic) hardening under *large* rotation** — exact only for the
>   isotropic spine; the §14.11 objectivity boundary (see §12).

---

## 3 · The stack at a glance

```
        ┌──────────────────────────────────────────────────────────────┐
        │  element LadrunoBrick … $matTag  -geom finite                 │   LEG 1
        │  • kinematics ledger: F = I + Σ u_a ⊗ ∇N_a   (per Gauss pt)    │   element
        │  • SolidTransformationFinite (updated-Lagrangian)             │   geometry
        │  • assembles ∫ Bᵀσ dv + ∫ Bᵀc B dv + ∫ Gᵀ Σ G dv (K_geo)      │   layer
        └───────────────┬──────────────────────────────────────────────┘
                        │  setTrialF(F)            ▲ Cauchy σ, spatial tangent c
                        ▼                          │
        ┌──────────────────────────────────────────────────────────────┐
        │  nDMaterial LogStrain $t $inner   (FiniteStrainNDMaterial)     │   LEG 2
        │  • bᵉ_tr = F_Δ bᵉ_n F_Δᵀ  →  εᵉ_tr = ½ ln bᵉ_tr  (spectral)   │   material
        │  • owns committed bᵉ_n per GP; Box 14.3 MATISU                 │   wrapper
        │  • τ → σ = τ / det F ;  spatial tangent (1/2J)[D:L:B]          │   (Hencky)
        └───────────────┬──────────────────────────────────────────────┘
                        │  setTrialStrain(εᵉ_tr)   ▲ Kirchhoff τ, small-strain D
                        ▼                          │
        ┌──────────────────────────────────────────────────────────────┐
        │  nDMaterial LadrunoJ2 $K $G -iso voce … -kin 3 …              │   LEG 3
        │  • UNCHANGED small-strain return map (radial/AF return)        │   constitutive
        │  • header-only LadrunoJ2Kernel.h — one verified core           │   law
        └──────────────────────────────────────────────────────────────┘
```

Three commands, three layers, one data path. The **only** non-standard call in
the whole chain is `setTrialF(F)` (Leg 1 → Leg 2); everything below it is the
ordinary OpenSees `setTrialStrain → getStress/getTangent` cycle, which is *why*
Leg 3 needs no finite-strain awareness.

---

## 4 · The one idea — three seams, one contract

`LadrunoBrick` is carved into **three seams** so that geometric nonlinearity is
*additive*, not a rewrite. The trifecta lives across these seams:

| Seam | Isolates | `linear` | `corot` | `finite` |
|---|---|---|---|---|
| **1 · Kinematics ledger** | strain measure from $(\mathbf X,\mathbf u,\nabla N)$ at a GP | $\boldsymbol\varepsilon=\mathbf B\mathbf u$ | small $\boldsymbol\varepsilon$ in de-rotated frame | $\mathbf F=\mathbf I+\sum_a\mathbf u_a\otimes\nabla N_a$ |
| **2 · Geometry method** | de-rotate-in / re-rotate-out around an oblivious kernel | identity | EICR corotational + $\mathbf K_{\text{geo}}$ | updated-Lagrangian + spatial $\mathbf K_{\text{geo}}$ |
| **3 · Material adaptor** | element ↔ material boundary | `setTrialStrain(ε)` | `setTrialStrain(ε)` | **`setTrialF(F)`** → LogStrain |

The constitutive/integration **core** (the `std`/`bbar`/`uri`/`ssp`/`eas`
formulation logic) lives *inside* the seams and never learns which geometry
method wraps it. That orthogonality is the load-bearing design decision: it is
why `-geom corot -formulation bbar` needs no new kernel, and why `finite` reuses
every small-strain material through one adaptor.

> [!note] The formulation axis is orthogonal to the geometry axis
> `-formulation {std|bbar|uri|ssp|eas}` (anti-locking treatment) is independent of
> `-geom {linear|corot|finite}` (kinematic regime) — with **two genuine couplings
> under `finite`**: `bbar+finite = F-bar` (a seam-1 modification of **F**, not an
> oblivious wrapper), and `eas+finite` is stability-fragile (reserved). See
> [[LadrunoBrick_reference]] §8 and [[solid_transformation_wrapper]].

---

## 5 · Leg 1 — the element geometry layer (`SolidTransformation`)

### 5.1 What it is

A small interface — the **solid-element analogue of OpenSees' beam-column
`CrdTransf`/`CorotCrdTransf`**, which upstream has for *frames only*, never for
solids. It decouples *kinematics* from the *constitutive + integration core* and
does three jobs per residual/tangent evaluation: (1) build the strain measure the
core consumes from $(\mathbf X,\mathbf u,\nabla N)$; (2) de-rotate/localize global
kinematics into the core frame; (3) re-rotate/globalize the core's force/stiffness
back, **adding the geometric stiffness** the method implies.

| Method | Handles | Strain to material | Adds $\mathbf K_{\text{geo}}$ | Material library |
|---|---|---|---|---|
| `linear` | small disp + small rot | $\boldsymbol\varepsilon=\mathbf B\mathbf u$ | no | full small-strain ✓ |
| `corot` | **large rot**, small strain | small $\boldsymbol\varepsilon$ in de-rotated frame | yes (from $\partial\mathbf R/\partial\mathbf u$) | full small-strain ✓ |
| `finite` | **large strain** | $\mathbf F\to\varepsilon^e=\tfrac12\ln\mathbf b^e$ (via Leg 2) | yes (spatial) | small-strain *promoted by* Leg 2 |

### 5.2 Theory — the three methods

**`linear`** — identity. `localizeDisp`/`globalizeForce`/`globalizeStiff` are
pass-throughs; the element *is* the verified small-strain hex (regression-bit-
identical to upstream `Brick`/`bbarBrick`/`SSPbrick`).

**`corot`** — Element-Independent Corotational (EICR). A rotation $\mathbf R$
(iterative Higham polar of a centroidal **F** — chosen over eigen-of-**U** so
repeated stretches never bite) strips the rigid-body rotation; the unchanged
small-strain kernel runs in the corotated frame; forces/stiffness rotate back out
with the corotational geometric stiffness $\mathbf K_{\text{geo}}$ from
$\partial\mathbf R/\partial\mathbf u$:

$$
\mathbf f_{\text{global}}=\mathbf T^{\mathsf T}\mathbf f_{\text{local}},\qquad
\mathbf K_{\text{global}}=\mathbf T^{\mathsf T}\mathbf K_{\text{local}}\mathbf T+\mathbf K_{\text{geo}}.
$$

The **defining patch test**: zero stress under arbitrary rigid-body rotation.
Valid for large displacement + large rotation, *small* material strain — and it
reuses the **entire** small-strain material library and every core formulation.

**`finite`** — spatial multiplicative (updated-Lagrangian), **not**
total-Lagrangian. The earlier draft mixed a TL story (Green–Lagrange **E** → PK2
**S**) with a spatial log-strain adaptor; those don't compose. The ledger emits
the **full deformation gradient** at each GP,
$\mathbf F=\mathbf I+\sum_a\mathbf u_a\otimes\nabla_X N_a$, hands it to Leg 2 via
`setTrialF(F)`, receives **Cauchy σ** + the **spatial material tangent**, and
assembles spatially on the *current* configuration:

$$
\mathbf f_{\text{int}}=\int \mathbf B^{\mathsf T}\boldsymbol\sigma\,\mathrm dv,
\qquad
\mathbf K=\underbrace{\int \mathbf B^{\mathsf T}\mathbf c\,\mathbf B\,\mathrm dv}_{\text{material}}
+\underbrace{\int \mathbf G^{\mathsf T}\boldsymbol\Sigma\,\mathbf G\,\mathrm dv}_{\text{geometric }(\mathbf K_{\text{geo}})}.
$$

> [!important] Why `finite` is spatial, and why the small-strain map survives
> The route is Box 14.3 (`MATISU`):
> $\mathbf b^e_{\text{tr}}=\mathbf F_\Delta\mathbf b^e_n\mathbf F_\Delta^{\mathsf T}
> \to \varepsilon^e_{\text{tr}}=\tfrac12\ln\mathbf b^e_{\text{tr}}
> \to[\text{unchanged small-strain return map}]\to\boldsymbol\tau\to\boldsymbol\sigma=J^{-1}\boldsymbol\tau$.
> The **exponential-map** plastic update (dSNPO 14.72) is the load-bearing
> ingredient: it makes the small-strain map reusable *verbatim* (Remark 14.5) and
> keeps plastic incompressibility **exact** (Remark 14.4). Do **not** stack
> `corot` on `finite` — `finite` handles rotation intrinsically (spectral **F** is
> frame-indifferent); stacking double-counts it.

**F-bar.** `bbar + finite = F-bar`, the large-strain volumetric-locking cure
(dSNPO §15.1). It is **not** "the bbar kernel run obliviously" — it is a seam-1
modification of the gradient fed to the material:
$\bar{\mathbf F}=(J_0/J)^{1/3}\mathbf F$ with $J_0=\det\mathbf F_0$ at the element
centroid. The residual is unchanged (standard **B** on the actual config); the
tangent gains the coupling term (dSNPO eq. 15.10), which is **generally
unsymmetric** ⇒ needs an unsymmetric solver.

### 5.3 Architectural decisions

- **`SolidTransformation` as a class family**, not booleans:
  `SolidTransformation{,Linear,Corot,Finite}`, selected by id (0/1/2), rebuilt
  from the id in `recvSelf`. Element-agnostic — keyed only on nodal kinematics +
  ∇N + the core's local f/K — so *any* later Ladruno solid reuses it.
- **The geometric stiffness lives on the element, not the material.** dSNPO's
  "all-in-one spatial modulus **a**" is split: the material returns the
  *constitutive* part $\mathbf c=(1/2J)[D:L:B]$; the element owns
  $\int\mathbf G^{\mathsf T}\boldsymbol\Sigma\mathbf G\,\mathrm dv$ from the Cauchy
  stress and its *own* shape-function gradients (which the material cannot know).
  This keeps `LogStrainNDMaterial` element-agnostic and matches the conventional
  updated-Lagrangian split (Bonet & Wood). The two conventions are equivalent;
  the arbiter is the element's FD-tangent patch test on a rotated+stretched
  element.
- **F-bar lives in the element ledger**, not the transformation:
  `SolidTransformationFinite` stays a pure identity, and `LadrunoBrick` branches
  on `formulation==BBAR` inside the finite path (`centroidFbar()`).

### 5.4 OpenSees implementation (Leg 1)

- Files: `SRC/element/.../SolidTransformation{,Linear,Corot,Finite}.{h,cpp}`;
  consumed by `LadrunoBrick` (`ELE_TAG 33002`).
- Finite path: `deformationGradient()` (`LadrunoBrick.cpp:1094`),
  `centroidFbar()` (`:1131`), `updateFinite()` (`:1160`),
  `formResidAndTangentFinite()` (`:1233`).
- Parser-enforced combos: `corot`+`{std,bbar}`, `finite`+`{std (plain F), bbar
  (F-bar)}`; `finite` additionally requires the material be a
  `FiniteStrainNDMaterial` (else a construction-time error).
- Full element detail: [[LadrunoBrick_reference]] §8; design rationale and the
  seam math: [[solid_transformation_wrapper]].

---

## 6 · Leg 2 — the material wrapper (`LogStrain` / Hencky)

### 6.1 What it is

A **material-side adaptor** that lifts an *unchanged* small-strain 3D `NDMaterial`
to a genuine finite-strain (large rotation **and** large strain) material by the
logarithmic (Hencky) strain-space technique. The inner return map is reused
verbatim; the wrapper supplies the kinematic pre/post-processing around it
(dSNPO Box 14.3 `MATISU`). It is `nDMaterial LogStrain $tag $innerTag`
(`ND_TAG_LogStrainNDMaterial 33010`), a concrete subclass of the fork-local
abstract base `FiniteStrainNDMaterial` (which adds `setTrialF(F)`).

### 6.2 Theory — the Hencky lift

**Kinematics.** Multiplicative split (Lee 1969) $\mathbf F=\mathbf F^e\mathbf F^p$;
the Eulerian logarithmic (Hencky) elastic strain from the elastic left Cauchy–Green
tensor $\mathbf b^e=\mathbf F^e\mathbf F^{e\mathsf T}$:

$$
\boldsymbol\varepsilon^e=\tfrac12\ln\mathbf b^e
=\sum_a\ln(\lambda^e_a)\,\mathbf N_a\otimes\mathbf N_a .
$$

**Work-conjugacy — the fact that drives the whole design.** The stress conjugate
to the Hencky strain is the **Kirchhoff stress** $\boldsymbol\tau$, and the
elastic law $\boldsymbol\tau=\mathbf D^e:\boldsymbol\varepsilon^e$ uses the
*identical* matrix as infinitesimal isotropic elasticity. So feeding the inner
small-strain material the Hencky strain "as if it were total strain" returns
$\boldsymbol\tau$; Cauchy stress is recovered at the end as
$\boldsymbol\sigma=\boldsymbol\tau/\det\mathbf F$.

**Algorithm (Box 14.3, principal-stress form §14.6, what is actually coded):**

```
(i)   F update:        F_{n+1} = F_Δ · F_n            (element supplies full F)
(ii)  elastic predictor (push stored Hencky strain forward):
        bᵉ_tr   = F_Δ · bᵉ_n · F_Δᵀ                   spectral → eigenvalues b_i, vectors e_i
        εᵉ_tr_i = ½ ln b_i                            principal Hencky strains
(iii) UNCHANGED small-strain return map  (inner->setTrialStrain(εᵉ_tr)):
        εᵉ_{n+1} = εᵉ_tr − Δγ ∂Ψ/∂τ ,  return τ_i
(iv)  σ = (1/det F) Σ_i τ_i e_i⊗e_i                   Cauchy stress for internal force
```

**Spatial consistent tangent (§14.5):** the wrapper returns the *constitutive*
part $\mathbf c=(1/2J)[\,\mathbf D:\mathbf L:\mathbf B\,]$, where $\mathbf D$ is
the **reused small-strain consistent tangent**, $\mathbf L=\partial\ln[\mathbf
b^e_{\text{tr}}]/\partial\mathbf b^e_{\text{tr}}$ (the tensor-log derivative), and
$\mathbf B$ is the push-forward kinematic term. The geometric term
$-\sigma_{il}\delta_{jk}$ is the element's (see §5.3).

### 6.3 Architectural decisions

- **Fork-local base, additive registration only.** `setTrialF(F)` lives on
  `FiniteStrainNDMaterial` — **zero edit to vanilla `NDMaterial.h`**; inner
  materials untouched. The only vanilla touches are the unavoidable registration
  points (`classTags.h`, broker, command table, CMake) as strictly-additive
  `// Ladruno` lines.
- **The hidden state is $\mathbf b^e_n$**, not a plastic multiplier and not total
  strain (Box 14.3: $\mathbf b^e_n=\exp[2\varepsilon^e_n]$). The adaptor **owns**
  committed $\mathbf b^e_n$ + $\mathbf F_n$ per GP and **wraps** the inner
  material's `commitState`/`revertToLastCommit`/`sendSelf`/`recvSelf` — *this
  state-ownership coupling, not the stress formula, is the genuinely hard part*.
- **Repeated-eigenvalue degeneracy is the adaptor's responsibility, shipped from
  day one.** The log-map tangent contains
  $(\ln b_i-\ln b_j)/(b_i-b_j)$ off-diagonals that are 0/0 when principal stretches
  coincide — and that strikes the *simplest* states (uniaxial = 2-fold, hydrostatic
  = 3-fold, axisymmetric). The fix branches on multiplicity with the analytic
  Taylor limit (Miehe 1998 / dSNPO App. A.5 `ISO2`/`DISO2`), **not**
  ε-perturbation. The element's A2/A2′ tests (hydrostatic/equibiaxial) are the
  external net.
- **The J-seam is centralized and flagged.** All $1/J$ scaling happens at one
  `toInner`/`fromInner` seam, so a future `-cauchyYield` flag (divide by J before
  the inner yield check, for pressure-sensitive soils at large $J$) is ~10 lines,
  not a rebuild. Default: feed the inner law Kirchhoff $\boldsymbol\tau$
  (metals/moderate strain: $J\approx1$).

### 6.4 The material matrix — what lifts cleanly

The binding requirement: the inner law must be a genuine **3D order-6 symmetric
$\boldsymbol\sigma(\boldsymbol\varepsilon)$** (`getType()=="ThreeDimensional"`).

- 🟢 **GREEN (v1 targets):** `ElasticIsotropicThreeDimensional` (becomes Hencky
  hyperelasticity), `J2ThreeDimensional`/**`LadrunoJ2`** (the ideal case — exact
  plastic incompressibility), `DruckerPrager3D`, the `ASDPlasticMaterial3D`
  family, the 3D soils (`PressureDependMultiYield…`, `ManzariDafalias3D`, …),
  `ASDConcrete3D`.
- 🟡 **YELLOW (wraps, but a physics decision):** pressure-sensitive plasticity
  (Kirchhoff-vs-Cauchy pressure at large $J$), damage/softening (lch
  regularization keyed to the small-strain measure), rate/viscoplastic.
- 🔴 **RED (won't lift from the material side):** plane-stress-only (FSAM, MCFT
  panels), fiber/condensation adaptors (BeamFiber/PlateFiber), non-continuum
  (Acoustic, Contact), axisymmetric, `*Thermal`.

Full matrix: [[09_finite_strain_material_wrapper]].

### 6.5 OpenSees implementation (Leg 2)

- Files: `SRC/material/nD/FiniteStrainNDMaterial.{h,cpp}` (base),
  `SRC/material/nD/LogStrainNDMaterial.{h,cpp}` (adaptor), `LogStrainKernel`-style
  spectral math.
- Pattern copied from `MinMaxNDMaterial`/`InitStrainNDMaterial` (the
  NDMaterial-wraps-NDMaterial idiom: deep-copy inner, serialize inner classTag +
  state).
- `nDMaterial LogStrain $tag $innerTag` (Tcl + Python).
- Independent numpy oracle `tests/logstrain_reference.py` (spectral, the three
  eigenvalue branches, the spatial tangent vs FD) is the contract the C++ must
  match.

---

## 7 · Leg 3 — the constitutive law (`LadrunoJ2`)

`LadrunoJ2` is the **flagship small-strain inner law** of the trifecta — a single
self-contained **rate-independent von Mises (J2)** `nDMaterial` unifying nonlinear
**isotropic** (Voce + linear) and nonlinear **kinematic** (Chaboche /
Armstrong–Frederick) hardening. It is the OpenSees analogue of Abaqus `*PLASTIC,
COMBINED`. classTag `33011`. Full reference: [[LadrunoJ2_guide]].

### 7.1 Why it is the natural inner law

- **It presents the *exact* `J2ThreeDimensional` inner contract** the Hencky
  wrapper expects: engineering-shear strain in / true stress out, constant elastic
  $\mathbf C^e$, stateful $\boldsymbol\varepsilon^p$ subtraction, constant
  `getInitialTangent`. So `nDMaterial LogStrain $t $j2` lifts it to finite strain
  with **zero change to `LogStrainNDMaterial`**.
- **Exact plastic incompressibility** ($\det\mathbf F^p=1$) under the
  exponential-map update — the canonical clean finite-strain plasticity case.
- **One verified kernel.** The return map + tangent live in a header-only,
  OpenSees-free `LadrunoJ2Kernel.h`; the small-strain material, every dimensional
  view, and the finite path all call the *same* code — there is never a second
  implementation to drift.

### 7.2 Theory (capsule — full in [[LadrunoJ2_guide]])

Von Mises yield on the backstress-shifted deviator
$\boldsymbol\xi=\mathbf s-\boldsymbol\alpha$:
$f=\lVert\boldsymbol\xi\rVert-\sqrt{2/3}\,\sigma_y(\bar\varepsilon^p)$. Voce+linear
isotropic law $\sigma_y=\sigma_0+Q_\infty(1-e^{-b\bar\varepsilon^p})+H_{\text{iso}}\bar\varepsilon^p$;
Chaboche backstress
$\dot{\boldsymbol\alpha}_k=\tfrac23 C_k\dot{\boldsymbol\varepsilon}^p-\gamma_k\boldsymbol\alpha_k\dot{\bar\varepsilon}^p$.
The backward-Euler return map collapses to a **single scalar Newton** on $\Delta\gamma$
(Kobayashi–Ohno reduction $\mathbf n=\mathbf M/\lVert\mathbf M\rVert$); the
consistent tangent is assembled exactly (with the AF cross-term). Parameter limits
recover every simpler model: `-kin 0` ⇒ exactly upstream `J2Plasticity`.

### 7.3 The finite-strain objectivity boundary (load-bearing)

> [!warning] Combined hardening is *not* objective under large rotation
> The lift is **exact for the isotropic spine**: the elastic state co-rotates and
> isotropic yield depends only on $\lVert\mathbf s\rVert,\bar\varepsilon^p$, so a
> superposed rotation rotates Cauchy stress rigidly (verified to ~$10^{-9}$). But
> the **backstress $\boldsymbol\alpha$ is stored in the inner's *fixed* frame and
> does not co-rotate**, so $\lVert\mathbf s-\boldsymbol\alpha\rVert$ is not
> rotation-invariant. This is the dSNPO §14.11 framework boundary, pinned as a
> strict `xfail`. **Practical rule:** for `-geom finite`, use **isotropic
> hardening** (`-kin 0`) and the lift is exact at any deformation; combined
> hardening is valid in small-rotation / `corot` regimes. The planned fix is a
> finite-strain-native J2 that co-rotates $\boldsymbol\alpha$ each step — the
> kernel extraction is precisely its enabler ([[16_finite_native_j2_adr]]).

---

## 8 · The seam-3 contract (what binds them)

The entire trifecta hinges on **one interface call** at the element↔material
boundary. Getting this contract right (locked 2026-06-01, jointly by the element
and material efforts) is what lets the three legs ship independently.

```cpp
// Leg 1 (element, per Gauss point, finite path) calls:
mat->setTrialF(F);             // FULL deformation gradient F (3×3) in
sigma = mat->getStress();      // Cauchy σ = J⁻¹ τ   (6, Voigt)  → internal force ∫Bᵀσ dv
c     = mat->getTangent();     // SPATIAL CONSTITUTIVE tangent c = (1/2J)[D:L:B]
```

The five contract clauses, each a deliberate decision:

1. **Element passes the *full* F, not $\mathbf F_\Delta$.** The adaptor stores
   $\mathbf F_{\text{committed}}$, derives $\mathbf F_\Delta=\mathbf F\mathbf
   F_{\text{committed}}^{-1}$ itself, and holds $J=\det\mathbf F$ — so it owns the
   Kirchhoff↔Cauchy conversion at a single seam.
2. **`getStress()` returns Cauchy σ**, so the element assembles $\int\mathbf
   B^{\mathsf T}\boldsymbol\sigma\,\mathrm dv$ spatially with no pull-back, treating
   the adaptor as an ordinary `NDMaterial`.
3. **`getTangent()` returns the *constitutive* part $\mathbf c$ only.** The
   geometric/initial-stress term is the element's job (§5.3). Arbiter: the
   element's FD-tangent patch test.
4. **The adaptor owns committed $\mathbf b^e_n$ per GP** and wraps the inner
   material's commit/revert/serialize. The element commits the material as usual;
   the wrapper's `commitState` sets $\mathbf b^e_n\leftarrow\exp[2\varepsilon^e_{n+1}]$.
5. **The adaptor guarantees finite, non-NaN σ + tangent under repeated stretches**
   (degeneracy handling, §6.3). The element's A2/A2′ degeneracy battery is the net.

> [!note] Neither half is "finite strain" alone
> The element supplies the full **F** and the geometric stiffness; the adaptor
> supplies Cauchy stress + the spatial constitutive tangent and owns the
> multiplicative hidden state. They ship together but were *designed* separately
> against this contract — that is the whole architectural point.

---

## 9 · Putting it together — modeling recipes

### 9.1 The canonical large-strain elastoplastic stack

```python
import openseespy.opensees as ops

E, nu = 206900.0, 0.29
K, G = E/(3*(1-2*nu)), E/(2*(1+nu))

# Leg 3 — the small-strain law (ISOTROPIC for an objective finite-strain lift)
#   Simo necking law: sig_y = 450 + 265(1 - e^{-16.93 p}) + 129.24 p
ops.nDMaterial('LadrunoJ2', 1, K, G, '-iso', 'voce', 450.0, 265.0, 16.93, 129.24)

# Leg 2 — lift it to finite strain (Hencky / log-strain)
ops.nDMaterial('LogStrain', 2, 1)

# Leg 1 — F-delivering element; bbar ⇒ F-bar (necking is isochoric-plastic → std locks)
ops.element('LadrunoBrick', 101, *nodes, 2, '-formulation', 'bbar', '-geom', 'finite')

# F-bar tangent is UNSYMMETRIC → unsymmetric solver
ops.system('UmfPack')          # or 'FullGeneral'
```

```tcl
# Tcl mirror
nDMaterial LadrunoJ2 1 $K $G -iso voce 450.0 265.0 16.93 129.24
nDMaterial LogStrain 2 1
element LadrunoBrick 101 $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 2 -formulation bbar -geom finite
system UmfPack
```

### 9.2 Decision guide — which legs do I need?

| Problem | Legs | Element command |
|---|---|---|
| Ordinary stress analysis | Leg 3 only | `LadrunoBrick … $j2` (default `-geom linear`) |
| Column/member buckling, P-Δ, snap-through | Legs 1+3 (**no** wrapper) | `LadrunoBrick … $j2 -geom corot` |
| Metal forming / necking / Taylor bar | **all three** | `LadrunoBrick … $logstrain -geom finite -formulation bbar` |
| Hyperelastic rubber block (large strain, elastic) | Legs 1+2 over `ElasticIsotropic3D` | `LadrunoBrick … $logstrain -geom finite -formulation bbar` |
| Large-deformation soft soil | all three (inner = a 3D soil) | as above; consider future `-cauchyYield` |

> [!tip] Material-choice rules for `-geom finite`
> - **Isotropic hardening** (`-kin 0`) ⇒ fully objective at any deformation.
> - **Combined hardening** ⇒ only for small-rotation / `corot`; xfail at large
>   rotation (§7.3).
> - **Near-incompressible inner** (J2 at yield, ν→0.5) ⇒ use `-formulation bbar`
>   (F-bar) or the standard hex *will* lock — verified in B2/B3/E4.

### 9.3 Geometry / formulation cheat-sheet

- `-geom linear` — default; small everything. Any formulation.
- `-geom corot` — large rotation, small strain; `std`/`bbar` only; no wrapper
  needed.
- `-geom finite` — large strain; `std` (plain **F**) or `bbar` (**F-bar**); needs
  a `FiniteStrainNDMaterial` (`LogStrain`); `uri`/`ssp`/`eas` + finite reserved.

---

## 10 · Solver & analysis guidance

These are the hard-won, test-verified operational rules (logged in
[[LEDGER_quirks]] and [[18_finite_strain_validation_report]]):

> [!warning] Implicit bending into plasticity needs Krylov + EnergyIncr
> `-geom finite` + `LogStrain(LadrunoJ2)` bending into yield does **not** converge
> under plain `Newton` or `NewtonLineSearch` (residual grows from step 0).
> **`KrylovNewton` + `test EnergyIncr`** is robust. Also use a refined (≥2×2)
> cross-section — a 1-element-wide column bends too poorly for stable plastic
> Newton. Elastic, load-controlled finite runs are fine with plain `Newton`.

> [!danger] F-bar ⇒ unsymmetric tangent ⇒ unsymmetric solver
> `bbar+finite` (F-bar, dSNPO eq. 15.10) is generally unsymmetric. A symmetric
> solver (`BandSPD`/`ProfileSPD`) silently drops the coupling and breaks Newton
> convergence. Use `UmfPack`, `FullGeneral`, or `SparseGEN`. The factory emits a
> one-time advisory.

> [!warning] Explicit `criticalTimeStep()` is reference-config
> `criticalTimeStep()` does **not** shrink as elements compress (it is computed on
> the reference configuration — bit-identical before/after 33 % compression +
> 2× mushrooming in the Taylor-bar test). So explicit `-geom finite` must carry a
> safety factor < 1: 0.3·dt_cr is safe through full mushrooming; 0.5 is fine for
> short transits but risks instability through large strain. Pair with `-lumped`
> mass, `system Diagonal`, and `integrator CentralDifferenceLadruno`.

> [!note] EnergyBalance recorder sign on finite elements
> The `EnergyBalance` recorder reports internal energy with a **flipped sign** for
> the finite-strain element (KE column is correct). Compare `|IE|` to the kinetic
> energy change until `EnergyBalanceRecorder.cpp` internal-energy accumulation vs
> `LadrunoBrick -geom finite` resisting-force sign is audited.

---

## 11 · Validation summary

The trifecta is validated end-to-end in [[18_finite_strain_validation_report]] —
**41 Zone-A tests, all PASS** (Phases P1–P5), with independent numpy oracles
(`logstrain_reference.py`, `ladrunoj2_reference.py`) so an indexing bug cannot be
mirrored. Highlights:

| Phase | What it pins | Result |
|---|---|---|
| **P1 · L1 analytical (A1–A5)** | uniaxial/hydrostatic/equibiaxial/shear/dilatation vs closed form; degeneracy branches; exact plastic incompressibility ($\mathrm{tr}\,\boldsymbol\tau=3K\ln J=0$) | PASS (≤1e-5, A5 1e-7) |
| **P1 · L2 locking (B1–B3)** | h-convergence; near-incompressible block ν=0.4999 (std locks 0.10×, F-bar ν-insensitive); isochoric-J2 cantilever | PASS |
| **P1 · C1 necking** | Simo–Armero bar: neck localizes at imperfection, Considère geometric softening, ε̄ᵖ_neck ≫ ε̄ᵖ_end | PASS (physics; tight quantitative → Zone-B) |
| **P2 · geom-NL (A7,C4,C5,D4)** | corot elastica roll-up, large-rotation cantilever (≤2.3 %), Euler buckling, corot↔finite consistency | PASS |
| **P3 · locking (B4,E4)** | Cook's membrane (bbar/std > 1.3), rubber block (std/bbar ≈ 9×) | PASS |
| **P4 · explicit (Taylor bar)** | copper at 227 m/s: L_f/L₀=0.670 (~1.5 %), mushroom r/r₀=2.15 (~2 %), energy balance | PASS |
| **P5 · cross-val (D1,D2,D3,D5)** | `LadrunoBrick(linear,std)`↔`stdBrick` bit-identical; `LadrunoJ2`↔`J2Plasticity` bit-identical; std↔bbar on homogeneous F; hex↔tet bracket | PASS |

Plus the per-leg batteries: `LadrunoJ2` (`tests/test_ladrunoJ2_material.py`,
+ kernel oracle, + 33-agent adversarial review, **no correctness bugs**); the
element bit-for-bit anchors vs `Brick`/`bbarBrick`/`SSPbrick`; the wrapper's
spectral/degeneracy/tangent oracle. The full trifecta also passed a 56-agent
deep review (2026-06-02) that fixed COROT-1 (rotated body/gravity load), SEAM-1
(recorder returned Kirchhoff not Cauchy), PLUMB-1 (an `exit(-1)` on a bad inner),
and GEOM-1 (finite-mat + linear/corot guard) — see [[project_finite_strain_trifecta_review]].

---

## 12 · Limitations & boundaries

> [!caution] Carried-in boundaries (per leg)
> - **Combined hardening + large rotation** — the §14.11 objectivity boundary
>   (§7.3). Use isotropic hardening with `-geom finite`; combined hardening is
>   exact only in small-rotation / `corot`. Pinned xfail; v2 native-J2 fix planned.
> - **Pressure-sensitive yield at large volumetric strain** — the wrapper feeds
>   the inner law Kirchhoff $\boldsymbol\tau$ (J≈1 for metals). Soils needing the
>   yield read on true Cauchy pressure want the deferred `-cauchyYield` J-seam flag.
> - **F-bar tangent is unsymmetric** — needs an unsymmetric solver throughout
>   locking/necking studies.
> - **`uri`/`ssp`/`eas` under `corot`/`finite` are reserved** — large-strain
>   enhanced-strain (`eas+finite`) hourglasses in compression (dSNPO §15.2.5).
> - **Plane-stress / beam-shell finite kinematics, anisotropic & single-crystal
>   finite plasticity** — out of scope (the wrapper needs a 3D order-6 inner law
>   driven by a full **F**).
> - **Explicit dt_cr is reference-config** — carry a safety factor < 1 (§10).
> - **C1 tight quantitative neck-ratio** match to Simo needs a mesh-converged fine
>   model (Zone-B); the coarse structured lattice reproduces the physics, not the
>   3 % radius tolerance.

---

## 13 · References

**Theory**
- de Souza Neto, Perić & Owen (2008), *Computational Methods for Plasticity* —
  **Ch. 14** (finite-strain multiplicative plasticity: §14.3 log strain, §14.4
  exponential map, Box 14.3 `MATISU`, §14.5 spatial tangent, §14.6 principal-stress
  form, §14.11 kinematic-hardening boundary) and **Ch. 15** (§15.1 F-bar eq. 15.5/
  15.10, §15.2 EAS). The primary framework reference.
- Belytschko, Liu, Moran & Elkhodary, *Nonlinear Finite Elements for Continua and
  Structures* — corotational (§4.6), hourglass/assumed strain (§8.6–8.7).
- Bonet & Wood, *Nonlinear Continuum Mechanics for FEA* — the updated-Lagrangian
  material/geometric tangent split.
- Lee (1969) — multiplicative split $\mathbf F=\mathbf F^e\mathbf F^p$.
- Miehe (1998) — closed-form derivative of isotropic tensor functions / repeated
  eigenvalue treatment for the log-map tangent.
- Simo & Miehe (1992); Eterovic & Bathe (1990); Perić & Owen (1991) — log-strain /
  exponential-map return mapping.
- Armstrong & Frederick (1966); Chaboche (1986); Kobayashi & Ohno (2002) — the
  nonlinear-kinematic hardening law and its scalar return-map reduction.

**Within this repo**
- [[LadrunoBrick_reference]] — Leg 1, the element (geometry seam §8).
- [[solid_transformation_wrapper]] — Leg 1 design rationale, the three methods.
- [[09_finite_strain_material_wrapper]] — Leg 2, the Hencky adaptor (full material
  matrix, decisions D1–D4).
- [[LadrunoJ2_guide]] — Leg 3, the constitutive law (full theory + algorithm).
- [[16_finite_native_j2_adr]] — the planned finite-strain-native J2 (the §14.11 v2).
- [[17_finite_strain_validation_plan]] / [[18_finite_strain_validation_report]] —
  the V&V plan and execution record.
- Source: `SRC/element/ladrunoBrick/` (33002), `SRC/element/.../SolidTransformation*`,
  `SRC/material/nD/{FiniteStrainNDMaterial,LogStrainNDMaterial}.{h,cpp}` (33010),
  `SRC/material/nD/LadrunoJ2.{h,cpp}` + `LadrunoJ2Kernel.h` (33011).

---

> [!info] Maintenance
> This guide is the **single entry point** for the finite-strain stack. When any
> leg changes (a new `-geom` combo, the §14.11 v2 native-J2, the `-cauchyYield`
> flag), update here *and* the per-leg doc, the [[LEDGER_implementations]] rows
> (33002/33010/33011), and the banner. Keep the bit-for-bit oracles
> (`stdBrick`/`J2Plasticity`) and the numpy mirrors as the correctness anchors.
