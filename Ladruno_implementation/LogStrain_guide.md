---
title: "LogStrain — Hencky (Logarithmic) Finite-Strain Material Wrapper"
project: Ladruno
type: reference / user guide
status: shipped (v1, elastic + GREEN plastic inner; isotropic spine objective)
classTag: 33010
material: nDMaterial LogStrain
base: FiniteStrainNDMaterial
related:
  - "[[finite_strain_trifecta_guide]]"
  - "[[09_finite_strain_material_wrapper]]"
  - "[[LadrunoJ2_guide]]"
  - "[[16_finite_native_j2_adr]]"
  - "[[LadrunoBrick_reference]]"
  - "[[solid_transformation_wrapper]]"
  - "[[Ladruno_materials_guide]]"
tags:
  - material
  - wrapper
  - finite-strain
  - large-deformation
  - hencky
  - log-strain
  - reference
---

# LogStrain — Hencky (Logarithmic) Finite-Strain Material Wrapper

> [!abstract] One-line summary
> `LogStrain` is a **material-side adaptor** that lifts an *unchanged* small-strain
> 3D `nDMaterial` to a genuine **finite-strain** model — large rotation **and**
> large strain — by the logarithmic (Hencky) strain-space technique. The inner
> material's return map is reused **verbatim** (de Souza Neto, Perić & Owen Box
> 14.3, *MATISU*); the wrapper supplies only the material-independent kinematic
> pre-/post-processing around it. It is the OpenSees analogue of how **Abaqus**
> obtains its finite-strain plasticity, and the constitutive-lift partner of
> [[LadrunoJ2_guide|LadrunoJ2]] and [[LadrunoBrick_reference|LadrunoBrick `-geom
> finite`]].

> [!note] Scope of this document — the material-side deep dive
> This guide is the **material-focused** reference for the `LogStrain` /
> `FiniteStrainNDMaterial` adaptor: the Hencky theory, the spectral algorithm
> exactly as coded, the spatial tangent, the GREEN/YELLOW/RED inner-material
> matrix, the hidden state, the degeneracy net, and the OpenSees wiring. For the
> **whole large-deformation stack** — the element geometry layer
> (`SolidTransformation`), this wrapper, and the constitutive law, plus the
> end-to-end validation record — start at the
> [[finite_strain_trifecta_guide|finite-strain trifecta guide]]. For the
> chronological design/decision log see [[09_finite_strain_material_wrapper]].

---

## Contents

- [[#1 Intended use cases|1. Intended use cases]]
- [[#2 Why it exists — what drove the design|2. Why it exists]]
- [[#3 Continuum theory — the Hencky lift|3. Continuum theory]]
- [[#4 Numerical formulation — Box 14.3 MATISU|4. Numerical formulation]]
- [[#5 The consistent spatial tangent|5. Spatial tangent]]
- [[#6 The hidden state and commit cycle|6. Hidden state]]
- [[#7 Repeated-eigenvalue degeneracy|7. Degeneracy]]
- [[#8 The material matrix — what lifts cleanly|8. Material matrix]]
- [[#9 OpenSees implementation|9. OpenSees implementation]]
- [[#10 Command syntax and usage|10. Command syntax]]
- [[#11 Verification and validation|11. Verification]]
- [[#12 Limitations and boundaries|12. Limitations]]
- [[#13 References|13. References]]

---

## 1. Intended use cases

`LogStrain` exists to give the **entire 3D small-strain material library**
large-deformation reach without touching a single return map. It is the right
tool whenever **strain (not just rotation) gets large** and the constitutive law
you want is already a verified small-strain 3D `nDMaterial`.

| Use case | Why LogStrain |
|---|---|
| **Metal forming, upsetting, ductile necking, Taylor-bar impact** | Wraps `J2ThreeDimensional` / [[LadrunoJ2_guide\|LadrunoJ2]] into multiplicative finite-strain J2 with **exact** plastic incompressibility — the canonical large-strain plasticity case. |
| **Large-strain elasticity (rubber blocks, elastomer pads)** | Wrapping `ElasticIsotropicThreeDimensional` yields a **Hencky hyperelastic** law — a real finite-strain elastic model, well-behaved in compression (not the spurious-softening StVK). |
| **Large-deformation / post-failure soft soils** | Lifts the 3D soil family (`PressureDependMultiYield…`, `ManzariDafalias3D`, Drucker–Prager) to finite strain — with a Kirchhoff-vs-Cauchy-pressure caveat at large $J$ (§8, §12). |
| **Concrete continuum at large strain / localization** | Wraps `ASDConcrete3D` (a YELLOW case — `lch` regularization is keyed to the small-strain measure). |
| **Any future GREEN 3D `nDMaterial`** | The lift is per-*class*-free: a new 3D order-6 law gets finite strain the moment `LogStrain` wraps it — no per-material code. |

> [!tip] When you do **not** need `LogStrain`
> If your problem is **large rotation but small strain** (column buckling, P-Δ,
> snap-through), you do **not** need this wrapper — use `LadrunoBrick -geom corot`
> with the *unchanged* small-strain material directly. `LogStrain` is specifically
> the leg that the **`-geom finite`** (large-strain) regime requires. And do
> **not** stack `corot` on `finite`: the log-strain path handles rotation
> intrinsically (spectral $\mathbf F$ is frame-indifferent), so stacking
> double-counts it.

---

## 2. Why it exists — what drove the design

Before `LogStrain` the fork had nonlinear *elements* and explicit dynamics, but
**every continuum material was small-strain**: large-deformation soil, metal
forming, necking, rubber compression, and finite-rotation column buckling *with
yielding* were all out of reach. There are three classical wrapper families for
closing that gap, and only one is right:

| Family | Route | Verdict |
|---|---|---|
| **(a) St. Venant–Kirchhoff** | feed Green–Lagrange $\mathbf E$ to the elastic $\mathbf D$, $\mathbf S=\mathbf D{:}\mathbf E$ | degenerate elastic-only; **spurious softening in compression**; useless for plasticity |
| **(b) Hypoelastic / Jaumann** | corotational rate form | documented shear-stress **oscillation** under large simple shear (dSNPO §14.10.3); only a small-elastic-strain corotational engine, not a true finite-strain solution |
| **(c) Logarithmic / Hencky** | spectral $\tfrac12\ln\mathbf b^e$ + exponential-map plastic update | **the standard, validated technique** — reuses every small-strain return map verbatim |

`LogStrain` is family **(c)**. The reason it is the highest-leverage choice is
structural: it **reuses the entire existing 3D material library** (J2,
Drucker–Prager, the soils, ASDPlastic, concrete-3D) without rewriting a single
return map. Abaqus does exactly this for its finite-strain plasticity; the
technique is validated across the literature (Eterovic & Bathe 1990; Weber &
Anand 1990; Perić & Owen 1991; Simo & Miehe 1992; Cuitiño & Ortiz 1992).

> [!important] The formal license — dSNPO **Remark 14.5**
> With the logarithmic (Hencky) elastic strain **and** the exponential-map plastic
> update, *every* small-strain stress-update **and** tangent subroutine is reused
> **without modification**. Only the kinematic pre/post-processing of Box 14.3
> steps (i), (ii), (iv) is new — and that processing is **material-independent**.
> This is the single fact that makes the whole "wrapper" framing literal.

> [!note] Design philosophy — one adaptor, every material
> `LogStrain` is a pure spatial map $\mathbf F \mapsto (\text{Cauchy }
> \boldsymbol\sigma,\ \text{spatial tangent } \mathbf c)$ wrapping an **unchanged**
> small-strain `NDMaterial`. The inner law never learns it is being driven at
> finite strain — it only ever sees `setTrialStrain → getStress/getTangent`. That
> orthogonality is the load-bearing decision: it is why **any** GREEN 3D material
> gets large-strain capability for free, and why the adaptor is element-agnostic
> and reusable.

---

## 3. Continuum theory — the Hencky lift

> [!info] Notation
> Direct (bold), index, and Voigt forms are shown side by side. Upper-case indices
> $I,J,K,L$ are reference, lower-case $i,j,k,l$ are spatial. Voigt ordering is the
> OpenSees convention $\{11,22,33,12,23,31\}$ with **engineering** shear
> $\gamma_{ij}=2\varepsilon_{ij}$ on the off-diagonal slots (the conversion to/from
> true-tensor shear happens at the seam to the inner material).

### 3.1 Kinematics — multiplicative split and the Hencky strain

The deformation gradient splits multiplicatively (Lee 1969) into elastic and
plastic parts:

- **Direct:** $\mathbf F = \mathbf F^e\,\mathbf F^p$
- **Index:** $F_{iJ} = F^e_{iK}\,F^p_{KJ}$

The wrapper works with the **Eulerian logarithmic (Hencky) elastic strain**, built
from the **elastic left Cauchy–Green tensor** $\mathbf b^e = \mathbf F^e\mathbf
F^{e\mathsf T} = (\mathbf V^e)^2$:

$$
\boldsymbol\varepsilon^e
= \ln\mathbf V^e
= \tfrac12\ln\mathbf b^e
= \sum_{a=1}^{3}\ln(\lambda^e_a)\,\mathbf N_a\otimes\mathbf N_a,
\qquad
\mathbf b^e = \sum_{a=1}^{3}(\lambda^e_a)^2\,\mathbf N_a\otimes\mathbf N_a,
$$

where $\lambda^e_a$ are the elastic principal stretches and $\mathbf N_a$ the
spatial principal directions. The volumetric part is
$\varepsilon^e_v = \operatorname{tr}\boldsymbol\varepsilon^e = \ln J^e$,
$J^e = \det\mathbf F^e$.

> [!important] Why J2 is the clean case — exact plastic incompressibility
> A **traceless plastic flow vector ⇒ isochoric plastic flow, exactly** — just as
> in small strain (dSNPO 14.64–14.65). Combined with the **exponential-map**
> plastic update this keeps $\det\mathbf F^p = 1$ to machine precision, not just to
> integration order. That is why von Mises (J2) is the canonical, exact
> finite-strain target: the deviatoric small-strain return map carries over with
> its incompressibility intact.

### 3.2 Constitutive law — Kirchhoff stress is conjugate to Hencky strain

The free energy is written on the log strain, $\psi(\boldsymbol\varepsilon^e,
\alpha)$. The stress that falls out is the **Kirchhoff stress** $\boldsymbol\tau$
— *this is the work-conjugacy fact that drives the entire design*:

- **Direct:** $\boldsymbol\tau = \bar\rho\,\partial\psi/\partial\boldsymbol\varepsilon^e$
- **Hencky elastic law (the reused small-strain $\mathbf D$):**
  $\boldsymbol\tau = \mathbf D^e:\boldsymbol\varepsilon^e$ — the **identical**
  matrix as infinitesimal isotropic elasticity.

So feeding the inner small-strain material the Hencky strain "as if it were total
strain" returns the **Kirchhoff** stress. Yield and flow are posed in Kirchhoff
stress, $\Phi(\boldsymbol\tau,A)\le 0$, $\tilde{\mathbf D}^p =
\dot\gamma\,\partial\Psi/\partial\boldsymbol\tau$ — again format-identical to the
small-strain scheme. The **Cauchy** stress (what the element needs for the spatial
internal force) is recovered only at the very end:

$$
\boxed{\ \boldsymbol\sigma = \dfrac{\boldsymbol\tau}{\det\mathbf F}\ }
$$

> [!note] This is *literally* "the small-strain material, reused"
> Because $\boldsymbol\tau = \mathbf D^e:\boldsymbol\varepsilon^e$ is the same
> matrix algebra as small strain, the inner material's `setTrialStrain →
> getStress` cycle — fed the Hencky strain — returns exactly the Kirchhoff stress
> the finite theory needs. The discrete return map (dSNPO 14.89–14.90) is
> functionally identical to the infinitesimal scheme (7.25). Nothing inside the
> inner law changes.

---

## 4. Numerical formulation — Box 14.3 MATISU

The integration is **spatial multiplicative** (updated-Lagrangian), *not*
total-Lagrangian — Green–Lagrange $\mathbf E$ / PK2 $\mathbf S$ (TL) and a spatial
log-strain adaptor do not compose. The coded path is dSNPO Box 14.3 (`MATISU`) in
its **principal-stress form** (§14.6), which sidesteps assembling 4th-order
log-map tensors during the stress update.

Given the full deformation gradient $\mathbf F_{n+1}$ from the element (the
adaptor derives the relative gradient $\mathbf F_\Delta = \mathbf F_{n+1}\mathbf
F_n^{-1}$ itself):

```
(i)   Kinematics — F update:
        F_Δ      = F_{n+1} · F_n⁻¹                       (element supplies full F)

(ii)  Elastic predictor — push the stored Hencky strain forward:
        bᵉ_n          = exp[2 εᵉ_n]                       (stored hidden state)
        bᵉ_tr         = F_Δ · bᵉ_n · F_Δᵀ                 (14.94)
        spectral decompose bᵉ_tr → eigenvalues b_i, eigenvectors e_i
        εᵉ_tr,i       = ½ ln b_i                          (14.104, principal Hencky strain)

(iii) GOTO the UNCHANGED small-strain return map (inside the inner material):
        inner->setTrialStrain( assembleVoigt(εᵉ_tr, e_i) )   → τ_i  (Kirchhoff)
        εᵉ_{n+1} = εᵉ_tr − Δγ ∂Ψ/∂τ|_{n+1}               (14.89/14.90)
      ** This is exactly OpenSees' existing setTrialStrain → getStress. **

(iv)  Cauchy stress for the internal force:
        σ = (1/det F_{n+1}) · Σ_i τ_i e_i⊗e_i             (14.105 / Box 14.3 iv)
```

The **exponential-map** discretisation of the flow rule (dSNPO 14.72–14.73) is the
load-bearing ingredient at step (iii): it is what makes the inner return map
format-identical to small strain **and** keeps $\det\mathbf F^p = 1$ exactly for
traceless flow (a plain backward-Euler, 14.74, loses volumetric accuracy). The
plastic update never needs $\mathbf F^p$ explicitly; the wrapper stores the
elastic state $\mathbf b^e_n$ (equivalently $\boldsymbol\varepsilon^e_n$) as
history (§6).

> [!note] The principal-stress form is what is coded
> Spectrally decomposing $\mathbf b^e_{\text{tr}} = \sum_i b_i\,\mathbf
> e_i\otimes\mathbf e_i$, taking principal Hencky strains $\varepsilon_i^{\text{tr}}
> = \tfrac12\ln b_i$, running the small-strain return map in principal space to get
> $\tau_i$, and reassembling $\boldsymbol\sigma = J^{-1}\sum_i\tau_i\,\mathbf
> e_i\otimes\mathbf e_i$ avoids the 4th-order tensor-log machinery in the *stress*
> path. It reappears only in the *tangent* (§5).

---

## 5. The consistent spatial tangent

For implicit (Newton) global solvers the wrapper must return a tangent consistent
with the discrete update. The full spatial modulus is (dSNPO 14.99):

$$
a_{ijkl} = \underbrace{\frac{1}{2J}\,[\,\mathbf D:\mathbf L:\mathbf B\,]_{ijkl}}_{\text{constitutive }\mathbf c}
\;-\; \underbrace{\sigma_{il}\,\delta_{jk}}_{\text{geometric / initial-stress}},
$$

with the three constitutive factors:

| Factor | Definition (dSNPO) | Role |
|---|---|---|
| $\mathbf D = \partial\boldsymbol\tau/\partial\boldsymbol\varepsilon^e_{\text{tr}}$ | 14.100 | the **reused small-strain consistent tangent** — the *only* material-dependent term |
| $\mathbf L = \partial\ln[\mathbf b^e_{\text{tr}}]/\partial\mathbf b^e_{\text{tr}}$ | 14.101 | the tensor-log derivative (the "intricate transformation"; **care at repeated eigenvalues**, §7) |
| $\mathbf B$, $B_{ijkl}=\delta_{ik}(\mathbf b^e_{\text{tr}})_{jl}+\delta_{jk}(\mathbf b^e_{\text{tr}})_{il}$ | 14.102 | the push-forward kinematic term |

> [!important] Seam contract — the wrapper returns $\mathbf c$, **not** the full $\mathbf a$
> The split of 14.99 across the element↔material seam is **locked** (2026-06-01,
> element + material teams agreed): **`LogStrain::getTangent()` returns the
> constitutive part $\mathbf c = (1/2J)[\,\mathbf D:\mathbf L:\mathbf B\,]$ only.**
> The **element** owns the geometric/initial-stress term $-\sigma_{il}\delta_{jk}$,
> assembling it as $\int \mathbf G^{\mathsf T}\boldsymbol\Sigma\,\mathbf G\,\mathrm
> dv$ from the Cauchy stress and its *own* shape-function gradients (which the
> material cannot know). This keeps `LogStrainNDMaterial` element-agnostic and
> matches the conventional updated-Lagrangian split (Bonet & Wood); it is
> equivalent to dSNPO's "all-in-one $\mathbf a$" recipe. A separate accessor,
> `getSpatialTangentTensor`, returns the full $\mathbf a$ for callers that want it.
> **Arbiter:** the element's finite-difference tangent patch test on a single
> rotated + stretched element.

> [!note] One-line summary (Remark 14.5, again)
> The only material-dependent object in the whole tangent is $\mathbf D$, the
> *small-strain consistent tangent the inner material already computes*. The
> log-map factors $\mathbf L$ and $\mathbf B$ are material-independent kinematics.

---

## 6. The hidden state and commit cycle

This is the genuinely hard part of the adaptor — not the stress formula, but the
**state ownership**.

> [!important] The wrapper owns $\mathbf b^e_n$ — not a plastic multiplier, not total strain
> The committed history `LogStrain` stores per Gauss point is the **elastic left
> Cauchy–Green tensor** $\mathbf b^e_n = \exp[2\,\boldsymbol\varepsilon^e_n]$ plus
> the committed gradient $\mathbf F_n$. It is **not** a plastic multiplier and
> **not** the total strain. The wrapper **wraps** the inner material's own state:
> it delegates `commitState` / `revertToLastCommit` / `sendSelf` / `recvSelf` to
> the inner law and manages $\mathbf b^e_n,\mathbf F_n$ around it. This copies the
> NDMaterial-wraps-NDMaterial idiom of `MinMaxNDMaterial` and
> `InitStrainNDMaterial` (deep-copy inner, serialize inner classTag + state).

The commit cycle:

| Call | What `LogStrain` does |
|---|---|
| `setTrialF(F)` | run Box 14.3 (i),(ii); call `inner->setTrialStrain(εᵉ_tr)` (Box 14.4 runs *inside* the unchanged inner) |
| `getStress()` | $\boldsymbol\sigma = J^{-1}\boldsymbol\tau$ (Cauchy, from the inner's Kirchhoff $\boldsymbol\tau$ at the J-seam) |
| `getTangent()` | the constitutive spatial part $\mathbf c$ (§5) |
| `commitState()` | $\mathbf b^e_n \leftarrow \exp[2\,\boldsymbol\varepsilon^e_{n+1}]$, $\mathbf F_n\leftarrow\mathbf F_{n+1}$; then `inner->commitState()` |
| `revertToLastCommit()` | restore $\mathbf b^e_n,\mathbf F_n$; then `inner->revertToLastCommit()` |
| `sendSelf` / `recvSelf` | round-trip the inner pair **plus** committed $\mathbf b^e_n$ per GP (the parallel/database path must carry the wrapper's own history) |

> [!warning] The J-seam is centralized (decision D2 — deferred but architected)
> **All** $1/J$ scaling happens at a single `toInner` / `fromInner` seam — never
> scattered across the stress and tangent paths. The default feeds the inner law
> the **Kirchhoff** stress $\boldsymbol\tau$ (for metals / moderate strain $J\approx
> 1$, this is a non-issue). Because $\mathbf F$ is stored as wrapper state, $J=\det
> \mathbf F$ is always in hand, so a future `-cauchyYield` flag — divide by $J$
> before the inner yield check, for pressure-sensitive soils at large volumetric
> $J$ — is roughly **10 lines**, not a rebuild. See §8 (YELLOW) and §12.

---

## 7. Repeated-eigenvalue degeneracy

> [!warning] The #1 trap — and it strikes the *simplest* states
> The tensor-log derivative $\mathbf L$ (14.101) and the spectral $\ln\mathbf b^e$
> contain off-diagonal terms
> $$\frac{\ln b_i - \ln b_j}{b_i - b_j},\qquad i\ne j,$$
> which are $0/0$ when principal stretches coincide ($b_i\to b_j$). The analytic
> limit is $1/b_i$, but naïve code divides by $\sim 0$ and returns **NaN**. The
> cruel part: degeneracy hits exactly the *simplest, most common* deformation
> states — **uniaxial** (2-fold repeat), **hydrostatic / equibiaxial** (3-fold),
> and any **axisymmetric** strain. It is load-bearing from the very first run.

**Ownership and fix.** Degeneracy handling lives **inside the wrapper** (it is part
of the seam-3 `MATISU` contract), and it ships from **day one**. The fix branches
on eigenvalue multiplicity — all-distinct / two-equal / three-equal — and uses the
**analytic Taylor limit** in each branch (Miehe 1998; dSNPO Appendix A.5,
subroutines `ISO2` / `DISO2`). It does **not** use $\varepsilon$-perturbation,
which is noisier near the threshold. Because the principal-stress form (§4) bounds
the problem to **three** eigenvalues, this is a closed, finite set of branches.

> [!check] The contract and its external net
> The wrapper **guarantees finite, non-NaN $\boldsymbol\sigma$ + tangent under
> repeated stretches**. The element's degeneracy battery (A2/A2′:
> hydrostatic/equibiaxial) is the external net for this guarantee, and the numpy
> oracle's three eigenvalue branches (`logstrain_reference.py`, with an analytic
> tangent vs finite-difference probe) is the contract the C++ must match.

---

## 8. The material matrix — what lifts cleanly

The binding requirement is structural: the inner material must expose a genuine
**3D, order-6, symmetric-tensor $\boldsymbol\sigma(\boldsymbol\varepsilon)$** law
(`getType()=="ThreeDimensional"`, `getOrder()==6`). The wrapper feeds it the
6-vector Hencky strain (engineering shear) and reads the returned 6-vector as
**Kirchhoff** stress.

### 🟢 GREEN — lifts cleanly (v1 targets)

| Inner material | Note |
|---|---|
| `ElasticIsotropicThreeDimensional` | becomes the **Hencky hyperelastic** law — a real finite-strain elastic model, well-behaved in compression (not StVK's spurious softening) |
| `J2ThreeDimensional`, [[LadrunoJ2_guide\|LadrunoJ2]] | **the ideal target** — isochoric plastic flow ⇒ exp-map preserves $\det\mathbf F^p=1$ *exactly*; the canonical demo |
| `DruckerPrager3D` | works (pressure-sensitive caveat ↓) |
| `ASDPlasticMaterial3D` family | 3D-only by design; cleanest inner laws |
| 3D soils — `PressureDependMultiYield`/`02`/`03`, `PressureIndependMultiYield`, `ManzariDafalias3D`, `SAniSandMS3D`, `BoundingCamClay3D`, … | 3D cyclic/soil laws (pressure-sensitive caveat ↓) |
| `ASDConcrete3DMaterial`, `PlasticDamageConcrete3d` | 3D damage / plastic-damage (regularization caveat ↓) |

### 🟡 YELLOW — wraps mechanically, but the *physics* needs a decision

- **Pressure-sensitive plasticity** (Drucker–Prager, the soils, Cam-Clay): the
  yield surface is evaluated on what the inner law *thinks* is Cauchy pressure, but
  the wrapper hands it **Kirchhoff** $\boldsymbol\tau = J\boldsymbol\sigma$. At
  large volumetric strain $J$ departs from 1, so the confining pressure the model
  sees is off by $J$. This is a *modeling choice* (and exactly the framework de
  Souza Neto & Perić 1996 used for soils) — resolve it per-material via the latent
  `-cauchyYield` J-seam flag (§6, D2).
- **Damage / softening** (ASDConcrete3D, PlasticDamageConcrete3d):
  characteristic-length / fracture-energy regularization is keyed to the
  *small-strain* measure; the IMPL-EX tangent must survive the log-map transform.
- **Rate / viscoplastic** (J2 with $\eta$, viscoplastic soils): the rate measure
  shifts $\dot\varepsilon\to$ Hencky rate (dSNPO §14.8 — format-identical).

### 🔴 RED — will **not** lift from the material side

| Group | Why |
|---|---|
| Plane-stress-only (`FSAM`, MCFT/RC panels, `PlaneStressUserMaterial`, …) | plane stress + finite strain must *solve* the out-of-plane stretch to keep $S_{33}=0$ (dSNPO §14.7 — separate plan) |
| Fiber/condensation adaptors (`BeamFiber*`, `PlateFiber*`, `J2PlateFibre`, …) | static-condensation wrappers for beam/shell elements that don't deliver $\mathbf F$ |
| Non-continuum (`AcousticMedium`, `ContactMaterial3D`) | not a strain-driven $\boldsymbol\sigma(\boldsymbol\varepsilon)$ at all |
| Axisymmetric (`ElasticIsotropicAxiSymm`, `J2AxiSymm`, …) | constrained-3D with a hoop term; needs an axisymmetric-aware $\mathbf F$ |
| `*Thermal` variants | temperature-coupled; out of scope |

### ⚪ Pass-through wrappers — compose *around* the inner law

`InitStrainNDMaterial`, `InitStressNDMaterial`, `MinMaxNDMaterial`,
`Parallel3DMaterial`, `Series3DMaterial`, `OrthotropicMaterial`,
`InitialStateAnalysisWrapper` — these wrap the *inner* 3D law; `LogStrain` lifts
the composed result.

> [!success] Capability statement
> `LogStrain` turns **any GREEN material** into a finite-strain model supporting
> arbitrary rotations and (for J2) arbitrary strains, with multiplicative
> plasticity and **exact** plastic incompressibility — *with zero change to the
> material's return map*. Elastic-wrapped becomes Hencky hyperelasticity. The one
> intrinsic assumption is **moderate elastic strains** (total strain and rotation
> are unrestricted) — automatically satisfied by metals and soils.

---

## 9. OpenSees implementation

### 9.1 Class and files

| Item | Value |
|---|---|
| Class | `LogStrainNDMaterial : public FiniteStrainNDMaterial` (concrete subclass) |
| Base | `FiniteStrainNDMaterial : public NDMaterial` (fork-local; adds `setTrialF(F)` + Cauchy/spatial-tangent accessors) |
| classTag | **`ND_TAG_LogStrainNDMaterial = 33010`** (`SRC/classTags.h`) |
| Source | `SRC/material/nD/LogStrainNDMaterial.{h,cpp}`, `SRC/material/nD/FiniteStrainNDMaterial.{h,cpp}` |
| Spectral / log-map math | `LogStrainKernel`-style (Jacobi 3×3 spectral, A.52/A.53 degeneracy-safe tensor-log derivative) |
| Parser | `nDMaterial LogStrain` in the NDMaterial command table |
| Broker | `FEM_ObjectBroker` case for parallel/database |
| Reference patterns | `MinMaxNDMaterial`, `InitStrainNDMaterial` (NDMaterial-wraps-NDMaterial), `J2ThreeDimensional` (inner contract + Voigt order), `OrthotropicMaterial` (3D tensor rotation) |

> [!note] Architecture — the `setTrialF` interface (decision D4)
> `setTrialF(F)` lives on the **fork-local** `FiniteStrainNDMaterial` base, **not**
> on vanilla `NDMaterial.h`. The element holds that base type and calls `setTrialF`
> through it (via a `dynamic_cast`) — so **inner materials are untouched and
> vanilla `NDMaterial.h` is never edited.** The only vanilla touches are the
> unavoidable, strictly-additive `// Ladruno` registration lines: the `classTags.h`
> `#define`, the `FEM_ObjectBroker` case, the NDMaterial command-table entry, and
> `CMakeLists.txt`. Promoting `setTrialF` onto base `NDMaterial` would be a
> separate, deliberately-surfaced decision; not now. Sibling members of the base:
> `LadrunoJ2Finite` (33012, see §12) and `InitDefGrad` (33013).

### 9.2 The seam-3 contract (element ↔ material)

The whole stack hinges on one interface call at the element↔material boundary:

```cpp
// LadrunoBrick (element, per Gauss point, finite path) calls:
mat->setTrialF(F);            // FULL deformation gradient F (3×3) in
sigma = mat->getStress();     // Cauchy σ = J⁻¹ τ   (6, Voigt) → ∫ Bᵀσ dv
c     = mat->getTangent();    // SPATIAL CONSTITUTIVE tangent c = (1/2J)[D:L:B]
```

The contract clauses (each a deliberate decision):

1. The element passes the **full** $\mathbf F$, not $\mathbf F_\Delta$; the adaptor
   derives $\mathbf F_\Delta = \mathbf F\mathbf F_n^{-1}$ itself and holds $J=\det
   \mathbf F$ — so it owns the Kirchhoff↔Cauchy conversion at a single seam.
2. `getStress()` returns **Cauchy** $\boldsymbol\sigma$, so the element assembles
   $\int\mathbf B^{\mathsf T}\boldsymbol\sigma\,\mathrm dv$ spatially with no
   pull-back, treating the adaptor as an ordinary `NDMaterial`.
3. `getTangent()` returns the **constitutive** part $\mathbf c$ only; the geometric
   term is the element's (§5).
4. The adaptor **owns committed $\mathbf b^e_n$ per GP** and wraps the inner
   material's commit/revert/serialize (§6).
5. The adaptor **guarantees finite, non-NaN $\boldsymbol\sigma$ + tangent under
   repeated stretches** (degeneracy, §7).

---

## 10. Command syntax and usage

### 10.1 Grammar

```tcl
nDMaterial LogStrain $tag $innerTag
```

`$innerTag` must reference an already-defined GREEN 3D `nDMaterial` (§8). The
wrapper deep-copies the inner law, so the inner tag may be reused or redefined
afterward without affecting the wrapped instance.

```python
ops.nDMaterial('LogStrain', tag, innerTag)
```

### 10.2 The canonical large-strain elastoplastic stack

```tcl
# Leg 3 — the small-strain law (ISOTROPIC ⇒ objective finite-strain lift, see §12)
#   Simo necking law: sig_y = 450 + 265(1 - e^{-16.93 p}) + 129.24 p
nDMaterial LadrunoJ2 1 $K $G -iso voce 450.0 265.0 16.93 129.24

# Leg 2 — lift it to finite strain (this wrapper)
nDMaterial LogStrain 2 1

# Leg 1 — an F-delivering element; bbar ⇒ F-bar (necking is isochoric-plastic → std locks)
element LadrunoBrick 101 $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8 2 -formulation bbar -geom finite

# F-bar tangent is UNSYMMETRIC → unsymmetric solver
system UmfPack
```

```python
import openseespy.opensees as ops
E, nu = 206900.0, 0.29
K, G = E/(3*(1-2*nu)), E/(2*(1+nu))

ops.nDMaterial('LadrunoJ2', 1, K, G, '-iso', 'voce', 450.0, 265.0, 16.93, 129.24)
ops.nDMaterial('LogStrain', 2, 1)                       # lift mat 1 to finite strain
ops.element('LadrunoBrick', 101, *nodes, 2,
            '-formulation', 'bbar', '-geom', 'finite')
ops.system('UmfPack')                                   # F-bar ⇒ unsymmetric
```

```tcl
# Large-strain ELASTIC (Hencky hyperelasticity) — rubber/elastomer block
nDMaterial ElasticIsotropic 1 $E $nu          ;# 3D order-6 inner
nDMaterial LogStrain 2 1
element LadrunoBrick 101 $n1 ... $n8 2 -formulation bbar -geom finite
```

> [!tip] Solver guidance (brief — full stack in the trifecta guide)
> - **Implicit bending into plasticity** under `-geom finite` does not converge
>   under plain `Newton`; use **`KrylovNewton` + `test EnergyIncr`** and a refined
>   (≥2×2) cross-section.
> - **`bbar+finite` (F-bar) has an unsymmetric tangent** — use an **unsymmetric
>   solver** (`UmfPack`, `FullGeneral`, `SparseGEN`); a symmetric solver silently
>   drops the coupling.
> See [[finite_strain_trifecta_guide]] §10 for the complete operational rules.

---

## 11. Verification and validation

### 11.1 The independent oracle (the correctness contract)

The C++ is pinned against an **independent numpy oracle**,
`tests/logstrain_reference.py`, which implements the spectral decomposition, the
**three eigenvalue branches** (the degeneracy net), and the spatial tangent
checked against a finite-difference probe (dSNPO 14.95). Because the oracle is a
separate implementation, an indexing or branch bug in the C++ cannot be mirrored
into the reference.

### 11.2 Single-point, prescribed-`F` battery (ships with the material)

Driving `setTrialF(F)` directly — no element needed:

| Test | What it pins |
|---|---|
| **T0 — elastic round-trip** | wrap `ElasticIsotropic3D`; pure rotation ⇒ zero stress (frame indifference); uniaxial stretch $\lambda$ ⇒ Hencky-hyperelastic closed form |
| **T1 — J2 uniaxial finite stretch** | wrap `J2ThreeDimensional`/`LadrunoJ2`; $\sigma(\text{log strain})$ vs the 1-D finite-J2 closed form; **exact** $\det\mathbf F^p=1$ |
| **T3 — simple-shear negative control** | large simple shear must **not** show the spurious Jaumann oscillation (dSNPO §14.10.3) — a direct discriminator of family (c) vs (b) |
| **TD — degeneracy battery** | uniaxial (2-fold), hydrostatic (3-fold), equibiaxial, axisymmetric → finite, non-NaN $\boldsymbol\sigma$ + tangent; analytic tangent vs FD agrees (§7) |

### 11.3 Element-level (joint — through `LadrunoBrick -geom finite`)

Necking of a circular bar (the canonical finite-strain J2 benchmark),
volumetric-locking check (near-incompressible J2 block, `std` vs F-bar), and the
finite-stretch patch test run through the element. These are recorded in the
end-to-end validation report; see [[finite_strain_trifecta_guide]] §11 (41 Zone-A
tests, all PASS, with the `logstrain_reference.py` and `ladrunoj2_reference.py`
oracles as anchors).

> [!note] Implementation milestones
> Seam-3 contract shipped in **PR #68** (the `FiniteStrainNDMaterial` base);
> the numpy oracle + elastic `LogStrainNDMaterial` v1 followed; the first GREEN
> **plastic** inner ([[LadrunoJ2_guide|LadrunoJ2]], via the extracted
> header-only `LadrunoJ2Kernel.h`) wired in **PR #81/#97** — `nDMaterial LogStrain
> $t $j2` with **zero change to `LogStrainNDMaterial`**. Full log:
> [[09_finite_strain_material_wrapper]].

---

## 12. Limitations and boundaries

> [!caution] Known boundaries
> - **Combined (kinematic) hardening is not objective under large rotation** — the
>   §14.11 boundary. The lift is **exact for the isotropic spine** (the elastic
>   state co-rotates and isotropic yield depends only on $\lVert\mathbf
>   s\rVert,\bar\varepsilon^p$, so a superposed rotation rotates Cauchy stress
>   rigidly — verified to ~$10^{-9}$). But an inner **backstress**
>   $\boldsymbol\alpha$ is stored in the inner's *fixed* frame and does **not**
>   co-rotate, so $\lVert\mathbf s-\boldsymbol\alpha\rVert$ is not
>   rotation-invariant. Pinned as a strict `xfail`. **Practical rule: under `-geom
>   finite`, use isotropic hardening (`-kin 0`) for an exact lift.** The fix is the
>   finite-strain-**native** material `LadrunoJ2Finite` (co-rotates
>   $\boldsymbol\alpha$ each step) — see [[16_finite_native_j2_adr]].
> - **Pressure-sensitive yield at large volumetric strain** — the wrapper feeds the
>   inner law **Kirchhoff** $\boldsymbol\tau$ ($J\approx1$ for metals/moderate
>   strain). Soils needing the yield read on true Cauchy pressure want the deferred
>   `-cauchyYield` J-seam flag (§6, D2).
> - **The inner must be a genuine 3D order-6 $\boldsymbol\sigma(\boldsymbol
>   \varepsilon)$.** Plane-stress-only, fiber/condensation, non-continuum,
>   axisymmetric, and `*Thermal` materials do **not** lift from the material side
>   (§8 RED).
> - **F-bar (`bbar+finite`) yields an unsymmetric element tangent** — needs an
>   unsymmetric solver throughout near-incompressible / locking studies (an element
>   concern, surfaced here because it gates the canonical stack).
> - The intrinsic continuum assumption is **moderate elastic strains** (total
>   strain and rotation are unrestricted) — satisfied by metals and soils.

---

## 13. References

**Theory**
- de Souza Neto, Perić & Owen (2008), *Computational Methods for Plasticity:
  Theory and Applications* (Wiley) — **Ch. 14** (finite-strain multiplicative
  plasticity: §14.3 log strain, §14.4 exponential map, Box 14.3 `MATISU`, §14.5
  spatial tangent eq. 14.99, §14.6 principal-stress form, §14.7 plane stress,
  §14.10 hypoelastic/Jaumann, §14.11 kinematic-hardening boundary), **Ch. 15**
  (F-bar §15.1), and **App. A.5** (`ISO2`/`DISO2`, repeated-eigenvalue treatment).
  Hencky hyperelasticity §13.2.3. The primary framework reference.
- Lee (1969) — the multiplicative split $\mathbf F=\mathbf F^e\mathbf F^p$.
- Miehe (1998) — closed-form derivative of isotropic tensor functions; the
  repeated-eigenvalue treatment for the log-map tangent.
- Eterovic & Bathe (1990); Weber & Anand (1990); Perić & Owen (1991); Simo & Miehe
  (1992); Cuitiño & Ortiz (1992) — log-strain / exponential-map return mapping.
- de Souza Neto & Perić (1996) — the framework applied to soil plasticity.
- Bonet & Wood, *Nonlinear Continuum Mechanics for FEA* — the updated-Lagrangian
  material/geometric tangent split.
- Abaqus Theory Manual — the finite-strain plasticity that uses exactly this
  technique.

**Within this repo**
- [[finite_strain_trifecta_guide]] — the **whole-stack entry point** (element +
  wrapper + material, with the end-to-end validation summary).
- [[09_finite_strain_material_wrapper]] — the chronological design/decision log
  (D1–D4, the full material matrix, the implementation log).
- [[LadrunoJ2_guide]] — the flagship GREEN inner law (combined-hardening J2).
- [[16_finite_native_j2_adr]] — `LadrunoJ2Finite`, the finite-strain-native fix for
  the §14.11 combined-hardening objectivity boundary.
- [[LadrunoBrick_reference]] — the F-delivering element (`-geom finite`, geometry
  seam §8).
- [[solid_transformation_wrapper]] — the element geometry layer (linear/corot/finite).
- [[Ladruno_materials_guide]] — the material catalog (§4.1 seeds this guide).
- Source: `SRC/material/nD/{LogStrainNDMaterial,FiniteStrainNDMaterial}.{h,cpp}`,
  spectral/log-map kernel; oracle `tests/logstrain_reference.py`. classTag **33010**.
```
