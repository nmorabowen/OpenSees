---
title: Finite-strain material wrapper (logarithmic / Hencky strain-space)
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - implementation
  - material
  - element
  - finite-strain
---

# Finite-strain material wrapper — logarithmic (Hencky) strain-space

> **Primary source.** de Souza Neto, Perić & Owen, *Computational Methods for
> Plasticity: Theory and Applications* (Wiley, 2008), **Ch. 14** (finite-strain
> elastoplasticity) and **Ch. 15** (large-strain incompressible elements). Local
> copy: `C:\Users\nmora\Desktop\FEM expert\Computational Methods for Plasticity (Neto, Peric and Owen 2008).pdf`.
> **Page offset: PDF page = book page + 25.** Equation numbers below are the
> book's; page cites are given as `book p. N (PDF M)`.

## What

A **material-side adaptor** that lifts an *unchanged* small-strain 3D
`NDMaterial` to a genuine finite-strain (large rotation **and** large strain)
material, by the logarithmic (Hencky) strain-space technique — option **(c)** of
the three wrapper families (StVK / hypoelastic-Jaumann / log-strain). The
material's return-mapping algorithm is reused **verbatim**; the wrapper supplies
the kinematic pre-/post-processing around it (Box 14.3, *MATISU*).

Because every one of these wrappers needs the **deformation gradient F**, not a
strain vector — and *no* OpenSees solid element currently delivers F (audit
below) — this is a **two-axis co-design split across two teams** that meet at
**seam 3** (the constitutive seam):

1. **Material axis — THIS PLAN (seam 3).** A new fork-local class
   `LogStrainNDMaterial` (the Hencky adaptor) plus a thin abstract base
   `FiniteStrainNDMaterial` that adds a `setTrialF(F)` entry point. The adaptor
   owns the per-GP hidden state `bᵉ` and the spectral/log-map machinery; it is a
   pure spatial map `F → (Cauchy σ, spatial tangent a)` wrapping an **unchanged**
   small-strain `NDMaterial`.
2. **Element axis — COMPANION PLAN (not this one).** Owned by the
   *solid-transformation-wrapper* effort: `SolidTransformation::finite` on
   `LadrunoBrick` (classTag 33002) computes **F / F_Δ** per GP, calls seam 3, and
   assembles the spatial internal force `∫ bᵀσ dv` + tangent (with F-bar, Ch. 15
   §15.1). See `Ladruno_implementation/solid_transformation_wrapper.md` and
   [[project_solid_transformation_wrapper]]. **We do NOT build a separate
   element** — that would duplicate `LadrunoBrick + finite`.

> **Cross-team reconciliation (settled 2026-06-01).** The `finite` route is
> **SPATIAL multiplicative (Box 14.3 *MATISU*), NOT total-Lagrangian** — Green–
> Lagrange/PK2 (TL) and the spatial log-strain adaptor don't compose. This
> *supersedes* the earlier D1=TL decision below. The adaptor owns committed `bᵉ`
> per GP; repeated-eigenvalue degeneracy is **our** responsibility (D3).

**In scope (v1):** isotropic, isothermal, 3D. The seam-3 adaptor
(`LogStrainNDMaterial` + `FiniteStrainNDMaterial` base) around any *GREEN*
material, its single-point (prescribed-`F`) verification battery, Tcl + Python
`nDMaterial LogStrain` commands. **Phase 0 (do first): lock the seam-3 signature**
— `setTrialF`/state-ownership contract — with the companion team before v3 element
wiring. Element-level benchmarks (necking, locking, patch) are *joint*, run
through `LadrunoBrick + finite`.

**Not in scope (v1, → follow-ups):** plane-stress finite plasticity (needs the
§14.7 nested-iteration path — separate plan), axisymmetric elements, anisotropic
/ single-crystal finite plasticity (Ch. 16), thermomechanical coupling,
kinematic hardening at finite strain (§14.11 — works but is a v2 once the
isotropic spine is proven), and the hypoelastic/Jaumann engine (option **b**,
§14.10 — only worth it as the corotational small-elastic-strain path).

## Why

The fork already has nonlinear *elements* (Bézier simplices) and explicit
dynamics, but every continuum material is small-strain: large-deformation soil,
metal forming, necking, rubber-like compression, and finite-rotation column
buckling with yielding are all out of reach. The log-strain wrapper is the
*highest-leverage* way to close that gap because it **reuses the entire existing
3D material library** (J2, Drucker–Prager, the soils, ASDPlastic, concrete-3D)
without touching a single return map. Abaqus does exactly this for its
finite-strain plasticity; it is the standard, validated technique
(`book p. 574 (PDF 599)`, §14.1, which cites Eterovic & Bathe 1990, Weber &
Anand 1990, Perić & Owen 1991, Simo 1992, Simo & Miehe 1992, Cuitiño & Ortiz
1992 — and notes it works for soil plasticity, de Souza Neto & Perić 1996, and
damage).

---

## What this can do — the material matrix

The binding requirement is structural: the inner material must expose a genuine
**3D, order-6, symmetric-tensor σ(ε)** law (`getType()=="ThreeDimensional"`,
`getOrder()==6`). The wrapper feeds it the 6-vector Hencky strain in the
canonical OpenSees order `{ε₀₀,ε₁₁,ε₂₂, 2ε₀₁,2ε₁₂,2ε₂₀}` (engineering shear) and
interprets the returned 6-vector as **Kirchhoff** stress (the work-conjugate to
Hencky strain — see Theory). Audit of `SRC/material/nD/` ⇒:

### 🟢 GREEN — lifts cleanly (v1 targets)

| Material | File | Note |
|---|---|---|
| `ElasticIsotropicThreeDimensional` | `ElasticIsotropicThreeDimensional.cpp` | becomes the **Hencky hyperelastic** law (§13.2.3, `book p. 528`) — a real finite-strain elastic model, well-behaved in compression, *not* the spurious-softening StVK of option (a) |
| `J2ThreeDimensional`, `SimplifiedJ2` | `J2ThreeDimensional.cpp`, `SimplifiedJ2.cpp` | **the ideal target** — isochoric plastic flow ⇒ exp-map preserves `det Fᵖ=1` *exactly* (14.65, `book p. 589`); the canonical Miehe/dSNPO demo |
| `DruckerPrager3D` | `UWmaterials/DruckerPrager3D.cpp` | works; pressure-sensitive caveat ↓ |
| `ASDPlasticMaterial3D` family | `ASDPlasticMaterial3D/` | 3D-only by design; cleanest inner laws |
| `ManzariDafalias3D(/RO)`, `SAniSandMS3D`, `BoundingCamClay3D`, `J2CyclicBoundingSurface3D`, `CycLiqCP3D`, `CycLiqCPSP3D`, `MultiaxialCyclicPlasticity3D`, `StressDensityModel3D` | `UWmaterials/`, `UANDESmaterials/`, … | 3D soils/cyclic; pressure-sensitive caveat ↓ |
| `PressureDependMultiYield`/`02`/`03`, `PressureIndependMultiYield`, `MultiYieldSurfaceClay` | `soil/` | **3D via the `ndm==3` switch** (verified — `getType()` returns `"ThreeDimensional"`) |
| `ASDConcrete3DMaterial`, `PlasticDamageConcrete3d` | `ASDConcrete3DMaterial.cpp`, … | 3D damage/plastic-damage; regularization caveat ↓ |

### 🟡 YELLOW — wraps mechanically, but the *physics* needs a decision

- **Pressure-sensitive plasticity** (Drucker–Prager, the soils, Cam-Clay,
  StressDensity): the yield surface is evaluated on what the unchanged material
  *thinks* is Cauchy pressure, but the wrapper hands it **Kirchhoff** stress
  `τ = Jσ` (14.33). At large volumetric strain `J` departs from 1, so the
  confining pressure the model sees is off by `J`. This is a *modeling choice*,
  not a bug — and it is exactly the framework de Souza Neto & Perić (1996) used
  for soils. Decide per-material whether to (a) accept τ-pressure, or (b) divide
  by `J` before the yield check (needs a flag). Also: reference/confining-stress
  state and `InitialStateAnalysisWrapper` / `updateMaterialStage` staging is set
  in small-strain terms.
- **Damage / softening** (ASDConcrete3D, PlasticDamageConcrete3d): characteristic
  -length / fracture-energy regularization is keyed to the small-strain measure;
  the IMPL-EX tangent must survive the log-map tangent transform (§14.5).
- **Rate / viscosity** (J2 η, viscoplastic soils): rate measure shifts ε̇ → Hencky
  rate; §14.8 covers finite viscoplasticity, format-identical.

### 🔴 RED — will NOT play along from the material side

| Group | Why | Where the real path lives |
|---|---|---|
| Plane-stress-only (`FSAM`, `ConcreteS`, `PlaneStressUserMaterial`, `PlasticDamageConcretePlaneStress`, `ElasticIsotropicPlaneStress2D`, RC-panel/MCFT) | plane stress + finite strain must *solve* the out-of-plane stretch to keep `S₃₃=0` | §14.7 `book p. 601 (PDF 626)` — nested iteration / plane-stress-projected; **separate plan** |
| Fiber/condensation adaptors (`BeamFiber*`, `PlateFiber*`, `J2BeamFiber*`, `J2PlateFibre`, `ConcreteMcftNonLinear*`) | static-condensation wrappers for beam/shell elements that don't deliver F | beam/shell finite-strain kinematics, out of scope |
| Non-continuum (`AcousticMedium` order-1 pressure; `ContactMaterial3D` traction–separation) | not a strain-driven σ(ε) at all | n/a |
| Externs (`FeapMaterial`, `ExternalNDMaterial`, `WrapperNDMaterial`, UMAT) | type/order at runtime; a UMAT may *already* be finite-strain | n/a |
| Axisymmetric (`ElasticIsotropicAxiSymm`, `J2AxiSymm`, …) | constrained-3D with a hoop term; needs an axisymmetric-aware F | v2 element |
| `*Thermal` variants | temperature-coupled; out of scope | v2 |

### ⚪ Pass-through wrappers — wrap the *inner* 3D law, not these

`InitStrainNDMaterial`, `InitStressNDMaterial`, `MinMaxNDMaterial`,
`Parallel3DMaterial`, `Series3DMaterial`, `OrthotropicMaterial`,
`InitialStateAnalysisWrapper` — compose around them.

**Capability statement.** v1 turns any GREEN material into a finite-strain model
supporting arbitrary rotations and (for J2) arbitrary strains, with multiplicative
plasticity and *exact* plastic incompressibility, *with zero change to the
material's return map*. Elastic-wrapped becomes Hencky hyperelasticity. The one
intrinsic assumption is **moderate elastic strains** (total strain/rotation
unrestricted) — automatically satisfied by metals and soils (`§14.1`).

---

## Theory — the wrapper math (three notations)

Conventions: bold upper = 2nd-order tensor; `I,J,K,L` reference indices, `i,j,k,l`
spatial; Voigt order `{11,22,33,12,23,13}` with **engineering** shear `γ=2ε`.

### 1. Kinematics: multiplicative split & the Hencky strain

Multiplicative decomposition (Lee 1969), `book p. 578 (PDF 603)`, eq (14.20):

- **Direct:**  **F** = **Fᵉ Fᵖ**
- **Index:**   `F_iJ = Fᵉ_iK Fᵖ_KJ`

The **Eulerian logarithmic (Hencky) elastic strain**, `book p. 582 (PDF 607)`,
eq (14.30), built from the elastic left Cauchy–Green tensor **Bᵉ = Fᵉ Fᵉᵀ = (Vᵉ)²**:

- **Direct:**  **εᵉ** = ln **Vᵉ** = ½ ln **Bᵉ**
- **Spectral (index):**  `εᵉ = Σₐ ln(λᵃₑ) Nₐ⊗Nₐ`, where `Bᵉ = Σₐ (λᵃₑ)² Nₐ⊗Nₐ`

Volumetric/deviatoric split (14.31): `εᵉ_v = tr[εᵉ] = ln Jᵉ`, `Jᵉ = det Fᵉ`.
Key consequence (14.64–14.65, `book p. 588–589`): a **traceless flow vector
⇒ isochoric plastic flow, exactly**, just as in small strain. This is why J2 is
the clean case.

### 2. Constitutive law: Kirchhoff stress is conjugate to Hencky strain

The free energy is written on the log strain, `ψ(εᵉ,α)` (14.32). The stress that
falls out is the **Kirchhoff stress** (14.33, `book p. 583 (PDF 608)`) — *this is
the work-conjugacy fact that drives the whole design*:

- **Direct:**  **τ** = ρ̄ ∂ψ/∂**εᵉ**
- **Hencky elastic law (the reused small-strain D):**  **τ** = **Dᵉ : εᵉ**
  — *identical* matrix to the infinitesimal isotropic elasticity tensor.

Yield and flow are posed in Kirchhoff stress (14.35, 14.37): `Φ(τ,A)≤0`,
`D̃ᵖ = γ̇ ∂Ψ/∂τ`. Cauchy stress is recovered only at the end:
**σ = τ / det F**.

### 3. The integration algorithm (Box 14.3 *MATISU* + Box 14.4)

`book p. 596 (PDF 621)`. Given the incremental displacement `Δu` over a step:

```
(i)   Kinematics — update F:
        F_Δ      = I + ∇ₙ[Δu]                         (14.76)
        F_{n+1}  = F_Δ · F_n

(ii)  Elastic predictor (push the stored Hencky strain forward):
        Bᵉ_n          = exp[2 εᵉ_n]                    (14.93)
        Bᵉᵗʳ_{n+1}    = F_Δ · Bᵉ_n · F_Δᵀ              (14.94)
        εᵉᵗʳ_{n+1}    = ½ ln[Bᵉᵗʳ_{n+1}]               (spectral, 14.104)
        αᵗʳ_{n+1}     = α_n

(iii) GOTO Box 14.4 — the UNCHANGED small-strain return map:
        εᵉ_{n+1} = εᵉᵗʳ_{n+1} − Δγ ∂Ψ/∂τ|_{n+1}        (14.89/14.90)
        α_{n+1}  = α_n + Δγ H_{n+1}
        Φ(τ_{n+1}, A_{n+1}) = 0
      → identical functional format to the infinitesimal scheme (7.25).
      ** This is exactly OpenSees' existing setTrialStrain→getStress. **

(iv)  Cauchy stress for the internal force:
        σ_{n+1} = (1/det F_{n+1}) · τ_{n+1}             (Box 14.3 iv)
```

The **exponential-map** discretisation of the flow rule (14.72–14.73,
`book p. 591 (PDF 616)`) is what makes step (iii) format-identical to small
strain *and* keeps `det Fᵖ=1` exactly for traceless flow (a plain backward-Euler,
14.74, would not — it loses volumetric accuracy). The plastic update never needs
`Fᵖ` explicitly; we store `εᵉ_n` as the state (Box 14.3 stores the Hencky strain).

**Principal-stress form (§14.6, `book p. 599 (PDF 624)`)** — the form we will
actually code, because it sidesteps assembling 4th-order log-map tensors:
spectrally decompose `Bᵉᵗʳ = Σᵢ bᵢ eᵢ⊗eᵢ` (14.103), take principal Hencky strains
`εᵢᵗʳ = ½ ln bᵢ` (14.104), run the principal-space small-strain return map to get
`τᵢ`, reassemble `σ = (1/J) Σᵢ τᵢ eᵢ⊗eᵢ` (14.105).

### 4. The consistent spatial tangent (§14.5)

`book p. 597 (PDF 622)`, eq (14.99):

- `a_ijkl = (1/2J) [ D : L : B ]_ijkl − σ_il δ_jk`
  - **D** = ∂τ/∂εᵉᵗʳ (14.100) — *the small-strain consistent tangent, reused
    verbatim*: the only material-dependent term.
  - **L** = ∂ ln[Bᵉᵗʳ]/∂Bᵉᵗʳ (14.101) — derivative of the tensor log (the
    “intricate transformation”; needs **care at repeated eigenvalues**, Miehe
    1998 / dSNPO Appendix A.5, subroutine `DISO2`/`ISO2`).
  - **B** (14.102), `B_ijkl = δ_ik (Bᵉᵗʳ)_jl + δ_jk (Bᵉᵗʳ)_il`.
  - `−σ_il δ_jk` is the geometric (initial-stress) term.

> **The one-line summary (Remark 14.5):** with the Hencky elastic law, *all*
> small-strain integration **and** tangent subroutines are reused without
> modification; only the kinematic pre/post-processing of Box 14.3 (i),(ii),(iv)
> is new, and it is **material-independent**.

### Why not (a) or (b)

- **(a) StVK** — degenerate elastic-only case: feed Green–Lagrange **E** into
  `D`. `S = D:E`. Spurious softening in compression; useless for plasticity.
- **(b) Hypoelastic/Jaumann** — §14.10, `book p. 615 (PDF 640)`; documented
  shear-stress oscillation under large simple shear (§14.10.3). Only the
  corotational small-elastic-strain engine, not a true finite-strain solution.

---

## The element axis — owned by the companion (not this plan)

**Audit result:** *no* OpenSees solid element computes F. `Brick`, `BbarBrick`,
`Twenty_Node_Brick`, `FourNodeTetrahedron`, `TenNodeTetrahedron`, the quads — all
build a linear `B` and call `mat->setTrialStrain(ε)` with a 6- (3D) or 3-vector
(2D). The only large-deformation elements are beam/shell. `BbarBrick` is the
B-bar (mean-dilatation) precedent (`SRC/element/brick/BbarBrick.cpp`,
`computeBbar()` ~L1054–1147), the small-strain ancestor of the finite-strain
**F-bar**.

The element that closes this gap is **`LadrunoBrick` + `SolidTransformation::finite`**
(companion plan), *not* a new element here. That layer: computes `F` / `F_Δ` per
GP, calls **seam 3** (our adaptor), receives Cauchy σ + spatial tangent **a**, and
assembles the **spatial** internal force `∫ bᵀσ dv` + tangent. F-bar
(`F̄ = (J̄/J)^{1/3} F`, §15.1 `book p. 648 (PDF 673)`; eq. 15.5) is a *seam-1
kinematics-ledger* modification on their side — `bbar + finite = F-bar` (not an
oblivious wrapper); `eas + finite` is stability-fragile (§15.2.5). §15.1.9
(`book p. 665`) F-bar-Patch for simplex ties into BezierTri6/Tet10.

**Why F-bar matters (their concern, our awareness):** wrapped J2 and the soils
are (near-)isochoric, so a plain displacement hex locks volumetrically at finite
strain. Our single-point tests don't see locking (one GP); the *joint* element
tests (T4) do.

---

## Where (files)

- **New — material (THIS plan):**
  - `SRC/material/nD/FiniteStrainNDMaterial.{h,cpp}` — thin abstract base on
    `NDMaterial` adding `virtual int setTrialF(const Matrix &F)` + Cauchy/spatial-
    tangent accessors. **No vanilla `NDMaterial.h` edit** — new subclass.
  - `SRC/material/nD/LogStrainNDMaterial.{h,cpp}` — the Hencky adaptor (seam 3).
    Holds an inner `NDMaterial*` (any GREEN 3D law) + committed `bᵉ` per instance;
    implements Box 14.3 (i),(ii),(iv) + §14.6 + §14.5 + the degeneracy branch.
- **Element side:** none here — `LadrunoBrick + finite` is the companion's
  (`SRC/element/...`, classTag 33002).
- **Modify (additive `// Ladruno` hooks only):**
  - `SRC/classTags.h` — `ND_TAG_LogStrainNDMaterial 33010` (material tag only;
    the element tag 33002 belongs to the companion).
  - `SRC/api/FEM_ObjectBroker.cpp` — register the material for parallel/broker.
  - the NDMaterial command file (`OpenSeesNDMaterialCommands.cpp` / interpreter) —
    `nDMaterial LogStrain` factory registration.
  - `SRC/material/nD/CMakeLists.txt`.
  - `Ladruno_scripts/banner_features.txt` (+ `patch_banner.py`).
- **Reference (copy patterns):** `MinMaxNDMaterial.cpp` / `InitStrainNDMaterial.cpp`
  (the NDMaterial-wraps-NDMaterial pattern, incl. `getCopy`, `sendSelf`/`recvSelf`
  of a sub-material — the adaptor wraps the inner material's commit/revert/serialize),
  `J2ThreeDimensional.{h,cpp}` (inner-material interface + Voigt order),
  `OrthotropicMaterial.cpp` (3D tensor rotation bookkeeping).
- **Ledgers:** add a `draft` row to [[LEDGER_implementations]] (done); update
  [[LEDGER_vanilla_files]] when the additive hooks land; log gotchas in
  [[LEDGER_quirks]] (eigenvalue degeneracy in the log-map tangent).

## How — API & data flow

**User-facing (the "wrapper" framing made literal):**

```tcl
# 1. the unchanged small-strain 3D law
nDMaterial J2ThreeDimensional 1  $K $G $sig0 $sigInf $delta $H
# 2. lift it to finite strain (seam 3 — this plan)
nDMaterial LogStrain 2 1
# 3. an F-delivering element (companion plan: LadrunoBrick + finite geometry)
element LadrunoBrick 1  $n1 ... $n8  2  -geom finite     ;# companion syntax, TBD
```

```python
ops.nDMaterial('J2ThreeDimensional', 1, K, G, sig0, sigInf, delta, H)
ops.nDMaterial('LogStrain', 2, 1)                       # wrap mat 1 (seam 3)
# element via the companion's LadrunoBrick + finite
```

**Seam-3 contract & commit cycle (Phase 0 — lock this with the companion).** The
element layer calls `mat->setTrialF(F)` (or `setTrialFincr(F_Δ)`); the adaptor
does Box 14.3 (i),(ii), calls `inner->setTrialStrain(εᵉᵗʳ)` [Box 14.4 runs *inside
the unchanged material*], then exposes `getStress()` = **Cauchy σ** and
`getTangent()` = **spatial tangent a** (14.99). **State ownership:** the adaptor
owns committed `bᵉ_n` per GP (elastic left Cauchy–Green — *not* a plastic
multiplier, *not* total strain) and **wraps** the inner material's
`commitState`/`revertToLastCommit`/`sendSelf`/`recvSelf`. `commitState` sets
`bᵉ_n ← exp[2 εᵉ_{n+1}]` from the inner material's committed Hencky strain.

```cpp
// LogStrainNDMaterial::setTrialF  (sketch — principal-stress form, §14.6; SPATIAL)
int LogStrainNDMaterial::setTrialF(const Matrix &F) {
  Ftrial = F;                                   // (i)
  Matrix Be_tr = pushForward(Be_committed, F, Fcommitted);   // FΔ bᵉ FΔᵀ  (14.94)
  // (ii) spectral decomposition of Be_tr  → eigvals b_i, eigvecs e_i
  spectral3x3_degenSafe(Be_tr, bvals, evecs);   // branch on repeated eigenvalues (D3)
  for (int a=0;a<3;a++) eps_tr(a) = 0.5*log(bvals[a]);   // (14.104)
  // (iii) hand Hencky strain to the UNCHANGED small-strain material as "total strain"
  inner->setTrialStrain( assembleVoigt(eps_tr, evecs) );  // 6-vector, eng. shear
  return 0;
}
const Vector& LogStrainNDMaterial::getStress() {           // Cauchy (spatial)
  const Vector &tau = fromInner( inner->getStress() );     // τ (Kirchhoff); J-seam (D2)
  J = det(Ftrial);
  sigma = (1.0/J) * tau;                          // Box 14.3 (iv)
  return sigma;
}
```

## Testing — analytical references

**Single-point, prescribed-`F` (OURS — no element needed; the first deliverable).**
Drive `LogStrainNDMaterial::setTrialF(F)` directly from a unit harness:

- **T0 elastic round-trip.** Wrap `ElasticIsotropicThreeDimensional`; pure
  rotation ⇒ zero stress (frame indifference); uniaxial stretch λ ⇒ Hencky
  hyperelastic closed form `τ = D : (ln λ …)`. Compare to §13.2.3 / §13.4.3.
- **T1 J2 uniaxial finite stretch.** Wrap `J2ThreeDimensional`; monotonic
  uniaxial; check `σ(logarithmic strain)` against the 1-D finite J2 closed form
  (§14.2, `book p. 575`) and exact volume preservation `det Fᵖ=1`.
- **T3 simple-shear sanity (negative control).** Large simple shear with the
  log-strain adaptor must **not** show the spurious Jaumann stress oscillation of
  §14.10.3 — a direct discriminator of (c) vs (b).
- **TD degeneracy battery (D3, load-bearing).** Uniaxial (2-fold), hydrostatic
  (3-fold), equibiaxial, axisymmetric → adaptor returns **finite, non-NaN** σ +
  tangent; analytic-tangent vs finite-difference probe agrees. *This is the net
  for our own degeneracy bug and the companion's contractual guarantee.*

**Element-level (JOINT — run through the companion `LadrunoBrick + finite`):**

- **T2 necking of a circular bar** (§14.9.2, `book p. 607`) — canonical
  finite-strain J2 benchmark; neck radius / load-displacement.
- **T4 volumetric-locking check.** Near-incompressible J2 block, `finite` *without*
  vs *with* F-bar; F-bar must remove pressure checkerboarding (§15.1.6).
- **T5 patch test** (rigid-body + constant-F) at finite stretch.

Wire into the Zone-A pytest battery alongside `tests/` (numpy-based). The
single-point battery is C++/Python-harness only and ships with **this** PR's code;
element-level is gated on the companion element landing.

## Decisions (resolved) & open questions

### ✅ D1 — Formulation: **SPATIAL multiplicative (Box 14.3 MATISU)** — *revised 2026-06-01*

**Supersedes the earlier D1=TL.** Cross-team settlement: the `finite` route is
**spatial multiplicative**, not total-Lagrangian. Green–Lagrange **E** / PK2 **S**
(TL) and a spatial log-strain adaptor **don't compose** — the canonical, only-
self-consistent route is Box 14.3: `bᵉᵗʳ = F_Δ bᵉ_n F_Δᵀ → εᵉ = ½ ln bᵉ`
(Eulerian, 14.30) → **unchanged** small-strain return map → Kirchhoff τ → Cauchy
**σ = J⁻¹ τ**. This is where "small-strain material reused verbatim" is *literally*
true (14.89–90 = infinitesimal 7.25). The **exponential-map** integrator (14.72) is
the load-bearing ingredient (Remark 14.5; makes plastic incompressibility exact,
Remark 14.4).

> **Governing principle — the MATERIAL emits SPATIAL quantities.**
> `LogStrainNDMaterial` is a pure map `F → (Cauchy σ, spatial tangent a)`. The
> companion element assembles spatially `∫ bᵀσ dv` directly — no pull-back. (A
> hypothetical TL caller could still pull back to S/C, but that's not our route.)
> The log-strain path handles large rotation *and* strain intrinsically (spectral
> F is frame-indifferent) — it is the option-(c) *alternative* to a corotational
> (option-b) wrapper, **not** a client of it; do **not** stack `corot` on top of
> `finite` (double-counts rotation).

### ✅ D4 — Interface extension: **fork-local base, additive registration only** (resolved 2026-06-01)

`setTrialF(F)` lives on the fork-local `FiniteStrainNDMaterial` base; the companion
element holds that type and calls it via that pointer — **zero `NDMaterial.h`
edit**, inner materials untouched. The only vanilla touches are the unavoidable
registration points (`classTags.h` `#define`, `FEM_ObjectBroker` case, the
NDMaterial command table, `CMakeLists.txt`) — strictly-additive `// Ladruno`
lines, the same kind every shipped fork feature already adds. Promoting `setTrialF`
onto base `NDMaterial` would be a separate, deliberately-surfaced decision; not now.

### 🟡 D2 — Kirchhoff vs Cauchy pressure: **deferred, but architected** (open)

Don't resolve the physics now; make it a flag, not a rebuild. Three constraints
guarantee that: (1) **store F as wrapper state** ⇒ `J = det F` always in hand at
the inner call; (2) **centralize all J-scaling at one seam** (`toInner`/`fromInner`)
— never scatter `/J` across the stress and tangent paths; (3) **the tangent uses
that same seam**, so the convention can't drift. Default: feed the inner law
Kirchhoff τ (metals/moderate strain: J≈1, non-issue). Latent `-cauchyYield` flag
divides by J at the seam when a soil at large volumetric strain needs its
friction/cohesion read on Cauchy pressure. Deferral cost ≈ 10 lines.

### 🟡 D3 — Repeated-eigenvalue handling: **OURS, ship from day 1** (open, #1 trap)

**Ownership: this team** (it lives inside seam-3 `MATISU`). The log-map (14.101)
and spectral ln **B** contain off-diagonal terms `(ln bᵢ − ln bⱼ)/(bᵢ − bⱼ)`, i≠j.
When principal stretches coincide (`bᵢ→bⱼ`) this is 0/0; the limit is `1/bᵢ` but
naive code divides by ~0 → NaN. It strikes the *simplest* cases — **uniaxial**
(2-fold), **hydrostatic/equibiaxial** (3-fold), any axisymmetric strain — i.e.
exactly the T0/T1/TD set, load-bearing from the first run. Fix (Miehe 1998; dSNPO
App. A.5, `ISO2`/`DISO2`): branch on multiplicity (all-distinct / two-equal /
three-equal) with the analytic Taylor limit, **not** ε-perturbation (noisier near
the threshold). Bounded **3-eigenvalue** problem (principal-stress form §14.6).
Contract with the companion: our `update()` returns finite, non-NaN σ + tangent
under repeats; **their** finite T0/T1 battery is the external net for this bug.

### Remaining open

> [!question]
> **`getCopy`, `sendSelf/recvSelf` of the wrapped inner material.** Follow
> `MinMaxNDMaterial`/`InitStrainNDMaterial` exactly (deep-copy inner, serialize
> inner's classTag + state). Parallel (`OpenSeesMP`) must round-trip the pair —
> **plus** our committed `bᵉ_n` per GP.

> [!question]
> **Shared util is thinner than first thought (settled 2026-06-01).** The companion
> `corot` extracts R via **iterative Higham polar — not eigen-of-U** — so the two
> efforts do **not** share a polar/spectral kernel. The genuine shared asset shrinks
> to (a) the degeneracy-safe symmetric 3×3 eigensolver (ours, inside seam 3) and
> (b) the **degenerate-stretch test fixtures** — may end up being just the fixtures.

> [!question]
> **Phase 0 — lock the seam-3 signature** (`setTrialF`/`setTrialFincr`, who owns
> `bᵉ_n`, σ+a return convention) jointly with the companion **before** v3 element
> wiring. This is the first action item, not a parallel one.

- **Numerical:** conditioning of `exp[2εᵉ]`/`ln Bᵉ` at large strain; consistency
  of the analytic tangent vs a finite-difference probe (gate in CI).
- **Backwards compat:** purely additive — new classes + additive hooks; no
  existing model changes behavior.

## References

- **de Souza Neto, Perić, Owen (2008)**, *Computational Methods for Plasticity*,
  Wiley. **Ch. 14** (Boxes 14.2 `p.587`, 14.3 `p.596`, 14.4 `p.597`; §14.5 tangent
  `p.597`; §14.6 principal-stress `p.599`; §14.7 plane-stress `p.601`; §14.8
  viscoplastic; §14.10 hypoelastic/Jaumann `p.615`; §14.11 kinematic hardening
  `p.633`) and **Ch. 15** (F-bar §15.1 `p.648`, F-bar-Patch simplex §15.1.9
  `p.665`, EAS §15.2, u/p §15.3). Hencky hyperelasticity §13.2.3 `p.528`.
- Lee (1969); Lee & Liu (1967) — multiplicative split **F = Fᵉ Fᵖ**.
- Eterovic & Bathe (1990); Weber & Anand (1990); Perić & Owen (1991);
  Simo (1992); Simo & Miehe (1992); Cuitiño & Ortiz (1992) — log-strain /
  exponential-map return mapping.
- **Miehe (1998)** — closed-form derivative of isotropic tensor functions;
  repeated-eigenvalue treatment for the log-map tangent.
- de Souza Neto & Perić (1996) — the framework applied to soil plasticity.
- OpenSees: `SRC/element/brick/BbarBrick.cpp` (B-bar precedent),
  `SRC/material/nD/{J2ThreeDimensional,MinMaxNDMaterial,InitStrainNDMaterial,OrthotropicMaterial}.cpp`.

## Implementation log

- **2026-06-01 — seam-3 contract shipped (PR #68, merged).** `FiniteStrainNDMaterial.h`
  abstract base on `ladruno`: `setTrialF(F)`, Cauchy σ / spatial-tangent accessors,
  `setTrialStrain` disabled. Plan PR #67 merged.
- **2026-06-01 — numpy reference oracle (branch `guppi/logstrain-material`, 21/21).**
  `tests/logstrain_reference.py` + 2 test files. Box A.5 spectral, (A.51) iso
  function, (A.52/A.53) derivative with all three eigenvalue branches (D3 net vs
  FD), J2 radial-return + consistent tangent, spatial tangent a (14.99) vs the FD
  definition (14.95) — elastic + plastic, exact plastic incompressibility, simple-
  shear no-Jaumann-oscillation. The oracle the C++ must match.
- **2026-06-01 — C++ `LogStrainNDMaterial` (elastic v1) written + wired.**
  `SRC/material/nD/LogStrainNDMaterial.{h,cpp}`; Jacobi 3×3 spectral, A.52/A.53
  degeneracy-safe tensor-log derivative, spatial MATERIAL tangent `(1/2J)[D:L:B]`
  (geometric term = element K_geo — **contract clarification to propagate to the
  companion**: `getTangent()` is the symmetric material part, not the full
  non-minor-symmetric `a`). Owns committed `bᵉ_n`+`F_n`; wraps inner state.
  classTag `ND_TAG_LogStrainNDMaterial 33010`; `nDMaterial LogStrain $tag $inner`;
  CMake. **g++ 15.2 `-fsyntax-only` clean** (concrete class ⇒ all pure virtuals
  implemented). v1 correct for ELASTIC inner; plastic-inner state protocol = next
  decision (avoid εᵖ double-count when reusing a stateful return map).
- **NEXT:** full Intel build of the worktree, then the openseespy acceptance test
  (drive `nDMaterial LogStrain` over `ElasticIsotropicThreeDimensional` with
  prescribed `F`, assert match to the oracle ~1e-10), then open the PR. Then the
  plastic-inner protocol.

*(move to `Ladruno_internal/implemented_finite_strain_wrapper.md` when shipped)*
