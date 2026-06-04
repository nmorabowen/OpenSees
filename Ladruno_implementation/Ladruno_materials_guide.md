---
title: "Ladruno Materials Guide — the fork-authored constitutive library"
project: Ladruno
type: reference / user guide (catalog)
status: shipped (all listed materials on `ladruno`)
related:
  - "[[LadrunoJ2_guide]]"
  - "[[LadrunoUniaxialJ2_guide]]"
  - "[[LadrunoRebarBuckling_guide]]"
  - "[[finite_strain_trifecta_guide]]"
  - "[[09_finite_strain_material_wrapper]]"
  - "[[15_lemaitre_ductile_damage_adr]]"
  - "[[16_finite_native_j2_adr]]"
  - "[[staged_deformation_gradiend]]"
  - "[[20_ladruno_embedded_reinforcement_adr]]"
tags:
  - material
  - catalog
  - reference
  - guide
  - plasticity
  - finite-strain
  - steel
---

# Ladruno Materials Guide — the fork-authored constitutive library

> [!abstract] What this is
> The **single catalog** of every constitutive material the Ladruno fork adds on
> top of vanilla OpenSees — the plasticity core, the finite-strain & staged
> wrappers, and the steel/rebar overlays. For each: a theory capsule, the
> architectural decision behind it, the OpenSees command, and the intended use
> case. Materials with a dedicated deep-dive (`LadrunoJ2`, `LadrunoUniaxialJ2`,
> `LadrunoRebarBuckling`) are summarized here and linked; the rest are documented
> in full. The companion docs are the element catalog [[LadrunoBrick_reference]]
> and the large-deformation stack [[finite_strain_trifecta_guide]].

This guide is organized by **family**, because the fork's materials are not a flat
list — they share kernels, compose by nesting, and follow a small set of recurring
architectural patterns. Read §1–§2 first (the patterns); then jump to the family
you need.

---

## Contents

1. [[#1 · The material map|The material map]]
2. [[#2 · Architectural patterns (read once)|Architectural patterns]]
3. [[#3 · Family A — the J2 plasticity core|Family A: J2 plasticity core]]
4. [[#4 · Family B — finite-strain & staged wrappers|Family B: finite-strain & staged wrappers]]
5. [[#5 · Family C — steel & rebar overlays|Family C: steel & rebar overlays]]
6. [[#6 · The Lemaitre ductile-damage mode|Lemaitre damage]]
7. [[#7 · Composition — how the wrappers nest|Composition]]
8. [[#8 · Class-tag registry|Class-tag registry]]
9. [[#9 · Use-case decision guide|Use-case decision guide]]
10. [[#10 · References|References]]

---

## 1 · The material map

| Material | Kind | classTag | Role | Deep-dive |
|---|---|---|---|---|
| **`LadrunoJ2`** | nDMaterial | 33011 | combined-hardening (Voce+Chaboche) von Mises — the flagship 3D plastic law | [[LadrunoJ2_guide]] |
| **`LadrunoUniaxialJ2`** | uniaxialMaterial | 33000 | uniaxial twin (fiber/truss/zeroLength); the calibration oracle | [[LadrunoUniaxialJ2_guide]] |
| **`LadrunoJ2Finite`** | nDMaterial | 33012 | finite-strain-**native** combined J2 (co-rotating backstress) | [[16_finite_native_j2_adr]] |
| **`LogStrain`** | nDMaterial (wrapper) | 33010 | Hencky log-strain lift: any GREEN 3D law → finite strain | [[finite_strain_trifecta_guide]], [[09_finite_strain_material_wrapper]] |
| **`InitDefGrad`** (`StagedDefGrad`) | nDMaterial (wrapper) | 33013 | stress-free **finite** staged-construction birth (`F₀`) | [[staged_deformation_gradiend]] |
| **`StagedStrain`** | nDMaterial (wrapper) | 33014 | stress-free **small-strain** staged birth (`ε₀`) | [[staged_deformation_gradiend]] |
| **`LadrunoRebarBuckling`** | uniaxialMaterial (wrapper) | 33001 | reinforcing-bar buckling-average degradation overlay | [[LadrunoRebarBuckling_guide]] |
| **`LadrunoBondSlip`** | uniaxialMaterial | 33002 | 1D bond-slip τ–s (MC2010 backbone) for embedded rebar | [[20_ladruno_embedded_reinforcement_adr]] |

Plus the **Lemaitre ductile-damage mode** (`-damage lemaitre …`) — not a class, a
mode on both J2 materials (§6). And three **shared header-only kernels** that are
the reason the J2 family stays consistent (§2.1).

---

## 2 · Architectural patterns (read once)

Four patterns recur across the whole library. Knowing them makes every material
below predictable.

### 2.1 One verified kernel, many faces

The J2 family never re-implements its return map. The stress-update + tangent live
in **header-only, OpenSees-free kernels** over plain `double[]`:

| Header | Owns | Consumed by |
|---|---|---|
| `LadrunoJ2Kernel.h` | the von Mises return map + consistent tangent (scalar Newton on Δγ) | `LadrunoJ2`, all its dimensional views, `LadrunoJ2Finite` |
| `LadrunoHardening.h` | the Voce+linear yield law `σ_y(p)` (the **1e-12 oracle contract**) | `LadrunoJ2`, `LadrunoUniaxialJ2`, `LadrunoJ2Finite` |
| `LadrunoDamage.h` | the Lemaitre damage helpers (`Y`, `D` update) | `LadrunoJ2`, `LadrunoUniaxialJ2` (`-damage` mode) |

The payoff: the *same* verified code serves the small-strain material, every
dimensional view, the finite-strain path, and the uniaxial twin — there is never a
second implementation to drift. The 1D and 3D classes are even pinned to agree to
~$10^{-12}$ under uniaxial stress (that is `LadrunoUniaxialJ2`'s job as an oracle).

### 2.2 The wrapper / composition pattern

Many materials are **wrappers** — they hold an inner `*Material*` and modify its
input or output, copying the `MinMaxNDMaterial`/`InitStrainNDMaterial` idiom
(deep-copy inner, serialize inner classTag + state, delegate commit/revert). This
is how the fork gets combinatorial reach without per-material rewrites:

- `LogStrain` wraps a small-strain 3D law → lifts it to finite strain.
- `InitDefGrad`/`StagedStrain` wrap a material → give it a stress-free birth config.
- `LadrunoRebarBuckling` wraps a steel `UniaxialMaterial` → adds buckling degradation.

Wrappers **nest** (§7) — that is a deliberate design choice, not an accident.

### 2.3 The `FiniteStrainNDMaterial` base (the F-interface)

Finite-strain materials are driven by the **deformation gradient**, not a strain
vector. The fork adds one abstract base, `FiniteStrainNDMaterial`, with
`setTrialF(F)` + Cauchy-stress / spatial-tangent accessors — **without editing
vanilla `NDMaterial.h`**. `LadrunoBrick -geom finite` drives any
`FiniteStrainNDMaterial` through a single `dynamic_cast`. Members of this base:
`LogStrain` (33010), `LadrunoJ2Finite` (33012), `InitDefGrad` (33013). See
[[finite_strain_trifecta_guide]] §8 for the full seam contract.

### 2.4 The hygiene rules (non-negotiable, fork-wide)

Every fork material obeys rules that several upstream materials violate:

- **No `exit()` — ever.** A non-converged/singular solve returns a **status code**
  → an `opserr` warning + a step-cut request, never a libc `exit(-1)` that kills the
  Python kernel. (`SimplifiedJ2`, `InitStrain` on non-3D inners both `exit(-1)`.)
- **Real `getInitialTangent`** = the *elastic* tangent (correct for modal /
  initial-stiffness analyses). `SimplifiedJ2` returns the plastic tangent — a bug.
- **Real `revertToLastCommit`/`revertToStart`** — state genuinely rolls back, so
  adaptive step-cutting works. (`SimplifiedJ2`'s are empty stubs.)
- **Real `setParameter`/`Print`/`sendSelf`/`recvSelf`** — parallel- and
  database-safe, not stubs.

---

## 3 · Family A — the J2 plasticity core

The deviatoric (pressure-independent) von Mises plasticity family. All three share
`LadrunoJ2Kernel.h` + `LadrunoHardening.h`.

### 3.1 `LadrunoJ2` (nDMaterial, 33011) — the flagship

**What.** A single self-contained **rate-independent von Mises** `nDMaterial`
unifying nonlinear **isotropic** (Voce + linear) and nonlinear **kinematic**
(Chaboche / Armstrong–Frederick) hardening — the OpenSees analogue of Abaqus
`*PLASTIC, COMBINED`. Supersedes both upstream J2 materials (`J2Plasticity`: no
kinematic; `SimplifiedJ2`: linear-only + real defects).

**Theory capsule.** Yield on the backstress-shifted deviator
$\boldsymbol\xi=\mathbf s-\boldsymbol\alpha$:
$f=\lVert\boldsymbol\xi\rVert-\sqrt{2/3}\,\sigma_y(\bar\varepsilon^p)$; Voce law
$\sigma_y=\sigma_0+Q_\infty(1-e^{-b\bar\varepsilon^p})+H_{\text{iso}}\bar\varepsilon^p$;
Chaboche backstress
$\dot{\boldsymbol\alpha}_k=\tfrac23C_k\dot{\boldsymbol\varepsilon}^p-\gamma_k\boldsymbol\alpha_k\dot{\bar\varepsilon}^p$.
The backward-Euler return map collapses to a **single scalar Newton** on $\Delta\gamma$
(Kobayashi–Ohno: $\mathbf n=\mathbf M/\lVert\mathbf M\rVert$); the consistent tangent
is assembled exactly with the AF cross-term. **Parameter limits recover every
simpler model** — `-kin 0` ⇒ bit-identical to upstream `J2Plasticity`.

**Architecture.** Header-only kernel (§2.1); 5 dimensional views (3D / PlaneStrain /
AxiSymm / PlateFiber / PlaneStress) from one class via a `vmap` index table +
nested-Newton condensation for the plane-stress σ₂₂=0 constraint.

**OpenSees.**
```tcl
nDMaterial LadrunoJ2 $tag $K $G  -iso voce $sig0 $Qinf $b $Hiso  -kin $N $C1 $g1 $C2 $g2 …
```
**Use case.** Cyclic/seismic steel continuum (panel zones, plate buckling, BRBs,
shear links) where the Bauschinger effect and correct ratcheting rate matter;
monotonic ductile metal forming; the inner law of the finite-strain stack.
**Full reference: [[LadrunoJ2_guide]].**

### 3.2 `LadrunoUniaxialJ2` (uniaxialMaterial, 33000) — the uniaxial twin

**What.** The 1-D twin of `LadrunoJ2`: a rate-independent uniaxial J2
`UniaxialMaterial` with the same Voce+linear isotropic and Chaboche kinematic
hardening, consumable by **fiber sections, trusses, and zeroLength** elements.

**Theory capsule.** Because the 1D relative stress $\sigma_{\text{tr}}-X$ cannot
*rotate*, the return map is simpler than the 3D one: a **fixed direction**
$s=\operatorname{sign}(\sigma_{\text{tr}}-X_n)$ + a **scalar Newton on $\Delta\bar p$**;
analytic tangent $E_{\text{alg}}=Eh/(E+h)$,
$h=\sigma_y'+\sum_k a_k/(1+\gamma_k\Delta\bar p)^2$. Shares `LadrunoHardening.h`
verbatim with the 3D class — the **1e-12 oracle contract**.

**Architecture / intent.** Deliberately **not** a "better Steel02". Its three jobs:
(1) the **verification oracle** that pins the 3D Chaboche calibration under uniaxial
stress; (2) a fiber material with **true multi-backstress ratcheting** (what
Menegotto–Pinto cannot do); (3) the clean **core** for future steel layers
(fracture/LCF, local-buckling, weld/HAZ). Hygiene rules of §2.4 all satisfied (fixes
the `SimplifiedJ2` modal bug).

**OpenSees.**
```tcl
uniaxialMaterial LadrunoUniaxialJ2 $tag $E $sig0  -iso voce $Qinf $b $Hiso  -kin $N $C1 $g1 …
```
**Use case.** Steel fibers in fiber sections needing genuine ratcheting; truss/
zeroLength steel; the calibration oracle for a 3D `LadrunoJ2` model.
**Full reference: [[LadrunoUniaxialJ2_guide]].**

### 3.3 `LadrunoJ2Finite` (nDMaterial, 33012) — finite-strain-native

**What.** A `FiniteStrainNDMaterial` subclass that does combined-hardening J2 at
finite strain **natively** — it owns all the finite-strain state (the elastic left
Cauchy–Green $\mathbf b^e_n$, $\bar\varepsilon^p_n$, and the backstresses
$\boldsymbol\alpha_{n,k}$ **as spatial tensors**) and **co-rotates the backstress**
each step.

**Why it exists — the §14.11 boundary.** The wrapper path
(`LogStrain → LadrunoJ2`) is exact for the **isotropic spine** but *not*
frame-indifferent for kinematic hardening under large rotation: the backstress
$\boldsymbol\alpha$ lives in the inner material's *fixed* frame and the wrapper has
no channel to rotate it, so $\lVert\mathbf s-\boldsymbol\alpha\rVert$ is not
rotation-invariant (pinned as a strict `xfail`). This is the dSNPO §14.11
limit — fundamental to the "unchanged inner" architecture, not a wiring bug.

**Theory — the one new step.** Before the return map, push the backstress forward by
the incremental material rotation $\mathbf R_\Delta=\operatorname{polar}(\mathbf
f_\Delta)$:
$$
\tilde{\boldsymbol\alpha}_{n,k}=\mathbf R_\Delta\,\boldsymbol\alpha_{n,k}\,\mathbf R_\Delta^{\mathsf T},
\qquad \mathbf R_\Delta=\operatorname{polar}(\mathbf F\mathbf F_n^{-1}).
$$
Then the **verified `LadrunoJ2Kernel.h` runs verbatim** (the kernel is frame-agnostic;
this is the payoff of having extracted it). For a rigid increment $\mathbf
f_\Delta=\mathbf Q$ the push gives $\tilde{\boldsymbol\alpha}=\mathbf Q\boldsymbol\alpha\mathbf Q^{\mathsf T}$,
exactly cancelling the rotation of $\mathbf s$; as $\mathbf f_\Delta\to\mathbf I$ it
reduces to the small-strain case. Numpy oracle: superposed-rotation objectivity
$\sim$3e-13 (v1 wrapper fails at 6.17).

**Consistent tangent — two channels.** (A) the standard log-strain dependence
(reused spatial tangent); (B) the **co-rotated backstress** $\tilde{\boldsymbol\alpha}=\mathbf R(\mathbf f_\Delta)\boldsymbol\alpha\mathbf R^{\mathsf T}$, new to this
material. Channel B crosses the tangent-test tolerance ($\sim$1e-4…2e-3) so it is
included; shipped first as a numeric R-perturbation (PR #127), then as the **analytic**
$\partial\mathbf R/\partial\mathbf F$ chain (PR #139). Also ships `-implex` (Δγ-extrapolation,
PR #134): a constant SPD elastic tangent for explicit/quasi-static use and the
enabler when finite J2 is later paired with softening.

**OpenSees.**
```tcl
nDMaterial LadrunoJ2Finite $tag $K $G  -iso voce $sig0 $Qinf $b $Hiso  -kin $N … <-implex>
element LadrunoBrick … $tag -geom finite      ;# the only consumer (3D, F-interface)
```
**When to use it vs the wrapper.** Use `LadrunoJ2Finite` when you need **combined
(kinematic) hardening AND large rotation** — finite cyclic / buckling-brace-like
loops. For **isotropic** hardening at finite strain, the wrapper path
(`LogStrain → LadrunoJ2`, `-kin 0`) is already exact and simpler. 3D-only; default
(no `-implex`) is bit-identical to a fully-implicit run. **ADR: [[16_finite_native_j2_adr]].**

---

## 4 · Family B — finite-strain & staged wrappers

`FiniteStrainNDMaterial`-based (F-interface) and small-strain adaptors that change
*how the element's kinematics reach the real constitutive law*.

### 4.1 `LogStrain` (nDMaterial, 33010) — the Hencky finite-strain lift

**What.** The material-side adaptor that lifts an *unchanged* small-strain 3D
`NDMaterial` to a genuine finite-strain (large rotation + large strain) material by
the logarithmic (Hencky) strain-space technique (dSNPO Box 14.3 `MATISU`). The inner
return map is reused **verbatim**.

**Theory capsule.** Multiplicative split; Eulerian Hencky strain
$\varepsilon^e=\tfrac12\ln\mathbf b^e$; the stress conjugate to it is **Kirchhoff**
$\boldsymbol\tau=\mathbf D^e:\varepsilon^e$ with the *identical* small-strain elastic
matrix; recover Cauchy $\boldsymbol\sigma=\boldsymbol\tau/\det\mathbf F$. The wrapper
does the spectral pre/post-processing (incl. the degeneracy branches at repeated
stretches) and returns the **constitutive** spatial tangent; the element owns the
geometric stiffness.

**Use case.** Any large-strain problem over a GREEN inner law — metal forming /
necking, rubber-block compression, soft-soil large deformation. The detailed treatment
(material matrix, the seam contract, modeling recipes) is in
**[[finite_strain_trifecta_guide]]** and **[[09_finite_strain_material_wrapper]]**.

```tcl
nDMaterial LadrunoJ2 1 $K $G -iso voce 450 265 16.93 129.24   ;# isotropic ⇒ objective
nDMaterial LogStrain 2 1
element LadrunoBrick … 2 -formulation bbar -geom finite
```

### 4.2 `InitDefGrad` / `StagedDefGrad` (nDMaterial, 33013) — finite staged birth

**What.** A `FiniteStrainNDMaterial` wrapper that makes an element **born stress-free
at the current deformed geometry** in a staged analysis (staged construction: a new
member, a concrete lift, a backfill layer). `nDMaterial InitDefGrad $tag $innerTag
<-noInitF> <-F0 f11..f33>`.

**The problem it solves.** OpenSees stores displacement from the original mesh
`X_ref`, which **never re-zeros**. A continuum element appended mid-stage computes
$\varepsilon=\mathbf B\mathbf u$ against `X_ref` and is born carrying the full
accumulated strain. Line elements dodge this with `initialDisp`; continuum elements
had **no such hook**, and the obvious patches fail (additive `InitStrain` is
non-objective at large strain; moving node coords corrupts *shared* nodes).

**Theory — the per-GP birth deformation gradient.** At birth the wrapper captures
$\mathbf F_0=\mathbf I+\operatorname{Grad}_X(\mathbf u_0)$ at each Gauss point, then
feeds the inner the **relative** gradient
$$
\mathbf F_{\text{rel}}=\mathbf F\,\mathbf F_0^{-1},\qquad \mathbf F=\mathbf F_{\text{rel}}\cdot\mathbf F_0,
$$
the standard "stress-free intermediate configuration" multiplicative split (like
$\mathbf F=\mathbf F^e\mathbf F^p$). **Objective** by construction: a post-birth rigid
rotation $\mathbf Q$ gives $\mathbf F\to\mathbf Q\mathbf F$ so
$\mathbf C_{\text{rel}}=\mathbf F_{\text{rel}}^{\mathsf T}\mathbf F_{\text{rel}}$ is
unchanged ⇒ zero spurious stress. Degenerates to identity at $t=0$ ($\mathbf F_0=\mathbf I$).

**The architectural decision (and the subtlety).** It rides on **seam 2** (a reusable
F-interface wrapper) — because it **IS-A** `FiniteStrainNDMaterial`, `LadrunoBrick
-geom finite` consumes it through the *same* `dynamic_cast` with **zero element
edits**. Crucially, in the *spatial* seam the PK1 "push-back" bookkeeping
($\mathbf P=\mathbf P_{\text{rel}}\mathbf F_0^{-\mathsf T}$, weight by $J_0$)
**disappears**: (1) the incremental gradient $\mathbf F_\Delta=\mathbf F\mathbf
F_n^{-1}$ is invariant to $\mathbf F_0$, and (2) Cauchy stress is reference-independent
— so the element integrates $\int\mathbf b^{\mathsf T}\boldsymbol\sigma\,\mathrm dv$ on
the current config with its own $\det\mathbf F$, no $J_0$ weight. *Do not add a $J_0$
weight or $\mathbf F_0^{-\mathsf T}$ push-back to the spatial code — it is already
correct.*

**Trigger = auto-capture** on the first `setTrialF` (which evaluates at the committed
deformed state). Each mid-stage append gets a fresh wrapper instance ⇒ **multi-lift
staged construction falls out for free**, no `InitialStateAnalysis` global toggle.
`-noInitF` opts out (≡ bare inner); `-F0` sets a known birth gradient;
`setResponse("F0")` reports it.

**Use case.** Staged construction of large-deformation continuum models — appending
finite-strain bricks onto an already-deformed mesh and having them born neutral.
**Note: [[staged_deformation_gradiend]].**

### 4.3 `StagedStrain` (nDMaterial, 33014) — small-strain staged birth

**What.** The **small-strain** member of the `Staged*` family (the additive analog of
§4.2). `nDMaterial StagedStrain $tag $innerTag <-noInit> <-eps0 …>`. Captures the
birth strain $\varepsilon_0$ at the first `setTrialStrain` and feeds the inner
$\varepsilon_{\text{rel}}=\varepsilon-\varepsilon_0$, so at birth the element is
**genuinely virgin** (zero stress *and* zero plastic history).

**Why a new class, not upstream `InitStrain`.** `InitStrainNDMaterial` (Petracca) is a
fine *fixed additive prestrain* but bites staged use three ways: (1) **3D-only** —
`exit(-1)`s on a non-3D inner (dead for quad / plane-strain / axisymmetric);
(2) **fixed `ε0`** — no auto-capture, so one tag can't birth a field of elements with
different birth strains; (3) it **adds** a prestrain rather than **subtracting** the
captured birth strain. `StagedStrain` fixes all three: **order-general** (`ε0` sized
to the inner's order; `getCopy(type)` adapts to the element's view), **auto-capturing**,
and **graceful** (no `exit`). The tangent is an **exact passthrough** ($\varepsilon_0$
constant ⇒ $\partial\sigma/\partial\varepsilon=\partial\sigma/\partial(\varepsilon-\varepsilon_0)$),
so no FD-tangent gate is needed. Additive ⇒ small strain only.

**Use case.** Staged construction of ordinary (small-strain) continuum models in 2D
*or* 3D — the everyday staged-build case. `InitStrain` is left untouched for fixed
prestrain (and they **compose**, see §7). **Note: [[staged_deformation_gradiend]].**

---

## 5 · Family C — steel & rebar overlays

Uniaxial materials for reinforced-concrete and steel detailing — a stress-modifying
wrapper and a bond law.

### 5.1 `LadrunoRebarBuckling` (uniaxialMaterial, 33001) — buckling overlay

**What.** A **stress-modifying wrapper** around *any* tension–compression
`UniaxialMaterial`. In compression past a slenderness-dependent onset it applies a
reinforcing-bar **buckling-average** degradation — the section-average compressive
stress of a bar buckling between ties drops well below the bare-bar stress — while
leaving the wrapped material **byte-untouched**:
$$
\sigma_{\text{buckled}}=r(e,\lambda)\cdot\sigma_{\text{bare}},
$$
an opt-in geometric overlay (`-lsr 0` is an identity gate). Two backbone models:
**Dhakal–Maekawa** (`-model dm`) and **Gomes–Appleton** (`-model ga`); v2 adds
cyclic **re-straightening** (a {PASS, BUCKLING, RESTRAIGHTEN} state machine that
raises σ back toward the bare curve on load reversal).

**Architecture.** A wrapper (§2.2) with an analytic tangent and a committed-clone
backbone probe; designed for `LadrunoUniaxialJ2`, composes over `Steel02`/`Steel4`.

**OpenSees.**
```tcl
uniaxialMaterial LadrunoRebarBuckling $tag $steelTag  -model dm -lsr $L_over_D …
```
**Use case.** Longitudinal reinforcing bars in RC columns/walls where post-spalling
bar buckling governs the compression response. **Full reference:
[[LadrunoRebarBuckling_guide]].**

### 5.2 `LadrunoBondSlip` (uniaxialMaterial, 33002) — 1D bond-slip τ–s

**What.** A 1-D **bond-slip** `UniaxialMaterial`: bond stress τ vs slip s, the axial
constitutive slot for embedded-reinforcement modeling (drives `LadrunoEmbeddedRebar`
Mode P). v1 = **MC2010 (Model Code 2010) monotonic backbone** (cyclic → v2).

**Theory.** The classic four-regime CEB-FIP/MC2010 backbone with the ascending
power-law $\tau=\tau_{\max}(s/s_1)^\alpha$, a plateau, a linear descent to the
residual $\tau_f$, exposing $\{\tau_{\max}, s_1, s_2, s_3, \tau_f, \alpha\}$ as inputs
(the regime is **not** hardcoded) plus a **7th explicit input** — the bond fracture
energy $G_F^{\text{bond}}$ — for mesh-objective regularization. Unit contract:
$F=\tau\cdot(\pi d)\cdot L_{\text{trib}}$ (perimeter × tributary length).

**Three must-fix subtleties (all handled in v1):**
- **Initial-slip regularization.** $\tau=\tau_{\max}(s/s_1)^\alpha$ has
  $\mathrm d\tau/\mathrm ds\to\infty$ at $s\to0$ ⇒ singular first tangent. Fixed by a
  **linear segment for $s<s_0\approx0.1 s_1$** (or a tangent cap $K_{\max}$).
- **Softening ⇒ negative tangent.** Past $\tau_{\max}$ the descending branch needs
  **DisplacementControl / ArcLength / IMPLEX** — load control diverges.
- **Mesh objectivity.** Node-lumped bond softening is **not** objective; the post-peak
  slope is regularized by $G_F^{\text{bond}}$ scaled by characteristic length
  ($\text{lch}_{\text{ref}}/\text{lch}$, mirroring ASDConcrete3D). v1 is first-order
  O(h²) (≳6–8 segments per development length).

**Architecture / scope.** v1 = well-confined pull-out; cover-controlled **splitting**
bond failure → v1.1 (a small-cover splice test is xfail'd, to avoid silently
over-predicting bond). Exposes `dissipatedEnergy`/`energy` responses. classTags
33005 (`LadrunoEmbeddedRebar`) / 33002 are collision-free (per-registry namespaces).

**OpenSees.**
```tcl
uniaxialMaterial LadrunoBondSlip $tag $tauMax $s1 $s2 $s3 $tauf $alpha $GFbond
# consumed by the embedded-rebar element's Mode P axial slot:
element LadrunoEmbeddedRebar $et $rebarNode $hostEle $xi $eta $zeta -mode P -bond $tag -perimeter $p -ltrib $L
```
**Use case.** Embedded reinforcement with **explicit bond-slip** (own-DOF rebar + 1D
τ–s, the DIANA-style model) — pull-out, development-length, anchorage studies.
**ADR: [[20_ladruno_embedded_reinforcement_adr]] §D4.**

---

## 6 · The Lemaitre ductile-damage mode

Not a class — an **optional, OFF-by-default mode** on both J2 materials
(`-damage lemaitre $r $s $pD $Dc`), sharing `LadrunoDamage.h` (dSNPO §12.3/§12.4
strain-equivalence / effective-stress ductile damage). Used for ductile fracture
and low-cycle-fatigue studies.

**The key property — exact decoupling.** The return map runs on the **undamaged
effective stress** $\tilde{\boldsymbol\sigma}$ (which carries no damage), so the state
update **decouples exactly**: with damage OFF the behaviour is **byte-identical** to
the undamaged material (no new classTag, no cost). All the $(\Delta\gamma, D)$
coupling lives in the consistent tangent,
$\mathbf D=(1-D)\mathbf D^{\text{alg}}-\tilde{\boldsymbol\sigma}\otimes\partial D/\partial\boldsymbol\varepsilon$.

- **3D (`LadrunoJ2`):** triaxiality-dependent energy release rate $Y$; validated
  A–E + perf, 26/26 (gmsh DEN-bar mesh-objectivity, Coffin-Manson LCF).
- **1D (`LadrunoUniaxialJ2`):** triaxiality ≡ 1/3 ⇒ $R_v\equiv1$ ⇒
  $Y=\tilde\sigma^2/(2E)$; the V4 3D↔1D 1e-6 oracle at triaxiality 1/3.
- **IMPL-EX (`-implex`):** committed extrapolated $\tilde D$ ⇒ SPD tangent for
  softening robustness.

**Full theory + validation: [[15_lemaitre_ductile_damage_adr]]** and the validation
bundle `lemaitre_validation/`.

---

## 7 · Composition — how the wrappers nest

The wrappers are designed to **stack**, and **order matters**. The load-bearing rules:

```
# Finite-strain plasticity (isotropic ⇒ objective)
element → LadrunoBrick(-geom finite) → LogStrain( LadrunoJ2 -kin 0 )

# Staged construction, finite strain
element → LadrunoBrick(-geom finite) → InitDefGrad( LogStrain( LadrunoJ2 ) )

# Staged construction WITH prestress, small strain  (StagedStrain OUTERMOST — hard rule)
element → StagedStrain( InitStrain( realMaterial, ε_pre ) )

# Buckling-prone reinforcing bar
fiber  → LadrunoRebarBuckling( LadrunoUniaxialJ2 )

# Ductile fracture at finite strain (isotropic)
element → LadrunoBrick(-geom finite) → LogStrain( LadrunoJ2 -damage lemaitre … -kin 0 )
```

> [!important] Two ordering rules you must not get wrong
> - **`StagedStrain` outermost over `InitStrain`.** `StagedStrain(InitStrain(mat,
>   ε_pre))`: at birth `StagedStrain` feeds `InitStrain` `ε_rel=0`, `InitStrain` adds
>   `ε_pre`, the real material sees `ε_pre` ⇒ born carrying **exactly** the prestress
>   `σ(ε_pre)` with no inherited geometric stress. Nest the other way and the
>   prestrain is captured into `ε0` and **cancels**. This is *why* staged-birth and
>   prestrain stay separate composable wrappers, not one `-prestrain` mode.
> - **Do not stack `corot` on `finite`, and use isotropic hardening under `finite`.**
>   `LogStrain` handles rotation intrinsically; combined hardening through the wrapper
>   is non-objective at large rotation (use `LadrunoJ2Finite` instead, §3.3).

> [!note] Finite-strain prestrain is NOT a freebie
> Prestrain composes for free in *small* strain (`StagedStrain ∘ InitStrain`, both
> `setTrialStrain`-driven). In *finite* strain there is **no free path** —
> `InitStress`/`InitStrain` are not `FiniteStrainNDMaterial`s and cannot nest in the
> F-chain. A finite prestrain needs a dedicated F-driven wrapper (future work;
> pre-deformation `F_pre` is the clean objective route). See
> [[staged_deformation_gradiend]] §"Deferred".

---

## 8 · Class-tag registry

Per-registry 33000-bands (Element / nDMaterial / uniaxial / Integrator / Recorder
each have their own space, so the same number across registries is **not** a
collision). The fork's material rows:

| classTag | Registry | Material | Define |
|---|---|---|---|
| 33000 | uniaxial | `LadrunoUniaxialJ2` | `MAT_TAG_LadrunoUniaxialJ2` |
| 33001 | uniaxial | `LadrunoRebarBuckling` | `MAT_TAG_LadrunoRebarBuckling` |
| 33002 | uniaxial | `LadrunoBondSlip` | `MAT_TAG_LadrunoBondSlip` |
| 33010 | nDMaterial | `LogStrain` | `ND_TAG_LogStrainNDMaterial` |
| 33011 | nDMaterial | `LadrunoJ2` | `ND_TAG_LadrunoJ2` |
| 33012 | nDMaterial | `LadrunoJ2Finite` | `ND_TAG_LadrunoJ2Finite` |
| 33013 | nDMaterial | `InitDefGrad` (`StagedDefGrad`) | `ND_TAG_InitDefGradNDMaterial` |
| 33014 | nDMaterial | `StagedStrain` | `ND_TAG_StagedStrainNDMaterial` |

(Authoritative source: `SRC/classTags.h` + [[LEDGER_implementations]]. The 1D J2
backstress has no `LadrunoJ2Finite` analogue — a scalar has no frame to co-rotate.)

---

## 9 · Use-case decision guide

| You are modeling… | Reach for | Notes |
|---|---|---|
| Cyclic/seismic steel continuum (panel zone, plate buckling) | `LadrunoJ2` | combined hardening; Bauschinger + ratcheting |
| Steel **fibers** with true ratcheting | `LadrunoUniaxialJ2` | fiber/truss/zeroLength; the ratcheting Menegotto–Pinto can't do |
| Monotonic ductile metal, isotropic | `LadrunoJ2` (`-kin 0`) | ≡ upstream `J2Plasticity`, with hygiene |
| **Large-strain** plasticity, isotropic | `LogStrain(LadrunoJ2 -kin 0)` + `-geom finite` | exact & objective; the trifecta |
| **Large-strain** plasticity, **combined** hardening + large rotation | `LadrunoJ2Finite` + `-geom finite` | the §14.11-native material |
| Large-strain **elastic** (rubber block) | `LogStrain(ElasticIsotropic3D)` | Hencky hyperelasticity |
| Ductile **fracture / LCF** | `LadrunoJ2`/`LadrunoUniaxialJ2` `-damage lemaitre` | OFF by default ⇒ free when unused |
| **Staged construction**, small strain | `StagedStrain(realMat)` | order-general 2D+3D; born virgin |
| **Staged construction**, finite strain | `InitDefGrad(LogStrain(mat))` | objective birth `F₀`; multi-lift for free |
| Reinforcing **bar buckling** | `LadrunoRebarBuckling(LadrunoUniaxialJ2)` | opt-in geometric overlay |
| **Bond-slip / pull-out** of embedded rebar | `LadrunoBondSlip` + `LadrunoEmbeddedRebar` | DisplacementControl/ArcLength past peak |
| Soils, concrete, pressure-sensitive | (upstream / ASDConcrete / ASDPlastic) | the J2 family is **deviatoric-only** |

> [!tip] The deviatoric boundary
> The whole J2 family is **pressure-independent**. For soils/concrete/rock use a
> Drucker–Prager / Mohr–Coulomb / cap / ASDConcrete / ASDPlastic law — and note any
> 3D one of those can itself be lifted to finite strain by `LogStrain` (the
> [[09_finite_strain_material_wrapper]] material matrix lists the GREEN set).

---

## 10 · References

**Per-material deep-dives (in this repo)**
- [[LadrunoJ2_guide]] — `LadrunoJ2` full reference (theory, return map, tangent, views).
- [[LadrunoUniaxialJ2_guide]] — the uniaxial twin.
- [[LadrunoRebarBuckling_guide]] — the buckling overlay (DM/GA, cyclic re-straightening).
- [[finite_strain_trifecta_guide]] — `LogStrain` + the element/wrapper stack.
- [[09_finite_strain_material_wrapper]] — the Hencky adaptor design + GREEN/YELLOW/RED matrix.
- [[16_finite_native_j2_adr]] — `LadrunoJ2Finite` (co-rotating backstress, channels A/B, IMPL-EX).
- [[staged_deformation_gradiend]] — `InitDefGrad`/`StagedStrain` (staged stress-free birth).
- [[15_lemaitre_ductile_damage_adr]] — the Lemaitre damage mode + validation bundle.
- [[20_ladruno_embedded_reinforcement_adr]] — `LadrunoBondSlip` (§D4) + the embedded-rebar element.
- [[LEDGER_implementations]] — authoritative class-tag + PR registry.

**Theory**
- de Souza Neto, Perić & Owen (2008), *Computational Methods for Plasticity* —
  Ch. 6/7 (return mapping, consistent tangent), §9.2/§9.4 (plane stress),
  §12.3/§12.4 (Lemaitre), Ch. 14 (finite-strain multiplicative plasticity, §14.11
  kinematic hardening). The primary framework reference.
- Armstrong & Frederick (1966); Chaboche (1986); Kobayashi & Ohno (2002) — nonlinear
  kinematic hardening and the scalar return-map reduction.
- Dhakal & Maekawa (2002); Gomes & Appleton (1997) — reinforcing-bar buckling.
- *fib* Model Code 2010 / CEB-FIP — the bond-slip τ–s backbone.
- Lemaitre (1985); Lee (1969) — ductile damage; the multiplicative split.

**Source**
- `SRC/material/nD/{LadrunoJ2,LadrunoJ2Finite,LogStrainNDMaterial,FiniteStrainNDMaterial,InitDefGradNDMaterial,StagedStrainNDMaterial}.{h,cpp}`
- `SRC/material/nD/{LadrunoJ2Kernel,LadrunoHardening,LadrunoDamage,LogStrainKernel}.h`
- `SRC/material/uniaxial/{LadrunoUniaxialJ2,LadrunoRebarBuckling,LadrunoBondSlip}.{h,cpp}`

---

> [!info] Maintenance
> This is the **catalog entry point** for fork-authored materials. When a new
> material ships (or a mode is added), add a row to §1 + §8 here, write or link its
> deep-dive, and update [[LEDGER_implementations]] + the banner. Keep the shared
> kernels (`LadrunoJ2Kernel`/`LadrunoHardening`/`LadrunoDamage`) as the single source
> of truth for the J2 family — never fork the return map.
