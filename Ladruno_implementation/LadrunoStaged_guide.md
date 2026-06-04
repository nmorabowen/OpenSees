---
title: "Staged* — Stress-Free Staged-Construction Birth (InitDefGrad / StagedStrain)"
project: Ladruno
type: reference / user guide
status: shipped (v1 — InitDefGrad/StagedDefGrad finite 33013, StagedStrain small-strain 33014)
classTag:
  - 33013
  - 33014
material:
  - nDMaterial InitDefGrad
  - nDMaterial StagedStrain
related:
  - "[[staged_deformation_gradiend]]"
  - "[[constraints_reference_position]]"
  - "[[LogStrain_guide]]"
  - "[[finite_strain_trifecta_guide]]"
  - "[[LadrunoJ2_guide]]"
  - "[[LadrunoBrick_reference]]"
  - "[[solid_transformation_wrapper]]"
  - "[[Ladruno_materials_guide]]"
tags:
  - material
  - staged-analysis
  - finite-strain
  - wrapper
  - reference
  - guide
---

# Staged* — Stress-Free Staged-Construction Birth

> [!abstract] One-line summary
> The `Staged*` family makes a **continuum element born genuinely neutral at the
> current deformed geometry** when it is appended mid-analysis (staged construction:
> a new member, a concrete lift, a backfill layer). Two members, one idea:
> **`StagedStrain`** (`nDMaterial`, classTag 33014) is the **small-strain** member —
> it subtracts the captured birth strain $\boldsymbol\varepsilon_0$; **`InitDefGrad`**
> (alias `StagedDefGrad`, `nDMaterial`, classTag 33013) is the **finite-strain**
> member — it captures the birth deformation gradient $\mathbf F_0$ and feeds the
> inner the *relative* gradient $\mathbf F_{\text{rel}}=\mathbf F\,\mathbf F_0^{-1}$.
> Both auto-capture at birth, so **multi-lift staged construction falls out for free**
> with no global `InitialStateAnalysis` toggle.

This document is the **descriptive reference** for the shipped family: the shared
problem, the two mechanics (additive small-strain vs. multiplicative finite-strain),
the spatial-seam reconciliation that makes the finite wrapper a near pass-through, the
OpenSees wiring, the composition rules, and the boundaries. For the chronological
design log and the theory derivation see [[staged_deformation_gradiend]]; for the
staged-analysis frame context (why `X_ref` never re-zeros) see
[[constraints_reference_position]].

---

## Contents

- [[#1 Intended use cases|1. Intended use cases]]
- [[#2 The shared problem — `X_ref` never re-zeros|2. The shared problem]]
- [[#3 The two mechanics at a glance|3. The two mechanics]]
- [[#4 InitDefGrad / StagedDefGrad — the finite-strain member|4. InitDefGrad (finite)]]
- [[#5 StagedStrain — the small-strain member|5. StagedStrain (small)]]
- [[#6 Composition — order is a hard rule|6. Composition]]
- [[#7 OpenSees implementation|7. OpenSees implementation]]
- [[#8 Command syntax and recipes|8. Command syntax]]
- [[#9 Verification and validation|9. Verification]]
- [[#10 Limitations and deferred work|10. Limitations]]
- [[#11 References|11. References]]

---

## 1. Intended use cases

The `Staged*` family exists for one job: **staged construction / sequential analysis
of continuum models**, where part of the model is built, loaded, and held while new
elements are added onto the already-deformed mesh. The appended element must be **born
neutral** — zero stress at birth — and start picking up load only from the next
increment.

| You are modeling… | Reach for | Why |
|---|---|---|
| **Everyday small-strain staged build** (2D or 3D) — adding a slab, a wall lift, a backfill layer to a model that has settled under its own weight | `StagedStrain( realMat )` | order-general (quad / plane-strain / plane-stress / axisymmetric / brick); born virgin (zero stress *and* zero plastic history) |
| **Large-deformation staged construction** — appending finite-strain bricks onto a mesh that has already deformed substantially | `InitDefGrad( LogStrain( realMat ) )` + `-geom finite` | objective birth gradient $\mathbf F_0$; survives large rotation/stretch at birth |
| **Staged build *with* a prestress** (small strain) | `StagedStrain( InitStrain( realMat, ε_pre ) )` | composes — born carrying exactly $\sigma(\varepsilon_{\text{pre}})$, no inherited geometric stress (§6) |
| **Multi-lift sequence** (several appends at different stages) | either member | each append is a fresh wrapper instance ⇒ each captures its *own* birth state automatically |

> [!tip] Pick the member by **strain regime**, not by element
> `StagedStrain` is **additive** ($\boldsymbol\varepsilon_{\text{rel}}=\boldsymbol\varepsilon-\boldsymbol\varepsilon_0$) and therefore correct for **small strain
> only**. `InitDefGrad` is **multiplicative** ($\mathbf F_{\text{rel}}=\mathbf F\mathbf
> F_0^{-1}$) and **objective**, so it is the one to use whenever the birth (or
> post-birth) deformation is large enough that a rigid rotation would otherwise inject
> spurious stress. When in doubt at finite strain, use `InitDefGrad`.

> [!note] What this is **not**
> This family does not *deactivate* or stress-relieve an existing element — OpenSees
> has no stress-relieving element deactivation (see [[constraints_reference_position]]
> §7). It only controls the **birth** state of a newly-added element. It is also not a
> prestress/prestrain facility on its own — prestress comes from *composing* with
> `InitStrain` (§6).

---

## 2. The shared problem — `X_ref` never re-zeros

OpenSees stores nodal displacement $\mathbf u$ measured from the **original mesh
$X_{\text{ref}}$ at $t=0$**, and that frame **never re-zeros** — not at a stage
boundary, not after a `wipeAnalysis`, not ever (see
[[constraints_reference_position]] §1). Position is always

$$\text{position} = X_{\text{ref}} + \mathbf u .$$

So when an analysis deforms under stage 1 and a new continuum element is then appended
to part of the model, that element computes its strain against $X_{\text{ref}}$ —

$$\boldsymbol\varepsilon = \mathbf B\,\mathbf u \quad(\text{small strain}),
\qquad
\mathbf F = \mathbf I + \operatorname{Grad}_X(\mathbf u) \quad(\text{finite}),$$

— and is therefore **born carrying the full accumulated strain**. That is exactly the
opposite of what staged construction wants.

**Line elements dodge this; continuum elements do not.** A `Truss` / beam-column folds
an `initialDisp` / `nodeI*InitialDisp` snapshot into its length / corotational
transform, so it is born stress-free at the deformed shape by default. **Continuum
elements have no such hook.** And the two obvious patches both fail:

| Naive patch | Why it fails |
|---|---|
| **Additive material offset** (`InitStrainNDMaterial`, `InitialStateAnalysisWrapper`) | Exact for *small* strain, but at large deformation it is **non-objective** — a post-birth rigid rotation injects spurious strain — and **additive** where finite strain demands multiplicative composition. |
| **Move the node coordinates** $X_{\text{ref}}\to X_{\text{ref}}+\mathbf u_0$ | **Wrong** — the nodes are *shared* with the already-active elements, which still need the original $X_{\text{ref}}$. Corrupts the neighbors. |

The fix must therefore be **local to the appended element**, must keep the single
$X_{\text{ref}}$ frame, and — at finite strain — must be **multiplicative** and
**objective**. The `Staged*` family delivers exactly that, as a **material wrapper**
so that no element needs editing.

---

## 3. The two mechanics at a glance

Both members are wrappers that hold an inner material, snapshot its birth kinematics on
the first evaluation, and thereafter feed the inner the deformation **relative to
birth**. They differ only in the kinematic measure they subtract/divide out:

| | **`StagedStrain`** (33014) | **`InitDefGrad` / `StagedDefGrad`** (33013) |
|---|---|---|
| Strain regime | **small** strain | **finite** strain (large rotation + stretch) |
| Base class | `NDMaterial` (strain-driven) | `FiniteStrainNDMaterial` (F-driven) |
| Birth snapshot | $\boldsymbol\varepsilon_0$ (the birth strain) | $\mathbf F_0$ (the birth deformation gradient) |
| Fed to inner | $\boldsymbol\varepsilon_{\text{rel}}=\boldsymbol\varepsilon-\boldsymbol\varepsilon_0$ (**additive**) | $\mathbf F_{\text{rel}}=\mathbf F\,\mathbf F_0^{-1}$ (**multiplicative**) |
| Objective under post-birth rotation? | n/a (small strain) | **yes** — by construction |
| Born | **virgin** (zero stress *and* zero plastic history) | stress-free at birth config |
| Tangent | exact passthrough — **no FD-tangent gate** | spatial seam — validated by the element's **FD-tangent gate** |
| Opt-out flag | `-noInit` | `-noInitF` |
| Explicit birth value | `-eps0 …` | `-F0 f11..f33` |
| Consumed by | strain-based elements: `FourNodeQuad`, `stdBrick`, … (2D + 3D) | `LadrunoBrick -geom finite` (the F-interface) |

Both **degenerate to identity at $t=0$** (birth at the undeformed state ⇒
$\boldsymbol\varepsilon_0=\mathbf 0$ / $\mathbf F_0=\mathbf I$ ⇒ zero overhead), exactly
like a `Truss`'s nonzero-disp guard on `initialDisp`.

---

## 4. InitDefGrad / StagedDefGrad — the finite-strain member

### 4.1 The clean idea: a per-Gauss-point birth deformation gradient

At birth the wrapper captures, **at each integration point**, the deformation gradient
of the already-accumulated (committed) displacement field:

$$\mathbf F_0 = \mathbf I + \operatorname{Grad}_X(\mathbf u_0)
\qquad(\mathbf u_0 = \text{committed disp at birth, gradient w.r.t. } X_{\text{ref}}).$$

Thereafter the element computes the **total** $\mathbf F$ from $X_{\text{ref}}$ exactly
as always (nothing about the framework changes), and the wrapper feeds the inner the
**relative** deformation gradient

$$\boxed{\;\mathbf F_{\text{rel}} = \mathbf F\,\mathbf F_0^{-1},
\qquad \mathbf F = \mathbf F_{\text{rel}}\cdot\mathbf F_0.\;}$$

This is the **multiplicative split** with $\mathbf F_0$ mapping
$X_{\text{ref}}\to\text{birth}$ (the stress-free intermediate configuration) and
$\mathbf F_{\text{rel}}$ mapping $\text{birth}\to\text{current}$ (the stress-driving
part) — structurally identical to $\mathbf F=\mathbf F^e\mathbf F^p$ in finite
plasticity or $\mathbf F=\mathbf F^e\mathbf F_g$ in growth mechanics. It is standard
"stress-free intermediate configuration" theory, not a hack.

> [!check] Objective by construction — the property the additive offset lacks
> Under a post-birth rigid rotation $\mathbf Q$: $\mathbf F\to\mathbf Q\mathbf F$, so
> $\mathbf F_{\text{rel}}\to\mathbf Q\,\mathbf F_{\text{rel}}$, and the relative
> right Cauchy–Green tensor
> $$\mathbf C_{\text{rel}} = \mathbf F_{\text{rel}}^{\mathsf T}\mathbf Q^{\mathsf T}\mathbf Q\,\mathbf F_{\text{rel}} = \mathbf F_{\text{rel}}^{\mathsf T}\mathbf F_{\text{rel}}$$
> is **unchanged** ⇒ **zero spurious stress** under rigid rotation. This is exactly the
> failure mode of additive `InitStrain` at large strain, cured.

### 4.2 The spatial reconciliation — the PK1 push-back disappears

The design note ([[staged_deformation_gradiend]]) develops the energy-consistency
accounting in **first-PK / total-Lagrangian** terms: push the stress back through
$\mathbf F_0$ ($\mathbf P=\mathbf P_{\text{rel}}\,\mathbf F_0^{-\mathsf T}$) and weight
the integration by $J_0=\det\mathbf F_0$. **That bookkeeping is correct for the PK1
dialect — but the seam that actually shipped is spatial, and in the spatial dialect it
collapses to nothing.** Two facts do it:

1. **The incremental gradient is invariant to $\mathbf F_0$.**
   $$\mathbf F_\Delta = \mathbf F_{\text{rel}}\,\mathbf F_{\text{rel},n}^{-1}
   = \bigl(\mathbf F\mathbf F_0^{-1}\bigr)\bigl(\mathbf F_n\mathbf F_0^{-1}\bigr)^{-1}
   = \mathbf F\,\mathbf F_n^{-1}.$$
   So $\mathbf F_0$ changes **only** the stress-free *initial* elastic state at birth;
   the step-to-step response (and hence the inner's $\mathbf b^e$ evolution) is
   identical with or without it.
2. **Cauchy stress is reference-independent.** The inner returns
   $\boldsymbol\sigma=\boldsymbol\tau_{\text{rel}}/J_{\text{rel}}$, which **is** the
   true physical Cauchy stress (the Kirchhoff $\boldsymbol\tau_{\text{rel}}$ is relative
   to the birth config, $J_{\text{rel}}=\det\mathbf F_{\text{rel}}=\det\mathbf F/\det\mathbf F_0$). The element integrates $\int\mathbf b^{\mathsf T}\boldsymbol\sigma\,\mathrm dv$ on the **current** configuration with its **own** $\det\mathbf F$.

Net result: there is **no $\mathbf P=\mathbf P_{\text{rel}}\mathbf F_0^{-\mathsf T}$
push-back** (a PK1 artifact) and **no $J_0$ integration weight** (the current-config
integration already carries the right measure). The inner's `getJ()` returns
$J_{\text{rel}}$, which the element does **not** use for assembly. The wrapper is a
near pass-through.

> [!warning] Do not "fix" the spatial code with a PK1 correction
> It is tempting to read the design note's §"The subtlety" and add a $J_0$ integration
> weight or an $\mathbf F_0^{-\mathsf T}$ push-back to the spatial code. **Do not.** In
> the spatial seam the code is **already correct**; adding either term would
> double-count and break it. The PK1 accounting is kept in the design note for
> intuition only.

The one genuine consistency point — the $J_{\text{rel}}$ inside the spatial tangent
$\mathbf c=\tfrac{1}{2J_{\text{rel}}}[\mathbb D:\mathbb L:\mathbb B]$ matching the
element's $\mathrm dv=\det\mathbf F\,\mathrm dV$ — is enforced by the element's
**finite-difference tangent gate**, the same arbiter that validated the base finite
seam and F-bar.

### 4.3 Trigger: auto-capture on the first `setTrialF`

The wrapper captures $\mathbf F_0$ on its **first `setTrialF(F)` call**, and that first
call is precisely the birth gradient: `Domain::addElement → setDomain → update`
evaluates the appended element at the **committed deformed state**, so the wrapper's
first $\mathbf F$ is $\mathbf F_{\text{birth}}$. After capture, every call forwards
$\mathbf F_{\text{rel}}=\mathbf F\,\mathbf F_0^{-1}$ to the inner.

> [!important] Multi-lift falls out for free
> Each mid-stage append gets a **fresh wrapper instance**, so each captures its **own**
> birth gradient at its own activation time. Multi-lift staged construction needs **no
> `InitialStateAnalysis`-style global toggle** — it is automatic.

Controls:
- **`-noInitF`** — opt out; the wrapper becomes a transparent pass-through ($\equiv$ the
  bare inner).
- **`-F0 f11..f33`** — supply a known birth gradient explicitly (the nine row-major
  components). It is **kept across `revertToStart`**; auto-capture re-arms instead.
- **`setResponse("F0")`** — reports the captured (or supplied) birth gradient.

### 4.4 Why a material wrapper (seam-2), not an element hook (seam-1)

Two seams were possible: (1) the element captures $\mathbf F_0$ per GP and handles the
push-back/weight itself, or (2) a reusable F-interface **material wrapper**. The fork
shipped **seam-2** because the element already exposes the F-interface
(`FiniteStrainNDMaterial::setTrialF`): the wrapper needs **zero element edits** and is
**reusable across any future F-based solid**. Because `InitDefGrad` **IS-A**
`FiniteStrainNDMaterial`, `LadrunoBrick -geom finite` consumes it through the *same*
`dynamic_cast<FiniteStrainNDMaterial*>` it already does. Seam-1 would have edited the
brick's finite path for no gain.

`InitDefGrad` rides **on top of** the existing `finite` (total-Lagrangian) kinematics
method of the [[solid_transformation_wrapper]] — it is the capture-at-birth hook on
that `finite` ledger, not a fourth kinematics method.

---

## 5. StagedStrain — the small-strain member

### 5.1 Mechanic: subtract the captured birth strain

`StagedStrain` captures the birth strain $\boldsymbol\varepsilon_0$ at the **first
`setTrialStrain`** (= birth), then feeds the inner the relative strain

$$\boxed{\;\boldsymbol\varepsilon_{\text{rel}} = \boldsymbol\varepsilon - \boldsymbol\varepsilon_0.\;}$$

At birth $\boldsymbol\varepsilon_{\text{rel}}=\mathbf 0$, so the element is born
**genuinely virgin** — zero stress *and* zero plastic history (a J2 inner inherits no
plastic strain). The tangent is an **exact passthrough**, because $\boldsymbol\varepsilon_0$ is constant:

$$\frac{\partial\boldsymbol\sigma}{\partial\boldsymbol\varepsilon}
= \frac{\partial\boldsymbol\sigma}{\partial(\boldsymbol\varepsilon-\boldsymbol\varepsilon_0)} .$$

So — unlike the finite wrapper — **no FD-tangent gate is needed**. Being additive,
`StagedStrain` is valid for **small strain only**; finite strain uses `StagedDefGrad`.

### 5.2 Why a new class, not the upstream `InitStrain`

`InitStrainNDMaterial` (Petracca 2024) is a fine **fixed, additive prestrain** material,
but it bites *staged* use three ways. `StagedStrain` fixes all three:

| `InitStrainNDMaterial` limitation | Consequence | `StagedStrain` fix |
|---|---|---|
| **3D-only** — `getOrder()` hardcoded to 6, and it `exit(-1)`s (kernel kill) on a non-3D inner | Dead for `FourNodeQuad` / plane-strain / plane-stress / axisymmetric | **Order-general** — $\boldsymbol\varepsilon_0$ sized to `inner->getOrder()`, `getCopy(type)` adapts the inner to the element's view; 2D *and* 3D |
| **Fixed $\varepsilon_0$** — no auto-capture | One tag can't birth a *field* of elements that each have a different birth strain | **Auto-capturing** — each instance snapshots its own birth strain at first `setTrialStrain` |
| **Adds** a prestrain | Wrong sign for staged birth (you want to *remove* the inherited strain) | **Subtracts** the captured birth strain |
| `exit(-1)` on a bad inner | Kills the Python kernel | **Graceful** — no `exit` |

`InitStrain` is left **untouched** for its intended job (fixed prestrain) — and the two
**compose** (§6).

---

## 6. Composition — order is a hard rule

The `Staged*` wrappers are designed to **nest**, and **order matters**. This is why
staged-birth and prestrain stay *separate composable wrappers* instead of one
`-prestrain`-mode class — staged-birth, prestrain, and any future offset all come from
nesting.

### 6.1 Small-strain staged build **with** a prestress

Nest with **`StagedStrain` OUTERMOST** (element-facing):

```text
element → StagedStrain( InitStrain( realMaterial, ε_pre ) )
```

At birth `StagedStrain` feeds `InitStrain` $\boldsymbol\varepsilon_{\text{rel}}=\mathbf 0$, `InitStrain` adds $\boldsymbol\varepsilon_{\text{pre}}$, and the real
material sees $\boldsymbol\varepsilon_{\text{pre}}$ ⇒ born carrying **exactly the
prestress** $\boldsymbol\sigma(\boldsymbol\varepsilon_{\text{pre}})$, with **no
inherited geometric stress**.

> [!warning] Nest the other way and the prestrain cancels
> If you put `InitStrain` outermost — `InitStrain( StagedStrain( realMaterial ),
> ε_pre )` — the prestrain is captured into $\boldsymbol\varepsilon_0$ at birth and
> **cancels**. The order `StagedStrain( InitStrain( … ) )` is a **hard rule**.

### 6.2 Finite-strain staged construction

```text
element → LadrunoBrick(-geom finite) → InitDefGrad( LogStrain( realMaterial ) )
```

`InitDefGrad` wraps the `LogStrain` Hencky lift (which itself wraps any GREEN
small-strain 3D law), giving an objective stress-free birth at large deformation.

### 6.3 Finite-strain prestrain is **not** a freebie

Prestrain composes for free in **small** strain (`StagedStrain ∘ InitStrain`, both
`setTrialStrain`-driven). In **finite** strain there is **no free path**, because the
finite chain is `setTrialF`-driven and `InitStress` / `InitStrain` are **not**
`FiniteStrainNDMaterial`s — they cannot nest in the F-chain. A finite prestrain
therefore needs a **dedicated F-driven wrapper** (future work), with a genuine design
fork:

1. **Pre-deformation $\mathbf F_{\text{pre}}$** (recommended) — feed the inner
   $\mathbf F_{\text{rel}}\cdot\mathbf F_{\text{pre}}$, so at stress-free birth it sees
   $\mathbf F_{\text{pre}}\to\boldsymbol\sigma(\mathbf F_{\text{pre}})$. The clean
   multiplicative analog of `InitStrain`'s additive $\boldsymbol\varepsilon_{\text{pre}}$; **objective by construction**; would compose as
   `StagedDefGrad( FinitePrestrain( … ) )`.
2. **Initial Cauchy stress $\boldsymbol\sigma_{\text{pre}}$** (the `InitStress` analog)
   — additive on stress, but a *fixed* $\boldsymbol\sigma_{\text{pre}}$ does **not**
   co-rotate ⇒ **non-objective under large rotation** (the exact failure mode that
   motivated `InitDefGrad`). A correct version must co-rotate the prestress. Harder;
   only if a consumer needs it.

> [!note] The unvetted partial path
> A partial path exists — `LogStrain( InitStress( real ) )`, the small-strain
> `InitStress` at the bottom of the Hencky stack — but its **large-rotation
> objectivity is unverified**, which is precisely why finite prestrain is its own task,
> not a freebie.

---

## 7. OpenSees implementation

### 7.1 Classes and tags

| Item | `InitDefGrad` (finite) | `StagedStrain` (small) |
|---|---|---|
| Command | `nDMaterial InitDefGrad` (alias `StagedDefGrad`) | `nDMaterial StagedStrain` |
| Class | `InitDefGradNDMaterial` | `StagedStrainNDMaterial` |
| Base | `FiniteStrainNDMaterial` (F-interface) | `NDMaterial` (strain-interface) |
| classTag | **`ND_TAG_InitDefGradNDMaterial = 33013`** | **`ND_TAG_StagedStrainNDMaterial = 33014`** |
| Birth state | $\mathbf F_0$ per GP | $\boldsymbol\varepsilon_0$ per GP (order-general) |
| Source | `SRC/material/nD/InitDefGradNDMaterial.{h,cpp}` | `SRC/material/nD/StagedStrainNDMaterial.{h,cpp}` |

Both follow the fork's **wrapper / composition pattern** (deep-copy the inner,
serialize the inner classTag + state, delegate `commitState` / `revertToLastCommit` /
`revertToStart` / `sendSelf` / `recvSelf` to the inner), and both obey the fork hygiene
rules — **no `exit()`**, real `revert*`, real serialization (see
[[Ladruno_materials_guide]] §2.2, §2.4).

> The class-tag bands are per-registry; 33013 / 33014 sit in the `nDMaterial` band.
> Authoritative source: `SRC/classTags.h` + [[LEDGER_implementations]].

### 7.2 The birth-capture cycle

| Phase | `InitDefGrad` | `StagedStrain` |
|---|---|---|
| First evaluation (birth) | first `setTrialF(F)` → capture $\mathbf F_0$, forward $\mathbf F_{\text{rel}}=\mathbf F\mathbf F_0^{-1}$ | first `setTrialStrain(ε)` → capture $\boldsymbol\varepsilon_0$, forward $\boldsymbol\varepsilon_{\text{rel}}=\varepsilon-\varepsilon_0$ |
| Subsequent | forward $\mathbf F\mathbf F_0^{-1}$ to the inner | forward $\varepsilon-\varepsilon_0$ to the inner |
| `commitState` | delegate to inner ($\mathbf F_0$ frozen) | delegate to inner ($\varepsilon_0$ frozen) |
| `revertToStart` | re-arm auto-capture (keep an explicit `-F0`/`-eps0` if supplied) | re-arm auto-capture |
| Stress / tangent out | inner's Cauchy $\boldsymbol\sigma$ + spatial tangent (no push-back) | inner's $\boldsymbol\sigma$ + tangent (exact passthrough) |

---

## 8. Command syntax and recipes

### 8.1 Grammar

```tcl
# Finite-strain staged birth (alias StagedDefGrad)
nDMaterial InitDefGrad  $tag  $innerTag  <-noInitF>  <-F0 $f11 $f12 $f13 $f21 $f22 $f23 $f31 $f32 $f33>

# Small-strain staged birth
nDMaterial StagedStrain $tag  $innerTag  <-noInit>   <-eps0 $e0 $e1 ...>
```

| Option | Member | Meaning |
|---|---|---|
| `-noInitF` | `InitDefGrad` | opt out — transparent pass-through ($\equiv$ bare inner) |
| `-F0 f11..f33` | `InitDefGrad` | supply a known birth gradient (9 row-major components); kept across `revertToStart` |
| `-noInit` | `StagedStrain` | opt out — transparent pass-through ($\equiv$ bare inner) |
| `-eps0 …` | `StagedStrain` | supply an explicit birth strain (sized to the inner's order) |

With no option the wrapper **auto-captures at birth** (the staged-construction default).

### 8.2 Small-strain staged build (the everyday case)

```tcl
# Phase 1: a base material settles under self-weight (the "holder" region).
nDMaterial ElasticIsotropic 1 200000.0 0.3
element FourNodeQuad ...

# Phase 2: append a NEW quad region with a StagedStrain-wrapped material.
#          It auto-captures its birth strain → born neutral at the deformed shape.
nDMaterial StagedStrain 10 1            ;# wrap material 1
# element FourNodeQuad <new> ... 10     ;# appended mid-analysis
```

```python
import openseespy.opensees as ops
ops.nDMaterial("ElasticIsotropic", 1, 200000.0, 0.3)
# ... phase-1 analysis runs and commits, the mesh deforms ...
ops.nDMaterial("StagedStrain", 10, 1)          # auto-capture birth strain
ops.element("FourNodeQuad", 99, *new_nodes, 1.0, "PlaneStrain", 10)
# the appended quad is born stress-free at the current deformed geometry
```

### 8.3 Small-strain staged build **with** a prestress (composition)

```tcl
# StagedStrain OUTERMOST over InitStrain (hard rule — see §6.1)
nDMaterial ElasticIsotropic 1 200000.0 0.3
nDMaterial InitStrain   20 1 $eps_pre     ;# fixed prestrain ε_pre
nDMaterial StagedStrain 30 20             ;# staged birth, neutral except for σ(ε_pre)
# element ... 30   →  born carrying exactly σ(ε_pre), no inherited geometric stress
```

### 8.4 Finite-strain staged construction

```tcl
# element → LadrunoBrick(-geom finite) → InitDefGrad( LogStrain( realMaterial ) )
nDMaterial LadrunoJ2   1 $K $G -iso voce 350.0 60.0 10.0 0.0   ;# isotropic ⇒ objective inner
nDMaterial LogStrain   2 1                                     ;# Hencky finite lift
nDMaterial InitDefGrad 3 2                                     ;# auto-capture F0 at birth
element LadrunoBrick $eleTag $n1 $n2 $n3 $n4 $n5 $n6 $n7 $n8  3  -geom finite
```

```python
ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", 350.0, 60.0, 10.0, 0.0)
ops.nDMaterial("LogStrain", 2, 1)
ops.nDMaterial("InitDefGrad", 3, 2)            # or "StagedDefGrad"
ops.element("LadrunoBrick", 1, *node_tags, 3, "-geom", "finite")
# inspect the captured birth gradient:
F0 = ops.eleResponse(1, "material", 1, "F0")
```

> [!tip] Prefer an **isotropic** inner under finite staged construction
> Wrapping `LogStrain( LadrunoJ2 -kin 0 )` keeps the chain objective. Combined
> (kinematic) hardening through the `LogStrain` wrapper is **not** exactly objective
> under large rotation (the dSNPO §14.11 boundary — see [[LadrunoJ2_guide]] §10.1);
> for finite kinematic hardening use `LadrunoJ2Finite` instead.

---

## 9. Verification and validation

Both materials have **no Python `setTrialF`/`setTrialStrain` entry point**, so they are
validated **through real elements** in Zone-A pytest batteries.

### 9.1 `InitDefGrad` — `tests/test_initDefGrad_material.py`

All cases run through `LadrunoBrick -geom finite`:

| Test intent | What it pins |
|---|---|
| `-noInitF` pass-through | wrapper $\equiv$ bare inner (stress *and* resisting force) |
| explicit `-F0=G`, impose $\mathbf F=G$ | **stress-free birth** |
| **post-birth rigid rotation** $\mathbf F=\mathbf Q\cdot G$ | **zero stress** (the additive-`InitStrain` failure mode, cured) |
| reduce-to-relative | `InitDefGrad(-F0=G)` at $\mathbf F=M$ Cauchy $\boldsymbol\sigma$ == bare inner at $\mathbf F=M\,G^{-1}$ |
| wrapped consistent tangent | == FD of the resisting force (the wrapper **preserves Newton**) |
| **genuine two-phase staged construction** | deform + commit a holder brick, append an auto-capture `InitDefGrad` brick on the same nodes → born with **~zero force**, holder **unchanged**, reported $\mathbf F_0$ == the deformed gradient |

### 9.2 `StagedStrain` — `tests/test_stagedStrain_material.py`

Run through `FourNodeQuad` (2D) and `stdBrick` (3D):

| Test intent | What it pins |
|---|---|
| `-noInit` | $\equiv$ bare (2D *and* 3D) |
| explicit-`eps0` | reduce-to-relative |
| **two-phase staged stress-free birth, 2D *and* 3D** | born neutral; also the **dimensional-generality proof** — order-3 capture for the quad, order-6 for the brick |
| `PlaneStress` birth | the plane-stress view births correctly |
| **J2 inner born virgin** | no inherited plastic strain |
| **composition `StagedStrain∘InitStrain`** | born stress independent of the birth deformation and equal to $\boldsymbol\sigma(\varepsilon_{\text{pre}})$ |
| **graceful failure** | succeeds on a non-3D inner exactly where `InitStrain` `exit(-1)`s |

---

## 10. Limitations and deferred work

> [!caution] Known boundaries
> - **`StagedStrain` is small-strain only** — it is additive. For large deformation use
>   `InitDefGrad`, which is multiplicative and objective.
> - **No element deactivation / stress relief.** The family controls only the *birth*
>   state of a *newly-added* element; OpenSees has no stress-relieving deactivation of
>   an existing element (see [[constraints_reference_position]] §7).
> - **Finite-strain combined hardening** inside the inner is not exactly objective
>   under large rotation when wrapping `LogStrain( LadrunoJ2 )` — that is the inner's
>   §14.11 boundary, *not* a `Staged*` defect; the staged wrapper itself is objective.
>   Use an isotropic inner, or `LadrunoJ2Finite`.
> - **Finite-strain prestrain is not yet shipped** — there is no F-driven prestrain
>   wrapper; the small-strain `StagedStrain ∘ InitStrain` route has no finite analog
>   (§6.3). Deferred (v2) for `InitDefGrad`: a `totalStrain` recorder channel (Hencky
>   from $X_{\text{ref}}$, on top of the default *relative* Hencky channel).

> [!note] Do not re-add the PK1 push-back
> The single most likely correctness mistake when touching `InitDefGrad` is to add a
> $J_0$ integration weight or an $\mathbf F_0^{-\mathsf T}$ push-back to the spatial
> code. **The spatial seam is already correct** (§4.2); both terms would double-count.

---

## 11. References

**Within this repo**
- [[staged_deformation_gradiend]] — the **design log**: the $\mathbf F_0$ theory, the
  PK1↔spatial reconciliation, the small-strain companion, the composition rules, and
  the deferred finite-prestrain analysis. (This guide is the how-to; that is the why.)
- [[constraints_reference_position]] — the single-frame ($X_{\text{ref}}$ never
  re-zeros) and SP/MP/element trichotomy that motivates the whole family.
- [[LogStrain_guide]] / [[finite_strain_trifecta_guide]] — the Hencky finite-strain
  lift that `InitDefGrad` wraps.
- [[LadrunoJ2_guide]] — the flagship inner law; §10.1 is the §14.11 objectivity boundary
  referenced above.
- [[LadrunoBrick_reference]] / [[solid_transformation_wrapper]] — the `-geom finite`
  element and the `linear`/`corot`/`finite` kinematics layer `InitDefGrad` rides on.
- [[Ladruno_materials_guide]] — the material catalog (§4.2/§4.3 staged wrappers, §7
  composition rules).
- [[LEDGER_implementations]] — authoritative class-tag + PR registry.

**Theory**
- de Souza Neto, Perić & Owen (2008), *Computational Methods for Plasticity* — Ch. 14
  (finite-strain multiplicative split, the stress-free intermediate configuration);
  §14.11 (kinematic hardening at finite strain). The primary framework reference.
- Lee (1969) — the multiplicative decomposition $\mathbf F=\mathbf F^e\mathbf F^p$ that
  $\mathbf F=\mathbf F_{\text{rel}}\mathbf F_0$ mirrors.

**Source**
- `SRC/material/nD/InitDefGradNDMaterial.{h,cpp}` (classTag 33013),
  `SRC/material/nD/StagedStrainNDMaterial.{h,cpp}` (classTag 33014).
- Tests: `tests/test_initDefGrad_material.py`, `tests/test_stagedStrain_material.py`.
