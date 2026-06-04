---
title: "Staged Activation — Stress-Free Birth of Mid-Stage Elements (StagedStrain / StagedDefGrad)"
project: Ladruno
type: reference / user guide
status: shipped (StagedDefGrad/InitDefGrad finite #163; StagedStrain small-strain 2D+3D #182; review-hardened #183)
classTags:
  - 33013  # ND_TAG_InitDefGradNDMaterial  (StagedDefGrad / InitDefGrad, finite)
  - 33014  # ND_TAG_StagedStrainNDMaterial (StagedStrain, small strain)
materials:
  - nDMaterial StagedStrain
  - nDMaterial StagedDefGrad   # alias of InitDefGrad
related:
  - "[[staged_deformation_gradiend]]"
  - "[[constraints_reference_position]]"
  - "[[solid_transformation_wrapper]]"
  - "[[09_ladruno_brick]]"
  - "[[09_finite_strain_material_wrapper]]"
  - "[[LEDGER_implementations]]"
tags:
  - material
  - staged-construction
  - staged-analysis
  - finite-strain
  - small-strain
  - reference
---

# Staged Activation — Stress-Free Birth of Mid-Stage Elements

> [!abstract] One-line summary
> The **`Staged*` family** is two `nDMaterial` wrappers that make a continuum element
> **born stress-free at the current deformed geometry** when it is appended *after*
> the model has already deformed (staged construction: a new member, a concrete lift,
> a backfill layer). **`StagedStrain`** does it for **small-strain** elements
> (additive re-reference `ε_rel = ε − ε₀`); **`StagedDefGrad`** (= **`InitDefGrad`**)
> does it for **finite-strain** elements (multiplicative, objective re-reference
> `F_rel = F·F₀⁻¹`). Both **auto-capture** the birth state on the element's first
> evaluation — one material tag births a whole lift, no global toggle.

This document is the **descriptive reference**: the problem, the continuum theory, the
architecture exactly as coded, the OpenSees wiring, and the intended use. For the
chronological design log see [[staged_deformation_gradiend]]; for the reference-frame
behaviour of SP/MP constraints and line vs continuum elements see
[[constraints_reference_position]].

---

## Contents

1. [Intended use — and what it is *not*](#1-intended-use--and-what-it-is-not)
2. [The problem — the never-rezeroed displacement frame](#2-the-problem--the-never-rezeroed-displacement-frame)
3. [The family at a glance](#3-the-family-at-a-glance)
4. [Theory](#4-theory)
5. [Architecture](#5-architecture)
6. [Command syntax](#6-command-syntax)
7. [How to use — recipes](#7-how-to-use--recipes)
8. [Prestrain by composition (and the finite gap)](#8-prestrain-by-composition-and-the-finite-gap)
9. [Recording, responses & parameters](#9-recording-responses--parameters)
10. [Gotchas, limitations & boundaries](#10-gotchas-limitations--boundaries)
11. [Verification & validation](#11-verification--validation)
12. [References & related](#12-references--related)

---

## 1. Intended use — and what it is *not*

**Use it when** you build part of a model, load it, commit, and then **add elements**
that should not inherit the accumulated deformation as stress:

- staged construction / sequential erection (frames, walls, bridges);
- soil/rock placement (backfill lifts, embankment construction);
- concrete pours / lifts cast against an already-deformed structure;
- any "activate this element now, stress-free, at wherever the nodes currently are".

**It is not:**

- **Element *deactivation*.** OpenSees has no general stress-relieving element removal;
  `remove element` simply deletes. The `Staged*` family is **activation only**.
- **A prestrain/prestress facility.** It births an element with **zero** stress. For an
  imposed prestrain you *compose* with the upstream `InitStrain` — see §8.
- **A new element or analysis mode.** It is a constitutive wrapper; the analysis driver,
  integrator, and element are unchanged. Staging itself is your normal Python/Tcl loop
  (`analyze → addElement → analyze`).

---

## 2. The problem — the never-rezeroed displacement frame

OpenSees stores nodal displacement `u` measured from the **original mesh `X_ref` at
t = 0**, and that frame **never re-zeros** — not at a stage boundary, not after
`wipeAnalysis` ([[constraints_reference_position]] §1):

```
position = X_ref + u
```

So an element appended after the model has drifted to `u = u₀` computes its strain from
`X_ref` and is **born carrying the full accumulated deformation as spurious stress**.

Line elements (`Truss`, beam-columns) dodge this: they snapshot the committed nodal
displacement at birth and subtract it (`Truss`'s `initialDisp`, folded into the length
`L`). **Continuum elements have no such hook** — `FourNodeQuad`, `stdBrick`, tets compute
`ε = B·u` straight from the node pointers. The `Staged*` family is the missing hook,
implemented at the **material** level so it needs **zero element changes**.

> [!info] Why the material level (seam-2)?
> The re-referencing rides inside the material the element already holds. Each element
> gets its own per-Gauss-point material copies (`getCopy`), so a single tag births any
> number of elements, each snapshotting its **own** birth state. No `CrdTransf`-style
> element seam is needed, and the element is oblivious. This is the small/finite analogue
> of `Truss`'s `initialDisp`, lifted to a full strain/deformation-gradient capture.

---

## 3. The family at a glance

| | **`StagedStrain`** | **`StagedDefGrad`** ( = `InitDefGrad`) |
|---|---|---|
| Regime | **small strain** | **finite strain** |
| classTag | `33014` | `33013` |
| Base class | `NDMaterial` | `FiniteStrainNDMaterial` |
| Driven by | `setTrialStrain(ε)` | `setTrialF(F)` |
| Re-reference | additive `ε_rel = ε − ε₀` | multiplicative `F_rel = F·F₀⁻¹` |
| Captures | birth strain `ε₀` | birth deformation gradient `F₀` |
| Dimensions | **2D + 3D** (PlaneStrain/Stress/AxiSymm/3D) | 3D |
| Inner material | any small-strain `NDMaterial` | a `FiniteStrainNDMaterial` (e.g. `LogStrain`) |
| Element | `FourNodeQuad`, `stdBrick`, `LadrunoBrick -geom linear/corot`, tris/tets | `LadrunoBrick -geom finite` |
| Objective under post-birth rotation | n/a (rotations infinitesimal) | **yes** |
| Tangent | exact passthrough (no FD gate) | spatial, FD-gated |
| Opt-out flag | `-noInit` | `-noInitF` |
| Explicit-birth flag | `-eps0 e1 e2 …` | `-F0 f11 … f33` |

**Pick by the element's strain regime.** A small-strain element → `StagedStrain`; a
`-geom finite` element → `StagedDefGrad`. The type system enforces the match
(`StagedDefGrad` rejects a non-finite inner; `StagedStrain` wraps a strain-driven one).

---

## 4. Theory

### 4.1 Small strain — additive re-reference (`StagedStrain`)

At activation, capture the birth strain `ε₀ = ε(birth) = B·u₀`. Thereafter feed the
inner material the **relative** strain:

```
ε_rel = ε − ε₀          (deformation measured from the stress-free birth state)
```

At birth `ε_rel = 0`, so the inner (started pristine) returns **zero stress and zero
plastic history** — the element is born **genuinely virgin** (§4.5). Because `ε₀` is a
constant, the tangent is an **exact passthrough**: `∂σ/∂ε = ∂σ/∂(ε − ε₀)` — no
finite-difference tangent gate is needed. This is **exact for small strain** (the strain
measure is linear in `u`, rotations are infinitesimal).

### 4.2 Finite strain — multiplicative re-reference (`StagedDefGrad`)

At activation, capture the birth **deformation gradient** `F₀ = I + Grad_X(u₀)`.
Thereafter feed the inner finite-strain material the **relative** gradient:

```
F_rel = F · F₀⁻¹        (the multiplicative split F = F_rel · F₀, Lee)
```

This is the standard *stress-free intermediate configuration* idea — structurally
identical to `F = Fₑ·Fₚ` in finite plasticity or `F = Fₑ·F_g` in growth mechanics. At
birth `F_rel = I` → zero stress, virgin history. It **degenerates correctly**: appended
at t = 0 (`u₀ = 0`) ⇒ `F₀ = I` ⇒ `F_rel = F`, ordinary behaviour with zero overhead.

### 4.3 Why additive fails for finite strain (the rotation proof)

The reason finite strain needs the **multiplicative** form (and `InitStrain`'s additive
offset is *wrong* there) is **objectivity under large rotation**. Subtracting `u₀`
builds the gradient additively:

```
F_naive = I + Grad_X(u − u₀) = F − F₀ + I
```

Append the element at the deformed shape, then apply a **pure rigid-body rotation `Q`**
(no new strain — it must stay stress-free). Under that rotation `F → Q·F₀`, so:

```
F_naive = Q·F₀ − F₀ + I  ≠  I    for any Q ≠ I   →  spurious strain → spurious stress  ✗
```

The multiplicative version does not manufacture stress:

```
F_rel = F·F₀⁻¹ → (Q·F₀)·F₀⁻¹ = Q   ⇒   C_rel = F_relᵀF_rel = QᵀQ = I   →  zero stress  ✓
```

Deformation gradients **compose multiplicatively** (`F = F_rel·F₀`), they do not subtract.
Additive subtraction coincides with the truth only in the infinitesimal limit — which is
exactly why `StagedStrain` (additive) is correct for small strain and `StagedDefGrad`
(multiplicative) is required for finite strain.

### 4.4 The spatial seam (why the PK1 bookkeeping vanishes)

The shipped finite seam (`FiniteStrainNDMaterial`: `setTrialF(F)` in, **Cauchy σ** out;
the element integrates `∫bᵀσ dv` on the *current* configuration with its own `det F`) is
**spatial**, not PK1 / total-Lagrangian. Two facts collapse the `F₀` bookkeeping that a
PK1 derivation would carry (`P = P_rel·F₀⁻ᵀ`, weight by `J₀`):

1. **The incremental gradient is invariant to `F₀`:**
   `F_Δ = F_rel·F_rel,n⁻¹ = (F·F₀⁻¹)(F_n·F₀⁻¹)⁻¹ = F·F_n⁻¹`. So `F₀` changes **only** the
   stress-free initial state, never the step-to-step response.
2. **Cauchy stress is reference-independent.** The inner returns `σ = τ_rel / J_rel`,
   which **is** the true physical Cauchy stress (`J_rel = det F_rel = det F / det F₀`).
   The element integrates on the current config with `det F`, so **no `F₀⁻ᵀ` push-back
   and no `J₀` integration weight** are needed.

So the wrapper is a near pass-through: capture `F₀`, forward `F_rel`, delegate stress /
tangent / state. The one genuine consistency point — the `J_rel` inside the spatial
tangent `c = (1/2J_rel)[D:L:B]` matching the element's `dv` — is enforced by the
element's **FD-tangent gate** (the same arbiter that validated the base finite seam).

### 4.5 Born genuinely virgin

In both regimes the appended element is born **virgin**, not "stress zeroed but history
carried". The inner material starts pristine (`bᵉ = I` / `εᵖ = 0`), and the relative
measure at birth is the identity (`F_rel = I` / `ε_rel = 0`), so the constitutive law
returns zero stress **and** zero internal variables. The element then accrues stress and
plastic history only from deformation **beyond** birth — exactly the physics of a member
cast/placed into a deformed structure.

---

## 5. Architecture

### 5.1 Material-level (seam-2) wrappers

Both are thin `NDMaterial` wrappers around an inner material. Neither touches any element
code: `LadrunoBrick -geom finite` already obtains its material via
`dynamic_cast<FiniteStrainNDMaterial*>` and calls `setTrialF`; `FourNodeQuad`/`stdBrick`
call `setTrialStrain`. The wrappers slot in transparently because they *are* the
respective base class.

### 5.2 `StagedStrain` internals

- **Dimension-general.** `ε₀`/buffers are sized to `inner->getOrder()` (3 for
  PlaneStrain/PlaneStress, 4 for AxiSymmetric, 6 for 3D). `getType()`/`getOrder()`
  delegate to the inner, so the element sees its own view.
- **The `getCopy` adoption pattern (a real OpenSees subtlety).** The generic material
  bases (`ElasticIsotropicMaterial`, …) implement the no-arg `getCopy(void)` as
  `exit(-1)` "subclass responsibility". So the constructor builds a canonical
  `"ThreeDimensional"` template, and `getCopy(type)` re-derives the element's actual
  view (PlaneStrain/Stress/AxiSymm/3D) via the inner's **inherited `getCopy(type)`
  switch** (which rebuilds from `E, ν, ρ`). The copy then **adopts** that typed inner
  directly — never re-deriving back to 3D (which would feed a 2D element an order-6
  material). One adopting constructor + one `adoptCopy()` helper serve the factory and
  both `getCopy` paths; construction is graceful (`isValid()`), with **no `exit(-1)`**.
- **State.** `{active, captured, eps0Explicit, ε₀}` + `totalStrain` (for the recorder).
  Capture latches on the first `setTrialStrain`, survives `commit`/`revertToLastCommit`,
  re-arms only on `revertToStart` (auto mode). All serialized in `sendSelf`/`recvSelf`.

### 5.3 `StagedDefGrad` / `InitDefGrad` internals

- **Wraps an inner `FiniteStrainNDMaterial`** — typically `LogStrain` over a GREEN 3D
  small-strain `NDMaterial` (Elastic, `LadrunoJ2`, …). **IS-A** `FiniteStrainNDMaterial`,
  so the finite element consumes it with zero edits.
- On `setTrialF(F)`: capture `F₀` (first call), forward `setTrialF(F·F₀⁻¹)`; delegate
  `getStress` (Cauchy σ), `getTangent`, `getSpatialTangentTensor`, `getJ`, and the full
  state cycle to the inner. `F₀` + `captured` are serialized.

### 5.4 The auto-capture lifecycle (per-GP, committed birth)

```
ops.element(...)                ──► getCopy(type) per Gauss point
                                     each copy: captured=false  (the COPY, once, at creation)
   (prior stage already committed → nodes deformed)
analyze() step = birth:
   first update → setTrialStrain/F at the COMMITTED deformed state
                                     ──► capture ε₀ / F₀  (the REBASE, once, latched)
   thereafter: forward ε−ε₀ / F·F₀⁻¹ on every step, every iteration
   commitState / revertToLastCommit  ──► ε₀/F₀ untouched (persist)
```

Two distinct moments: the **copy** is at `ops.element(...)`; the **rebase** is on the
first material evaluation, which for a freshly-appended element lands at the **committed
deformed state** (the integrator predictor does not move displacement). Each appended
element gets fresh copies ⇒ **multi-lift falls out for free**, no `InitialStateAnalysis`
global toggle. `F₀`/`ε₀` then ride every subsequent step for the element's whole life.

### 5.5 Files, classTags, dispatch

| | `StagedStrain` | `StagedDefGrad` / `InitDefGrad` |
|---|---|---|
| Source | `SRC/material/nD/StagedStrainNDMaterial.{cpp,h}` | `SRC/material/nD/InitDefGradNDMaterial.{cpp,h}` |
| classTag | `ND_TAG_StagedStrainNDMaterial = 33014` | `ND_TAG_InitDefGradNDMaterial = 33013` |
| Dispatch names | `StagedStrain`, `StagedStrainNDMaterial` | `InitDefGrad`, `StagedDefGrad` (+ `…NDMaterial`) |
| Tests | `tests/test_stagedStrain_material.py` | `tests/test_initDefGrad_material.py` |
| Inner seam | any small-strain `NDMaterial` | `FiniteStrainNDMaterial` ([[09_finite_strain_material_wrapper]]) |

---

## 6. Command syntax

### 6.1 `StagedStrain` (small strain)

```
nDMaterial StagedStrain $tag $innerTag <-noInit> <-eps0 e1 e2 ...>
```

- `$innerTag` — any small-strain 3D-capable `NDMaterial` (it must expose a
  `"ThreeDimensional"` view; the 2D views are re-derived from it).
- `-noInit` — pass-through (`ε₀ ≡ 0`); bit-identical to the bare inner.
- `-eps0 e1 e2 …` — explicit birth strain (engineering-shear Voigt; component count =
  the element order). Skips auto-capture. The **subtract** convention (opposite sign to
  `InitStrain`'s additive prestrain). If the size doesn't match the element's order it is
  ignored with a warning and the element auto-captures instead.

### 6.2 `StagedDefGrad` / `InitDefGrad` (finite strain)

```
nDMaterial StagedDefGrad $tag $innerTag <-noInitF> <-F0 f11 f12 f13 f21 f22 f23 f31 f32 f33>
nDMaterial InitDefGrad   $tag $innerTag <-noInitF> <-F0 ...>     ;# identical (alias)
```

- `$innerTag` — a `FiniteStrainNDMaterial` (e.g. a `LogStrain` tag).
- `-noInitF` — pass-through (`F₀ ≡ I`).
- `-F0 …` — explicit birth deformation gradient (9 row-major components).

> [!warning] Flag spelling differs by member
> Small strain uses `-noInit`; finite uses `-noInitF`. (Historical: `InitDefGrad` shipped
> first with `-noInitF`/`-F0`.) The opt-out semantics are identical.

---

## 7. How to use — recipes

### 7.1 The rule

> **If your element runs `-geom finite` (drives `setTrialF`), wrap with `StagedDefGrad`
> over a `LogStrain`-lifted material. Otherwise (small-strain / linear / corot element)
> wrap with `StagedStrain` over the plain material.** The constitutive model rides
> *inside* either path unchanged.

### 7.2 Finite-strain staged block

```python
# ---- Stage 1: existing structure deforms and commits ----
# ... build nodes, stage-1 elements, apply load ...
ops.analyze(1)                                  # nodes move to the deformed shape, committed

# ---- Stage 2: append a finite-strain block, born stress-free ----
ops.nDMaterial("ElasticIsotropic", 11, E, nu)   # (or LadrunoJ2, etc.)
ops.nDMaterial("LogStrain",        12, 11)       # lift to finite strain (seam-3)
ops.nDMaterial("StagedDefGrad",    13, 12)       # <-- staged-activation wrapper

ops.element("LadrunoBrick", 99, *conn, 13, "-formulation", "std", "-geom", "finite")

ops.integrator("LoadControl", 0.0)              # optional zero-load equilibration
ops.analyze(1)                                  # element 99 auto-captures F0 at birth -> ~0 stress
# subsequent staged loading: element 99 responds only to NEW deformation
```

### 7.3 Small-strain staged block

```python
# ---- Stage 1: holder deforms + commits (nodes are now displaced) ----
ops.analyze(1)

# ---- Stage 2: append a small-strain block, born stress-free ----
ops.nDMaterial("ElasticIsotropic", 21, E, nu)   # or J2Plasticity, ASDConcrete3D, ...
ops.nDMaterial("StagedStrain",     23, 21)       # <-- one tag births the whole lift

# 3D brick, or a 2D quad — same wrapper tag, dimension-general:
ops.element("stdBrick", 99, *conn3d, 23)
# ops.element("quad", 99, *conn2d, 1.0, "PlaneStrain", 23)

ops.integrator("LoadControl", 0.0)
ops.analyze(1)                                  # auto-captures eps0 at birth -> ~0 stress
```

### 7.4 Multi-lift

Because capture is **per-Gauss-point copy**, a single wrapper tag births any number of
elements, each snapshotting its own birth state:

```python
ops.nDMaterial("ElasticIsotropic", 11, E, nu)
ops.nDMaterial("LogStrain", 12, 11)
ops.nDMaterial("StagedDefGrad", 13, 12)          # ONE wrapper, defined once

for lift in lifts:
    # ... prior structure already deformed & committed ...
    for e, conn in lift.elements:
        ops.element("LadrunoBrick", e, *conn, 13, "-geom", "finite")  # same tag 13
    ops.analyze(1)        # every new element captures its OWN F0 here, independently
```

Use a *different* tag per lift only if the lifts have different **constitutive**
properties (e.g. concrete strength) — not for the `F₀`/`ε₀` capture, which is per-copy.

### 7.5 The lifecycle on the pseudo-time axis

- `ε₀`/`F₀` are captured **once**, at the element's first evaluation after it is added
  (= the committed deformed birth state), then **frozen** and reused on every step for
  the element's whole life.
- The latch survives `commitState` and a non-converged `revertToLastCommit`; a genuine
  `revertToStart` (full restart) re-arms auto-capture.
- The inner material chains its own state (`bᵉ`, plastic strain) **relative to birth**,
  exactly as a t = 0 element would relative to `X_ref`.

---

## 8. Prestrain by composition (and the finite gap)

The `Staged*` family births **stress-free**. For an imposed **prestrain/prestress** in
**small strain**, *compose* with the upstream `InitStrain`, nesting `StagedStrain`
**outermost** (element-facing):

```python
ops.nDMaterial("ElasticIsotropic", 1, E, nu)
ops.nDMaterial("InitStrain",       2, 1, *eps_pre)   # additive prestrain (inner)
ops.nDMaterial("StagedStrain",     3, 2)             # staged birth (OUTER)
# element -> StagedStrain( InitStrain( real, eps_pre ) )  =>  born carrying sigma(eps_pre)
```

At birth `StagedStrain` feeds `InitStrain` `ε_rel = 0`; `InitStrain` adds `ε_pre`; the
real material sees `ε_pre` → born carrying **exactly the prestress `σ(ε_pre)`**, with no
inherited geometric stress.

> [!important] Order is load-bearing
> `StagedStrain` must be the **outer** wrapper. Nest the other way
> (`InitStrain(StagedStrain(...))`) and the prestrain is captured into `ε₀` and
> **cancels** — born fully stress-free, prestrain ignored.

**Finite-strain prestrain is out of scope (deferred).** The finite chain is
`setTrialF`-driven and `InitStress`/`InitStrain` are `setTrialStrain`-driven, so they
**cannot nest** in the F-chain — there is no free composition. A finite prestrain needs a
dedicated F-driven wrapper; the recommended route is a **pre-deformation `F_pre`** (feed
the inner `F_rel·F_pre` → `σ(F_pre)`, objective), with a co-rotating initial Cauchy
stress `σ_pre` as the harder variant. See [[staged_deformation_gradiend]] §"Deferred".

---

## 9. Recording, responses & parameters

**`StagedStrain`** (`eleResponse … material <i> <chan>` or recorder `material`):

| channel | returns |
|---|---|
| `eps0` / `birthStrain` | the captured birth strain `ε₀` |
| `totalStrain` | the total strain from `X_ref` |
| `stress` / `strain` / `tangent` | delegated to the inner — `strain` is the **relative** `ε − ε₀` (the element's physical state) |
| any inner-specific channel | delegated (e.g. `plasticStrain`, `backStress`) |

**`StagedDefGrad` / `InitDefGrad`:**

| channel | returns |
|---|---|
| `F0` / `birthF` / `F0capture` | the captured birth deformation gradient `F₀` (9, row-major) |
| `stress` / `tangent` / … | delegated to the inner (`LogStrain` reports Cauchy σ, relative Hencky strain, spatial tangent) |

Both delegate `setParameter` to the inner, and `StagedStrain` keeps `getStressSensitivity`
/`commitSensitivity` (DDM) parity with `InitStrain`.

---

## 10. Gotchas, limitations & boundaries

> [!warning] Pick the wrapper to match the element regime
> A small-strain element wrapped in `StagedDefGrad` (or vice-versa) is a type error the
> wrappers reject. `StagedDefGrad` needs a `FiniteStrainNDMaterial` inner (`LogStrain`);
> `StagedStrain` needs a strain-driven inner.

- **Why not the upstream `InitStrain` for staged birth?** `InitStrain` is **3D-only**
  (and `exit(-1)`-kills on a non-3D inner), has a **fixed** `ε₀` (no per-element
  auto-capture), and **adds** a prestrain rather than subtracting the birth strain.
  `StagedStrain` fixes all three for staged use; `InitStrain` remains the right tool for
  a *known imposed prestrain* (compose them, §8).
- **Capture timing.** `ε₀`/`F₀` is taken on the **first** material evaluation. For a
  freshly-appended element that is the committed deformed birth state (the integrator
  predictor does not move displacement). A pathological garbage first-trial that later
  reverts would latch a wrong birth reference — the same first-touch assumption the whole
  family (and `Truss`'s `initialDisp`) makes.
- **`getCopy(type)` invariant.** Dimension generality assumes
  `inner.getCopy("3D").getCopy(type) == inner.getCopy(type)` — true for materials that
  route `getCopy(type)` through `E, ν, ρ` (ElasticIsotropic, J2Plasticity, …). A material
  that bakes different state into its 3D-vs-2D copy paths would violate it. The factory
  validates only the 3D view, so a 2D-incapable inner fails (graceful `null`) at
  `getCopy(type)`, not at registration.
- **`setTrialStrainIncr`.** Supported, but reconstructing the absolute birth strain from
  increments alone is fundamentally limited; the common continuum elements drive
  `setTrialStrain(total)`. (`-noInit` is a faithful incremental pass-through as of #183.)
- **No deactivation.** Activation only (§1).

---

## 11. Verification & validation

**`StagedStrain`** — `tests/test_stagedStrain_material.py`, driven through `FourNodeQuad`
(2D) + `stdBrick` (3D):

- `-noInit` ≡ bare inner (stress + force), 2D and 3D;
- explicit `-eps0` reduce-to-relative;
- **two-phase staged stress-free birth in 2D *and* 3D** (also the dimensional-generality
  proof — order-3 capture for the quad, order-6 for the brick);
- PlaneStress staged birth;
- **J2 inner born virgin** (no inherited plastic strain after a yielded holder);
- **composition `StagedStrain∘InitStrain`** → born stress independent of the birth
  deformation and equal to `σ(ε_pre)`;
- mismatched-`eps0` graceful fallback; missing-inner graceful failure.

**`StagedDefGrad` / `InitDefGrad`** — `tests/test_initDefGrad_material.py`, through
`LadrunoBrick -geom finite`:

- `-noInitF` ≡ bare inner; explicit-`F0` stress-free birth;
- **post-birth rigid-rotation objectivity** (the additive-offset failure mode);
- reduce-to-relative (`-F0=G` at `F=M` Cauchy σ == bare inner at `F=M·G⁻¹`);
- **wrapped consistent tangent == FD of the resisting force** (the wrapper preserves
  Newton);
- **genuine two-phase staged construction** (deform+commit holder → append auto-capture
  brick → born ~zero force, holder unchanged).

Both batteries are Zone-A (CI-gated). The implementation was hardened by a high-effort
code review (#183) — all `exit(-1)` paths removed, `getCopy` unified, `setTrialStrainIncr`
passthrough fixed.

---

## 12. References & related

- [[staged_deformation_gradiend]] — the design note / ADR: the `F₀` theory, the
  spatial-vs-PK1 reconciliation, the small-strain companion, and the deferred
  finite-prestrain plan.
- [[constraints_reference_position]] — the SP / MP / line-vs-continuum reference-frame
  trichotomy that motivates this (why `fix` snaps to reference, `equalDOF` doesn't, and
  continuum elements inherit strain).
- [[solid_transformation_wrapper]] — the `linear`/`corot`/`finite` geometry-method layer
  the finite wrapper rides on.
- [[09_ladruno_brick]] — `LadrunoBrick`, the first solid consumer (`-geom finite`).
- [[09_finite_strain_material_wrapper]] — `LogStrain` / `FiniteStrainNDMaterial`, the
  seam-3 inner the finite wrapper consumes.
- [[LEDGER_implementations]] — the build-control rows (classTags 33013 / 33014, PRs).
- de Souza Neto, Perić & Owen, *Computational Methods for Plasticity* (2008) — Ch. 14
  (the spatial multiplicative / log-strain framework the finite inner implements).
