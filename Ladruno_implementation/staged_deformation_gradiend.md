---
title: Staged deformation gradient — stress-free activation of large-deformation continuum elements
project: Ladruno
status: IMPLEMENTED v1 (InitDefGradNDMaterial, seam-2, classTag 33013)
tags:
  - design
  - staged-analysis
  - finite-strain
  - elements
  - opensees-internals
related:
  - "[[constraints_reference_position]]"
  - "[[solid_transformation_wrapper]]"
  - "[[09_ladruno_brick]]"
  - "[[LEDGER_quirks]]"
---

# Staged deformation gradient (`F₀`)

> **Status:** IMPLEMENTED as `InitDefGradNDMaterial` (`nDMaterial InitDefGrad`,
> classTag 33013) — the **seam-2** material-wrapper route below. This note keeps the
> theory; the **Implementation (v1)** section at the bottom records what shipped and,
> importantly, **reconciles the PK1 dialect of the "energy consistency" section with
> the spatial seam that actually landed** (the `J₀`/push-back accounting is subsumed —
> read that before re-deriving anything from §"The subtlety"). See
> [[constraints_reference_position]] §4 for the motivating context.

## Problem

In a staged analysis, the model deforms under stage 1, then we add new elements to
part of the model (staged construction: a new member, a lift of concrete, a backfill
layer). We want the appended element **born stress-free at the current deformed
geometry** — it should carry zero stress at birth and start picking up load only from
the next increment.

OpenSees stores displacement `u` from the original mesh `X_ref`, and that frame
**never re-zeros** (see [[constraints_reference_position]] §1). So an appended element
computes `ε = B·u` against `X_ref` and is born carrying the full accumulated strain.
The line elements (truss, beam-column) dodge this with an `initialDisp` /
`nodeI*InitialDisp` snapshot, but **continuum elements have no such hook**, and the
obvious patches don't generalize:

- **Additive material offset** (`InitStrainNDMaterial`, `InitialStateAnalysisWrapper`):
  exact for *small* strain, but for a large-deformation element it is **non-objective**
  (a post-birth rigid rotation injects spurious strain) and **additive** where finite
  strain demands multiplicative composition. → trouble.
- **Move the node coordinates** `X_ref → X_ref + u₀`: **wrong** — nodes are *shared*
  with the already-active elements, which still need the original `X_ref`.

The fix must therefore be **local to the appended element**, **multiplicative**, and
**objective** — and it must stay inside the single `X_ref` frame.

## The clean idea: per-Gauss-point birth deformation gradient

At activation, the element captures, at **each integration point**, the deformation
gradient of the already-accumulated (committed) displacement field:

```
F₀ = I + Grad_X(u₀)          # gradient w.r.t. X_ref, u₀ = committed disp at birth
```

Thereafter the element computes the **total** `F` from `X_ref` exactly as always
(nothing about the framework changes) and feeds the material the **relative**
deformation gradient:

```
F      = I + Grad_X(u)        # still measured from X_ref — the never-rezeroed frame
F_rel  = F · F₀⁻¹             # deformation from the birth (stress-free) config
```

The material's strain measure is built from `F_rel` (`C_rel = F_relᵀ F_rel`,
`E_rel = ½(C_rel − I)`, log strain, etc.). This is the **multiplicative split**

```
F = F_rel · F₀
```

with `F₀` mapping `X_ref → birth` (stress-free) and `F_rel` mapping `birth → current`
(stress-driving) — structurally identical to `F = Fₑ·Fₚ` in finite plasticity or
`F = Fₑ·F_g` in growth mechanics. It is standard "stress-free intermediate
configuration" theory, not a hack.

## Why it is the *right* object within `X_ref`

- **Single frame preserved.** `F` is still computed from the displacement field
  measured from `X_ref`. We only **post-multiply by a stored local tensor**. `X_ref`
  remains the one global reference and never re-zeros; `F₀` is a per-GP "where was I
  born" map.
- **Objective** — the property the additive offset lacked. Under a post-birth rigid
  rotation `Q`: `F → QF`, so `F_rel → Q·F_rel`, and
  `C_rel = F_relᵀ Qᵀ Q F_rel = C_rel` is **unchanged** → zero spurious stress under
  rigid rotation.
- **Multiplicative**, as finite strain requires (`E_rel ≠ E_total − E₀` except in the
  infinitesimal limit).
- **History-variable hygiene.** `F₀` is frozen at birth and lives in the GP state
  alongside the material's own history — commits/reverts with everything else, no
  global bookkeeping, no analysis-rebuild quirk.
- **Degenerates correctly.** Built at `t = 0` (`u₀ = 0`) ⇒ `F₀ = I` ⇒ `F_rel = F`,
  i.e. ordinary behavior with zero overhead — exactly like `Truss`'s nonzero-disp
  guard on `initialDisp`.

## The subtlety that must be respected (energy consistency)

> [!warning] PK1 dialect — superseded by the spatial seam in practice.
> The accounting below is written in **first-PK / total-Lagrangian** terms
> (`P = P_rel·F₀⁻ᵀ`, weight by `J₀`). The seam that actually shipped
> (`FiniteStrainNDMaterial`: `setTrialF(F)` in, **Cauchy σ** out, element integrates
> `∫bᵀσ dv` on the *current* config with its own `det F`) is **spatial**, and in that
> dialect this bookkeeping **disappears** — see **Implementation (v1) → the spatial
> reconciliation**. Keep this section for the PK1 intuition, but do **not** add a
> `J₀` weight or an `F₀⁻ᵀ` push-back to the spatial code: it is already correct.

The element still integrates over the **original reference volume** `V₀`, but the
material's stress-free state is now the *birth* config. To keep the internal
force/energy and the tangent consistent, push the stress back through `F₀` and scale
the integration weight by `J₀ = det F₀`:

- material returns the stress work-conjugate to `F_rel` (e.g. first PK `P_rel = ∂ψ/∂F_rel`);
- assemble with `P = P_rel · F₀⁻ᵀ` (chain rule through `F_rel = F·F₀⁻¹`);
- weight by `J₀` so the integral is effectively over the birth volume, not `V₀`.

Get this mapping right and **both** the material and geometric tangents stay
consistent → quadratic convergence preserved. This is the classic prestressed /
pre-deformed reference formulation; the only "new" storage is `F₀` (a 3×3, or 2×2 in
2-D) per Gauss point.

## Comparison with the alternatives

| Approach | Stays in `X_ref`? | Large-def correct? | Cost |
|---|---|---|---|
| Additive material offset (`InitStrain*`, ISA wrapper) | yes | **no** (non-objective, additive) | trivial |
| Move node coords `X_ref → X_ref + u₀` | **no** — corrupts *shared* nodes | n/a | breaks neighbors |
| Proxy-node element wrapper | yes | only small-strain | heavy plumbing |
| **Per-GP `F₀`, `F_rel = F·F₀⁻¹`** | **yes** | **yes** | per-GP `F₀` + stress map |

`F₀` is the finite-strain generalization of what already exists: `Truss`'s
`initialDisp` folded into `L` is a 1-D `F₀`; the beam `CrdTransf`'s `nodeIInitialDisp`
is the rotation-free version. We are lifting that to the full deformation gradient so
it survives large rotation and stretch.

## Relationship to the solid transformation wrapper

`F₀` is **not** a fourth kinematics method — it rides *on top of* the existing
`finite` (total-Lagrangian) method of the [[solid_transformation_wrapper]]. That
layer answers *"how do I compute strain/stress under large rotation and strain"*
(its `linear`/`corot`/`finite` axis); `F₀` answers the orthogonal staged-construction
question *"what is this element's stress-free reference when it is born mid-stage."*
Concretely: the wrapper's `finite` ledger already produces `F = I + Grad_X(u)`; the
staged-activation feature is the **capture of `F₀` at birth + the post-multiply
`F_rel = F·F₀⁻¹`** before that `F` reaches the material/seam-3 adaptor. So this is a
capture-at-activation hook on the `finite` ledger of [[solid_transformation_wrapper]],
and its first consumer would be the same element, [[09_ladruno_brick]] (`-geom finite`).

## Where it should live

Two candidate seams (decide at implementation time):

1. **Element/GP kinematics layer** — the appended finite-strain continuum element
   stores `F₀` per GP, computes `F_rel`, and handles the `P = P_rel·F₀⁻ᵀ`, `J₀`
   mapping. Most direct; element-specific.
2. **An `F`-interface nD material wrapper** — a *multiplicative* analog of
   `InitStrainNDMaterial`: stores `F₀`, receives `F` from the element, computes
   `F_rel = F·F₀⁻¹`, and returns the mapped stress/tangent. Reusable across any
   `F`-based finite-strain element. (Distinct from `InitStrainNDMaterial`, which is
   additive on `ε` and therefore only small-strain-correct.)

Seam 2 is the more reusable; it depends on whether the target elements expose an
`F`-based material interface or a strain-based one.

## Open questions for the implementation pass

- **Capture trigger.** Snapshot `F₀` in `setDomain`/first `update` when committed disp
  is nonzero (mirroring `Truss`), or via an explicit "activate" command analogous to
  `InitialStateAnalysis on/off`? An explicit toggle is clearer for true staged
  construction with multiple lifts.
- **Which elements first.** Pick the finite-strain continuum element(s) actually in
  use in the fork's models as the pilot; small-strain elements can keep using the
  additive material offset (it's exact for them).
- **Strain measure & material API.** Confirm the target material's conjugacy
  (`P`/`F`, `S`/`E`, or rate form) so the `F₀` push-back factor is derived once,
  correctly.
- **Recorder semantics.** Decide whether recorders report total strain (from `X_ref`)
  or relative strain (from birth). Likely expose both; default to relative for the
  appended element so post-processing matches its physical state.
- **Class tags / ledger.** New element or material ⇒ allocate a class tag and add rows
  to [[LEDGER_implementations]] when code lands.

## Implementation (v1) — `InitDefGradNDMaterial`

Shipped 2026-06-03 as the **seam-2** route (the reusable F-interface material
wrapper), because the finite-strain seam it needed already exists in-tree.

**What it is.** `nDMaterial InitDefGrad $tag $innerTag <-noInitF> <-F0 f11..f33>`
(classTag `ND_TAG_InitDefGradNDMaterial = 33013`). A `FiniteStrainNDMaterial`
subclass wrapping an inner `FiniteStrainNDMaterial` (e.g. `LogStrain`). Because it
**IS-A** `FiniteStrainNDMaterial`, `LadrunoBrick -geom finite` consumes it through
the same `dynamic_cast<FiniteStrainNDMaterial*>` it already does — **zero element
edits**. On `setTrialF(F)`: capture `F₀` (first call), then forward
`setTrialF(F·F₀⁻¹)` to the inner; `getStress`/`getTangent`/`getSpatialTangentTensor`/
`getJ`/`getStrain`/commit/revert/sendSelf delegate down.

**The spatial reconciliation (read this before re-deriving §"The subtlety").** The
note above is PK1; the shipped seam is spatial, and two facts collapse the `F₀`
bookkeeping to nothing:

1. **The incremental gradient is invariant to `F₀`.**
   `F_Δ = F_rel·F_rel,n⁻¹ = (F·F₀⁻¹)(F_n·F₀⁻¹)⁻¹ = F·F_n⁻¹`. So `F₀` changes **only**
   the stress-free initial elastic state at birth; the step-to-step response (and
   hence the inner's `bᵉ` evolution) is identical with or without it.
2. **Cauchy stress is reference-independent.** The inner returns `σ = τ_rel/J_rel`,
   which **is** the true physical Cauchy stress (Kirchhoff `τ_rel` is relative to
   the birth config, `J_rel = det F_rel = det F/det F₀`). The element integrates
   `∫bᵀσ dv` on the **current** configuration with its own `det F`. So there is **no
   `P = P_rel·F₀⁻ᵀ` push-back** (a PK1 artifact) and **no `J₀` integration weight**
   (the current-config integration already carries the right measure). The inner's
   `getJ()` returns `J_rel`, which the element does **not** use for assembly.

Net: the wrapper is a near pass-through. The one genuine consistency point — the
`J_rel` inside the spatial tangent `c = (1/2J_rel)[D:L:B]` matching the element's
`dv = det F · dV` — is enforced by the element's **FD-tangent gate**, the same
arbiter that validated the base finite seam and F-bar.

**Trigger = auto-capture on the first `setTrialF`.** That first call is the birth
gradient: `Domain::addElement → setDomain → update` evaluates the element at the
**committed** deformed state, so the wrapper's first `F` is `F_birth`. Each mid-stage
append gets a fresh wrapper instance ⇒ **multi-lift staged construction falls out for
free** with no `InitialStateAnalysis`-style global toggle. Degenerates to identity at
`t = 0` (`F₀ = I`). Opt-out `-noInitF` (pass-through ≡ bare inner); explicit `-F0` for
a known birth gradient (kept across `revertToStart`; auto-capture re-arms instead);
`setResponse("F0")` reports the captured gradient.

**Why seam-2 over seam-1 (element captures `F₀`).** The element already exposes the
F-interface (`FiniteStrainNDMaterial::setTrialF`), so the wrapper needs no element
change and is reusable across any future F-based solid. Seam-1 would have edited
`LadrunoBrick`'s finite path for no gain.

**Validation** (`tests/test_initDefGrad_material.py`, all through `LadrunoBrick
-geom finite` since the material has no Python `setTrialF`): `-noInitF` pass-through ≡
bare inner (stress + force); explicit `-F0=G`, impose `F=G` → stress-free birth;
**post-birth rigid rotation** `F=Q·G` → zero stress (the additive-`InitStrain`
failure mode); reduce-to-relative (`InitDefGrad(-F0=G)` at `F=M` Cauchy σ == bare
inner at `F=M·G⁻¹`); **wrapped consistent tangent == FD of resisting force** (the
wrapper preserves Newton); **genuine two-phase staged construction** (deform+commit a
holder brick, append an auto-capture `InitDefGrad` brick on the same nodes → born with
~zero force, holder unchanged, reported `F₀` == the deformed gradient).

**Resolved open questions** (from the list below): *capture trigger* → auto-capture
(no explicit command); *which elements first* → `LadrunoBrick -geom finite` via the
material seam (element-agnostic); *strain measure / material API* → the inner's
`FiniteStrainNDMaterial` contract (Cauchy σ + spatial tangent, no push-back);
*recorder semantics* → default **relative** Hencky (the appended element's physical
state) via the inner, with the birth `F₀` exposed; total-strain channel deferred.
*class tags / ledger* → 33013, row added to [[LEDGER_implementations]].

**Deferred (v2):** a `totalStrain` recorder channel (Hencky from `X_ref`); the small-
strain additive companion is unchanged (`InitStrain` stays correct there); an explicit
re-birth/deactivation command (out of scope — OpenSees has no stress-relieving element
deactivation anyway, §7 of [[constraints_reference_position]]).

## Small-strain companion — `StagedStrainNDMaterial` (`StagedStrain`)

Shipped alongside the finite wrapper as the **small-strain** member of the `Staged*`
family (`StagedStrain` small / `StagedDefGrad` = `InitDefGrad` finite). classTag
`ND_TAG_StagedStrainNDMaterial = 33014`.

**Why a new class rather than the upstream `InitStrain`.** `InitStrainNDMaterial`
(Petracca 2024) is a fine *fixed, additive prestrain* but bites staged use three ways:
(1) **3D-only** — `getOrder()` is hardcoded to 6 and it `exit(-1)`s (kernel kill) on a
non-3D inner, so it is dead for `FourNodeQuad` / plane-strain / plane-stress /
axisymmetric; (2) **fixed `ε0`** — no auto-capture, so one tag can't birth a field of
elements with different birth strains; (3) it **adds** a prestrain rather than
**subtracting** the captured birth strain. `StagedStrain` fixes all three: order-general
(`ε0` sized to `inner->getOrder()`, `getCopy(type)` adapts the inner to the element's
view), auto-capturing, and graceful. `InitStrain` is left untouched for prestrain.

**Mechanic (additive — exact for small strain).** Capture `ε0` at the first
`setTrialStrain` (= birth), then feed the inner `ε_rel = ε − ε0`. At birth `ε_rel = 0`
⇒ born **genuinely virgin** (zero stress *and* zero plastic history). The tangent is an
**exact passthrough** (`ε0` constant ⇒ `∂σ/∂ε = ∂σ/∂(ε−ε0)`), so — unlike the finite
wrapper — **no FD-tangent gate is needed**. Additive ⇒ valid for small strain only;
finite uses `StagedDefGrad`.

**Composability (the load-bearing design choice).** `StagedStrain` and the upstream
`InitStrain` compose by nesting, with **`StagedStrain` OUTERMOST** (element-facing):

```
element → StagedStrain( InitStrain( realMaterial, ε_pre ) )
```

At birth `StagedStrain` feeds `InitStrain` `ε_rel = 0`, `InitStrain` adds `ε_pre`, the
real material sees `ε_pre` ⇒ born carrying **exactly the prestress `σ(ε_pre)`**, no
inherited geometric stress. Nest the other way and the prestrain is captured into `ε0`
and **cancels** — so the order is a hard rule. This is *why* the two stay separate
composable wrappers instead of one `-prestrain`-mode class: staged-birth, prestrain, and
any future offset (`InitStress`, thermal) all come from nesting.

**Validation** (`tests/test_stagedStrain_material.py`, through `FourNodeQuad` 2D +
`stdBrick` 3D, since the material has no Python `setTrialStrain`): `-noInit` ≡ bare
(2D+3D); explicit-`eps0` reduce-to-relative; **two-phase staged stress-free birth in
2D *and* 3D** (also the dimensional-generality proof — order-3 capture for the quad,
order-6 for the brick); PlaneStress birth; **J2 inner born virgin** (no inherited
plastic strain); **composition `StagedStrain∘InitStrain`** → born stress independent of
the birth deformation and equal to `σ(ε_pre)`; **graceful failure** where `InitStrain`
`exit(-1)`s.

### Deferred — finite-strain prestrain (NOT in scope)

Prestrain/prestress *composes for free in small strain* (`StagedStrain ∘ InitStrain`,
both `setTrialStrain`-driven), but **there is no free path in finite strain**, because
the finite chain is `setTrialF`-driven and `InitStress`/`InitStrain` are not
`FiniteStrainNDMaterial`s — they cannot nest in the F-chain. A finite prestrain
therefore needs a **dedicated F-driven wrapper** (future work), with a genuine design
fork:

1. **Pre-deformation `F_pre`** (recommended) — feed the inner `F_rel·F_pre`, so at
   stress-free birth it sees `F_pre` → `σ(F_pre)`. The clean multiplicative analog of
   `InitStrain`'s additive `ε_pre`; **objective by construction**; composes as
   `StagedDefGrad(FinitePrestrain(...))`.
2. **Initial Cauchy stress `σ_pre`** (the `InitStress` analog) — additive on stress, but
   a *fixed* `σ_pre` does **not** co-rotate ⇒ **non-objective under large rotation**
   (the exact failure mode that motivated `InitDefGrad`); a correct version must
   co-rotate the prestress. Harder; only if a consumer needs it.

(An unvetted partial path exists — `LogStrain(InitStress(real))`, the small-strain
`InitStress` at the bottom of the Hencky stack — but its large-rotation objectivity is
unverified, which is precisely why finite prestrain is its own task, not a freebie.)

## See also

- [[constraints_reference_position]] — the SP / MP / element trichotomy and the
  small-vs-finite-strain caveat that leads here.
- [[solid_transformation_wrapper]] — the `linear`/`corot`/`finite` kinematics layer
  this feature rides on (`F₀` is a capture-at-birth hook on its `finite` ledger).
- [[09_ladruno_brick]] — the first solid element that would consume both.
- [[LEDGER_quirks]] — the mid-stage append gotcha for continuum elements.
