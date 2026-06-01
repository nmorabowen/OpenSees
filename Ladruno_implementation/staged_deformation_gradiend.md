---
title: Staged deformation gradient — stress-free activation of large-deformation continuum elements
project: Ladruno
status: design note (forward-looking; implementation TBD)
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

> **Status:** design note. No code yet — this captures the theory and the intended
> shape of a future Ladruno feature. See [[constraints_reference_position]] §4 for the
> context (what happens when you append elements mid-stage) that motivates it.

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

## See also

- [[constraints_reference_position]] — the SP / MP / element trichotomy and the
  small-vs-finite-strain caveat that leads here.
- [[solid_transformation_wrapper]] — the `linear`/`corot`/`finite` kinematics layer
  this feature rides on (`F₀` is a capture-at-birth hook on its `finite` ledger).
- [[09_ladruno_brick]] — the first solid element that would consume both.
- [[LEDGER_quirks]] — the mid-stage append gotcha for continuum elements.
