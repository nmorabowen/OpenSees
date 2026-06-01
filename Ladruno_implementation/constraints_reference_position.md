---
title: Constraints, elements & the reference position in staged analysis
project: Ladruno
tags:
  - reference
  - constraints
  - elements
  - staged-analysis
  - opensees-internals
related:
  - "[[LEDGER_quirks]]"
  - "[[staged_deformation_gradiend]]"
---

# Constraints, elements & the reference position

**What configuration does a model object act in when you add it *after* the model has
already deformed?** This is the central question for staged / sequential analysis
(deform under stage 1, hold the load, then change the BCs *or* add material to part of
the model). The answer differs sharply between single-point constraints, multi-point
constraints, and elements — and getting it wrong produces spurious forces, unexpected
snap-back, or a member born pre-stressed.

> **One-liner:** `fix`/`sp` are **absolute** (reference-frame): adding one mid-stage
> *moves* the node. `equalDOF`/`rigid*` ties are **incremental** (deformed-frame):
> adding one mid-stage *doesn't* move anything. **Elements are a mixed bag** — line
> elements (truss, beam-column) are born stress-free at the deformed shape by default;
> most continuum elements inherit the accumulated strain. Same "capture state at
> add-time" idea underneath everywhere, but **the defaults differ per object.**

All file/line references below are into this fork's `SRC/` tree.

---

## 1. The one displacement frame

OpenSees stores nodal **displacement `u` measured from the original mesh `X_ref` at
t = 0**. There is exactly one displacement frame and **it never re-zeros** — not at a
stage boundary, not after a `wipeAnalysis`, not ever. Everything downstream is
consistent with this:

```
position = X_ref + u
```

Element internal forces are computed from `X_ref + u` (including corotational and
P-Δ geometric transforms), so the element and the constraint live in the *same*
frame. The subtlety is never a frame mismatch — it is **what value the constraint
pins, and relative to which configuration.**

A useful gotcha-phrasing: *"a constraint fixes deformation, not position"* is true but
misleading, because deformation is itself measured from the undeformed mesh. For a
single-point constraint, `u = 0` and "back at the original position" are the **same
statement** (`u = 0 ⟺ position = X_ref`).

---

## 2. Single-point constraints (`fix`, `sp`) — absolute / reference-frame

### Behavior

An `SP_Constraint` prescribes the **total value** of one displacement DOF:
`u = value`. `fix` is just `value = 0`. Adding it mid-stage to a node that has
already drifted to `u = d` drives that DOF back toward the prescribed value on the
next `analyze` → the node is pulled toward its reference position, and the attached
elements pick up spurious strains/forces.

### Source trail

- **`fix` builds a homogeneous SP with value 0.**
  `OPS_HomogeneousBC` → `new SP_Constraint(node, dof, 0.0, true)`
  ([SP_Constraint.cpp:80](../SRC/domain/constraints/SP_Constraint.cpp)).

- **The current displacement *is* captured at add-time** — but only stored:
  `SP_Constraint::setDomain()` records `initialValue = U(dofNumber)`
  ([SP_Constraint.cpp:380](../SRC/domain/constraints/SP_Constraint.cpp)).

- **Whether that captured value is used is gated by `retZeroInitValue`**, which
  defaults to `true` → `getInitialValue()` returns **0**, i.e. the captured value is
  ignored ([SP_Constraint.cpp:317](../SRC/domain/constraints/SP_Constraint.cpp);
  flag declared with its semantics at
  [SP_Constraint.h:81](../SRC/domain/constraints/SP_Constraint.h)).

- **All three handlers subtract `getInitialValue()`** — but since it is 0 by default,
  the subtraction is a no-op and the *absolute* value is enforced:
  - Penalty: `resid = alpha * (constraint - (nodeDisp - initialValue))`
    ([PenaltySP_FE.cpp:139](../SRC/analysis/fe_ele/penalty/PenaltySP_FE.cpp)).
  - Lagrange: same `initialValue` read
    ([LagrangeSP_FE.cpp:143](../SRC/analysis/fe_ele/lagrange/LagrangeSP_FE.cpp)).
  - Transformation: `setTrialDisp(value + initial_value, i)` under
    `#define TRANSF_INCREMENTAL_SP`
    ([TransformationDOF_Group.cpp:1055-1056](../SRC/analysis/dof_grp/TransformationDOF_Group.cpp);
    macro at [TransformationDOF_Group.h:44](../SRC/analysis/dof_grp/TransformationDOF_Group.h)).

- **The "stay-in-place" path is not cleanly reachable from script.** Both `fix` and
  `sp` default `retZeroInitValue = true`, and even the `sp -subtractInit` flag sets it
  to `true` again ([OpenSeesPatternCommands.cpp:1056,1065,1092](../SRC/interpreter/OpenSeesPatternCommands.cpp)).
  So the incremental SP branch is compiled in but effectively dormant.

**Net:** with default commands, a mid-stage `fix`/`sp` enforces the *absolute*
displacement → **moves the node to the reference frame.**

---

## 3. Multi-point constraints (`equalDOF`, `rigidLink`, `rigidDiaphragm`) — incremental / deformed-frame

### Behavior

An `MP_Constraint` prescribes a **linear relation** between a constrained node's DOFs
and a retained node's DOFs (`equalDOF` ⇒ identity). Added mid-stage, it does **not**
snap the constrained node onto the retained node. It captures the current relative
configuration and ties only the **future increments**, so the offset present at
tie-time is preserved and the initial constraint force is **zero**.

This is the "install at the current deformed state" behavior that `fix` lacks — and
for MP it is the default, with no flag required.

### Source trail

- **Both nodes' current displacements are captured at add-time**, unconditionally:
  `MP_Constraint::setDomain()` stores `Uc0` (constrained) and `Ur0` (retained)
  ([MP_Constraint.cpp:294-313](../SRC/domain/constraints/MP_Constraint.cpp); getters at
  [MP_Constraint.cpp:390-398](../SRC/domain/constraints/MP_Constraint.cpp)).

- **Every MP-capable handler enforces the relation on offset-removed displacements**
  — there is **no** `retZeroInitValue`-style flag; the subtraction is always on:
  - Penalty: `UU(i) = Uc(cdof) - Uc0(i)`, `UU(...) = Ur(rdof) - Ur0(i)`
    ([PenaltyMP_FE.cpp:222-238](../SRC/analysis/fe_ele/penalty/PenaltyMP_FE.cpp)).
    Equilibrium is therefore `(Uc - Uc0) = C·(Ur - Ur0)`, **not** `Uc = C·Ur`.
  - Lagrange: identical `Uc - Uc0` / `Ur - Ur0`
    ([LagrangeMP_FE.cpp:243-260](../SRC/analysis/fe_ele/lagrange/LagrangeMP_FE.cpp)).
  - Transformation: under `#define TRANSF_INCREMENTAL_MP` it transforms only the
    **increment** (`modUnbalance -= modTrialDispOld`) and applies it with
    `incrTrialDisp`, never overwriting the standing offset
    ([TransformationDOF_Group.cpp:518-525,552-560](../SRC/analysis/dof_grp/TransformationDOF_Group.cpp);
    macro at [TransformationDOF_Group.h:45](../SRC/analysis/dof_grp/TransformationDOF_Group.h)).

- **At tie-time** `Uc = Uc0` and `Ur = Ur0`, so the constraint is satisfied with zero
  residual → **zero force, zero jump.**

**Net:** a mid-stage `equalDOF`/rigid tie enforces the *increment* → **preserves the
deformed offset, no snap.**

---

## 4. Appending elements mid-stage — element-dependent

### Behavior

When you append an element after the model has deformed, **what it treats as its
stress-free reference is decided per-element**, not globally:

- **Line elements** — `Truss`, `CorotTruss`, and beam-columns (via their geometric
  transform) — are born **stress-free at the deformed geometry** by default. They
  snapshot the committed nodal displacement at birth and subtract it.
- **Most continuum elements** — `FourNodeQuad`, `stdBrick`, tets, etc. — have **no**
  such hook. They compute `ε = B·u` straight from the node pointers, so an element
  added at `u = d` is born referenced to `t = 0` and **instantly carries the
  accumulated displacement as spurious strain.**

There is no single global answer the way there is for SP and MP.

### Mechanism

`Domain::addElement` initializes the element *immediately* at add-time, against the
current committed state ([Domain.cpp:473-485](../SRC/domain/domain/Domain.cpp)):

```cpp
element->setDomain(this);   // element captures its reference geometry NOW
element->update();          // ...and computes its first state NOW
...
this->domainChange();       // analysis rebuilds the FE/DOF model next analyze()
```

So whatever the element snapshots in `setDomain()` is taken against the domain's
**committed** state at the moment of the call.

### Source trail — the capture path (line elements)

- **`Truss`** reads committed nodal disp (`getDisp()`) in `setDomain` and, when
  `useInitialDisp` is on and the disp is nonzero, stores it as `initialDisp` *and folds
  it into the reference length* `L` ([Truss.cpp:407-419](../SRC/element/truss/Truss.cpp)),
  then subtracts it in `computeCurrentStrain`
  ([Truss.cpp:1132-1137](../SRC/element/truss/Truss.cpp)). **`useInitialDisp` defaults to
  `true`** ([Truss.cpp:92](../SRC/element/truss/Truss.cpp) — *"default is the previous
  behavior"*). The capture only fires when the committed disp is nonzero, so t = 0
  models pay nothing.
- **Beam-columns** do the same one layer down, in the geometric transform:
  `LinearCrdTransf3d` (and PDelta/Corot) capture `nodeIInitialDisp`/`nodeJInitialDisp`
  from `getDisp()` when nonzero
  ([LinearCrdTransf3d.cpp:211-227](../SRC/coordTransformation/LinearCrdTransf3d.cpp), guarded
  by `initialDispChecked`) and subtract them when forming the basic deformations.

### Source trail — the no-capture path (continuum)

- **`FourNodeQuad`** builds `B` from `getCrds()` (reference coords,
  [FourNodeQuad.cpp:1647-1650](../SRC/element/fourNodeQuad/FourNodeQuad.cpp)) and strain
  from `getTrialDisp()` directly
  ([FourNodeQuad.cpp:578-581](../SRC/element/fourNodeQuad/FourNodeQuad.cpp)) — no
  initial-disp path. Appended mid-stage, it inherits the full accumulated strain.

### Fixing the continuum case — two wrapper seams

There's no `CrdTransf`-style seam to bolt onto a continuum element (continuum elements
have no global→basic transform layer). The re-referencing has to sit elsewhere:

1. **Material-level (idiomatic; already in-tree).** Wrap the `nDMaterial` and offset in
   strain space. `InitialStateAnalysisWrapper` does `mStrain = strain_from_element +
   mEpsilon_o` ([InitialStateAnalysisWrapper.cpp:184](../SRC/material/nD/UWmaterials/InitialStateAnalysisWrapper.cpp));
   `InitStrainNDMaterial` / `InitStressNDMaterial` offer fixed offsets. For stress-free
   birth, set the material's zero-stress point to the birth strain `ε₀ = B·u₀` so
   `σ(ε₀) = 0`. The element still *computes* the full strain; only the **stress** is
   neutralized.
2. **Element-level proxy-node wrapper.** A generic `Element` that captures `u₀` at
   `setDomain` and feeds the inner element proxy `Node`s reporting `u − u₀`. Truer to
   the "transformation wrapper" idea, but heavier: you must forward mass, damping,
   recorders, response, and the inner element pulls disp from its `Node*` pointers
   directly (so it needs *proxy nodes*, not a thin shell). OpenSees ships no generic
   element-wrapper base.

### The small-strain vs finite-strain caveat

| | Small-strain continuum (`FourNodeQuad`, `stdBrick`) | Finite-strain / large-displacement |
|---|---|---|
| `ε = B·u` is… | **linear** in `u`, `B` fixed from ref coords | **nonlinear** (`F = I + ∂u/∂X`) |
| Subtract `u₀` / offset `ε₀` | **exact** stress-free birth ✓ | **wrong** — additive, fixed-frame; not objective, doesn't compose multiplicatively |
| Correct fix | (offset is fine) | re-reference the configuration (`F_rel = F·F₀⁻¹`) — see **[[staged_deformation_gradiend]]** |

The additive material/displacement offset is **exact for small strain** but breaks for
a large-deformation element: it is non-objective (a post-birth rigid rotation produces
spurious strain) and it adds in strain space where finite strain requires multiplicative
composition of deformation gradients. The clean within-`X_ref` answer — store a per-GP
birth deformation gradient `F₀` and feed the material `F_rel = F·F₀⁻¹` — is developed in
its own design note: **[[staged_deformation_gradiend]]**.

---

## 5. The asymmetry at a glance

| | **SP** (`fix`, `sp`) | **MP** (`equalDOF`, `rigid*`) | **Element** (line: truss / beam-col) | **Element** (continuum) |
|---|---|---|---|---|
| What it sets | a **value** (`u = value`) | a **relation** (`Uc = C·Ur`) | its stress-free reference geometry | its stress-free reference geometry |
| State captured at add-time | `initialValue` | `Uc0`, `Ur0` | `initialDisp` / `nodeI*InitialDisp` | **none** |
| Used by default? | **No** (gated by `retZeroInitValue`) | **Yes** (always) | **Yes** (`useInitialDisp` default true) | **n/a** |
| Frame of effect | **reference** (absolute) | **deformed at tie-time** | **deformed at birth** | **reference** (t = 0) |
| Adding it mid-stage | **moves** node → spurious force | **no move** → zero force | **stress-free birth** | **inherits strain** → spurious stress |

Same underlying idea everywhere — *"snapshot the state when the object enters the
domain"* — but the **defaults differ per object class.** MP and the line elements were
hardened for staged construction; SP kept its legacy absolute-value default, and the
continuum elements never got the hook at all.

---

## 6. Practical recipes

**Install a ground support that holds the current deformed position (zero initial
force), via SP:** prescribe the *current* displacement rather than `fix`:

```python
d = ops.nodeDisp(n, dof)        # current total displacement
ops.sp(n, dof, d, '-const')     # target = current → no jump
# needs an active pattern, or pass -pattern <tag>
```

**Genuinely return a DOF to its t = 0 position:** `fix` is correct and the resulting
forces are physical, not spurious.

**Catch the deformed structure against another node with zero initial force:** just
`equalDOF` (or a rigid tie) to that node — the offset-preservation is automatic.

**Any sudden BC change in a transient stage injects an impulse** (energy spike) even
when the value is right; ramp the prescribed displacement with a `timeSeries` rather
than stepping it.

**Appending an element and wanting stress-free birth:** a `Truss`/beam-column gives it
for free (commit the deformation first). For a **continuum** element, neutralize the
stress at the constitutive level (`InitStrainNDMaterial` / `InitialStateAnalysisWrapper`
with the birth strain) — valid for small strain only; for large deformation see
**[[staged_deformation_gradiend]]**.

---

## 7. Caveats

- **Handler scope.** The offset-preservation for MP holds for the **MP-capable
  handlers**: Transformation, Penalty, Lagrange. The **Plain** handler is not the one
  to use for nontrivial MP staged ties.
- **Vintage.** The MP `Uc0`/`Ur0` machinery (and `TRANSF_INCREMENTAL_*`) is relatively
  recent upstream (staged-construction hardening). On **older** OpenSees binaries
  `equalDOF` mid-analysis *did* snap the constrained node onto the retained one — so
  old forum lore describing a jump reflects the legacy behavior, not this build.
- **`isConstant` / load patterns.** A non-constant `sp` inside a `pattern` ramps as
  `valueC = loadFactor * valueR` ([SP_Constraint.cpp:331-338](../SRC/domain/constraints/SP_Constraint.cpp));
  `-const` pins it. This is orthogonal to the reference-vs-deformed question but
  matters when scripting imposed support motion.
- **Element capture is committed disp, not trial.** The `initialDisp`/`nodeI*InitialDisp`
  snapshot is taken from `getDisp()` (committed) — commit the deformation before
  appending, or the birth reference is wrong. `-useInitialDisp 0` turns the capture off
  for `Truss` (forces reference-frame birth → inherits strain).
- **Dynamics.** A newly appended element's mass/stiffness enters the system
  discontinuously — same impulse caveat as a sudden BC change. And this is stress-free
  *activation* only; OpenSees has no general element *deactivation* with stress relief
  (`remove element` just deletes it).

---

## 8. See also

- [[LEDGER_quirks]] — the condensed gotcha entries:
  *"A `fix`/`sp` added mid-analysis snaps the node to the REFERENCE frame…"* and its
  companion *"…but `equalDOF`/MP constraints are the OPPOSITE…"*.
- [[staged_deformation_gradiend]] — the clean within-`X_ref` design for true
  stress-free activation of large-deformation continuum elements (per-GP `F₀`,
  `F_rel = F·F₀⁻¹`). Forward-looking; implementation TBD.
