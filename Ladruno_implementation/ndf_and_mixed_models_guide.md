---
title: NDF and Mixed-DOF Models — How degrees-of-freedom-per-node work in OpenSees
project: Ladruno
status: living-reference
priority: high
owner: nmora
tags:
  - reference
  - guide
  - ndf
  - dof
  - mixed-models
  - opensees-internals
  - constraints
  - mass
  - loads
  - explicit
related:
  - "[[LEDGER_quirks]]"
  - "[[23_ladruno_embedded_node_adr]]"
  - "[[24_ladruno_coupling_constraints_adr]]"
  - "[[constraints_reference_position]]"
  - "[[ladruno_continuum_elements_guide]]"
aliases:
  - "ndf guide"
  - "mixed ndf"
  - "degrees of freedom per node"
---

# NDF and Mixed-DOF Models in OpenSees

> [!abstract] What this document is
> A single reference for **how `ndf` (degrees-of-freedom-per-node) is defined,
> stored, and consumed** across OpenSees — and exactly how far you can push
> **mixed-ndf models** (a single domain holding nodes with different dof counts:
> e.g. 3-dof solid nodes next to 6-dof shell/beam nodes, or u-p nodes with an
> extra pressure dof). It traces ndf through every layer that touches it:
> the model builder, the `Node`, the solver (DOF_Group / numberer), **element
> formulations**, **boundary conditions & constraints**, **loads**, and **mass**.
>
> All file/line references are into this fork's `SRC/` tree. Where a fork-authored
> class is the right tool (e.g. `LadrunoEmbeddedNode`), it is linked, but the
> body is upstream-faithful OpenSees behavior unless marked otherwise.

> [!tip] One-liner
> **`ndf` is a per-`Node` property, not a global model property.** The builder's
> `-ndf` is only a *default stamped onto new nodes at creation*. Everything
> downstream is sized off each node's own `numberDOF`. As a result, mixed-ndf is
> **native and safe at the node/solver layer**, and **brittle-to-unsupported at
> the element/rigid-constraint layer**. That asymmetry is the entire story —
> whatever breaks in a mixed model breaks at an *element or rigid tie that spans
> two dissimilar nodes*, never in the numbering or the solve.

---

## 0. The mental model (read this first)

```
  model -ndf N        ──>  a DEFAULT integer held in the interpreter command object
        │                  (OpenSeesCommands::ndf / PythonModelBuilder::ndf).
        │                  Mutable. Re-issuing `model` overwrites it. No domain check.
        ▼
  node tag x y z      ──>  stamps the CURRENT default onto Node::numberDOF
  node tag x y z -ndf K     (… unless -ndf K overrides it for THIS node)
        │
        ▼
  Node::numberDOF     ──>  THE source of truth. Sizes trial/commit disp,vel,accel,
        │                  mass (ndf×ndf), unbalanced load (ndf), and the node's
        │                  DOF_Group (ndf equation slots).
        ▼
  DOF_Group(node)     ──>  sized to node->getNumberDOF(). Numberer assigns eqn #s
  PlainNumberer            SEQUENTIALLY across DOF_Groups of DIFFERING sizes — it
        │                  does not know or care that the model is mixed.
        ▼
  SOLVE               ──>  fully ndf-agnostic. Never the failure point.

  ── the danger is ELSEWHERE ──
  Element::setDomain  ──>  each element ASSUMES a per-node dof count. Reaction to a
                           node whose ndf differs ranges from clean exit(-1) to
                           SILENT memory corruption (see §3).
  rigid* / equalDOF   ──>  most REQUIRE ndf parity on both nodes (see §5).
```

**Takeaway:** to reason about a mixed model, forget "the model's ndf." Ask, for each
node, *what is its `numberDOF`*, and for each element/tie, *what does it assume about
the nodes it spans*.

---

## 1. Where ndf is born and stored

### 1.1 The global default (mutable, per-builder)

`model basic -ndm M -ndf N` parses and stores the default:

- `OPS_model()` in `SRC/interpreter/OpenSeesCommands.cpp:1350-1411` parses `-ndm`/`-ndf`,
  then calls `cmds->setNDF(ndf)` / `cmds->setNDM(ndm)`.
- Default when `-ndf` omitted: `ndm 1→ndf 1`, `2→3`, `3→6`
  (`OpenSeesCommands.cpp:1380-1402`).
- Python path: `ops_Model()` in `SRC/interpreter/PythonModelBuilder.cpp:124-174`
  → `pBuilder.set(ndm,ndf)`.
- **It is freely re-settable mid-model.** `setNDF` is an unconditional setter
  (`OpenSeesCommands.h:80-81`); `PythonModelBuilder::set` overwrites unconditionally
  (`PythonModelBuilder.h:64`). Re-issuing `model` changes the default for *subsequent*
  nodes; **nodes already built keep their original ndf.** No validation, no warning.

> [!warning] This re-settability is one legitimate way to build a mixed model
> (set `-ndf 3`, create the solid nodes, set `-ndf 6`, create the shell nodes), but
> it is also a silent footgun: a stray second `model` call retroactively changes the
> default and you get a mix you didn't intend, with zero diagnostic.

### 1.2 The per-node value (the real source of truth)

- Tcl: `TclCommand_addNode` reads the builder default
  (`SRC/modelbuilder/tcl/TclModelBuilder.cpp:1008-1009`) and honors a per-node
  `-ndf` override (`TclModelBuilder.cpp:1081-1099`).
- Python / OPS API: `OPS_Node()` in `SRC/domain/node/Node.cpp:66-221` reads
  `OPS_GetNDF()` and honors `-ndf`/`-NDF` (`Node.cpp:168-189`).
- Stored in `Node::numberDOF` (`SRC/domain/node/Node.h:187`), set in every Node
  constructor (`Node.cpp:272-370`), returned by `getNumberDOF()` (`Node.cpp:548-552`).
- All state vectors are sized off it: `createDisp()` allocates `4*numberDOF` and
  builds `trial/commit/incr/incrDelta` Vectors of length `numberDOF`
  (`Node.cpp:1815-1839`); vel/accel similarly `2*numberDOF`.

### 1.3 `ndm` vs `ndf` — independent

`ndm` (spatial dimension) sizes the coordinate vector `Crd`; `ndf` sizes the state
vectors. They are stored and validated separately. ndm=2 with ndf=2 (plane solid),
ndf=3 (plane frame), and ndf=4 (plane u-p) are all legal and common.

**Per-node ndf override syntax:**

```python
ops.model('basic', '-ndm', 3, '-ndf', 3)      # default for solids
ops.node(1, 0.0, 0.0, 0.0)                     # ndf = 3
ops.node(2, 1.0, 0.0, 0.0, '-ndf', 6)          # ndf = 6 for THIS node (e.g. a shell node)
ops.node(3, 0.0, 1.0, 0.0, '-ndf', 4)          # ndf = 4 (u-p)
```

```tcl
model basic -ndm 3 -ndf 3
node 1 0 0 0
node 2 1 0 0 -ndf 6
```

### 1.4 ndf is **immutable after node creation** — it is a build-time decision

A node's ndf is fixed the instant the node is constructed. `Node::numberDOF` is
assigned in exactly two places: the constructors (`Node.cpp:272-370`, from the
`ndf` you pass) and `Node::recvSelf` (`Node.cpp:1568`) — and the latter is
*deserialization*, i.e. building a **brand-new** Node from a channel/database, not
mutating a live one. There is **no public `setNumberDOF`, no Tcl/Python command,
and no analysis hook** that changes a live node's ndf (`Node.h:70` exposes only the
`getNumberDOF` getter; the field is private). It *has* to be this way — at
construction the node sizes its disp/vel/accel vectors, its `ndf×ndf` mass matrix,
and its `DOF_Group` all off `numberDOF`; changing it live would leave every one of
those mis-sized.

The only way to "change" a node's ndf is to **destroy and rebuild**:
`remove node <tag>` then re-create it with the new ndf. But `Domain::removeNode`
(`Domain.cpp:1147`) just pulls the node out and hands back the pointer — it does
**not** delete the elements/constraints/loads attached to it, which then dangle. So
in practice "change a node's ndf" = tear down and rebuild that node *and everything
connected to it*, then let `domainChange()` re-number.

> [!important] Consequence for mixed models
> You must decide each node's ndf **up front**. You cannot promote an ndf=2 solid
> node to ndf=3 later to attach a frame to it (see the staged-construction worked
> example, §8). This is *the* reason cross-formulation staging needs a separate
> node + a tie rather than an in-place upgrade.

---

## 2. The solver layer is genuinely ndf-agnostic (the good news)

- Each node gets exactly **one** `DOF_Group`, sized to *that node's*
  `getNumberDOF()` (`SRC/analysis/dof_grp/DOF_Group.cpp:60-64`; it `exit(-1)`s only
  if a node has `numDOF <= 0`).
- `PlainHandler::handle()` creates one DOF_Group per node with **no uniformity
  check** and loops over `id.Size()` = that node's ndf
  (`SRC/analysis/handler/PlainHandler.cpp:101-141`).
- `TransformationConstraintHandler` sizes from `nodPtr->getNumberDOF()`
  (`TransformationConstraintHandler.cpp:183`).
- `PlainNumberer::numberDOF()` assigns equation numbers **sequentially across
  DOF_Groups of differing sizes** (`SRC/analysis/numberer/PlainNumberer.cpp:98-115`).
  It iterates `theID.Size()` per group; heterogeneity is invisible to it.

> [!success] Consequence
> A domain with 3-dof solid nodes interleaved with 6-dof shell nodes numbers,
> assembles, and solves with no complaint. **Mixing ndf never breaks DOF numbering
> or the linear solve.** Whatever breaks, breaks at an element or constraint that
> *spans* dissimilar nodes (§3, §5).

---

## 3. Element ↔ node ndf — the danger zone

### 3.1 How the mapping actually works

An element advertises only a **total** `getNumDOF()` (sum over its nodes), *not* a
per-node breakdown (`SRC/element/Element.h:57-61`; e.g. `Truss::getNumDOF` returns a
stored `numDOF`). At assembly, `FE_Element::setID()` walks each attached node's
DOF_Group and **greedily pulls ALL of that node's equation numbers, in order,** into
the element's `myID` (`SRC/analysis/fe_ele/FE_Element.cpp:252-284`):

```cpp
for (int i=0; i<numGrps; i++) {            // each node
  const ID &theDOFid = dofPtr->getID();
  for (int j=0; j<theDOFid.Size(); j++)    // ALL of the node's dofs, positionally
    if (current < numDOF) myID(current++) = theDOFid(j);
}
```

The element's local K/M/R is then mapped **positionally** onto those equations. There
is **no reconciliation** between what the element assumed per node and what the node
actually carries. If the element's `numDOF` and the sum of the nodes' real ndf
disagree, you get truncation (too-small numDOF) or out-of-bounds reads/writes
(too-large), depending on the element.

### 3.2 The four contracts elements implement

| Contract | Elements (examples) | On node-ndf mismatch | File:line |
|---|---|---|---|
| **Hard `exit(-1)`** | `ElasticBeam2d`/`3d` (require ndf=3 / 6) | Kills the process with a message | `ElasticBeam3d.cpp:567-580` |
| **Warn-but-continue** | `ShellMITC4` (needs 6) | Prints *"NEEDS 6 dof — GARBAGE RESULTS or SEGMENTATION FAULT WILL FOLLOW"* then proceeds and reads OOB | `ShellMITC4.cpp:359-375` |
| **Silent bail** | `FourNodeQuad` (needs 2) | `return`s *without* calling `setDomain` → element half-initialized, later assembly corrupts | `FourNodeQuad.cpp:474-485` |
| **No check at all** | `Brick` (hardcodes `getNumDOF()=24`, assumes 3/node) | If nodes are ndf=6, `setID` writes 48 eqns into a 24-slot map → overflow / garbage | `Brick.cpp:251-271, 319-322` |
| **Adaptive (the exception)** | `Truss`, `ZeroLength`, `ZeroLengthSection` | Dispatch on `(ndm, ndf)`; **require both nodes equal**; size to `2*ndf` | `Truss.cpp:299-386`, `ZeroLengthSection.cpp:234-247` |

Two structurally important facts:

1. **Elements consume node dofs positionally and greedily.** A 2D truss on ndf=3
   nodes silently absorbs the unused rotational dof into a 6×6 element with zero
   rotational stiffness — harmless here. The *same* greediness is what makes a `Brick`
   on ndf=6 nodes catastrophic. The behavior isn't "use the first K dofs I need"; it's
   "take everything the node has and hope the count matches what I assumed."

2. **`ZeroLength` / `ZeroLengthSection` are the only formulation-flexible family.**
   `-dir` maps each 1D material onto a chosen dof index via a transformation matrix
   `t1d`, decoupling *material count* from *node ndf* (`ZeroLength.cpp:142-164` parses
   `-dir`; `setTran1d` builds the map). This is why zeroLength is the canonical tool
   for stitching dof *meanings* together at a coincident point. Note even here the two
   nodes must have **equal** ndf — the flexibility is in *which* dofs the materials
   act on, not in bridging dissimilar nodes.

> [!danger] The core hazard
> "Mixed ndf" is safe *within a node's own element family*. It becomes unsafe the
> moment a **single element straddles two nodes of different ndf** — and the failure
> mode ranges from a clean crash (beam) to **silent memory corruption** (brick). No
> standard element auto-adapts to a per-node-varying ndf except by the
> equal-both-ends dispatch of Truss/ZeroLength.

### 3.3 Per-family quick table

| Element family | dofs/node assumed | Behavior when a node's ndf differs |
|---|---|---|
| Truss | 1/2/3 by (ndm,ndf), both ends equal | warn + fallback numDOF=2 |
| ZeroLength / ZeroLengthSection | any, both ends equal; scales `2*ndf` | warn + bail if ends differ |
| ElasticBeam2d/3d | 3 / 6 | `exit(-1)` |
| FourNodeQuad (and plane solids) | 2 | silent `return`, no setDomain |
| Brick / 3D solids | 3 (hardcoded 24 total) | **no check** → corruption if ndf≠3 |
| ShellMITC4 / shells | 6 | warn loudly, continue (UB) |

### 3.4 The CST ndf=3 *trap* — a permissive check that doesn't mean what it says

The 2D CST (`Tri31`) looks like it tolerates 3-dof nodes: `setDomain` accepts both
(`dofNd != 2 && dofNd != 3` → it only rejects if *neither*), and `getNumDOF()`
returns the **sum** of the nodes' dofs (`Tri31.cpp:450-460, 491-501`). **It does
not actually work with ndf=3 nodes.** Its stiffness is a hardwired static **6×6**
matrix (`Tri31::K(matrixData, 6, 6)`, `Tri31.cpp:51`) filled with **stride-2**
indexing (`ia += 2`, `Tri31.cpp:620-635`). So with three ndf=3 nodes,
`getNumDOF()` returns 9 while `getTangentStiff()` returns a 6×6 → the `FE_Element`
9×9 tangent **rejects** the `addMatrix` on size mismatch and the element
**contributes zero stiffness** (silent, no message). `FourNodeQuad` is blunter — it
flat-refuses anything but ndf=2 and silently bails (`FourNodeQuad.cpp:474-485`).

> [!warning] Don't build 2D continuum nodes as ndf=3 "to leave room for a frame."
> The CST's ndf=3 branch is a trap, the quad rejects it outright, and even if a
> solid *did* accept the extra dof, that rotational dof would have **zero stiffness**
> in a soil-only stage → singular `K`. Build soil nodes ndf=2; attach frames via a
> tie to *separate* ndf=3 nodes (§8).

---

## 4. The physics gap nobody auto-fixes: rotational compatibility

Solid elements (ndf=3) have **no rotational dofs and no rotational stiffness**.
Shells/beams (ndf=6) carry 3 rotations. When you connect them:

- **No code anywhere** synthesizes rotational stiffness for a solid node, and none
  auto-couples a shell's rotations to a solid's translations.
- Tie only the 3 translations and the shell/beam node **spins freely** at the
  junction — the moment leaks. OpenSees solves this *wrong* model without complaint;
  it is a *modeling* incompleteness, not a software error.

This is exactly the problem the fork's **`LadrunoEmbeddedNode`** ([[23_ladruno_embedded_node_adr]])
solves the right way. At `setDomain` it reads **each node's own ndf**, lays out
per-node dof offsets, always couples the first `ndm` translational dofs, and
*optionally* reconstructs:

- **rotation** (`-rot`): `θ_c = skew(∇u_host)` from the host shape-function gradient —
  enabled only if the constrained node actually has `ndf ≥ ndm+nrot`;
- **pressure / u-p** (`-pressure`): enabled only if *all* nodes have `ndf ≥ ndm+1`.

If a participating node lacks the required dofs it **degrades to translation-only with
a warning** rather than corrupting (`SRC/element/ladrunoEmbeddedNode/LadrunoEmbeddedNode.cpp:139-210`).
Honest limit (ADR 23, UR-4): on a CST / TET4 host, `∂N/∂x` is element-wise constant, so
the rotation tie collapses to a single element-wide rigid spin.

> [!tip] Pattern to replicate
> Every robust mixed-ndf feature in the fork wins by **re-reading `getNumberDOF()`
> per node and coupling only the dofs that physically exist.** Any future
> cross-formulation tie should do the same — never assume a uniform model ndf.

See also [[24_ladruno_coupling_constraints_adr]] for the broader coupling-constraint
strategy.

---

## 5. Boundary conditions & constraints

### 5.1 Single-point (`fix`, `fixX/Y/Z`, `sp`) — two-stage validation

- `fix` infers count from its argument list and loops only `0..ndf-1`
  (`SRC/domain/constraints/SP_Constraint.cpp:45-93`); **extra values beyond ndf are
  silently dropped**, too few → error.
- `Domain::addSP_Constraint` rejects `dof ≥ ndf` **non-fatally** (returns false, logs)
  (`SRC/domain/domain/Domain.cpp:565-629`, check at `:582`).
- A **second, fatal** check at `SP_Constraint::setDomain` does `exit(-1)` if
  `dofNumber >= U.Size()` (`SP_Constraint.cpp:362-387`, `:376`).

### 5.2 Multi-point — `equalDOF` vs `equalDOFmixed`

- **`equalDOF`** with no explicit dof list **assumes both nodes share the builder
  ndf** (`MP_Constraint.cpp:49-101`, `ndf = OPS_GetNDF()` at `:75`) and couples
  `0..ndf-1` on both. Across dissimilar nodes this lands on the fatal setDomain bounds
  check (`MP_Constraint.cpp:281-320`, `exit(-1)` at `:300/:308`).
- **`equalDOFmixed`** is the *intended* mixed-ndf bridge: explicit
  `(retainedDof, constrainedDof)` index pairs, validated against each node's real ndf
  at setDomain (`MP_Constraint.cpp:103-198`). This is the sanctioned way to tie a
  3-dof node's translations to a 6-dof node's translations.

```tcl
# Tie translations only: constrained node 2 (ndf=3 solid) to retained node 1 (ndf=6 beam)
equalDOF_Mixed 1 2 3   1 1   2 2   3 3
# Shell/beam rotations (dofs 4,5,6 on node 1) are intentionally left uncoupled.
```

### 5.3 Rigid ties **require ndf parity** — they cannot bridge formulations

- `RigidBeam` / `RigidRod`: **require identical ndf** on both nodes; on mismatch they
  **silently abort** (return, log) (`RigidBeam.cpp:104-117`, `RigidRod.cpp:71-78`).
- `RigidDiaphragm`: hardcoded to `(6 dof, 3D)` or `(3 dof, 2D)` patterns; nodes that
  don't match are **silently ignored** with a warning
  (`RigidDiaphragm.cpp:94-99, 218-220`).

> [!warning] Across dissimilar-ndf nodes, your only sanctioned tools are
> `equalDOF_Mixed`, `zeroLength -dir`, or an ndf-aware embedding element
> (`LadrunoEmbeddedNode`). **Every `rigid*` helper assumes ndf parity** and will
> quietly skip or abort otherwise.

### 5.4 Constraint handlers

`Plain`, `Transformation`, `Lagrange`, `Penalty` handlers all size their DOF_Groups
off each node's own ndf and rely on the constraint `setDomain` bounds checks above;
none enforce model-wide ndf uniformity (`PlainHandler.cpp:101-141`;
`TransformationConstraintHandler.cpp:180-249`). Out-of-range dof indices that slip
past the earlier checks become fatal in `LagrangeMP_FE`/`PenaltyMP_FE::getResidual`.

---

## 6. Loads

- **Nodal `load`**: vector size is inferred from the argument count with *no*
  pre-check against the node (`SRC/interpreter/OpenSeesPatternCommands.cpp:140-210`),
  but `Node::addUnbalancedLoad` **validates `add.Size() == numberDOF`** and rejects on
  mismatch (returns −1, logs; no padding, no truncation)
  (`SRC/domain/node/Node.cpp:937-962`). The load then flows to the residual via
  `DOF_Group::addPtoUnbalance` (`DOF_Group.cpp:401-413`).
- **Element / beam loads** (`eleLoad`, `Beam*UniformLoad`, `SelfWeight`) live in
  *element-local* space; the **element** distributes them to node dofs (self-weight →
  translational dofs only). ndf correctness is the element's responsibility, not the
  load's (`ElementalLoad.cpp:81-92`).
- **`UniformExcitation`** writes the influence vector to dof index `theDof` with **no
  bounds check against the node's ndf** (`UniformExcitation.cpp:303-365`). In a mixed
  model this is a real footgun — exciting `dof 2` is fine on ndf≥3 nodes but undefined
  on ndf=2 nodes. (Logged in [[LEDGER_quirks]].)

---

## 7. Mass — and the explicit zero-mass trap

- Nodal `mass` must supply **exactly ndf values** → stored as a diagonal `ndf×ndf`
  matrix; `Node::setMass` rejects any non-`ndf×ndf` matrix (`Node.cpp:1269-1291`;
  parsing `Node.cpp:130-144`, `OpenSeesMiscCommands.cpp:444-483`). For ndf=6, entries
  3–5 are **rotational inertias** (`Ixx, Iyy, Izz`).
- Element mass is `numDOF×numDOF`, mapped exactly like stiffness via
  `FE_Element::addMtoTang` / `DOF_Group::addMtoTang`
  (`FE_Element.cpp:389-403`, `DOF_Group.cpp:351-363`). **Lumped mass puts ZERO on
  rotational dofs** — beam lumped leaves dofs 3–5/9–11 empty
  (`ElasticBeam3d.cpp:785-834`), and `ASDShellQ4` *explicitly* comments
  *"Rotational mass neglected for the moment"* (`ASDShellQ4.cpp:1152`). Consistent
  mass *does* populate rotational terms (beam via `Jx/A`).

> [!danger] Zero-mass-dof trap in EXPLICIT dynamics
> An ndf=6 node with translational-only mass (or a lumped shell/beam) leaves
> **zero mass on rotational dofs** → singular/ill-conditioned `M`. `CriticalTimeStep`
> **silently filters** those modes (DGGEV β-threshold, `CriticalTimeStep.cpp:166-179`)
> rather than warning, so `dt_cr` comes back **unconservatively large** and the run
> can go unstable **with no diagnostic.** In a mixed-ndf explicit model the binding
> constraint is *mass, not stiffness*: every active dof — including rotations on ndf=6
> nodes — needs nonzero mass. Relevant to [[central_difference_ladruno_guide]].
> (Logged in [[LEDGER_quirks]].)

---

## 8. Staged construction across formulations — the frame-on-soil problem

This is the canonical mixed-ndf workflow and it ties §1.4 (immutability), §3
(element contracts), §4 (rotational gap) and §5 (constraints) together.

**Scenario.** A 2D plane-strain **soil domain meshed with CST** (`Tri31`, ndf=2:
ux, uy), possibly with an excavation. **Stage 1:** geostatic / settle the soil.
**Stage 2:** introduce a **frame** (beam-columns, which in 2D need ndf=3: ux, uy,
θz).

### 8.1 Why you cannot just attach the frame to the soil nodes

- Soil CST nodes are **ndf=2** — a continuum has no rotational dof.
- 2D frame elements need **ndf=3** and `exit(-1)` / fail otherwise (§3.2).
- ndf is **fixed at creation** (§1.4) — you cannot promote an existing ndf=2 soil
  node to ndf=3 in stage 2.
- The "build the shared nodes ndf=3 from the start" escape **doesn't work**: the
  CST's ndf=3 branch is a trap (§3.4), the quad refuses it, and even if it didn't,
  the rotational dof would be **zero-stiffness/singular** during the soil-only
  stage 1.

So introducing the frame **requires a kinematic compatibility** — a constraint or
an interface/embedding element. That is not a software workaround: it is the
**physically correct** statement. The CST continuum only *has* translational dofs,
so the mechanically meaningful coupling between frame and soil is **translational**;
the frame's rotational dof is carried by the frame itself (its bending stiffness /
foundation detail). The tie *is* the physics (§4).

### 8.2 The correct construction

Give the frame its **own ndf=3 nodes**, coincident with (or embedded in) the soil,
and couple **translations only**, added in stage 2:

| Situation | Tool | Notes |
|---|---|---|
| Frame node coincident with a single soil node | `equalDOF_Mixed 1 2 3  1 1  2 2  3 3`* (couple dofs 1↔1, 2↔2) | hard constraint; frame's dof 3 (θz) left free |
| Frame interpolated into a soil **element face** (non-matching mesh) | embedding element (`ASDEmbeddedNodeElement` / [[23_ladruno_embedded_node_adr|LadrunoEmbeddedNode]]) | penalty, explicit-friendly, interpolates host shape fns |
| Frame *inside* the soil (piles) | embedding element | same |

\*example couples translations; never tie the frame's rotation to the soil — there
is nothing there to tie it to.

### 8.3 Does adding the tie **yank** the model, or is it born at the deformed shape?

The decisive property is whether the tie **captures the relative offset at
add-time** (ties *increments* → no jolt) or enforces an **absolute** displacement
equality (ties total disp → jolts the new node to chase the host's accumulated
displacement). The three tools differ — and one fork element currently gets it
wrong:

| Tie | Gap definition | Added in stage 2 onto a deformed host |
|---|---|---|
| **`equalDOF_Mixed`** (MP_Constraint) | increment tie; **`Uc0`/`Ur0` captured at `setDomain`** (`MP_Constraint.cpp:300-309`) | **No jolt** — born stress-free at deformed shape ✅ |
| **`ASDEmbeddedNodeElement`** (upstream) | `(u − U0)`; **`m_U0` captured at `setDomain`** (`ASDEmbeddedNodeElement.cpp:468-472, 795-796`), serialized in send/recvSelf | **No jolt** ✅ |
| **`LadrunoEmbeddedNode` v1** | `g = u_slave − Σ N·u_host`; **no offset capture** (`LadrunoEmbeddedKernel.cpp computeGap`) | **Jolts** — drives `u_slave → N·u_host`, force spike ❌ |

> [!danger] LadrunoEmbeddedNode v1 dropped the parent's offset capture
> The gap is built from raw trial displacements with **no `g0`** (verified:
> `LadrunoEmbeddedKernel.cpp computeGap` = `u_slave(k) − Σ N_i u_host,i(k)`, and the
> element has no `U0`/`g0` member). So a frame node created in stage 2 (`u_slave=0`)
> tied onto already-settled soil (`u_host≠0`) sees an activation gap `−N·u_host` and
> the penalty **yanks it by the full settlement**, with a spurious force spike. The
> fix is the parent `ASDEmbeddedNodeElement`'s `m_U0` pattern: capture `g0` at
> `setDomain`, use `(g − g0)`, serialize it. Until then: for staged frame-on-soil
> use **`equalDOF_Mixed`** (coincident node) or **`ASDEmbeddedNodeElement`**
> (non-matching mesh), or keep `LadrunoEmbeddedNode` **present from stage 1** so the
> slave tracks the host from zero and the gap never opens. (Logged in
> [[LEDGER_quirks]].)

> [!note] Caveat on `ASDEmbeddedNodeElement`
> It is a **small-deformation** penalty — the host interpolation/lever uses
> undeformed `getCrds()`, not the current position. The `U0` capture removes the
> initial-state jolt cleanly, but the embedding map itself does not follow large
> host rotations. For settling-soil + frame SSI (small strains after geostatic)
> that is the right regime.

### 8.4 General staging caveats (see [[constraints_reference_position]])

- **`fix`/`sp` are absolute (reference-frame):** adding one mid-stage *moves* the
  node to satisfy it. Releasing a fixed rotational dof mid-stage is delicate — one
  more reason to avoid the "fix-then-release the soil rotation" route of §8.1.
- **`equalDOF`/embedding ties are incremental (deformed-frame)** *when they capture
  the offset* — born stress-free at the current shape (the whole point of §8.3).
- **Line elements (truss/beam) default to stress-free birth** at the deformed
  length when added mid-stage.

---

## 9. Practical envelope for mixed-ndf modeling

**Supported / safe**

- Heterogeneous ndf across the domain — via per-node `node ... -ndf K` or a re-issued
  `model -ndf`. Numbering and solve are fully ndf-agnostic.
- Each element family living on its own consistent-ndf node set (all the solid nodes
  ndf=3, all the shell nodes ndf=6, etc.).
- Crossing formulations via **`equalDOF_Mixed`**, **`zeroLength -dir`**,
  **`ASDEmbeddedNodeElement`**, or **`LadrunoEmbeddedNode`** (`-rot` / `-pressure`
  only when the dofs exist). For **staged** addition onto an already-deformed host,
  prefer `equalDOF_Mixed` / `ASDEmbeddedNodeElement` (offset-captured, no jolt) —
  `LadrunoEmbeddedNode` v1 jolts unless present from stage 1 (§8.3).

**Brittle / unsupported**

- A single element spanning nodes of *different* ndf → crash (beam) at best, **silent
  corruption** (brick) at worst.
- `rigidLink` / `rigidDiaphragm` / `RigidBeam` / `RigidRod` across dissimilar ndf →
  silently skipped or aborted.
- Tying a 6-dof shell/beam to a 3-dof solid and expecting **moment transfer** →
  physically incomplete; use an embedding/transition strategy.
- `UniformExcitation` dof index and explicit zero-mass rotational dofs are
  **unchecked** — your responsibility.

**Checklist when debugging a mixed model**

1. For each node: what is its `numberDOF`? (Don't think "the model's ndf".)
2. For each element: does any element span nodes whose ndf differ, or differ from what
   that element assumes? (That's where it breaks.)
3. For each tie: is it `equalDOF_Mixed`/`zeroLength`/embedding (ok) or a `rigid*`
   (needs parity)?
4. Does every loaded dof exist on its node? Does every *active* dof have mass
   (especially under explicit)?
5. Are shell/beam rotational dofs at a solid junction actually carrying stiffness, or
   spinning free?

6. **Staged cross-formulation tie?** Is it offset-captured (no jolt:
   `equalDOF_Mixed`, `ASDEmbeddedNodeElement`) or absolute (jolts:
   `LadrunoEmbeddedNode` v1 added post-deformation)? (§8.3)

---

## 10. Source map (quick index)

| Concern | Primary file:line |
|---|---|
| `model -ndf` parse + default | `SRC/interpreter/OpenSeesCommands.cpp:1350-1411` |
| Python `model` | `SRC/interpreter/PythonModelBuilder.cpp:124-174` |
| `node ... -ndf` override (Tcl) | `SRC/modelbuilder/tcl/TclModelBuilder.cpp:1081-1099` |
| `node ... -ndf` override (OPS) | `SRC/domain/node/Node.cpp:168-189` |
| `Node::numberDOF` storage | `SRC/domain/node/Node.h:187`; ctors `Node.cpp:272-370` |
| ndf immutable (only reassign = recvSelf) | `SRC/domain/node/Node.cpp:1568`; `removeNode` `Domain.cpp:1147` |
| State-vector sizing | `SRC/domain/node/Node.cpp:1815-1839` |
| DOF_Group sized per node | `SRC/analysis/dof_grp/DOF_Group.cpp:60-64` |
| Sequential numbering | `SRC/analysis/numberer/PlainNumberer.cpp:98-115` |
| Element local→global map | `SRC/analysis/fe_ele/FE_Element.cpp:252-284` |
| Beam ndf check (exit) | `SRC/element/elasticBeamColumn/ElasticBeam3d.cpp:567-580` |
| Shell ndf warn | `SRC/element/shell/ShellMITC4.cpp:359-375` |
| Quad ndf bail | `SRC/element/fourNodeQuad/FourNodeQuad.cpp:474-485` |
| Brick (no check) | `SRC/element/brick/Brick.cpp:251-271, 319-322` |
| CST ndf=3 trap (6×6 K vs `getNumDOF`=9) | `SRC/element/triangle/Tri31.cpp:51, 450-460, 491-501, 620-635` |
| ZeroLength `-dir` | `SRC/element/zeroLength/ZeroLength.cpp:142-164` |
| `fix` / SP two-stage check | `SRC/domain/constraints/SP_Constraint.cpp:45-93, 362-387`; `Domain.cpp:582` |
| `equalDOF` / `equalDOF_Mixed` | `SRC/domain/constraints/MP_Constraint.cpp:49-101, 103-198` |
| Rigid ndf-parity requirement | `RigidBeam.cpp:104-117`, `RigidRod.cpp:71-78`, `RigidDiaphragm.cpp:94-99` |
| Nodal load size check | `SRC/domain/node/Node.cpp:937-962` |
| `UniformExcitation` (no dof check) | `SRC/domain/pattern/UniformExcitation.cpp:303-365` |
| Nodal mass `ndf×ndf` check | `SRC/domain/node/Node.cpp:1269-1291` |
| Lumped mass zero rotational | `ElasticBeam3d.cpp:785-834`, `ASDShellQ4.cpp:1152` |
| Explicit zero-mass filtering | `SRC/analysis/integrator/CriticalTimeStep.cpp:166-179` |
| Staged tie offset-capture (MP) | `SRC/domain/constraints/MP_Constraint.cpp:300-309` |
| ASDEmbeddedNode `m_U0` capture (no jolt) | `SRC/element/CEqElement/ASDEmbeddedNodeElement.cpp:468-472, 795-796` |
| Mixed-ndf embedding (fork) | `SRC/element/ladrunoEmbeddedNode/LadrunoEmbeddedNode.cpp:139-210` |
| LadrunoEmbeddedNode gap (no `g0` → jolts) | `SRC/element/ladrunoEmbeddedRebar/LadrunoEmbeddedKernel.cpp` (`computeGap`) |

---

*Verification note: §1–§7 are source-read against this fork's `SRC/` tree (multi-agent
sweep, 2026-06-07). The runtime claims most worth a guard rail — `UniformExcitation`
dof bounds and the explicit zero-mass-rotational-dof trap — are flagged in
[[LEDGER_quirks]]. Behaviors here are upstream OpenSees unless attributed to a
`Ladruno*` class.*
