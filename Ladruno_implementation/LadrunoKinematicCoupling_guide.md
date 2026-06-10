---
title: "LadrunoKinematicCoupling — RBE2 / Kinematic-Coupling Element"
project: Ladruno
type: reference / user guide
status: v1 shipped (rigid driver, U + rotation, -dof component select, penalty/AL, derived K_r, bipenalty default-off, g0 birth, 2D+3D); built + 16/16 Zone-A + 4-lens code review; MERGED to ladruno (PR #221, 2026-06-09)
element: element LadrunoKinematicCoupling
classTag: 33012 (element)
related:
  - "[[29_ladruno_kinematic_coupling_rbe2_adr]]"
  - "[[24_ladruno_coupling_constraints_adr]]"
  - "[[LadrunoDistributingCoupling_guide]]"
  - "[[LadrunoEmbeddedNode_guide]]"
  - "[[ndf_and_mixed_models_guide]]"
  - "[[LEDGER_implementations]]"
  - "[[LEDGER_vanilla_files]]"
tags:
  - guide
  - element
  - rbe2
  - kinematic-coupling
  - rigid-body
  - rigid-link
  - ndf-mismatch
  - moment-transfer
  - penalty-method
  - augmented-lagrangian
  - explicit-dynamics
updated: 2026-06-10
---

# LadrunoKinematicCoupling — RBE2 / kinematic coupling

`element LadrunoKinematicCoupling` (classTag **33012**) makes a **reference node R**
**rigidly drive** a set of **slave nodes** `{S_i}`: every tied slave DOF follows R's
rigid-body motion,

```
u_i = u_R + θ_R × d_i          (slave translations follow R, with the moment arm d_i)
θ_i = θ_R                       (slave rotations follow R, where the slave carries them and the DOF is tied)
```

with `d_i = x_i − x_R`. A load or prescribed displacement at R is transmitted to the set
**with rigid stiffness** — the region moves as one rigid body. It is the OpenSees
realization of **Nastran RBE2** / **Abaqus `*COUPLING, kinematic`** / **LS-DYNA
`*CONSTRAINED_NODAL_RIGID_BODY`**, and it generalizes the classic `rigidLink`/`rigidDiaphragm`
to an arbitrary node set, an offset reference, selectable tied components, finite penalty
(or augmented-Lagrangian) enforcement, and explicit-dynamics safety.

This is the single user/developer reference. It covers the **theory** (the rigid tie, the
gaps, the transport sign, the force transfer), the **capabilities** (the option grammar),
the **explicit-safety** machinery, the **implementation** (where the math lives,
serialization, responses), the **validation** battery, and **use cases**. The design record
is [[29_ladruno_kinematic_coupling_rbe2_adr|ADR 29]]; the family scoping is
[[24_ladruno_coupling_constraints_adr|ADR 24]].

> [!important] RBE2 vs RBE3 — the single most common wiring mistake
> In **RBE2 (this element)** the *single* node R is the **master**: it rigidly drives the
> *many* slaves and **adds stiffness** to them (the region becomes rigid). In **RBE3
> ([[LadrunoDistributingCoupling_guide|distributing coupling]])** the *single* node is the
> **dependent**, interpolated from the *many* independents, and the set is **not** stiffened.
> The "one vs many" roles **flip**. Reach for **RBE2** when the region must genuinely move as
> a **rigid body** (a loading platen, a rigid offset, a rigid connection block); reach for
> **RBE3** to **introduce/transmit a load or BC at a point while the region stays flexible**.
> Using RBE2 where RBE3 is wanted **over-stiffens** the patch (artificial local rigidity);
> using RBE3 where RBE2 is wanted leaves the region too soft.

---

## 1. Theory — the rigid coupling

### 1.1 The kinematic constraint (rigid-body driver)

R is the **master**. Each tied slave DOF is constrained to R's rigid-body field. With
`d_i = x_i − x_R` the moment arm from R to slave `i`, the per-slave constraint gaps are

```
g_i^t = u_i − ( u_R + θ_R × d_i )        (translation rows — the transport block)
g_i^r = θ_i − θ_R                         (rotation rows — identity blocks, where tied)
```

`g_i^t = 0` says the slave translation is R's translation **plus** the rigid lever
`θ_R × d_i`; `g_i^r = 0` says the slave rotation equals R's rotation. Unlike RBE3 there is
**no weighted centroid, no position-inertia `I_c`, no eigensolve, no least-squares fit** —
the tie is a *direct* per-slave rigid link. This makes RBE2 substantially simpler than its
distributing sibling (it reuses RBE3's DOF-scatter / AL / bipenalty / `g0` scaffold but
replaces the geometry resolve and the `B` build).

> [!important] Transport sign flip vs RBE3 (don't copy the line verbatim)
> RBE2 carries `−θ_R × d_i` in `g_i^t`, so the transport block is
> `∂g_i^t/∂θ_R = +[d_i]_×`, the **negative** of the operator RBE3 uses (RBE3 carries
> `+θ_R × (x_c − x_R)`). The `B` assembly is `B += −transOp`. Copying RBE3's `buildB` line
> as-is reverses every offset moment; the validation suite asserts the **sign** of the
> transport couple (test 2), not just its magnitude.

### 1.2 Selectable tied components (`-dof`) and the ragged gap layout

`-dof c1 … cK` chooses which slave components are dependent (Nastran's RBE2 dependent-DOF
list): components `1..ndm` are translations, `ndm+1..ndm+nrot` are rotations. Default = every
DOF the slave possesses. One element therefore subsumes:

- **translation-only** tie (`-dof 1 2 3` in 3D) — slaves follow R's translation + transport, but spin freely;
- **translation + transport** of an offset reference;
- **full rigid** tie (default on 6-DOF slaves) — translations *and* rotations driven.

Because slaves can be a **mix of 3-DOF and 6-DOF** nodes, the number of tied DOFs varies
per slave: the gap vector is **ragged**. The element resolves the layout **once** at
`setDomain` into three parallel index arrays — `gapNode` (which slave), `gapDof`
(node-local component), `gapIsRot` (translation vs rotation row) — so a flat uniform stride
can never mis-index the moment reaction silently.

### 1.3 Penalty formulation, residual, and consistent tangent

Penalize the gaps with a per-row penalty `D_i = diag(K_t I, K_r I)` (translation rows get
`K_t`, rotation rows get `K_r`):

```
Π_p = ½ Σ ( K_t |g_i^t|² + K_r |g_i^r|² )         (penalty potential)
t    = D g                                          (coupling traction, per row)
P    = Bᵀ t                                         (resisting force)
K    = Σ Bᵢᵀ D_i Bᵢ = BᵀDB                          (consistent tangent)
```

`g = B u − g0` is **linear** in the DOFs and `B` is **constant** (geometry only), so
`K = BᵀDB` is the **exact** penalty Hessian — symmetric, PSD, state-independent (hence
`getInitialStiff ≡ getTangentStiff`). `B` (size `nGap × nDOF`) is built once at `setDomain`.
Its non-zero blocks for tied slave `i`:

| ∂/∂ | `g_i^t` (translation rows) | `g_i^r` (rotation rows) |
|---|---|---|
| `u_R`  | `−I`              | `0` |
| `θ_R`  | `+[d_i]_×` (transport, sign-flipped vs RBE3) | `−I` |
| `u_i`  | `+I`              | `0` |
| `θ_i`  | `0`               | `+I` |

The force transfer is the work-conjugate transpose: a force at R balances against the
slave reactions (`Σ reaction = −F`), and a **moment** at R enters even a **3-DOF slave
face** as a **self-equilibrated force couple** (`Σ f_i = 0`, `Σ d_i × f_i = −M`) — the same
ndf-mismatch moment-transfer driver as RBE3, but here the patch is held **rigid** rather
than left free to deform.

> [!note] Reference configuration only (finite-rotation boundary)
> In the reference configuration `B` is constant and `K` is exact. Under a **finite** rigid
> rotation the lever `θ_R × d_i` is the first-order (small-rotation) form; a finite-rotation
> update of `d_i` (with its geometric-stiffness term) is deferred (ADR 29), shared with
> RBE3. Fine for the small-relative-rotation regime of rigid offsets, platens, and
> connection blocks.

### 1.4 Degeneracy & refusals

RBE2 has no `I_c` to go rank-deficient, so its degeneracy handling is about **ill-posed
ties**, not rotation-axis dropping:

| Situation | Handling |
|---|---|
| **Self-tie** (a slave tag `==` the reference node) | **refused** — but *inertly*: `valid=false`, `nGap=0`, matrices allocated zero (it does **not** early-return before `allocate`, which would null-`K` crash) |
| **Duplicate slave** | refused inertly (same path) |
| **All-coincident slaves** (`d_i = 0`) with a rotation tied | the default `K_r = K_t·ℓ²` is **floored** strictly positive (a zero `ℓ²` can't leave a tied rotation unpenalized); translation-only + bipenalty reports a **finite** `dtcr` (no `2√(I_p/0) = +Inf` trap) |
| Reference node with too few DOFs for the tied rotations | refused at `setDomain` |

Health check: `eleResponse $tag tiedDOFs` returns the total tied-DOF count (`Σ` over slaves
of their tied components). `0` means the element went inert (a refusal fired); a full rigid
3D tie of `N` 6-DOF slaves should read `6N`.

### 1.5 No damping; mass via bipenalty (default OFF)

A pure coupling carries **no physical damping** (`getDamp ≡ 0`; Rayleigh factors are refused
so a `βK` can't shrink the explicit step). Unlike RBE3 the reference node is **often a real
massed node** (a platen, a footing, an equipment block), so the mass-penalty machinery
**defaults OFF** — see §5.

---

## 2. Capabilities — the option grammar

```
element LadrunoKinematicCoupling $tag $refNode $N $s1 ... $sN
        [-dof $c1 ... $cK]               # dependent components on each slave (default: all the slave has)
        [-k {$Kt | auto}] [-kAlpha $a] [-host $eleTag]   # auto needs a representative -host
        [-kr $Kr]                        # rotational-tie penalty (default DERIVED K_t·ℓ²)
        [-enforce {penalty | al}]        # default penalty
        [-bipenalty {-dtcr $dt | -wcap $beta}]   # default OFF (R is often massed)
        [-absolute]                      # opt out of g0 initial-gap (offset) capture
```

| Token | Meaning | Default |
|---|---|---|
| `$refNode` | the **master** reference node (ndf ≥ ndm+nrot: **6** in 3D, **3** in 2D) | — |
| `$N $s1..sN` | count + tags of the **slave** nodes (ndf ≥ ndm; may mix 3-/6-DOF) | — |
| `-dof $c1..cK` | dependent components per slave: `1..ndm` trans, `ndm+1..ndm+nrot` rot | all DOFs the slave has |
| `-k $Kt` / `-k auto` | translational penalty; `auto` = `kAlpha·max\|K_host(i,i)\|` (needs `-host`) | `1e12` |
| `-kAlpha $a` | multiplier for `-k auto` | `1e3` |
| `-host $eleTag` | one **representative** slave-side element, used only to scale `-k auto` / `-wcap` | none |
| `-kr $Kr` | rotational-tie penalty; omit to **derive** `K_r = K_t · ℓ²` (floored) | derived |
| `-enforce` | `penalty` or `al` (augmented Lagrangian, implicit) | `penalty` |
| `-bipenalty` | explicit critical-step control (see §5); **needs** `-dtcr`/`-wcap` budget | **off** |
| `-absolute` | keep the **absolute** tie (skip stress-free `g0` capture); `-noInitGap` alias | off (capture on) |

Build the reference node with the rotational DOFs it needs — `ndf 6` in 3D, `ndf 3` in 2D
(per-node `-ndf`, see [[ndf_and_mixed_models_guide]] §1.3); the element **refuses** a
reference node with too few DOFs. Slave nodes need at least translations (`ndf ≥ ndm`); a
slave needs `ndf ≥ ndm+nrot` for its rotations to be tieable.

---

## 3. Penalty & auto-scaling (`-k`, `-kr`)

`-k $Kt` sets the translational penalty directly (default `1e12`). `-k auto` scales it from a
representative slave-side host element's initial-stiffness diagonal,
`K_t = kAlpha · max|K_host(i,i)|` — mesh/material-independent conditioning — but it requires
**`-host`** because an RBE2 set has **no single host element**. Without `-host`, use a numeric
`-k`.

The **rotational** penalty `K_r` is, by default, **derived** from `K_t` and the geometry:

```
K_r = K_t · ℓ²,   ℓ² = (lever-weighted length scale, floored to max_i |d_i|²)
```

This is mandatory, not cosmetic: a translation gap is a length and a rotation gap is a
rotation (radians), so a single penalty would be off by `ℓ²` and wreck the conditioning (and
the explicit step). For the flat-face fixture (`a = 1`, all `|d_i|² = 2`), `ℓ² = 2` so
`K_r = 2·K_t` (validation test 14). Override with `-kr $Kr` only with a specific reason.

---

## 4. Enforcement strategies (`-enforce`)

### 4.1 `penalty` (default)
The gaps are driven toward zero by `K_t`/`K_r`; the **force/moment transfer is exact for any
penalty** (§1.3), and the residual *kinematic* gap is `O(1/K)`. For most rigid-tie use the
default `1e12` is stiff enough that the region is rigid to round-off.

### 4.2 `al` (augmented Lagrangian)
Adds per-gap-row multipliers `λ` (size `nGap`) updated by a per-step Uzawa recursion in
`commitState` (`λ += D g`), recovering near-exact rigidity at **moderate** stiffness — use
when a very high `K_t` is hurting conditioning but you still need the tie tight.
**Implicit only**: under an explicit integrator the per-sub-step Uzawa update has no
equilibrium iteration to converge against; combining `-enforce al` with `-bipenalty` is
**refused at parse** (the bipenalty flag is dropped with a warning).

---

## 5. Explicit stability — bipenalty (`-bipenalty`, default OFF)

A penalty tie to a **massless** tied DOF has an unbounded frequency → zero stable step in
explicit central difference. RBE2's twist vs RBE3: **the master R is frequently a real,
massed node** (platen / footing / equipment), so `-bipenalty` **defaults OFF** and only fires
where you ask for it.

When on, the element scans **every tied DOF of R *and* every slave** and lumps a penalty mass
only on those that are **actually massless** — a **massless slave** is the RBE2-specific
hazard that RBE3's R-centric lumping would miss (validation test 9). Per-DOF stiffness for
sizing is the **Gershgorin row-sum** of the assembled penalty tangent (`≥ λ_max` ⇒
conservative; no eigensolve — a deliberate simplification vs the ADR's `jacobi3`):

```
-dtcr $dt :  m_p(DOF) = k_dof · (dt/2)²   on each massless tied DOF
             ⇒ self-reported critical step = min over lumped DOFs of 2√(m_p/k_dof) = dt
-wcap $β  :  m_p(DOF) = k_dof / (β·ω_host)²                 (needs -host for ω_host)
```

A `k_dof ≤ 0` guard skips untied/zero-stiffness DOFs, which is what prevents the
all-coincident `2√(I_p/0) = +Inf` trap. Penalty mass is lumped **diagonally** (`getMass`
stays diagonal → `DiagonalSOE`-safe). If **every** tied DOF already carries mass, the element
lumps **nothing** (R's own mass is not double-counted) and the self-report returns `0` —
"no opinion" (validation test 10). The bound is exposed via
`Element::getExplicitCriticalTimeStep` (folded into `ops.criticalTimeStep` / `-cflAbort`), and
Rayleigh factors are refused so a spurious `βK` can't shrink the step.

> [!note] `-dtcr` is a user-asserted target
> For a host-less node set there is no single host to derive an exact reduced-mass bound from;
> `-dtcr` sets the step you want the coupling to permit. Prefer `-dtcr` over `-wcap` for RBE2.

---

## 6. Deferred / not-yet-built

| Feature | Status |
|---|---|
| Finite-rotation geometric stiffness (large relative rotation) | deferred (§1.3) — shared with RBE3 |
| LS-DYNA CNRB mass-condensation integrator (true rigid-body inertia on R) | deferred (ADR 24 D3b) |
| General N-node linear-equation primitive `Σ cₖuₖ = 0` | separate element (ADR 24 D4) |
| RBE3 / distributing coupling (the flexible sibling) | shipped — [[LadrunoDistributingCoupling_guide]] |

---

## 7. Diagnostics & responses

`eleResponse $tag <name>` / `recorder Element -ele $tag -<name>`:

| Response | Aliases | Size | Meaning |
|---|---|---|---|
| `force` | `couplingForce` | nGap | the coupling traction `t = D g` (incl. `λ` under AL) |
| `gap` | — | nGap | the constraint gap `g` (relative to `g0` if captured) |
| `kt` | `k` | 1 | resolved translational penalty `K_t` |
| `kr` | — | 1 | resolved (derived or user) rotational penalty `K_r` |
| `lambda` | `augLambda` | nGap | AL multiplier vector |
| `dtcr` | `dtCritical` | 1 | self-reported explicit critical step (`0` = no opinion; `−1` = bipenalty off) |
| `tiedDOFs` | `nGap` | 1 | total tied DOF count (`0` ⇒ the element went inert — a refusal fired) |

`tiedDOFs` is the quick health check: it should equal `Σ` over slaves of their tied
components (e.g. `6N` for a full rigid 3D tie of `N` 6-DOF slaves, `3N` for translation-only).

---

## 8. Implementation map

| Concern | Where |
|---|---|
| Element | `SRC/element/ladrunoKinematicCoupling/LadrunoKinematicCoupling.{h,cpp}` |
| Parser | `…/OPS_LadrunoKinematicCoupling.cpp` |
| Geometry + ragged layout | `resolveGeometry()` (`d_i = x_i − x_R`, `gapNode`/`gapDof`/`gapIsRot`, floored `ℓ²`, self-tie/duplicate refusals) — resolved **once** at `setDomain` from `getCrds()` |
| Constant gap operator | `buildB()` → `Matrix* B` (`nGap × nDOF`); `g = B·u − g0`; transport block `+[d_i]_×` (**sign-flipped** vs RBE3) |
| Residual / tangent | `getResistingForce` (`Bᵀt`), `getTangentStiff` (`BᵀDB`); `getInitialStiff ≡` tangent; per-row penalty via `rowPenalty(row) = gapIsRot(row) ? Kr : Kt` |
| Auto / derived penalties | `resolveAutoKt()` (`-k auto` off `-host`, derive `K_r = K_t·ℓ²`) |
| Damping bypass | `getDamp` / `getRayleighDampingForces` return element-owned **zeroed** `C0`/`dampF` (a no-op `setRayleighDampingFactors` alone would crash transient — see [[LEDGER_quirks]]) |
| Explicit | `getMass` (`M0` diagonal, massless-scan lumps), `resolveBipenalty()` (Gershgorin row-sum over R **and** slaves), `getResistingForceIncInertia`, `getExplicitCriticalTimeStep` |
| Serialization | `sendSelf`/`recvSelf` carry `dofSel`/weights flags/`λ`/`g0` + options; geometry, ragged layout & `B` **recomputed** on recv from coords |
| Plumbing | `classTags.h` (`ELE_TAG_LadrunoKinematicCoupling = 33012`), `FEM_ObjectBrokerAllClasses.cpp`, `OpenSeesElementCommands.cpp`, CMake — see [[LEDGER_vanilla_files]] |

### 8.1 Build-bug lessons folded in (don't relearn)
- **`getDamp` override is mandatory**, not optional: the no-op `setRayleighDampingFactors`
  means the base `Element` never allocates its lazy damping slot, so the base `getDamp` would
  deref `theMatrices[−1]` the first time a **transient** integrator forms the C-tangent → hard
  crash in the element's *primary* (dynamic) use. The 16 quasi-static-leaning tests can't catch
  it without the transient-Newmark smoke (test 11). Any element with a no-op
  `setRayleighDampingFactors` **must** override `getDamp()` + `getRayleighDampingForces()`.
- **Refuse inertly, after `allocate`** — self-tie / duplicate-slave set `valid=false` but still
  build zeroed matrices; an early-return-before-`allocate` would leave a null `K` → crash.
- **`getMass` double-zero guard** — `resolveBipenalty` is `bpResolved`-guarded so a second call
  doesn't wipe `M0`.
- **`lambdaAL` / `g0` first-call detection before resize** — resizing destroys the size-0
  fresh-vs-recv signal; check the size *before* resizing.
- **Test gotcha:** `ops.mass(node, …)` sizes the mass vector by the **model** `ndf`, not the
  node's — to mass a 6-DOF reference node, use `model -ndf 6` (massing a `-ndf 6` node in a
  `-ndf 3` model → "Node::setMass incompatible matrices").

---

## 9. Validation

Zone-A battery `tests/test_ladrunoKinematicCoupling_element.py` — **16/16**, plus full
Zone-A 633-pass no-regression and a 4-lens adversarial **code** review (2 CRITICAL + 6 MAJOR
folded in; ADR 29). The kinematic tests exploit a clean fact: with the slaves otherwise free,
their only stiffness is the penalty tie, so equilibrium drives each slave **exactly** onto R's
rigid prediction (gap → 0, independent of `K`).

- **Kinematics:** rigid translation (slaves follow, zero gap); **rotation + transport with the
  asserted SIGN** on an offset reference (the transOp sign-flip guard); full-rigid drives slave
  rotation (`θ_i = θ_R`, `tiedDOFs = 24`); translation-only `-dof 1 2 3` leaves slave rotation
  free (`tiedDOFs = 12`); single rigid link `N = 1` (generalized `rigidLink` with a moment arm).
- **Statics:** force balance (`Σ reaction = −F`); **moment into a 3-DOF face** as a
  self-equilibrated couple (`Σ f = 0`, `Σ d_i × f_i = −M`).
- **Robustness:** AL solves at finite `K`; derived `K_r = K_t·ℓ²` (`= 2·K_t` on the flat face);
  all-coincident → finite `dtcr` + floored `K_r > 0`; self-tie refused inertly.
- **Explicit / dynamics:** bipenalty massless-reference self-report (`dtcr == dt`);
  **massless-slave scanned** (the RBE2-specific hazard); **massed-R not double-counted**
  (`dtcr == 0`); **transient-Newmark smoke** (the `getDamp` regression guard).
- **Serialization:** FE_Datastore `sendSelf`/`recvSelf` round-trip (probes geometry-derived
  `kr` + `tiedDOFs` → fails if recv didn't reconstruct the element from coords).

---

## 10. Use cases & recipes

### 10.1 Loading platen / rigid bearing on a deformable face (the driver)
A rigid platen pressing on a solid specimen face: drive the face nodes rigidly from one
control node, prescribe its motion or load it, and the patch translates+rotates as one rigid
body (contrast RBE3, which would let the face deform under the load):

```python
ops.model('basic', '-ndm', 3, '-ndf', 3)
# ... solid face nodes 101..104 (ndf 3) ...
ops.node(1, xc, yc, zc, '-ndf', 6)                         # platen control node (master, 6-DOF)
ops.element('LadrunoKinematicCoupling', 1, 1, 4, 101, 102, 103, 104,
            '-dof', 1, 2, 3,                               # tie translations (3-DOF slaves)
            '-k', 'auto', '-host', 5001)                   # scale penalty off solid element 5001
ops.sp(1, 3, -0.01)                                        # press the platen down 10 mm
```

### 10.2 Rigid offset / rigid link (generalized rigidLink)
Tie a node to a master through a rigid arm — a beam end to an offset working point, a
sensor node rigidly carried by a member, a rigid connection block. `N = 1` reduces exactly to
`rigidLink` but with full moment-arm transport and an offset reference:

```python
ops.element('LadrunoKinematicCoupling', 7, masterNode, 1, slaveNode, '-k', 1.0e10)
```

### 10.3 Rigid body / rigid diaphragm over an arbitrary node set
Make a group of nodes move as one rigid body driven by a single reference (a rigid footing, a
rigid cap, an equipment skid, an arbitrary-shape diaphragm that `rigidDiaphragm` can't
express). Use `-dof` to tie only the in-plane DOFs for a diaphragm, or all 6 for a full rigid
body.

### 10.4 Beam/shell-to-solid rigid moment transfer (ndf-mismatch, rigid variant)
The rigid counterpart of RBE3's signature transfer: a 6-DOF beam/shell node framing into a
3-DOF solid face, delivering the member's moment into the continuum as a force couple **while
holding the patch rigid**. Choose RBE2 here when the joint block is genuinely stiff (a thick
gusset, an embedded plate); choose [[LadrunoDistributingCoupling_guide|RBE3]] when the face
must stay flexible.

```python
ops.node(1, xb, yb, zb, '-ndf', 6)                         # beam node (6-DOF master)
ops.element('LadrunoKinematicCoupling', 1, 1, 4, 101, 102, 103, 104, '-dof', 1, 2, 3)
ops.load(1, 0, 0, -P, 0, Mx, 0)                            # force + moment at the master
```

### 10.5 Explicit dynamics (impact / SSI / blast)
R is **often massed** here, so bipenalty stays off unless a tied DOF is massless. If R or a
slave is massless, add `-bipenalty -dtcr <dt>` at/below your explicit step and query the bound:

```python
ops.mass(1, m, m, m, J, J, J)                              # massed master -> no penalty mass lumped
ops.element('LadrunoKinematicCoupling', 1, 1, 4, 101, 102, 103, 104,
            '-k', 1.0e8, '-bipenalty', '-dtcr', 2.0e-6)    # fires only on massless tied DOFs
dt_cr = ops.eleResponse(1, 'dtcr')[0]                      # honored by ops.criticalTimeStep / -cflAbort
```

### 10.6 Staged construction (stress-free birth)
`g0` capture is **on by default**: a coupling added to an already-deformed model is born
force-free (the gap is measured relative to the activation state, including any pre-existing
offset between R's predicted position and the slaves). Add `-absolute` only for the legacy
absolute snap-to-master tie.

### 10.7 apeGmsh integration — generator contract (recommended)
For a rigid platen/connection-on-solid interface, apeGmsh should emit one
`LadrunoKinematicCoupling` per master: resolve the host face the master drives, pass those
nodes as the slave set, choose `-dof` by the slave ndf (translations for a 3-DOF face), name a
representative face element via `-host` for `-k auto`, and add `-bipenalty -dtcr` only when a
tied node is massless. A `g.couple(ref, set, mode='kinematic', dof=...)` wrapper is the natural
API surface (cf. ADR 24 §6 and the RBE3 `mode='distributing'` sibling).

### 10.8 Pitfalls
- **Expecting RBE2 to leave the patch flexible** → it doesn't; it makes the set **rigid**.
  Use [[LadrunoDistributingCoupling_guide|RBE3]] for a load-introduction-without-stiffening.
- **`tiedDOFs == 0`** → the element went inert: a **self-tie** (slave tag equals the reference)
  or a **duplicate slave**. Check the slave list.
- **Reference node built with too few DOFs** (e.g. `ndf 3` in 3D) → refused at `setDomain`.
  Build it `ndf 6` (3D) / `ndf 3` (2D).
- **`-k auto` without `-host`** → refused (no single host for a node set); use a numeric `-k`.
- **Single penalty across translation and rotation** → don't override `-kr` with `K_t`; the
  derived `K_t·ℓ²` is there for a reason (§3).
- **Bipenalty default surprise** → unlike RBE3, RBE2 bipenalty is **off** by default (R is
  usually massed). Turn it on explicitly if any tied node is massless in an explicit run.

---

## 11. References & related
- **[[29_ladruno_kinematic_coupling_rbe2_adr|ADR 29]]** — the design record (formulation, the
  transport sign flip, the ragged layout, explicit safety, the 4-lens review must-fixes).
- **[[24_ladruno_coupling_constraints_adr|ADR 24]]** — the coupling-constraint family (RBE2 /
  RBE3 / linear equation) and the Abaqus/LS-DYNA/Nastran alignment.
- **[[LadrunoDistributingCoupling_guide]]** — the **flexible sibling (RBE3)**: R is the
  dependent weighted-average of the set, adds no stiffness. Reach for it when RBE2 would
  over-stiffen.
- **[[LadrunoEmbeddedNode_guide]]** — the node-into-host **embedment** tie (shares the
  penalty/AL/bipenalty kernel).
- **[[ndf_and_mixed_models_guide]]** — mixed-ndf models; §4 the rotational-compatibility gap;
  §7 the explicit zero-mass-rotational trap.
- **Nastran** RBE2 · **Abaqus** `*COUPLING, kinematic` · **LS-DYNA**
  `*CONSTRAINED_NODAL_RIGID_BODY`.
</content>
</invoke>
