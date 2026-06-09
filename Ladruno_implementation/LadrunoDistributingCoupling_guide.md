---
title: "LadrunoDistributingCoupling — RBE3 / Distributing-Coupling Element"
project: Ladruno
type: reference / user guide
status: v1 shipped (U + rotation, penalty/AL, derived K_r, bipenalty, g0 birth, 2D+3D); built + 15/15 Zone-A + 4-lens code review
element: element LadrunoDistributingCoupling
classTag: 33011 (element)
related:
  - "[[28_ladruno_distributing_coupling_rbe3_adr]]"
  - "[[24_ladruno_coupling_constraints_adr]]"
  - "[[23_ladruno_embedded_node_adr]]"
  - "[[LadrunoEmbeddedNode_guide]]"
  - "[[ndf_and_mixed_models_guide]]"
  - "[[LEDGER_implementations]]"
  - "[[LEDGER_vanilla_files]]"
tags:
  - guide
  - element
  - rbe3
  - distributing-coupling
  - kinematic-coupling
  - ndf-mismatch
  - moment-transfer
  - penalty-method
  - augmented-lagrangian
  - explicit-dynamics
updated: 2026-06-09
---

# LadrunoDistributingCoupling — RBE3 / distributing coupling

`element LadrunoDistributingCoupling` (classTag **33011**) couples a **reference node R**
to a set of **independent nodes** `{I_i}` so that R's motion is the **weighted
least-squares rigid-body fit** of the set — and, dually, a load applied at R is
**distributed** to the set as a statically-equivalent force pattern, **adding no
stiffness** to the independents. It is the OpenSees realization of **Nastran RBE3** /
**Abaqus `*COUPLING, distributing`** / **LS-DYNA `*CONSTRAINED_INTERPOLATION`**.

The fork driver is the **ndf-mismatch moment-transfer interface**: a **6-DOF beam/shell
reference node** framing into a face/ring of **3-DOF solid nodes**. RBE3 delivers the
beam's **moment into the continuum as a distributed force couple** while leaving the
solid face **free to deform** and injecting **no rotational stiffness** into it — the
mechanically faithful answer where a translation-only tie lets the moment leak
([[ndf_and_mixed_models_guide]] §4) and a rigid (RBE2) tie over-stiffens the patch.

This is the single user/developer reference. It covers the **theory** (the fit, the
gaps, the force split, the degeneracy handling), the **capabilities** (the option
grammar), the **explicit-safety** machinery, the **implementation** (where the math
lives, serialization, responses), the **validation** battery, and **use cases**. The
design record is [[28_ladruno_distributing_coupling_rbe3_adr|ADR 28]]; the family scoping
is [[24_ladruno_coupling_constraints_adr|ADR 24]].

> [!important] RBE3 vs RBE2 — the single most common wiring mistake
> In **RBE3 (this element)** the *single* node is the **dependent** node, interpolated
> from the *many* independents; the set is **not** stiffened. In **RBE2 / kinematic
> coupling** the *single* node is the **master**, rigidly driving the *many* slaves
> (adds rigid stiffness). The "one vs many" roles **flip**. Reach for RBE3 to **apply or
> transmit a load/BC at a point while the region stays flexible**; reach for RBE2 when the
> region must genuinely move as a **rigid body**. RBE2 is not yet built (ADR 24 D3).

---

## 1. Theory — the coupling

### 1.1 Geometry: the weighted centroid and the position-inertia

With independent positions `x_i`, user weights `w_i > 0`, `W = Σ w_i`:

```
weighted centroid   x_c = (1/W) Σ w_i x_i          (the fit origin — NEVER x_R)
relative positions  r_i = x_i − x_c                ⇒ Σ w_i r_i = 0
position-inertia    I_c = Σ w_i [ (r_i·r_i) 1 − r_i ⊗ r_i ]   (3×3 in 3D; scalar Σw|r|² in 2D)
transport lever     d   = x_c − x_R
```

`I_c` is a **geometric** tensor (the weighted second-moment of the point cloud), **not**
a thin-plate inertia and **not** a penalty stiffness. The property `Σ w_i r_i = 0` (true
only because `x_c` is the *weighted* centroid) is what decouples the translation fit from
the rotation fit; measuring `r_i` from any other origin would couple them.

### 1.2 The kinematic constraint (the least-squares fit)

R's motion is the rigid-body field `(a, ω)` that best fits the independent translations,
`min Σ w_i |u_i − (a + ω × r_i)|²`. The minimiser is

```
translation   a = (1/W) Σ w_i u_i                  (weighted centroid translation)
rotation      ω = I_c⁻¹ Σ w_i (r_i × u_i)          (moment-balanced LSQ rotation)
```

The element enforces this via **constraint gaps** (kinematically complete for a
reference node that is **not** at the centroid — the `θ_R × d` transport term):

```
g_t = ( u_R + θ_R × (x_c − x_R) ) − (1/W) Σ w_i u_i          (size ndm)
g_r =   P θ_R − I_c⁺ Σ w_i (r_i × u_i)                        (size nrot: 3 in 3D, 1 in 2D)
```

`I_c⁺` is the spectral **pseudo-inverse** and `P` its range **projector** (§1.5). Dropping
the transport term (the naïve `g_t = u_R − Σ(w_i/W)u_i`) is silently wrong whenever R is
offset from the centroid and force + moment are transmitted together.

### 1.3 The work-conjugate force distribution

Penalise the gaps: `t_t = K_t g_t`, `t_r = K_r g_r`. The resisting force is `P = Bᵀ t`
where `B = ∂g/∂u` (constant — geometry only). The force delivered to independent node `i`
is then the **textbook RBE3 split** — the work-conjugate transpose of the kinematics:

```
f_i = (w_i/W) t_t            ← force part:  Σ f_i = t_t,    Σ r_i × f_i = 0
    + w_i (α × r_i),  α = I_c⁺ t_r   ← moment part: Σ r_i × f_i = t_r,  Σ f_i = 0
```

A **moment** at R becomes a **self-equilibrated swirl** of forces `w_i(α × r_i)` over the
independents — exactly how a continuum carries a moment (as distributed traction), with no
nodal couple and no rotational DOF required on the independents. **Net force and moment
balance hold for any penalty `K`** — distribution exactness does not depend on stiffness
(that is the defining RBE3 property; a finite `K` only relaxes the *kinematic* tie, which
augmented-Lagrangian restores — §4.2).

### 1.4 Penalty formulation, residual, and consistent tangent

```
Π_p = ½ K_t |g_t|² + ½ K_r |g_r|²        (penalty potential)
P   = Bᵀ t,   t = diag(K_t I, K_r I) · g          (resisting force)
K   = Bᵀ diag(K_t I, K_r I) B                     (consistent tangent)
```

Because `g = B u − g0` is **linear** in the DOFs and `B` is constant, `K = BᵀDB` is the
**exact** penalty Hessian — symmetric, PSD, and state-independent (so `getInitialStiff ≡
getTangentStiff`). `B` is built **once** at `setDomain` from the node coordinates; its rows
are:

| ∂/∂ | `g_t` (translation) | `g_r` (rotation) |
|---|---|---|
| `u_R`  | `I`                  | `0` |
| `θ_R`  | `[d]_×ᵀ` (transport) | `P` |
| `u_i`  | `−(w_i/W) I`         | `−w_i I_c⁺ [r_i]_×` |

> [!note] Reference configuration only (finite-rotation boundary)
> In the reference configuration `B` is constant and `K` is exact. Under a **finite**
> rigid rotation the gap vanishes only to **first order** (`O(φ_rel²)`, with `φ_rel` the
> rotation of R *relative to the set's frame*) — fine for small-relative-rotation
> beam/shell-on-solid transfer. A finite-rotation update of `r_i` (with its geometric
> stiffness term) is deferred (ADR 28 §4.2); it is **not** a parser flag.

### 1.5 Degeneracy — flat face is full-rank; collinear sets drop an axis

`I_c` can be rank-deficient. The element eigendecomposes it (a 3×3 Jacobi solve; scalar in
2D) at `setDomain` and forms the pseudo-inverse / projector over the **kept** subspace,
**dropping** axes below a unit-invariant floor:

```
λ_floor = max( 1e-8 · λ_max ,  1e-10 · W · L² ),   L² = max_i |r_i|²
keep eigenpair k ⇔ λ_k > λ_floor   ⇒   I_c⁺ = Σ_keep λ_k⁻¹ v_k v_kᵀ,  P = Σ_keep v_k v_kᵀ
```

| Independent set | `I_c` | reference rotation |
|---|---|---|
| **Coplanar, 2-D spread (a flat face — the main case)** | **full rank** | all resolvable ✓ |
| Collinear (an edge of nodes, or N = 2) | one zero eigenvalue (the line axis) | that axis **dropped**, left free, warned |
| N = 1 | `I_c = 0` | all rotation dropped (translation tie only) |
| Near-collinear strip | one tiny eigenvalue | dropped (spectral drop, **not** Tikhonov) |

> [!tip] The flat face is NOT degenerate
> A common worry — "a flat solid face is coplanar, so `I_c` is singular" — is **false**.
> `I_c` is a *position* inertia: for a square face its eigenvalues are `(4a², 4a², 8a²)` —
> the normal/drilling axis is the **largest**, not vanishing. All three reference rotations
> are resolvable. Regularization fires only for genuinely collinear / single-node sets.

A dropped rotation axis leaves the corresponding reference rotation DOF **unconstrained**
(physically correct — that rotation is genuinely indeterminate from the set). In practice
R carries that rotation through its own attached member (a beam), so the global system is
non-singular; a standalone model with a dropped axis and a free reference rotation will be
singular — restrain it or accept it is unconstrained. Read `eleResponse $tag nKept` (§7).

### 1.6 No damping; mass via bipenalty

A pure coupling carries **no physical damping** (`getDamp ≡ 0`; Rayleigh factors are
refused so a `βK` can't shrink the explicit step). The reference node is **massless** by
construction — the explicit-dynamics hazard handled in §5.

---

## 2. Capabilities — the option grammar

```
element LadrunoDistributingCoupling $tag $refNode $N $i1 ... $iN
        [-w $w1 ... $wN]                 # weights (default equal); -area deferred (warns)
        [-k {$Kt | auto}] [-kAlpha $a] [-host $eleTag]   # auto needs a representative -host
        [-kr $Kr]                        # rotational penalty (default DERIVED K_t·ℓ²)
        [-enforce {penalty | al}]        # default penalty
        [-bipenalty {-dtcr $dt | -wcap $beta}]
        [-absolute]                      # opt out of g0 initial-gap capture
```

| Token | Meaning | Default |
|---|---|---|
| `$refNode` | the **dependent** reference node (ndf ≥ ndm+nrot: **6** in 3D, **3** in 2D) | — |
| `$N $i1..iN` | count + tags of the **independent** nodes (ndf ≥ ndm) | — |
| `-w` | per-independent weights `w_i ≥ 0` (e.g. tributary areas) | all `1.0` |
| `-k $Kt` / `-k auto` | translational penalty; `auto` = `kAlpha·max\|K_host(i,i)\|` (needs `-host`) | `1e12` |
| `-kAlpha $a` | multiplier for `-k auto` | `1e3` |
| `-host $eleTag` | one **representative** independent-side element, used only to scale `-k auto` / `-wcap` | none |
| `-kr $Kr` | rotational penalty; omit to **derive** `K_r = K_t · (Σ w_i\|r_i\|²/W)` | derived |
| `-enforce` | `penalty` or `al` (augmented Lagrangian, implicit) | `penalty` |
| `-bipenalty` | explicit critical-step control (see §5) | off |
| `-absolute` | keep the **absolute** tie (skip stress-free `g0` capture) | off (capture on) |

Build the reference node with the rotational DOFs it needs — `ndf 6` in 3D, `ndf 3` in 2D
(per-node `-ndf`, see [[ndf_and_mixed_models_guide]] §1.3); the element **refuses** a
reference node with too few DOFs. Independent nodes need only translations (`ndf ≥ ndm`).

---

## 3. Penalty & auto-scaling (`-k`, `-kr`)

`-k $Kt` sets the translational penalty directly (default `1e12`). `-k auto` scales it from
a representative host element's initial-stiffness diagonal, `K_t = kAlpha · max|K_host(i,i)|`
— mesh/material-independent conditioning, but it requires **`-host`** because an RBE3 set
has **no single host element** (unlike the embedded-node tie). Without `-host`, use a
numeric `-k`.

The **rotational** penalty `K_r` is, by default, **derived** from `K_t` and the geometry:

```
K_r = K_t · ℓ²,   ℓ² = Σ w_i |r_i|² / W      (the lever-weighted scale)
```

This is mandatory, not cosmetic: `g_t` is a length and `g_r` is a rotation (radians), so a
single penalty would be off by a factor `ℓ²` and wreck the conditioning (and the explicit
step). The derived `K_r` makes a unit rotation gap cost energy comparable to a unit
translation gap. Override with `-kr $Kr` only if you have a specific reason.

---

## 4. Enforcement strategies (`-enforce`)

### 4.1 `penalty` (default)
The gap is driven toward zero by `K_t`/`K_r`; the **force distribution is exact for any
penalty** (§1.3), so for the dominant use case (introducing/transmitting a load) plain
penalty is mechanically sufficient. Residual kinematic gap is `O(1/K)`.

### 4.2 `al` (augmented Lagrangian)
Adds per-element multipliers `λ`, `λ_r` updated by a per-step Uzawa recursion in
`commitState` (`λ += K_t g_t`, `λ_r += K_r g_r`), recovering near-exact kinematic
satisfaction at **moderate** stiffness — use when you need R's motion to *be* the weighted
fit (not just the force distributed). **Implicit only**: under an explicit integrator the
per-sub-step Uzawa update has no equilibrium iteration to converge against, so it provides
no benefit; combining `-enforce al` with `-bipenalty` is refused at parse.

---

## 5. Explicit stability — bipenalty (`-bipenalty`)

The reference node is a **massless interpolation point**, so a penalty spring to it has an
unbounded frequency → zero stable step in explicit central difference. `-bipenalty` lumps a
**mass penalty** on R so the coupled modes self-bound, **per DOF class**:

```
-dtcr $dt :  m_p = K_t (dt/2)²   (translational),   I_p = K_r (dt/2)²   (rotational)
             ⇒ self-reported critical step = min( 2√(m_p/K_t), 2√(I_p/K_r) ) = dt
-wcap $β  :  m_p = K_t /(β·ω_host)²,  I_p = K_r/(β·ω_host)²   (needs -host for ω_host)
```

Both `m_p` and `I_p` are lumped on R only (`getMass` stays diagonal → `DiagonalSOE`-safe).
The derived `K_r = K_t·ℓ²` is what lets the rotation class self-bound to the **same** `dt`
— without a real `K_r` the massless rotation DOF would fall into the zero-mass-rotational
`CriticalTimeStep` trap ([[ndf_and_mixed_models_guide]] §7) and the run would go unstable
with no diagnostic. The element reports its bound through `Element::getExplicitCriticalTimeStep`
(folded into `ops.criticalTimeStep` / `-cflAbort`), and refuses Rayleigh factors so a
spurious `βK` can't shrink the step.

> [!note] `-dtcr` is a user-asserted target
> For a host-less node set there is no single host to derive an exact reduced-mass bound
> from; `-dtcr` sets the step you want the coupling to permit (the independents carry real
> assembled mass). Prefer `-dtcr` over `-wcap` for RBE3. Mass-redistribution to the
> independents (LS-DYNA `*CONSTRAINED_INTERPOLATION` style) is deferred (ADR 28 §6.4).

---

## 6. Deferred / not-yet-built

| Feature | Status |
|---|---|
| `-area` tributary-area auto-weights | parser **warns** and falls back to equal; supply `-w` explicitly |
| Finite-rotation geometric stiffness (large relative rotation) | deferred (§1.4) |
| Mass redistribution to independents (vs bipenalty) | deferred (ADR 28 §6.4) |
| RBE2 / kinematic coupling (rigid driver) | separate element (ADR 24 D3) |
| General N-node linear equation `Σ cₖuₖ=0` | separate (ADR 24 D4) |

---

## 7. Diagnostics & responses

`eleResponse $tag <name>` / `recorder Element -ele $tag -<name>`:

| Response | Size | Meaning |
|---|---|---|
| `force` (`couplingForce`) | ndm+nrot | the coupling traction `[t_t ; t_r]` (incl. `λ` under AL) |
| `gap` | ndm+nrot | the constraint gap `[g_t ; g_r]` (relative to `g0` if captured) |
| `kt` (`k`) | 1 | resolved translational penalty `K_t` |
| `kr` | 1 | resolved (derived or user) rotational penalty `K_r` |
| `lambda` (`augLambda`) | ndm | AL translational multiplier |
| `lambdaR` (`augLambdaR`) | nrot | AL rotational multiplier |
| `mpenalty` (`massPenalty`) | 1 | bipenalty translational mass `m_p` |
| `dtcr` (`dtCritical`) | 1 | self-reported explicit critical step |
| `nKept` (`rotationAxes`) | 1 | resolvable reference-rotation axes (`< nrot` ⇒ a degenerate set) |

`nKept` is the quick health check: it should equal `nrot` (3 in 3D, 1 in 2D) for a
well-posed set; a smaller value means the geometry dropped a rotation axis (§1.5).

---

## 8. Implementation map

| Concern | Where |
|---|---|
| Element | `SRC/element/ladrunoDistributingCoupling/LadrunoDistributingCoupling.{h,cpp}` |
| Parser | `…/OPS_LadrunoDistributingCoupling.cpp` |
| Geometry + operators | `resolveGeometry()` (x_c, r_i, d, ℓ², `jacobi3` eigensolve → `I_c⁺`/`P`, nKept) — resolved **once** at `setDomain` from `getCrds()` |
| Constant gap operator | `buildB()` → `Matrix* B` (ndm+nrot × nDOF); `g = B·u − g0` |
| Residual / tangent | `getResistingForce` (`Bᵀt`), `getTangentStiff` (`BᵀDB`); `getInitialStiff ≡` tangent |
| Translation reuse | `LadrunoEmbeddedKernel` math (`maxAbsDiagonal`, `massPenalty*`, `criticalTimeStep`) — the rotation+transport blocks are RBE3-specific |
| Damping bypass | `getDamp` / `getRayleighDampingForces` return element-owned **zeroed** `C0`/`dampF` (a no-op `setRayleighDampingFactors` alone would crash transient — see [[LEDGER_quirks]]) |
| Explicit | `getMass` (`m_p`/`I_p` on R), `getResistingForceIncInertia`, `getExplicitCriticalTimeStep` |
| Serialization | `sendSelf`/`recvSelf` carry weights/λ/λ_r/`g0`/`gr0` + flags; geometry & `B` **recomputed** on recv from coords |
| Plumbing | `classTags.h` (33011), `FEM_ObjectBrokerAllClasses.cpp`, `OpenSeesElementCommands.cpp`, CMake — see [[LEDGER_vanilla_files]] |

---

## 9. Validation

Zone-A battery `tests/test_ladrunoDistributingCoupling_element.py` — **15/15**, plus a
4-lens adversarial **code** review (math / serialization / numerics / integration; ADR 28
§8b). The battery exploits a clean fact: with the reference node otherwise free, its only
force is the penalty, so equilibrium drives the gap to **zero exactly** (independent of K)
and R lands exactly on the LSQ fit — making the kinematic assertions exact, not
penalty-approximate.

- **Kinematics:** weighted-average translation; LSQ rotation fit on a flat face; **offset-
  reference transport** (`u_R = −θ_R×d`, the term that distinguishes the complete
  formulation); rigid-translation zero-gap; 2D drilling fit; unequal-weight weighted mean.
- **Statics:** force balance (`Σ f_i = F`, equal split); **moment into a 3-DOF set**
  (`Σ f_i = 0`, `Σ r_i×f_i = −M` — the ndf-mismatch driver against a 3-DOF independent set).
- **Degeneracy:** flat-face full rank (`nKept = 3`); collinear drops one axis (`nKept = 2`,
  no NaN); derived `K_r = K_t·ℓ²`.
- **Dynamics / robustness:** bipenalty `dtcr == dt`; AL solves; **transient Newmark smoke**
  (the regression guard for the `getDamp` fix); FE_Datastore `sendSelf`/`recvSelf` round-trip
  (probes geometry-derived element state).

---

## 10. Use cases & recipes

### 10.1 Beam/shell node into a solid face — moment transfer (the driver)
A 6-DOF beam end-node `1` framing into the 4 nodes of a solid face; the beam's moment
spreads over the face, the face deforms, no rotational stiffness is injected:

```python
ops.model('basic', '-ndm', 3, '-ndf', 3)
# ... solid nodes 101..104 on the face (ndf 3) ...
ops.node(1, xb, yb, zb, '-ndf', 6)                 # beam node (6-DOF)
# ... the beam element attaches to node 1 ...
ops.element('LadrunoDistributingCoupling', 1, 1, 4, 101, 102, 103, 104,
            '-k', 'auto', '-host', 5001)           # scale penalty off solid element 5001
```

### 10.2 Control-point load / BC over a surface (no over-stiffening)
Apply a concentrated force/moment (or a prescribed displacement) at one reference point and
let it distribute over a bearing surface — a column base on a footing, an actuator head on a
specimen face, a prestress-anchor — without the artificial stiffness a rigid tie injects.
Use **tributary-area weights** so the distribution matches a uniform traction:

```python
ops.element('LadrunoDistributingCoupling', 7, 7, 9, *ring_nodes,
            '-w', *tributary_areas)
ops.load(7, 0, 0, -P, 0, 0, M)                     # force + moment at the control point
```

### 10.3 Pile head ↔ cap / reaction distribution
Spread a pile-head reaction across a cap node set, or introduce a column reaction into a
mat foundation, through a weighted average — moment-carrying, set free to deform.

### 10.4 Explicit dynamics (SSI / impact)
The reference node is massless, so add `-bipenalty -dtcr <dt>` and pick `dt` at/below your
explicit step; query the bound and abort if violated:

```python
ops.element('LadrunoDistributingCoupling', 1, 1, 4, 101, 102, 103, 104,
            '-k', 1.0e8, '-bipenalty', '-dtcr', 2.0e-6)
dt_cr = ops.eleResponse(1, 'dtcr')[0]              # honored by ops.criticalTimeStep / -cflAbort
```

### 10.5 Staged construction (stress-free birth)
`g0` capture is **on by default**: a coupling added to an already-deformed model is born
force- and stress-free (the gap is measured relative to the activation state). Add `-absolute`
only if you deliberately want the legacy absolute tie (a snap-to-set).

### 10.6 apeGmsh integration — generator contract (recommended)
For a beam/shell-on-solid interface, apeGmsh should emit one `LadrunoDistributingCoupling`
per reference node: resolve the host solid face the reference projects onto, pass that
face's nodes as the independent set with **tributary-area `-w`**, name a representative
face element via `-host` for `-k auto`, and add `-bipenalty -dtcr` in explicit jobs. The
reference node keeps its 6 DOFs; the solid face nodes stay 3-DOF. (A `g.couple(ref, set,
mode='distributing', weights=...)` wrapper is the natural API surface — cf. ADR 24 §6.)

### 10.7 Pitfalls
- **`nKept < nrot`** → the set dropped a rotation axis (collinear/near-collinear/single
  node). That reference rotation is unconstrained; restrain it or attach a member that
  carries it. Check `eleResponse $tag nKept`.
- **Reference node built with too few DOFs** (e.g. `ndf 3` in 3D) → refused at `setDomain`.
  Build it `ndf 6` (3D) / `ndf 3` (2D).
- **`-k auto` without `-host`** → refused (no single host for a node set); use a numeric `-k`.
- **Single penalty across translation and rotation** → don't override `-kr` with `K_t`; the
  derived `K_t·ℓ²` is there for a reason (§3).
- **Expecting RBE3 to make the set rigid** → it doesn't; that's RBE2 (not yet built).

---

## 11. References & related
- **[[28_ladruno_distributing_coupling_rbe3_adr|ADR 28]]** — the design record (formulation,
  degeneracy algorithm, explicit safety, 4-lens reviews).
- **[[24_ladruno_coupling_constraints_adr|ADR 24]]** — the coupling-constraint family (RBE2 /
  RBE3 / linear equation) and the Abaqus/LS-DYNA/Nastran alignment.
- **[[LadrunoEmbeddedNode_guide]]** — the sibling node-into-host **embedment** tie (shares the
  penalty/AL/bipenalty kernel; UR there is a continuum-spin tie, *not* the RBE3 LSQ fit).
- **[[ndf_and_mixed_models_guide]]** — how mixed-ndf models work; §4 the rotational-
  compatibility gap RBE3 fills; §7 the explicit zero-mass-rotational trap.
- **Nastran** RBE3 · **Abaqus** `*COUPLING, distributing` · **LS-DYNA**
  `*CONSTRAINED_INTERPOLATION`.
