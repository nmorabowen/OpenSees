---
title: Ladruno Contact — 2D plane-model contact engine
project: Ladruno
status: draft (rev 2 — post adversarial panel 2026-08-18)
priority: medium
owner: nmora
tags:
  - implementation
  - contact
  - 2d
  - plane
  - design-adr
---

# Ladruno Contact — 2D plane-model contact engine (ADR-85)

> **Design ADR, not code.** This scopes bringing the fork's contact engine
> (rigid plane, NTS penalty, mortar/ALM, friction, viscous stabilization, explicit
> SOFT) to **`-ndm 2` plane models**. It is subordinate to the
> [[48_ladruno_contact_capstone_adr]] capstone (architecture + contracts) and
> follows the lane discipline of [[39_ladruno_contact_domain_adr]] /
> [[41_ladruno_mortar_alm_contact_adr]] / [[57_ladruno_edge_edge_contact_adr]]:
> oracle-first numpy gates, loud refusal over silent wrongness, one shared
> geometry per primitive. Implementation is a separate, later session; this
> document fixes the kinematics, the API/dimension-routing decisions, the
> thickness convention, the phase gates, and the scope fence.
>
> **Rev 2 (2026-08-18):** a 3-reviewer adversarial panel (source-evidence /
> mechanics / plan lenses, all FIX-FIRST) corrected rev 1 in nine load-bearing
> places. See §Review record at the end for what changed and why.

## What

The shipped **pair (NTS) and mortar contact lanes are 3D-only and enforce it
loudly**: they pre-flight their surfaces and **abort** if any node has fewer
than 3 coordinates (`ladrunoSurfaceNodesOk`,
`SRC/analysis/handler/LadrunoContactHandler.cpp:100` — "contact needs -ndm 3",
ADR-78 P1; called only from the NTS/mortar/tie blocks at `:705-706`,
`:1128-1129`, `:1322-1323`). The **rigid-plane lane skips that pre-flight** and
is already reachable from a 2D deck today with a zero-padded 9-arg
`contactPlane` (§How/8) — untested, undeclared, and frictionless, but wired. A
2D plane model — the workhorse for plane-strain geomechanics, wall-panel
pounding, masonry joint studies, and cheap parameter sweeps on the fork's own
plane-element family ([[25_ladruno_plane_elements_adr]],
[[26_ladruno_plane_frontier_adr]], [[70_ladruno_plane_finite_triangles_adr]]) —
cannot declare any *surface-to-surface* contact.

This ADR designs the 2D lane:

- **2D NTS** — slave node vs 2-node master **line segment**: closed-form
  point-on-segment projection (no Newton), perpendicular normal with an
  explicitly *new* orientation-vote design (§How/2 — the shipped 3D default is
  a per-pair datum with a recorded ridge-flip hazard; 2D does not inherit it
  blindly), a **vertex policy** for corner contact (§How/1 — the 2D analogue of
  3D's edge-edge gap), penalty enforcement, consistent tangent, auto-kn from
  2-DOF/node plane elements.
- **2D friction** — the tangent space is **one-dimensional**, so both the NTS
  return map (ADR-39 P3) and the mortar return map (ADR-41 C3) collapse to a
  scalar stick/slip test (same displacement-not-position slip crux, same
  commit-cycle state). NTS friction lands in T2; **mortar friction lands in T3
  with its own gates** (it is promised here, so it is built and gated here).
- **2D mortar** — the 3D facet **polygon clip collapses to a 1D interval
  overlap**, integrated on the **slave-trace measure** with the explicit
  alignment Jacobian (§How/3 — the 2D analogue of the 3D kernel's
  `A_phys = A_aux/cos_t`): D/M segment matrices, weighted nodal gaps,
  commit-cycle ALM `λ_N`, and the mesh-tying variant.
- **Parity features** — explicit SOFT=1 Courant-stable penalty, `-visc` normal
  dashpot, the ADR-78 removal lane (pruning), the warn-latch and abort
  discipline. The rigid-plane dashpot and the mass helpers carry over for free
  (size-guarded, `LadrunoContactFE.cpp:435-460`); the **SEGMENT/MORTAR
  dashpots and the SOFT tangential sizing are stride-3 hardcoded and are
  ported, not inherited** (§How/6).
- **Rigid plane in 2D** — the `RIGID_PLANE` FE mode is **already
  ndm-parameterized** (`LadrunoContactFE.cpp:63-68`; gap/residual loops run
  `d < ndm`; the handler derives ndm per node,
  `LadrunoContactHandler.cpp:1446-1489`). T0 turns "appears wired to work"
  into a tested, declared capability.

**Taxonomy note (corrected in rev 2):** in 2D the narrow-phase taxonomy is
**two primitives plus a vertex policy**. The 3D trichotomy is vertex-face
(NTS), face-face (mortar), edge-edge (ADR-57); in 2D there is no edge-edge
analogue, but **vertex-vertex/corner contact is not dismissible**: a
penetrating slave node whose closest master feature is a shared endpoint
projects out-of-bounds on *both* adjacent segments and would be refused by
both — a corner Voronoi wedge of finite measure at finite penetration and
finite dt (explicit pounding of building corners, masonry joints — this ADR's
own headline applications). §How/1 designs the vertex policy; G-T1a gates it.

**Not in scope (fenced):**

- **Axisymmetric contact** — the integration measure (`r dθ` weighting) and the
  hoop coupling are a different formulation, not a parameterization. Deferral
  row, demand-driven.
- **Parallel/DDM** — the contact subsystem is serial-only (handler
  `sendSelf/recvSelf` are P1a stubs; [[78_ladruno_parallel_contact_adr]] owns
  that lane). 2D changes nothing here.
- **2D runtime body discovery / re-emission** — ADR-55/60 stay 3D-only until the
  2D substrate exists; their seams are dimension-agnostic by design so nothing
  here forecloses them.
- **2D nodal-normal smoothing** (ADR-63 analogue: node-averaged segment
  normals) — the shipped `-smoothNormal` machinery (opt-in,
  `LadrunoContactHandler.cpp:844+`) stays 3D-only; deferral row. The 2D vertex
  policy (§How/1) covers the corner-*contact* need without it.
- **SOFT=2 segment-based explicit penalty** (B2 analogue) — in 2D it would ride
  the T3 mortar interval exactly as B2 rides the 3D clip; deferred until T3
  lands and an explicit-2D demand exists.

## Why

- **The workaround is expensive and partial.** Today the only surface-contact
  route is a one-element-thick 3D slice with fixed out-of-plane DOFs. That
  reproduces plane strain but (i) multiplies DOFs ~2-4x and forces 3D
  solvers/elements on a 2D problem, (ii) **cannot** reproduce plane stress (the
  slice is kinematically plane strain; plane-stress condensation lives in the
  2D element family — and because this is a headline reason, a plane-stress
  deck is **gated** in T1b, not just claimed), and (iii) makes every quick 2D
  study — the main reason 2D exists — pay 3D cost.
- **The fork invested in 2D elements and left them contact-blind.** ADR-25/26/70
  built a serious plane-element family (linear + finite-strain triangles,
  quads). Foundations-on-soil, retaining walls, masonry joints, and pounding
  studies in 2D all want contact; the abort at
  `LadrunoContactHandler.cpp:100` is currently the end of the road for
  surface-to-surface.
- **Everything is simpler in 2D — this is a low-risk lane on a proven seam.**
  The projection is closed-form (the 3D Gauss-Newton, its detK guard, and the
  tolR/tolP escape hatch all evaporate); the mortar clip is interval
  intersection; friction is scalar. The substrate — bucket-sort broad phase,
  `LadrunoContactDomain` committed state + GC + warn latches,
  `LadrunoContactFE` adapter modes, the handler's abort discipline — is
  reused, not rebuilt.
- **The 3D-only guard exists for a reason we must keep honoring.** The pre-3D
  guard incident (`LadrunoContactHandler.cpp:75` comment): a 2D node in a
  `-ndm 2 -ndf 3` frame passed the ndf checks and `getCrds()(2)` read out of
  bounds — nondeterministic garbage geometry. The fix is NOT to delete the
  guard but to replace "must be 3D" with **"must be dimension-consistent, and
  the dimension routes the lane"** — with the guard's **exact-equality
  contract kept** (§How/8): the shipped code refuses `ndf != 3` at six sites
  (`LadrunoContactHandler.cpp:951/983/1149/1172/1337/1349`) because
  `FE_Element::setID()` packs each node's *full* DOF_Group ID sequentially — a
  node with extra DOFs would push rotation equations into the contact's myID
  slots and shift every later slot: silent mis-assembly. 2D NTS/mortar
  therefore requires **`ndf == 2` exactly**, not `ndf ≥ 2`.

## Where

- **New code:**
  - `SRC/domain/contact/LadrunoContact2DKernel.h` — header-only, OpenSees-free:
    closed-form point-on-segment projection, perp normal, **vertex-pair
    evaluation** (§How/1), gap, 2D B-operator, and the 2D mortar
    interval-overlap integrator (D/M/ḡ on the slave-trace measure), plus a 2D
    `normalOriented` sibling for auto-kn (§below). Mirrors the
    `LadrunoContactProjection.h` / `LadrunoMortarKernel.h` discipline (raw
    `double[]` + `<cmath>`, numpy-oracle testable build-free). Kept separate
    from the 3D headers rather than templated — the 3D code is shipped and
    gate-protected; we do not churn it for symmetry.
  - `Ladruno_implementation/contact_prototypes/proto_t*_contact2d.py` — the
    numpy oracles, one per phase gate.
  - `tests/_testbed/contact_dump.py` — **the T0 byte-identity harness** (§Phase
    plan): runs 7 canonical 3D contact decks (rigid plane, NTS, NTS auto-kn,
    NTS friction, mortar ALM, mortar tie, edge-edge) and writes full-precision
    (hexfloat) nodal displacements + the per-lane queries to a physics-only
    text artifact for bit-for-bit diffing across builds (the existing battery
    is pass/fail-with-tolerance and cannot certify byte-identity — panel
    finding). Note `ladrunoContactForce` is **NTS-only** (fed from the SEGMENT
    branch via `setNtsForce`); the other lanes use their own queries, and the
    rigid-plane lane has none — its force checks go through the penalty-depth
    identity `kn·Σ⟨−g⟩₊ = P` (T0-A erratum).
  - `tests/test_adr85_contact2d_*.py` — the OpenSeesPy battery additions.
- **Modify (extend, do not duplicate):**
  - `SRC/analysis/handler/LadrunoContactHandler.cpp` —
    `ladrunoSurfaceNodesOk` becomes dimension-consistency (§How/8); the
    NTS/mortar loops grow a 2D branch; the six `getNumberDOF() != 3` guards
    become `!= ndm` (equality preserved); **auto-kn**
    (`ladrunoResolveAutoKn`, `:140-189`) needs more than the advertised
    `3*nn → ndm*nn` generalization: its coordinate gather is an unguarded
    3-loop (`:148-149` — the same `getCrds()(2)` OOB class as the §Why
    incident) and it calls the tri/quad-only
    `LadrunoContactProjection::normalOriented` (`:156`) — both get 2D
    siblings, since `-kn auto`'s failure path is FATAL (`:1049-1055`) and is
    the first thing a 2D deck would hit; the SEGMENT/MORTAR `-visc` dashpot
    loops (stride-3, `:700-703` and the `addCtoTang` twins) are ported in
    T1b/T3.
  - **Broad-phase staging stays stride-3**: `LadrunoContactBucketSort.h`
    indexes `P[i*3+{0,1,2}]` unconditionally (`:73-80`, `:123-125`), and its
    degenerate-axis handling is safe (`if (hi <= lo)` collapses NZ to 1,
    `:180`) — so in 2D the handler fills `segCoords` with an **explicit
    z = 0 pad** (making the fill at `:783` dimension-conditional) and the Grid
    is reused untouched. Rev 1's "parameterize staging on ndm" applies to the
    FE-facing arrays only, NOT to `segCoords`. (`segDef`/`segRef`/`slaveRef`
    live inside the opt-in `-smoothNormal` block, `:844+`, which 2D fences
    out — untouched.)
  - `SRC/analysis/handler/LadrunoContactFE.{h,cpp}` — `SEGMENT` and `MORTAR`
    modes gain 2D construction paths. The load-bearing hardcodes are the
    **SEGMENT ctor** (`:91-93`: `FE_Element(..., 3*(1+nps))`,
    `resid/tang(3*(1+nps))`, `ndm(3)` fixed at construction) and the fixed
    stack buffers (`:674` — `B[15]`); `softKt` gets an explicit 2D branch
    (§How/6). `RIGID_PLANE` needs no change. No new Mode enum values and no
    new class tag — 2D is a *parameterization* of the existing modes, recorded
    in `LEDGER_implementations.md`.
  - `SRC/domain/contact/LadrunoContactSurface.{h,cpp}` — accept
    `nodesPerSeg = 2`; topology storage is already dimension-agnostic.
  - `SRC/interpreter/OpenSeesOutputCommands.cpp` — `contactSurface` accepts
    `-master 2` / `-slave-segments 2` **only when the referenced nodes are 2D**
    (declaration-time check; nodes already must exist at declaration, ADR-78);
    `contactPlane` gains the 2D form (§How/8); queries document their 2D
    meaning. This is an upstream file with existing `LEDGER_vanilla_files`
    rows — each phase that touches it owes its ledger row + `// Ladruno`
    marker **in that phase's PR**.
  - `SRC/domain/contact/LadrunoContactDomain.{h,cpp}` — `FrictionState` and
    the tie slots are reused **as 3-wide arrays with the third component held
    at zero** (they are hardcoded 3-vectors: `lambdaTie[3]` / `rtGlobal[3]`,
    `LadrunoContactDomain.h:339-340`, and `accumulateMortarTie(...,
    const double rFacet[3], ...)`, `:373-374` — NOT "dimension-blind" as rev 1
    claimed); `MortarNormalState` scalar λ_N reused as-is; `RigidPlane` struct
    unchanged (`p0[3], n[3]`, zero third component in 2D).
  - Banner (`Ladruno_scripts/banner_features.txt` + `patch_banner.py`) and the
    three ledgers — **in the same PR as each shipping change** (repo rule):
    the banner's existing contact line gains a `2D` clause in the first phase
    that ships a user-visible capability (T0), and T4 only finalizes wording.
- **Reference (copy patterns from):**
  - `LadrunoContactProjection.h` — the in-bounds tolerance slack (1e-9), the
    refusal-status idiom, the header discipline.
  - `LadrunoMortarKernel.h` — `PairResult` + degenerate-overlap refusal
    (`status = -1`) + the slave-measure Jacobian precedent
    (`A_phys = A_aux/cos_t`, `:460`) + the min-area sliver guard
    (`A_poly < 1e-12·A_min`, `:444`).
  - `LadrunoContactFE.cpp` SEGMENT/MORTAR — stateless-adapter assembly, softKn /
    gapModeInvMass (size-guarded for 2-DOF nodes), viscousActive, warn-latch.
  - ADR-57 E7's **ordered ownership protocol** (NTS ≻ edge ≻ mortar stand-down,
    handler `:1179-1184`) — the routing precedent §How/3 adapts.
- **Build:** no new CMake target, no external dep.

## How — governing kinematics and decisions

### 1. 2D NTS: closed-form projection, vertex policy, B-operator

Master segment `x(ξ) = (1−ξ)·X₀ + ξ·X₁`, `ξ ∈ [0,1]`, slave point `x_s`
(current = reference + trial displacement, as in 3D):

```
t = X₁ − X₀,   L² = t·t          (refuse L < τ_seg·L_min_ref — zero-length segment;
ξ* = (x_s − X₀)·t / L²            L_min_ref = min initial segment length of the surface,
                                  computed once at handle(); τ_seg default 1e-8)
in-bounds iff ξ* ∈ [−δ, 1+δ],  δ = 1e-9 (parametric slack, as 3D)
OUT-OF-BOUNDS ⇒ REFUSE the pair (as 3D evalSegment) — never clamp
n  = perp(t)/L = (−t_y, t_x)/L · σ     (σ = ±1, the orientation vote — §2)
g  = (x_s − x(ξ*)) · n                  (penetration ⇔ g < 0, fork convention)
```

The parametric coordinate is dimensionless by construction, so the 3D
projection's Newton/detK/tolR history is structurally absent — **but the gap
is not scale-free**: `g` is a difference of large coordinates and keeps the
`ε·|x|` absolute noise floor (≈1e-11 at coordinates ~5e4). Oracle tolerances
are therefore **relative** (`·L_seg` for parametric quantities, `·|x|_scale`
for gaps), never bare 1e-12 absolutes; the far-from-origin falsifier verifies
the *relative* bound.

**Vertex policy (new in rev 2 — the 2D corner primitive).** Refusing
out-of-bounds projections leaves the corner Voronoi wedge uncovered: a slave
node penetrating past a convex master corner projects out-of-bounds on *both*
adjacent segments. The kernel therefore also evaluates **vertex pairs**: when
the slave's closest master feature on a surface is a shared endpoint `X_v`
(both adjacent projections out-of-bounds on the vertex side), activate

```
n_v = (x_s − X_v)/‖x_s − X_v‖ · σ_v     σ_v = committed side sign (see below)
g_v = ‖x_s − X_v‖ · σ_side               (sign from the committed side state)
B_v = [ n_v | −n_v ]                     (slave + the one vertex node; ndof = 4)
```

with (a) the **sign committed at first capture** from the wedge of the two
adjacent segment normals — never re-derived from the current position (the
ADR-57 committed-sign lesson: a position-derived sign flips exactly while
interpenetrating); and (b) **shared-vertex ownership dedup**: a slave node is
owned by *either* a segment pair *or* a vertex pair per surface per step,
never both. G-T1a gates corner impact, slide-around-corner normal continuity,
and the no-double-count property.

**T1a oracle corrections (2026-08-18) — the shipped kernel contract, which
T1b wires verbatim:**

1. **The rev-2 premise was inverted.** A **convex** corner's exterior wedge is
   strictly non-penetrating (20000/20000 random wedge points, `gap > 0`,
   cross-checked against a half-plane material test), so a convex vertex pair
   **never carries force**; the convex corner's *interior* is covered twice by
   the adjacent segment strips — a double-count, handled by ownership, not a
   hole. The real hole is **concave (re-entrant)**: no exterior wedge exists,
   and the material-side region behind the corner is refused by *both*
   adjacent segments while being entirely penetrating — an unowned arc of
   length `R·|turn|` at penetration depth R (19.7% of a sweep path at
   R = 0.25, scaling linearly). **The handler instantiates vertex pairs at
   concave corners only.**
2. **Ownership is an ordered rule, not the slacked XOR.** "Segment in-bounds
   XOR wedge" is algebraically right but not floating-point closable:
   `xi_prev − 1` resolves an overshoot only to ulp(1) while the wedge
   predicate resolves it exactly (measured disagreement 1.15e-15 parametric,
   5.2 ulp; a 16400-station seam scan gave the XOR 80 double-owners and 1
   zero-owner; no tolerance closes a resolution mismatch). Shipped rule
   (ADR-57 E7 style): **segment precedence in surface order, then the
   unslacked vertex claim (`aPrev > 0 && aNext < 0`), else no owner** — total
   and unique over the full scan; in the seam overlap both candidates are the
   same contact (relative gap difference ≤ 1.7e-14, normal disagreement
   within its analytic `atan(tolIn·L/|gap|)` bound).
3. **The vertex side sign** comes from the adjacent-normal **bisector**, with
   the conditioning gate on `|r|` alone — NOT `|r|·|bis|`, because the
   bisector is a cancelling sum of unit vectors carrying ~eps *absolute*
   noise (the scaled gate produced 707 wrong signs in 48902 wedge hits; the
   `|r|` gate produced 0). `sign((xs−Xv)·n_prev)` is valid only for turns
   ≤ 90° and is not used.
4. The 2D B-operators FD-gate as the **exact** first variation of the gap
   (the `∂n/∂u`, `∂ξ/∂u` contributions cancel identically in 2D) — this is
   the proof behind §How/5 keeping `kn·BᵀB` as the main term.
5. **T1b gate correction (2026-08-18): the ownership contract had an
   open-chain-end hole.** Step 4 ("else nobody") was oracled on *interior*
   seams only; at a surface's **terminal ends** a hard refusal on a pair
   carrying O(P) force is a discontinuity no consistent-tangent Newton
   survives (probe-traced: master tilt pushes a boundary slave's ξ to
   −1.25e-7, both end pairs refuse, the block goes singular or limit-cycles —
   the measured patch/slice-twin stall). T1b interim: a **terminal-side
   acceptance window** (`NTS2D_END_SLACK = 1e-3`, parametric, on unchained
   ends only; gap = signed distance to the infinite line, B stays the exact
   first variation, so residual/tangent are consistent and C0 across the
   band; the cliff moves to 1+1e-3, disclosed). Interior seams keep the
   strict ordered rule untouched. The honest permanent treatment is a
   **radial end-cap vertex pair** at open terminals (C0, no window) — named
   as a T2+ item in LEDGER_quirks. Also from the same gate round: the vertex
   coincidence floor uses the **conditioning gauge** (`τ_perp·Lref`), not
   `τ_seg·Lref` — the two relative gauges collided with the standard 1e-8
   seed penetration (the concave singular-matrix failure).
6. **T1b panel correction (2026-08-18): segment chaining is an ENFORCED input
   contract, not a convention.** `prevFar`/`nextFar` arming assumes the
   master's 2-node segments are declared head-to-tail in surface order; the
   panel probe showed a *permuted* listing silently produced the exact
   double-owner answer (2× stiffness in the end-window band) and a *reversed*
   listing disarms the stand-down and drops every concave vertex pair — the
   T1a hole reopened silently. handle() therefore runs an O(nSeg)
   chain-integrity scan: any node shared by two segments must be a
   consecutive head-to-tail link (or the closed-loop wrap); anything else is
   a named FATAL (`ADR-85 chain`) with re-order/re-wind guidance. This is a
   user-facing declaration contract and goes in the T4 user guide.

Shape functions `N = [1−ξ*, ξ*]`; segment B-operator over `[slave | X₀ | X₁]`:

```
B = [ n | −N₀·n | −N₁·n ]        ndof = 2·(1+2) = 6
r = Bᵀ·tn,  tn = kn·⟨−g⟩₊        K_c = kn·BᵀB  (+ geometric ∂n/∂u term — §5)
```

Active-pair bookkeeping, per-(slave, segment) committed state, and the B3
`setNtsForce` snapshot reuse the 3D structures.

### 2. Normal orientation: a NEW 2D default, not a reuse (rev 2 correction)

Rev 1 claimed the 2D lane "reuses the handler's reference-configuration
centroid vote". The panel refuted this: that machinery
(`LadrunoContactHandler.cpp:845-900`) lives inside the **opt-in
`-smoothNormal` block** (ADR-63, off by default). The shipped engine actually
has **three** sign mechanisms:

1. **NTS default** — a *per-pair* datum `orientDir = slave_ref −
   segment_ref_centroid` (`:1017-1027`), annotated in-source as the "R3-buggy
   datum (it flips on a ridge)"; it also **degenerates silently** for a slave
   node lying ON its master segment (flush, zero-gap interface — exactly the
   masonry-joint / footing-on-soil decks this lane is for): the datum is
   tangential, `normalOriented`'s fail-safe returns false, and the pair is
   silently inactive.
2. **`-smoothNormal` global seed** — the interface-level reference vote
   (`:865-886`), opt-in, 3D-only (built on Newell facet normals).
3. **Mortar** — its own per-facet-pair centroid direction (`:1188-1194`) plus
   the WINDING-LUCK fix + `WARN_MORTAR_ORIENT` latch (`:1199-1226`).

**2D decision:** the 2D NTS and mortar lanes get a **new, interface-level
reference vote** (computed once per contact at handle() from reference-
configuration surface centroids — the *idea* of mechanism 2 without its 3D
Newell machinery), with:

- **named abort** when the vote is degenerate and no `-outward ox oy` is given
  (flush interfaces with coincident centroids are genuinely ambiguous — refuse,
  never guess; this is NEW behavior relative to the 3D default's silent drop,
  and is called out as such). **T1b panel correction (2026-08-18): degeneracy
  needs a MAGNITUDE gate, not only an angle gate.** The centroid datum must
  satisfy `|refDir| > 1e-6·Lref` — a flush deck seeded the fork-standard 1e-8
  into penetration gives `refDir = (0, −1e-8)`, perfectly conditioned in
  angle, and the vote then resolves *inverted off the seed's sign*: every
  pair inert, force ≡ 0, run converges (probe-confirmed — the ADR-78 P0
  signature through the lane's own safety mechanism). Below the magnitude
  gate the vote falls into the named refusal. A user `-outward` keeps the
  angle-only gate (an explicit direction is a deliberate input);
- `-outward` as the explicit override (2 components in 2D; the parser currently
  reads exactly 3 doubles, `:782` — a 2D branch);
- σ **fixed at pairing**, re-evaluated only on re-pairing (ADR-57 lesson);
- **flush 2D interfaces effectively require `-outward`** unless the centroid
  vote resolves — stated in the user guide, not discovered at runtime.

### 3. 2D mortar: interval clip on the slave-trace measure + explicit routing

Slave segment S and master segment M, both 2-node. Project S's endpoints onto
M's parametrization (closed form, §1); intersect with `[0,1]`:

```
[a,b] = [ξ_M(s₀), ξ_M(s₁)] ∩ [0,1]  (sorted)
b−a ≤ τ_ov ⇒ status = −1 (empty)     τ_ov relative to min(L_s, L_m)/L_m —
                                      the 2D analogue of the 3D A_poly < 1e-12·A_min
                                      guard; a fine-slave/coarse-master pair has
                                      b−a ≈ L_s/L_m legitimately small (panel catch)
```

**Integration measure (the load-bearing sentence rev 1 omitted):** the
integrals run on the **slave-trace measure**, `dΓ_s = (L_m/|t̂_s·t̂_m|)·dξ_m` —
the 2D analogue of the shipped kernel's `A_phys = A_aux/cos_t`
(`LadrunoMortarKernel.h:460`). Omitting the `1/|t̂_s·t̂_m|` Jacobian would bias
D/M/g̃ by cosθ for every non-parallel pair — and the flat patch test could
never catch it (cosθ = 1 there), so the T3 oracle includes a tilted
non-parallel-overlap case. With the measure stated, the exactness claim
holds: the map ξ_s ↦ ξ_m is affine for straight segments, n and the Jacobian
are constant per pair, so all integrands (D_ij, M_ij, g̃_i) are degree ≤ 2 and
**2-point Gauss is exact**.

**Routing (rev 2 correction — rev 1's "refused to the NTS lane, which covers
it" assumed a fallback that does not exist).** In the shipped engine `-mortar`
and NTS are separate declared lanes; a mortar `status = −1` is silently
skipped. 3D closed its equivalent hole (ADR-57 E7) with an ordered ownership
protocol — NTS ≻ edge-edge ≻ face-mortar with stand-down
(`LadrunoContactHandler.cpp:1179-1184`). The 2D lane **adopts that protocol
minus the edge tier**: NTS ≻ mortar with stand-down (no double-count when both
are declared on the same surfaces), and an alignment-refused (near-
perpendicular) mortar pair is covered iff an NTS contact is declared on the
same surfaces. If it is NOT, the refusal is **loud**: a new
`WARN_2D_PERP_NO_NTS` engine latch fires naming the pair — never a silent
zero. G-T3(d) gates the protocol, not a fictional fallback.

**Friction (ADR-41 C3, in scope here — rev 2):** the mortar return map over
the soft/ALM pressure collapses to the same scalar as §4, with the unified
cone `min(μN+c, τ_max·h·ℓ)` clamps h-scaled per §7. Gated in G-T3 (mortar
block-on-incline + μ=0 byte-identity), not left orphaned in the What section.

Downstream is structure-preserving reuse of the ADR-41 machinery: nodal
weighted gaps on the global accumulators, KKT active set, commit-cycle ALM
`λ_N` via `MortarNormalState` Uzawa-on-commit, `-tie` on the `lambdaTie[3]` /
`rtGlobal[3]` global accumulators (3-wide, third slot held at zero — §Where).
`analyze_augmented` / `ladrunoBeginAugment` work untouched.

### 4. 2D friction: the return map is a scalar

The tangent space at a 2D contact is `span{t̂}`, `t̂ = perp(n)` — one
dimension. The ADR-41 covariant-metric machinery (`R = g·r`, 2×2 metric)
collapses: metric `g₁₁ = L²`, one slip variable. Decision: a scalar path
(`returnMap1D`) beside the 3D kernel in the same header, sharing the constants
and the τ_max/cohesion clamps. The two cruxes carry over verbatim:
displacement-not-position slip; committed stick state in `FrictionState`
(second tangential slot ≡ 0).

**Oracle role (rev 2 correction):** rev 1 called the plane-constrained parity
test "the referee" for scalar-vs-degenerate-3D. The panel showed it has **zero
discriminating power** — on a plane-constrained configuration the second
covariant direction vanishes, so a correct scalar path, a degenerate 3D path,
*and a wrong-metric 3D path* all pass. The decision is therefore made **before
T2** by a recorded 30-minute numpy experiment (evaluate the 3D kernel
degenerately; keep whichever is simpler to maintain if both are exact), and
the plane-constrained test is demoted to what it actually is: an
equivalence/regression gate.

### 5. Consistent tangent (rev 2: correct flag names; T1b correction)

`K_c = kn·BᵀB` — and in 2D this main term is **gradient-exact**: the T1a FD
gate proved the B-operator is the exact first variation of the gap (the
`∂n/∂u`, `∂ξ/∂u` contributions cancel identically in 2D), so `kn·BᵀB` is the
consistent force-linearization on a *fixed* geometry.

**T1b correction (2026-08-18):** the phase-plan inference "therefore the
`-geomtan` geometric term is zero in 2D" was an over-read — first variation
exact ≠ second variation zero. The Newton **curvature block** (master-segment
rotation, `g` nonlinear in the master DOFs) is nonzero — the same class the 3D
B3 term assembles. T1b ships `-geomtan` on a 2D pair as **accepted,
documented no-op** (identical to the 3D *default* path: never wrong forces,
only Newton rate on rotating masters); the true 2D curvature term is
**deferred, demand-driven** (a Hertz-style curved-master 2D benchmark is the
trigger). G-T1b(e) is therefore a **byte-identity gate** (± flag identical),
not an iteration-count gate.

Flag mapping, per the shipped parser:

- geometric **normal** term → **`-geomtan`** (`OpenSeesOutputCommands.cpp:809`)
  — symmetric, any solver in 3D; accepted no-op in 2D (above);
- the **non-symmetric consistent friction** tangent → **`-consistanttan`**
  (`:798`) — keeps its shipped unsymmetric-solver warning.

(Rev 1 conflated the two and would have gated the symmetric normal term behind
the unsymmetric-solver flag — a capability regression the panel caught.)
Explicit path never forms K.

### 6. Explicit SOFT, viscous, removal: what is free vs what is a port (rev 2)

- **Free (verify, don't port):** `softKn` normal sizing — `gapModeInvMass` /
  `ladrunoNodeMass` / `ladrunoInvMassProj` are size-guarded for 2-DOF nodes
  and, under the §How/8 `ndf == ndm` guard, slot 2 is structurally zero.
  Rigid-plane `-visc` (`d < ndm` loops, `LadrunoContactFE.cpp:660-670`).
  Removal-lane pruning (`pruneMissingNodes` whole-segment drop holds for
  nps = 2).
- **Ports (rev 1 wrongly claimed these free):**
  - SEGMENT/MORTAR `-visc` dashpots: stride-3 hardcoded (`:700-703` + the
    `addCtoTang` twins) — ported in T1b (SEGMENT) and T3 (MORTAR).
  - `softKt` (B1-kt tangential sizing): the shipped code always builds two 3D
    tangents from the least-aligned coordinate axis (`:553-566`) — for an
    in-plane n it selects the OUT-OF-PLANE direction as t1. It needs an
    explicit 2D branch with the single tangent `t = perp(n)`. The rev 1 claim
    "one tangent direction makes it exact" is the *goal* of that branch, and
    only per-mode: for anisotropic nodal mass (m_x ≠ m_y) the coupled kn+kt
    system's ω²_max can approach the sum of the per-mode targets — the 3D B2
    coupled-K_c caveat (SOFSCL > 0.25 warn) transfers to 2D verbatim rather
    than vanishing.
- `k_soft = SOFSCL·4·m_eff/dt²` itself is a dimension-independent Courant
  argument — unchanged.

### 7. Thickness convention (rev 2: full enumeration + no double-scaling)

2D plane elements assemble **actual forces** (thickness h is inside their
stiffness). Contact must match:

- **NTS penalty**: kn is force/length on a nodal gap — as in 3D, where kn was
  never a pressure. Auto-kn from the element stiffness absorbs h automatically
  (the K nodal blocks contain it). No thickness parameter.
- **Mortar**: interval integrals produce force **per unit thickness**.
  Declaration gains `-thickness h` (default **1.0**, loudly documented),
  applied **once, where each density enters the residual**, to: `ε_N`, `ε_T`,
  the viscous `μ_c`, the friction stress clamps `-tauMax` / `-cohesion` (these
  are the h-sensitive ones — an unscaled clamp makes the stick/slip threshold
  wrong by exactly h against the h-carrying pressure), and the tie stiffness.
  **`λ_N` is NOT separately scaled** — it accumulates from `ε_N·ḡ` in the
  Uzawa update and inherits h; scaling it again would be an h² error (panel
  catch). Mismatched slave/master element thicknesses: the declared h is
  authoritative; warn-latch if a discoverable element thickness disagrees.
- Queries: `ladrunoContactForce` returns force **(NTS lane only** — it is fed
  exclusively from the SEGMENT branch; rigid-plane/mortar/tie/edge decks must
  use their own queries or the penalty-depth identity); `ladrunoMortarPenetration`
  is geometric (length) — dimension-blind.

### 8. API and dimension routing (rev 2: one oracle, not two)

**One command surface; the dimension oracle is the referenced surface's node
coordinates — everywhere, including `contactPlane`.** Rev 1 proposed branching
`contactPlane` on `OPS_GetNDM()`; the panel refuted it: interpreter-state ndm
is mutable (repeat `model` commands don't wipe the domain), has no
`ndm ∈ {0,1}` story, and on the classic Tcl path (`SRC/tcl/commands.cpp`, the
OpenSeesMP route) resolves through a static `TclModelBuilder*` with **no null
guard** (`elementAPI_TCL.cpp:1355`) — a builder dependency the contact parser
has never had. Instead:

- `contactPlane` derives its dimension from `getSurface(slaveSurfTag)`'s node
  coordinates — the same oracle `contactSurface` uses. The **9-arg zero-padded
  form stays permanently valid** (it is what already works today); a 7-arg 2D
  form (`contactPlane tag ssTag nx ny px py kn <opts>`) is accepted as pure
  ergonomics when the surface is 2D, and cannot disagree with anything.
- `contactSurface tag -master 2 ...` — `nodesPerSeg = 2` accepted **iff every
  referenced node has 2D coordinates** (declaration-time; nodes must already
  exist). `nps ∈ {3,4}` iff 3D. The current `>= 3` parser check
  (`OpenSeesOutputCommands.cpp:371`) doubles as a 3D sanity gate — its
  replacement must keep a 3D deck from declaring 2-node segments.
- `-outward` takes 2 components in 2D (parser branch; currently reads 3).
- **Guards** (`ladrunoSurfaceNodesOk` successor): every node of both surfaces
  must have the **same** `getCrds().Size() ∈ {2,3}`; that size must match the
  declared segment arity (2 ⇔ nps 2; 3 ⇔ nps 3/4); and on the NTS/mortar/tie
  lanes **`ndf == ndm` exactly** (the shipped equality contract — §Why). Only
  the RIGID_PLANE lane keeps `ndf ≥ ndm` (its adapter couples the first ndm
  DOFs by construction). Any violation is a **named abort**.
- **Guard-contract precision (T0 review):** the `ndf == ndm` equality is a
  **multi-node-FE packing requirement** — `FE_Element::setID()` packs each
  node's full DOF_Group sequentially, so extra DOFs shift *later nodes'*
  slots; that shifting needs ≥2 DOF_Groups. The single-node RIGID_PLANE
  adapter deliberately keeps `ndf ≥ ndm`: 3D shells-on-plane (`ndf 6`) rely
  on it today, `setID` truncates with an upstream warning, and the retained
  translational equations are correct. A 2D `ndf 3` deck on the plane lane is
  therefore accepted-with-warning, exactly like its shipped 3D analogue.
- **Disclosed over-refusal (T0 review):** the declaration-time
  dimension-consistency rule applies to **all** surface kinds, including
  `-slave` node sets whose only consumer is the per-node-capable rigid-plane
  lane. A mixed 2D/3D `-slave` set that previously reached `contactPlane` is
  now refused at declaration — deliberate: the set's future consumers are
  unknown at declaration time, and a mixed set on any pair lane is the
  historical OOB reproducer. Loud over-refusal beats a consumer-dependent
  contract (the ADR-78 P2 disclosed-over-refusal precedent).

## Phase plan (oracle-first, gated; T = two-D lane)

Sessions follow [[ladruno-session-pacing]]: one phase per session, C++ batched
before one build (`OpenSees OpenSeesPy`, both mtimes + the `ladrunoBuild` stamp
checked; oracles under `python3.12`). **Every phase gate includes:** ledger
rows + `// Ladruno` markers for the files that PR touches, **in that PR**
(repo rule — rev 1 wrongly parked these at T4).

Gate language (rev 2): "battery green, **N passed with N unchanged** from the
pre-PR count, recorded in the PR body" is the checkable regression clause (the
shipped ADR-39/41/57 idiom); **byte-identity** claims run through the
`contact_dump.py` harness — full-precision observables from 7 canonical 3D
decks, diffed bit-for-bit against the pre-change binary (both `ladrunoBuild`
stamps recorded). Refusal falsifiers assert the **named message substring**,
not just a nonzero exit: a declaration refusal raises in-process in serial
(`return −1` → `OpenSeesError`) and tears the job down only under MPI
(`LadrunoContactAbort`), and a handler-time refusal surfaces as a failed
`analyze` — the child-process test pattern covers all three shapes uniformly
(T0-A erratum).

- **T0 — guards, declaration, rigid-plane acceptance, dump harness.**
  Dimension-consistency pre-flight (with the `ndf == ndm` equality contract);
  `nps = 2` declaration path; `contactPlane` surface-derived dimension + 7-arg
  form; **build `contact_dump.py` and capture the pre-change baseline**;
  banner clause + vanilla-ledger row for the parser file.
  **Gate G-T0:** (a) 3D battery N-unchanged + `contact_dump` bit-identical
  across the T0 build; (b) 2D block-on-rigid-plane statics + explicit drop:
  converges, holds `g ≥ 0`, reactions balance; SOFT + visc variants;
  (c) falsifier matrix: mixed-dim abort; `-ndm 2 -ndf 3` node on an NTS
  surface refused (equality guard — the historical OOB reproducer); 3D deck
  declaring `nps = 2` refused; 2D pair/mortar declarations get a **named
  not-yet-supported refusal** (T1b/T3) rather than silent acceptance;
  9-arg `contactPlane` on a 2D surface still works (back-compat). (The
  flush-interface degenerate-vote falsifier was misfiled here in rev 2 — the
  vote is built in T1b and its falsifier lives in G-T1b.)
- **T1a — 2D NTS kernel + oracle (header-only PR).** `LadrunoContact2DKernel.h`
  projection/normal/vertex-policy/B + `proto_t1_nts2d.py`. The header is
  included by **no TU** ⇒ the build is byte-identical **by construction** (the
  ADR-57 E0 precedent).
  **Gate G-T1a:** oracle parity at *relative* tolerances over randomized
  segments — clamped/out-of-bounds refusal, far-from-origin (relative-bound
  check), zero-length, reversed-winding, **corner falsifiers**: corner impact
  activates a vertex pair; slide-around-convex-corner hands off
  segment→vertex→segment with continuous normal; shared-vertex no-double-count.
- **T1b — 2D NTS wiring.** FE SEGMENT 2D path (ctor, buffers, dashpot port);
  handler 2D branch (z-padded `segCoords`, `ndm*nn` auto-kn + 2D
  `normalOriented` sibling); the new interface-level orientation vote;
  `softKt` 2D branch; the `contact -outward` **2-component 2D branch** (T0
  left the 3-double read — a 2D deck writing `-outward ox oy` today gets an
  unnamed 3D-signature error; T1b gives it the named 2D form) and the
  `-soft`-pre-scan refusal-ordering cleanup (T0 documents that a 2D `-soft`
  deck is refused by the `-soft` subsystem's message before the ADR-85 named
  refusal — loud but wrongly attributed; fix when this loop is rebuilt).
  **Gate G-T1b:** (a) 3D battery N-unchanged + `contact_dump` bit-identical;
  (b) 2D compression patch across a non-matching interface transmits the load
  (NTS tolerance); (c) **slice-twin gate, honestly stated**: first measure the
  element-formulation floor (2D plane vs one-layer brick with contact replaced
  by a tied/merged mesh — the floor is a recorded artifact), then require
  contact **observables** (summed normal force, per-node gaps) to agree within
  floor + 2%; (d) **plane-stress acceptance** (a plane-stress patch under
  contact vs its analytic through-thickness force — the §Why headline, now
  gated); (e) `-geomtan` gate — **restated at T1b** (§How/5 correction):
  byte-identity of results with and without the flag (accepted no-op; the
  iteration-count form was unsatisfiable once the no-op decision landed);
  (f) auto-kn sizes within the 3D-twin's factor; (g) flush-interface
  degenerate orientation vote → named abort (moved here from G-T0 — the vote
  lands in this phase). **T1b vote policy addition:** a per-surface sigma
  requires a unanimous nonzero vote across the master segments — a SPLIT vote
  (strongly curved or inconsistently wound master) is a named refusal with
  re-wind / split-the-surface / `-outward` guidance, never a guess; a
  >90°-curved master surface is declared as separate contacts in T1b.
- **T2 — NTS friction + explicit parity.** Pre-T2: the recorded numpy
  reuse-vs-scalar experiment (§How/4). Scalar return map; SOFT/visc/removal
  gate coverage.
  **Gate G-T2:** (a) block-on-incline stick/slip threshold exact at
  `tanφ = μ`; (b) sliding energy balance under explicit (dissipation =
  μ·N·slip); (c) SOFT pounding stable at `dt = 0.9·dt_cr`, SOFSCL default;
  (d) visc statics-inert byte-identity (via dump harness); (e) μ = 0
  byte-identical to T1b; (f) **removal gate**: `remove element` mid-run drops
  a 2D segment, contact continues on the survivors (rev 1 promised this and
  gated nothing).
- **T3 — 2D mortar + ALM + tie + mortar friction + thickness.** Interval
  integrator on the slave-trace measure + `proto_t3_mortar2d.py` (including
  the tilted non-parallel case that catches a missing Jacobian); ALM; `-tie`;
  the ownership protocol + `WARN_2D_PERP_NO_NTS`; C3 scalar friction; MORTAR
  dashpot port; `-thickness`.
  **Gate G-T3:** (a) the mortar headline: machine-precision uniform pressure
  transfer across a non-matching 2D interface; (b) `analyze_augmented` drives
  penetration below an ε_N-independent augTol; (c) tie patch exact on
  non-matching meshes incl. a shared-node T-junction (the ADR-41 C4 crux);
  (d) **routing protocol**: perpendicular pair + NTS declared ⇒ NTS owns it,
  no double-count; perpendicular pair + mortar-only ⇒ `WARN_2D_PERP_NO_NTS`
  fires (never silent); (e) mortar block-on-incline at `tanφ = μ` + μ = 0
  byte-identity; (f) a `-thickness 2` twin scales tractions by exactly 2 and
  leaves geometry (gaps) untouched.
- **T4 — battery, benchmark, ship.** 2D Hertz cylinder-on-plane **with the
  convention pinned**: `b² = 4P′R/(πE*)`, `p_max = 2P′/(πb)` where **P′ =
  applied force per unit thickness** (the deck uses unit-thickness elements so
  the applied F equals P′ numerically) and **E\* is the plane-strain
  combination** `1/E* = (1−ν₁²)/E₁ + (1−ν₂²)/E₂` — both written in the test
  docstring; full 2D battery; **one full Zone-A sweep**; final banner wording;
  user-guide section (incl. the 2D units/thickness table and the flush-
  interface `-outward` requirement).
  **Gate G-T4:** Hertz within the 3D lane's tolerance band; slice-twin re-run
  green; full-sweep N recorded; docs merged.

Estimated shape: **6 PRs** (T0, T1a, T1b, T2, T3, T4). T0 and T1a are
deliberately small; T1b is the only PR whose diff touches shipped 3D
assembly, and it contains nothing else.

## Execution plan — sessions, models, effort (rev 2: staffing re-argued)

Staffing principle, sharpened by the panel: tier follows **failure-mode
visibility**, not difficulty. The fork's booked worst-case is the silent
plausible-wrong run (ADR-78 P0: converged, balanced, wrong by 2×). A phase
whose bugs are *loud* (a machine-precision gate that fails in your face) can
run a tier lower than a phase whose bugs are *quiet* (a mis-parameterized
stride in shipped 3D staging that still converges).

| Phase | Model | Effort | Why this tier |
|---|---|---|---|
| T0 guards + declaration + harness | Opus | high | Grew in rev 2 (dump harness, equality-contract guard rewrite, banner/ledger). Parser + pre-flight edits on the path every 3D deck traverses; the refusal matrix is the contract. |
| T1a kernel + oracle | Opus | medium | Self-contained header + numpy oracle, byte-identical by construction; the vertex-policy math is the only subtlety. |
| T1b NTS wiring | **Fable** | high | The lane's highest *silent*-failure surface: shipped stride-3 staging, dashpot port, auto-kn siblings, the new orientation vote — every bug here is the quiet kind. (Rev 1 had Opus here and Fable on T3; the panel's inversion argument is accepted.) |
| T2 friction + explicit parity | Sonnet | high | Well-specified once the pre-T2 experiment (recorded artifact) closes the design question; every gate has an analytic answer. |
| T3 mortar + ALM + tie + friction | Opus | high | Crux-dense but **loud**: new code in an isolated header, and the machine-precision patch/tie/routing gates self-diagnose. (Rev 1's Fable slot moved to T1b.) |
| T4 battery + Hertz + ship | Sonnet | medium | Mechanical breadth plus one bounded calibration. |

**Escalation rules (mechanical, not self-policed — rev 2):**

- T2: the **second failing build** of the parity test — builds are stamped;
  the stamp of each failing build is recorded in `_adr85_t2_design.md` — hands
  the session to Opus high.
- T4: a Hertz miss escalates the *diagnosis* to Fable — a benchmark miss at
  the end of the lane means an earlier gate lied; finding which is a
  cross-phase reasoning task.
- Any phase: a `contact_dump` diff or a battery-count drop freezes the phase
  until root-caused; never "re-ran and it passed".

**Review gates:**

- Every PR: `/code-review` at **high** minimum (`max` where warranted;
  `xhigh` is not a code-review level — rev 1 error).
- T1b and T3: an adversarial multi-reviewer pass at **max** effort (or
  `/code-review ultra` for the cloud panel) before merge, **against a named
  contract table**: for T1b — the dump-harness diff, the six equality guards,
  the z-pad fill, the dashpot port, the orientation-vote degenerate abort; for
  T3 — the slave-measure Jacobian (tilted-case oracle), the ownership
  protocol matrix, the h-scaling list of §How/7, μ=0/tie byte-identity. At
  least one reviewer is prompted to REFUTE the byte-identity and patch-test
  claims.
- T0's review checks the refusal matrix (2D/3D × nps 2/3/4 × **ndf == ndm** ×
  vote-degenerate) line by line against §How/8.

## Risks / open questions

> [!question]
> **Thickness discovery.** Can the mortar lane cheaply discover the adjacent
> element thickness to *validate* the declared `-thickness` (warn-latch on
> mismatch), or is that coupling not worth it? T3 decides; the declared value
> is authoritative either way.

> [!question]
> ~~**Vertex-policy sign robustness / concave activation.**~~ **ANSWERED by
> the T1a oracle (2026-08-18), with the rev-2 premise inverted** — see the
> T1a corrections block in §How/1: convex vertex pairs never carry force,
> concave corners are the real hole and the only vertex pairs the handler
> instantiates; ownership is the ordered rule; the side sign is the
> `|r|`-gated bisector.

> [!question]
> **Does the T0 rigid-plane acceptance find real 2D bugs?** The lane appears
> wired (ndm-generic FE, per-node ndm derivation) but has never run in 2D. If
> T0 uncovers structural surprises, T0 grows; the gate stays.

- **Byte-identity risk in 3D** remains the top regression concern, now with a
  real observable: the `contact_dump` harness. The T1b panel checks its diff,
  not a prose claim.
- **`nodesPerSeg = 2` reaching 3D code**: any path computing `3*(1+nps)` or
  indexing `X[4][3]` with nps = 2 must be structurally unreachable (dimension
  routing precedes every kernel call); T1b falsifier decks probe it.
- **Recorder/query semantics** in 2D (force vs force-per-thickness) are
  documented in the T4 user guide — the silent-unit trap is the 2D analogue
  of the fork's "plausible but wrong by a factor" incidents.
- **Numerical:** `τ_seg` (physical zero-length refusal, relative to the
  surface's min initial segment length) and `τ_ov` (parametric overlap sliver,
  relative to `min(L_s,L_m)/L_m`) are **distinct named thresholds** (rev 1
  conflated them as one symbol); defaults set by the T1a/T3 oracles.
- **Backwards compatibility:** pure addition — no shipped 3D deck changes
  meaning; the 9-arg `contactPlane` form keeps working.

## Review record

**2026-08-18 — adversarial panel (3 reviewers: source-evidence, mechanics,
plan/process; all returned FIX-FIRST).** Sustained findings applied in rev 2:

1. Guard contract inverted: `ndf ≥ ndm` → **`ndf == ndm` exactly** on
   NTS/mortar (six shipped `!= 3` sites; `FE_Element::setID` packing).
   *(source lens; independently re-verified in-session)*
2. Orientation-vote "reuse" actually cited the opt-in `-smoothNormal` block;
   the default per-pair datum is ridge-flip-hazardous and silently degenerate
   on flush interfaces → §How/2 rewritten as a NEW interface-level vote with
   named aborts. *(source + mechanics lenses, convergent)*
3. Mortar→NTS fallback does not exist → §How/3 adopts the ADR-57 E7 ownership
   protocol minus the edge tier + `WARN_2D_PERP_NO_NTS`. *(mechanics blocker)*
4. 2D corner/vertex contact is a real hole, not measure-zero → vertex policy
   in §How/1 + G-T1a corner falsifiers; taxonomy claim softened. *(mechanics)*
5. Mortar integration measure unstated → slave-trace Jacobian
   `L_m/|t̂_s·t̂_m|` written down + tilted oracle case. *(mechanics)*
6. "Byte-identical battery" was not a checkable observable → `contact_dump.py`
   harness (T0) + N-unchanged clause; slice-twin gate restated with a measured
   formulation floor + contact observables; plane-stress gate added.
   *(plan blockers)*
7. `contactPlane` on `OPS_GetNDM()` rejected → surface-derived dimension;
   9-arg form permanently valid. *(plan)*
8. Flag conflation fixed: `-geomtan` (symmetric normal) vs `-consistanttan`
   (non-symmetric friction). *(plan; re-verified in-session)*
9. Staffing inverted where failure modes are silent vs loud: T1 split into
   T1a/T1b, **Fable moved to T1b**, T3 → Opus; escalation counters made
   mechanical; panel effort levels corrected; ledger/banner obligations moved
   into every phase PR. *(plan)*

Also corrected: SEGMENT/MORTAR `-visc` and `softKt` are ports, not inherited;
tie accumulators are 3-vectors reused with a zero slot; auto-kn needs a 2D
`normalOriented` sibling (its FATAL failure path is the first thing a 2D deck
hits); `segCoords` stays stride-3 z-padded for bucket-sort reuse; tolerance
semantics made relative; τ thresholds split; T2 parity oracle demoted from
referee to regression gate; thickness list completed (clamps, μ_c, ε_T) with
λ_N explicitly not double-scaled; Hertz P′/E\* pinned.

## Implementation log

- **T0 shipped** (PR #749, 2026-08-18): declarations + guards + rigid-plane
  acceptance + `contact_dump.py`. Gates: dump bit-identical pre/post (twice),
  battery 121/121, T0 suite 13/13 (8-row refusal matrix + control + 5
  rigid-plane tests). 8-angle review: 10 findings — 7 fixed in-PR, 1
  disclosed over-refusal, 2 deferred to T1b by name. Review also caught (and
  the test matrix now locks) the cross-surface dimension hole in the first
  C++ cut.
- **T1a shipped** (2026-08-18): `LadrunoContact2DKernel.h` (header-only, no
  TU includes it ⇒ build byte-identical by construction) +
  `proto_t1_nts2d.py` (95/95) + `t1a_cpp_check.cpp` (44/44); 190600-record
  C++↔numpy parity **bit-for-bit** on every continuous field; gap validated
  against 60-digit arithmetic (≤ 8·eps·|x|). Design corrections recorded in
  §How/1 (concave-only vertex pairs; ordered ownership; bisector side sign).
  Owed to T1b: add the header to `stamp_headers.py` GLOBS.
- **T1b panel deferrals (2026-08-18 — recorded per the FIX-FIRST verdicts;
  each with a phase owner):**
  - **D1 (T4):** G-T1b(f) auto-kn gates positivity, boundedness, `∝ E`, and
    equilibrium — NOT the absolute magnitude; the 2D↔brick-twin `kn = P/pen`
    comparison lands with the T4 benchmark battery.
  - **D2 (closed in T1b):** the segment-chaining input contract — enforced by
    the handle()-time scan (§How/1 item 6); T4 documents it in the user
    guide.
  - **D3 (T2):** the `-initial` (`addKiToTang`) 2D arm gains its first test
    in this PR (ModifiedNewton `-initial` patch rerun); the full
    implicit-transient 2D lane (`addCtoTang` under HHT/GN) is T2's gate.
  - **D4 (T2+):** the `NTS2D_END_SLACK = 1e-3` width is engineering margin
    (~3.9 orders above the measured tilt drift), not a measured optimum; the
    radial end-cap vertex pair replaces the window and its cliff.
  - **D5 (T2):** G-T2(f)'s removal gate must name the **vertex case**: a
    removed concave-vertex node whose adapter holds `prevFar`/`nextFar`
    pointers outside its connectivity. Also T2: the cross-type caveat — a
    reference-concave corner deformed convex keeps its committed side sign
    (the new side-flip warn latch makes it loud; collapse decks are the
    reachable case).
  - **D6 (T2):** the `-soft` pre-scan refusal-ordering cleanup, listed in
    rev 2's T1b bullet, was NOT done in T1b (with the lane live it is now
    correctly attributed in serial anyway); T2 owns it alongside the loop
    rebuild for friction.

- **T2 shipped** (2026-08-18): 2D NTS friction live -- the scalar return map
  (`returnMap1D`/`tangentBlock1D` appended to `LadrunoFrictionKernel.h`, a pure
  append so the 3D kernel is byte-identical by construction), the `-mu` FATAL
  retired into the live lane, `-consistanttan` wired, `softKt`'s 2D branch
  reached for the first time, and D3/D5/D6 discharged (D6 partially -- see
  below). Design note: `_adr85_t2_design.md`.
  - **The pre-T2 experiment (How/4) settled REUSE-vs-SCALAR on measurement, not
    taste**: evaluating the shipped 3D kernel degenerately takes the WRONG
    stick/slip branch in **89.2%** of disagreements adjudicated in exact
    rational arithmetic (vs 10.8% scalar), and -- the dominant finding -- a
    GLOBAL-FRAME VECTOR slip store loses the slip increment to catastrophic
    cancellation (median **1.5e-2**, max 100% relative error at
    `|ds|/|s| ~ 1e-15`, the regime a converging Newton lives in) where a
    tangential SCALAR store is exactly zero-error. This is the ADR-41 C3
    displacement-not-position crux resurfacing in a new frame; it is now a
    LEDGER_quirk. Artifacts:
    `contact_prototypes/t2_reuse_vs_scalar_check.cpp` (includes and calls the
    SHIPPED header) + `proto_t2_{reuse_vs_scalar,slip_cancellation}.py`.
  - **Gates MEASURED** at stamp `bce163ccf`: `contact_dump`
    `da73b6f8...7782` bit-identical to the pre-change baseline (184 lines, 0
    decks failed, captured twice); 3D battery **142 passed, N unchanged**; 2D
    suite **38 -> 54**; **G-T2(a) exact** -- last-STICK r = 0.99999999999999989,
    first-SLIP r = 1, a **0.5 ulp** bracket around T/N == mu.
  - **Three ADR corrections found during T2 and owed back to this document:**
    1. **The 3D battery is 142, not 121.** No subset of the adr39/41/57 globs
       (70/50/22) sums to 121, and no adr39/41/57 test changed since T0
       (`git diff 87a882464..HEAD -- tests/`), so 142 was the T0 count too and
       "121" is a stale/unreconstructible figure. There is also **no
       `tests/test_adr78_*.py`** at all, though the gate language names an
       "adr78 glob". T2 pins the glob explicitly and records N = 142.
    2. **There was no `-mu` refusal row to invert.** The phase plan describes
       "the retired `-mu` row" inverting into an acceptance test, but no test
       in `tests/test_adr85_contact2d_*.py` ever asserted the `-mu` 2D FATAL.
       T2 therefore ADDS a friction acceptance row (on the `pair_2d_now_live`
       precedent, asserting force TRANSFER); nothing was retired.
    3. **D6 is only partially discharged.** The ordering fix is scoped to the
       Contact/NTS `-soft` prescan (the lane the T0 note named); the
       RigidPlane and Mortar `-soft` mass-check loops were not re-audited and
       keep the original mis-attributing order.
  - **NOT gated (carried forward):** the D5 cross-type caveat -- a
    reference-concave corner deformed convex keeping its committed side sign.
    The `WARN_VTX2D_SIDE_FLIP` latch is wired
    (`LadrunoContactFE.cpp:1032`, shipped T1b) but no rig reproduced it across
    four attempts (kinematic drive -> singular; tethered -> slave leaves the
    wedge; untethered -> coincidence-floor guard / unstable). It is
    pre-existing T1b code, so this is inherited coverage debt, not a T2
    regression. **T2+ follow-up**, alongside D4's radial end-cap vertex pair.

*(per-phase design notes go to `_adr85_t*_design.md` if a phase needs one)*
