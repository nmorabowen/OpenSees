---
title: Ladruno Contact — perpendicular edge-edge contact (the cos_t→0 degeneracy)
project: Ladruno
status: draft (design-only — implementation is a later session)
priority: medium
owner: nmora
tags:
  - implementation
  - contact
  - edge-edge
  - narrow-phase
  - design-adr
---

# Ladruno Contact — perpendicular edge-edge contact (ADR-57)

> **Design ADR, not code.** This scopes the new narrow-phase algorithm for the
> `cos_t→0` case the shipped face-mortar lane (B2) documented as deferred. It is
> subordinate to the [[48_ladruno_contact_capstone_adr]] capstone (architecture +
> contracts + status-of-record) and graduates deferral row **#10 / the B2 follow-up
> note** out of [[47_ladruno_contact_deferrals_adr]]. Implementation is a separate,
> later session; this document fixes the kinematics, the substrate seam, the
> oracle-first phase gates, and the scope fence.

## What

Two faceted bodies can touch **edge-on** — a slab corner resting on a beam's top
edge, two cross-stacked bars, an L- or T-junction of walls, a fragment edge raking
across another fragment's edge after element removal. In that configuration the two
contacting facets are mutually **perpendicular** (their normals are ~orthogonal), so:

- the **face-mortar clip degenerates** — projecting the slave facet onto the master
  facet's plane yields a line, not an area, so there is **no in-plane overlap to
  integrate** (`A_clip → 0`); the shipped `LadrunoMortarKernel` returns the empty/
  degenerate-overlap status (`PairResult::status = −1`), contributing zero force; and
- **node-to-segment (NTS) misses it** — neither edge has a slave *node* that projects
  in-bounds onto the other facet's interior; the contact is between the *interiors of
  two line features*, not a point-on-face.

This ADR designs the missing primitive: a **segment-to-segment (edge-edge) penalty
contact** that activates exactly where face-mortar and NTS both go blind. The contact
between two non-parallel edges is a **single closest-point pair** (point-like, like
NTS — not distributed, like mortar), with the **common-perpendicular normal**
`n = (ê_s × ê_m)/‖ê_s × ê_m‖`.

**In scope (this ADR):** the governing kinematics (segment-segment closest point, the
edge normal + edge gap, first variation / B-operator); the **routing detector** that
sends a candidate facet-pair to edge-edge vs face-mortar vs NTS (and refuses the
*parallel*-edge degeneracy); a new OpenSees-free header kernel
`LadrunoEdgeKernel.h`; the seam into the existing substrate
(`LadrunoContactBucketSort` broad phase → `LadrunoContactDomain` path state →
`LadrunoContactFE` adapter, **new `EDGE_EDGE` mode**); penalty enforcement (MVP) with
friction (reusing `LadrunoFrictionKernel`), the explicit Courant-stable SOFT variant,
and an optional commit-cycle ALM `λ_N`; per-phase oracle-first gates.

**Not in scope (fenced to ADR-47 / ADR-55):** self-contact edge-edge (a surface vs
itself); **runtime body-discovery** edge-edge (debris/collapse — the re-emitting,
arbitrary-pair broad phase of [[55_ladruno_contact_runtime_discovery_adr]]); edge
**smoothing** (rounded/blended edges, Munjiza distributed-potential continuity across a
feature transition); and the consistent geometric `∂n/∂u, ∂s/∂u` tangent beyond the
gated MVP. The scope fence (§Fence) draws the lines precisely.

## Why

- **B2 named it as the one remaining gap in the explicit-contact narrow phase.** The
  SOFT=2 segment-based penalty closed the corner/edge/T-intersection cases NTS misses
  *when the facets still overlap in-plane* (cos_t away from 0). The capstone B2 row and
  the handoff both record the residual: *"the perpendicular edge-edge case (cos_t→0 ⇒
  the clip degenerates) needs a dedicated edge-edge treatment (deferred)."* This ADR is
  that treatment. With it, the narrow-phase taxonomy is **complete** for faceted
  bodies: vertex-face (NTS), face-face (mortar), **edge-edge (here)** — the three
  irreducible polyhedral contact cases (Wriggers 2006 §4; the skill's
  `broad_phase_collision.md §5`).
- **Edge-edge is genuinely irreducible.** It is the classic textbook gap: an edge-edge
  contact is **not reducible to vertex-face** (Ericson 2005 §5.1; SAT's edge×edge
  separating axes exist precisely because no face normal captures it). A model that
  only does vertex-face + face-face silently transmits **zero force** through an
  edge-on contact — a *correctness* hole, not an accuracy one (a slab corner would pass
  straight through a beam edge).
- **It is the enabling primitive for the collapse direction.** Once bodies fragment
  (ADR-51 element removal → ADR-54 FDEM-lite), arbitrary edges rake across arbitrary
  edges. The *predefined-surface* edge-edge designed here is the narrow-phase building
  block that ADR-55's runtime discovery will later feed with discovered edge pairs — so
  this ADR deliberately keeps the narrow phase **independent of how the pair was
  found** (a stateless adapter over 4 nodes), exactly so ADR-55 can reuse it.
- **The substrate is ready.** Every piece this needs already exists in a reusable form:
  the mortar lane's facet pairing already enumerates near facet-pairs (brute force — §7);
  `LadrunoContactProjection` already has the vec3 primitives; `LadrunoFrictionKernel`
  already is the one return map; the
  `LadrunoContactFE` Mode enum already has the extension point; `LadrunoContactDomain`
  already has the committed-state + GC + ALM patterns to copy. This is **net-new
  geometry on a proven seam**, not a new subsystem.

## Where

- **New code:**
  - `SRC/domain/contact/LadrunoEdgeKernel.h` — header-only, OpenSees-free
    segment-segment closest-point + edge normal/gap + the B-operator (mirrors the
    `LadrunoMortarKernel.h` / `LadrunoContactProjection.h` discipline; numpy-oracle
    testable build-free). The **only genuinely new geometry**.
  - `Ladruno_implementation/contact_prototypes/proto_e*_edge_edge.py` — the numpy
    oracles (one per phase gate, the fork's oracle-first standard).
  - `tests/test_adr57_edge_edge_*.py` — the OpenSeesPy battery additions.
- **Modify (extend, do not duplicate):**
  - `SRC/domain/contact/LadrunoContactDomain.{h,cpp}` — a new `EdgeEdgeState`
    committed-state struct + `getOrCreateEdgeEdgeState` + GC (`edgeGC{Begin,Mark,End}`)
    + extend `commit()`/`revertToLastCommit()` to iterate it (capstone contract #1, the
    same explicit obligation the mortar lane had); a new `EdgeContact` definition struct
    + `addEdgeContact` (or fold a `-edgeedge` flag onto `MortarContact` — see §How/API);
    optional `λ_N` ALM per edge-pair (reuse the `MortarNormalState` Uzawa-on-commit
    pattern).
  - `SRC/analysis/handler/LadrunoContactFE.{h,cpp}` — **new `Mode EDGE_EDGE = 4`** with
    `getResidual`/`getTangent`/`addCtoTang` for the 4-node edge-pair; stateless view.
  - `SRC/analysis/handler/LadrunoContactHandler.{h,cpp}` — the pairing loop: for each
    near facet-pair the **mortar brute-force enumeration** returns (§7 — *not* the bucket
    sort, which is NTS point-vs-segment only), run the **routing detector**; on an
    edge-edge route, enumerate the candidate edge pairs and inject one `EDGE_EDGE`
    adapter per active edge pair.
  - `SRC/classTags.h` — no new class tag needed if `EDGE_EDGE` is a *mode* of the
    existing `LadrunoContactFE` (preferred; mirrors RIGID_PLANE/SEGMENT/MORTAR sharing
    one runtime tag). Record the decision in `LEDGER_implementations.md`.
  - Command surface in the interpreter (`OPS_LadrunoContact*`), banner feature line,
    and the three ledgers — per the fork's build-control rules.
- **Reference (copy patterns from):**
  - `LadrunoMortarKernel.h` — the clip→Gauss header discipline + the `PairResult`
    return-struct idiom + the convexity/degeneracy refusal pattern.
  - `LadrunoContactFE.cpp` SEGMENT/MORTAR modes — the stateless-adapter
    residual/tangent assembly + the B1/B2 explicit-SOFT fast path.
  - `LadrunoContactDomain` `FrictionState` + `MortarNormalState` — the committed-state /
    GC / Uzawa-on-commit templates to mirror for `EdgeEdgeState`.
- **Build:** no new CMake target or external dep — header-only kernel + existing
  contact translation units. No `Ladruno_internal/01_compilation_journal` entry needed.

## How — governing kinematics

### 1. Segment-segment closest point (the new primitive)

Two edges, each a parametric segment over its two end-nodes:

```
slave  edge  A(s) = a₀ + s·ê_s·L_s,   ê_s = (a₁−a₀)/‖a₁−a₀‖,   s ∈ [0,1]
master edge  B(t) = b₀ + t·ê_m·L_m,   ê_m = (b₁−b₀)/‖b₁−b₀‖,   t ∈ [0,1]
```

The closest-point pair `(s*, t*)` minimizes `f(s,t) = ‖A(s) − B(t)‖²`. With
`d₀ = a₀ − b₀`, `r₁ = a₁−a₀`, `r₂ = b₁−b₀`, the stationarity `∂f/∂s = ∂f/∂t = 0` is the
2×2 linear system (Ericson 2005, *Real-Time Collision Detection* §5.1.9
`ClosestPtSegmentSegment`; Wriggers 2006 §4.3 beam/edge contact):

```
[ r₁·r₁   −r₁·r₂ ] [ s ]   [ −d₀·r₁ ]
[ r₁·r₂   −r₂·r₂ ] [ t ] = [ −d₀·r₂ ]

denom = (r₁·r₁)(r₂·r₂) − (r₁·r₂)²  =  ‖r₁‖²‖r₂‖²·sin²θ_edge
```

`denom` is the **dual-purpose degeneracy gauge** (Lens-A finding): `denom → 0` ⇔
`sinθ_edge → 0` ⇔ the edges are collinear/parallel, so (i) the closest-point pair is
**not unique** (a whole interval) **and** (ii) the cross-product normal `n = ê_s×ê_m`
(§2) is **ill-conditioned** — normalizing `n = m/‖m‖` divides by `sinθ_edge`, and
`∂n/∂u ∝ 1/sin²θ_edge`, so the normal *direction* blows up well before the closest
point becomes non-unique. The threshold `τ_∥` therefore governs **both** uniqueness and
normal conditioning, and is set by the conditioning bound (the stricter of the two): at
`θ_edge = 5°`, `1/sin²θ_edge ≈ 130×` amplification of nodal perturbation into `n` —
arguably already too permissive, so the kernel uses a **conditioning-justified**
`τ_∥` (tracked as an open E0-gated number, default `sin²(15°) ≈ 0.067`, **not** 5°) and
records the amplification at the chosen angle. Wriggers 2006 §4.3 uses the same
`ê_s×ê_m` beam-to-beam normal and notes it degenerates as the edges align.

The kernel branches:

- **`denom ≥ τ_∥·‖r₁‖²‖r₂‖²`** (non-parallel, well-conditioned normal): solve the
  unclamped line-line `(s,t)`, then recover the segment-segment closest point with
  **Ericson's exact branch ladder** (RTCD §5.1.9 `ClosestPtSegmentSegment`, verbatim —
  *not* the single back-substitution the draft sketched, which gives the wrong vertex-
  vertex sub-region): clamp `s`; recompute `t = (d₀+s·r₁)·r₂/(r₂·r₂)`; **if** `t ∉
  [0,1]`, clamp `t` and recompute `s = (t·r₂ − d₀)·r₁/(r₁·r₁)`; clamp `s`. This
  conditional ladder is exact in all 8 Voronoi sub-regions (interior, 4 edge, 4
  vertex); a single unconditional pass is not. (E0 must include the vertex-vertex
  double-clamp falsifier.)
- **`denom < τ_∥·…`** (parallel / ill-conditioned): **refuse** — return `status =
  PARALLEL`. Parallel edges are *not* an edge-edge contact; routing (§2) never sends a
  parallel pair here, but — per Lens-C — the refusal must **cross-check the clip area
  before dropping the pair** (a near-parallel, edge-on X-crossing has no clip area
  *and* is refused here *and* has no in-bounds NTS node: a triple-miss lost-contact
  hole). The handler treats a `PARALLEL` return in a sub-`τ_A` clip neighborhood as a
  routing escalation (E1 falsifier), not a silent zero.
- **Zero-length / coincident-node edge** (`‖a₁−a₀‖ < τ_len·L_ref` or likewise master —
  a real case once elements are removed in the collapse motivation): return `status =
  DEGENERATE`; the `ê` normalization is undefined. (E0 zero-length falsifier.)

Outputs: closest points `P_s = A(s*)`, `P_m = B(t*)`, and the **margin'd interior flag**
`s*,t* ∈ [δ, 1−δ]` (δ ∝ edge-length scale). Only a **strictly margin-interior** pair is
a true edge-edge crossing where the cross-product normal is the line connecting the
closest points (§2); a pair at/near an endpoint (`s*≈1`) is a near-vertex contact where
`g_N = w·n ≠ ±‖w‖` — routing **demotes** it to NTS/vertex, and an E0 falsifier
(`|w·n| == ‖w‖` to tolerance across the admitted band) enforces the demotion. The open
interval `(0,1)` of the draft was insufficient — `s*=0.999` is numerically a vertex.

### 2. Edge normal, edge gap, and the degenerate-cos_t detector

**Edge normal.** For two crossing edges the contact normal is the **common
perpendicular** — the unique direction orthogonal to both edge tangents:

```
m = ê_s × ê_m,    n = m / ‖m‖           (‖m‖ = sinθ_edge — defined ⟺ non-parallel)
```

**Sign — body-fixed, committed (NOT position-derived).** The draft's
`n ← sign((P_s−P_m)·n)·n` is **self-referential and masks interpenetration** (Lens-A
BLOCKER-class catch): `g_N = w·n` is the very quantity whose sign encodes penetration,
so as bodies pass *through* contact, `w` reverses, the rule flips `n`, and the
recomputed `g_N` flips back to positive — the contact reads **open while
interpenetrating**. Unlike NTS (whose master *facet* normal is a body-fixed geometric
orientation independent of the slave position), the edge cross-product normal has no
intrinsic anchor, so a position-derived sign is genuinely unstable through `g_N=0`.
**Fix:** orient `n` **once per persistent pair** from a body-fixed reference available
at routing — the two facet normals `{n̂_s, n̂_m}` (or the master facet outward) — and
**hold the sign committed in `EdgeEdgeState`** (§5), re-evaluating only on re-pairing,
exactly as friction history is carried. The optional user `outward` is the tie-breaker
only at first capture. An E2 falsifier drives a pair from `g_N>0` through `0` into
penetration and asserts `g_N` **decreases monotonically** (no flip).

**Edge gap (signed scalar).** With `w = P_s − P_m` (separation vector at the closest
points) and the committed sign:

```
g_N = w · n           (penetration ⇔ g_N < 0, the fork's convention)
```

For the **strictly margin-interior** closest point of two straight segments, `w ⊥ ê_s`
and `w ⊥ ê_m` **exactly** (not "to first order"), so `w ∥ n` exactly and `g_N = ±‖w‖` —
the minimum segment-segment distance, signed. This identity holds **only** in the
margin-interior set (§1); at a clamped/near-vertex pair `g_N = w·n` is the projection,
not the magnitude, which is precisely why §1 demotes those pairs rather than firing
edge-edge on them.

**The routing detector.** A candidate (slave facet, master facet) pair is pre-classified
by a single cheap scalar, the **facet-normal alignment** — named `align_fn` here to
**avoid colliding with the kernel's existing `cos_t`** (Lens-C catch): in the shipped
`LadrunoMortarKernel`, `cos_t = |n_s·n0|` is the area-Jacobian denominator
`A_phys = A_aux/cos_t`, guarded by a hard `cos_t < 1e-12 ⇒ return −1` (see
`LadrunoMortarKernel.h:301`). The router scalar has the **same operands** but a
**different role** (a 50°/70° classifier band, not a 1e-12 divide guard):

```
align_fn  =  | n̂_s · n̂_m |        ∈ [0,1]   (facet normals; same operands as the kernel's cos_t)
```

| `align_fn` regime | geometry | candidate owner |
|---|---|---|
| `→ 1` | facets **parallel / facing** | face-mortar (clip has area) or NTS |
| `→ 0` | facets **perpendicular / edge-on** | **edge-edge (this ADR)** — clip has no area |
| in between | oblique | resolved by **ordered precedence** below |

> **Resolving the B2 "cos_t→0" shorthand.** The capstone/handoff wrote "cos_t→0 ⇒ the
> clip degenerates." Verified against the shipped oracles, that `cos_t` is the
> facet-normal alignment `|n̂_s·n̂_m|` — but it conflates **two** distinct in-tree
> failure points as `cos_t→0`: the **zero-area clip** (a perpendicular slave facet
> projects to a line ⇒ `signed_area < area_tol` ⇒ empty) **and** the **Jacobian
> blow-up** `A_aux/cos_t`. Either way the mortar pair contributes nothing — that is the
> degeneracy edge-edge fills. `align_fn` must **not** be confused with
> `sinθ_edge = ‖ê_s×ê_m‖` (the *edge-direction* gauge of §1): edge-edge wants the
> *facets* perpendicular (`align_fn→0`) **and** the *edges* non-parallel
> (`sinθ_edge` away from 0) — opposite roles, distinct symbols.

**Ordered precedence with a guaranteed owner (NOT a clean partition).** The draft
claimed "exactly one lane, no superposed force — a partition." That over-claims:
(a) the clip-area arbiter only adjudicates edge-edge vs face-mortar — **NTS is a third
lane** the partition ignored, and a slave node near a master edge can fire NTS *while*
the boundary edges fire edge-edge, double-counting; (b) the ADR's own oblique-band
question concedes a pair can present *both* a small clip and a valid edge-edge pair.
So the detector is specified as **ordered precedence with hysteresis and a
guaranteed-owner fall-through**, not a partition, and the residual double-count is
*bounded and gated*, not *zero by construction*:

0. **Proximity gate.** Facets farther apart than `d_band` own **no** lane (`NONE`) —
   genuinely separated, no contact primitive. Within `d_band` an owner is **guaranteed**.
1. **Regime split on `align_fn` (refined oracle-first in E1 — the flat precedence of the
   draft was wrong).** The E1 oracle showed that a flat "NTS ≻ edge-edge ≻ mortar"
   *steals facing-facet contact from the mortar lane*: a facing overlap has slave nodes
   projecting in-bounds (the NTS condition), so flat NTS-precedence would mis-route the
   whole mortar regime to NTS. Fix: **`align_fn ≥ τ_face` (≈cos 50°) ⇒ FACE-MORTAR owns
   the pair outright** — its weighted-gap integral already subsumes the in-bounds slave
   nodes; a separate NTS there would double-count. NTS-vs-edge precedence applies **only
   in the non-facing regime**.
2. **Non-facing regime (`align_fn < τ_face`): ordered precedence NTS ≻ edge-edge ≻
   sliver-mortar.**
   - **NTS** if a slave **vertex** projects in-bounds on the master face with `g_N ∈
     [−d_pen, +d_band]` — the more local primitive owns it, so edge-edge cannot
     double-count a vertex-on-face (the Lens-C third-lane fix).
   - else **EDGE-EDGE** if ≥1 of the (≤4)×(≤4) edge pairs is a *margin-interior,
     well-conditioned* (E0-`EE_OK`) crossing with `g_N ≤ d_band`.
   - else **FACE-MORTAR** — the **guaranteed owner** for everything else in proximity
     (even a sliver clip). Face-mortar stands down *only if* edge-edge succeeds, so
     `A_clip < τ_A` is necessary-but-not-sufficient to vacate it: **no coverage hole**.

E1 gates the lot (oracle **16/16**): the regime split, the NTS-vs-edge double-count
falsifier, the oblique-band no-coverage-hole falsifier, the parallel-but-offset facing
case (caught by mortar, not dropped), and the hysteresis no-chatter sweep (a monotone
tilt through 45° produced a single EDGE_EDGE→FACE_MORTAR transition). `align_fn` is the
cheap regime pre-classify, **not** the arbiter — correctness comes from the
proximity-gated precedence. If precedence ever chatters in the band, a partition-of-force
*blend* is the documented fallback; the MVP ships hard precedence with the E1 falsifiers.

### 3. First variation and the B-operator (penalty force + main tangent)

The 4-node connectivity is `[a₀, a₁, b₀, b₁]` (12 DOFs in 3D). Holding the closest-point
parameters and normal frozen (the standard NTS/mortar first-order simplification; the
curvature terms are the gated consistent-tangent refinement, §How/E4):

```
δg_N = n·(δP_s − δP_m),   δP_s = (1−s)δa₀ + s δa₁,   δP_m = (1−t)δb₀ + t δb₁
⇒  B = [ (1−s)·nᵀ | s·nᵀ | −(1−t)·nᵀ | −t·nᵀ ]      (1×12)
```

This is the **NTS B-operator with the master shape functions `Nᵢ` replaced by the edge
linear weights `(1−t, t)`** and the slave node replaced by the slave-edge weights
`(1−s, s)` — a clean reuse of the existing assembly shape.

- **Normal penalty force:** `t_N = ε_N⟨−g_N⟩` (Macaulay; active only in penetration);
  residual `f = t_N·Bᵀ` distributes the concentrated edge force to the 4 nodes with the
  edge weights (self-equilibrating: `Σf = 0` since the slave and master weights each sum
  to 1 with opposite sign).
- **Main tangent:** `K = ε_N·BᵀB` (12×12, symmetric, rank-1, PSD) — solver-safe on any
  system, exactly as the NTS/mortar main blocks.

### 4. Friction (reuse `LadrunoFrictionKernel` verbatim)

The contact tangent plane is spanned by two unit tangents orthogonal to `n`. Since
`n = ê_s×ê_m` is orthogonal to `ê_s` **by construction**, the projection in the draft's
`ê_s⊥` is a no-op (Lens-A): the frame is simply `t₁ = ê_s`, `t₂ = n × ê_s`. Note
`‖t₂‖ = ‖(ê_s×ê_m)×ê_s‖ = sinθ_edge`, so the **tangent frame inherits the same
`1/sinθ_edge` conditioning** as the normal (§1) — the `τ_∥` gate covers it too. The
tangential slip
increment is the relative motion of the two closest points projected onto `{t₁,t₂}`,
**built from displacement increments, not positions** (the C3.1 / NTS-SEGMENT lesson:
the closest-point construction makes the relative *position* purely normal — a
[[LEDGER_quirks]] trap already documented). Feed `Δg_T`, the committed elastic slip
`gpT`, and the engagement origin `gT0` into the shipped `frictionReturnMap` with the
unified cone `cap = min(μ|t_N| + c, τ_max)`; assemble `tFric` via `Bᵀ` per tangent and
the symmetric `frictionTangentBlock` (`Csl` non-symmetric Coulomb branch behind
`-consistanttan`, like everywhere else). **No new friction code** — edge-edge is just a
fourth consumer of the one return map (NTS, mortar, SOFT=2, now edge-edge).

### 5. Path state — `EdgeEdgeState` (Domain-owned, committed-only)

Keyed by **(contactTag, slaveEdgeId, masterEdgeId)** — a 3-int key reusing the existing
`PairKey` machinery. The edge ids must be **rebuild-stable** (the slot carries committed
friction + λ + the sign anchor across the adapter rebuild every `handle()`); to avoid a
lossy hash collision silently merging two edge states (Lens-B), use a **deterministic
global edge ordinal** (matching the existing `segIndex`/`slaveFacetIndex` discipline) or
the **ordered node-tag pair as a composite key** — never a lossy hash:

```
struct EdgeEdgeState {
    // friction (mirrors FrictionState exactly)
    double gpT[3], gpTtrial[3], gT0[3];   bool engaged;
    double gT0committed[3];               bool engagedCommitted;  // implicit-revert double-buffer
    // optional ALM (mirrors MortarNormalState.lambdaN — ONE scalar per edge pair, point-like)
    double lambdaN;                       // committed normal multiplier (≤0); updated ONLY in commit()
    double gN_committed;                  // the committed gap for the Uzawa update + query
    // committed normal SIGN anchor (§2 Lens-A fix): the body-fixed orientation captured at
    // first engagement, so g_N cannot flip through contact. +1/−1; 0 = not yet captured.
    int    signN;
};
```

The `signN` field is the §2 fix: the normal's orientation is captured **once** from the
facet-normal reference at first engagement and committed; the adapter applies it every
`getResidual` rather than re-deriving the sign from `w·n` — this is what makes `g_N`
monotone through contact (the E2 falsifier).

Lifecycle obligations (capstone contract #1, #2): `commit()` promotes `gpT←gpTtrial`,
`gT0committed←gT0`, and (if ALM) `lambdaN ← min(0, lambdaN + ε_N·g_N_committed)` once per
`Domain::commit()`; `revertToLastCommit()` restores the trials from committed. A new GC
trio (`edgeGCBegin/Mark/End`) drops dead slots each `handle()` (the engine survives
`domainChanged` — the `theEQs` leak class). **Extend the existing `commit`/`revert`
loops to iterate `theEdgeEdgeStates`** — the shipped hooks iterate only friction +
mortar slots, so this is an *explicit obligation on the edge-edge PR*, not inherited
(the same trap the mortar lane had, capstone contract #1).

### 6. Enforcement tiers (mirror the shipped lanes)

- **Penalty (MVP):** `t_N = ε_N⟨−g_N⟩`, `ε_N` from `-kn`/`-kn auto` (sized from the
  owning solid element, the existing `knAuto` path). Single-pass, implicit + explicit.
- **Explicit Courant-stable SOFT:** under `CentralDifferenceLadruno`, replace `ε_N` with
  `k_soft = SOFSCL·4·m_eff/dt²`, `m_eff = 1/(B M⁻¹ Bᵀ)` the **gap-mode generalized mass
  of the 4-node edge operator** using the assembled nodal-mass cache (the B1/B2
  `ladrunoBuildNodalMass` cache, reused verbatim — `B = [(1−s)n, s n, −(1−t)n, −t n]`).
  So the edge contact's own `ω·dt = 2√SOFSCL ≤ 2` and never throttles `dt_cr` — the
  progressive-collapse enabler, now for edge-on impact. Friction `softKt` follows the
  B1-kt n→t rule. Explicit-only (`dynamic_cast` to CDL); implicit / SOFT-absent
  byte-identical.
- **Commit-cycle ALM (optional `λ_N`):** one scalar Uzawa per edge pair (point-like ⇒
  *much* simpler than the mortar per-node global accumulator — no shared-node
  variational-consistency problem, because an edge pair's `g_N` is one number). Drives
  penetration → an `ε_N`-independent tol; the held-load `analyze_augmented` proc (D1)
  works unchanged. Implicit-only (the mass-only explicit tangent degenerates Uzawa to
  single-pass penalty — the disclosed ADR-47 limitation, inherited).

### 7. Broad-phase seam — inherit the mortar lane's facet pairing, add NO spatial structure

**Correction (Lens-B catch): the mortar lane does NOT use the bucket sort.** The shipped
mortar pairing (`LadrunoContactHandler.cpp:521-528`) is an explicit `O(nSegS·nSegM)`
**brute force** (its own comment defers a "slave-aware mortar broad phase / mortar
P2.5"); `LadrunoContactBucketSort` is a **point-vs-master-segment** query keyed on a
single slave-node coordinate, used by the NTS lane only, and its superset contract is
proven for *point* slaves. So edge-edge **cannot** "reuse the bucket sort to iterate
facet-pairs" — there is no such facility. Edge-edge instead **inherits the mortar lane's
brute-force facet enumeration** (over which it is a trivial superset):

1. Iterate the mortar lane's near (slave facet, master facet) pairs as today (brute
   force; or the deferred mortar-P2.5 broad phase if/when that lands — an *optional
   optimization*, not an edge-edge prerequisite).
2. For each facet-pair, run the **routing detector** (§2). On an edge-edge route,
   enumerate the `(≤4)×(≤4)` boundary-edge pairs, run the `LadrunoEdgeKernel`
   closest-point on each, **keep** the margin-interior + well-conditioned + in-band
   pairs, and inject one `EDGE_EDGE` adapter per kept pair.
3. **De-dup shared edges** (a facet edge is shared by two facets): key the injected
   adapter by the **ordered edge-node-tag pair** in a handler-side `std::set` — new
   handler code that *mirrors* the bucket sort's sparse-set stamp idiom, not a reuse of
   that structure (which is internal to `Grid`).

The cost is **zero new broad-phase memory** (the enumerator is the existing brute-force
loop) and the superset property is trivial (every near facet-pair is visited). The
reference-config scoping limit is inherited: the facet geometry is built at `handle()`
from reference coords, so **large sliding** needs re-emission — the same epoch-re-emit
follow-on the subsystem shares, and the gateway to ADR-55's runtime discovery.

**Required small plumbing (Lens-B):** the `align_fn` pre-classify needs the two facet
normals, which `PairResult` does **not** currently expose (`LadrunoMortarKernel.h:82-89`
returns status/area/D/M/g only; `facetNormal()` exists but is internal). The edge-edge
PR adds a cheap facet-normal computation at the router and/or surfaces `align_fn` from
the kernel — a few lines, flagged here so it is not mistaken for a free read. The
**authoritative** arbiter (clip-area degeneracy via `status==−1` / `area≈0`) *is*
directly observable today.

### Public API (proposed)

The least-friction surface: edge-edge is an **opt-in modifier on a mortar contact**
(the surfaces are already declared as faceted master/slave), gated **off by default** so
the null build is byte-identical:

```tcl
constraints LadrunoContact
contactSurface 1 -kind masterSegments -faces $master
contactSurface 2 -kind slaveSegments  -faces $slave

# face-mortar with the edge-edge fallback enabled (cos_t→0 pairs routed to edge-edge):
contact 11 -master 1 -slave 2 -mortar -enforce alm -epsN auto -augTol 1e-8 \
        -edgeedge [-edgeKn auto] [-edgeMu 0.3] [-edgeBand <d>]
#   -edgeedge      : enable the cos_t→0 edge-edge fallback (DEFAULT OFF ⇒ byte-identical)
#   -edgeKn        : the edge-edge penalty (default = the mortar epsN)
#   -edgeMu        : edge-edge friction (default = the mortar mu)
#   -edgeBand      : the gap activation band d_band (default from -kn auto length scale)
```

`-soft` (explicit SOFT), `-visc` (viscous), and `-consistanttan` compose with `-edgeedge`
exactly as on the mortar lane. A standalone `-formulation edge-edge` (no face-mortar
companion) is a later convenience, not the MVP.

## Per-phase gates (oracle-first — numpy proto before C++, the fork discipline)

| Phase | Delivers | Oracle (`proto_e*`) + falsifier | C++ gate |
|---|---|---|---|
| **E0** ✅ **SHIPPED** | `LadrunoEdgeKernel.h` — segment-segment closest point + edge normal/gap | `proto_e0_closest_point.py` **28/28**: vs an independent scan reference on all **8 Voronoi sub-regions** incl. the **vertex-vertex double-clamp**; **parallel refusal** (denom<τ_∥, τ_∥=sin²15° conditioning-justified) + near-parallel; **zero-length** refusal; conditioning bounded across the band; `|w·n|==‖w‖` only margin-interior (+ its failure clamped); FD `∂g_N/∂coord` == analytic B to **1.1e-9**; committed-sign monotone-through-contact + the buggy-`w·n`-rule penetration-mask demo (the A-4 fix, shown empirically). | `e0_cpp_check.cpp` **13/13** bit-for-bit (B==FD identical 1.13e-9); header compiles warning-free; header-only ⇒ no TU includes it yet ⇒ build byte-identical |
| **E1** ✅ **oracle SHIPPED** (C++ wiring folds into E2) | the **routing detector** | `proto_e1_router.py` **16/16**: the **align_fn regime split** (facing→mortar; non-facing→NTS≻edge≻sliver — oracle-first caught that flat NTS-precedence steals facing-facet contact), X-crossing→EDGE_EDGE, an **NTS-vs-edge double-count falsifier** (non-facing vertex-poke→NTS), the **oblique-band no-coverage-hole falsifier** (in-proximity pair ⇒ non-NONE owner), the **parallel-but-offset facing** case caught by mortar (not dropped), the proximity NONE gate, and the hysteresis no-chatter sweep (1 clean transition through 45°). | handler unit (lands with the E2 wiring — routing without an adapter to inject is inert): edge-on→EDGE_EDGE, facing→MORTAR, slave-vertex-in-bounds→NTS; **byte-identity when `-edgeedge` absent** |
| **E2** ✅ **SHIPPED** | penalty normal force + main tangent | `proto_e2_penalty.py` **23/23**: two-bar edge-on test, force ⟂ both edges (along `n`), self-equilibrium `Σf=0`, `K=ε_N BᵀB` FD-checked symmetric PSD rank-1, penetration `δ=P/ε_N`; **monotone `g_N` through contact** driven through a real penalty step (gap>0 → 0 → penetration, no sign flip — the committed-sign-anchor falsifier + the buggy-`w·n`-rule penetration-mask demo). | `test_adr57_edge_edge_1` **3/3**: two crossed bars restrained at `δ≈P/ε_N`, Newton converges (the `K_c` check); **NTS-passes-through vs edge-edge-restrained** (the headline falsifier — same geometry NTS contact transmits ZERO force ⇒ the edge sinks ~10³× deeper); **the EE-1 oblique-band regression** (a 60° pair, non-degenerate clip, edge-edge owns ⇒ single-primitive δ≈P/ε_N — guards the face-mortar stand-down); byte-identity when no pair routes / `-edgeedge` absent. Full contact battery **145 passed** (no regression). 3-reviewer adversarial gate PASS (2 findings folded). |
| **E3** | friction (reuse `LadrunoFrictionKernel`) | `proto_e3_friction.py`: stick/slip on the edge tangent plane, `a = g(sinθ−μcosθ)` incline sign, Tresca cap, slip-from-displacement (the C3.1 trap), `μ=0` byte-identical. | `test_adr57_edge_edge_2`: explicit raking bar (Coulomb opposes motion), implicit stick converges; symmetric tangent solver-safe |
| **E4** | consistent geometric tangent (`-edgegeomtan`, **gated off**) | `proto_e4_geomtan.py`: `∂n/∂u, ∂s/∂u, ∂t/∂u` analytic == FD on a skew large-sliding pair; symmetric; flat ⇒ byte-identical. | local Newton iter-count improvement on a curved-rake case; ProfileSPD-safe |
| **E5** | explicit Courant-stable SOFT | `proto_e5_soft.py`: 4-node edge `m_eff = 1/(B M⁻¹ Bᵀ)` closed form; `ω·dt = 2√SOFSCL`; stiff diverges vs soft bounded; energy restitution. | `test_adr57_edge_edge_3`: edge-on impact at the STRUCTURAL dt; implicit byte-identical |
| **E6** | optional ALM `λ_N` (one scalar/pair) | `proto_e6_alm.py`: penetration → ε_N-independent tol; release → `λ_N→0`, F=0; eqn count constant across augmentations. | held-load `analyze_augmented` drives `‖g_N‖→augTol`; opt-in byte-identity |
| **E7** | integration + regression | — | full battery (a slab corner on a beam edge; cross-stacked bars; an L-junction), **byte-identical when no pair routes edge-edge**, 3-reviewer adversarial gate |

Every `proto_e*` is numpy-only (build-free); Zone-A CI is the C++ no-regression gate;
keep the capstone status-of-record row + the three ledgers current **in the same PR** as
each phase (the build-control rule).

## Fence — what this ADR is and is NOT

**vs the shipped face-mortar lane (the complement, not a competitor).** Face-mortar owns
`align_fn → 1` (facing facets, the clip has area); edge-edge owns `align_fn → 0`
(perpendicular facets, the clip degenerates). The routing detector is **ordered
precedence with a guaranteed owner** (§2) — NTS ≻ edge-edge/face-mortar, with the
clip-area test the edge-edge↔face-mortar arbiter and a fall-through so no band pair drops
to zero force; the residual double-count is bounded by `τ_A` and E1-gated (it is *not* a
clean partition — that over-claim was the Lens-C MAJOR). Edge-edge **reuses** mortar's
kernel (for the clip-area degeneracy signal), the projection's vec3 primitives, the
friction kernel, the SOFT mass cache, and
the ALM Uzawa pattern; it **duplicates none** of them.

**vs ADR-47 (the deferral ledger).** This ADR **promotes** the edge-edge follow-up out of
the B2 deferral note into a committed design. It does **not** pull in the rest of the
ADR-47 set: the **parallel-edge** degeneracy stays NTS/vertex/face-mortar (refused by the
kernel, §1); **edge smoothing** (rounded/blended edges, Munjiza distributed-potential
continuity across a corner-tumbles-off-an-edge transition) stays deferred under ADR-47
#4/#4a (slide-line smoothing); **dual-basis / true-LM** stay ADR-47 #1/#2 (edge-edge is
point-like penalty/Uzawa, no mortar basis); **anisotropic friction** stays #5.

**vs self-contact + runtime discovery (ADR-47 #3 / ADR-55) — KERNEL OWNERSHIP RESOLVED.**
This ADR designs **predefined-surface** edge-edge: a *declared* master facet-surface
against a *declared* slave facet-surface, enumerated by the *existing reference-config*
mortar facet pairing. **Self-contact** (a surface raking its own other edge) and
**runtime body-discovery** (debris fragments, arbitrary edge-vs-edge, no declared
pairing) are **out of scope** → ADR-55's re-emitting, self-pair-excluding, current-config
broad phase + self-contact's penalty-only over-constraint handling (ADR-47 #3).
**Critical fence fix (Lens-C):** ADR-55 was written the same day and its **P3 row +
subsystem-5 + title** independently claim to *build* "a new edge-edge narrow-phase
kernel" — a direct **kernel-ownership collision** with this ADR. Resolution: **ADR-57
owns the edge-edge *kinematics + kernel* (`LadrunoEdgeKernel` + the `EDGE_EDGE`
adapter); ADR-55's P3 is demoted to *consume* that kernel over runtime-discovered
pairs**, not re-build it. The adapter is deliberately kept **discovery-agnostic** (a
stateless 4-node view) so ADR-55 feeds it discovered pairs without re-designing the
kinematics. ADR-55 is edited reciprocally when this ADR lands (its P3 row + subsystem-5
bullet now point here).

**Re-open triggers for what this ADR still defers:** **edge-end handoff smoothing**
re-opens if the E7 sliding-rake gate shows force chatter as a contact point tumbles off
an edge *end* (the edge-edge→vertex transition) that the gap-band + `-visc` cannot damp.
*Citation correction (Lens-C):* the right reference is the **Munjiza distributed
potential** (a *new* ADR-47 row), **not** ADR-47 #4a — #4a is *averaged nodal-normal*
smoothing for *facet*-normal chatter across facet junctions, which does nothing for the
edge-end primitive handoff. Self-contact/runtime edge-edge re-opens with ADR-55 / the
FDEM-lite collapse line (ADR-54).

**vs dimensionality.** Edge-edge is a **3D-only** phenomenon — `n = ê_s×ê_m` requires
3D, and in 2D two non-parallel crossing segments *are* a vertex/segment NTS penetration,
not edge-edge. A 2D contact surface **never routes here**.

## Decision log (resolved design questions)

| Decision | Resolution | Rationale |
|---|---|---|
| Distributed vs point contact | **Point-like** (single closest-point pair, NTS-shaped B-operator) | Two non-parallel edges touch at one point; a distributed (mortar) treatment has no overlap measure to integrate — that *is* the degeneracy |
| Normal definition | **Common perpendicular** `n = (ê_s×ê_m)/‖ê_s×ê_m‖`; sign **body-fixed + committed** (`EdgeEdgeState.signN`), NOT position-derived | The unique direction ⟂ both edges; a `w·n`-derived sign masks interpenetration through contact (Lens-A) — anchor it once per pair |
| The two gauges | **facet-normal `align_fn`=\|n̂_s·n̂_m\|** (route to edge-edge as →0) vs **edge-direction `sinθ_edge`=‖ê_s×ê_m‖** (refuse as →0) — distinct symbols | `align_fn` renamed off the kernel's existing `cos_t` (an area-Jacobian guard, Lens-C); perpendicular *facets* + non-parallel *edges* |
| Routing | **Ordered precedence with a guaranteed owner** (NTS ≻ then align_fn-classified edge-edge/face-mortar), hysteresis band, residual double-count *bounded by τ_A & E1-gated* | NOT a clean partition — NTS is a third lane + the oblique band can present both; precedence + fall-through closes the coverage hole (Lens-C) |
| `τ_∥` parallel gate | Set by **normal conditioning** (`1/sin²θ_edge`), default `sin²(15°)`, E0-gated — NOT the draft's 5° | 5° → ~130× normal amplification; conditioning bites before non-uniqueness (Lens-A) |
| New kernel vs reuse projection | **New header `LadrunoEdgeKernel.h`** (segment-segment is new geometry) reusing the projection vec3 primitives; Ericson §5.1.9 exact clamp ladder | Point-surface ≠ line-line; but no duplicated vec math |
| State key | **(contactTag, slaveEdgeId, masterEdgeId)** in `PairKey`; ids = deterministic global ordinal OR ordered node-tag composite — **never a lossy hash** | Stable across `handle()` rebuilds; a hash could silently merge two edge states (Lens-B) |
| Enforcement | Penalty MVP + explicit-SOFT + **one-scalar ALM** (point-like) | Inherits all shipped tiers; ALM is *simpler* here (no shared-node accumulator) |
| API | **`-edgeedge` opt-in modifier on `-mortar`**, off by default | Surfaces already declared; null/off ⇒ byte-identical |
| Broad phase | **Inherit the mortar lane's brute-force facet pairing** — NOT the bucket sort (which is NTS point-vs-segment only) | Lens-B factual fix; zero new spatial structure; the deferred mortar-P2.5 broad phase is an optional optimization |
| Dimensionality | **3D-only**; a 2D surface never routes here | `n=ê_s×ê_m` needs 3D; 2D crossing = vertex/NTS |
| Discovery generalization | **Out of scope → ADR-55**; ADR-57 owns the *kernel*, ADR-55 *consumes* it | Kernel-ownership collision resolved (Lens-C); adapter discovery-agnostic |

## Risks / open questions

> [!question]
> **Force discontinuity as the contact point migrates off an edge end.** When `s*` or
> `t*` reaches an endpoint, the contact transitions from edge-edge to vertex/NTS; the
> normal can swing. Mitigation: the gap-band + the routing hysteresis + viscous `-visc`
> damping; the E7 sliding-rake gate is the falsifier. If it chatters, that is the
> smoothing re-open trigger (ADR-47 #4a / Munjiza distributed potential) — *not* a defect
> in this design, a documented hand-off.

> [!question]
> **Oblique band routing (align_fn ≈ 0.5).** A 45° facet pair can present *both* a small
> clip and a valid edge-edge pair. The §2 design ships **ordered precedence with a
> guaranteed-owner fall-through** (not a clean partition); the residual double-count is
> bounded by `τ_A` and gated by the E1 NTS-double-count + no-coverage-hole falsifiers.
> Open only: does precedence chatter across steps in the band, forcing the documented
> partition-of-force *blend* fallback? E1's hysteresis sweep decides — not a correctness
> risk (force is always owned), a robustness one.

> [!question]
> **Edge-end handoff chatter.** As `s*`/`t*` reach an edge end, edge-edge hands off to a
> vertex/NTS primitive and `n` can swing. Mitigation: the gap-band + routing hysteresis +
> `-visc`; the E7 sliding-rake gate is the falsifier. If it chatters, that is the
> *Munjiza-distributed-potential* re-open trigger (a **new** ADR-47 row — **not** #4a; see
> §Fence) — a documented hand-off, not a defect here.

- **Near-term use case is narrow — primary value is de-risking the ADR-55 kernel
  (Lens-C).** The headline motivators (fragment edges after element removal, tumbling
  polyhedra) are *fenced to ADR-55*; the in-scope predefined-surface cases (slab corner
  on a beam edge, cross-stacked bars, L/T wall junctions) are real but niche, and several
  are arguably better modeled by mesh-tying/NTS. The honest framing: the predefined-
  surface MVP's chief payoff is **proving the `LadrunoEdgeKernel` kinematics** that ADR-55
  will reuse, with the structural cases as secondary validation — not broad immediate use.
- **No new build/ABI/dep risk** — header-only kernel, existing translation units, serial
  (the whole contact subsystem is serial-only; parallel `m_eff` is BLOCKED on the contact
  handler parallelization, the *other* queued ADR — edge-edge inherits that fence, does
  not widen it).
- **Backwards compatibility:** `-edgeedge` off by default ⇒ every existing model is
  byte-identical; an active edge-edge adapter declares real 4-node connectivity ⇒ the
  numberer permutation differs (the standard active-adapter caveat, capstone contract #6).

## Adversarial review log

**Gate 1 — 2026-06-24 (3 independent reviewers, refute-by-default, code- + primary-source-
grounded). Verdict: SALVAGEABLE → all dispositions folded; now PASS-equivalent.** Three
lenses — (A) kinematics/numerics, (B) substrate-fit/contracts, (C) routing/scope/cross-ADR
— raised **8 MAJOR + 7 MINOR + 4 NIT**. No fatal flaw in the core framing; the
segment-segment skeleton, the substrate seam (mode dispatch, commit/GC hooks, friction/
SOFT/mortar/projection kernel reuse — Lens-B verified each at file:line), and the API all
held. The MAJORs, all folded:

| ID | Lens | Issue | Fix (folded) |
|---|---|---|---|
| A-1 | kin | Endpoint-clamp recovery wasn't Ericson's exact ladder (wrong vertex-vertex sub-region) | §1 now quotes the conditional clamp ladder (RTCD §5.1.9); E0 vertex-vertex falsifier |
| A-2 | kin | `g_N=±‖w‖` true only for strict interior; open-interval flag admits `s*=0.999` vertices | Margin'd interior `[δ,1−δ]`; "to first order" deleted; E0 `\|w·n\|==‖w‖` falsifier |
| A-3 | kin | Cross-product normal `1/sin²θ_edge` ill-conditioned; `τ_∥`=5° too permissive | `τ_∥` set by conditioning (default `sin²15°`), E0-gated; dual-purpose gauge documented |
| A-4 | kin | **Position-derived sign rule masks interpenetration through contact** (self-referential `w·n`) | Body-fixed **committed** sign `EdgeEdgeState.signN`; E2 monotone-through-contact falsifier |
| B-1 | sub | **Broad-phase misattribution** — mortar lane is BRUTE FORCE, not the bucket sort (NTS-only) | §7 rewritten to inherit brute-force facet pairing; mortar-P2.5 noted as optional |
| C-1 | scope | `cos_t` symbol overloaded — in-tree it's the area-Jacobian guard, not a router band | Router scalar renamed `align_fn`; same operands, different role, documented |
| C-2 | scope | **Kernel-ownership collision with ADR-55** (its P3 builds "a new edge-edge kernel") | Fence fix: ADR-57 owns the kernel, ADR-55 P3 demoted to *consume* it (ADR-55 edited) |
| C-3/4 | scope | "Partition, no superposed force" over-claimed — ignores NTS 3rd lane + oblique band | Re-specified as **ordered precedence + guaranteed-owner fall-through**; E1 falsifiers |

`fold` (MINOR/NIT, applied): friction-frame no-op projection → `t₁=ê_s,t₂=n×ê_s` (A);
`d_band` tied to edge-length fraction (A); zero-length edge guard + E0 falsifier (A);
2D-only fence (A); `align_fn`/facet-normal not in `PairResult` → flagged required plumbing
(B); edge-id key = ordinal/composite, never a lossy hash (B); de-dup = *new* set mirroring
the stamp idiom, not a reuse (B); near-parallel-edge-on triple-miss → hard E1 gate +
clip cross-check (C); edge-end smoothing citation corrected #4a → Munjiza distributed
potential, new ADR-47 row (C); capstone component row annotated "routing open E1 gates"
(C); near-term use-case reframed as "de-risk the ADR-55 kernel" (C). `confirmed-sound`
(Lens-B, no change): EDGE_EDGE mode extensibility, commit-hook obligation accuracy,
friction-kernel 4th-consumer, SOFT mass-cache, clip-degeneracy probe safety
(`cos_t<1e-12` guard before the `A_aux/cos_t` division), null-build byte-identity.

## Implementation log

- **E0 (geometry) — SHIPPED** (this PR). `SRC/domain/contact/LadrunoEdgeKernel.h`
  (header-only, OpenSees-free): `closestPtSegSeg` (Ericson RTCD §5.1.9 exact clamp
  ladder + the conditioning-justified `τ_∥=sin²15°` parallel refusal + zero-length
  DEGENERATE guard), `edgeNormal` (common-perpendicular), `edgeGap` (body-fixed
  **committed** sign — the A-4 fix), `bOperator` (envelope-theorem-exact first variation).
  Oracle `proto_e0_closest_point.py` **28/28** + standalone `e0_cpp_check.cpp` **13/13**
  bit-for-bit (B==FD identical to 1.13e-9). Header-only ⇒ not yet included by any TU ⇒
  build byte-identical (the kernel is inert until E1/E2 wire it into the handler/adapter).
  Oracle caught one ADR wording imprecision to fold later: the *normal direction* `∂n/∂u`
  conditioning is `1/sinθ_edge`; the closest-point *parameter* `∂(s,t)/∂u` conditioning is
  `1/sin²θ_edge` (the dominant gauge the §1 "1/sin²" refers to) — both bounded at the band.
- **E1 (routing logic) — oracle SHIPPED** (this PR). `proto_e1_router.py` **16/16** pins
  the routing decision tree. **Design refinement caught oracle-first:** the draft's flat
  "NTS ≻ edge-edge ≻ mortar" precedence steals facing-facet contact from the mortar lane
  (a facing overlap has slave nodes in-bounds — the NTS condition); fixed to an
  **`align_fn` regime split** (facing `align_fn≥τ_face` ⇒ mortar owns; non-facing ⇒
  NTS≻edge≻sliver), folded back into §2. The proximity-gated precedence is a *guaranteed
  owner* (no coverage hole) and chatter-free (1 transition through a 45° sweep). The C++
  realization (the handler routing loop) **folds into E2** — routing without an adapter
  to inject is inert, so it ships with the `EDGE_EDGE` adapter, not before.
- **E2 (penalty normal force + main tangent) — SHIPPED** (this PR). The first live C++/build
  phase: the inert E0 kernel + E1 routing oracle become a real force. **`LadrunoContactFE` gains
  `Mode EDGE_EDGE = 4`** — a stateless 4-node adapter `[a₀,a₁,b₀,b₁]` (`FE_Element(tag, 4, 12)`,
  **no new class tag** — a mode of the one adapter like RIGID_PLANE/SEGMENT/MORTAR): `getResidual`
  runs `LadrunoEdgeKernel::closestPtSegSeg` at the trial config, fetches the body-fixed committed
  sign from the Domain-owned `EdgeEdgeState`, and assembles `f = t_N·B`, `t_N = ε_N⟨−g_N⟩`;
  `addKtToTang`/`addKiToTang` assemble `K = ε_N·BᵀB` (symmetric, rank-1, PSD); `addCtoTang` is a
  no-op (viscous is E5+). **`LadrunoContactDomain::EdgeEdgeState`** (the §5 struct — `signN` + its
  revert double-buffer carrying; friction E3 + one-scalar ALM E6 fields zeroed/inert) is keyed by a
  5-int **composite** `(contactTag, ordered slave-edge nodes, ordered master-edge nodes)` — never a
  lossy hash — with `getOrCreateEdgeEdgeState` + the `edgeGC{Begin,Mark,End}` dead-slot trio; and
  `commit()`/`revertToLastCommit()` were **extended to iterate `theEdgeEdgeStates`** (the explicit
  obligation — the shipped hooks iterated only friction + mortar, exactly the §5 trap). The
  **handler routing loop** (the E1 `align_fn` regime-split precedence, now non-inert) runs over the
  mortar lane's **brute-force** facet pairing (the §7 fix — *not* the bucket sort): on an EDGE_EDGE
  route it enumerates the `(≤4)×(≤4)` boundary-edge pairs, keeps the margin-interior/well-conditioned
  /in-band crossings (the E0 kernel decides), de-dups by the ordered edge-node-tag tuple, and injects
  one `EDGE_EDGE` adapter per kept pair. The **`-edgeedge`** opt-in modifier on `-mortar`
  (`-edgeKn auto|val`, `-edgeBand <d>`) is **OFF by default ⇒ byte-identical**; on a perpendicular
  pair the mortar adapter is geometrically inert (the clip degenerates, the very gap edge-edge fills)
  so the two never double-count. Oracle `proto_e2_penalty.py` **23/23**; `test_adr57_edge_edge_1`
  **4/4** (crossed-bars restrained + the NTS-passes-through headline falsifier + the oblique-band
  EE-1 regression + byte-identity); full contact battery **145 passed** (no regression).
  **3-reviewer adversarial gate (Gate 2 — refute-by-default) → PASS, 2 findings folded:** **EE-1
  (MAJOR)** — the mortar lane injected its adapter unconditionally, so in the *oblique* band (50°–90°,
  the clip is non-degenerate, not just at exactly 90°) the mortar pressure and the edge-edge force
  superposed (a double-count the "mortar is inert on a perpendicular pair" note over-claimed). Folded:
  a `ladrunoEdgeEdgeOwns` **single-source-of-truth** routing helper the *mortar* loop now consults to
  **stand down** when edge-edge owns the pair (the §2 "face-mortar vacates only if edge-edge succeeds",
  now implemented), + the `test_e2_oblique_no_double_count` regression. **E2-2 (MINOR)** — the default
  sign-anchor `orientDir` used a centroid difference that can vanish near-coincident centroids; folded
  to the always-non-degenerate slave **facet normal** (the §2 "orient once from {n̂_s, n̂_m}"). Rejected
  (sound / out-of-scope): the per-pair rank-deficiency (inherent to penalty contact — NTS/mortar share
  it; the host elements or a regularizer supply the differential mode) and an unreachable null-deref.
- **E3→E7** — pending (design-only above; each lands oracle-first → C++ → gate → PR,
  updating the capstone row + ledgers in the same PR).
