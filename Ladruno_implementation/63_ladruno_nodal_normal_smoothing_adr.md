# ADR-63 — Averaged Nodal-Normal Smoothing for NTS Contact (smooth `N(X)` master-normal field)

- **Status:** **P0+P1 SHIPPED (#457); P2.1 facet-ownership SHIPPED (#460); P2.2+P2.3 friction/re-emit VALIDATED (local, 2026-07-01)** — the
  design gate (below) is complete and the symmetric-first MVP is implemented and green. **Q-TANGENT
  RESOLVED:** ship the **symmetric frozen-field `kn·BᵀB`** tangent first; the full `∂n_smooth/∂u` is a
  conditionally-required, evidence-gated P3 (see §"Gate decision (RESOLVED)"). **P2.1 closes the one real P1
  limitation** — the spurious adjacent-facet activation at a sharp convex ridge (see the "NEW finding at
  bring-up" §, now RESOLVED, and the "P2.1 implementation notes" § below): a per-segment shared-edge
  ownership guard, **gap-aware** (reject a smoothed projection on a SHARED interior edge only when its
  penetration is large vs the facet size ⇒ keeps the owner + the at-apex contact; a blunt near-edge reject
  reproduced the reverted pass-through). Validation: build-free oracle `proto_nodal_normal_selfcheck.cpp`
  **35/35**; in-solver `tests/test_adr63_smoothnormal_p1.py` **5/5** + `tests/test_adr63_smoothnormal_p2.py`
  **2/2** (quad + tri convex-ridge press-into: no spurious ejection); full contact battery
  (ADR-39/41/57/60/63) **152** with the feature OFF byte-identical. **P2.2+P2.3 (friction + re-emit over a
  curved master) — VALIDATED, NO engine code**: `segmentActive` already threads the smoothed `n` into the
  gap operator AND the friction slip (`tangentPart(drel, n, gTvec)`), and `-reemit`/`-smoothNormal`/`-mu`
  compose with no mutual refusal — so a frictional block dragged across a convex curved master with
  `-reemit -smoothNormal -mu -outward` sustains contact across the crest, vs pass-through without `-reemit`
  (the ADR-60 "exposed combo" closed); flat-master friction is smooth==faceted byte-identical
  (`tests/test_adr63_smoothnormal_p23.py`, 3/3; battery **157**). Documented caveats (not new bugs): the
  AUTO global sign stays ill-conditioned for edge-grazing slave clouds (F2/F3/F5 → `-outward`); mild
  near-apex over-stiffness on sharp crests (the P2.1 gap-aware guard keeps the small-gap non-owner;
  non-destabilizing on realistic arcs; full single-owner = ADR-57 #4b); traction-continuity asserted
  qualitatively via sustained contact. **Q-IMPLICIT-NEWTON RESOLVED — outcome (a), CONVERGES (local,
  2026-07-01):** the implicit quasi-static rig (`tests/test_adr63_smoothnormal_p24_implicit.py`, 3/3;
  StaticAnalysis + DisplacementControl dragging a slave across a FLAT middle facet whose *smoothed* normal
  rotates ~15–22° from steep-neighbour nodal-normal tilt) shows the frozen-field symmetric `kn·BᵀB`
  Newton converges in **2 iterations every step**, *independently of load-step coarseness* (whole facet
  in a SINGLE step = maximal within-step rotation ⇒ still 2 iters), and up to steep facets / kn≤1e6 /
  press≤3000. The dropped `∂n_smooth/∂u` is `O(kn·gN)=O(press)` — scaled by the small penalty penetration
  `gN≈press/kn` ⇒ sub-dominant to `kn·BᵀB`; it only strains at LARGE penetration (`gN`≳15% of the facet),
  where the FACETED `-geomtan` consistent tangent diverges too (penalty NTS misused, NOT a smoothed-normal
  defect ⇒ P3 wouldn't rescue it). So **steelman-B's stall concern is empirically REFUTED and P3 stays
  deferred WITH EVIDENCE.** The `-smoothNormal`+`-geomtan` silent-downgrade warning is present and fires.
  Battery **160** OFF byte-identical; test-only, no engine code, no classTag. **Remaining:** P3 (full
  `∂n_smooth/∂u`, gated — now a genuinely-optional follow-up), CI port; separate follow-ups: an auto-sign
  robust/per-component vote for edge-grazing slaves, and a quantitative junction-traction-continuity gate.
- **Owner:** Mora Bowen · Palacios · Abell · Guppi
- **Priority:** high — **resolves the open [[60_ladruno_finite_sliding_reemission_adr]] R3 item** (curved-master
  `orientDir` flip → silent pass-through, today's only `-outward`-caveated combo) and the in-fork
  [[41_ladruno_mortar_alm_contact_adr]] **Q-NORMAL** faceted-normal chatter for the NTS lane.
- **Parents / prior art:** [[39_ladruno_contact_domain_adr]] (NTS broad+narrow phase, the `normalOriented`
  derived-normal sign machinery), [[60_ladruno_finite_sliding_reemission_adr]] (finite-sliding re-emit; **R3**
  is the load-bearing trigger for this ADR), [[41_ladruno_mortar_alm_contact_adr]] (Q-NORMAL statement of the
  chatter problem), [[47_ladruno_contact_deferrals_adr]] (rows **#4 / #4a / #4b** — this ADR is **#4a pulled
  forward**), [[48_ladruno_contact_capstone_adr]] (status-of-record).
- **Pull-forward:** ADR-47 ledger row **#4a** ("averaged nodal-normal smoothing, smooth `N(X)` field") — the
  cheap interim **preferred over** full slide-line Hermite smoothing (#4). Mark #4a `pulled forward → ADR-63`
  in ADR-47 in the shipping PR.
- **Class tags:** none new. Behavior of `LadrunoContactHandler`
  (`HANDLER_TAG_LadrunoContactHandler = 33002`) + `LadrunoContactDomain` + a new header-only geometry module.

---

## What

Replace the NTS contact normal — today a **per-segment faceted** normal whose **outward sign** is chosen
**per (slave, segment) pair** from the slave's reference side (auto `orientDir`) or a fixed `-outward`
vector — with a **smooth nodal-normal field** `N(X)` on the master surface, opted in per contact with
`-smoothNormal`:

1. **Per master node** `a`, average the **consistently-wound** adjacent segment normals (area- or
   angle-weighted) into a single nodal normal `n_a` (Abaqus TG §5.1.1, verbatim: *"Unit normals are computed
   at every master node (averaging adjacent segment normals)… giving a smoothly varying normal field `N(X)`
   over the master surface"*).
2. **At a query point** with master parametric coords `(ξ,η)` on segment `s`, interpolate the segment's nodal
   normals with the **segment shape functions**: `n_smooth(ξ,η) = normalize(Σ_i N_i(ξ,η) · n_{a_i})`.

This single change kills **both** ADR-60 R3 failure modes at once:

- **(a) the ridge flip (R3, HIGH).** The outward sense becomes a **global property of the surface**
  (its coherent winding + an `-outward` / slave-centroid seed), captured **once**, not a per-slave datum.
  A slave that slides over a convex ridge can no longer flip *its* pair's sign → the silent pass-through is
  **structurally impossible**, not merely persisted-around (the ADR-60 BLOCKER-ORIENTDIR persistence was a
  *partial* fix; a global surface sign is the *complete* one — see §How D2).
- **(b) the junction chatter (ADR-41 Q-NORMAL).** Adjacent segments share the **same** nodal normals at
  their shared edge nodes, so `n_smooth` is **continuous (C0) across segment junctions** — the faceted
  `n`-jump that makes the contact point oscillate (Abaqus TG §5.1.2) is removed. The ADR-60 P1 continuity
  gate can tighten from its current LOOSE resultant tol (forced by the faceted jump) toward a continuous-force
  assertion.

**Honest continuity claim:** nodal-normal interpolation gives a **C0-continuous normal field** (the normal
itself is continuous; its tangential derivative still jumps at a junction). It is **not** the C1 *surface*
of full Hermite slide-line smoothing (#4, TG §5.1.2). #4a removes the **jump** (the chatter driver); the
residual slope discontinuity is the gap #4 would close and is the documented limit of this ADR.

**Scope (NTS lane only).** This ADR builds the reusable nodal-normal machinery and wires it to the **NTS**
penalty path (the lane with the silent-pass-through *correctness* hole). See §Fence for the mortar /
edge-edge boundary.

OUT: full slide-line **Hermite** C1 smoothing (#4 — heavier, only if the P3 chatter gate proves #4a
insufficient); **mortar** GP-normal smoothing (the same header is reusable, but mortar already finite-sliding
*correct* — its Q-NORMAL is a quality, not correctness, gap; a future pass / amendment); **edge-edge** normal
(its own body-fixed sign `signN`, ADR-57); closed / self-enclosing / non-manifold master auto-orientation
(refused, fall back to faceted `-outward` — see BLOCKER-WINDING).

---

## Why

`LadrunoContactProjection::normalOriented()` (`LadrunoContactProjection.h:139-157`) derives a **single faceted
normal** `n = normalize(g1 × g2)` for a segment and flips it to satisfy `n·refDir > 0`, where `refDir`
(`orientDir`) is supplied by the handler. Two structural defects follow:

1. **The sign is per-slave (R3).** `LadrunoContactHandler.cpp:663-674` computes `orientDir` **per (slave,
   segment) pair**, each `handle()`, from the **reference** config: explicit `ct.outward` if given, else
   `auto = (slave ref coords) − (segment ref centroid)`. On a **sharp convex ridge**, a slave that slides
   from one face to the other has its reference-side vector point to the **opposite half-space** after the
   crossing → the derived normal **flips sign** → the repulsive normal becomes attractive → **silent
   pass-through** (zero/inverted force, no warning). ADR-60 R3 (HIGH) documents this; the interim is
   `-outward` (a fixed global vector — but a single vector cannot point outward on *both* faces of a ridge,
   so `-outward` only papers over the case where the two faces are near-coplanar). The convex-ridge gate
   test is **absent** today.
2. **The normal is C0-discontinuous across junctions (Q-NORMAL).** Each segment owns an independent facet
   normal, so `n` **jumps** at every master inter-element boundary. Under sliding the active contact point
   oscillates between segments (Abaqus TG §5.1.2) → **force chatter**. Today the ADR-60 P1 continuity gate
   asserts only a **LOOSE** resultant tolerance precisely because of this documented jump
   (`60_…adr.md` P1 oracle (b), citing Q-NORMAL).

ADR-60's deferred BLOCKER-ORIENTDIR (persist `orientDir` per `(contactTag, slaveTag)` from the first
reference handle) was correctly diagnosed in the ADR-60 post-code review as a **partial** fix: *"Persistence
of `orientDir` alone is only a PARTIAL fix — it stops the penetrating-rederive inversion but NOT a legitimate
ridge crossing (the persisted reference direction is for the wrong half-space after the slave crosses). Full
fix needs consistent-winding normals or nodal smoothing (ADR-47)."* This ADR is that full fix.

---

## Where

- **New** `SRC/domain/contact/LadrunoContactNormalField.h` — **header-only, OpenSees-free** (mirrors
  `LadrunoContactProjection.h` / `LadrunoContactReemit.h` / `LadrunoContactKernel.h`: raw `double[]` + STL,
  no OpenSees types, so the engine **and** a build-free oracle share **one** logic). Three pure stages:
  - `buildAdjacency(mTags, nSeg, nps, …)` — the **segment-edge adjacency graph** (segments sharing an
    unordered node-pair edge) + a **non-manifold / connectedness** report.
  - `propagateOrientation(adjacency, seed, …)` — **flood-fill manifold orientation**: a per-segment winding
    sign `σ_s ∈ {+1,−1}` giving globally-coherent winding from a seed, **topological only** (depends on
    connectivity, not coordinates → computed **once per membership fingerprint**).
  - `nodalNormals(segCoords, σ_s, weight, globalSign, …)` → `n_a[node][3]`: area- or angle-weighted average
    of the consistently-wound adjacent facet normals, plus `smoothNormal(seg, ξ, η, n_a) → n` (the shape-fn
    interpolation + renormalize). Oracle: `contact_prototypes/proto_nodal_normal.py` + a standalone main.
- **Modified** `SRC/domain/contact/LadrunoContactProjection.h` — a **new** entry point
  `evalSegmentSmooth(nps, Xseg, xs, nodalNorm[4][3], gap, n, N, …)` paralleling `evalSegment()`: facet
  closest-point projection for `x̄,(ξ,η)` (unchanged), then `n = smoothNormal(ξ,η)` and `gap = n·(xs − x̄)`.
  The existing `evalSegment()`/`normalOriented()` stays **byte-verbatim** for the OFF path.
- **Modified** `SRC/domain/contact/LadrunoContactDomain.{h,cpp}` — per **master surface** (keyed by surface
  tag): the cached topological winding `σ_s` + the per-handle nodal-normal field `n_a[…]` + an
  orientation **status** (oriented OK / refused-non-manifold / refused-closed) + the **captured global sign**.
  Built **once per handle** from current coords (like the bucket grid), **strictly under** `if(smoothNormal)`,
  **guarded by the R1 membership fingerprint** (`theReemitFp` machinery, [[60_…]] R1) so a re-mesh recomputes.
  Config field `smoothNormal` with a **default member initializer** (Lens-B; stack garbage must not flip it on).
- **Modified** `SRC/analysis/handler/LadrunoContactHandler.cpp` — NTS loop only: when `ct.smoothNormal`, hand
  the engine the master `segCoords`, request the field, and construct each `LadrunoContactFE` with a pointer
  to its segment's four nodal normals (and the OFF path unchanged). Capture the **global sign** once (D2).
- **Modified** `SRC/analysis/handler/LadrunoContactFE.cpp` — `segmentActive()` calls `evalSegmentSmooth()`
  with the segment's nodal normals when smoothing is on (else `evalSegment()` verbatim); the `B`-operator and
  gap use the smoothed `n`. Adapter stores the four nodal normals (or a fallback `orientDir` if refused).
- **Vanilla footprint:** **none.** No `Domain.cpp` touch (unlike ADR-60); this is entirely within the
  contact subsystem + the handler.
- **Reference sources (faithful):**
  - **Abaqus TG §5.1.1** (small-sliding) — the averaged-nodal-normal `N(X)` field is *verbatim* this ADR's
    construction. (TG §5.1.2 is the heavier Hermite C1 #4.) Read recipe in the `abaqus-theory-contact-loading`
    skill; `references/contact-kinematics.md` §"Small-sliding".
  - **LS-DYNA Theory §26.11** (bucket-sort re-search + surrounding-segment local tracking, pp.545–550) — the
    global-research-paired-with-smoothed-local-surface architecture, already cited by ADR-60. *NB:* LS-DYNA
    Theory **§17.1.4 "Surface Smoothing" is mesh (ALE) smoothing, NOT contact** — do **not** cite it for this;
    LS-DYNA's contact-surface smoothing lives in its segment-based contact, which this ADR does not need to
    cite beyond the §26.11 architecture parallel. (Marked to avoid the miscitation trap.)

---

## How

### Mechanism (one field, computed once per handle, frozen within the step)

1. **Topological winding (once per fingerprint).** From `mTags` build the segment-edge adjacency graph; pick
   a seed segment; flood-fill a per-segment winding sign `σ_s` so all segments are coherently wound. `σ_s`
   depends only on connectivity → cached and **invalidated only when the [[60_…]] R1 membership fingerprint
   changes** (re-mesh / element removal). Refuse (named error) on a **non-manifold edge** (>2 segments share
   an edge), a **closed / self-enclosing** surface, or a **disconnected multi-shell without `-outward`**
   (BLOCKER-WINDING) → fall back to the faceted `-outward` path for that surface.
2. **Global outward sign (captured once).** Fix the surface's single global sign from **`-outward` if given**,
   else a **one-time global vote**: `sign = sgn(Σ_a n_a · (slave-surface reference centroid − master surface
   centroid))`. The sign is a property of the **surface (+ slave cloud)**, captured at the first handle and
   held — it is **never** a per-slave datum and **never** re-derived from a penetrating current config
   (BLOCKER-GLOBAL-SIGN). This is what makes the ridge flip structurally impossible.
3. **Nodal normals (each handle, current coords).** `n_a = normalize(Σ_{s∈adj(a)} w_s σ_s (g1×g2)_s)`,
   `w_s` = facet area (default) or incident angle. Computed from `getCrds()+getDisp()` so the field rotates
   with the deformed master handle-to-handle (the deformed feed re-emit already uses).
4. **Query (per residual eval).** `segmentActive()` projects the slave onto its segment (facet closest point,
   unchanged), then `n_smooth = normalize(Σ_i N_i(ξ,η) n_{a_i})`, `gap = n_smooth·(xs−x̄)`,
   `B = [n_smoothᵀ | −N_i n_smoothᵀ]`. The nodal-normal field is **frozen within the step** (recomputed only
   at the next handle), so the within-step tangent has a well-defined posture (D4).

### Decisions

- **D1 — Geometry = averaged nodal normal + shape-fn interpolation (Abaqus §5.1.1).** `n_a` from
  consistently-wound, area-weighted adjacent facet normals; `n_smooth` from the segment shape functions at the
  projection. **C0-continuous normal** (no jump across junctions); **not** C1 (the #4 Hermite gap, deferred).
- **D2 — Outward sign is GLOBAL, captured once (the R3 fix).** Surface-level sign from `-outward` or a
  one-time slave-centroid vote — never per-slave, never from a penetrating config. **Subsumes** ADR-60
  BLOCKER-ORIENTDIR (the persistent per-(contactTag,slaveTag) orientDir): a single global surface sign is
  strictly stronger than a persisted per-slave one (a persisted per-slave direction is still wrong for the
  far half-space after a legitimate crossing). `-outward` keeps working as the explicit seed **and** the
  refused-surface fallback.
- **D3 — Winding by flood-fill manifold orientation (the hard part).** `σ_s` from edge-adjacency BFS; a shared
  edge traversed in the **same** node order by both segments ⇒ opposite winding ⇒ flip. Topological → cached
  per fingerprint. Refuse non-manifold / closed / multi-shell-without-outward (BLOCKER-WINDING).
- **D4 — Tangent: ship the SYMMETRIC `kn·BᵀB` with a frozen field first (RESOLVED — see §Gate decision).**
  With the nodal-normal field **frozen within the step**, the consistent tangent of the penalty normal is
  exactly `kn·BᵀB` evaluated with the smoothed `B`. Differentiating the *frozen-geometry residual*
  (`r = kn·gN·B(n̄)`, `gN = n̄·(xs−x̄)`) treating `n̄` as the constant it literally is gives `∂r/∂u =
  kn·B(n̄)ᵀB(n̄)` exactly — so this is an **exact modified-Newton on per-step-frozen geometry, NOT an
  inconsistent-Newton with a residual/tangent mismatch**: the residual normal *is* the same frozen `n̄` the
  tangent assumes, the converged solution is **tangent-independent**, and quadratic convergence to the
  frozen-geometry fixed point is preserved (the only dropped piece is the `∂N_i/∂u` facet-projection block the
  faceted default drops too). The dropped configuration dependence is the **same *kind*** of approximation the
  shipped faceted default already makes (the per-segment `normalGeomTangent` B3 block is itself opt-in via
  `-consistentNormal`) — **but its magnitude is generically *larger*** for the smoothed field on curved
  masters: the faceted tri-3 facet block is `≡ 0` (a planar facet, `κ ≡ 0`,
  `LadrunoContactProjection.h:230,238`), whereas the smooth field deliberately reintroduces a non-zero
  normal-rotation (the very term that fixes chatter is the term being dropped), so a curved-master Newton
  **stall is more probable here** than on the piecewise-flat faceted path, and is the documented trigger for
  P3 (Q-IMPLICIT-NEWTON).
  - **Symmetry advantage is FRICTIONLESS-IMPLICIT-scoped.** The `kn·BᵀB` Gram form is symmetric and keeps the
    symmetric solver — but **only matters for frictionless implicit**: under `-mu>0` sliding the friction
    tangent (`LadrunoFrictionKernel`, Abaqus TG §5.2.3) is **already non-symmetric**, so an implicit frictional
    run **already** requires an unsymmetric solver (UmfPack/FullGeneral) regardless of the normal-tangent
    posture. The **load-bearing** arguments for symmetric-first are therefore: **(i)** the primary NTS lane is
    **explicit** (`CentralDifferenceLadruno`), where the contact FE returns the integrator-routed (**zero**)
    tangent and contributes a **residual only** — the full `∂n_smooth/∂u` would **never be evaluated** in the
    load-bearing lane; and **(ii)** it minimizes the P1/P2 adversarial + assembly surface.
  - **The full `∂n_smooth/∂u`** (differentiating the nodal-normal average **and** the interpolation through
    the mesh motion) expands the contact `FE_Element` connectivity from `{slave + segment}` (**≤5 nodes /
    15 DOF** — the current `normalGeomTangent` footprint, `Kgeom[15][15]`, `LadrunoContactProjection.h:245-249`)
    to `{slave + segment + master 1-ring}` (**~10–17 nodes**), grows the dense-block area **~4–12×**, **re-colors
    the sparse graph on every active-set / contacted-segment change under sliding**, and is **non-symmetric
    under slip** (Abaqus §5.1.2 second variation) ⇒ unsymmetric solver. It contributes **nothing** to the
    explicit lane (no tangent assembled) and is therefore **deferred to a conditionally-required, gated P3**
    (`-consistentNormalSmooth`), **implicit-only** (see BLOCKER-TANGENT-POSTURE).
- **D5 — Opt-in, OFF byte-identical.** `-smoothNormal` (default off, default-member-initialized). OFF ⇒ no
  field built, `evalSegment()`/`normalOriented()` verbatim ⇒ **byte-identical** (the full contact battery is
  the gate, exactly as ADR-60 P0). ON over a **flat** master ⇒ all facet normals equal ⇒ `n_smooth ≡ n_facet`
  to round-off (a **consistency** assertion `‖n_smooth − n_facet‖ < 1e-14`, *not* a byte-identity claim).
- **D6 — Scope = NTS lane.** The header is lane-agnostic and reusable, but only the NTS adapter is wired here
  (the lane with the *correctness* hole). Mortar / edge-edge fenced (§Fence).
- **D7 — Composition with re-emit (ADR-60).** The field is built from the **same deformed feed** re-emit uses
  and **invalidated by the same R1 fingerprint**; the global sign capture composes with re-emit's anchor
  capture. With D2 the **R3 `-outward` caveat is lifted**: `-reemit + -mu + curved master` becomes safe
  (§"Interaction with ADR-60").

### Why not just persist `orientDir` (the ADR-60 deferred fix)?

Persisting the per-(contactTag,slaveTag) `orientDir` from the first reference handle (ADR-60
BLOCKER-ORIENTDIR) stops the *penetrating-rederive* inversion but **not a legitimate ridge crossing**: the
persisted reference direction points into the **wrong half-space** once the slave is genuinely on the far
face. A **global surface sign** (D2) has no per-slave direction to be wrong — the same outward field serves
every slave on the surface. That is why this ADR resolves R3 and the persistence patch does not. (Recorded so
the cheaper patch is not re-proposed.)

---

## Design-gate BLOCKERs (status after the gate)

1. **BLOCKER-WINDING (the hard part)** — auto-orientation needs a **coherent manifold winding**.
   `propagateOrientation` must **refuse** (named error, fall back to faceted `-outward`, do **not** silently
   smooth) on: a **non-manifold edge** (>2 segments share an edge), a **closed / self-enclosing** surface
   (a single global `-outward` cannot point outward on all faces), or a **disconnected multi-shell without
   `-outward`** (each component's global sign is independently ambiguous). The supported MVP target is a
   **single open connected manifold patch** — the curved-master R3 case (ramp / cylinder section / dome cap).
   P1 gate includes a non-manifold refuse-and-fallback test. **Folded:** a `-outward` that resolves *some*
   shells but not others (a partially-aligned multi-shell) refuses the unaligned components only.
2. **BLOCKER-GLOBAL-SIGN** — the outward sense **must** be a single surface (+slave-cloud) datum captured
   once (D2); a per-slave or per-segment sign re-introduces R3. The vote uses the **slave-surface reference
   centroid**, not any individual slave's live position. P1 gate = convex-ridge crossing keeps a
   **non-flipped, repulsive** normal across ≥3 crossings.
3. **BLOCKER-IDENTITY** — OFF must be **bitwise identical** AND assert (a) no field allocated/built when
   `smoothNormal` off, (b) `evalSegment` path byte-verbatim, (c) the `smoothNormal` struct field has a
   **default member initializer** = `false` (Lens-B; uninitialized POD could flip it true), (d) the full
   ADR-39/41/57/60 contact battery byte-identical. P0 gate.
4. **BLOCKER-FALLBACK** — never emit a **zero/garbage** normal, at **either** of two distinct degeneracies:
   **(a)** a **degenerate nodal normal** (`‖Σ w_s σ_s n_s‖ → 0`, e.g. a sharp fold where two facets nearly
   cancel, or a node with one sliver facet), and **(b)** a **degenerate interpolated normal at the query
   point** (`‖Σ N_i(ξ,η) n_{a_i}‖ → 0` — antiparallel *corner* normals blended across a folded element, a
   failure the per-node guard (a) does **not** catch because each `n_a` can be individually unit yet their
   shape-fn blend collapses). Both — and a **refused surface** — fall back to the per-segment
   `normalOriented()` for that node/segment/query with a **one-time named warning**. The smoothed path is an
   *upgrade*, never a *regression* to undefined. *(Note: the frozen-field tangent silently rides through (b)
   with a slightly-wrong direction; the full `∂n_smooth/∂u` would instead expose it as a near-singular
   `(I−n⊗n)/‖·‖` block — another reason the frozen MVP is the more robust first ship.)* P1 oracle adds a
   folded-element query where corner normals nearly cancel → must fall back, not renormalize garbage.
5. **BLOCKER-PROJECTION-CONSISTENCY** — keep the **facet** closest-point projection for `x̄,(ξ,η)`
   (the slave-to-projection vector is `∥` the facet normal, not `n_smooth`); evaluate `n_smooth` at `(ξ,η)`
   and form `gap = n_smooth·(xs−x̄)`. The resulting normal/projection inconsistency is **second-order in the
   facet inter-segment angle** (Abaqus §5.1.1 instead anchors the projection *on* the smooth field — the
   heavier construction, deferred). **Critically: the in-bounds + penetration gate still uses the facet
   projection**, so [[60_…]] R5 **slide-off** behavior (force → 0 cleanly off the edge, no phantom traction)
   is **unchanged** — smoothing alters the *direction*, not the *active set*.
6. **BLOCKER-TANGENT-POSTURE** — the MVP ships the **symmetric `kn·BᵀB`, frozen-field** tangent (D4),
   explicit-OK; the field is frozen within the step so the within-step `∂n/∂u = 0` (**modified-Newton, not
   inconsistent-Newton** — the residual normal is the byte-same frozen field the tangent assumes; the
   converged answer is tangent-independent). Three gated invariants:
   - **(frozen-field consistency, P0/P1 assert)** the residual `n̄` MUST be the **byte-same** handle-frozen
     field the tangent assumes — `segmentActive()` must **not** re-derive `n_smooth` from live coords mid-residual.
     P1 asserts `residual-n̄ ≡ tangent-n̄` (one frozen normal per step, threaded through `segmentActive`).
     Without this the consistency the whole verdict rests on is unenforced.
   - **(explicit prohibition, P3)** the full `∂n_smooth/∂u` (1-ring stencil, non-symmetric) is **implicit-only**
     and MUST NOT be assembled under an explicit integrator — returning a non-zero tangent corrupts the
     mass-only effective system (the P2a `getTangent`-returning-K trap, [[39_…]]). Under explicit the contact
     FE returns the integrator-routed (zero) tangent and the smoothed-normal change is **residual-direction-only**.
   - **(convergence tripwire, P2)** the P2 Q-IMPLICIT-NEWTON gate MUST exercise a **genuinely curved** (not
     piecewise-flat) implicit master; if the frozen-field tangent fails to **converge** (divergence/stall, not
     merely slow), **P3 is promoted from optional to required** and curved-implicit smoothed-NTS is documented
     **explicit-first until P3 ships**.
   *Gate decision RESOLVED (§below): symmetric-first.*

---

## Phased rollout (oracle-first)

| Phase | Deliverable | Validation oracle |
|---|---|---|
| **P0** | `LadrunoContactNormalField.h` (adjacency / flood-fill orientation / weighted nodal normals / `smoothNormal` query) + engine field store keyed by surface tag + R1-fingerprint-guarded build + `-smoothNormal` parser (default off, default member initializer) + `evalSegmentSmooth` entry. **ZERO behavior change.** | numpy/standalone: edge-adjacency + flood-fill winding on a mixed-winding patch (all `σ_s` coherent); non-manifold / multi-shell **refuse** reported; nodal-normal average vs hand-computed; `smoothNormal` interp continuity at a shared edge. OpenSees: full contact battery **byte-identical**; no field built with the flag off (BLOCKER-IDENTITY). |
| **P1** | Wire the NTS adapter to `evalSegmentSmooth` (**frictionless**); **global-sign capture** (D2); **refuse + fallback** (BLOCKER-WINDING/FALLBACK); symmetric frozen-field tangent. | numpy: slave dragged over a piecewise-flat **convex ridge** — `n_smooth` stays repulsive, **no sign flip** (the no-pass-through, tight). OpenSees: **convex-ridge gate** — two master facets at a sharp angle (e.g. 90°), a slave dragged over the ridge: **(a) FAILS today** (auto `orientDir` flip → pass-through, `min_z` dives), **PASSES on** (`min_z` bounded, contact maintained across ≥3 crossings); **(b) chatter gate** — resultant slave force **continuous** across the ridge to a **tighter** tol than the faceted LOOSE tol (cite the removed Q-NORMAL jump); **(c) flat-master consistency** `‖n_smooth−n_facet‖<1e-14`; **(d) non-manifold surface → refuse + faceted fallback** (no crash, warning surfaced); **(e) slide-off unchanged** (R5 regression still green). |
| **P2** | NTS **friction** under the smoothed normal; compose with **re-emit** (lift the R3 `-outward` caveat); membership-fingerprint recompute on change. **+ the BLOCKER-TANGENT-POSTURE convergence tripwire.** | OpenSees: dragged block **with friction** over the ridge — traction continuous, no spike at the junction; **`-reemit + -mu + curved master`** now sustains contact across crossings (the ADR-60 "exposed combo" closed); energy balance closed; OFF still byte-identical. **Q-IMPLICIT-NEWTON convergence gate** — a **genuinely curved** (NOT piecewise-flat) implicit master (quasi-static curved bearing / dome cap under Newmark, `-smoothNormal -consistentNormal`): record Newton iteration count vs the faceted `-consistentNormal` baseline; if frozen-field converges comparably → B's concern empirically refuted, P3 stays deferred *with evidence*; if it **stalls** → trip the P3-promotion. *(Friction note: under `-mu>0` the implicit tangent is already non-symmetric (Abaqus §5.2.3), so this gate measures iteration count, not solver symmetry.)* **Emit a one-time warning when `-smoothNormal` + `-consistentNormal` are both set pre-P3** ("smoothed-normal consistent tangent not yet available; using frozen-field modified-Newton — ADR-63 P3") so the silent-downgrade trap is a *documented limitation*, not a surprise. |
| **P3** *(optional, gated — promoted to **required** if the P2 convergence tripwire trips)* | Full `∂n_smooth/∂u` consistent tangent (`-consistentNormalSmooth`, 1-ring stencil ~10–17 nodes, non-symmetric, **implicit-only — never assembled under an explicit integrator**) for quadratic Newton on strongly-curved implicit masters; ledgers/banner/quirks; mark ADR-47 #4a `pulled forward → ADR-63`; flip **ADR-60 R3 → RESOLVED**. | numpy: analytic `∂n_smooth/∂u` == FD ~1e-8 on a warped patch (mirrors the proto_b3 discipline). OpenSees: Newton iteration-count drop vs the frozen-field tangent on a curved master; **converged answer unchanged** (tangent-independence is the invariant); explicit lane byte-unaffected (no tangent assembled). |
| **P4** *(future / amendment)* | Reuse `LadrunoContactNormalField.h` in the **mortar GP loop** (the original ADR-41 Q-NORMAL #4a target) — separate ADR or ADR-41 amendment, not this scope. | mortar sliding-patch chatter gate (ADR-41 P3). |

---

## Fence (scope vs the contact program)

- **NTS lane only.** The mortar lane derives its GP normal from `LadrunoMortarKernel::facetNormal`
  (its own clip); that is a **separate** Q-NORMAL instance and mortar is **already finite-sliding correct**
  (brute force), so its chatter is a *quality* gap, not the *correctness* hole this ADR closes. The header is
  built lane-agnostic so a future mortar pass (P4 / ADR-41 amendment) reuses it — but **no mortar wiring
  here**.
- **Edge-edge (ADR-57)** carries its own **body-fixed normal sign** `signN` (`EdgeEdgeState`); it has no
  shared nodal-normal field and is **out of scope**.
- **#4 full Hermite slide-line smoothing** stays deferred in ADR-47; re-open **only if** the P1/P3 chatter
  gate shows the C0 nodal-normal field (this ADR) leaves residual oscillation the C0 field cannot damp.
- **ADR-60 R3** is **resolved by P1** (the global-sign field) — its row flips to RESOLVED in the shipping PR;
  the `-outward` interim becomes the *fallback-on-refuse* path, not the recommended curved-master route.

---

## Interaction with ADR-60 (re-emit) — the R3 caveat lifted

ADR-60's "Exposed combo today" is `-reemit + -mu` on a **curved** master (R3: `orientDir` re-derivation can
flip a convex-ridge normal). With this ADR's **global surface sign** (D2):

- **R3 (a) ridge flip — eliminated.** No per-slave sign exists to flip; a slave crossing the ridge sees the
  same global outward field. `-reemit + friction + curved` becomes safe.
- **R1 membership fingerprint composes.** The winding `σ_s` and nodal-normal topology are keyed by the **same**
  `membershipFingerprint(mTags,n,nps)` ([[60_…]] R1); a re-mesh that reorders `mTags` invalidates **both** the
  friction slots **and** the normal field, recomputing orientation (and re-checking the refuse conditions —
  a re-mesh could change manifold-ness).
- **Deformed feed composes.** The field is built from `getCrds()+getDisp()` each handle — the same committed
  deformed config the re-emit grid uses — so the slave slides on the **deformed smoothed** master.
- **Persistence subsumed.** ADR-60's deferred BLOCKER-ORIENTDIR (persist per-slave orientDir) is **retired**
  in favor of D2's global sign (strictly stronger). The ADR-60 R3 row should cite this ADR as the resolution.

---

## Ledger / classTag / banner bookkeeping (for the shipping PR)

- **classTags:** none new (handler stays 33002 — behavior-only, like ADR-60's re-emit and the EDGE_EDGE mode).
- **`LEDGER_implementations.md`:** row — *Averaged nodal-normal smoothing (NTS contact `N(X)` field)*, files =
  `LadrunoContactNormalField.h` (new) + `LadrunoContactProjection.h` + `LadrunoContactDomain.{h,cpp}` +
  `LadrunoContactHandler.cpp` + `LadrunoContactFE.cpp`, status per phase, PR.
- **`LEDGER_vanilla_files.md`:** **no new row** — zero vanilla touch (all changes within the contact subsystem
  + handler, both already fork-authored).
- **`LEDGER_quirks.md`:** (a) `-smoothNormal` gives a **C0** normal (no jump), **not** C1 (#4 Hermite is C1);
  (b) the smoothed normal uses the **facet** projection (`x̄` unchanged) — a documented second-order
  inconsistency; (c) auto-orientation **refuses** non-manifold / closed / multi-shell masters → faceted
  `-outward` fallback; (d) the global outward sign is a **surface** datum (kills R3), not per-slave; (e) the
  field is **frozen within the step** ⇒ the smoothed-NTS tangent is **modified-Newton (frozen `∂n/∂u`), NOT
  inconsistent-Newton** — the residual normal is the same frozen field, so residual⊥tangent are mutually
  consistent and the **converged answer is tangent-independent**; the full `∂n_smooth/∂u` (P3) changes only
  Newton *rate* on curved implicit masters, never the converged answer or the explicit lane; (f) **frictional
  NTS implicit is unsymmetric regardless of normal-tangent posture** (friction sliding tangent is non-symmetric,
  Abaqus §5.2.3) — do **not** cite "symmetric solver preserved" as a reason for the frozen-field normal tangent
  in the `-mu>0` case; (g) `-smoothNormal` + `-consistentNormal` together (pre-P3) **silently use the
  frozen-field tangent** (the smoothed-normal consistent tangent is P3) — a one-time warning surfaces it.
- **`47_…deferrals_adr.md`:** mark row **#4a** `pulled forward → ADR-63` (in the same PR).
- **`60_…adr.md`:** flip the **R3** row to RESOLVED (cite ADR-63); update "Exposed combo today" — the curved
  `-reemit + -mu` combo is closed by `-smoothNormal`.
- **Banner:** `shipped`-gated line in `Ladruno_scripts/banner_features.txt` when P1 lands, then
  `patch_banner.py`.
- **Header stamp:** run `Ladruno_scripts/stamp_headers.py` on the new `LadrunoContactNormalField.h` (add to
  the GLOBS) — non-optional for a new fork-authored file.

---

## Risks / open questions

- **Q-TANGENT — RESOLVED (§Gate decision):** symmetric frozen-field `kn·BᵀB` first; full `∂n_smooth/∂u` is a
  conditionally-required gated P3. The convex-ridge **correctness** fix is geometry/sign only (tangent-independent
  converged answer); the frozen-field tangent is an exact modified-Newton, symmetric (frictionless-implicit),
  explicit-OK, and matches the shipped default's posture. The curved-implicit Newton-convergence risk is gated
  by the P2 tripwire + the silent-downgrade warning rather than pre-empted by building B up front.
- **Q-PROJECTION-ANCHOR:** the MVP keeps the **facet** projection and smooths only the normal (second-order
  inconsistency). The Abaqus §5.1.1 **anchor-on-the-smooth-field** construction (project *along* `N(X)`) is
  the exact form — heavier (Newton on the smooth field), deferred unless the inconsistency shows up in the
  chatter gate.
- **Q-WEIGHTING:** area-weight (default) vs angle-weight vs unweighted nodal normals. Area-weight is the
  robust default (a sliver facet contributes little); angle-weight is the textbook crease-robust choice.
  Settle empirically at P1; expose `-normalWeight area|angle` only if the gate shows a difference.
- **Q-MULTISHELL:** the MVP refuses a disconnected multi-shell master without `-outward`. A
  per-component sign-vote (each connected component votes against its own nearest slaves) is a clean
  follow-up if a multi-shell curved master use case lands.
- **Q-RESTART:** the field + winding + global sign are **reconstructed at the next handle** (re-derivable from
  `mTags` + coords + the captured-once vote), so a database/`sendSelf` round-trip that drops them costs at most
  one conservative recompute (mirrors ADR-60 Q-RESTART). No serialization added.
- **Q-IMPLICIT-NEWTON:** does Newton converge with the frozen-field tangent on a strongly-curved master?
  P2/P3 gate; if not, the P3 full `∂n_smooth/∂u` is the remedy (or declare smoothed-NTS explicit-first, like
  ADR-39/60).

---

## Gate decision (RESOLVED — symmetric-first)

**Decision (2026-06-30, after a 2-iteration tangent-posture review loop — 3 independent lenses + a steelman-B
advocate + a coherence check, all converging):** #4a ships with the **symmetric frozen-field `kn·BᵀB`** tangent
first (P1/P2); the full **`∂n_smooth/∂u`** is a **conditionally-required, evidence-gated P3** (promoted to
required only if the P2 curved-implicit convergence tripwire trips). The decision was **stable across 2
iterations with no decision-flipping (BLOCKER/HIGH-that-refutes) finding** — the steelman for the alternative
**conceded** to symmetric-first.

**Why symmetric-first wins (rationale ordered by weight):**

1. **The primary lane is explicit — it assembles no tangent.** Under `CentralDifferenceLadruno` the contact FE
   returns the integrator-routed (zero) tangent and contributes a **residual only**; the full `∂n_smooth/∂u`
   would **never be evaluated** in the load-bearing lane. Building a wider non-symmetric stencil the primary
   lane cannot use is pure cost.
2. **The correctness content is tangent-independent.** The R3 ridge-flip fix (global surface sign, D2) and the
   chatter fix (C0 normal) are **residual / sign / active-set** effects; the **converged answer is identical**
   regardless of tangent. P1's convex-ridge + chatter gates pass or fail on the residual and the active set,
   neither of which the consistent tangent touches. Shipping the fix without the full tangent is the *whole*
   fix for the stated problem, not a half-fix.
3. **Smaller adversarial + assembly surface.** Full `∂n_smooth/∂u` widens the FE connectivity from ≤5 nodes to
   the ~10–17-node 1-ring, re-colors the graph on every active-set/segment change, and is non-symmetric ⇒
   unsymmetric solver — a global cost paid by *all* implicit users to serve the strongly-curved-implicit-Newton
   subset.
4. **Symmetry parity with the incumbent.** The frozen-field `kn·BᵀB` posture is **exactly the shipped faceted
   default's posture** (its B3 `∂n/∂u` block is itself opt-in via `-consistentNormal`); demanding
   consistency-from-day-one of the *new* feature while the *incumbent default* ships without it is an
   inconsistent bar. *(The symmetric-solver advantage itself is decisive only for the frictionless-implicit
   lane — under `-mu>0` the friction tangent is already non-symmetric, Abaqus §5.2.3.)*

**The conceded cost (folded as guards, not a flip):** the steelman's strongest point — a user on a curved
master setting `-smoothNormal -consistentNormal` gets the consistent tangent **silently downgraded** to
frozen-field, in the feature's own headline use case, with the "P3 vaporware" risk landing hardest *because*
R3 is about curved masters — is real but **HIGH, not BLOCKER** (no wrong answer, only slow/stalled Newton in a
deferrable, recoverable subset). It is mitigated **without** building B up front, via three folded guards: the
**P2 convergence tripwire** (binds P3-promotion to a *measured* condition, killing the discretionary-backlog
vaporware risk), the **evidence-gated P2 oracle** (refutes or confirms the concern with data), and the
**one-time `-smoothNormal`+`-consistentNormal` warning** (turns the silent downgrade into a documented limit).

---

## Implementation notes (P0+P1 built, 2026-06-30)

Built together in one local build (the phases are tightly coupled and the build cycle is ~20 min):
- **Files:** `LadrunoContactNormalField.h` (new, header-only — `propagateOrientation` flood-fill winding,
  `newellAreaNormal` area-weighted facet normal, `nodalNormals` with the global sign, `smoothNormal` query +
  refusal `Status`); `evalSegmentSmooth` in `LadrunoContactProjection.h` (facet projection + in-bounds gate
  unchanged, smoothed normal + degenerate-blend fallback to `normalOriented`); `LadrunoContactDomain.{h,cpp}`
  (per-master-surface field store, fingerprint-cached σ, `setNormalField`/`getSegNodalNorm`/`clearNormalFields`,
  `Contact.smoothNormal` default-init); `LadrunoContactFE.{h,cpp}` (`setSmoothNormals` + frozen `nodalNorm[4][3]`;
  `segmentActive` branches to `evalSegmentSmooth`; faceted B3 `∂n/∂u` suppressed under smoothing — D4);
  `LadrunoContactHandler.cpp` (per-contact field build from the deformed config + the once-captured global seed,
  refuse + silent-downgrade warnings, installs each segment's normals); `OpenSeesOutputCommands.cpp` (`-smoothNormal`
  parse + `-mortar` refusal). No new classTag; one vanilla file (the command parser, ledgered).
- **Bug caught on the convex-ridge gate (folded):** the auto global-sign seed `slave_centroid − master_centroid`
  used the **flat `mTags` connectivity** for the master centroid, which **double-counts shared/ridge nodes**
  and biased the centroid toward the high-valence ridge → on a 90° tent the biased master-centroid `z` ≈ the
  slave `z`, so the seed's `z` went slightly negative → the global sign `G` flipped → the whole field pointed
  INWARD → `gap > 0` ("not penetrating") → silent pass-through *that looked exactly like the R3 bug the feature
  fixes*. Fix: dedupe master tags (`std::set<int>`) before averaging → unique-node centroid (the ADR D2 intent).
  → `LEDGER_quirks.md`. The auto sign stays fragile when the slave cloud straddles the master centroid plane —
  `-outward` is the documented escape (BLOCKER-GLOBAL-SIGN already names it).
- **Gates green:** oracle 18/18; `tests/test_adr63_smoothnormal_p1.py` 3/3; ADR-39/41/57/60 battery (150) OFF
  byte-identical. Local build recipe = the ADR-60 one (junction MUMPS + `LADRUNO_OPENSEES_QUIET=1` +
  `-S` runner `Ladruno_scripts/_run_adr63_tests.py`).

## P1 adversarial review (4-lens, 2026-06-30 — findings + dispositions folded)

A scoped post-code review (4 lenses) on the built P1, targeting the novel/fragile surface the tests
under-cover. Two lenses came back CLEAN; two surfaced real issues, all fixed or documented below.

**Lens — frozen-field tangent consistency: CLEAN.** Verified `residual n̄ ≡ tangent n̄` (both derive `B`
from one pure `segmentActive`→`evalSegmentSmooth` reading the frozen `nodalNorm`); the tangent is the
exact modified-Newton of the frozen-geometry residual; the explicit lane assembles no tangent; the
faceted B3 suppression under smoothing is total. No correctness bug. The dropped `∂n̄/∂u` is the
documented, gated-P3 convergence-rate approximation.

**Lens — OFF byte-identity + shared-header blast radius: CLEAN.** The existing projection functions are
byte-verbatim; the new header leaks nothing (namespaced); the `Contact` struct grows one default-init
bool with no serialization/`sizeof` consumer; the mortar/edge loops are untouched (field build is
`ct.smoothNormal`-gated, NTS-only); the parser wiring is correct. Verified green: the full ADR-39/41/57/60
battery (147) OFF byte-identical.

**Lens — global-sign robustness (the demonstrated-fragile area): fixed.**
- **F1 (BLOCKER) — the sign was re-decided every handle from the DEFORMED vote, not captured once (D2
  violation): a deformable master under `-reemit` could flip the whole field mid-run.** FIXED — the
  engine captures the global sign ONCE (first OK build) into `NormalField.globalSign`/`signCaptured` and
  freezes it; `nodalNormals` now takes the frozen `G`. `voteSign()` decides it once. (F6 — ref-seed vs
  deformed-vote config mismatch — mooted by the freeze.)
- **F2/F3/F5 (HIGH/MED, silent) — `vote·seed ≈ 0` is a coin-flip (two-sided master, saddle, straddle);
  `-outward` fed the vote rather than overriding it.** FIXED — `voteSign` reports a **sign confidence**
  `|cos∠(vote,seed)|`; the handler warns ONCE (loud) when the auto sign is ill-conditioned (`conf < 0.1`)
  and recommends `-outward`. Converts silent-wrong → loud. (`-outward` given ⇒ the warning is suppressed;
  the user owns the direction.)
- **F4 (MED) — the ADR claimed closed surfaces are refused but the code only checked manifold-ness.**
  FIXED — `propagateOrientation` now counts boundary edges and returns `CLOSED` when there are none (a
  box/tube), refusing → faceted fallback. Oracle: a tetrahedron refuses.
- **F7 (MED, silent) — on refusal the fallback silently reverts to the R3-prone faceted path.** FIXED —
  the refuse warning now states the R3 protection is NOT active on that surface and to prefer `-outward`.

**Lens — coverage gaps: fixed.**
- **GAP-2 (HIGH, silent) — the degenerate-blend fallback used the per-slave auto `orientDir` (the
  R3-buggy datum).** FIXED — for a `-smoothNormal` contact the adapter's fallback `orientDir` is the
  captured GLOBAL seed (`smoothSeed`), so the degenerate-query fallback is also global-sign-correct.
- **GAP-4 (MED, silent) — a missing master node's backfill polluted the SHARED nodal normals (no
  per-pair skip, unlike the bucket grid).** FIXED — a missing master node refuses the field (→ faceted
  fallback), rather than backfilling a garbage facet normal.
- **GAP-5 (LOW) — the `1e-300` absolute degeneracy floor normalized a near-cancelling node to noise.**
  FIXED — a RELATIVE test (`|m_a| < 1e-9·max-incident-contribution`) leaves a cancelling node zero → the
  query falls back.
- **GAP-1 (coverage) — tri-3 masters were solver-unverified.** FIXED — a flat 2-triangle patch in-solver
  test exercises the full `nps=3` chain (`test_adr63_smoothnormal_p1.py`).

**NEW finding at bring-up (documented limitation) — spurious adjacent-facet activation at a sharp
ridge.** A slave sitting AT a sharp convex ridge has its closest-point land on the shared edge of the
*adjacent* facet (barycentric ≈ 0, marginally in-bounds), which then reads a large spurious penetration
and injects a big ejecting force. This is a **pre-existing NTS facet-ownership issue** (the faceted path
only *incidentally* prunes it — at a 90° ridge the neighbor's normal is ~⟂ the per-pair `orientDir`, so
the `normalOriented` fail-safe kills it). The smoothed path inherits it and has no such incidental prune;
a blunt interior-margin fix was tried and rejected (it removed a spuriously-*helpful* neighbor force and
destabilized other geometries). **Disposition:** ~~DEFERRED~~ **RESOLVED at P2.1 (2026-07-01)** — a
**gap-aware** shared-edge ownership guard (see "P2.1 implementation notes" below): reject a smoothed
projection on a SHARED interior edge only when its penetration is large vs the facet size, which drops the
spurious large-gap non-owner but keeps the true owner AND the at-apex contact. (The literal "reject on a
shared edge" was tried FIRST and reproduced exactly the reverted blunt fix — apex pass-through — confirming
the gap-magnitude gate is the load-bearing ingredient.) Full closest-facet / interior selection under
sliding remains a separate concern related to [[57_ladruno_edge_edge_contact_adr]] #4b edge-handoff.
→ `LEDGER_quirks.md`.

## P2.1 implementation notes (facet ownership at a sharp convex ridge, 2026-07-01)

**The one real P1 limitation, closed.** A per-segment shared-edge ownership guard, smoothed-path-only; OFF
and the faceted path stay byte-identical (battery **152**). Files: `LadrunoContactNormalField.h`
(`segmentSharedEdges` = the per-segment interior/shared-edge mask, reuses `propagateOrientation`'s
undirected-edge tally, TOPOLOGICAL ⇒ cached with σ; `onSharedInteriorEdge` = the parametric edge test whose
local-edge→(ξ,η)-boundary map matches `LadrunoContactProjection::shape()`'s node order); `…Projection.h`
(`evalSegmentSmooth` computes the gap, then if the projection is on a SHARED edge **and**
`−gap > edgeGapFrac·(|g1|+|g2|)` returns false); `…Domain.{h,cpp}` (`NormalField.sharedEdge` cached beside
σ; `getSegSharedEdge`); `LadrunoContactFE.{h,cpp}` (`sharedEdge[4]` member default all-zero ⇒ guard inert;
`setSmoothNormals(nn, se)`); `LadrunoContactHandler.cpp` (installs the mask). No new classTag; no vanilla
touch.

- **Why gap-aware, not a blunt reject.** A frictionless slave slides UP-slope to the ridge apex (the smoothed
  normal is more vertical than the facet normal, so a facet-perpendicular load tilts up-ridge). At the apex
  BOTH facets project onto the shared edge, so a blunt near-edge reject drops both ⇒ **pass-through** (it
  regressed the P1 ridge gate to min_d=−338 — the same reverted blunt fix). The spurious non-owner is
  distinguished by its penetration ∝ its LATERAL distance from the ridge, whereas the owner and a genuine
  at-apex contact read a SMALL gap; so the guard rejects only a LARGE-gap shared-edge projection. Constants
  `edgeGapFrac=0.05`, `edgeTol=0.1` are dimensionless (gap relative to facet size; parametric edge band) ⇒
  robust across ridge angle + model scale.

- **P1 R3 gate corrected (ill-posed rig).** `test_p1_smoothnormal_holds_the_ridge_facet` had been passing for
  the wrong reason — the spurious ejection (the very P1 bug) pinned `min_d`. A frictionless, laterally-free
  slave under a lateral load on a convex ridge has NO equilibrium (it slides up-slope and launches), so a
  laterally-free rig cannot test "held". The slave's x,y DOFs are now FIXED (only z free), isolating the R3
  **sign**: with the global-sign field the normal stays repulsive (min_d≈−1e-3, held); with the faceted auto
  sign it FLIPS (min_d≈−318, driven through). Both assertions hold; the R3 fix is intact.

- **New gates.** `tests/test_adr63_smoothnormal_p2.py` (quad + tri convex ridge, slave near the ridge with
  x,y fixed, pressed straight down: NO upward ejection, held). Oracle P2.1 checks: the shared-edge mask on the
  quad tent + tri pair; reject-large-gap-non-owner / keep-owner / keep-apex / keep-free-edge via
  `evalSegmentSmooth`.

- **Residual limitation (P2.3 note).** Near the apex, within the small-gap band, the non-owner is still kept
  (a harmless small force) — mild over-stiffness at the exact ridge, not a pass-through. Full closest-facet
  selection (one owner even while sliding across an interior edge) is deferred to ADR-57 #4b. This matters for
  P2.3 (sliding under `-reemit` over a curved master) and should be re-examined there.

## P2.2 + P2.3 notes (friction + re-emit over a curved master, 2026-07-01)

**No engine code — the composition already works.** `LadrunoContactFE::segmentActive` builds the friction
slip `gTvec` via `LadrunoFrictionKernel::tangentPart(drel, n, gTvec)` with the SAME `n` the gap operator
uses — the smoothed normal when `useSmoothNormal`. So friction is projected against `n_smooth` for free.
The parser refuses `-reemit`/`-smoothNormal` only with `-mortar`; they compose with each other and with
`-mu` on an NTS contact, and the handler fires the field-build (`ct.smoothNormal`) and re-emit
(`ct.enableReemit`) blocks independently off the SAME deformed feed (D7). Gates
(`tests/test_adr63_smoothnormal_p23.py`, 3/3):

- **P2.3 — curved crossing sustained.** A frictional block dragged across an 8-facet convex arc with
  `-reemit -smoothNormal -mu -outward` stays in contact across the crest (maxpen<0.05, reaches the far
  side); WITHOUT `-reemit` it slides off its reference candidates past the crest and the press drives it
  through. This closes the ADR-60 "exposed combo" (`-reemit + -mu` on a curved master): smoothing fixes the
  SIGN (a global datum, no per-slave flip at the ridge — R3), re-emit fixes the SEARCH (deformed re-pair).
- **P2.2 — flat-friction consistency.** On a flat master `n_smooth == n_facet`, so a frictional slide with
  `-smoothNormal` is byte-identical to the faceted one (D5, now for friction).

**Caveats folded (documented, not new bugs):**

- **Auto-sign conditioning is UNCHANGED.** The global outward vote is ill-conditioned when the slave cloud
  grazes the master edge-on (seed ~⟂ the field) — the pre-existing F2/F3/F5 warning. A single slave
  *starting to the side* of a curved arc votes a near-horizontal seed whose tiny z-component can flip the
  sign INWARD ⇒ pass-through even with `-smoothNormal`. So the ADR-60 R3 `-outward` caveat is **lifted only
  when the sign is well-conditioned** (slave cloud over the master) or `-outward` is given — it is NOT a
  blanket lift. The gates pass `-outward` for determinism.
- **Near-apex over-stiffness (the P2.1-guard interaction under sliding, re-examined).** On a SHARP crest the
  gap-aware guard keeps the small-gap non-owner near the apex ⇒ two facets briefly co-active with one
  smoothed normal ⇒ a mild extra upward bump as a block crests at speed (measured: it crosses slightly less
  far, never diverges). NEGLIGIBLE on realistic shallow arcs (maxpen 0.0096 frictionless vs 0.0122
  frictional). It is a *quality* effect, not a pass-through; the proper single-owner fix (closest-facet
  selection while sliding) is [[57_ladruno_edge_edge_contact_adr]] #4b, deferred. No code change this phase.
- **Traction continuity** (the C0-normal chatter fix) is a quality property; asserted here only qualitatively
  via sustained contact — a quantitative junction-traction gate is a follow-up (a clean measurement needs a
  quasi-static/implicit rig, which overlaps the Q-IMPLICIT-NEWTON tripwire below).
- **Q-IMPLICIT-NEWTON tripwire NOT yet run.** The curved-IMPLICIT convergence gate (frozen-field Newton on a
  genuinely curved master, gating P3-promotion) needs an implicit quasi-static rig and is a separate focused
  study — explicit point-mass-on-curve is too ballistic to double as the implicit convergence probe.
  *(RESOLVED at P2.4 below, 2026-07-01.)*

## P2.4 notes — Q-IMPLICIT-NEWTON convergence tripwire (2026-07-01)

**Outcome (a): CONVERGES — steelman-B's stall concern empirically REFUTED, P3 stays a deferred option.**
Test-only, NO engine code, no classTag. File: `tests/test_adr63_smoothnormal_p24_implicit.py` (3/3);
battery **160** OFF byte-identical.

**The rig (isolating `∂n_smooth/∂u` cleanly).** A symmetric convex bump of THREE quad facets: two steep
outer facets and a FLAT horizontal middle facet. All master nodes FIXED ⇒ the nodal-normal field `n_a` is
genuinely CONSTANT (a property of the fixed master geometry) ⇒ **no `-reemit` and no facet migration
needed** — which matters because there is NO validated static+reemit path (every ADR-60 reemit test is
explicit/CDL). The middle facet is flat, so its *faceted* normal is a constant vertical (0,0,1), but its
two shared-edge nodal normals inherit the steep outer facets' ~±22° tilt; so dragging a slave ACROSS the
flat middle facet rotates the *smoothed* normal `n_smooth(ξ)=ΣN_i(ξ)n_i` from vertical toward the tilt
while the facet geometry stays flat. The dropped `∂n_smooth/∂u` is therefore ENTIRELY the nodal-blend
rotation through the moving projection — exactly the P3 term, isolated from facet-projection curvature.

**Drive = StaticAnalysis + DisplacementControl on the slave x-DOF.** This is the only sanctioned
displacement driver — `LadrunoContactHandler::handle()` REFUSES a non-homogeneous (imposed-displacement)
SP and names DisplacementControl as the way ("it drives a DOF via the load factor, no imposed SP needed").
Implicit is also *required*: the explicit lane assembles NO tangent (the contact FE returns the
integrator-routed zero tangent), so only an implicit integrator's `addKtToTang` exercises the frozen-field
`kn·BᵀB`. The slave is seated at the flat-facet centre (x=0, where `n_smooth` is vertical by symmetry ⇒ no
frictionless seat-slide) under a constant press decoupled via `loadConst` (a single combined-load
DisplacementControl is DEGENERATE for a frictionless slave — the two equilibrium equations share one load
factor ⇒ equilibrium pins to a single slope), then dragged out. A weak lateral spring (ks≪kn) removes the
crest's zero-lateral-stiffness singularity at the seat step only (during the drag x is displacement-
controlled ⇒ non-singular anyway); its force is <0.5% of the contact force.

**Result.** The frozen-field symmetric `kn·BᵀB` Newton converges in **2 iterations every step**, and — the
load-bearing evidence — *independently of load-step coarseness*: the whole facet crossed in a SINGLE step
(maximal within-step normal rotation) still converges in 2 iterations. It stays 2 up to steep facets
(hz≤4), stiff penalty (kn≤1e6), heavy press (≤3000). **Why it never stalls:** the dropped `∂n_smooth/∂u`
term is `O(kn·gN)=O(press)` — scaled by the small penalty penetration `gN≈press/kn` ⇒ sub-dominant to the
kept `kn·BᵀB` in any well-posed penalty contact (the same reason the shipped faceted default drops its own
B3 by default and converges fine). The ONLY regime where iterations climb (a swept `gN≳15%` of the facet,
i.e. a soft-penalty misuse) is where the FACETED `-geomtan` consistent tangent ALSO diverges (and even the
seat fails) — a penalty-NTS breakdown, not a smoothed-normal defect, so P3's full `∂n_smooth/∂u` would not
rescue it either. → `LEDGER_quirks.md`.

**Rig lessons (folded).** (i) A lone frictionless slave on a convex master has zero lateral contact
stiffness at the crest (vertical n ⇒ `kn·nx²=0`) and no lateral equilibrium off-crest ⇒ the seat solve is
singular/runaway; the weak lateral spring + seat-at-the-symmetric-centre resolves it. (ii) The
combined-load DisplacementControl degeneracy (one λ scaling both press and drag) forces the `loadConst`
two-pattern split. (iii) No static+reemit coverage existed ⇒ the fixed-master constant-field rig sidesteps
it entirely. These are why the rig looks the way it does — recorded so the next agent doesn't re-derive them.

## Decision log

- 2026-06-30 — ADR opened as **#4a pulled forward** (ADR-47) to resolve **ADR-60 R3** (curved-master
  `orientDir` flip) + the NTS slice of ADR-41 Q-NORMAL. Next free ADR number was 63 (verified against
  `Ladruno_implementation/`; highest was 62).
- 2026-06-30 — Grounding: **Abaqus TG §5.1.1** averaged-nodal-normal `N(X)` is *verbatim* the construction
  (vs §5.1.2 Hermite C1 = the deferred #4); LS-DYNA §26.11 global-research + local-tracking architecture
  parallel (and the explicit note that LS-DYNA §17.1.4 is *mesh*, not contact, smoothing).
- 2026-06-30 — **5-lens adversarial gate (below); dispositions folded.** Key outcomes: the **global surface
  sign** (D2) replaces ADR-60's per-slave persisted orientDir as the *complete* R3 fix; **winding by
  flood-fill** with non-manifold/closed/multi-shell **refusal + faceted fallback** (BLOCKER-WINDING/FALLBACK);
  **facet projection retained** so R5 slide-off is untouched (BLOCKER-PROJECTION-CONSISTENCY); **frozen-field
  symmetric tangent** recommended, full `∂n/∂u` deferred to a gated P3 (BLOCKER-TANGENT-POSTURE);
  scope **NTS-only** with the header built reusable for a future mortar pass.
- 2026-06-30 — **Q-TANGENT RESOLVED via a dedicated 2-iteration tangent-posture review loop** (3 independent
  lenses [numerics/Newton, explicit/symmetry, Abaqus-precedent/cost] all favored symmetric-first at ~0.85; a
  steelman-B advocate **conceded** to A; a coherence check returned GO). Decision: **symmetric frozen-field
  `kn·BᵀB` first; full `∂n_smooth/∂u` is a conditionally-required gated P3.** Position stable across 2
  iterations, no decision-flipping finding. Dispositions folded: parity claim softened (dropped term is
  same-*kind*-but-larger-magnitude than the faceted ≡0 facet block); symmetric advantage scoped to
  **frictionless**-implicit (friction already non-symmetric, Abaqus §5.2.3) with explicit-lane irrelevance +
  smaller adversarial surface promoted to the load-bearing rationale; **frozen-field invariant** asserted
  (residual-`n̄` ≡ tangent-`n̄`); **BLOCKER-FALLBACK** extended to the **query-point blend** norm; **P3
  explicit prohibition** + **P2 curved-implicit convergence tripwire** + **`-smoothNormal`+`-consistentNormal`
  silent-downgrade warning** added; stencil cost stated precisely (≤5 → ~10–17 nodes, ~4–12× block area).

## Adversarial review log (5-lens, 2026-06-30 — findings + dispositions)

**Lens A — winding / manifold orientation (the hard part).** (A1, BLOCKER) a declared `MASTER_SEGMENTS`
surface can have **mixed segment winding**; averaging facet normals *without* coherent orientation cancels or
garbles `n_a` → **adopted: flood-fill `σ_s` from edge-adjacency (D3), BLOCKER-WINDING.** (A2, HIGH) a
**non-manifold edge** (>2 segments) makes orientation propagation ill-defined → **adopted: refuse + faceted
fallback (named error), do NOT guess.** (A3, HIGH) a **closed/self-enclosing** master (a box, a tube wrapping
on itself) has no single global `-outward` → **adopted: refuse closed surfaces; the supported target is an
open manifold patch (the R3 case).** (A4, med) a **disconnected multi-shell** has per-component sign ambiguity
→ **adopted: refuse without `-outward`; per-component vote is a Q-MULTISHELL follow-up.** (A5, med) seed
choice affects nothing topological but the **global sign** must not ride the seed's arbitrary geometric
normal → **adopted: global sign from `-outward` or the slave-centroid vote (D2/BLOCKER-GLOBAL-SIGN), not the
seed.**

**Lens B — byte-identity OFF + flag plumbing.** (B1) reuse the `if(ct.smoothNormal)` gate so OFF builds no
field, allocates nothing, and runs `evalSegment` verbatim → adopted (BLOCKER-IDENTITY). (B2, HIGH) the new
`smoothNormal` struct field is POD; uninitialized stack garbage could flip it true → **adopted: default member
initializer `= false`** (mirrors ADR-60 Finding-B4). (B3) flat-master ON is **not** byte-identical to OFF
(interp path differs) but **is** numerically identical → **adopted: assert `‖n_smooth−n_facet‖<1e-14`
consistency, not byte-identity (D5).** (B4) the field build must sit **strictly under** the flag and be keyed
by surface tag so a mixed model (one smoothed contact, one faceted) is correct per-contact → adopted.

**Lens C — projection / gap / active-set consistency.** (C1, HIGH) the closest-point projection makes
`(xs−x̄) ∥` the **facet** normal, not `n_smooth`; naively re-projecting along `n_smooth` is the heavy §5.1.1
anchor construction → **adopted: keep the facet projection, smooth only the normal; document the second-order
inconsistency (BLOCKER-PROJECTION-CONSISTENCY / Q-PROJECTION-ANCHOR).** (C2, HIGH) does smoothing perturb the
**active set** and thus [[60_…]] R5 slide-off? **No** — the in-bounds + penetration gate uses the **facet**
projection (unchanged); smoothing changes only the force *direction* → **R5 regression must stay green
(P1 oracle (e)).** (C3, med) `gap = n_smooth·(xs−x̄)` can differ slightly in sign near a grazing crossing from
the facet gap → bounded by the inter-segment angle; the penetration gate still uses the facet `evalSegment`
status → benign, documented. (C4) the `B`-operator must use the **same** `n_smooth` in the residual and the
friction slip `gTvec` projection → adopted (one normal per eval, threaded through `segmentActive`).

**Lens D — sign correctness (the R3 fix itself).** (D1, BLOCKER) if the outward sign stays **per-slave** the
ridge flip survives the smoothing → **adopted: global surface sign captured once (D2/BLOCKER-GLOBAL-SIGN);
this is THE fix.** (D2, HIGH) a **degenerate nodal normal** (`‖Σ w_s σ_s n_s‖→0` at a sharp fold, or a
node owning one sliver facet) yields a zero/garbage normal → **adopted: BLOCKER-FALLBACK to faceted
`normalOriented()` + one-time warning; never emit garbage.** (D3, med) the existing `normalOriented`
**fail-safe** (refuse when `refDir ⟂ n`, `LadrunoContactProjection.h:149-155`) must be preserved on the
fallback path → adopted (fallback *is* that code). (D4) the global-sign **vote** must use the **reference**
slave centroid (config-independent), not a live position, else a deformed slave cloud could vote the sign
mid-analysis → **adopted: reference centroid, captured once.**

**Lens E — composition + scope + completeness.** (E1) compose with ADR-60: same deformed feed, **same R1
fingerprint** invalidates the field (a re-mesh re-checks manifold-ness) → adopted (D7). (E2) the field is the
**complete** replacement for ADR-60 BLOCKER-ORIENTDIR persistence (global ⊃ per-slave) → adopted; ADR-60 R3
row resolves here. (E3, med) **tangent**: frozen-field `kn·BᵀB` is symmetric; the full `∂n_smooth/∂u` widens
to a non-symmetric 1-ring stencil → **adopted: symmetric-first (D4/BLOCKER-TANGENT-POSTURE), full term gated
P3; escalated to the dedicated tangent-posture loop below.** (E4, major) **mortar/edge-edge**: ADR-41 Q-NORMAL
is filed under the *mortar* GP normal, but the load-bearing *correctness* hole is the NTS R3 → **adopted: scope
NTS-only, build the header reusable, fence mortar/edge-edge to a future pass (D6/§Fence).** (E5) **recorder
continuity**: `getNtsForce` per-pair snapshots change direction with the smoothed normal — disclose in the P3
doc (summed force is the continuous quantity), mirrors ADR-60 P3. (E6) **restart**: field reconstructed at next
handle → Q-RESTART, no serialization (mirrors ADR-60). (E7) **C0 vs C1 honesty**: the gate must not
over-claim — #4a is a continuous (C0) normal, not the C1 Hermite surface; the chatter gate asserts
*continuity*, not *smoothness* of slope → **adopted: explicit C0 claim in §What + the P1 (b) tol wording.**

### Tangent-posture review loop (Q-TANGENT, 2-iteration convergence, 2026-06-30)

A focused adversarial loop on the one open gate decision (symmetric frozen-field `kn·BᵀB` first vs full
`∂n_smooth/∂u` up front). **Converged: symmetric-first, stable across 2 iterations, no decision-flipping
finding.**

**Iteration 1 (3 independent lenses, all favored A ~0.85):**
- **T-A (numerics / Newton, vs the shipped B3 path).** Verified the frozen-field tangent is an **exact
  modified-Newton on per-step-frozen geometry** (`∂r/∂u = kn·B(n̄)ᵀB(n̄)` with `n̄` the literal constant),
  **not** an inconsistent-Newton — converged answer tangent-independent. (T-A1, MED) the parity claim was
  **overstated**: the dropped `∂n_smooth/∂u` is the same *kind* but **larger magnitude** than the faceted
  default's (tri-3 facet block `≡0`), so curved-master stall is more probable → **folded** (D4, P3 likelihood).
  (T-A2, MED) **BLOCKER-FALLBACK** guarded only the *nodal-average* norm, not the *query-point blend* norm →
  **folded** (BLOCKER-FALLBACK (b)). (T-A3) frozen-field **invariant** (residual-`n̄` ≡ tangent-`n̄`) must be a
  P0/P1 assertion → **folded** (BLOCKER-TANGENT-POSTURE).
- **T-B (explicit-dynamics / symmetry / solver).** (T-B1) the **primary lane is explicit** — the contact FE
  assembles **zero** tangent (the P2a `getTangent`-routes-through-integrator lesson), so the full term is never
  evaluated in the load-bearing lane → **promoted to the lead rationale.** (T-B2, LOW) the symmetric-solver
  advantage is **moot under `-mu>0`** (friction tangent already non-symmetric, Abaqus §5.2.3) → **folded:
  symmetric advantage scoped frictionless-implicit; explicit prohibition on the P3 full term added.**
- **T-C (Abaqus §5.1.1 frozen-affine precedent / stencil cost).** The frozen-field-refreshed-each-handle
  posture is **more faithful to finite sliding than §5.1.1** (which freezes for the whole analysis); the only
  omitted §5.1.2 ingredient is the curvature second-variation = the P3 term. **Stencil verified vs source:**
  `normalGeomTangent` couples ≤5 nodes / 15 DOF; the full term couples the **1-ring** (~10–17 nodes, ~4–12×
  block area, re-colored on active-set change) → **folded** (D4 precise cost). No correctness dependence on the
  tangent found (R3 + chatter are residual/sign/active-set).

**Iteration 2 (steelman-B + coherence check — confirmed stability):**
- **T-D (steelman for the alternative).** Strongest case for building B up front: a user on a curved master
  with `-smoothNormal -consistentNormal` gets the consistent tangent **silently downgraded** to frozen-field —
  in the feature's headline use case — plus the **"P3 vaporware"** risk (R3 *is* curved masters, and P3 is the
  tangent that makes them converge). **Severity HIGH, not BLOCKER** (no wrong answer, only slow/stalled Newton
  in a deferrable, recoverable subset). **Verdict: CONCEDE to A** — the correctness fix is tangent-independent,
  the incumbent default already ships frozen-normal, and B forces a global non-symmetric/1-ring cost to serve a
  narrow subset. **Mitigations folded (not a flip):** P2 convergence **tripwire** binding P3-promotion to a
  *measured* condition (kills the vaporware risk), an **evidence-gated** P2 oracle, and a **one-time
  `-smoothNormal`+`-consistentNormal` warning** (silent downgrade → documented limit).
- **T-E (coherence check).** GO to fold, with reconciliations applied: "P3 more likely" / "promote to required"
  worded as a *probabilistic* / *conditional tripwire* (P3 row keeps `optional, gated` + a back-ref); the
  §Gate-decision "symmetric" rationale **edited in place** (not appended) so symmetry isn't claimed at two
  scopes; D-ii/D-viii folded as **refinements** of existing rows (no ledger duplication). No
  contradictions, no smuggled code; net-new = query-point fallback, P3 explicit prohibition, P2 tripwire.
