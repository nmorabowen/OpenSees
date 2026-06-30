# ADR-60a — PFEM Cross-Pollination Amendment (finite-sliding contact re-emission)

- **Status:** AMENDMENT to [[60_ladruno_finite_sliding_reemission_adr]]. Strategy + grounding only — no
  new scope, no new classTag. Folds the OpenSees PFEM remesh study into ADR-60's references and P4 roadmap.
- **Owner:** Mora Bowen · Palacios · Abell · Guppi
- **Date:** 2026-06-30
- **Parents:** [[60_ladruno_finite_sliding_reemission_adr]] (the host), [[39_ladruno_contact_domain_adr]]
  (the engine), [[55_ladruno_contact_runtime_discovery_adr]] (open-universe successor).
- **Prior art (studied):** OpenSees PFEM subsystem — `SRC/element/PFEMElement/` (mesher + BackgroundMesh),
  `SRC/analysis/analysis/PFEMAnalysis.{h,cpp}`, `OTHER/Triangle` (Shewchuk) + `OTHER/Tetgen` (Si).

---

## Why this amendment

We scoped the PFEM (Particle Finite Element Method) remesh machinery to test a hypothesis: **contact
re-emission is PFEM remesh restricted to the interface.** The hypothesis holds. PFEM rebuilds a Lagrangian
mesh every step by re-triangulating a moving particle cloud (Delaunay + alpha-shape) — the same
"rebuild connectivity from the current deformed configuration" move ADR-60 makes for the NTS broad phase.
The study delivers three things to ADR-60: it **independently validates** two core design choices, it
supplies a **concrete reference implementation** for the re-attach loop and a cost-optimization roadmap, and
it sharply isolates the **one thing PFEM cannot help with** — friction history transfer — as the load-bearing
novelty that is genuinely ours to design.

This is grounding, not new scope. ADR-60's phases, BLOCKERs, and fences are unchanged; this amendment adds
references and a P4 reference implementation, and reinforces where the adversarial scrutiny must concentrate.

---

## The paradigm map (why it transfers)

PFEM's defining decomposition: **nodes (particles) are persistent and carry all state; elements are
disposable connectivity, rebuilt every step.** This is *already* the ADR-39/60 architecture — the engine
(`LadrunoContactDomain`) owns persistent path state keyed to entities, and the assembly-side
`LadrunoContactFE` adapters are stateless views rebuilt every `handle()` (`LadrunoContactDomain.h:22-26`).

| PFEM | Ladruno contact | State carried? |
|---|---|---|
| Background bucket grid (`BackgroundGrid`/`BCell`/`BNode`) | `LadrunoContactBucketSort::Grid` | none — pure geometry, rebuild free |
| Delaunay/clip → which triangle a node joins | NTS narrow phase → which segment a slave pairs with | transient pairing; flips under sliding |
| Persistent fluid `Node` (vel, pressure) | persistent slave `Node` + `FrictionState`/`MortarNormalState` | yes — entity-keyed, survives the rebuild |
| `Pressure_Constraint` re-link after remesh | `frictionGC*` + `getOrCreateFrictionState(c,s,segIndex)` re-attach | the re-attach contract |

The reason this decomposition makes contact the *safe* re-emit target (vs. solid remesh, which hits the
Gauss-point state-transfer wall) is that almost all contact state is **node/segment-keyed, not
Gauss-point-keyed** — it survives a connectivity rebuild by construction. The lone exception is friction
(see "the gap PFEM does not close").

---

## What PFEM independently validates

1. **Commit-boundary cadence — keep it.** PFEM remeshes *between* steps, never inside the Newton/equilibrium
   iteration: `PFEMAnalysis::analyze()` re-handles at the top of the next step and `identify()` runs at the
   commit boundary, riding the `domainChange()` flag the standard drivers poll. Remeshing mid-iteration would
   destroy Jacobian continuity and chatter. **ADR-60's trigger placement — in `Domain::commit()` →
   `domainChange()`, consumed by the re-handle on the next step (`60_…adr.md` §How, `Domain.cpp:2201-2208`) —
   is the identical discipline.** PFEM is external corroboration that the commit-boundary is the correct (only)
   safe re-emit point. This directly answers Q-IMPLICIT-NEWTON's worry: the fix is not "make Newton converge
   through a re-pair mid-iteration," it is "never re-pair mid-iteration" — exactly what we do.

2. **Re-attach by stable entity key.** PFEM's `identify()` is: disconnect every `Pressure_Constraint` → walk
   the new element connectivity → reconnect by **node key**. ADR-60's `frictionGCBegin/Mark/End` +
   `getOrCreateFrictionState(contactTag, slaveTag, segIndex)` (`LadrunoContactDomain.h:405-412`) is the same
   disconnect / re-walk / reattach-by-stable-key pattern. PFEM's `identify()` is therefore the **reference
   implementation** for our re-attach loop — including its **isolated-node handling** (particles that fall
   outside the remeshed boundary are driven by gravity alone), which is the direct analog of
   `BLOCKER-SLIDE-OFF` (a slave that loses all candidate segments → drop the friction slot, force → 0).

---

## Concrete reference: the PFEM file map

For the P3 doc deliverable and any future implementer, the PFEM sources to read alongside ADR-60:

- `SRC/analysis/analysis/PFEMAnalysis.cpp` — `analyze()` (commit-boundary re-handle loop, adaptive-dt) and
  `identify()` (disconnect/reconnect-by-node-key; isolated-particle handling). **The re-attach reference.**
- `SRC/element/PFEMElement/PFEMMesher2D.cpp` — the classic full-remesh: collect deformed coords → triangulate
  → alpha-shape filter → delete old / create new → reconnect constraints. The "rebuild from current config"
  template (but per-step + global, which we deliberately do not copy — see fences).
- `SRC/element/PFEMElement/BackgroundMesh.cpp` (`remesh()`) + `BCell`/`BNode` — the **fixed-grid** alternative
  that re-buckets moving particles **without** a full re-triangulation. The P4 reference (below).

---

## The cost path: BackgroundMesh ⇒ the P4 patch adapter

ADR-60 D6 is explicit that baseline re-emit rides a **full `domainChanged`** (renumber + graph + SOE rebuild)
at frequency ≈ 2D/h — the deliberate, cheap-to-implement subset, default OFF. PFEM has *two* meshing
architectures, and the second is precisely the optimization ADR-60 already names as P4:

- **Classic full-remesh** (`PFEMMesher2D`): re-triangulate the whole cloud every step. ≈ ADR-60 baseline P1
  (full `domainChanged` on every re-emit).
- **BackgroundMesh** (`BackgroundMesh::remesh`): a fixed cell grid; only particles that crossed a cell are
  re-bucketed; no global re-triangulation. ≈ ADR-60 **P4 patch / local-tracking** — each slave keeps a
  neighborhood patch of segments; within-patch sliding is narrow-phase-only (no `domainChanged`); re-emit
  fires only on patch-exit (`60_…adr.md` §"The patch / local-tracking alternative (P4)").

**Action for P4:** cite `BackgroundMesh::remesh` + the `BCell`/`BNode` re-bucket as the working reference for
"re-bucket without a full rebuild." It is the same data structure as `LadrunoContactBucketSort` and the same
move, already proven in-tree under a transient analysis. This de-risks P4's connectivity refactor.

---

## The gap PFEM does NOT close (the load-bearing novelty)

**PFEM fluids carry no tangential history.** Their stress is a pure function of the *current* strain rate and
pressure (incompressible Navier–Stokes with bubble stabilization), so a PFEM remesh is genuinely state-free —
and offers **no template** for the single hard part of contact re-emission: continuing stick/slip friction
history when a slave **crosses a segment boundary** under sustained sliding.

ADR-60 D4 (re-engage on a fresh slot, do *not* transfer) is correct for the crossing **instant** — the fresh
slot self-captures `gT0 = current gTvec` ⇒ `gTeff − gpT = 0` ⇒ zero stick force ⇒ traction-continuous at the
crossing (`60_…adr.md` D4; `LadrunoContactDomain.h:228-238`). But **sustained-drag continuity across a
crossing** — a block sliding under steady friction whose contact point walks from facet A to facet B without
the traction blinking — is P4's job, and the frame transform it needs has **no PFEM analog**:

$$g_T^{(B)} = R_{A\to B}\,g_T^{(A)}, \qquad \text{re-evaluate } \lVert\tau\rVert \lessgtr \mu p_N \text{ in B's local frame.}$$

Get this wrong and every facet crossing injects a spurious tangential impulse on what should be smooth
sliding. This is the place to concentrate adversarial scrutiny in the P2/P4 gates: **PFEM provides no safety
net here.** (Mortar buys slack — its `λ_N` is scalar/node-keyed/frame-free and survives finite sliding for
free; the mortar `gT0` staleness on dominant-facet change is the pre-existing, out-of-scope Q-MORTAR-GT0.)

---

## Fences / non-goals (what we are NOT importing from PFEM)

- **No runtime gmsh/Delaunay re-triangulation of the contact surfaces.** Surfaces are *declared* (closed
  universe). Runtime mesh generation of new bodies/surfaces is ADR-55's open universe, not ADR-60.
- **No alpha-shape.** PFEM detects the free surface by which triangles survive a circumradius test; the
  contact active set is **gap-driven** (`gap < 0`), not circumradius-driven. The analogy is real but weak and
  is not lifted.
- **apeGmsh stays a SEEDER, not a runtime remesher.** If apeGmsh ever touches this, its role is to *generate*
  the declared surfaces / contact pairs / initial configuration (and, separately, to seed an actual OpenSees
  PFEM fluid model — particle cloud + `Pressure_Constraint`s + the `analysis PFEM` stack). Runtime re-emission
  stays in C++, behind the `FE_Element`-only-in-`handle()` wall ADR-60 §How documents. This amendment does not
  propose any apeGmsh code.

---

## Actions (fold into ADR-60 phases)

- **P3 (doc):** add the PFEM file map above to ADR-60 §Where "Reference sources," and a `LEDGER_quirks` note:
  *PFEM `identify()` is the re-attach reference; `BackgroundMesh::remesh` is the P4 re-bucket reference; PFEM
  validates the commit-boundary cadence but offers no friction-transfer template.*
- **P4:** adopt `BackgroundMesh::remesh` + `BCell`/`BNode` as the named reference for within-patch
  re-bucket-without-rebuild; this is the concrete prior art for the local-tracking adapter.
- **P2/P4 gates:** flag the friction frame-transform as the highest-scrutiny item — it is the one piece with
  no in-tree or PFEM precedent.
