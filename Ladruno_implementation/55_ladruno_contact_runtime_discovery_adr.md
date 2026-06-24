---
title: Contact runtime discovery (pillar-2) — the open-universe broad phase + self-contact + edge-edge
project: Ladruno
status: draft
priority: medium
owner: nmora
tags:
  - implementation
  - subsystem
  - contact
  - broad-phase
  - self-contact
  - collapse
  - explicit
---

# Contact runtime discovery (pillar-2)

> **ADR 55.** Promotes the **runtime contact-pair discovery** subsystem out of the
> [[47_ladruno_contact_deferrals_adr]] deferral ledger into a **scoped design ADR**,
> **billed to the contact lane** ([[48_ladruno_contact_capstone_adr]]). Created
> 2026-06-24 as the load-bearing conclusion of the E-FEM/FDEM-lite adversarial review:
> **this — not the fracture front-ends — is the program's long pole.**
>
> **Why it's its own ADR, on the contact lane:** it has **independent, named users
> that need no crack at all** — pounding between adjacent buildings, foundation
> rocking/uplift re-contact, sliding/impact — i.e. the same dense-urban-retrofit
> research the contact lane already serves. Building it here de-risks **every**
> collapse route ([[51_ladruno_element_removal_adr]] recontact,
> [[54_ladruno_fdem_lite_adr]], [[53_ladruno_embedded_discontinuity_adr]]) **for
> free**, and is justified on its own merits regardless of whether continuum fracture
> is ever built.

## What

The subsystem that turns the fork's **predefined-pair, reference-config** contact
into **open-universe, current-config** contact: at runtime, **discover which body /
facet pairs are near**, among bodies whose adjacency is **not known a priori** and
**changes every step** — including **self-contact** and **edge-edge**.

> `[GATE — review finding]` **This is the genuinely hard ~60% of contact, and it is
> ~0% built.** The fork owns the **narrow phase** (NTS + mortar + friction + tie +
> viscous, ADR-39/41) + an **intra-pair culler** (`LadrunoContactBucketSort`,
> reference-config). The contact literature treats the narrow phase as the
> *additive, well-understood* layer that discovery is built *on top of*. So "we own
> the expensive 60%" (asserted by the fracture ADRs) **inverts the difficulty**:
> owned = the tractable ~40%; **this ADR is the hard ~60%.**

The six subsystems (per `opensees-contact/broad_phase_collision.md §6`), plus the two
structural gaps the review surfaced:

1. **Current-config geometric proxies** — per-body/facet AABB from *current* coords,
   **refit each step**, inflated by a Verlet skin so the candidate set survives a few
   steps. (Today's bucket sort is **reference-config** — a `domainChanged`-built grid.)
2. **A re-emitting broad phase** — rebuilt from current coords every step/epoch
   (LS-DYNA re-sorts every 10–15 cycles; the bucket-sort header calls this "the
   follow-on rung").
3. **A size-dispersity broad phase** — **BVH / loose octree / hierarchical hash**.
   `[GATE]` This **REPLACES** the uniform bucket sort, not extends it: the bucket
   sort sizes cells to the *median segment*, and fragmentation/debris **creates**
   size dispersity (slab → chips), which a single cell size handles *fatally* badly
   (`broad_phase_collision.md §3`).
4. **A candidate-pair manager** — add/persist/remove a hash-set of near pairs, with
   **friction-/λ-history carry across re-pairing** (else `LadrunoFrictionKernel` slip
   history is lost) and self-adjacency exclusion.
5. **An edge-edge narrow phase** `[GATE]` — the shipped narrow phase is
   node-to-segment / facet-mortar **only**; edge-edge is **irreducible** (NOT a
   vertex-face special case — `§5`) and tumbling polyhedral fragments hit edge-first
   constantly. A **new narrow-phase kernel**, currently invisible in the program.
6. **A runtime adapter pool + runtime surface registration** `[GATE]` — the deepest
   gap: `addSurface`/`addContact` are **parse-time only**; `LadrunoContactDomain`
   state is keyed by pre-declared `contactTag` and has **no free-body/fragment
   struct** (ADR-51's finding, verbatim). Creating/destroying `LadrunoContactFE`
   adapters and *registering new surfaces mid-analysis* is an unbuilt subsystem.

**Self-contact** `[GATE]` is treated as a **distinct sub-problem**, not a bullet: a
folding/piling body must be tested against *itself*, which the master/slave
predefined-pair model **cannot express**. It needs a **single-surface formulation**
fed by a BVH/hash, with **topological-adjacency exclusion** and a coplanar angle
filter. It is the *typical* debris case, and is currently fenced to
[[47_ladruno_contact_deferrals_adr]].

**NOT in this ADR:** the fracture front-ends (53/54) and element removal (51) — this
is the *contact* subsystem they all consume. MPI contact decomposition stays deferred
([[47_ladruno_contact_deferrals_adr]]).

## Why (independent of fracture)

- **Named contact-lane users need it with no crack:** pounding (adjacent buildings,
  dense-urban retrofit), foundation **rocking/uplift re-contact**, sliding bearings
  post-uplift, member-to-member post-buckling contact in frames. All require
  **runtime discovery of pairs whose contact is not predefined** — exactly this
  subsystem. The contact capstone (48) explicitly parks these.
- **It is the shared long pole for the whole collapse program.** Removal-recontact
  (51 pillar-2), FDEM-lite (54), and E-FEM (53) **all** terminate in the same
  discrete stage and **all** gate on this. Build it **once, here**, and the rest
  inherit it.
- **It plays to the fork's edge.** The fork's genuine contact head-start (NTS +
  mortar + viscous) is the *narrow phase*; this completes it into a real contact
  engine the way LS-DYNA/Abaqus/Code_Aster have, rather than re-deriving the
  constitutive layer that's already done.

## Where

- **`SRC/domain/contact/LadrunoContactBucketSort.h`** — to be **superseded** by a
  current-config, re-emitting BVH/octree/hash for size-dispersity (keep the bucket
  sort for the static/predefined small-motion path).
- **`SRC/domain/contact/LadrunoContactDomain.{h,cpp}`** — needs a **free-body /
  fragment representation** + a **runtime `addSurface`/`addContact` path** (today
  parse-time only) + the candidate-pair manager + GC that tolerates topology change.
- **`SRC/analysis/handler/LadrunoContactFE.*`** + `LadrunoContactHandler` — a
  **runtime adapter pool** (create/destroy adapters as pairs appear/separate), within
  the forced-`domainChanged`/re-handle model (ADR-51).
- **New:** the edge-edge narrow-phase kernel; the single-surface self-contact path.
- Class tags / banner at implementation; [[LEDGER_vanilla_files]] for any Domain seam.

## Risks / open questions

> [!question] Q-SCALE — this is the largest single item in the collapse program
> Six subsystems + a broad-phase **rewrite** + a new narrow-phase kernel + a
> different-class self-contact formulation. Scope and resource it as the **dominant
> cost**, ahead of any fracture front-end — not as a "P3 line item."

> [!question] Q-EXPLICIT-CHURN — the active set churns every step
> Open-universe discovery + per-step re-pairing argues for **explicit penalty**
> (forces to RHS, no re-factorization) — matching the contact lane's explicit-first
> decision (ADR-39). Implicit would need an open-close active-set iteration (DDA-style)
> around Newton; far more invasive.

> [!question] Q-SELF-CONTACT-DATA-MODEL — the master/slave model can't express it
> Self-contact needs single-surface + adjacency exclusion; how much of
> `LadrunoContactDomain`'s pair-keyed state must change to host it?

> [!question] Q-STATE-REKEYING — committed pair-state under topology change
> `FrictionState`/`MortarNormalState` are keyed by surface-relative ordinals;
> runtime surface creation/renumbering must not alias history onto the wrong pair
> (the ADR-51 review's silent-wrong-friction bug, at scale).

> [!question] Q-MPI — deferred, but don't design it out
> MPI contact decomposition stays in [[47_ladruno_contact_deferrals_adr]]; the data
> structures here should not assume a single `Domain` (mirror ADR-51's MP posture).

## Phasing

| Phase | Scope | Gate |
|---|---|---|
| **P1** | Current-config proxies + re-emitting broad phase (BVH/hash) over **existing predefined surfaces**; candidate-pair manager with history carry | Large-sliding contact (the bucket-sort "follow-on rung") works without analyst re-analyze; friction history survives re-pairing |
| **P2** | Runtime surface/adapter registration + free-body struct in `LadrunoContactDomain` | A body added/removed mid-analysis gets contact with no leak; pair-state not aliased |
| **P3** | Edge-edge narrow-phase kernel | Tumbling polyhedra contact edge-first correctly |
| **P4** | Single-surface **self-contact** (BVH/hash + adjacency exclusion) | A folding/piling body self-contacts without spurious adjacent-facet forces |
| → consumers | Removal-recontact (51), FDEM-lite (54), E-FEM (53) all consume P1–P4 | — |

## See also
- [[48_ladruno_contact_capstone_adr]] · [[39_ladruno_contact_domain_adr]] ·
  [[41_ladruno_mortar_alm_contact_adr]] — the contact lane this belongs to.
- [[47_ladruno_contact_deferrals_adr]] — where self-contact / discovery were fenced
  (this ADR promotes them out).
- Consumers: [[51_ladruno_element_removal_adr]] (recontact) ·
  [[54_ladruno_fdem_lite_adr]] · [[53_ladruno_embedded_discontinuity_adr]].
- Skill: `opensees-contact` (`broad_phase_collision.md` — the six-subsystem spec).
- Primary: Munjiza & Andrews (1998, NBS); Munjiza (2004); Ericson (2005, *Real-Time
  Collision Detection* — BVH/SAP/hashing); Teschner et al. (2003, spatial hashing);
  Shi (1992, DDA open-close).
