---
title: Ladruno Contact — deferred-features ledger (the ADR-47 hard-defer set)
project: Ladruno
status: draft
priority: low
owner: nmora
tags:
  - implementation
  - contact
  - deferral-ledger
---

# Ladruno Contact — deferred-features ledger (ADR-47)

> Stub created 2026-06-23 to give the contact subsystem's **hard-deferred** items a documented home
> (referenced by [[39_ladruno_contact_domain_adr]], [[41_ladruno_mortar_alm_contact_adr]], and the
> [[48_ladruno_contact_capstone_adr]] capstone). This is a **ledger, not a design** — each row records
> *why* the item is out of the committed P0→P4 / A→C scope and the **concrete trigger** that would
> re-open it into its own full ADR. Created before Track C1 lands, per the capstone's ADR-47 risk
> (C1/C2 gates lean on the dual-basis and faceted-normal deferral rationale).

## The deferred set

| # | Deferred item | Why deferred (rejection rationale) | Re-open trigger | Source |
|---|---|---|---|---|
| 1 | **Dual / biorthogonal mortar basis** (diagonal `D`, cheap nodal `λ`, LBB/inf-sup-optimal) | Standard basis + **overlap clipping** + finite-`epsN` ALM already passes the non-matched **patch test** (the MVP bar); the dual basis is an optimization, not a correctness need. | C1 patch-test gate reports residual interface-pressure oscillation that finite-`epsN` ALM cannot suppress, **or** an inf-sup/LBB-stable formulation is genuinely required. | ADR-41 Q-MORTARLITE |
| 2 | **True-LM / saddle-point enforcement + inf-sup stabilization** | Zero-diagonal pathology + active-set churn forces a full `domainChanged()` re-number; commit-cycle **Uzawa-over-penalty** gives the exact constraint at finite penalty with **zero new DOFs**. | Machine-precision constraint satisfaction is required **and** a stabilized saddle-point solver is in-tree. | ADR-41 Q-DOF |
| 3 | **Self-contact** (a surface vs itself) | 3D self-contact needs penalty-only over-constraint handling, per-node contact-element generation, and a surface-vs-itself broad phase; narrow structural-seismic payoff. | A post-buckling member-on-member / crushing / wrinkling use case lands. | ADR-39 v2/v3; ADR-41 |
| 4 | **Full slide-line Hermite smoothing** (C1-continuous junctions) | Heavy machinery (per-junction Hermite blends, smoothing factor). The cheaper **averaged nodal-normal field** (4a) is the right-sized interim for the faceted-normal chatter problem. | P3 sliding-patch **chatter gate** shows faceted-normal oscillation that 4a cannot fix. | ADR-41 Q-NORMAL (TG §5.1.2) |
| 4a | **Averaged nodal-normal smoothing** (smooth `N(X)` field) | The cheap interim for Q-NORMAL; still deferred from the MVP (faceted per-GP normal ships in P1–P4). | Pull **earlier than #4** if the P3 chatter gate fails — preferred over full slide-line smoothing. | ADR-41 Q-NORMAL (TG §5.1.1) |
| 4b | **Edge-end handoff smoothing** (Munjiza distributed potential — continuity as an edge-edge contact point tumbles off an edge *end* into a vertex/NTS primitive) | Distinct from #4/#4a, which smooth *facet*-normal chatter across facet junctions; #4b smooths the **edge-edge→vertex primitive transition** that ADR-57 hands off. Deferred from the ADR-57 MVP (gap-band + `-visc` are the interim damper). | ADR-57's E7 sliding-rake gate shows force chatter at an edge-end handoff that the gap-band + `-visc` cannot damp. | ADR-57 §Fence / Risks; Munjiza (2004) distributed potential |
| 5 | **Anisotropic / elliptic friction** (two principal `μ`, friction ellipse) | Niche for structural-seismic; doubles the return-map state + oracle surface. Clean extension (scaled-shear reuses the return map) when needed. | Orthotropic rock-joint / laminated-or-fabric-bearing project. *Read TG §5.2.3 Eq.5.2.3-10/11 directly first — the skill's ellipse RHS is schematic.* | ADR-41 Q-ANISO |
| 6 | **Pressure/velocity-dependent `μ(N,v)` wired into contact** | Machinery exists (`frictionModel/` `VelDependent`/`VelPressureDep`) but is unwired to contact surfaces. *Adopt-later, not deferred-forever.* | After the constant-`μ` mortar kernel ships (C3); velocity-weakening sliding-interface dissipation needed. | ADR-41 Q-MUDEP |
| 7 | **Custom global-Uzawa `LadrunoAugmentedNewton` `EquiSolnAlgo`** | Commit-cycle ALM on **stock `NewtonRaphson`** (EmbeddedRebar precedent) is the MVP; subclassing `NewtonRaphson` (header says "not expected … to be subclassed") is the fragile path. | A P2/P3 gate fails *specifically because within-step augmentation is required* **and** the held-load `analyzeAugmented` proc proves insufficient (e.g. a single monotonic step that must land at `augTol`). | ADR-41 Q-DRIVER |
| 8 | **NTN / NTS-via-mortar-weights unification** | Low priority; NTS (ADR-39) and mortar (ADR-41) already coexist behind the `-formulation` selector. | A unification need surfaces (e.g. a single code path simplification that pays for itself). | ADR-41 |
| 9 | **Coupled-multiphysics surfaces** (pore-pressure gap flow, thermal conduction/radiation, Joule, frictional heat, acoustic-structural) + **pressure penetration** | Out of scope for a structural-seismic fork — each is a multiphysics program (own constitutive surface law + field-DOF wiring + validation battery); the fork's targets (pounding/uplift/rocking/sliding) are purely mechanical. | A coupled u-p / thermal-contact / SSI-liquefaction project — likely its **own ADR**, not merely un-deferred here. | Abaqus-scope skeptic (TG §5.2.2/§5.2.4–5.2.7, §6.4) |
| 10 | **Perpendicular edge-edge contact** (cos_t→0; two facets meeting edge-on, the mortar face-clip degenerates, NTS misses the interior-line pair) | Was deferred at B2 (SOFT=2) as a net-new narrow-phase algorithm (segment-segment closest point + common-perpendicular normal) — not a tweak of the shipped face-mortar/NTS lanes. | **TRIGGER FIRED** — the B2 follow-up note made it the one remaining gap in the faceted narrow-phase taxonomy (zero force through an edge-on contact is a *correctness* hole). **PROMOTED → [[57_ladruno_edge_edge_contact_adr]]** (2026-06-24; design-only). | B2 follow-up note ([[48_ladruno_contact_capstone_adr]]) |

## Not deferred-pending — disclosed limitations (recorded so they aren't re-litigated)

- **Explicit ALM.** Full ALM mortar is **implicit-only**; under `CentralDifferenceLadruno` the tangent
  is mass-only so Uzawa degenerates to single-pass penalty. This is a permanent property of the
  accuracy-first lane (complementary to ADR-39's explicit-first NTS), not a deferral. (ADR-41 Q-EXPLICIT)
- **Exact (Lagrange) stick / rough friction.** Both lanes are penalty/Uzawa by construction; penalty
  stick with the `γ_crit = 0.5%·L` elastic-slip bound is the committed engineering trade. Re-open only
  with #2 (true-LM).

## Relationship

This ledger is subordinate to the [[48_ladruno_contact_capstone_adr]] capstone (which owns the
architecture, contracts, and committed roadmap). When a deferred item is re-opened, it graduates to
its **own** numbered ADR; this row then links to it and is marked promoted.
