---
title: LadrunoCohesive — a predefined-interface cohesive zone element (SKETCH, study later)
project: Ladruno
status: sketch
priority: medium
owner: nmora
tags:
  - implementation
  - element
  - cohesive
  - interface
  - masonry
  - fracture
  - sketch
---

# LadrunoCohesive — predefined-interface cohesive zone element

> **ADR 56 — SKETCH / study-later seed.** Incomplete by design (per the folder
> convention). Captures the idea + named users + reuse + the questions to resolve
> before committing. Born from the continuum-element-removal analysis: binary
> erosion is ruled out; the **energy-honest** way to separate a continuum is a
> cohesive `G_f` debond — and a cohesive element on a **known interface** delivers
> most of that value with **none** of the deferred fragmentation machinery.
>
> **Decoupled, deliberately, from:** continuum fragmentation
> ([[54_ladruno_fdem_lite_adr]], deferred), runtime contact discovery
> ([[55_ladruno_contact_runtime_discovery_adr]], the long pole), and E-FEM
> ([[53_ladruno_embedded_discontinuity_adr]], a spike). This element stands on its
> **own predefined-interface users** — same logic as billing pillar-2 to the contact
> lane: build the useful shared component for its own users, not for the speculative
> program that would later consume it.

## What

A **zero-thickness cohesive interface element** between two pre-meshed bulk faces,
carrying a **traction–separation law** (Mode I normal + Mode II/III shear, mixed-mode)
whose softening area = `G_f`. Holds the interface together pre-peak, dissipates `G_f`
through softening, transmits ~zero traction once fully open (then hands to
**predefined-surface** contact, which already exists). **Implicit / quasi-static
first** (the use cases are static-to-cyclic, which sidesteps the explicit-`dt`
pathology of CZM).

## Why — it has real, named, decision-grade users (no fragmentation needed)

- **Masonry mortar joints — the killer app.** Masonry cracks at the brick–mortar
  interface; a `G_f` mixed-mode cohesive joint is the standard model. (Why the AEM
  library is masonry-heavy — this gets most of AEM-for-masonry without building AEM.)
- **RC / precast interfaces** — construction / cold / lift joints, segmental joints,
  cover-spalling planes, shear-friction interfaces.
- **Bond & debonding** — rebar–concrete bond-slip, FRP-jacket / epoxy debonding
  (retrofit), composite delamination.
- **Rock joints / soil–structure interfaces.**

All **predefined** interfaces → no runtime discovery, no self-contact, no
fragmentation, no pillar-2. And it is the **energy-honest replacement for binary
continuum element removal** on a known plane.

## Reuse (most of it already exists — keep the build small)

- **The constitutive law:** `LadrunoCohesiveHingeBiaxial`
  ([[34_ladruno_cohesive_hinge_biaxial_adr]]) — `G_f`-exact, rigid-penalty pre-peak,
  exp/linear softening, irreversible secant-to-origin damage, **Benzeggagh–Kenane
  mixed-mode**. Port from rotation-jump/moment → displacement-jump/traction
  (a re-derivation, not a `getCopy`).
- **The element scaffold:** `ZeroLengthInterface2D` / `ZeroLengthND` — zero-thickness
  kinematics, auto-normal, local→global scatter.
- **Convergence through softening:** the fork's viscous (Duvaut–Lions) regularization
  (`-visc`) + arc-length / dissipation lanes.

## Scope / non-goals

- **In:** predefined interface, Mode I/II + mixed-mode, `G_f`, implicit/quasi-static,
  debond→**predefined** contact handoff.
- **Out (sibling ADRs):** runtime insertion / emergent fragmentation
  ([[54_ladruno_fdem_lite_adr]]); runtime contact discovery + self-contact
  ([[55_ladruno_contact_runtime_discovery_adr]]); path-objective cracks-through-elements
  ([[53_ladruno_embedded_discontinuity_adr]]); explicit-dynamics CZM (the `dt`/compliance
  bind — revisit only if a dynamic use case appears).

## Open questions — TO STUDY LATER

> [!question] Q-EXISTING — check the tree before building (don't rebuild)
> What cohesive/interface pieces already exist? Stock OpenSees has `Bond_SP01`,
> `ZeroLengthInterface2D/3D`, `ZeroLengthContact*`, and some cohesive/interface
> materials; the fork has the cohesive-hinge law. **Scope the build around the real
> gap** — a general `G_f` mixed-mode *surface* cohesive element with a clean contact
> handoff — not a rebuild of what's there. (First task.)

> [!question] Q-INTRINSIC — intrinsic vs extrinsic
> Intrinsic (stiff elastic pre-peak, pre-placed) adds compliance but is trivial to
> insert; extrinsic (rigid-until-onset) avoids compliance but needs activation logic.
> For *predefined static* interfaces, intrinsic with a stiff penalty is likely fine
> (the `dt` bind is an explicit problem). Confirm.

> [!question] Q-LAW-PORT — rotation/moment → displacement/traction
> How much of `LadrunoCohesiveHingeBiaxial` ports cleanly? The damage/softening/BK
> algebra should transfer; the kinematic mapping and the penalty/units change.

> [!question] Q-HANDOFF — debond → predefined contact
> When fully open, register/activate the (predefined) contact pair so the faces carry
> no tension + friction on closure. Predefined ⇒ this rides existing contact (cf.
> [[54_ladruno_fdem_lite_adr]] Q-V1-COHERENCE: pre-declare the pair while bonded,
> force-inert until debond).

> [!question] Q-VALIDATION — masonry-first
> Validate against a masonry wall / triplet shear / couplet bond test (mortar-joint
> `G_f`, friction) — the named user — plus a mesh-objectivity check on the dissipation.

> [!question] Q-DE-RISKS-54 — bonus
> A validated predefined-interface cohesive element is the de-risked front-end if
> continuum fragmentation ([[54_ladruno_fdem_lite_adr]]) is ever taken up — only the
> insertion + discovery parts would remain. Note but don't scope here.

## Phasing (sketch)

| Phase | Scope |
|---|---|
| P0 | Tree check (Q-EXISTING); decide reuse vs new; scope the gap |
| P1 | Port the law → `LadrunoCohesive` interface element (Mode I, predefined plane, implicit); single-element `∫T·dδ = G_f` exact |
| P2 | Mode II + Benzeggagh–Kenane mixed-mode; viscous-regularized convergence |
| P3 | Masonry-joint validation; debond→predefined-contact handoff |

## See also
- [[34_ladruno_cohesive_hinge_biaxial_adr]] — the cohesive law to port.
- [[54_ladruno_fdem_lite_adr]] — the *fragmentation* sibling (this is its P1, decoupled).
- [[55_ladruno_contact_runtime_discovery_adr]] · [[53_ladruno_embedded_discontinuity_adr]].
- [[51_ladruno_element_removal_adr]] — member-scale removal; this is the
  energy-honest continuum-interface counterpart.
- Primary: Dugdale (1960); Barenblatt (1962); Hillerborg et al. (1976, fictitious
  crack); Camacho & Ortiz (1996); Ortiz & Pandolfi (1999); Lourenço (1996, masonry
  interface CZM).
