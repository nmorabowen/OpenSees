---
title: FDEM-lite — cohesive-interface fragmentation (DEFERRED; build-vs-buy unresolved)
project: Ladruno
status: deferred
priority: low
owner: nmora
tags:
  - implementation
  - subsystem
  - fracture
  - cohesive
  - fdem
  - collapse
  - contact
  - explicit
---

# FDEM-lite — cohesive-interface fragmentation (DEFERRED)

> **ADR 54.** **Revised 2026-06-24 after a 4-lens adversarial review** that
> **corrected the reuse targets and refuted the cost framing**. **Status: DEFERRED** —
> the v1 "sidestep" delivers no fragmentation value, the real value gates on the
> shared long pole ([[55_ladruno_contact_runtime_discovery_adr]]), and **build-vs-buy
> is unresolved** (ELS / coupling may dominate building). Corrections folded in as
> `[GATE]`; per-lens verdicts under *Adversarial review log*. The thesis actually gets
> **stronger** once corrected — see the reuse fix.
>
> **Route A** (cohesive-at-boundaries, mesh-biased path) of continuum collapse; the
> path-objective sibling is **Route B**, [[53_ladruno_embedded_discontinuity_adr]].
>
> `[GATE]` **Honesty caveat carried from the parent (§2.4), previously dropped:**
> the program's verified finding is **workflow advantage, NOT accuracy parity**
> ("AEM/FDEM beats FEM at collapse" was **refuted 0-2**). FDEM-lite would land
> *south of* a validated tool, more conveniently — a weak reason to spend quarters.

## What

The Munjiza-faithful continuum-fragmentation route: bulk FE + **zero-thickness
cohesive interface elements** (`G_f` traction–separation) at element boundaries that
soften, fully debond, and become **contact** surfaces; integrated explicitly.

> `[GATE]` **Stop calling the v1 deliverable "fragmentation."** With **predefined**
> fracture+contact surfaces (the v1 sidestep), this ships *predefined separation* —
> breaking apart along surfaces you marked in advance. That is **not** emergent
> fragmentation and demonstrates ~nothing a careful predefined-contact model
> couldn't. Real (emergent) fragmentation needs [[55_ladruno_contact_runtime_discovery_adr]].

## The reuse fix (this is what makes Route A actually lower-risk)

> `[GATE]` **B1 — "reuse the contact kernels for the cohesive law" is REFUTED.** The
> contact kernels are **constitutively the opposite thing**: `LadrunoContactKernel`
> is reversible, compression-only, no-tension (`tn = kn⟨−gap⟩₊`, Macaulay-clamped);
> `LadrunoFrictionKernel` is a pressure-driven return map (`cap = min(μN+c, τmax)`),
> non-softening. A cohesive law needs **tension-carrying, irreversibly-softening,
> `G_f`-bounded, unload-to-origin damage** — the exact thing contact *forbids* by
> design. You cannot express cohesion as a "mode" on the contact kernel.
>
> `[GATE]` **B2 — the right reuse target ALREADY EXISTS:**
> **`LadrunoCohesiveHinge` / `LadrunoCohesiveHingeBiaxial`**
> ([[34_ladruno_cohesive_hinge_biaxial_adr]]) — `G_f`-exact area, rigid-penalty
> pre-peak, exp/linear softening, irreversible secant-to-origin damage, and
> **Benzeggagh–Kenane mixed-mode already coded** (the "mixed-mode criterion" the old
> draft wanted to build). The element scaffold to copy is
> **`ZeroLengthInterface2D` / `ZeroLengthND`** (zero-thickness kinematics + auto-normal
> + local→global scatter). So the new work is the **kinematics + a displacement-jump
> port of an existing, validated softening law** — *not* a green-field cohesive model.
> (The ASDConcrete3D crack-band the old draft cited is `lch`-bound to a continuum
> material and does **not** factor out to a surface law — dropped.) **This — not
> contact-kernel reuse — is why Route A is genuinely lower-risk than E-FEM.**

## Why (corrected, and the counter-weight)

- **Lower-risk than E-FEM — for the right reason.** The **cohesive softening law is
  already built and validated** in the fork (B2). E-FEM's make-or-break (stress-locking)
  has no analogue here; the constitutive heart is done.
- `[GATE]` **But it's a lateral trade, not a free win.** Intrinsic CZM (the v1 form)
  imports its **own** pathologies the old draft booked as "accepted": **mesh-density-
  dependent artificial compliance** (more pre-placed interfaces ⇒ softer bulk — the
  *very* mesh-dependence the continuum route is meant to reduce); the **penalty-
  stiffness ↔ stable-`dt` bind** (the stiffness you raise to cut compliance shrinks
  `dt`, and mass-scaling the cohesive zone to recover `dt` adds spurious inertia *at
  the fracture front*); and **wave reflection** off the pre-placed compliant
  interfaces. Extrinsic/adaptive insertion removes all three but needs runtime
  topology change.
- **Mesh-biased path** is accepted here — but `[GATE]` note the *fragment statistics*
  (count, size distribution — the actual P4 deliverable) are **mesh-controlled**, not
  just the trajectory; the FDEM/AEM literature is explicit. And this directly
  **contradicts the pair's "mesh-independence trifecta" thesis** (Route A throws away
  the path leg Route B says is the point) — a tension the program must resolve.

> `[GATE]` **Build-vs-buy (a co-equal option, not a footnote).** For the named
> research, weigh: (a) **use ELS** for the few continuum-fragmentation cases (>50
> validations, scriptless but mature); (b) **couple OpenSees-continuum ↔ an existing
> DEM/FDEM engine** for the discrete stage (Lu Xinzheng's FEM→DEM precedent) — plays
> to the fork's genuine edge (the *continuum* stage) instead of re-deriving discrete
> mechanics; (c) **just use removal** ([[51_ladruno_element_removal_adr]]) if the
> question is member-scale. Building FDEM-lite is justified only if a *named,
> decision-grade* question needs fragment-resolved answers no bought tool gives.

## Where (corrected)

- **Cohesive law:** `LadrunoCohesiveHingeBiaxial` ([[34_ladruno_cohesive_hinge_biaxial_adr]])
  — port its softening/damage/BK math from rotation-jump/moment to
  displacement-jump/traction (a re-derivation, not a `getCopy`).
- **Element scaffold:** `SRC/element/zeroLength/ZeroLengthInterface2D.{h,cpp}` and
  `ZeroLengthND` — copy the zero-thickness kinematics.
- **Contact (discrete stage):** [[39_ladruno_contact_domain_adr]] kernels — usable for
  the discrete stage on **predefined** surfaces; **NOT** the cohesive law, and **NOT**
  for runtime face creation (see the handoff gap).
- **Explicit substrate:** `CentralDifferenceSMS`, `LadrunoMassScaling`,
  `CriticalTimeStep` — real (subject to the `dt` bind above).
- New `ELE_TAG` for the cohesive interface element; record in
  [[LEDGER_implementations]] + banner at implementation.

## Risks / open questions

> [!question] Q-HANDOFF — the debond→contact transition has NO PLUMBING `[GATE]`
> `addSurface`/`addContact` are **parse-time only**; `LadrunoContactDomain` state is
> keyed by pre-declared `contactTag` and has **no free-body/fragment struct**
> (ADR-51's finding). "Register the faces as a contact pair" is described as a method
> call but is an **unbuilt subsystem** (part of pillar-2). The `G_f`-release energy
> argument is sound *constitutively*; the **registration/timing is the unsolved part**.

> [!question] Q-V1-COHERENCE — the v1 sidestep needs an unstated mechanism `[GATE]`
> v1 ("predefined surfaces") is coherent **only** if the cohesive faces are
> pre-declared as a contact pair **while still bonded** (force-inert, zero-effect
> until debond), so the handoff is a **force-activation of an already-numbered pair**,
> NOT a runtime registration. State this — and the byte-identity break
> (ADR-48 contract #6) + the cost of carrying N latent contact pairs. If v1 instead
> expects to *register* on failure, it is **secretly pillar-2-dependent**.

> [!question] Q-INTRINSIC — compliance / `dt` / mass-scaling / wave-reflection triangle
> Unbounded in the old draft (see Why). Bound it, or commit to **extrinsic/adaptive**
> insertion sooner (removes compliance + reflection).

> [!question] Q-SELFCONTACT — a different problem class, not a bullet `[GATE]`
> A folding/piling fragment is the *typical* case, and the master/slave predefined-pair
> data model **cannot express it** — self-contact needs a single-surface BVH/hash with
> adjacency exclusion ([[55_ladruno_contact_runtime_discovery_adr]] §self-contact).

> [!question] Q-FRAGMENT-STATS — P4 validation must accept mesh-controlled statistics
> Boundary-biased cracking gives mesh-controlled fragment count/size; a validation
> "vs ELS/FDEM" must state the metric tolerates this, or it will fail for the wrong
> reason.

## Phasing (deferred; if ever taken up)

| Phase | Scope | Gate |
|---|---|---|
| **P0 — decision** | Resolve build-vs-buy (ELS / couple / removal) + name a decision-grade use case | A named question no bought tool answers — **else stop here** |
| P1 | Port `LadrunoCohesiveHingeBiaxial` → displacement-jump cohesive interface element (`ZeroLength*` scaffold) + explicit; predefined plane | A bar splits with `∫T·dδ = G_f` exact; byte-identical when no interface fails |
| P2 | Force-activation handoff on **pre-declared-while-bonded** contact pairs (existing contact) | Predefined-surface *separation* demo (NOT fragmentation) with clean release |
| P3+ | **Requires [[55_ladruno_contact_runtime_discovery_adr]]** — emergent fragmentation + self-contact + validation | (gated on pillar-2) |

## Adversarial review log (2026-06-24, 4 lenses)
- **Lens B (cohesive):** REFUTED "reuse contact kernels" (constitutively opposite);
  retargeted reuse to `LadrunoCohesiveHinge[Biaxial]` (already built) + `ZeroLength*`
  scaffold; dropped the `lch`-bound ASDConcrete3D crack-band; reframed "lower-risk" on
  the correct ground; surfaced the intrinsic-CZM compliance/`dt`/reflection triangle.
- **Lens C (contact):** REFUTED "own the expensive 60%" — owned = narrow phase (~40%,
  tractable); unbuilt = discovery/re-emit/edge-edge/self-contact (~60%, the long pole,
  → ADR 55). Handoff has no runtime registration; v1 coherent only via pre-declared-
  while-bonded pairs; self-contact a different class.
- **Lens D (strategy):** scope-creep + build-vs-buy unresolved → **deferred**; §2.4
  caveat restored; v1 "fragmentation" renamed "predefined separation"; trifecta
  contradiction (Route A drops the path leg) flagged.

## See also
- [[55_ladruno_contact_runtime_discovery_adr]] — the shared **long pole** that real
  fragmentation (this) and Route B both require.
- [[53_ladruno_embedded_discontinuity_adr]] — **Route B** (E-FEM, path-objective spike).
- [[34_ladruno_cohesive_hinge_biaxial_adr]] — the cohesive softening law to port.
- [[51_ladruno_element_removal_adr]] — member-scale collapse (the near-term build).
- [[39_ladruno_contact_domain_adr]] · [[48_ladruno_contact_capstone_adr]] — contact.
- Primary: Munjiza (2004); Camacho & Ortiz (1996); Ortiz & Pandolfi (1999);
  Hillerborg et al. (1976); Bažant & Oh (1983).
