---
title: Embedded Strong Discontinuity (E-FEM) — a KILL-GATED research spike for mesh-objective fracture
project: Ladruno
status: spike (kill-gated)
priority: low
owner: nmora
tags:
  - implementation
  - element
  - fracture
  - discontinuity
  - localization
  - eas
  - spike
---

# Embedded Strong Discontinuity (E-FEM) — a kill-gated spike

> **ADR 53.** **Revised 2026-06-24 after a 4-lens adversarial review** that
> **falsified the framing** (not the instinct). Demoted from a product ADR to a
> **single-gate research spike**: one make-or-break question — **stress-locking on a
> fully-separated element** — decides go/no-go **before any phasing spend**. If it
> doesn't clear cheaply on the EAS brick, **kill it**. Corrections folded in below
> as `[GATE]`; per-lens verdicts under *Adversarial review log*.
>
> **The path-objectivity leg** of the mesh-independence trifecta. The lower-risk,
> ship-first sibling is **Route A**, [[54_ladruno_fdem_lite_adr]]. **Both** gate on
> the shared long pole, [[55_ladruno_contact_runtime_discovery_adr]] (contact
> pillar-2). Framing: [[50_aem_opensees_scoping_report]].
>
> `[GATE]` **Honesty caveat carried from the parent (§2.4), previously dropped:**
> the program's *verified* finding is a **workflow/automation advantage, NOT proven
> accuracy supremacy** — "AEM/FDEM beats FEM at collapse" was **refuted 0-2** in the
> deep-research verification. E-FEM does not change that; it improves *mesh-bias*, not
> *accuracy parity*.

## What

A **research spike** to test whether an **embedded strong discontinuity (E-FEM /
SDA)** — a displacement-jump mode letting a crack run *through* a `LadrunoBrick`
element, governed by a cohesive `G_f` law — is viable on the fork's **shipped EAS
brick** ([[19_ladruno_brick_eas_simo_rifai]]).

> `[GATE]` **The EAS-reuse claim, corrected (was the load-bearing oversell).** Only
> the **static-condensation block algebra** reuses (`LadrunoBrick.cpp:2733-2734` —
> bordering `[α]`→`[α;ζ]` is trivial). **Everything that makes E-FEM E-FEM is
> green-field:**
> - **Internal equation** differs *in kind*: EAS is a *volume* L2-orthogonality
>   `∫MᵀσdV=0` (`:2678`); the cohesive condition is a *surface* traction balance
>   `T(⟦u⟧)=σ·n` with a path-dependent softening `T`. Not the same residual family.
> - **M-operator** differs: `computeMenh` builds mean-zero odd bubbles for the EAS
>   patch test (`:2545-2563`); the SDA `G_d=(φ⊗n)ˢ` is a regularized *surface jump*
>   with an orientation — a new operator, not a parameter.
> - **Inner solver class change:** shipped is a **fixed-9, SPD, one-step-linear**
>   `Kaa.Solve` (`:2695`); a softening, `σ·n`-coupled, SKON-non-symmetric ζ makes the
>   inner system **variable-dim, indefinite, nonlinear**. That's a *rewrite* of the
>   inner Newton, not a widening.
> - **Lifecycle mismatch:** `α` is born at t=0 and committed every step; **ζ activates
>   mid-analysis on a criterion** with a frozen normal `n` + onset state that must
>   commit/revert/serialize. The shipped element has no "DOF activated at step k".
> - **Factual:** `condenseEAS` **does not exist** as a method — the condensation is
>   **inlined** (`:2650-2651`, `:2733-2734`); ADR-19 planned to extract it and never
>   did. A P0 would have to extract it first.
>
> Honest split: *reuse the condensation algebra; build a new internal-equation solver
> beside it.* "Build-on-investment" applies to the **plumbing**, not the **mechanics**.

## Why (and the honest counter-weight)

The **path-objectivity** leg `G_f` cannot fix: crack-band regularization (shipped in
the concrete models) makes *dissipation* mesh-objective, but the crack still follows
element boundaries. E-FEM is the variational route to a path that runs where the
physics wants — and it's the consistent limit of strain-softening (Simo–Oliver–Armero
1993), curing the localization ill-posedness on a `G_f` surface.

> `[GATE]` **"Mesh-objective path" downgraded — it was overclaimed.** v1 is
> **element-local with no global tracking**, so the crack surface need not be
> continuous across faces and **can zigzag** (Q-TRACKING). What v1 actually delivers
> is **objective *dissipation* (already shipped via `G_f`) + a *less* mesh-biased but
> NOT orientation-objective path.** Honest headline: *element-local embedded
> discontinuity (untracked)*, not "mesh-objective fracture."

> `[GATE]` **Build-vs-buy is unresolved and weighs against building.** Assembling
> `G_f` + E-FEM + contact ≈ Munjiza FDEM. The parent §2.4 says the payoff is
> *workflow, not accuracy*; ELS (>50 validations) and real FDEM exist. **For the
> fork's named research (retrofit / collapse-risk), member-scale removal
> ([[51_ladruno_element_removal_adr]]) plausibly answers the question, and
> coupling/ELS dominates building a research-grade fragmenter.** This spike is
> justified ONLY if a *named, decision-grade* question needs a mesh-objective crack
> path that no bought tool provides — name it before spending past P0.

## Where (corrected)

- **`SRC/element/ladrunoBrick/LadrunoBrick.cpp`** — the EAS seam: the **inlined**
  condensation (`:2650-2651`, `:2733-2734`), `computeMenh` (`:2545-2563`), the `α`
  state (`:266-273`, `:363`, `:392`) — the **template** for ζ's commit/serialize, NOT
  a container that hosts it. A P0 extracts `condenseEAS` from the inlined block.
- **New (all green-field):** the surface operator `G_d`, the onset criterion + frozen
  `n`, the cohesive surface law, the non-symmetric variable-dim inner solve, and the
  new committed state (ζ, `n`, activation flag, onset strain) + its serialization.
- **Cohesive law:** reuse [[34_ladruno_cohesive_hinge_biaxial_adr]]'s
  `LadrunoCohesiveHingeBiaxial` softening math (G_f-exact, Benzeggagh–Kenane) — a
  re-derivation to a displacement-jump/traction surface law, **not** the contact
  kernel (see [[54_ladruno_fdem_lite_adr]] for why).

## Risks / open questions

> [!question] Q-STRESSLOCK — **the single kill gate (P0)**
> Does a fully-separated element transmit ~zero traction, or does it lock (KOS/SOS)?
> SKON is the textbook cure but is **non-symmetric** ⇒ `FullGeneral`/`UmfPack` AND a
> non-symmetric inner solve. **Decide go/no-go here before any further spend.**

> [!question] Q-INNER-SOLVE `[GATE: new]`
> SKON ⇒ the element-internal system loses SPD and the shipped "converges in 1 step"
> linearity (ADR-20). The `Kaa.Solve` LU runs but its contract is gone — inner-Newton
> robustness under a softening law is unproven. Part of the P0 gate.

> [!question] Q-EAS-HOURGLASS — **mostly a non-issue; was based on a false premise**
> `[GATE]` ADR-20 is `status: rejected` — there is **no shipped EAS hourglass
> stabilization to survive**, and the small-strain brick has **no spurious hourglass
> mode** (ADR-19 rank test: 6 ZEMs). The real EAS hourglass risk is *finite-strain
> compressive*, which this spike **excludes** (small strain only). Downgraded.

> [!question] Q-ONSET — bifurcation onset is non-trivial `[GATE: hardened]`
> Loss-of-ellipticity onset (`det(n·Cep·n)=0` over `n`) needs the continuum tangent
> per GP + an eigen-search, none shipped; and the brick has **8 GPs** — which `n`
> does the single embedded surface take (8-GP→single-normal reduction undefined)? And
> onset timing in the solve cycle is undefined under `algorithm Linear` (the ADR-19
> `update()` staleness trap). Rankine is the brittle fallback.

> [!question] Q-TRACKING — element-local path can zigzag (not orientation-objective)
> No global tracking in v1; accepted, but it caps the "mesh-objective path" claim.

> [!question] Q-BUILD-VS-BUY — name the decision-grade question, or don't build
> See Why. If no named question needs a mesh-objective path that ELS/coupling can't
> give, **this spike ends at P0 regardless of the stress-locking result.**

## Phasing (spike, not product)

| Phase | Scope | Gate |
|---|---|---|
| **P0 — KILL GATE** | Extract `condenseEAS`; add a ζ jump mode + a `LadrunoCohesiveHingeBiaxial`-derived surface law on ONE brick; SKON; **stress-locking test on a fully-separated single element** + inner-solve robustness | Fully-separated element transmits ~0 traction AND the inner solve converges. **Fail ⇒ kill.** Also requires a *named* decision-grade use case (Q-BUILD-VS-BUY). |
| P1+ | Only if P0 clears: onset criterion + normal; notched-bar mesh-objectivity (dissipation); handoff to separation+contact | (re-scope after P0; depends on [[55_ladruno_contact_runtime_discovery_adr]]) |

## Adversarial review log (2026-06-24, 4 lenses)
- **Lens A (E-FEM mechanics):** REFUTED-by-overclaim the "condensed identically / exactly the EAS plumbing" framing (only the condensation algebra reuses; internal eqn, M-operator, inner solver, activation lifecycle all green-field); FACTUAL `condenseEAS` doesn't exist; REFUTED Q-EAS-HOURGLASS premise (ADR-20 rejected, no hourglass mode small-strain); WEAKENED "mesh-objective path" (zigzag); new risks: non-symmetric inner solve, onset timing/`update()`, 8-GP normal, onset-state serialization. The engineering instinct survives; the "ride the seam" sell does not.
- **Lens D (strategy):** §2.4 honesty caveat dropped → restored; build-vs-buy unresolved → E-FEM demoted to a kill-gated spike; the trifecta self-contradicts (Route A abandons the path leg this ADR says is essential) → flagged.
- Result: demoted to a single-gate spike; framing corrected; build-vs-buy and a named use case gate the spend.

## See also
- [[54_ladruno_fdem_lite_adr]] — **Route A** (cohesive interface, lower-risk, ship-first).
- [[55_ladruno_contact_runtime_discovery_adr]] — the shared **long pole** (contact
  pillar-2); both routes gate on it.
- [[19_ladruno_brick_eas_simo_rifai]] · [[34_ladruno_cohesive_hinge_biaxial_adr]] —
  the EAS condensation algebra + the cohesive softening law to reuse.
- [[51_ladruno_element_removal_adr]] — member-scale collapse (the near-term build).
- Primary: Simo & Rifai (1990); Simo, Oliver & Armero (1993); Oliver (1996); Jirásek
  (2000, SOS/KOS/SKON); Linder & Armero (2007); Wriggers & Reese (1996).
