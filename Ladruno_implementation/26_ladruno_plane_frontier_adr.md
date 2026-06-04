---
title: 2D continuum elements — frontier & contribution thesis
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - adr
  - element
  - solid
  - 2d
  - fracture
  - localization
  - frontier
  - research
---

# 2D continuum elements — frontier & contribution thesis

> Sibling: [[25_ladruno_plane_elements_adr]] — the Phase-0 substrate (brick
> technology ported to 2D, byte-verified against upstream). **This** document
> argues *why the 2D family is worth building beyond a port*: it is the platform
> on which the fork can ship genuinely novel element technology and rigorous
> localization limiters — developed and validated in 2D, then **lifted to the
> brick**. It is a decision document, not a build plan; each track that we
> green-lights gets its own ADR.

## Thesis

A 2D port of `LadrunoBrick` is *necessary substrate but not a contribution* — the
small-strain locking cures (std/bbar/SSP/EAS) already exist scattered across
upstream (`Tri31`, `FourNodeQuad`, `SSPquad`, `EnhancedQuad`). What elevates the
family from port to contribution are formulations that **fit the OpenSees
element black box** yet are at or beyond the 2025 state of the art:

1. **Coupled-field localization limiters** — implicit **gradient-enhanced
   damage** and **phase-field fracture** — assembled exactly like the existing
   `FourNodeQuadUP` u-p element (an extra nodal DOF, statically coupled, **zero
   analysis-core change**).
2. **Element-technology frontier that respects the black box** — the **Virtual
   Element Method** (arbitrary polygons) and the **Scaled Boundary FEM**
   super-element (analytic crack-tip singularities + unbounded domains).

And the strategic multiplier: **2D is the proving ground for 3D.** Every one of
these is dramatically cheaper to get right where meshes are cheap, a Newton step
is inspectable, and crack paths are visible. The plane family is the lab; the
brick is production. We prove a formulation in `LadrunoQuad`, then lift it.

## The governing constraint (the filter every idea is judged against)

An OpenSees element is a black box: given nodal displacements at *its own* nodes,
it returns **K** and **R**, consuming a material point through
`setTrialStrain`/`setTrialF`. The "analysis core" we will **not** touch = DOF
numbering, `FEM_ObjectBroker`, solver/integrator, constraint handler. So:

| An element legally **can** | An element **cannot** (without core surgery) |
|---|---|
| Choose its own integration / strain operator | Add/remove DOFs at runtime (XFEM enrichment) |
| Carry internal non-nodal state and **statically condense** it (EAS α, E-FEM jump) | See another element's Gauss points (true nonlocal integral) |
| **Declare its own `ndf` per node** — add a nodal scalar field | Drive a staggered/alternating global solve unaided |
| Have a **variable node count** (n-gon) | Remesh / h-adapt |
| Override `getCharacteristicLength` | |

Two precedents in *our* tree make this concrete and decide most of what follows:

- **`FourNodeQuadUP`** already adds a pore-pressure **nodal DOF** (ndf = 2+1) and
  assembles a coupled **K** with **zero** core changes. ⇒ *any micromorphic /
  coupled-field element is legal.* This is the key that unlocks gradient-damage
  and phase-field.
- **`LadrunoBrick`** already does **static condensation** of internal EAS modes,
  and `LadrunoArcLength`/`LadrunoDynamicRelaxation` show we can ship custom
  integrators without core surgery. ⇒ *E-FEM jumps and a staggered phase-field
  driver are in-character.*

Notation throughout: direct / index / Voigt where it illuminates; 2D Voigt
{11, 22, 12}, engineering shear γ₁₂ = 2ε₁₂.

## The CST verdict (why the triangle is honestly a dead end)

Plain CST is not merely "low order"; for the fork's target physics it is the
*wrong* substrate, for two compounding reasons:

- **Incompressible locking is catastrophic, not gradual.** The isochoric
  constraint tr(ε)=0 is *one constraint per element*; with ~½ as many elements as
  nodes, the constraint-to-DOF ratio chokes J2 flow at the perfectly-plastic
  limit. B-bar cannot save a 1-point triangle (nothing to average).
- **Localization is mesh-direction biased.** A constant-strain field cannot carry
  a shear band at an arbitrary angle — bands snap to mesh lines.

`Tri6` (LST) fixes the gradient but is quadratic (midside-node cost, distortion
fragility) and still has *no* locking cure. So in [[25_ladruno_plane_elements_adr]]
`LadrunoCST` ships only as the baseline / E-FEM carrier. The *real* triangle
frontier (strain smoothing, below) is paradigm-breaking and tracked separately.

## Frontier tracks

Each: theory in brief, OpenSees fit against the filter, port-vs-contribution
verdict, key references.

### T1 — Implicit gradient-enhanced damage (coupled u–ē) · **CONTRIBUTION, recommended first**

The nonlocal integral that regularizes softening (Pijaudier-Cabot & Bažant 1987)
is replaced by its PDE form (Peerlings, de Borst, Brekelmans, de Vree 1996): a
nonlocal equivalent strain ē governed by a Helmholtz screening equation.

- Direct:  ē − c ∇²ē = e_eq(ε),  c = ℓ_c²
- Index:   ē − c ē_{,ii} = e_eq
- Weak/Voigt: ∫(δē ē + c ∇δē·∇ē) dΩ = ∫ δē e_eq dΩ  → a standard
  **diffusion–reaction stiffness block** in ē, coupled to the displacement block
  through the damage D(ē) that degrades the stress σ = (1−D)·σ̄.

The element carries **u (2 dofs) + ē (1 dof) per node** and assembles a coupled
3-field-per-node **K** — *architecturally identical to `FourNodeQuadUP`. No core
change.* It introduces a true material length ℓ_c, gives mesh-objective **and**
bias-reducing localization, and converges to a finite-width band.

- **Fit:** ✅ exactly the u-p precedent. **Verdict:** highest rigor-per-cost
  upgrade we can make; the natural *first* coupled-field element and the proof
  step before phase-field. **Lift to 3D:** straightforward once 2D is proven.
- **Refs:** Peerlings et al., *IJNME* 39 (1996) 3391; Pijaudier-Cabot & Bažant,
  *J. Eng. Mech.* 113 (1987) 1512; de Borst & Mühlhaus (gradient plasticity).

### T2 — Phase-field fracture (coupled u–d) · **CONTRIBUTION, the headline**

A crack is a smeared scalar field d∈[0,1] with its own Ginzburg–Landau-type PDE
and length ℓ (Bourdin–Francfort–Marigo 2000; Miehe, Welschinger, Hofacker 2010):

- Energy: Ψ = ∫ g(d) ψ⁺(ε) + ψ⁻ dΩ + G_c ∫ ( d²/2ℓ + (ℓ/2)|∇d|² ) dΩ
- d-equation (Voigt/code): another diffusion–reaction block, *structurally the
  same shape as T1*, driven by the positive elastic energy ψ⁺ (tension/compression
  split) through a history field ℋ for irreversibility.

Nucleates, branches, and merges cracks with **no tracking, no enrichment**, and
consumes our existing concrete material's effective stress through g(d). Wu's
unified phase-field/cohesive theory (2017) ties it *back* to crack-band and
cohesive fracture — i.e. it is the rigorous generalization of the Tier-0 lₐ
regularization [[25_ladruno_plane_elements_adr]] already ships.

- **Fit:** ✅ coupled u–d element (add an ndf) **+** ideally a **staggered**
  (alternating-minimization) driver for robustness — and we already ship custom
  integrators (`LadrunoArcLength`, `LadrunoDynamicRelaxation`) without core
  surgery, so a staggered `LadrunoPhaseFieldSolver` is in-character. Monolithic
  Newton also works (non-convex → needs the line-search/arc-length we have).
- **Verdict:** the flagship fracture contribution; **2D-first is the correct
  order** (validation against Miehe's single-edge-notch tension/shear benchmarks
  is a 2D exercise). **Lift to 3D:** the brick gains it after 2D proves the
  staggering and the split.
- **Refs:** Bourdin, Francfort, Marigo, *JMPS* 48 (2000) 797; Miehe, Welschinger,
  Hofacker, *IJNME* 83 (2010) 1273; Miehe, Hofacker, Welschinger, *CMAME* 199
  (2010) 2765; Wu, *JMPS* 103 (2017) 72 (unified PF/cohesive).

### T3 — Virtual Element Method polygons · **CONTRIBUTION (element technology)**

First-order VEM on an arbitrary polygon. The stiffness splits into a
**consistency** part (exact for rigid + constant-strain fields, from a polynomial
L²-projection Π needing only boundary integrals — no interior quadrature) plus a
**stabilization** part (a dof-based algebraic term, morally identical to our
hourglass control):

- Direct: **K** = **K**_c + **K**_s = (projection of the exact constant-strain
  response) + α (I − Π)ᵀ(I − Π).

It talks to `LadrunoJ2`/`ASDConcrete` through the *same* `setTrialStrain` seam at
the projection point. Crucially, **a polygonal VEM element is still one element
owning its own nodes** — just a variable node count, which OpenSees supports.
Hanging nodes are trivial (a T-junction is a polygon with a collinear vertex) →
this is the natural substrate for quadtree adaptivity and image-/grain-based
meshing. Nonlinear/finite-strain VEM is *active research* (open problem: the
stabilization-energy choice in the nonlinear range — exactly the question our
damage-scaled hourglass already reasons about), so a finite-strain `LadrunoVEM`
plane element sits at the publishable edge.

- **Fit:** ✅ variable-node element, standard material seam. **Verdict:** the
  distinctive *element-technology* contribution that respects the black box.
  **Lift to 3D:** polyhedral VEM is much harder — this is a genuinely 2D-first
  capability (a reason the 2D family is not just a stepping stone).
- **Refs:** Beirão da Veiga, Brezzi, Cangiani, Manzini, Marini, Russo, *M3AS* 23
  (2013) 199; Wriggers, Reddy, Rust, Hudobivnik, *Comput. Mech.* 60 (2017) 253;
  Chi, Beirão da Veiga, Paulino, *CMAME* 318 (2017) 148.

### T4 — Scaled Boundary FEM super-element · **CONTRIBUTION (most fork-aligned)**

Semi-analytical: discretize only the boundary, solve **analytically in the radial
direction** (Wolf & Song 1996). Two payoffs that compound *our existing
investment*:

1. **Crack-tip singularities + analytic stress-intensity factors** with no tip
   refinement — the rigorous discrete-fracture counterpart to T1/T2's smeared
   cracks.
2. **Unbounded domains represented exactly** — a natural sibling of the fork's
   absorbing-boundary / PML / DRM soil-structure stack
   ([[absorbing_boundaries_and_pml_guide]]).

Architecturally a many-node **super-element** with a dense **K** → fits as an
element with no core change.

- **Fit:** ✅ super-element. **Verdict:** the most *differentiated* direction
  because it stacks on existing fork strengths (wave truncation, SSI). **Lift to
  3D:** 2D SBFEM is mature; 3D is research — 2D-first again.
- **Refs:** Wolf & Song, *Finite-Element Modelling of Unbounded Media* (1996);
  Song & Wolf, *CMAME* 147 (1997) 329; Song, *The Scaled Boundary FEM* (2018).

### T5 — Strain-smoothing triangles (ES-/NS-/bES-FEM) · contribution but **PARADIGM-BREAKING**

Replace the compatible strain with one smoothed over a smoothing domain via the
divergence theorem:

- Direct: ε̄ = (1/A_s)∮_{∂Ω_s} sym(u ⊗ n) dΓ  — needs only shape-function
  *values on the cell boundary*; **no Jacobian, no derivatives** ⇒ distortion-immune,
  extends to polygons.

**ES-FEM** (edge-based) gives near-`Tri6` accuracy at T3 cost (superconvergent);
**NS-FEM** (node-based) is volumetric-locking-free automatically (nodal B-bar);
**bES-FEM** (selective edge/node) is accurate *and* lock-free — the genuinely
better-than-CST-and-Tri6 triangle for plasticity.

- **Fit:** ⚠️ an edge/node smoothing domain spans *multiple elements' nodes* —
  S-FEM is fundamentally a **nodal/edge** assembly, not an element assembly,
  which is exactly what the OpenSees core assumes you won't do. Realizable only as
  a **patch/macro-element** owning a cluster of triangles → changes the
  user-facing model. **Verdict:** real contribution, but a *separate ambitious
  track*, not part of the plane-family port. Flag, don't schedule yet.
- **Refs:** Liu, Nguyen-Thoi, Lam, *J. Sound Vib.* 320 (2009) 1100 (ES-FEM);
  Liu et al., *Comput. Struct.* 87 (2009) 14 (NS-FEM); Liu & Nguyen-Thoi,
  *Smoothed Finite Element Methods* (CRC, 2010).

### T6 — Pian–Sumihara assumed-stress hybrid quad · cheap real win

Variationally a Hellinger–Reissner element; the most elegant locking-free Q4 and,
unlike EAS, it **does not hourglass in finite-strain compression** (the EAS
failure flagged in [[25_ladruno_plane_elements_adr]] Risks). A natural
`-formulation hybrid` once the finite path exists.

- **Fit:** ✅ standard element (assumed-stress modes condensed). **Verdict:**
  small, self-contained upgrade; the right *finite-strain* companion to / 
  replacement for `eas` in compression-dominated (concrete) regimes.
- **Refs:** Pian & Sumihara, *IJNME* 20 (1984) 1685; de Souza Neto Ch. 15
  (EAS-finite hourglass); Wriggers & Reese, *CMAME* 135 (1996) 201.

### T7 — Stabilized equal-order u–p (Dohrmann–Bochev) · cheap real win

Lets the *simplest* equal-order u and p interpolations satisfy near-incompressibility
without LBB, via a polynomial-pressure-projection stabilization (a small
∫(p−Πp)(q−Πq)/2μ block). Works identically on triangles and quads — the modern
mixed cure, cheaper than Taylor–Hood, more general than Q1P0.

- **Fit:** ✅ a few lines on top of a mixed element. **Verdict:** clean rigorous
  incompressibility for both `LadrunoQuad` and `LadrunoCST`.
- **Refs:** Dohrmann & Bochev, *IJNMF* 46 (2004) 183.

## Summary — port vs. contribution, against the filter

| Track | Rigor | Fits black box? | Verdict |
|---|---|---|---|
| std/bbar/SSP/EAS/corot/finite (ADR 25) | — | ✅ | **Port.** Necessary substrate. |
| T6 Pian–Sumihara hybrid | med | ✅ | Cheap win (non-hourglassing finite quad). |
| T7 Stabilized u–p | med-high | ✅ | Cheap win (rigorous incompressibility, tri+quad). |
| Crack-band + lₐ (ADR 25 P2) | low | ✅ | Floor. Ships now. |
| **T1 Gradient-enhanced damage** | high | ✅ (u-p precedent) | **Contribution — recommended FIRST.** |
| **T2 Phase-field fracture** | high | ✅ element + staggered driver | **Flagship contribution.** |
| **T3 VEM polygons** | high | ✅ (variable node count) | **Contribution (element tech, 2D-native).** |
| **T4 SBFEM super-element** | high | ✅ (super-element) | **Contribution (most fork-aligned).** |
| T5 ES-/NS-FEM triangles | high | ⚠️ patch/macro only | Contribution, **paradigm-breaking — separate track.** |
| Nonlocal *integral* damage | high | ⚠️ needs global GP registry | Skip in favor of its PDE form (T1). |
| XFEM/GFEM | high | ❌ core surgery | Out of scope (dynamic DOFs + level sets). |

## The 2D → 3D ladder (why the lab matters)

The intended flow for every contribution track:

1. **Prove in `LadrunoQuad`** — visible crack paths, cheap meshes, inspectable
   Newton steps, benchmarks that are 2D by construction (Miehe SENT/SENS for T2,
   notched-bar localization for T1).
2. **Lift the coupled-field / material machinery to `LadrunoBrick`** — the damage
   PDE block, the phase-field split, the staggering all generalize; the brick
   gains them as a new `-regularization`/coupled mode once 2D is green.
3. **VEM/SBFEM stay 2D-distinctive** (T3 polyhedral and T4 3D are research-grade)
   — these are the tracks that justify the 2D family *on its own*, not only as a
   stepping stone.

## Recommended decision

Build [[25_ladruno_plane_elements_adr]] (Phase 0 substrate) to completion, then
take **T1 (gradient-enhanced damage)** as the first contribution — it is the
cheapest rigorous step, it directly exercises the `FourNodeQuadUP` coupled-field
pattern we will reuse for everything else, and it de-risks **T2 (phase-field)**,
the headline. **T3 (VEM)** and **T4 (SBFEM)** are parallel, independently
valuable element-technology tracks; schedule by appetite. **T6/T7** are cheap
upgrades to fold into the substrate when the finite/mixed paths land. **T5** is
flagged for a future, separate macro-element ADR.

Each green-lit track gets its own ADR (proposed numbering: T1 → 27, T2 → 28,
T3 → 29, T4 → 30) with the full coupled-field weak form, the element DOF layout,
the staggered/monolithic solve decision, and the 2D benchmark battery.

## Reserved class tags (per [[LEDGER_implementations]] convention — RESERVED, not in classTags.h until built)

Element registry (band ≥ 33000; 33000–33008 used: …EmbeddedNode=33006,
LadrunoQuad=33007, LadrunoCST=33008 — ADR 25 shipped at 33007/33008 after a
sibling PR took the originally-reserved 33006 for LadrunoEmbeddedNode):

- `ELE_TAG_LadrunoVEM     = 33009` — T3 polygonal VEM (RESERVED)
- `ELE_TAG_LadrunoSBFEM   = 33010` — T4 scaled-boundary super-element (RESERVED)
- Coupled-field elements (T1/T2): decide at their ADRs whether they are **new
  element classes** (`ELE_TAG_LadrunoGradientQuad`, `…PhaseFieldQuad`) or a
  **coupled mode** of `LadrunoQuad` that raises ndf — the latter is cleaner but
  the brick's serialization packs no ndf-switch today. Reserve `33011`–`33012`
  in the element registry for them.

> Tags are per-registry; these element-registry numbers do not collide with the
> identical integrator/material numbers already in use.

## References (consolidated)

- de Souza Neto, Perić, Owen — *Computational Methods for Plasticity* (2008):
  Ch. 9 plane-stress projection, Ch. 12 damage, Ch. 15 F-bar/EAS/u-p + EAS-finite
  hourglass.
- Belytschko, Liu, Moran, Elkhodary — *Nonlinear FE for Continua and Structures*
  (2014): Ch. 8 element technology (Hu-Washizu, B-bar, patch tests).
- Peerlings et al. 1996; Pijaudier-Cabot & Bažant 1987 (gradient/nonlocal damage).
- Bourdin–Francfort–Marigo 2000; Miehe et al. 2010 (×2); Wu 2017 (phase-field).
- Beirão da Veiga et al. 2013; Wriggers et al. 2017; Chi et al. 2017 (VEM).
- Wolf & Song 1996; Song & Wolf 1997; Song 2018 (SBFEM).
- Liu, Nguyen-Thoi et al. 2009 (×2); Liu & Nguyen-Thoi 2010 (S-FEM).
- Pian & Sumihara 1984; Wriggers & Reese 1996 (hybrid stress / EAS instability).
- Dohrmann & Bochev 2004 (stabilized equal-order u-p).
- Bažant & Oh 1983 (crack band).
