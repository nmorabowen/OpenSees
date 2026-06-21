---
title: ContactDomain — broad-phase contact subsystem (NTS penalty, explicit-first, hybrid-capable)
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - implementation
  - subsystem
  - contact
  - explicit
  - implicit
---

# ContactDomain — broad-phase contact subsystem

> ADR 39. Supersedes the abandoned 2026-06-01 per-element `LadrunoContactNTS`
> scoping (plan doc + prototypes were never committed and are lost; the
> per-element value collapses into this subsystem's narrow phase). Decision
> reached 2026-06-21 after the ADR-30 gate cleared. See
> [[Ladruno_explicit_roadmap]] §5.7 + Q4.
>
> **Revised 2026-06-21 after a 4-lens adversarial design gate** (assembly
> architecture · contact mechanics · code reuse/footprint · scope/phasing). The
> gate confirmed the core architecture but falsified four ADR claims; all
> corrections are folded in below and flagged `[GATE]`. Findings + verdicts
> recorded at the bottom under *Adversarial review log*.
>
> **Strengthened 2026-06-21 against LS-DYNA Theory §26.7–26.12 + FEM first
> principles** (Ibrahimbegovic Ch 5, de Souza Neto Ch 6–7, Laursen/Wriggers): see
> *Theory grounding* under How — concrete `-kn auto`/SOFT formulas, the
> closest-point orthogonality result (force exact, only tangent lags), friction as
> a non-associative pressure-dependent radial return, the contact-energy
> diagnostic, segment-based bucket sort + the 10–15-cycle re-sort that validates
> epoch re-emit, and the Hertz/ironing/patch-test acceptance battery.

## What

A **contact engine parallel to `Domain`** that makes contact a *first-class
capability* — the user defines two surfaces, not every contact pair. Each step it
runs a **broad phase** (bucket sort) to find candidate pairs, then a **narrow
phase** (node-to-segment closest-point projection + penalty normal force +
Coulomb friction) that produces, per active pair, a contact **force** `F_c` and
(for the implicit leg) its **tangent** `K_c`.

**v1 ships on the EXPLICIT leg `[GATE]`:** bucket-sort broad phase + node-to-segment
(NTS) penalty contact + **IMPL-EX** Coulomb friction, injected as `F_c` into the
RHS. This is the robust-by-construction path and matches the research need
(pounding / uplift / sliding are all explicit per roadmap decision table). The
**hybrid plumbing is built and zero-force-validated in v1** (P1 runs the same
model under an explicit *and* an implicit integrator), but the **implicit
*frictional* Newton leg is a later gated phase, not a v1 ship gate** — see the
phasing table and Risk Q-FRIC. This honours "hybrid in v1" *architecturally*
without letting the least-de-risked part (frictional implicit convergence) hold
the ship.

**Why not "hybrid is free":** the force/tangent *contract* is genuinely shared
(explicit calls only `getResistingForce()`, implicit also calls
`getTangentStiff()` — verified in source, see Where). But **friction is not one
path** `[GATE]`: IMPL-EX is an explicit-flavoured linearisation; inside implicit
Newton it converges to a Δt-lagged, not the true, Coulomb equilibrium. The
reference element `ZeroLengthContactASDimplex` itself carries *two* friction
branches for exactly this reason. So the honest design is **one shared
kinematics + normal-penalty core, two friction constitutive branches**
(explicit IMPL-EX, implicit consistent return-map). The normal `K_c` is also only
*approximately* clean (the moving-normal `∂n/∂u` and projection `∂ξ/∂u` terms are
lagged — see Risk Q-TAN).

**NOT in v1 (later tiers / follow-up plans):**
- Implicit **frictional** Newton with a consistent (non-IMPL-EX) Coulomb tangent
  — phase P3.5, gated on the Δt-independent convergence test (Risk Q-FRIC).
- P4 soft-constraint penalty `SOFT=1` (Courant-stable stiffness) — but a *minimal
  auto-`k_n` default* is pulled INTO v1 (Risk Q-KN).
- P5 segment-based penalty `SOFT=2` (corner/edge/T-intersection).
- P6 **tied interfaces** (§26.9) — the one tier that rides on the ADR-30
  projection handler for its explicit leg.
- v2/v3: self-contact, eroding contact, mortar / segment-to-segment
  (Puso & Laursen 2004), MPI contact decomposition.
- **NTS does not pass the contact patch test** `[GATE]` (this is the textbook
  motivation for mortar/STS). v1 transmits resultants correctly but interface
  *stress fields* and large-sliding shear are mesh-sensitive. Quantitative
  interface-stress fidelity awaits mortar (v3). The v1 sell is "pounding + uplift
  kinematics + *qualitative* sliding," not quantitative sliding.

## Why

**The single largest capability gap between OpenSees and the rest of the
nonlinear-FEM world.** Without it none of the contact-driven research is
possible: pounding between adjacent buildings (dense-urban Latin-American
retrofit), foundation uplift and rocking with real toe contact, pile-soil gap
formation, soil-wall sliding, slab-on-soil separation, sliding bearings and
friction pendulums post-uplift, post-buckling member-to-member contact in
collapsing frames. The existing OpenSees contact (`ZeroLength`, `Contact2D/3D`,
`SimpleContact3D`) makes the **user pre-define every contact pair** — fine for
bench problems, broken for realistic geometry. A real `ContactDomain` with a
broad-phase search closes that gap the way LS-DYNA / Abaqus / Code_Aster do.
(Note `[GATE]`: those mature codes do NOT ship one hybrid contact engine —
LS-DYNA contact is explicit-only force-injection; Abaqus Standard vs Explicit are
*separate* contact implementations. Our hybrid is achievable because penalty
gives a clean residual in both, but it is sequenced explicit-first deliberately.)

**Why now (the gate cleared).** Contact was parked behind ADR 30
(`LadrunoProjectionHandler`) and the explicit-dynamics substrate. As of the
2026-06-21 rebase that substrate is **shipped**: ADR 30 P0–P6 (momentum-conserving
explicit constraint projection, tie-force query, whole explicit family); mass
scaling HRZ (ADR 35) + selective lumped/consistent/Olovsson + MPI (ADR 36/38) so
tiny contact-zone elements don't throttle global `dt`; the
`getExplicitCriticalTimeStep()` seam (#190).

**Scoping nuance that shaped the architecture:** the ADR-30 projection handler
serves only the **tied** tier (a *static* constraint set resolved at
`handle()`-time). Sliding penalty contact's active-set churns every step, so it
does **not** ride on the projection — it was never truly blocked. Waiting still
paid off (the tied tier P6 + the whole explicit / mass-scaling foundation
matured), but v1 penalty NTS depends on none of it.

## Where

- **New subsystem:** `OpenSees/SRC/domain/contact/` (parallel to
  `SRC/domain/domain/`):
  - `LadrunoContactDomain.{h,cpp}` — owns surfaces, broad-phase grid, active-pair
    list; queried each step; held by `Domain` (nullable → zero cost).
  - `LadrunoContactBucketSort.{h,cpp}` — spatial-hash broad phase (cell ≈ max
    element edge; §26.11). *Shared with SPH later (§5.8).* **Added AFTER the
    narrow phase is validated on brute-force pairs** `[GATE]` (see phasing).
  - `LadrunoContactSurface.{h,cpp}` — a meshed surface (master segments / slave
    nodes) from an element face-set or node-set.
  - `LadrunoContactPair` (pooled) — slave-node ↔ master-segment: gap, local frame
    `{n,t1,t2}`, friction state.
  - `LadrunoContactFE.{h,cpp}` — the **`FE_Element` adapter** with **conservative
    static connectivity** (see How): active pairs churn *inside* it (values only).
  - **`LadrunoContactHandler.{h,cpp}` — a custom `ConstraintHandler`** `[GATE]`
    that emits the contact `FE_Element`s during `handle()`. **This, not a direct
    `AnalysisModel::addFE_Element`, is the only supported injection path** (see
    How / Risk Q-WIRE). Mirrors the in-fork `LadrunoProjectionHandler` precedent.
  - `LadrunoContactFriction.{h,cpp}` — friction kernel with **two branches**
    `[GATE]`: `updateImplex()` (explicit) and `updateImplicit()` (consistent
    return-map). Built by re-deriving, NOT lifting (see reuse notes).
- **Reference (RE-DERIVE from, do not lift — quantified in the review log):**
  - `SRC/element/UWelements/SimpleContact3D.{h,cpp}` — the closest-point
    projection algebra is in `project()` (`:575-638`), `GetPoint()` (`:640-662`),
    `UpdateBase()` (`:665-684`), normal `n=g1×g2` (`:352`). **~110 of 1270 lines
    are reusable math** `[GATE]`; the rest is a 6-node/18-DOF Lagrange-multiplier
    assembly that does NOT transfer to penalty. Rewrite the projection as **pure
    functions** (it currently mutates element members).
  - `SRC/material/nD/UWmaterials/ContactMaterial3D.{h,cpp}` — its **consistent
    Coulomb return-map tangent** `C_ss = k·(g − r⊗r)`, `C_sl = μ·R` (`:303-369`)
    is the asset for the *implicit* friction branch `[GATE]`. But it is a
    Lagrange-multiplier law (no penalty normal of its own) → re-derive into a
    penalty form. Defects to avoid: `getInitialTangent` returns empty (`:372`),
    `mFrictFlag` is a `static` class global (`:47`).
  - `SRC/element/zeroLength/ZeroLengthContactASDimplex.{h,cpp}` — penalty-normal
    structure + IMPL-EX extrapolation for the *explicit* friction branch.
    **Confirmed live bugs to FIX (not "quarantine") `[GATE]`:** `recvSelf` reads
    `ddata(39)` from a size-31 Vector (`:613`, OOB → corrupts `gap0.y`); four
    `exit(-1)` in `setDomain` (`:269,280,296,308`); shadowed `dE1/dE2` dead-code
    in the analytic tangent (`:1085-1086` → friction-damage tangent terms silently
    zero). Also: it uses a process-wide `static std::map` for K/M/R buffers
    (`:80-84`) — the helper must NOT inherit that (re-entrancy under the
    per-segment adapter fan-out).
  - `SRC/element/ladrunoEmbeddedRebar/LadrunoEmbeddedKernel.{h,cpp}` — reuse the
    **scalar** helpers only `[GATE]`: `criticalTimeStep`, `massPenaltyDtcr`,
    `effectiveCouplingStiffness`, `maxAbsDiagonal` (`:66-80`). Its
    `computeGap`/`assemble*` are embedded-tie operators (`B=[I|−N_iI]`), NOT
    contact NTS kinematics — do not reuse those.
  - `SRC/analysis/handler/LadrunoProjectionHandler.{h,cpp}` — the in-fork custom
    ConstraintHandler precedent + the tied-tier (P6) enforcement path.
  - `SRC/analysis/fe_ele/FE_Element.{h,cpp}`, `SRC/analysis/model/AnalysisModel`,
    `SRC/analysis/handler/PlainHandler.cpp` — the assembly contract we adapt to.
- **Modify (vanilla — realistic footprint ~5–7 files / ~8–11 methods `[GATE]`,
  NOT "one hook + one pointer"; ledger + `// Ladruno`-tag each):**
  - `Domain.{h,cpp}` — own a `LadrunoContactDomain*` + accessors + dtor;
    participate in `hasDomainChanged()` when the contact topology epoch changes.
  - `Node.{h,cpp}` — a `contactForce` `Vector*` slot mirroring the ADR-30
    `projTieForce` precedent (`Node.h:236`, `Node.cpp:2201-2226`) for recording.
  - `AnalysisModel.{h,cpp}` — path for `LadrunoContactHandler` to add the contact
    FE_Elements (uses existing `addFE_Element`, `AnalysisModel.h:65`).
  - `DirectIntegrationAnalysis.cpp` (+ `StaticAnalysis.cpp`) — trigger
    `ContactDomain.update()` (broad/narrow phase) in the step loop.
  - Explicit + implicit integrator base(s) — one `newStep`/`domainChanged` hook
    for broad-phase refresh (watch the ExplicitBathe `firstStep`-reset caveat the
    ADR-30 work hit; subclasses inherit).
- **classTags:** `LadrunoContactFE` → **ELE_TAG 33015** `[GATE]` (next genuinely
  free ELE slot; **NOT 33018** — that number is live in the ND registry as
  `ND_TAG_LadrunoRCFiniteStrain`, `classTags.h:587`). Reserve 33016/33017 in the
  ELE band if surfaces/pairs ever need broker tags. The custom handler gets a
  HANDLER_TAG after `LadrunoProjectionHandler`=33001. Record in
  [[LEDGER_implementations.md]] at the first code PR.
- **Build:** new `SRC/domain/contact/CMakeLists.txt` + the new handler under
  `SRC/analysis/handler/`; no external dep. Touch
  [[../Ladruno_internal/01_compilation_journal]].

## How

### The hybrid architecture (what the gate confirmed, and what it corrected)

**Confirmed sound (verified in source):** one narrow phase computes `(F_c, K_c)`
per active pair, exposed through the standard `FE_Element` contract, consumed
differently per integrator family:

| Integrator family | Calls on the contact FE | Behaviour | Q4 |
|---|---|---|---|
| Explicit (CentralDifference*, ExplicitBathe*, LNVD, SMS) | `getResistingForce()` (+`addMtoTang`, which a contact pair returns as 0 mass) | `F_c` → RHS; `getTangentStiff()` **never called** | (a) force-injection |
| Implicit (Newmark, HHT, statics + Newton) | `getResistingForce()` **and** `getTangentStiff()` | `F_c` in residual, `K_c` in SOE, Newton iterates | (b) assembled |

The first row is **verified literally**: `CentralDifferenceLadruno::formEleTangent`
is mass-only (`addMtoTang` only, no `addKtToTang`), so the explicit LHS never
pulls `K_c`. The default `system Diagonal` cannot even store off-diagonal contact
coupling, so a wide-connectivity adapter costs nothing on the explicit LHS. The
SOE is built for fixed-pattern / varying-value assembly, so inactive pairs
(structural zeros in pre-allocated slots) are a no-op.

**Two corrections the gate forced:**

1. **Friction is two branches, not one** `[GATE]`. The shared core is
   *kinematics + normal penalty*. Friction splits:
   - **explicit:** `LadrunoContactFriction::updateImplex()` — IMPL-EX Coulomb
     (extrapolated state, well-conditioned, the v1 ship path).
   - **implicit:** `updateImplicit()` — a *true* consistent Coulomb return-map
     tangent re-derived from `ContactMaterial3D`'s `(g − r⊗r)` form (P3.5).
   IMPL-EX must NOT be the implicit path — it converges Δt-lagged (Risk Q-FRIC).

2. **Wiring + cost model** `[GATE]`. Contact FE_Elements are injected by a
   **custom `LadrunoContactHandler` (a `ConstraintHandler`)** during `handle()` —
   the ONLY supported re-emission point, because `domainChanged()` does
   `AnalysisModel::clearAll()` then `handle()` then `numberDOF()` then
   `setSize()`. A direct `addFE_Element` outside `handle()` is wiped by the next
   `clearAll()`. **Active-set churn is free (values only) ONLY IF adapter
   connectivity is a conservative static superset** of every node a pair could
   couple within that segment over a topology epoch. A genuine topology change
   (slave enters a new segment's territory) costs a full `domainChanged()`
   re-number — *not* a "cheap regroup." Minimise its frequency via epoch-based
   re-emission (broad-phase refresh cadence, Risk Q-EPOCH). If a pair activates
   whose DOFs are NOT in the adapter's frozen connectivity, its force is silently
   dropped (`addB` ignores out-of-`id` entries) — so the broad-phase neighbour
   bound MUST be conservative.

### The active-set / static-topology problem (corrected fix)

OpenSees's numberer/SOE assume a static element set fixed at `domainChanged()`.
The fix is **bucket-local `LadrunoContactFE` adapters with conservative static
connectivity**: one per master segment, connectivity = that segment's nodes ∪ all
slave nodes that *could* reach it within a topology epoch. Pairs going
active/inactive flip values inside a fixed `getID()` → free, no re-number. Epoch
boundaries (bodies slid far enough to change segment territory) trip a deliberate
`Domain::domainChange()` → one full rebuild, amortised. **Honest cost model**
`[GATE]`: re-number frequency scales with `v_slip·Δt / cell`; for sustained
sliding/impact this is not rare, so size epochs (and the per-segment-vs-per-bucket
granularity, Risk Q-GRAN) against the *sliding* worst case, not the static one.

### Per-step data flow (explicit ship path)

```
integrator.newStep()
  └─> ContactDomain.update(currentDisp)            // once per step
        ├─ bucketSort.refresh()                    // rebuild grid if moved > ½ cell
        ├─ for each candidate (slaveNode, masterSeg):
        │     project → gap g, local frame {n,t1,t2}     // re-derived SimpleContact3D math
        │     if g < 0: pair ACTIVE
        │         F_n = -k_n·g·n                          // penalty normal (auto k_n default)
        │         (F_t)   = friction.updateImplex(slip, F_n)  // IMPL-EX Coulomb
        │     else: INACTIVE (F=0)
        └─ (epoch boundary only) re-emit adapters via handler → domainChange()

// standard assembly, integrator-agnostic:
LadrunoContactFE.getResistingForce()  -> Σ active-pair F_c     // both families
LadrunoContactFE.getTangentStiff()    -> Σ active-pair K_c     // implicit only (P3.5)
LadrunoContactFE.getExplicitCriticalTimeStep() -> 2·√(m_slave/k_n)  // CFL seam (#190)
```

### Public API (Tcl / Python)

```tcl
contactSurface 1 -faces  $masterEleTags ?-side $f?    ;# master (segments)
contactSurface 2 -nodes  $slaveNodeTags               ;# slave  (nodes)

contact 1 -master 1 -slave 2 \
        ?-kn auto|$kNormal?  ?-kt $kTang?  ?-mu $friction? \
        ?-formulation nts?  ?-softness 0?            ;# -kn auto = v1 default rule
constraints LadrunoContact                           ;# the custom handler (implicit + explicit)

recorder Contact ... -surface 1 contactPressure|gap|slip|status   ;# own source (Node::contactForce)
ops.contactForce(pairOrSurface)
```

- `constraints LadrunoContact` is mandatory (it's the injection handler). It must
  compose with the structural model's other constraints (it wraps/delegates to a
  base handler for the non-contact constraints) — **test interaction with
  `rigidDiaphragm`/`equalDOF`** (Risk Q-CONSTR).
- `-kn auto` (Risk Q-KN): a minimal mass-and-`dt` based default so the user is not
  left to guess `k_n` with only a `dt_cr` warning.

### Theory grounding (LS-DYNA Theory §26 + FEM first principles)

Consulted the LS-DYNA Theory Manual §26.7–26.12 and FEM-mechanics references
(Ibrahimbegovic Ch 5, de Souza Neto Ch 6–7, Laursen/Wriggers) to pin the
formulas and the failure modes. Key results baked into the design:

**Penalty = regularised Signorini/KKT.** Exact unilateral contact is the KKT set
`gₙ ≥ 0` (no penetration), `pₙ ≤ 0` (compression only, no adhesion), `pₙ·gₙ = 0`
(complementarity). Penalty replaces it with `pₙ = εₙ⟨gₙ⟩₋` (a one-sided spring) —
penetration `gₙ = −pₙ/εₙ = O(load/εₙ)`. **The accuracy↔conditioning trade-off is
fundamental**: large `εₙ` → small penetration but stiffer SOE + smaller `dt_cr`;
small `εₙ` → soft, excessive penetration. Augmented Lagrangian (Uzawa) recovers
the *exact* constraint at *finite* `εₙ` — the principled v2 upgrade, and we
already have AL machinery in `LadrunoEmbeddedKernel`.

**`-kn auto` has a concrete formula `[LS-DYNA 26.14a]`:** for a solid master
segment, `kₙ = f_si · K·A²/V` (K = bulk modulus, A = segment face area, V =
element volume), `f_si` default **0.10** — larger needs `dt` scaled back. The
default takes `min(master, slave)` stiffness and **fails when one side is very
soft/low-density** → fall back to the SOFT rule. **SOFT=1 stability stiffness
`[26.15]`:** `k_cs = ½·SOFSCL·m*/Δt_c²` (Courant-stable by construction, `m*` from
slave+master mass, `Δt_c` reset if the solution step grows); ship `k = max(k_cs,
k_penalty)` as the v1 `-kn auto` default so the user is never left guessing.
SOFT=2 segment-based `[26.17]` uses the reduced mass `m₁m₂/(m₁+m₂)` and updates
`Δt_c` only on >5 % growth — P5.

**The contact FORCE is first-order exact with a frozen projection `[GATE/Q-TAN
refinement]`.** The normal gap `gₙ = (x_s − x̄(ξ̄))·n` with `x̄` the closest-point
projection. Because at the projection `(x_s − x̄)·g_i = 0` (orthogonality), the
first variation `δgₙ = n·(δu_s − Σφᵢ δuᵢ)` — **`∂ξ̄/∂u` drops out of the
residual**. So the explicit ship path (force injection) is *exact* even though ξ̄
is computed once per step. Only the **tangent** picks up the dropped geometric
terms: the consistent tangent is `kₙ Bᵀ(n⊗n)B` (main, always present) **plus**
curvature/rotation terms `O(gₙ·κ)` (segment curvature κ, `∂n/∂u`). Those vanish
for flat segments (κ=0) and small penetration → **the lagged tangent is nearly
exact for flat-ish, stiff-penalty, small-rotation contact, and degrades only on
curved/rocking/large-sliding**. This sharpens Q-TAN and *reinforces* shipping
explicit-first (the force path needs none of the dropped terms).

**Friction = pressure-dependent J2 radial return (non-associative) `[26.19–26.23
/ de Souza Neto §7]`.** Cone `f = ‖t_T‖ − μ|t_N| ≤ 0`. Stick: `t_T = ε_T·g_T`.
Slip: radial return `t_T = μ|t_N|·t*_T/‖t*_T‖` — *identical in structure to 1D
plasticity* (trial elastic tangential force → check cone → return). The
LS-DYNA "elastic-plastic spring" IS this return map. **Non-symmetry is rigorous,
not incidental**: the slip flow is purely tangential while the cone normal has a
`t_N` component → non-associative → the consistent slip tangent
`k_T = (μ|t_N|/‖t*_T‖)[I_T − t̂_T⊗t̂_T] + μ·sign(t_N)(t̂_T⊗n)` has a
**non-symmetric** coupling block (the last term). This is exactly
`ContactMaterial3D`'s `C_ss = k(g − r⊗r)`, `C_sl = μR` — confirming the gate's
reuse target for the implicit branch. Optional `[26.24]`: velocity-dependent
`μ = μ_d + (μ_s − μ_d)e^{−c|v|}`; shear cap `[26.26]` `f ≤ κ·A_master` (so soft
materials can't carry unphysical interface shear).

**Contact energy is the master diagnostic `[26.8.4]`:** track
`E_contact += Σ ΔF_slave·Δd_slave + Σ ΔF_master·Δd_master` incrementally.
Frictionless → `E_contact` bounded and equals stored spring energy (slave/master
opposite sign); friction → positive dissipation. **Large NEGATIVE `E_contact` =
undetected initial interpenetration** — the single most useful health signal.
Wire it into the `EnergyBalanceRecorder` and make it a gate.

**Bucket sort, the right way `[26.11]`:** cell size `LMAX` = fraction of the
segment-diagonal characteristic length; cap `NX·NY·NZ ≤ min(NSN, 5000)`. **Search
*segments* by centroid, not nearest *node*** — nearest-node fails on non-uniform
/ poor-aspect meshes (the nearest nodes need not even belong to the nearest
segment, Fig 26.24). **Guard runaway nodes**: an unstable node flying off inflates
the bbox → `LMAX` huge → buckets collapse → search degrades (the historic "1000
nodes in a bucket"). **LS-DYNA re-sorts only every 10–15 cycles** (stores the 3
nearest, re-searches a node when contact is lost) — this is the canonical
amortisation and **directly validates our epoch-based re-emit (Q-EPOCH)**: do NOT
re-search every step.

**Tied interfaces `[26.9]`:** constrain *slave only*, coarse side = master;
distribute slave force **and mass** to the master, then interpolate slave accel
`a_s = Σφⱼ aⱼ`; contact point computed **once** (static set; `|ξ| > 1.02` → leave
unconstrained). **Conflicting-constraint detection is mandatory** (a node in two
tied sets, or also in an MP/EQ constraint or spot weld → hard error). This is
exactly the static-set, momentum-distributing structure our P6 rides onto the
ADR-30 projection, and the conflict refusal aligns with ADR-30's chain refusal.

### Phasing (each phase = its own PR, ADR-30-style phased loop + adversarial gate)

| Phase | Deliverable | Gate / reference solution |
|---|---|---|
| **P0** | Falsify & baseline, **no SRC**: reproduce a pounding + an uplift problem with existing `ZeroLengthContactASDimplex` (pre-defined pairs); standalone bucket-sort prototype vs brute-force O(n²). | pair-list == brute force; penetration = P/k_n |
| **P1** | `LadrunoContactDomain` + `LadrunoContactSurface` + `LadrunoContactFE` adapter + `LadrunoContactHandler`, fed by a **brute-force pair list (NO bucket sort yet)** `[GATE]`; adapter returns **zero** force. Proves the hybrid plumbing + the handler injection. | same model runs under CentralDifferenceLadruno **and** Newmark, contact-defined-but-inactive, **byte-identical to no-contact**; **rigid-body translate both bodies → ΣF_c=0, ΔE=0** `[GATE]` |
| **P2** | NTS closest-point projection + **standard penalty** normal force (frictionless), `F_c` (+ `K_c` for the implicit assembly), **still on brute-force pairs**. First sub-rung = **slave-nodes vs a rigid analytical plane** `[GATE]`, then 2 deformable blocks. **Initial-interpenetration scan + warn/abort at setup** `[GATE]`. | penetration = P/k_n exact (explicit + implicit); **`K_c` vs FD on a ROTATED/curved config** `[GATE]`; **frictionless release/reopen → exact return to F=0** `[GATE]`; **Hertz cylinder/sphere-on-halfspace → analytic `p(r)=p₀√(1−(r/a)²)`, peak `p₀` + radius `a` converge under refinement** `[GATE]`; **contact patch test — report the NTS oscillation magnitude vs mesh ratio + confirm resultant equilibrium EXACT (known limitation, not pass/fail)** `[GATE]`; **`E_contact` bounded ≥ 0 frictionless, NEGATIVE ⇒ flag init-penetration** `[GATE]` |
| **P2.5** | Swap brute-force pair generation for **bucket sort** (cell=LMAX, cap `NX·NY·NZ ≤ min(NSN,5000)`, segment-centroid search, runaway-node guard). | active-pair set == brute force; **runaway-node bbox guard fires** `[GATE]`; regroup-cost profile under **sustained sliding** (>5% rule) `[GATE]` |
| **P3 — v1 SHIP (explicit)** | **IMPL-EX Coulomb friction** (`updateImplex`), tangential `F_t` (radial return on the `‖t_T‖≤μ|t_N|` cone; optional shear cap `κ·A`). | stick: no slip; **sliding-block-on-incline vs analytic `a=g(sinθ−μcosθ)` + static for `tanθ<μ`** `[GATE]`; **ironing/sliding-patch (block pressed + dragged across a non-matched master) — frictional shear resultant + moving-normal stability** `[GATE]`; **1D central impact, measured restitution `e(k_n,dt)` reported, elastic `e≈1`** `[GATE]`; **oblique frictional impact vs analytic Coulomb impulse** `[GATE]`; **`E_contact` = frictional dissipation (positive)** `[GATE]`; explicit 500-step stable; **friction-state db send/recv roundtrip** `[GATE]`; **contact on a rigidDiaphragm/equalDOF node + conflicting-constraint refusal** `[GATE]` |
| **P3.5** | Implicit **frictional** Newton: `updateImplicit()` consistent return-map tangent (re-derived from ContactMaterial3D). | Newton converges; **converged answer matches reference Coulomb to a tolerance INDEPENDENT of Δt** (IMPL-EX fails this — that's the point) `[GATE]` |
| P4 | soft-constraint penalty `SOFT=1` (full Courant-stable `k_n`). | `dt_cr` not throttled by `k_n` |
| P5 | segment-based penalty `SOFT=2`. | corner/edge/T-intersection battery |
| P6 | **tied interfaces** (§26.9) via ADR-30 projection (explicit) / `MP_Constraint` (implicit). | tie exact + momentum-conserving (reuse ADR-30 tie battery) |

### Minimum viable smoke test (P2)

Slave node-set vs a **rigid analytical plane**, pushed by load `P`: assert
penetration `= P/k_n` to ~1e-6 under **both** `CentralDifferenceLadruno` and
`Newmark` — the single test that proves the hybrid path on the cleanest
kinematics. *Then* graduate to two deformable `LadrunoBrick` blocks. Reference =
linear penalty spring (rigid surface removes master compliance, so the assertion
is unambiguous).

## Risks / open questions

> [!question] Q-FRIC (was "hybrid is free")
> IMPL-EX friction inside implicit Newton converges to a **Δt-lagged** Coulomb
> state, not the true equilibrium (can't be fixed by tightening Newton tol). The
> reference `ZeroLengthContactASDimplex` proves this — it uses a numerical
> *implicit* tangent for `!use_implex`, IMPL-EX only for explicit. **Resolution:**
> two friction branches (explicit IMPL-EX ships in P3; implicit consistent
> return-map is P3.5, gated on a Δt-independent convergence test).

> [!question] Q-TAN (normal tangent is not "clean" — but the FORCE is exact)
> Refined by the theory grounding: by the closest-point **orthogonality
> condition**, the contact *force* (residual) is first-order exact with a frozen
> projection ξ̄ — so the **explicit ship path needs none of the dropped terms**.
> Only `K_c` (implicit) drops the `O(gₙ·κ)` curvature/`∂n/∂u` terms. So v1 `K_c`
> is a lagged geometric tangent: nearly exact for flat-ish/stiff-penalty/small-
> rotation, degrades Newton on curved/rocking/large-sliding. The P2 `K_c`-vs-FD
> gate MUST run on a rotated/curved config (a flat crush passes even with `∂n/∂u`
> missing → false confidence). This *reinforces* explicit-first.

> [!question] Q-AL (exact constraint at finite penalty — v2)
> Penalty penetration is `O(load/εₙ)` and never zero; the accuracy↔conditioning
> trade-off is fundamental. Augmented Lagrangian (Uzawa) recovers the exact
> non-penetration at *finite* `εₙ` and we already have the AL machinery in
> `LadrunoEmbeddedKernel`. Out of v1 (penalty ships), but the planned upgrade for
> precision contact / quasi-static.

> [!question] Q-WIRE (injection mechanism)
> Injection is a custom `ConstraintHandler` (emit FE_Elements in `handle()`), NOT
> a direct `AnalysisModel::addFE_Element` (wiped by `clearAll()` in
> `domainChanged()`). Confirmed against `DirectIntegrationAnalysis.cpp:406-473`.

> [!question] Q-EPOCH / Q-GRAN (cost model + granularity)
> Adapter connectivity must be a conservative static superset so pair churn is
> values-only; topology change pays a full re-number, frequency `~ v_slip·Δt/cell`.
> Per-segment (fine sparsity, more adapters) vs per-bucket (coarser, denser
> blocks) — resolve at P2.5 by profiling regroup + implicit-fill cost under
> *sustained sliding*, not the static case.

> [!question] Q-KN (penalty stiffness selection)
> Shipping with no `k_n` guidance is a usability trap. v1 `-kn auto` =
> `f_si·K·A²/V` (LS-DYNA 26.14a, `f_si`=0.10) with the **SOFT=1 Courant floor**
> `k = max(½·SOFSCL·m*/Δt_c², f_si·K·A²/V)` (26.15) so a soft/low-density side
> can't collapse the stiffness — on top of the `getExplicitCriticalTimeStep` seam
> + a hard error if `dt_cr` drops below a floor. Full SOFT=1/2 controls are P4/P5.

> [!question] Q-CONSTR (constraint-handler interaction)
> The fork's flagship cases (building pounding, foundation rocking) run with
> `rigidDiaphragm`/`equalDOF`. A contact slave that is also constraint-controlled
> is a real configuration. `LadrunoContactHandler` must compose with the base
> constraint handling; P3 gates a contact force landing on a diaphragm/equalDOF
> node under the transformation handler.

> [!question] Q-PATCH (NTS limitation)
> NTS fails the contact patch test (textbook; the motivation for mortar). v1
> transmits resultants but interface stress / large-sliding shear are
> mesh-sensitive. Disclosed in scope; P2 *documents* the error magnitude rather
> than asserting a pass. Mortar = v3.

- **CFL:** `k_n` sets explicit `dt_cr`; self-report via `getExplicitCriticalTimeStep()`
  (#190) so `ops.criticalTimeStep` honours it.
- **Parallel (MPI):** OUT of v1 — refuse partition-crossing pairs with a named
  error (mirrors ADR-30's partition-crossing refusal).
- **Backwards compatibility:** null `LadrunoContactDomain*` ⇒ byte-identical to
  stock; opt-in via the `contact*` commands + `constraints LadrunoContact`.
- **Vanilla footprint:** ~5–7 files / ~8–11 methods (Domain, Node, AnalysisModel,
  Direct/Static analysis, integrator base(s), broker/registration). Every edit
  ledgered + `// Ladruno`-tagged. Keep hot-path edits to single hooks.

## Implementation log

*(filled in once execution starts; P0 first — no SRC change.)*

---

## Adversarial review log (2026-06-21, 4 parallel reviewers, grounded in source)

**Verdict: SALVAGEABLE-WITH-CHANGES — core architecture sound, four claims
corrected (all folded in above).**

**Confirmed sound:** (1) explicit integrators never form `K_c`
(`CentralDifferenceLadruno::formEleTangent` is mass-only) → the force/tangent
hybrid contract is real. (2) the SOE is built for fixed-connectivity/varying-value
assembly → inactive pairs are free. (3) strategic direction (contact = largest
gap, prerequisites cleared) confirmed. (4) `Node::contactForce` recorder reuse is
honest. (5) the bucket-local-adapter idea is the right active-set fix.

**Corrected (BLOCKER/MAJOR):**
- **B1 [Q-FRIC]** "one friction path both worlds" — FALSIFIED by the reference
  code itself (ASDimplex has two branches). → split friction; ship explicit-first;
  implicit frictional Newton → P3.5 with a Δt-independent gate.
- **B2 [Q-WIRE/Q-EPOCH]** "re-emit on membership change = cheap regroup, stable
  numbering" — FALSE; any connectivity change forces a full `domainChanged()`. →
  conservative static connectivity (churn = values-only) + epoch-amortised
  re-number; injection via a custom `ConstraintHandler`, not direct
  `addFE_Element`.
- **M1 [footprint]** "one hook + one pointer" — understated; realistic ~5–7
  files / ~8–11 methods. → honest enumeration in Where.
- **M2 [reuse]** "reuse SimpleContact3D/ASDimplex" — optimistic; it's re-derive
  (~110/1270 lines of projection math) + fix three confirmed bugs (`ddata(39)`
  OOB `:613`, `exit(-1)` ×4, shadowed `dE1/dE2` `:1085`). → reframed + bug list +
  lift the consistent tangent from ContactMaterial3D, not the buggy ASDimplex one.
- **M3 [mechanics]** NTS patch-test failure undisclosed (Q-PATCH); normal `K_c`
  lagged (Q-TAN); ship-without-SOFT usability trap (Q-KN); P3 energy gate
  hand-wavy → restitution-pinned impact + analytic Coulomb. All folded into scope
  + the phasing gates.
- **M4 [scope]** bucket-sort-first inverted → narrow phase on brute-force pairs
  first (P1/P2), bucket sort P2.5; add rigid-surface first rung; add missing gates
  (patch test, release/reopen, rigid-body zero-energy, sliding-incline,
  db-roundtrip, rigidDiaphragm interaction).
- **m1 [classTag]** 33018 collides in the ND registry → **ELE 33015**.
