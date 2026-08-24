---
title: "ADR 71 — LadrunoUP: Biot u-p saturated-porous continuum family (one element, shape providers, shared kernel)"
status: SHIPPED — P0–P4 complete; ADR-78 `-geom corot` (#677); `ELE_TAG_LadrunoUP` 33017 live since #557
---

# ADR 71 — LadrunoUP: the u-p (solid displacement + pore pressure) continuum family

**Status:** SHIPPED — P0–P4 complete, plus ADR-78 `-geom corot` (#677); `ELE_TAG_LadrunoUP`
33017 has been live since #557. See `[[LEDGER_implementations]]` for the per-phase PR trail.
Opens the saturated-porous-media
axis for the whole fork continuum family (LadrunoQuad/CST geometry class, LadrunoBrick,
BezierTri6/BezierTet10, future LadrunoLST/Q9/H20). Companion theory for this ADR:
`/quake-research` (liquefaction/SSI framing), `/abaqus-theory` (coupled pore-fluid
elements C*P), `/kratos` (U-Pw architecture cross-check), `/fem-mechanics-expert`
(mixed methods, inf-sup).

---

## 1. Context & problem

### 1.1 Why the fork needs this

The explicit roadmap names **liquefaction with discontinuities** (PM4Sand/PDMY-class
post-triggering, gap formation, sand boils, lateral spreading) as a strategic target,
and the DRM/absorbing-boundary + monitor stack already handles the wave-input side.
What is missing is the *element physics*: today the fork's continuum family is
single-phase. Every drainage-sensitive problem — liquefaction triggering (r_u),
consolidation settlement vs time, staged construction on soft ground, rapid drawdown,
partially-drained seismic response, saturated free-field/DRM columns — needs the
**Biot u-p mixed formulation**: solid displacement u plus a pore-pressure nodal field
p, with drained vs undrained behavior *emerging* from the race between loading rate
and permeability instead of being assumed per-analysis.

### 1.2 Upstream prior art (the review)

Upstream ships **eight** u-p elements in two families, all replicated per geometry:

| Family | Elements | Pattern |
|---|---|---|
| `UP-ucsd` (Elgamal/Yang) | FourNodeQuadUP (40), BBarFourNodeQuadUP, Nine_Four_Node_QuadUP, BrickUP (46), BBarBrickUP, Twenty_Eight_Node_BrickUP | full-integration twins of Quad/Brick; BBar twins for dilatant/undrained; 9-4 and 20-8 **mixed-order** (quadratic u, corner-linear p) |
| `UWelements` (McGann/Arduino) | SSPquadUP (120), SSPbrickUP (122) | single-point stabilized equal-order, physical stabilization of the pressure block |

The lesson upstream teaches is exactly the disease ADR-70 diagnosed for finite
strain: **N geometries × M formulations = N·M hand-copied classes** that drift.
Upstream needed 6+ classes for two geometries and two solid formulations; our family
has five-plus geometries. Detailed per-element findings (matrix-block placement,
argument/unit conventions, quirks worth keeping or fixing) are in §3.2–3.5.

### 1.3 Fork element facts that shape the choice

(From the seam review of this tree, branch state 2026-07-10.)

- **Hardcoded external DOF counts everywhere**: `getNumDOF()` returns literal
  8/6/24/12/30 (LadrunoQuad.cpp:162, LadrunoCST.cpp:134, LadrunoBrick.cpp:351,
  BezierTri6.cpp:280, BezierTet10.cpp:257); class-level **static scratch** K/P/mass
  sized to match; assembly loops interleave `2a+i`/`3a+i`.
- **Parser + setDomain ndf gates**: Quad/CST parsers hard-require `ndf==2`
  (OPS_LadrunoQuad.cpp:39–44); setDomain checks reject any node whose ndf ≠ the
  element's expectation — i.e. they would reject the very nodes carrying a pressure
  DOF.
- **No public Gauss/B API** on Quad/CST/Tri6 (only LadrunoBrick/BezierTet10 expose
  `getInterpolationWeights/Gradients`, reference-config only).
- **Formulation flags are if/else-in-method**, multiplying axes
  (`-formulation std|bbar|ssp|eas` × `-geom` × hourglass × mass) that are
  meaningless or unvalidated under u-p.
- The **liftable part is small and stateless**: shape functions + Jacobians
  (Quad ~57 lines cpp:1321–1377; CST cpp:655–688; brick via shared `shp3d.h`;
  Bernstein Tri6 cpp:851–896; Tet10 cpp:1177–1236). Mean-dilatation B-bar is ~30
  lines per dimension and already written twice in the family.
- **Heterogeneous ndf is framework-legal**: upstream `Nine_Four_Node_QuadUP` runs
  corner ndf=3 / interior ndf=2 in one element — the door to Taylor–Hood pairs on
  our quadratic Bézier shapes.
- The ADR-70 idiom — **pure header kernel (`namespace`, raw doubles, runtime node
  count) + numpy oracle + standalone g++ check** — is proven and element-count-
  agnostic (`LadrunoFiniteStrain2DKernel.h`).

### 1.4 Why this is NOT ADR-59's fork-the-element case

ADR-59 (gradient concrete) also adds an extra scalar nodal field, and it forked the
element (`LadrunoGradientBrick`) — but only because its coupling **contaminates the
constitutive chain** (the nonlocal field lives *inside* the damage law; `K_uē` had to
be re-derived from kernel internals). Biot u-p is the opposite: the **effective-stress
principle keeps the material seam untouched** (`σ = σ′ − α·m·p`; the material sees
only σ′), and the coupling blocks Q/H/S are pure geometry + scalar parameters. The
u-p mechanics are *identical for every geometry* — only shape functions and GP rules
differ. That is precisely the shared-kernel shape (ADR-70's argument), not the
fork-the-element shape.

---

## 2. Decision (summary)

1. **One new element class `LadrunoUP`** (`ELE_TAG_LadrunoUP = 33017`), covering the
   whole geometry family through internal, stateless **shape providers**:
   v1 ships T3, Q4, H8 (equal-order linear p) and Bézier T6, Tet10 (Taylor–Hood:
   quadratic Bernstein u, vertex-linear p). `(ndm, nodeCount)` dispatches the
   provider; a `-shape` keyword is reserved for the first genuine collision —
   which is already on our own roadmap: **(2,6) Bernstein T6 (v1 default) vs the
   future Lagrange-LST T6 provider** (§4.2) — existing scripts keep the Bernstein
   default, `-shape lst` selects the newcomer. (A wedge is (3,6) — *not*
   ambiguous with T6.) **Existing elements are not touched.**
2. **A pure, shared `LadrunoUPKernel.h`** (namespace `ladruno_up`, raw doubles,
   runtime node counts — the ADR-70 idiom), assembling the Biot blocks
   Q (coupling), H (permeability), S (storage), the stabilization block for
   equal-order pairs, and the fluid body-force vector, from the same per-GP tuple
   every element already produces (N, dN/dx, w·detJ) plus the pressure basis.
   The provider interface names **which nodes carry p** as an explicit field (the
   TH corner subset; = all nodes for equal-order) — the H20 corner-p promise in
   §4.2 rides on this being first-class, not implied.
   Unit-tested standalone (g++) against a numpy oracle **before** any OpenSees glue.
3. **Pressure-interpolation is its own axis** (not welded to a geometry):
   equal-order + stabilization on linear shapes (SSP lineage), Taylor–Hood
   vertex-linear p on quadratic shapes (9-4/20-8 lineage, Abaqus C3D20P-style),
   selected by `-p equal|linear` with per-shape defaults.
4. **Solid part v1 = `-formulation std|bbar`** (mean-dilatation B-bar — the
   BBar*UP precedent; essential near-undrained), `-geom linear` only. EAS/URI/
   hourglass/corot/finite stay out of v1 (finite-strain u-p is a research axis,
   reserved phase). *Post-v1 amendment: `-geom corot` OPENED by ADR 78
   (3D lanes; large rotation + geometric stiffness, small material strain —
   the TIMs bearing-mechanism driver). `finite` stays reserved (material-side
   blocker, ADR 78 §1.1).*
5. **Honest pressure-DOF contract** (the load-bearing divergence from upstream,
   §3.2): the nodal DOF *is* p — upstream's DOF is ∫p·dt with p recorded as a
   *velocity*. Ours records/constrains/initializes p directly, matches
   Abaqus/Kratos semantics, and makes **static analysis well-posed** (drained
   solid + steady seepage) — upstream is structurally transient-only. Price:
   unsymmetric effective tangent ⇒ general solver (UmfPack/SuperLU/MUMPS SYM=0).
6. **The kernel carries a physics-mode seam**: Biot u-p is **mode 1**;
   **incompressible hybrid u/p** (pressure as volumetric Lagrange multiplier — the
   ADR-26/70 frontier item) is **reserved mode 2**, sharing Q, the pressure-basis
   providers, stabilization, and DOF plumbing, differing only in the p-row law and
   the deviatoric material contract. Designed-for now, built later (§7 phase table).
7. **Out of scope**: LadrunoSolidShell (covariant ANS/EAS kinematics — through-
   thickness consolidation is a different problem), axisymmetric (future provider
   axis), explicit u-p time integration in v1 (honest research phase, §6/§8 —
   the zero-mass pressure rows break the current explicit stack's assumptions;
   near-term explicit answer for undrained problems remains material-level
   `FluidSolidPorousMaterial`).

---

## 3. Formulation (grounded in theory)

### 3.1 Governing equations and semi-discrete form

Biot's mixture, u-p reduction (Zienkiewicz & Shiomi 1984; Zienkiewicz et al. 1999
Ch. 2–3): neglect the relative fluid acceleration ẅ (valid in the earthquake
frequency band; see §3.6), keep solid displacement **u** and pore pressure **p** as
fields. With tension-positive solid stress and the effective-stress split
σ = σ′ − α·m·p (m = [1,1,1,0,0,0]ᵀ):

- momentum (mixture):  M·ü + C·u̇ + ∫BᵀσʹdV − Q·p = f_u
- mass balance (fluid): Qᵀ·u̇ + S·ṗ + H·p = f_p

with

| Block | Definition | Meaning |
|---|---|---|
| M | ∫ ρ NᵤᵀNᵤ dV, ρ = (1−n)ρ_s + n ρ_f (saturated, from the material) | mixture mass |
| Q | ∫ Bᵀ α m N_p dV | volumetric coupling |
| H | ∫ (∇N_p)ᵀ k̄ ∇N_p dV, k̄ = k_hydraulic/γ_w | Darcy permeability (Laplacian) |
| S | ∫ N_pᵀ (1/Q̄) N_p dV, 1/Q̄ = n/K_f + (α−n)/K_s | storage/compressibility |
| f_p | ∫ (∇N_p)ᵀ k̄ ρ_f (b − ü) dV + boundary flux | seepage drive (body force + **dynamic seepage**) |

α = 1 − K_skel/K_s is Biot–Willis (≈1 for soils, down to ~0.5 for rock; upstream
hardwires α=1, K_s→∞ — our default, overridable). Storage 1/Q̄ confirmed verbatim
against the primary sources (Book 2nd ed. eq 2.17 p.22; ZS84 eq 4a + Appendix II).
At full saturation the book's two coupling matrices coincide (Q̃ = Q, eqs
3.26/3.29) — partial saturation would split them (χ_w weighting + 1/Q* with S_w,
Book eq 3.33); noted for a future `-unsat` axis, not v1.

The **dynamic-seepage** −k̄ρ_f·ü part of f_p is retained in the continuum u-p form
(Book eq 2.21; ZCB80 eq 11) but **dropped from the LHS in the book's recommended
discretization** — "its inclusion … will render the final equation system
nonsymmetric (Leung 1984)"; Chan (1988) measured the omission and found it
**insignificant**. SWANDYNE-II keeps it in the *force* term (Book p.56).
**Upstream OpenSees drops the term entirely** — the "dynamic seepage force" in
FourNodeQuadUP/BBarFourNodeQuadUP is *commented-out dead code*
(FourNodeQuadUP.cpp:824–831; no live u-p element in-tree assembles it)
⟨adversarial review UP-3/FW-F2⟩. We adopt the SWANDYNE treatment as a
**deliberate divergence**: residual force yes (evaluated with the *trial/relative*
acceleration — under UniformExcitation the ρ_f·üg seepage drive is knowingly
omitted, the same Chan-1988 class), tangent no; a `-dynSeepage off` switch exists
for upstream-parity runs (the P1 equivalence gate uses it), and the P0 oracle
FD-checks the omitted tangent so it stays *measured*, not assumed. The term's own
physics is gated by the B2/B5 analytic benchmarks, not by upstream parity.

K is the ordinary solid stiffness on σ′ — **the material library is reused
verbatim**, including the UCSD effective-stress soil models (PressureDependMultiYield
et al., confirmed present in `SRC/material/nD/soil/` with `updateMaterialStage`).

### 3.2 Block placement — the upstream contract, and the one we choose

**Upstream (verified in code, all 8 elements):** a *symmetric-negative* arrangement
whose price is a semantic trick — the nodal "pressure" DOF is the **time-integral
of p**; the physical pore pressure is the DOF's **rate**:

| API | upstream content |
|---|---|
| `getTangentStiff()` | [ K , 0 ; 0 , **0** ] — pressure rows/cols structurally zero |
| `getDamp()` | [ C_Ray , **−Q** ; **−Qᵀ** , **−H** ] (FourNodeQuadUP.cpp:519–580; SSPquadUP.cpp:429–488) |
| `getMass()` | [ M , 0 ; 0 , **−S** ] (FourNodeQuadUP.cpp:582–638) |
| `getResistingForce()` | [ ∫Bᵀσ′dV − ρb-terms ; +fluid-gravity ] — **no −Q·p in the static residual**; the pore-pressure force on the skeleton enters *only* via damp·vel |

Because both Q blocks sit in one matrix, the effective transient tangent stays
symmetric — that was the design goal. Consequences upstream lives with: pore
pressure is recorded with the **`vel`** node-recorder option; fixing the DOF is a
drained BC but *prescribing a nonzero p means prescribing a velocity*; the DOF
itself is ∫p·dt and **grows without bound**; and a `Static` integrator yields a
**structurally singular system** (K's p-rows are zero) — everything, including
"static" consolidation/gravity stages, must run under a TransientIntegrator
(UCSD practice: Newmark γ ≥ 0.6, large Δt).

**LadrunoUP decision — honest p-DOF (deliberate divergence):** the nodal DOF **is
p**. Blocks split by the time-derivative they multiply:

| API | LadrunoUP content |
|---|---|
| `getTangentStiff()` | [ K , −Q ; 0 , H ] |
| `getInitialStiff()` | [ K₀ , −Q ; 0 , H ] — K₀ from material initial tangent; **p-rows must NOT be zero** (a family-idiom copy of upstream's solid-only Ki makes initial-stiffness statics structurally singular); cached, invalidated by `updateMaterialStage` / permeability `setParameter` |
| `getDamp()` | [ C_Ray , 0 ; Qᵀ , S (+ α·L stab, §3.3) ] |
| `getMass()` | [ M , 0 ; 0 , 0 ] |
| `getResistingForce()` | [ ∫Bᵀσ′dV − Q·p − body ; H·p − f_seep ] — **no rate terms here** (static-path correctness depends on it) |
| `getResistingForceIncInertia()` | getResistingForce() **+ M·ü + getDamp()·v** — i.e. + [M·ü + C_Ray·u̇ ; Qᵀ·u̇ + (S+H̃)·ṗ], rates from trial nodal vel/accel |
| `addInertiaLoadToUnbalance()` | −M·(R·üg) with the p-slots of R·üg zeroed (belt-and-braces; M's p-rows are already zero) |

**Element-contract obligations (adversarial-review hardened).** Three framework
facts make the extra rows load-bearing ⟨UP-1/UP-2/FW-F3/FW-F5/FW-F6⟩:

1. **Integrators never assemble damp·vel or mass·accel into the residual** —
   `getDamp()`/`getMass()` feed only the effective tangent (Newmark.cpp:284–309);
   the transient unbalance comes solely from `getResistingForceIncInertia()`
   (TransientIntegrator.cpp:166–171 → FE_Element.cpp:519–539). The element adds
   M·ü + C·v itself, exactly as every upstream UP element does. Omit this and
   Newton converges cleanly onto the wrong equation (p-row = steady seepage only).
2. **Rayleigh must be formed from SOLID-only blocks.** The base-class helpers
   (`Element::getResistingForceIncInertia`, `getRayleighDampingForces`, the base
   `getDamp`, and the `setRayleighDampingFactors` Kc snapshot) all compose
   βK·`getTangentStiff()` — which under this contract contains −Q and H, injecting
   spurious u-forces from ṗ and spurious ṗ-diffusion damping. LadrunoUP overrides
   `getDamp`, `getResistingForceIncInertia`, and `getRayleighDampingForces`;
   C_Ray = αM·M_s + βK·K_s(current) + βK0·K_s(initial) + βKc·K_s(committed),
   u-u block only, with the element keeping its own committed solid-K copy
   (porting reference: SSPquadUP.cpp:437–449 — noting it mishandles βK0/βKc,
   which we fix). Contract test: with βK-only damping, every p-row/col of
   `getDamp()` equals exactly the Qᵀ/S/H̃ blocks.
3. **`setParameter` is the updateMaterialStage transport** — the stage command
   reaches materials *through the element* (MaterialStageParameter.cpp:60–86
   calls `theEle->setParameter({"updateMaterialStage",…})`). LadrunoUP forwards
   the catch-all to its materials (upstream idiom) **and co-registers itself** so
   `updateParameter` can dirty the stabilization-α / skeleton-moduli / K₀ caches
   (§3.3); the same path carries xPerm/yPerm/zPerm.

`update()` maps trial nodal u to strains and calls `setTrialStrain` only — **p
never enters the constitutive update**; it acts solely through Q in assembly (the
effective-stress seam, §1.4).

**This is not an invention — it is the textbook's own arrangement.** The book's
recommended discretization is literally this block structure (Book eq 3.66, p.68:
mass [M,0;0,0], "damping" [0,0;Q̃ᵀ,S̃], "stiffness" [K,−Q̃;0,H]), and ZS84 (p.81,
re eq 34) states plainly that the system "is no longer symmetric and indeed shows
the typical form encountered in many coupled problems." Upstream OpenSees's
rate-DOF trick is the deviation; we are returning to the canonical form.

What this buys: p is directly the recorded / constrained quantity (prescribed-head
`sp` constraints under load patterns; classic `recorder Node -dof … disp` works
today, per-node guarded); **static analysis becomes well-posed** as drained solid +
steady-state seepage (H·p = f_seep), a genuine capability upstream cannot offer;
no unbounded DOFs on long runs; semantics match Abaqus/Kratos (their p DOF is p).
Two honesty riders on the headline ⟨UP-6/FW-F7⟩:

- **Static well-posedness has two necessary conditions**: k̄ > 0 and **≥1 drained
  (p-fixed) node per hydraulically connected region** — otherwise H is a
  pure-Neumann Laplacian (H·1 = 0) and the block-triangular static tangent
  det(K)·det(H) is structurally singular (an all-impervious sealed model is the
  first thing a user will build). Parser/guide warn; P1 adds an
  all-impervious-static smoke expecting the loud singularity. Conversely,
  statically p is *not* a Lagrange multiplier (H alone determines it in the
  triangular system), so equal-order pairs **cannot checkerboard on the static
  path** — H̃ living in the never-assembled damp block is immaterial there.
- **Initialization is a recipe, not a magic IC**: (a) a static steady-seepage
  stage (the genuinely new route), (b) a gravity-stage transient, or (c) scripted
  `setNodeDisp $node $pDof $val -commit` (no bulk `ic` command exists; it triggers
  no element update). Displacement-zeroing staged tricks (`InitialStateAnalysis
  off` → `revertToStart`) **also zero nodal p** under this contract — sequence
  them *before* establishing the p field (guide + P4 recipe).

What it costs: the effective tangent c₁K+c₂C+c₃M is **unsymmetric** (−c₁Q vs c₂Qᵀ
land in different slots — unavoidable once the two Q blocks live in different
matrices) ⇒ requires a general solver, and **the wrong-solver failure mode is
silent, not loud** ⟨FW-F1⟩: the no-`system`-command default is **ProfileSPD**,
whose assembly keeps only upper-triangle-in-profile entries — one Q block is
*silently discarded* and the solve returns plausible wrong pore pressures.
Mitigations: mandatory parser notice at element creation; the P1 "wrong-solver
smoke" is a *documented divergence comparison* (ProfileSPD vs UmfPack on the
Terzaghi column) + a LEDGER_quirks row — no framework hook lets an element reject
an SOE. Available general solvers: UmfPack / SuperLU / FullGeneral / BandGeneral
(serial, unconditional); **MUMPS SYM=0 in the MPI targets only** (serial `system
Mumps` is not compiled in this build). Upstream comparisons become
**response-level** equivalence gates rather than matrix-level reduce-to (§7).

> [!note] Future solver optimization (not v1)
> The book shows the coupled Newton Jacobian "can be made symmetric by a simple
> scalar multiplication of the second row (provided K_T is itself symmetric)"
> (Book eq 3.47, p.62) — the scaling depends on the integrator constants, so it
> belongs to an integrator/SOE-level hook, not the element. If unsymmetric-solve
> cost ever matters at scale, a LadrunoUP-aware row-scaling wrapper can restore
> symmetric solvers without touching the element contract.

Rejected: shipping both contracts behind a flag — dual DOF semantics doubles the
test matrix and poisons the recorder/BC story. Transient discretization of ṗ rides
the structural integrator: Newmark kinematics give a consistent first-order rule
that carries a parasitic history root −(1−γ)/γ on the p-row, **damped only for
γ > ½** — one more reason the γ ≥ 0.6 rider (§3.6) is load-bearing ⟨UP-4⟩. The
mass matrix stays singular on p-rows in *any* contract (→ explicit constraint, §6).

### 3.3 Stability: the inf-sup obligation and the two cures

In the undrained/incompressible limit (k→0, K_f/n ≫ K_skeleton) p becomes a
near-Lagrange-multiplier for Qᵀu̇ ≈ 0 and the discrete pair (Nᵤ, N_p) must satisfy
LBB/inf-sup, or pressure **checkerboarding** and volumetric locking appear. Two
standard cures, both adopted, chosen per shape:

1. **Taylor–Hood (mixed order)** — quadratic u, linear p: LBB-stable, no tuning.
   Upstream's 9-4 and 20-8 elements; Abaqus's C3D20P (quadratic bricks carry p at
   corners only). Our quadratic **Bézier shapes get this natively** — with the
   pressure basis pinned precisely ⟨UP-8⟩: **p uses the linear barycentric basis
   L_i attached to the vertex nodes** — NOT the vertex *subset* of the Bernstein
   basis ({L_i²} is not a partition of unity, cannot represent constant p, and
   gives H·1 ≠ 0). Bernstein-u spans full P2 on the same parametric simplex, and
   under the family's **straight-side guard the shared geometry map is affine**
   (∇L_i exact constants) — that guard is therefore a *TH-geometry requirement*
   here, not just a u-BC nicety. Vertex interpolatory-ness of Bernstein matters
   only for u Dirichlet BCs and for the p DOF sitting at a material point.
   Corner nodes run ndf = ndm+1, mid-edge nodes ndf = ndm (heterogeneous
   ndf — the Nine_Four_Node_QuadUP precedent).
2. **Equal-order + stabilization** — for T3/Q4/H8 where all nodes carry p:
   augment the storage block with the pressure-Laplacian stabilization
   **S\* = S + H̃**, H̃ = ∫(∇N_p)ᵀ α ∇N_p dV, acting on ṗ (McGann 2012 eqs 70+74;
   2015 eqs 70+73; mirrored in-tree at SSPquadUP.cpp:490–544 /
   SSPbrickUP.cpp:550–609; in our contract it lands in the damp p-p slot).

   **Provenance (paper-verified):** the *direct α-method* (Zienkiewicz–Huang–
   Pastor 1994; Huang–Yue–Tham–Zienkiewicz 2004) — add α·∇·(rate of mixture
   momentum) to mass balance, yielding the Laplacian plus a stress-rate term the
   authors **deliberately omit** (would unsymmetrize + need stress recovery),
   making the scheme non-residual, "similar to Brezzi–Pitkäranta" (2012 p.303).
   NOT polynomial pressure projection — Dohrmann–Bochev stays the documented
   fallback operator, and the book's fractional-step split (§5.5, eqs 5.58–5.78)
   is the third route (adds the same Laplacian via the operator split; conditional
   stability Δt ≤ h/c_s — noted for the P7 explicit runway, not v1).

   **α pinned from the papers** (2012 eq 73, 2015 eq 74, identical):
   **α = α₀·h²/(K_s + (4/3)G_s), α₀ ∈ [0.1, 0.5]** — K_s, G_s are the **drained
   skeleton** moduli (both papers verbatim; 2015 App. B derives them from the
   effective-stress material's *initial isotropic elastic tangent* — exactly how
   our `auto` computes them per element, per material). Neither paper defines h
   for irregular elements; back-calculating their own examples supports **h =
   largest element dimension** (2012 footing → α₀≈0.254 with h=3 m; site-response
   closes only with h = the larger dimension) — we pin that definition and
   document it. Surface: `-stab auto <$alpha0>` (default α₀ = 0.25, warn outside
   [0.1, 0.5]) | `off` | `$alpha` (raw, upstream-style). Their only guidance:
   "limited sensitivity studies suggest the range … is acceptable"; h² scaling
   makes H̃ vanish under refinement — no over-stabilization diagnostics exist, so
   the P1 gate includes an α₀-sweep on the checkerboard test.

   **Independence confirmed:** the stabilization is NOT tied to their single-point
   solid — their own *fully-integrated* Q1-P1/H1-P1 comparators checkerboard too
   (2012 Fig 8a; 2015 Fig 4c): the instability is inf-sup-driven, and the
   non-residual Laplacian rides any solid treatment. Our full-integration Q and S
   are strictly more consistent than their rank-deficient one-point versions
   (2012 eqs 64–67). Porting rule inherited: **full integration of H and H̃ in
   3D** — their single-point H "does not produce an acceptable element" in 3D
   (2015 §4, p.130); we integrate both exactly everywhere.

   Two upstream stabilization traps deliberately dodged: our H̃ needs skeleton
   moduli only at *setup* (recompute on `updateMaterialStage` — SSP freezes its
   assumed-strain Kstab at the initial elastic tangent, wrong after staging/
   plasticity), and our S/H̃ never silently vanish when material density is zero
   (SSPquadUP returns a zero p-p mass block if ρ==0).

   Book cross-check on pairs (Fig 3.2, pp.63–64): BB-safe = Q8P4, T6P3,
   biquadratic/bilinear, MINI-type — our Taylor–Hood picks; equal-order Q4P4/T3P3
   "not fully acceptable at incompressible-undrained limits" though usable at
   high permeability — exactly the regime split our `-p`/`-stab` axes encode.

### 3.4 Solid-part formulation axis

`-formulation std|bbar`. B-bar (mean dilatation) matters twice here: (i) dilatant
sand plasticity near-undrained is volumetric-constraint-dominated (why upstream
ships BBar*UP twins); (ii) it composes with *both* pressure axes. Not carried into
v1: EAS (interaction with the p field is genuinely non-trivial), URI/hourglass,
`-geom corot|finite` (finite-strain consolidation = reserved research phase; the
ADR-70 kernel is the natural partner when it opens). *Post-v1: `-geom corot`
shipped by ADR 78 (3D lanes, per-block frame decisions recorded there);
`finite` remains reserved.*

### 3.5 What we keep vs fix from upstream (code-review findings)

**Keep** (proven conventions):
- per-node DOF interleaving [u…, p] with p last (`fix … 1` on the p slot =
  drained boundary; with our §3.2 contract p records as a nodal *disp* component);
- SSP-lineage hygiene: per-instance matrices, material `setResponse` passthrough,
  **raw K_f + porosity** input (element computes n/K_f itself) rather than quadUP's
  pre-combined `kc` the user must derive;
- exact integration of H (SSPbrickUP does 2×2×2 even though its solid is
  single-point);
- the 9-4/20-8 heterogeneous-ndf DOF-mapping approach (our Taylor–Hood substrate) —
  but generated from one mapping utility, not hand-rolled per element;
- `-lumped` solid-mass switch, SelfWeight load-pattern hook (`loadFactor·b`),
  live `setParameter` permeabilities.

**Fix** (each item is a verified upstream defect that becomes a LadrunoUP test):
- the rate-DOF contract → honest p (§3.2);
- replication drift: the initial-displacement feature exists in **three mutations**
  (FourNodeQuadUP's is doubly buggy — writes only index [0] in setDomain:284–314
  and tests the wrong node's pointer at update:390; Nine_Four reimplements it
  correctly; BBar quad silently dropped it). We drop the feature entirely —
  staged/initial-state analysis is the fork's answer ([[project_initdefgrad_staged]]);
- **zero element-level pore-pressure/flux responses anywhere upstream** (only the
  nodal DOF-rate trick records p) → LadrunoUP ships GP `porePressure`, effective
  *and* total `stresses`, Darcy `flux` responses from day one;
- serialization holes: BrickUP/SSPbrickUP never send `massType` (and never init it
  in the null ctor — indeterminate after MP migration); SSPquadUP never sends its
  stabilization α → ours serializes *all* ctor state, MP round-trip gated in P1;
- parameter-map bugs: brickUP `zPerm` wired to the y-permeability id (case 103
  unreachable, acknowledged in-tree); 9_4's "pressure" parameter is dead
  (settable, never applied) → ours: xPerm/yPerm/zPerm wired + tested;
- response bugs: BrickUP `stress3D6`/`strain3D6` *assign* `σ·0.125` per GP instead
  of accumulating (returns last-GP÷8); FourNodeQuadUP's setResponse stamps
  `eleType="BrickUP"`;
- load-path bugs: SSPquadUP's side-traction block leaks load into a **pressure DOF**
  (index bug, .cpp:1318) — v1 therefore ships `-body`/SelfWeight only and defers
  surface tractions to standard patterns; SSPbrickUP's ctor sets `applyLoad=1`
  then zeroes `appliedB` (body forces silently dead until the first load cycle);
- five parsers advertise `type? dM? dK?` tokens no code reads; every UCSD element
  carries **two parallel parsers** (classic-Tcl + OPS_) → LadrunoUP: one OPS_
  routine, both dispatch tables point at it;
- `mixtureRho()` clones with a commented-out "real" mixture rule and hardcoded
  e=0.7 in six files → we pin the convention: solid/mixture ρ comes from the
  material (`getRho()` = saturated mixture density, documented), fluid ρ_f is the
  element arg used *only* in the seepage body-force term (upstream's actual
  behavior, now stated instead of implied);
- statics-as-scratch (`getDamp` clobbering the shared K buffer mid-computation,
  class-wide `theNodes`) → per-instance buffers (ADR-40 alignment).

### 3.6 Validity envelope (documented, not enforced)

u-p drops the fluid acceleration **relative to the solid** (both momentum
equations; solid ü retained — Book §2.2.2 p.24; ZS84 §7 p.79). The quantitative
envelope is Zienkiewicz–Chang–Bettess (Géotechnique 30(4):385–395, 1980, eq 26;
Book Fig 2.2 p.29 — paper in hand): for a layer of depth L under angular
frequency ω, with V_c² = (D + K_f/n)/ρ and β = ρ_f/ρ,

  Π₁ = k′·V_c² / (g·β·ω·L²)  (drainage/loading-rate ratio),  Π₂ = ω²·L² / V_c²  (rate/wave-speed ratio)

(k′ = hydraulic conductivity [m/s].) The chart is log-log (Π₁ ∈ [10⁻², 10²]
abscissa, Π₂ ∈ [10⁻³, 10²] ordinate); Π₁ < 10⁻² is the undrained strip,
Π₁ > 10² drained; two downward-sloping boundaries (drawn at ~3% solution
discrepancy — "imprecise" per the paper) split the band:

- **Zone I** (bottom, slow): consolidation equation suffices (ü *and* relative
  fluid inertia negligible) — **u-p valid here too**;
- **Zone II** (middle): **u-p valid**, consolidation invalid — the earthquake band;
- **Zone III** (top: high frequency and/or very permeable): full Biot only —
  **out of scope** (upstream's u-p-U bricks survive only as reserved class tags
  38/39; source removed from tree).

Concrete calibration from the paper's dam example (L=50 m, V_c=1000 m/s, T ∈
[0.05, 5] s): **u-p covers the complete earthquake frequency range for
k′ < 10⁻³ m/s**; full-Biot effects only intrude at k′ ≳ 10⁻³ m/s *and* very high
frequency (ZCB80 pp.394–395; Book p.30 — beware the book's misprint "less than
0.5 s" where the chart implies *more*). The guide redraws Fig 2.2 from the paper.

Two dynamic artifacts to inherit knowingly (ZS84 Example 1, pp.84–88):
(a) *physical* compressible-wave ringing of p (period ≈ the undrained column's
natural period), damped naturally; (b) *numerical* p-oscillation under
step/coarse-Δt trapezoidal integration. Remedy = GN22/GN11 algorithmic damping:
**γ=0.6, β=0.3025** (ZS84 p.87 verbatim; = Book's β₂=0.605/β₁=β̄₁=0.6 first set,
p.61) — with the book's two riders: Dewoolkar (1996) found the 0.6-set
*over-damps* vs centrifuge data and the milder **γ=0.51, β=0.2575** set "gave
very good comparisons"; and with realistic soil plasticity + seepage the physical
damping usually suffices (Book §3.2.6 p.69). Guide documents both sets; P4 gates
pin the step-load column under each.

---

## 4. Architecture

```
SRC/element/ladrunoUP/
  LadrunoUPKernel.h        # PURE (namespace ladruno_up, raw doubles, runtime nN/nP):
                           #   Q/H/S/stab block assembly, fluid body force,
                           #   physics-mode seam (mode 1 = Biot; mode 2 reserved)
  LadrunoUPShapes.h        # stateless shape providers: T3, Q4, H8, BezT6, BezTet10
                           #   (N, dN/dξ, GP rules, Jacobian; lifted from family code)
  LadrunoUP.{h,cpp}        # the element: DOF map (incl. heterogeneous ndf),
                           #   update/commit cycle, block→API placement (§3.2),
                           #   loads, responses, sendSelf
  OPS_LadrunoUP.cpp        # parser (fatal ndf gate, per-shape defaults)
  CMakeLists.txt
tests/
  ladruno_up_kernel_check.cpp   # standalone g++ vs numpy oracle (ADR-70 idiom)
  ladruno_up_reference.py       # the oracle: blocks, stabilization, 1-el consolidation ODE
```

### 4.1 Element command (draft surface — pinned at P1)

```tcl
element LadrunoUP $tag $n1 … $nk $matTag \
    <-thick $t>                 ;# 2D only (plane strain implied), default 1.0
    -Kf $Kf -poro $n -rhoF $rhof \
    -perm $k1 $k2 <$k3>         ;# per-axis, k_hydraulic/γ_w  [L³·T/M]
                                ;#   (upstream's verified convention: H balances
                                ;#    ρ_f·perm·b with b an acceleration)
    <-permH $k1 $k2 <$k3> -gammaW $gw>  ;# sugar: hydraulic conductivity ÷ γw internally
    <-alpha $biotAlpha>         ;# default 1.0  <-Ks $Ks>
    <-body $b1 $b2 <$b3>>       ;# accelerations, default 0; solid rows use material ρ (sat.)
    <-fluidBody $f1 $f2 <$f3>>  ;# defaults to -body; drives the seepage source
    <-formulation std|bbar> <-pOrder equal|linear> <-lumped> \
    <-stab auto <$alpha0> | off | $alpha>   ;# auto: α=α₀h²/(K_s+4G_s/3),
                                            ;#   α₀ default 0.25 (papers: 0.1–0.5),
                                            ;#   h = largest element dimension,
                                            ;#   K_s,G_s from material initial tangent
    <-dynSeepage on|off>        ;# default on; off = upstream-parity (§3.1)
    <-geom linear>              ;# only accepted value; reserves the axis

# fluid bulk entered RAW (-Kf) with -poro n; storage uses 1/Q̄ = n/K_f + (α−n)/K_s.
# (quadUP's `bulk` is the pre-combined Q̄ ≈ K_f/n — the gate's arg-mapping converts;
#  SSP takes raw fBulk + void ratio e → n = e/(1+e). We take K_f + n directly.)
```

Naming ⟨FW-F8⟩: `-pOrder` (not `-p`) and `-Kf` (not `-bulk`) avoid semantic
collisions with the family's `-pressure` (surface load) and `-bv`
(bulk viscosity). **Unknown flags are parser-FATAL** — a deliberate break from
the family's warn-and-continue idiom, because a mistyped u-p flag silently
changes physics.

`(ndm, k)` selects the provider: (2,3) T3 · (2,4) Q4 · (2,6) Bézier T6 ·
(3,8) H8 · (3,10) Bézier Tet10.

**v1 legality matrix** ⟨scope-F4⟩: quadratic shapes accept only `-pOrder linear`
(TH — equal-order quadratic has no theory or gate behind it: parser-fatal, message
names the reserved future axis). Linear shapes accept `equal`; `linear` there is
accepted as a documented synonym (all nodes are vertices — identical
interpolation). On equal-order pairs, omitted `-stab` ⇒ `auto` (α₀ = 0.25) with a
one-line notice, so B4's `-stab off` leg is always an explicit opt-out; `-stab`
on TH pairs is parser-fatal (they don't need it).

**Loads & defaults** ⟨scope-F6⟩: ctor `-body`/`-fluidBody` are **always-on**
(upstream behavior); a `LOAD_TAG_SelfWeight` pattern **replaces** the active
values with loadFactor-scaled ones for BOTH the solid rows and the seepage source,
and `zeroLoad()` restores the ctor values — this is the staging knob the P4
gravity recipe drives. All other elemental loads are rejected with a named error
(surface tractions come later via standard patterns — not the SSP side-pressure
route, whose index bug §3.5 documents).

### 4.2 The growth axis (why this stays cheap)

Adding a geometry = **one provider + tests**: Q9 (equivalence gate:
`Nine_Four_Node_QuadUP`; honesty clause — the 9-node Lagrange basis is *new code*
with no in-fork donor), H20 serendipity (equivalence: `Twenty_Eight_Node_BrickUP`;
this is Abaqus C3D20P territory — and the in-flight **ADR-72 `LadrunoBrick20`**
becomes the natural shape-code donor once its H20 basis lands; the corner-p DOF
map is already first-class via the provider's p-carrier field, §2.2), LadrunoLST
T6-Lagrange once ADR-70 P3 lands (`-shape lst`, §2.1), wedges. No new class, no
new tag, no new broker/dispatch rows; the existing LadrunoUP ledger/banner rows
are *amended*, not duplicated (one row per feature). The kernel is runtime-sized;
nothing structural changes at 68 DOF (20-8).

### 4.3 Forward-compatibility: hybrid u/p (incompressibility) as mode 2

Same two-field skeleton, same Q, same pressure providers, same stabilization
machinery; differs in the p-row law (algebraic compliance −(1/K)M_p instead of
S·ṗ + H·p — no rates, so a genuine **static** path) and in the constitutive
contract (element takes the deviatoric part of the material response; p̂ carries
mean stress — dSNPO §15.3). Kernel requirement now: block assembly parameterized by
a mode struct so mode 2 adds a p-row law, not a rewrite. The condensed-pressure
(Abaqus C3D8H-style, element-internal p̂) variant would reuse the EAS
static-condensation wiring precedent. Reserved phase, own gate (§7).

### 4.4 Modeling surface & ecosystem

- **Mixed-ndf regions**: dry (ndf=ndm) and saturated (ndf=ndm+1) zones coexist via
  separate node sets + `equalDOF` on shared u-DOFs at the interface (the ADR-59
  documented pattern; user-guide material, not element magic). Guide caveat
  ⟨FW-F10⟩: use the **explicit DOF-list form** of `equalDOF` — the bare form
  sizes its constraint from the model-builder NDF and mis-sizes against the
  smaller-ndf node. apeGmsh's per-node ndf support covers emission; the apeGmsh
  emitter + `ladruno_apegmsh_contract` row are a companion (apeGmsh-repo) runway
  item, not this ADR.
- **BCs/ICs**: p-fixity = drained boundary; no fixity = impervious (with the §3.2
  static rider: statics need ≥1 drained node per connected region). Initial p per
  the §3.2 recipe (steady-seepage stage / gravity transient / `setNodeDisp
  -commit`), sequenced AFTER any displacement-zeroing staged step
  ([[project_initdefgrad_staged]] interplay documented in the guide).
- **Recorders** — one confirmed work item ⟨FW-F4⟩: the fork's own
  `Ladruno_NodeResults` PRESSURE source **hardcodes the upstream vel-trick**
  (reads `getTrialVel()[ndm]` — under our contract that is ṗ, silently wrong),
  and its DISPLACEMENT source buffers only the first ndm components (p does NOT
  "ride the disp channel" there). P1 makes `PressureSource` contract-aware
  (disp-slot for LadrunoUP nodes, vel-slot for upstream UP elements) and P4 gates
  it. The frozen MPCORecorder PRESSURE result keeps the vel-assumption and stays
  upstream-only (change-budget: our-hooks-only). Classic `recorder Node -dof …
  disp` works unmodified. Element `setResponse` adds `stresses` (effective),
  `stressesTotal`, `porePressure` (GP), Darcy `flux`; LadrunoRecorder/Monitor get
  the new topology rows per [[06_quadrature_global_gp_plan]]; charLength
  handshake mirrors the family.
- **Registration checklist** (all three dispatch points + plumbing):
  `OpenSeesElementCommands.cpp`, `TclElementCommands.cpp`,
  `FEM_ObjectBrokerAllClasses.cpp`, `ladrunoUP/CMakeLists.txt`, LADRUNO header
  stamp (`stamp_headers.py` GLOBS), `LEDGER_implementations` row,
  `banner_features.txt` line (at ship).

---

## 5. Class tags & registration

- `ELE_TAG_LadrunoUP = 33017` (ELE registry; free-slot comment at
  **classTags.h:942** — the line-943 note carries a stale "33016–33019 free"
  that gets corrected when P1 touches the file; 33016 stays reserved for
  LadrunoLST per ADR-70). Cross-registry occupants to stamp in the define
  comment, per house precedent ⟨FW-F9⟩: ND 33017 = LadrunoConcrete3D,
  ND 33018 = LadrunoRCFiniteStrain, LADRUNO 33019 = ComplexEigen — numerically
  equal, different registries, not collisions.
- **Reservation mechanics** ⟨scope-F7⟩ — the fork already lost one tag race
  (StagedStrain took 33014): **the ADR PR itself lands the
  `LEDGER_implementations` row** `LadrunoUP — ELE 33017 — RESERVED (ADR-71), not
  yet built` (docs-only reservation, the LogStrain2D #203 precedent); the
  classTags.h `#define` + comment land at P1. The ADR-72/33018 pencil is
  **out-of-tree hearsay until ADR-72 lands its own reservation row** —
  first-recorded-in-tree wins.
- **Mode 2 (P6 hybrid) ships under the same ELE tag 33017 behind the mode flag**
  — consistent with §4.2's no-new-tag promise and §8's naming rationale; a
  separate class, if the P6 mini-ADR ever wants one, is a new reservation.
- Kernel + providers: **no classTag** (pure helpers, the `*Kernel.h` idiom).
- Vanilla touches (strictly additive, `// Ladruno`): `SRC/classTags.h`, broker,
  the two command registries, top-level element `CMakeLists` include —
  vanilla-ledger rows in the shipping PRs per house rules.

---

## 6. Explicit dynamics: the honest note

The pressure equation is **first-order (parabolic)** — central difference does not
apply to the p-row, and the p-rows carry zero mass (§3.2), so the existing explicit
stack (CentralDifferenceLadruno/SMS — global lumped M) cannot advance them.
Structural, not a bug. Near-term explicit answer for undrained problems stays
**material-level** (`FluidSolidPorousMaterial`, confirmed in-tree — no extra DOFs,
works today). True explicit u-p is a research phase with three established routes:

1. **Staggered partitioning** — Zienkiewicz, Paul & Chan, IJNME 26(5):1039–1055
   (1988): unconditionally stable implicit-p/explicit-u operator splitting;
2. **Fractional-step** — Huang & Zienkiewicz, IJNME 43(6):1029–1052 (1998) (+
   Huang, Wu & Zienkiewicz, SDEE 21(2):169–179, 2001): Chorin-type pressure
   correction that *also* buys equal-order stability;
3. **Fully explicit** — Xu, Feng, Song, Du & Zhao, SDEE 141:106452 (2021):
   diagonalize both M and S, central-difference u + explicit p update;
   second-order, conditionally stable.

Common requirements for any route: **lumped S ≠ 0** (finite K_f/n — the exact
incompressible limit is off the explicit table) and a CFL governed by the
**undrained** wave speed plus a diffusion limit as k̄ grows — both interacting with
the ADR-65 Δt machinery ([[65_ladruno_explicit_dt_strategies_adr]]). The book's own
verdict is sobering and worth quoting in the eventual ADR: explicit u-p (β₂=0,
diagonal M) exists but "its limitation is very serious … invariably the
unconditionally stable, implicit form is used" (Book §3.2.4); and staggered schemes
carry the Turska–Schrefler (1993) **lower bound on Δt/h²** — you cannot refine Δt
without refining h. Own ADR once P0–P4 are real; until then the parser documents
the constraint and the integrator scorecard
([[49a_integrator_scorecard_2026-06-23]]) gains a u-p row.

---

## 7. Phases & exit gates

| Phase | Scope | Gate |
|---|---|---|
| **P0** | `LadrunoUPKernel.h` + providers + numpy oracle (no OpenSees build) | oracle reproduces Q/H/S/L blocks vs numpy integration per shape; block pattern/symmetry checks; 1-element consolidation ODE (kernel blocks → closed-form decay) matches analytic; **inf-sup smoke, pinned** ⟨scope-F9⟩: 4×4 structured patch per pair, u fixed on the whole boundary, p free; eigen-decompose the pressure Schur complement S_p = Qᵀ·K⁻¹·Q (+ H̃ when stabilized); gate = count of eigenvalues < 1e-10·λ_max — equal-order-no-stab ≥ 2 (constant + checkerboard), TH and stabilized exactly 1 (the constant mode); optional Chapelle–Bathe inf-sup constant over N ∈ {2,4,8} (non-degrading for TH); FD check quantifying the dropped dynamic-seepage tangent (§3.1); g++ check ≤1e-9 vs oracle |
| **P1** | `LadrunoUP` Q4 lane end-to-end (equal-order, `std|bbar`, `-stab`) | **upstream-equivalence vs FourNodeQuadUP, two legs** ⟨scope-F1/UP-4⟩ — run with `-dynSeepage off` (upstream lacks the term, §3.1): **(leg 1, tight)** Newmark **γ=½, β=¼** — the unique pair where the p-as-disp and p-as-vel parameterizations produce the identical discrete one-step rule (θ=½, zero memory term): u and p histories ≤1e-6 rel, convergence tests unbalance-based (upstream's ∫p·dt DOF grows unboundedly; disp-increment tests mis-scale); **(leg 2, production)** γ=0.6/β=0.3025: Δt-halving study, gate on mutual convergence (observed order ≥ 1, difference ≤ 1% at finest Δt) — 1e-6 is reserved for identical-operator checks only; SSPquadUP pressure-block cross-check at matched α; **Terzaghi 1D consolidation vs series** + **B2 ZS84 column** (tol ~1e-3 rel vs analytic at pinned Δt); **B3 Boone–Ingraffea** (hard numbers, §7.1); **B4 checkerboard** via the CB index (§7.1) with α₀-sweep; **NEW-capability gate: static drained + steady-seepage** solves + **all-impervious-static smoke** (expects loud singularity, §3.2) + **wrong-solver divergence comparison** (ProfileSPD vs UmfPack, documented + LEDGER_quirks row); Rayleigh contract test (βK-only → p-rows/cols of getDamp = exactly Qᵀ/S/H̃, §3.2); `bbar` gated on *behavior* (near-undrained no-lock vs `std`), NOT on BBarFourNodeQuadUP equality (theirs B-bars a 4-component material via the parameter-20 hack); FD consistency of all blocks incl. getResistingForceIncInertia; PressureSource recorder made contract-aware (§4.4); sendSelf round-trip incl. every ctor arg (`-lumped`, α₀, dynSeepage…); MP smoke |
| **P2** | H8 + T3 lanes | upstream-equivalence vs **BrickUP** (two-leg methodology as P1, `-lumped` both ways); 3D Terzaghi; **B4-3D** (McGann 2015 cube-footing checkerboard analog, CB index); **B4-T3** ⟨scope-F12⟩: B4 re-meshed with crossed-diagonal T3 split, CB gate + α₀-sweep — pins that the §3.3 α transfer to simplices holds (h = largest element dimension); T3 equal-order documented as the honest baseline (its locking pathology inherited from CST is *worse* under undrained constraint — measured and pinned) |
| **P3** | Taylor–Hood on Bézier T6 + Tet10 (vertex-p, heterogeneous ndf) | no-stabilization stability demonstrated: **B1 ZCB80 periodic column** at the undrained limit (Q̄ = 10⁹ MPa case — the pressure-locking benchmark, §7.1) + mesh-refinement pressure L2 convergence, no checkerboard; quadratic-u rate preserved; heterogeneous-ndf plumbing (DOF map, sendSelf, recorders, numberers) proven; straight-sided guard carried; equalDOF mixed-ndf interface example; 9_4-style cross-check vs upstream on a shared quad mesh where geometry permits |
| **P4** | Dynamics validation + ecosystem | **B5 Simon–Zienkiewicz–Paul full-Biot closed form** (§7.1): hard gates on total-stress front position/amplitude and π̂ plateau **at stations ξ > 40**; û(0,τ) is *measured-first* ⟨UP-5⟩ — a 1D u-p semi-discrete oracle pre-run quantifies its full-Biot-vs-u-p discrepancy (17% of û(0,50) rides on drained-face outflow ŵ inside the deviation layer); if > 0.5% it demotes to a documented-comparison plot; step-load p-oscillation pinned under BOTH Newmark sets (γ=0.6/β=0.3025 and γ=0.51/β=0.2575, §3.6); free-field liquefaction column vs upstream quadUP reference (PDMY + `updateMaterialStage` staging — exercises the §3.2 setParameter forwarding + H̃/K₀ cache-dirty path; two-leg equivalence methodology); gravity/hydrostatic init recipe in the guide (sequenced per §3.2); PressureSource/Monitor channels gated; banner + ledgers + user guide (`LadrunoUP_guide.md`) |
| **P5** *(demand)* | growth-axis proof: Q9 and/or H20 provider | upstream-equivalence vs Nine_Four_Node_QuadUP / Twenty_Eight_Node_BrickUP |
| **P6** *(reserved)* | hybrid u/p mode 2 (incompressibility, static path, deviatoric contract) | own mini-ADR amendment + gates (ν→0.5 cantilever/Cook's membrane family: std locks, mode-2 doesn't; matches ADR-70 F-bar cross-checks) |
| **P7** *(research)* | explicit u-p strategy | own ADR (see §6) |

Adversarial-gate policy (per [[feedback_adversarial_gate_when]]): **full Opus gate
at P0** (novel coupled math: block placement, stabilization, TH pairing) and **P3**
(heterogeneous-ndf novelty on Bézier). P1/P2 are carried by the upstream-equivalence
+ analytic gates; P4 is validation, not new math.

Each phase is one PR off `ladruno`, ledgers updated in-PR.

### 7.1 Pinned benchmark dossier (paper-verified parameters)

**B1 — ZCB80 periodic-load layer** (analytic, frequency domain; Book Fig 5.25
instance): column L=30 m × 1 m, q = 100·e^{iωt}, ω = 3.379 rad/s; E = 7.492×10⁸ Pa,
ν = 0.2, n = 0.333, k′ = 10⁻⁷ m/s, ρ_s = 2000, ρ_w = 1000 kg/m³; drained top,
impermeable rigid base; run Q̄ = 10⁴ MPa (compressible) and 10⁹ MPa
(undrained-limit → the checkerboard stressor). Reference: ZCB80 closed form
(eqs 12–30). **FE realization** ⟨scope-F5⟩ (OpenSees is real-valued/time-domain):
drive q(t) = 100·sin(ωt) transient; ring up N cycles until per-cycle amplitude
drift < 0.5%; least-squares fit A·sin(ωt) + B·cos(ωt) per node over the last full
cycle to reconstruct the phasor; compare |p̂(z)|, arg p̂(z), |û(z)| vs the closed
form — L2-rel tolerance pinned at implementation (target ≤ 2% for the Q̄=10⁴
leg); the 10⁹ leg additionally gates CB ≈ 0 and mesh-refinement convergence.
Δt ≤ T/100, pinned by one Δt-halving check (γ=0.6 numerical damping biases
steady-state amplitude at coarse Δt).

**B2 — ZS84 Example-1 column** (consolidation limit + oscillation behavior):
30 m column; load 1.0 kN/m² ramped by sin(πt/2t_p), t_p = 0.1 s; E = 30 MN/m²,
ν = 0.2, ρ_s = 2, ρ_f = 1 Mg/m³, K_f = 100 MN/m², K_s = ∞, k′ = 10⁻² m/s, n = 0.3;
top drained, base impermeable/fixed. Gates the 1D consolidation solution AND the
γ = 0.5 vs 0.6 oscillation comparison (ZS84 Figs 2/5/7).

**B3 — Boone–Ingraffea poroelastic column** (quasi-static, closed-form hard
numbers; Book §6.4.1 eqs 6.7–6.9): 6 m × 1 m, side rollers; instantaneous
q = 1 MPa and p = 1 MPa at top; E = 15960, G = 6000, K_s = 36000, K_f = 3000 MPa;
ν = 0.2, ν_u = 0.33, B = 0.62, n = 0.19, α = 0.79, k = 2×10⁻⁵ m²/(MPa·s).
**Targets: u(0⁺) = 0.254 mm, u(∞) = 0.079 mm, p(0⁺) = 0.410 MPa.** Also the only
gate exercising α ≠ 1.

**B4 — McGann checkerboard footing** (stabilization gate): 2D — 30×30 m,
uniform 10×10 mesh (h = 3 m), half-model symmetry; 3 m strip load ramped
0→0.1 kPa over 0.1 s then constant; drainage top surface ONLY; E = 25000 kPa,
ν = 0.3, ρ_sat = 2.67 Mg/m³, ρ_f = 1.0, k′ = 10⁻⁷ m/s, **K_f = 2.2×10¹² kPa**
(engineered incompressible-impermeable limit); Rayleigh 0.05M + 0.02K (solid
rows). Check field at t = 1 s: `-stab off` checkerboards, α = 6.8×10⁻⁵
(α₀ ≈ 0.25, h = 3 m) clean (2012 Figs 7–8). 3D analog: 30 m cube, quarter
symmetry, 3×3 m patch, graded 12³ mesh, e = 0.7, Rayleigh 0.05M + 0.0003K
(2015 Figs 2–5). **Quantitative CI metric** ⟨scope-F3⟩: checkerboard index
CB = |⟨p, χ⟩| / (‖p‖·‖χ‖) with χ the alternating lattice mode over the
structured mesh's interior nodes (χ_ij = (−1)^{i+j}; (−1)^{i+j+k} in 3D — the
graded 3D mesh is still structured-indexed). Primary gate: CB(off)/CB(auto) ≥ R
with R and the absolute CB(auto) ceiling pinned at P1 from the first passing run
and recorded in the test header; the α₀-sweep gates CB(auto; α₀)
monotone-bounded across [0.1, 0.5].

**B5 — Simon–Zienkiewicz–Paul dynamic column** (exact full-Biot closed form —
numpy oracle, no numerical inversion): semi-infinite confined column, drained
top, σ(0,t) = F(t) ∈ {step, sine ω = 62.83, triangular spike}; materials must be
**constructed dynamically compatible** (ακ = β ⇒ K_f derived): base
E = 3000, ν = 0.2, ρ = 0.3060, ρ_f = 0.2977, n = 0.333, k = 0.004883 (their
k = k_hyd/(ρ_f g)); mats 1/2/3 = α ∈ {1.0, 0.667, 0.333} per Table II (Q =
1.201×10⁵ / 1.385×10⁴ / 1.441×10⁴). Hard values: σ̂(ξ,τ) = 1(τ−ξ) exactly;
π̂ plateau = β = 0.9730 and σ̂′ → 1−β (mat 2); v₂ = 1/√a ∈ {0.1153, 0.5092, 1.0};
ŵ(0,50) = 10.7, ŵ(0,150) = 18.7, û(0,τ) = −τ − β·ŵ(0,τ). FE mapping: column
length > V_c·T (no base reflections), lateral rollers, impervious sides.
**Paper errata the oracle MUST keep** (verified algebraically during
extraction — do not "fix back"): (i) eq 43's front term is **+f(τ−ξ)** (printed
−); (ii) eq 44 must read **σ′ = σ − π** (printed +); (iii) Table I spike falling
branch is σ₀(2Δ−t)/Δ (printed (t−2Δ)/Δ); (iv) p.390/Fig-3 labels swap materials
1↔3's Q — **Table II is the correct one**. Sanity pins: σ(0,τ) = F, π(0,τ) = 0,
deep π̂ → β. Expected u-p deviation: the second (slow) P-wave boundary layer
ξ ≲ 40 (u-p drops it by construction) — hard gates live at stations ξ > 40
(σ̂ front, π̂ plateau); û(0,τ) is measured-first per the P4 row ⟨UP-5⟩.
**FE discretization pinned at implementation and recorded in the test header**
⟨scope-F8⟩: suggest uniform mesh with ≥ 10 elements per v₂·τ_gate, Δt at CFL/2
of the fast wave; front position = 50%-of-plateau crossing, gated within ±2Δz;
π̂ plateau within ±2% over τ ∈ [50, 150]; one Δz,Δt-halving leg confirms
convergence toward the closed form.

*(de Boer–Ehlers–Liu 1993 is NOT reproduced in the book and stays on the optional
wishlist — B1/B2/B3/B5 already cover analytic static, consolidating, and dynamic
gates, so it is no longer blocking.)*

---

## 8. Risks / open questions

> [!note] **Unit/argument conventions — RESOLVED** (both reports in)
> Permeability input is **k̄ = k_hydraulic/γ_w** [L³·T/M] — the upstream docs are
> dangerously vague (they print a cm/s soil table with no conversion note!) but the
> code is unambiguous (H multiplies ∇p directly; seepage term balances ρ_f·b). We
> pin k̄ as the `-perm` meaning, and add parser sugar `-permH $k1… -gammaW $gw`
> that divides internally — killing the single most common quadUP user error.
> Bulk entry: raw K_f + porosity n (§4.1), converting quadUP's pre-combined kc in
> the gate mappings. Recorded p is compression-positive (σ = σ′ − αp·I,
> tension-positive σ). Also inherited-docs quirk worth a LEDGER_quirks row at P1:
> the new opensees.github.io SSPquadUP page describes α as "1 − K_s/K_f" — garbled;
> wiki + paper + code (stabilization parameter, §3.3) are authoritative.

> [!question]
> **Where does plane-stress go?** Nowhere — u-p is volumetric physics;
> `LadrunoUP` 2D is plane-strain only (parser-fatal otherwise). Axisymmetric is
> the real 2D sibling (dams, wells) — future provider axis + measure change,
> deliberately not v1.

> [!question]
> **Static-scratch break.** New element uses per-instance buffers (ADR-40
> alignment) while the rest of the family uses class statics. Accepted
> inconsistency — flagged so nobody "fixes" it backwards.

- **Blast radius**: zero on shipped elements (nothing existing is edited beyond
  additive registration rows) — that was the point of option (d).
- **Heterogeneous-ndf ecosystem friction** (P3): recorders/constraint handlers/
  numberers meeting ndf=3-and-2 mixes; Nine_Four proves the Domain handles it, but
  our Monitor/LadrunoRecorder paths have never seen it — P3 gate covers.
- **Materials**: RESOLVED — `SRC/material/nD/soil/` ships PressureDependMultiYield
  {,02,03}, PressureIndependMultiYield, FluidSolidPorousMaterial +
  `TclUpdateMaterialStageCommand`; P4 uses PDMY directly. One inherited contract to
  *avoid*: the BBar quad's `updateParameter(20)` "give me 4-component 2D stress"
  hack only works with UCSD materials and silently breaks others — our bbar stays
  plane-strain 3-component (family-consistent), so no material special-casing.
- **Solver requirement (not just preference)**: symmetric-storage solvers
  (ProfileSPD, SparseSYM, MUMPS SYM=2) are wrong, not merely slow — and
  **ProfileSPD is the no-`system`-command default, failing SILENTLY** (drops one
  Q block in assembly; plausible wrong pore pressures). Full story + mitigations
  in §3.2 ⟨FW-F1⟩: mandatory parser notice, P1 divergence-comparison smoke,
  LEDGER_quirks row. MUMPS SYM=0 is MPI-targets-only in this build. Assembly
  cost of the extra field is trivial vs solve.
- **PyLiq1/TzLiq1/QzLiq1 interop gap**: those materials read mean effective stress
  from solid elements via **friend access to FourNodeQuad/FourNodeQuadUP
  specifically** (FourNodeQuadUP.h:81–83) — they will not see LadrunoUP. Near-term:
  keep upstream quads under py/tz/qz interfaces; clean fix (a public
  `meanEffectiveStress` response instead of friendship) is a small follow-up if
  demand appears.
- **Naming**: `LadrunoUP` (command + class). Alternatives (`LadrunoPorous`,
  `LadrunoBiot`) rejected: the fork mirrors upstream's *UP suffix so geotech users
  find it, and mode 2 (hybrid) keeps "UP" honest ("u + p field", not "porous").

---

## 9. Alternatives considered (and why rejected)

| Option | Verdict | Killer fact |
|---|---|---|
| (a) *UP twin per element (upstream pattern) | ❌ | ≥5 new classes re-forking the hardest code (EAS inner Newton, finite paths); every future lane fix mirrored N ways — the documented upstream drift disease |
| (b) Wrapper/decorator owning an inner solid element | ❌ | inner elements' setDomain/parsers **reject pressure-carrying nodes** (ndf gates); no public GP/B API to build Q/H/S from — wrapper needs in-class edits anyway, forfeiting its selling point (same verdict pattern as ADR-70's SolidTransformation rejection) |
| (c) `-up` flag inside each element | ❌ | external DOF size changes per instance vs hardcoded sizes + static scratch + `2a+i` interleaving in ~6 methods × 5 classes; flag matrix multiplies against formulation/geom/hourglass axes on shipped, validated elements |
| (d) **one new unified element + providers + pure kernel** | ✅ | zero risk to shipped code; mechanics written once; geometry growth = provider-only; TH + equal-order + stab as orthogonal axes; ADR-70-style oracle testability |

---

## 10. References

**Formulation** *(papers IN HAND at `C:\Users\nmora\Desktop\Soil papers\` —
extracted 2026-07-10, all §3 claims verified against them)*
- Zienkiewicz OC, Shiomi T (1984). *Dynamic behaviour of saturated porous media;
  the generalized Biot formulation and its numerical solution.* IJNAMG 8:71–96.
  ✅ in hand (`zienkiewicz1984.pdf`).
- Chan, Pastor, Schrefler, Shiomi, Zienkiewicz. *Computational Geomechanics*,
  **2nd ed., Wiley 2022** ✅ in hand (`chan2022.pdf`, full book; 1st ed. 1999).
  Load-bearing anchors: eq 2.17 (storage), eq 2.21 (u-p), §3.2.3 GN22/GN11 +
  damping sets (p.61), eq 3.47 (row-scaling symmetrization), eq 3.66 (block
  arrangement = ours), Fig 2.2 (zones), Fig 3.2 (BB pairs), §5.5 (fractional
  step), §6.4 (Boone–Ingraffea, aquifer), p.56 (dynamic-seepage LHS drop,
  Chan 1988 insignificance).
- Zienkiewicz, Chang, Bettess (1980). *Drained, undrained, consolidating and
  dynamic behaviour assumptions in soils.* Géotechnique 30(4):385–395
  ✅ in hand (`zienkiewicz1980.pdf`): Π₁/Π₂ eq 26, zone chart Fig 3.
- Li X. *FE Formulation of Poro-Elasticity*, Stanford TR147 — open PDF: full u-p
  derivation used for the first-pass verification.
- Lotfian & Sivaselvan, arXiv:1506.06785; Zhao & Choo, arXiv:1905.00671 — open
  cross-checks (Biot system; equal-order checkerboarding survey).

**Stabilization / elements**
- McGann, Arduino, Mackenzie-Helnwein (2012). *Acta Geotechnica* 7(4):297–311
  ✅ in hand (`mcgann2012.pdf`; journal often miscited as SDEE) and (2015)
  *Computers & Geotechnics* 66:126–141 ✅ in hand (`mcgann2015.pdf`).
  α = α₀h²/(K_s+4G_s/3), α₀∈[0.1,0.5] (eqs 73/74); H̃ eqs 70+74/73; footing
  checkerboard §6.3 both papers.
- α-method lineage: Zienkiewicz, Huang, Pastor (1994); Huang, Yue, Tham,
  Zienkiewicz (2004); Brezzi & Pitkäranta (1984) (non-residual Stokes analogy);
  Truty (2001), Truty & Zimmermann (2006) (alternative α definitions).
- Bochev & Dohrmann (2004): polynomial pressure projection (documented fallback);
  de Pouplana & Oñate (2016): FIC-stabilized u-p in Kratos.
- Yang Z, Elgamal A — UCSD u-p manual (soilquake.net/opensees, OSManual 2008):
  the `vel`-records-p contract, quadUP arg conventions.
- Jeremić group, *Solution Verification … Fully Coupled Porous Media* — open PDF
  (CV-J33): u-p vs u-p-U, k̄ unit precedent.

**Benchmarks** *(pinned specs in §7.1)*
- B1 ZCB80 periodic layer (in-hand paper + Book Fig 5.25 instance); B2 ZS84
  Example 1; B3 Boone & Ingraffea 1990 poroelastic column (via Book §6.4.1);
  B4 McGann 2012/2015 footings; B5 Simon, Zienkiewicz & Paul (1984) IJNAMG
  8:381–398 ✅ in hand (`simon1984.pdf` — exact closed form; **four verified
  errata documented in §7.1**); Terzaghi series; VELACS (Arulanandan & Scott
  1993) for the liquefaction tier; Schrefler–Scotta 2001 column (Book §8.5,
  numerical-comparison only — optional).

**Explicit (P7 runway)**
- Zienkiewicz, Paul & Chan (1988) IJNME 26(5):1039–1055 (staggered);
  Huang & Zienkiewicz (1998) IJNME 43(6):1029–1052 + Huang, Wu & Zienkiewicz
  (2001) SDEE 21(2):169–179 (fractional step); Xu, Feng, Song, Du & Zhao (2021)
  SDEE 141:106452 (fully explicit, lumped M+S).

**Mode 2 (hybrid) theory**
- dSNPO (2008) §15.3 mixed u/p; Abaqus Analysis/Theory: C*P pore-fluid elements
  (corner-p quadratic precedent = our TH), C3D8H hybrids; Kratos GeoMechanics
  U-Pw family (one templated element over geometries + diff-order variants —
  independent convergence on option (d)'s shape).

**Fork**
- [[70_ladruno_plane_finite_triangles_adr]] (kernel idiom),
  [[59_ladruno_gradient_concrete_adr]] (contrasting extra-field decision + ndf
  pattern), [[26_ladruno_plane_frontier_adr]] (u-p precedent + PPP),
  [[65_ladruno_explicit_dt_strategies_adr]], [[49a_integrator_scorecard_2026-06-23]],
  upstream `SRC/element/UP-ucsd/` + `SRC/element/UWelements/` (equivalence-gate
  substrate).

> [!note] PDF status (updated 2026-07-10)
> **Received & extracted**: Zienkiewicz & Shiomi 1984 · *Computational
> Geomechanics* 2nd ed. 2022 (full book) · Zienkiewicz, Chang & Bettess 1980 ·
> McGann 2012 + 2015 · Simon, Zienkiewicz & Paul 1984. All §3 formulation claims
> verified against them; benchmark parameters pinned in §7.1.
> **Still optional** (nothing blocking): de Boer, Ehlers & Liu 1993 (B5 + B1–B3
> cover the analytic gates); the explicit trio (Zienkiewicz–Paul–Chan 1988,
> Huang & Zienkiewicz 1998, Huang–Wu–Zienkiewicz 2001) — only needed when P7
> opens; Bochev–Dohrmann 2004 (fallback operator only).

## 11. Adversarial review log (2026-07-10 — pre-code gate, PASSED after repairs)

Three-critic panel (formulation math / framework reality / plan & scope), every
claim verified against tree code or algebra. **The block algebra, sign
conventions, honest-p placement, H̃ operator position, Zone-I validity claim, TH
pairing, heterogeneous-ndf feasibility, registration checklist, and testbed
layout all SURVIVED attack** (refuted findings recorded in the critics'
transcripts). Confirmed defects, all repaired in-place (⟨tags⟩ mark the edits):

- **UP-1 (critical)**: §3.2 table lacked `getResistingForceIncInertia` — the only
  API that carries Qᵀu̇+Sṗ into the transient residual (integrators use
  damp/mass for tangent only). Fixed: three new contract rows + normative
  obligations block.
- **FW-F1 (critical)**: "wrong-solver fails loudly" was false — ProfileSPD (the
  silent default) *discards* a Q block. Fixed: silent-wrong story + mitigations.
- **scope-F1/UP-4 (critical)**: the 1e-6 upstream-equivalence gate was
  mathematically undeliverable (p-as-disp vs p-as-vel Newmark closures coincide
  only at γ=½, β=¼). Fixed: two-leg gate design.
- **UP-3/FW-F2**: the claimed upstream dynamic-seepage precedent is commented-out
  dead code. Fixed: §3.1 corrected; `-dynSeepage` switch; equivalence runs use
  `off`.
- **UP-2/FW-F3**: Rayleigh formation from `getTangentStiff()` would smear −Q/H
  into damping. Fixed: solid-only normative spec + contract test.
- **FW-F4**: the fork's own `Ladruno_NodeResults` PressureSource hardcodes the
  upstream vel-trick (would silently record ṗ). Fixed: §4.4 work item + P1/P4
  gates.
- **FW-F5**: `updateMaterialStage` reaches materials *through the element*;
  forwarding + self-co-registration added to the contract.
- **UP-6/FW-F7**: static well-posedness conditions (k̄>0, ≥1 drained node/region)
  + honest IC recipe replace the overclaims.
- Plus: `-lumped` restored to the surface (scope-F2); CB checkerboard index
  (scope-F3); `-pOrder`/`-Kf` renames + fatal-unknown-flags + legality matrix
  (scope-F4/FW-F8); B1 phasor methodology (scope-F5); loads & defaults paragraph
  (scope-F6); ledger tag-reservation mechanics + cross-registry notes
  (scope-F7/FW-F9); B5 FE spec + û(0,τ) measured-first (UP-5/scope-F8); P0
  inf-sup smoke pinned (scope-F9); `-shape` example corrected to the real (2,6)
  Bernstein-vs-LST collision (scope-F10); getInitialStiff + p-never-enters-
  material contract (scope-F11/FW-F6); B4-T3 leg (scope-F12); TH pressure-basis
  wording (UP-8); equalDOF explicit-list caveat (FW-F10).

## 12. Implementation log

- 2026-07-13 (P4) — **`-dynSeepage` default AMENDED to `off`** (the §4.1 draft
  said `on`; the P1 log pre-authorized this revisit at B5). Evidence, both
  regimes measured: quasi-static consolidation DIVERGES under Δt-refinement
  (P1 ZS84 sweep, err 1.8e-2 → 8.7e-1 as Δt 0.08 → 0.005) and genuine
  dynamics misbehaves on B5 (wandering post-front p level 1.7–2.0 vs β=0.973;
  unbounded shallow-station growth ~1.7e4 by τ=100) — trial-acceleration
  noise feeds f_seep. `-dynSeepage on` stays available (SWANDYNE-parity
  research axis; the +G residual term remains FD-gated in the batteries).
  Companion B5 finding: `-stab auto` injects ~10% spurious deep-station p
  ringing on fast-wave propagation — guide rule: wave-propagation runs use
  `-stab off` (stabilization targets the undrained/checkerboard limit, not
  wave physics).

- 2026-07-13 (P7 pre-study) — **meshless-p / staggered-seam spike**
  (`adr71_meshless_p_spike/`, numpy-oracle idiom, RESULTS.md + 6 plots).
  Question probed: replace the FE pressure interpolation with a meshless
  (MLS) field. Verdicts: **GP-cloud p REFUTED structurally** (pressure space
  4× richer than u ⇒ nP−nU_free = 158 spurious Schur modes at every support
  radius; skeleton locks to 9–27% of reference settlement); node-cloud MLS =
  expensive non-cure (checkerboard persists); centroid cloud = rediscovered
  weaker stabilization. **Genuine payoff isolated and measured**: a fluid
  measure decoupled from solid element life-cycle — after element removal,
  monolithic u-p traps p (0.35q at Tv=1.5) vs persistent fluid drains
  (0.025q, Tv90=0.99). **The staggered architecture was then tried
  end-to-end and works**: fixed-stress split (solid-first sweep + final
  momentum resolve; oedometric L = α²/(K_dr+4G/3)·M_p — 3 iters vs 11 for
  classic α²/K_dr, 0.5× oedometric diverges; naive drained split diverges
  in 4 steps at soil coupling) + persistent FEM pressure overlay ≡
  monolithic BE (crack curve 0.00%, Terzaghi 2.8e-8); both sub-solves SPD
  (symmetric solvers return — offsets the honest-p unsymmetric-tangent
  cost), factorization ×3.3 cheaper, per-step ~2× at 3 iters (3–9
  problem-dependent, degrades toward the incompressible-impermeable limit).
  Implementation trap pinned: a fluid-first sweep with a stale-u predictor
  self-"converges" to p≡0 on the first iterate — any P7 staggered
  integrator must gate against it. Feeds §6/P7; overlay realization beats
  meshless as the starting point (meshless stays the large-deformation
  upgrade path).

- 2026-07-13 — **P7 OPENED as [[73_ladruno_porous_overlay_adr]]**
  (LadrunoPorousOverlay, PATTERN 33022 reserved): fixed-stress staggered
  overlay carrying the §6 staggered route; fractional-step subsumed via L;
  fully-explicit-both-fields rejected for v1. This element stays the
  primary implicit u-p tool; the overlay owns the removal/explicit lanes.

*(filled as phases land; move to `Ladruno_internal/` when complete)*
