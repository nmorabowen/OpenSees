---
title: "ADR 73 — LadrunoPorousOverlay: persistent-fluid staggered u-p (the overlay architecture; opens ADR-71 P7)"
status: draft — planning, NO code
---

# ADR 73 — LadrunoPorousOverlay: a pore-fluid field with its own life-cycle

**Status:** draft. Planning only — no code yet. Opens the ADR-71 §6 / P7
runway (staggered & explicit u-p) with a measured numerical prototype behind
every load-bearing claim: the 2026-07-13 pre-study spike
(`Ladruno_implementation/adr71_meshless_p_spike/` — `meshless_p_toy.py`,
`RESULTS.md`, six plots, shipped in PR #567). Companion theory:
[[71_ladruno_up_family_adr]] (§3 Biot blocks, §6 explicit routes),
`/quake-research`, `/opensees-contact` (element-removal ecosystem),
`/explicit-dynamics` (Δt machinery).

---

## 1. Context & problem

### 1.1 The two things monolithic u-p cannot do

ADR-71 shipped the fork's Biot u-p family (LadrunoUP 33017, P0–P4 complete).
Two structural gaps remain, both on the roadmap's **liquefaction with
discontinuities** target:

1. **Element removal seals the flow path.** In every monolithic u-p element
   (upstream and ours), the element carries BOTH the skeleton and the fluid:
   `remove element` deletes the H/S/Q contributions with the stiffness, so
   the pore fluid loses its continuum the moment the skeleton cracks —
   physically backwards (a gap is water-filled; often it *drains*). Measured
   in the spike (E4, 80 %-width crack at Tv = 0.05): pressure below the crack
   stays trapped at **0.35 q forever** vs 0.025 q truth — a **14× error in
   the unsafe direction** (trapped excess pressure is what drives
   liquefaction failure). No removal-based workflow (sand boils, gap flow,
   progressive collapse below the water table, staged excavation) is honest
   today.
2. **Explicit u-p at fork scale is impossible today** (claim sharpened
   2026-07-13 after NMB's field observation — verified in-tree). Upstream UP
   elements DO run under upstream `CentralDifference`: the rate-DOF trick
   parks S in `getMass()` (p-rows have "mass") and Q/H in `getDamp()`, and
   that integrator factors the full coupled LHS c₂C + c₃M
   (CentralDifference.cpp:155–156) — so it marches, conditionally (needs
   S > 0; dies at the incompressible limit), at the price of a **factored
   non-diagonal matrix + back-solve every step**. That price is the point:
   the fork's explicit stack (CentralDifferenceLadruno / SMS, `system
   Diagonal`, criticalTimeStep, mass scaling — the 18.6 M-element machinery)
   is **matrix-free diagonal-solve** (CentralDifferenceLadruno.cpp:246), and
   no monolithic u-p element can ride it: upstream's S is consistent and its
   coupling lives in damp (never on the diagonal LHS); honest-p LadrunoUP
   has zero-mass p-rows outright — and under *upstream* CentralDifference
   its damp-resident p-row leapfrogs pure diffusion = Richardson-class
   unconditional instability (quirks-ledger row at P1: users will try it).

### 1.2 What the pre-study proved (and refuted)

The spike asked "solve p meshless?" and answered with measurements
(anchors reproduced the shipped ADR-71 P0 inf-sup gates before any claim):

| finding | number |
|---|---|
| GP-cloud meshless p (particles at Gauss points) | **REFUTED structurally**: pressure space 4× richer than u ⇒ nP−nU_free = 158 spurious Schur modes at *every* support radius; skeleton locks to 0.09–0.27× reference settlement |
| node-cloud MLS p | expensive non-cure: checkerboard persists (CB 0.068 ≈ unstabilized FEM); nnz ×4.8, factor ×4.9 |
| the *architecture* behind the idea — a fluid field whose life-cycle is decoupled from solid elements | **REAL**: persistent fluid drains the E4 crack fully (Tv90 = 0.99) where element-tied fluid never does |
| fixed-stress staggered split + persistent FEM overlay (E6) | **≡ monolithic backward Euler**: crack probe curve 0.00 %, Terzaghi 2.8e-8; both sub-solves SPD; factorization ×3.3 cheaper; per-step ~2× at 3 iterations |

The conclusion this ADR builds on: the payoff belongs to the **overlay
architecture** (a separately-discretized, persistent pressure field talking
to the skeleton through the effective-stress seam), realized in **FEM form**
— meshless is demoted to a possible large-deformation upgrade (reserved
phase), and monolithic meshless in any form is dead (do not re-propose;
`RESULTS.md` has the numbers).

### 1.3 Fork facts that shape the design

- **The effective-stress seam already exists and is element-free**: the
  material sees only σ′ (skeleton strain → constitutive); the pore-pressure
  force on the skeleton is exactly the nodal vector −Q·p, and Q is pure
  geometry (ADR-71 §1.4). A solid element never needs to know p exists if
  someone else applies that force.
- **LoadPattern precedent**: H5DRM (PATTERN_TAG 6) is a shipped pattern that
  computes and applies its own nodal force field each step — the exact
  mechanism the overlay needs on the solid side.
- **Domain-hook precedent**: the contact engine attaches to the Domain
  (`Domain::setLadrunoContactDomain`, ADR-39) and re-scans on commit
  (ADR-60) — the life-cycle/removal hook this ADR needs is the same shape.
- **The kernel is already written**: `LadrunoUPKernel.h` + `LadrunoUPShapes.h`
  (pure headers, runtime-sized) assemble Q/H/S/H̃ from geometry alone —
  the overlay reuses them verbatim, zero new block math.
- **Driver precedent**: ADR-45 P3c drives dmumps directly outside the
  standard analysis loop — owning a small solver inside a fork object is
  established practice.

---

## 2. Decision (summary)

1. **One new domain object `LadrunoPorousOverlay`** (LoadPattern subclass,
   `PATTERN_TAG_LadrunoPorousOverlay = 33022` — PATTERN registry; numerically
   = EigenSOE_TAGS_FeastEigenSOE 33022, per-registry, not a collision). It
   owns the pore-pressure field **outside the Domain's DOF graph**: its own
   p vector, its own H/S/H̃/L blocks (assembled once from
   `LadrunoUPKernel.h` over a **geometry snapshot** of a user-named region of
   solid elements), its own drained-BC set, its own small SPD solver.
2. **The region is element-type-agnostic.** The overlay reads only node
   coordinates and connectivity: any 3/4-node 2D or 8-node 3D solid cell
   qualifies — LadrunoQuad/CST/Brick, vanilla quad/stdBrick/SSP, all at
   ndf = ndm, **zero changes to any element anywhere**. (Bézier/quadratic
   cells: reserved, same provider dispatch as ADR-71 §4.2.)
3. **Coupling = fixed-stress staggered split, solid-first** (the spike-pinned
   recipe): per step, the solid analysis (standard OpenSees, any integrator)
   runs with the overlay's frozen nodal forces −Q·pⁿ; at commit the overlay
   advances its fluid system (S\*+L+ΔtH) pⁿ⁺¹ = S\* pⁿ + L pⁿ − Qᵀ Δu and
   refreshes the forces. **L is mandatory** (the drained split diverges in 4
   steps at soil coupling strength — measured). Default
   L = α²/K_dr·M_p (the classic Kim–Tchelepi–Juanes modulus — the variant
   with a general unconditional-stability proof); `-fsL oedometric`
   (α²/(K_dr+4G/3)) is the measured fast opt-in for constrained problems
   (3 vs 11 iters on the column) — sufficient in 1D theory and on both toy
   geometries, but with no general 2D/3D proof and a measured cliff at 0.5×
   (⟨panel A-2⟩: empirical evidence must not masquerade as theory in a
   default). Manual values floor at the oedometric with a loud warning.
4. **v1 lanes**: fs1 single-pass implicit (O(Δt) splitting error, measured
   0.09–0.4 % at production Δt) ships first; the **iterated driver** (fs-k,
   fixed point ≡ monolithic BE) is P2; the **explicit lane** is P3 — central
   difference on u untouched, overlay advances implicitly at commit
   (optionally subcycled): this **is** the Zienkiewicz–Paul–Chan 1988
   implicit-p/explicit-u split and resolves ADR-71 §6's blocker with zero
   changes to the explicit stack (p never enters the integrator).
5. **Fluid life-cycle is the point**: the overlay's fluid measure survives
   `remove element` (geometry snapshot, commit-time rescan). Per-overlay
   policy `-onRemoval keep | drain $kFactor` (water-filled gap at soil k,
   or crack-as-drain with scaled k); the Q coupling of dead cells drops
   automatically (no skeleton to push on).
6. **Division of labor with ADR-71**: monolithic LadrunoUP remains THE tool
   for standard implicit u-p (statics, Taylor–Hood, the B1–B5 battery). The
   overlay is the **discontinuity / explicit / companion** lane. ADR-71's P7
   row ("explicit u-p strategy — own ADR") is **opened by this ADR**.
7. **Out of scope v1**: MP/partitioned domains (the overlay is per-process;
   interface flux across partitions is honest research — flagged, §8),
   quadratic/TH pressure cells, finite strain, partial saturation, the
   meshless fluid realization (P5 reserved upgrade path), thermal analogs.

---

## 3. Formulation (spike-pinned)

### 3.1 The split and its fixed point

Quasi-static Biot (dynamics adds M·ü to the solid row only — the fluid row
is untouched by the split):

- solid:  K u − Q p = f          (material sees σ′; −Q·p applied as nodal load)
- fluid:  Qᵀ u̇ + S\* ṗ + H p = f_seep,  S\* = S + α_stab·H̃ (ADR-71 §3.3 stab retained)

Backward-Euler fluid + fixed-stress iteration k within step n→n+1
(solid-first sweep):

1. u⁽ᵏ⁾ from the standard solid analysis with nodal loads −Q·p⁽ᵏ⁻¹⁾
2. (S\* + L + ΔtH) p⁽ᵏ⁾ = S\*·pⁿ + L·p⁽ᵏ⁻¹⁾ − Qᵀ(u⁽ᵏ⁾ − uⁿ) + Δt·f_seep
3. converged when ‖p⁽ᵏ⁾ − p⁽ᵏ⁻¹⁾‖/(‖p⁽ᵏ⁾‖ + p_scale) < tol; final
   momentum resolve with p⁽ᵏ⁾.

At convergence L cancels: the fixed point **is** the monolithic
backward-Euler solution (spike E6: Terzaghi 2.8e-8, near-undrained 4.5e-7,
crack curve 0.00 %). fs1 = one sweep: stable, O(Δt) splitting error
(0.09 % compressible / 0.4 % near-undrained at the E6 Δt).

**Sequencing is normative, not stylistic** ⟨spike trap⟩: a fluid-first sweep
with a stale displacement predictor **self-"converges" to p ≡ 0 on the first
iterate** (Δu = 0 in ⇒ p⁽¹⁾ = p⁽⁰⁾ ⇒ the tolerance passes vacuously) — it
converges cleanly onto the wrong equation with no residual warning. P1 ships
a contract test that would have caught it (step-load column: p(0⁺) ≈ q, not
0).

### 3.2 The stabilization modulus L (mandatory, bounded)

- **Naive drained split (L = 0) diverges in 4 steps, 10 orders of
  magnitude, in BOTH regimes** (measured): soil coupling strength
  τ = (α²/K_dr)/storage ≈ 10³. There is no `-fsL off`.
- Auto pin ⟨amended by panel A-2⟩: **default L = α²/K_dr · M_p** — the
  classic Kim–Tchelepi–Juanes fixed-stress modulus, the variant carrying a
  general unconditional-stability proof. **`-fsL oedometric`** selects
  α²/(K_dr + 4G/3): measured 3.1–4.6 iters mean vs 11.0 for the classic on
  the constrained column, sufficient by 1D theory and on both toy
  geometries — but with no general 2D/3D proof and a measured cliff (**0.5×
  the oedometric value diverges**), so speed stays an informed opt-in, not
  the default. Moduli per element from the material **initial isotropic
  elastic tangent** (v1: explicit `-moduli`, see the plan) — the `-stab
  auto` machinery of ADR-71 §3.3, same `updateMaterialStage` dirty-cache
  treatment. Manual `-fsL $scale` floors at the oedometric value with a
  loud warning. E7 measures both variants on the footing geometry so the
  P1 default rests on numbers for BOTH constrained and unconstrained cases.
- Known degradation, inherited honestly: the iterated variant's convergence
  rate → 1 toward the incompressible-impermeable limit (measured drift 3 → 9
  iters footing-like; Turska–Schrefler 1993 lower bound on Δt/h² applies).
  fs1 stays stable there (unconditional, KTJ class); accuracy is the Δt
  question, gated in P1/P2.

### 3.3 What the staggered fluid solve is NOT

The overlay's fluid system (S\*+L+ΔtH) is **SPD diffusion — no saddle point,
no inf-sup constraint, no unsymmetric solver**. Two consequences:

- **Symmetric solvers return.** The solid analysis is whatever it was
  (ProfileSPD, MUMPS SYM=2 all legal); the fluid solve is the overlay's own
  small SPD factorization. This claws back exactly what ADR-71's honest-p
  unsymmetric coupled tangent costs (measured: factorization ×3.3 cheaper
  than the coupled unsymmetric solve at 40×40; splu doesn't even exploit
  SPD).
- **Checkerboard is inherited, in degrees** ⟨A-8 precision⟩: iterating to
  convergence reproduces the monolithic pair exactly (its inf-sup character
  included); fs1 is not guaranteed clean either — the QᵀΔu forcing carries
  the unstable mode's imprint, with L acting as strong artificial storage
  that damps it only transiently. So the ADR-71 §3.3 α-stab (H̃) stays in
  S\* for equal-order cells in the undrained limit — same default, same
  `-stab auto <α₀>` surface, both lanes.

### 3.4 Dynamics and the explicit lane

Solid row dynamics ride the solid analysis unchanged (Newmark, HHT, central
difference — the overlay neither sees nor perturbs M, C_Ray, or the
integrator). Fluid row uses Δu of the committed step (no velocity read
needed under BE). `-dynSeepage` semantics inherit ADR-71's P4-amended
default **off** (trial-acceleration noise feeds f_seep — measured in both
regimes, ADR-71 §12).

**Explicit lane (P3)**: central difference advances u with the overlay's
frozen forces; at each commit (or every N steps, `-subcycle $N`, fluid
Δt = N·Δt_explicit) the overlay solves its implicit SPD system. This is the
ZPC-1988 implicit-p/explicit-u class of operator split — with the honest
rider ⟨A-7⟩ that their stability proof covers *their* scheme, not our
frozen-force variant, so the prediction (drained-speed CFL suffices when
the fluid is implicit) is exactly that: a prediction, measured at P0/P3
rather than cited. **Which wave speed governs the CFL under the split
(drained vs undrained) is MEASURED at P3, not assumed** — the P0 toy
extension emulates the explicit lane first and pins the envelope cheaply.

Positioning vs the existing "explicit UP" practice (§1.1): the overlay lane
keeps the solid solve **matrix-free diagonal** (the fork stack untouched,
mass scaling and criticalTimeStep semantics intact) and adds one small SPD
fluid back-solve per commit (subcyclable); the upstream
CentralDifference-with-UP route factors the full coupled c₂C+c₃M and
back-solves the whole model every step, needs S > 0, and owns no removal or
staging story. P3's gate table includes a head-to-head cost/accuracy leg
against that route on the B2 column (it is the honest incumbent, not a straw
man).

**Explicit-lane upgrades (added 2026-07-13 review; P3/P3b below):**

1. **`-fluidUpdate explicit` — the fully matrix-free option (P3b).** The
   *diffusion* CFL of a lumped-S\* forward p-update is Δt ≲ h²/(2d·c_v)
   (d = dimension), and for real soils it is enormously slack: k′ = 10⁻⁷
   m/s, M_oed = 100 MPa, h = 0.5 m ⇒ Δt_diff ≈ 60 s (2D) vs a solid
   explicit Δt ~ 10⁻⁴ s — five-plus orders of margin (even clean gravel at
   k′ = 10⁻³ clears it by orders). The real price sits elsewhere: with the
   fluid explicit, the *coupled* staggered scheme's stability is governed
   by the **undrained** wave speed (the Q·S⁻¹·Qᵀ stiffening — the Xu et
   al. 2021 route ADR-71 §6 catalogued), so the solid Δt_cr shrinks by the
   factor **√(1 + K_f/(n·M_oed))** ⟨panel A-1 correction⟩: for the fork's
   own soft benchmark soils this is ~13× (E = 25 MPa, ν = 0.3) to ~21×
   (E = 10 MPa) — the earlier "5–8×" holds only for stiff/dense profiles
   (M_oed ≳ 150 MPa). Mass scaling composes against it; the per-element
   Δt_cr advisory (item 3) computes the true factor, so nothing downstream
   hardcodes a range. What the option buys: **no factorization anywhere in
   the run** — the fluid step is a local axpy + one halo exchange — making
   the 18.6 M-class saturated explicit run feasible and **dissolving most
   of the §8 MP risk on this lane**. Implicit stays the default; explicit
   is a gated option (measured envelope, mass-scaling composability, L not
   needed — no iteration).
2. **Pipelined implicit fluid (P3 design note).** fs1 tolerates O(Δt) lag by
   construction, so the implicit fluid solve may run on its own thread over
   the next subcycle window, forces updating one window late — the overlay
   owns all its data structures (no Domain locks), so this is natural
   concurrency. Combined with `-subcycle`, the fluid wall-clock approaches
   zero.
3. **Overlay-aware Δt_cr advisory (P3 gate item).** criticalTimeStep/SMS
   price the pencil from the *drained* skeleton and cannot see the overlay;
   if the split's stability is undrained-governed (P0 measures this), the
   advisory is √(1 + K_f/(n·M_oed)) optimistic — ~13–21× for soft soils
   ⟨A-1⟩ ⇒ mysterious blow-ups. The overlay exposes the per-element factor;
   the ADR-65 advisory machinery multiplies it in; SMS then scales against
   the correct pencil.
4. **`-subcycle auto`**: N computed at setup from the diffusion time scale
   (h²/c_v vs Δt_explicit) instead of asked of the user.

### 3.5 Initialization — better than the monolithic recipe

The overlay owns H, so it can solve its own steady state: `-pInit steady`
(H p = f_seep with the drained set), `-pInit hydrostatic $gw <$phreatic>`,
or explicit nodal values. No `InitialStateAnalysis` interplay, no
staged-`sp` Transformation trap (ADR-71 P1 quirk) — the solid analysis never
sees p DOFs at all.

---

## 4. Architecture

```
SRC/domain/pattern/ladrunoPorousOverlay/
  LadrunoPorousOverlay.{h,cpp}   # LoadPattern subclass: fluid state (p, blocks,
                                 #   drained set, removal mask), advance(), SPD solve,
                                 #   nodal-load management, sendSelf (serial DB v1)
  OPS_LadrunoPorousOverlay.cpp   # parser (pattern command), fatal unknown flags
  CMakeLists.txt
# reuses verbatim (no copies):
#   SRC/element/ladrunoUP/LadrunoUPKernel.h   (Q/H/S/H̃ blocks)
#   SRC/element/ladrunoUP/LadrunoUPShapes.h   (T3/Q4/H8 providers, GP rules)
tests/
  (P0) adr71_meshless_p_spike/meshless_p_toy.py grows E7: exact v1 sequencing
       (fs1-no-resolve error), explicit-lane emulation + CFL envelope,
       subcycle ratio study, removal-step jump consistency
```

### 4.1 Command surface (draft — pinned at P1)

```tcl
pattern LadrunoPorousOverlay $tag \
    -region {$ele1 $ele2 ...} | -regionAll   ;# geometry snapshot of these solids
    -Kf $Kf -poro $n -rhoF $rhof \
    -perm $k1 $k2 <$k3>                       ;# k̄ = k_hydraulic/γw (ADR-71 convention)
    <-permH $k1 ... -gammaW $gw>              ;# sugar, divides internally
    <-layer {$eleSet} -perm $k1 $k2 <$k3> <-poro $n>>
                                              ;# per-layer override, repeatable —
                                              ;#   k AND porosity (storage n/Kf)
                                              ;#   vary by layer; layered profiles
                                              ;#   are ONE overlay (assembly is
                                              ;#   per-cell; free). Kf/rhoF stay
                                              ;#   overlay-global (same water)
    <-alpha $biotAlpha> <-Ks $Ks> \
    <-thick $t>                               ;# 2D regions, default 1.0 — must match
                                              ;#   the solid elements' thickness ⟨A-4⟩
    -drained {$nd1 $nd2 ...}                  ;# p-fixed set (≥1 per connected region
                                              ;#   for statics — ADR-71 §3.2 rider)
    <-pInit steady | hydrostatic $gw <$z0> | {$nd $val ...}> \
    <-stab auto <$alpha0> | off | $alpha>     ;# H̃ in S*, ADR-71 §3.3 semantics
    <-fsL classic | oedometric | $scale>      ;# default classic α²/K_dr (proven,
                                              ;#   ⟨A-2⟩); oedometric = fast opt-in;
                                              ;#   $scale floors at oedometric; no off
    <-staticMode hold | steady>               ;# STATIC-analysis commits: domain
                                              ;#   "time" is the LOAD FACTOR there
                                              ;#   ⟨A-3, §8⟩ — hold p (default) or
                                              ;#   re-solve steady seepage; transient
                                              ;#   fluid marching only under
                                              ;#   transient analyses
    <-onRemoval keep | drain $kFactor>        ;# fluid life-cycle policy
    <-fluidUpdate implicit | explicit>        ;# default implicit (SPD solve);
                                              ;#   explicit = lumped-S* forward step,
                                              ;#   matrix-free (P3b, §3.4 — undrained
                                              ;#   CFL governs the SOLID Δt then)
    <-subcycle auto | $N>                     ;# fluid advance every N commits
                                              ;#   (auto: N from h²/c_v vs Δt)
    <-fluidBody $b1 $b2 <$b3>> <-dynSeepage on|off>   ;# default off (ADR-71 §12)
```

- The pattern's load factor is **ignored by design** (the overlay manages its
  own force amplitudes; a TimeSeries on pore physics is meaningless) — parser
  says so once.
- **Unknown flags are parser-FATAL** (ADR-71 house rule — a mistyped porous
  flag silently changes physics).
- Iterated driver (P2): `LadrunoStaggeredAnalyze $nSteps $dt <-tol $t> <-maxIter $k>`
  — wraps revertToLastCommit + re-analyze + overlay advance until the §3.1
  test passes; plain `analyze` gives fs1.

### 4.2 The seam contract (what crosses, when)

| direction | payload | mechanism | timing |
|---|---|---|---|
| overlay → solid | nodal forces +Q·p (compression-positive p; solid gravity stays user-side with mixture ρ, ADR-71 §3.5 convention) | `node->addUnbalancedLoad()` inside the pattern's `applyLoad` override — the VERIFIED H5DRM mechanism (H5DRMLoadPattern.cpp:897) | applied **once per step at newStep** ⟨A-5 mechanism correction⟩: transient integrators call `Domain::update(newTime,dT)` → `applyLoad` (Domain.cpp:2381); iteration-time `updateDomain()` skips load application, and `formUnbalance` only READS nodal unbalance — so forces are constant within the step by construction. Values refreshed at commit are picked up at the NEXT newStep |
| solid → overlay | Δu = uⁿ⁺¹ − uⁿ of region nodes (committed) | direct Node reads | at advance(), post-commit |
| framework → overlay | advance() trigger; removal rescan | Domain::commit hook at the ADR-39 contact-hook slot (Domain.cpp:2233 — after node/element commitState, **before** the `dT = 0` reset at :2247, so committed disps are final AND Δt is still readable ⟨A-5⟩); one strictly-additive `// Ladruno` block — fallback: applyLoad new-step detection if the hook is refused at review | each commit (or every N with -subcycle) |
| element death → overlay | region cell k with no live element | commit-time rescan of region tags | Q-mask update; H/S per -onRemoval |

Everything else is private: p never appears in the DOF graph, recorders see
it only through the overlay's own response channels (P4), and no shipped
element, integrator, handler, or numberer changes.

### 4.3 Fluid solver

Region-sized SPD system, factored once per (Δt, stage, removal-state) and
reused. v1 pins the implementation at P1 from: (a) in-object banded
Cholesky (regions are meshes of one soil body — bandwidth fine), (b) MKL
PARDISO (present in all fork builds), (c) reuse of an OpenSees symmetric
SOE standalone. Bias: (a) for zero new dependencies; measured spike sizes
(≤ few 10⁴ p-unknowns) make this a non-critical choice — revisited only if a
profile says so (ADR-40 discipline).

### 4.4 Registration checklist

`SRC/classTags.h` (PATTERN_TAG 33022 + comment), pattern dispatch in the two
command registries, `FEM_ObjectBrokerAllClasses` pattern factory,
`ladrunoPorousOverlay/CMakeLists.txt` + parent include, LADRUNO header stamp
(`stamp_headers.py` GLOBS), `LEDGER_implementations` amendment,
`banner_features.txt` line (at ship), vanilla-ledger rows for the Domain
commit hook if taken.

---

## 5. Class tag & reservation

- **`PATTERN_TAG_LadrunoPorousOverlay = 33022`** (PATTERN registry).
  Cross-registry occupants of 33022, stamped per house precedent: EigenSOE
  FeastEigenSOE — numerically equal, different registry, not a collision.
  ELE-registry neighbors for orientation: 33021 = LadrunoCSTPair (last ELE),
  33019 pencilled H27.
- **Reservation mechanics** (ADR-71 ⟨scope-F7⟩ precedent): **this ADR's PR
  lands the `LEDGER_implementations` row** `LadrunoPorousOverlay —
  PATTERN 33022 — RESERVED (ADR-73), not yet built`; the classTags.h
  `#define` lands at P1.

---

## 6. Relation to ADR-71

- **Opens P7**: ADR-71 §7's P7 row ("explicit u-p strategy — own ADR, see
  §6") is this document. The three §6 routes map: staggered partitioning →
  §3.1/§3.4 (adopted, ZPC-1988 realization); fractional-step → subsumed (the
  overlay's split adds the same stabilizing Laplacian through L);
  fully-explicit-both-fields (Xu 2021) → **adopted as the gated P3b option**
  `-fluidUpdate explicit` (§3.4: diffusion CFL measured slack by orders for
  real soils; the honest price is the undrained-speed CFL on the solid,
  priced into the ADR-65 advisory) — not the v1 default.
- **Monolithic LadrunoUP is not deprecated by this**: implicit statics,
  steady seepage, Taylor–Hood, and the validated B1–B5 battery stay on the
  element. The overlay's gates cross-check against it (it is the reference
  operator — spike E6 methodology).
- The kernel/providers are shared headers; any future ADR-71 P5 shape
  (Q9/H20) becomes overlay-legal for free through the same provider dispatch.

---

## 7. Phases & exit gates

| Phase | Scope | Gate |
|---|---|---|
| **P0** | toy extension E7 (numpy, no OpenSees build) | exact v1 sequencing pinned: fs1-without-final-resolve error vs monolithic measured over Δt-halving (expect O(Δt); pins the P1 default); explicit-lane emulation (CD solid + implicit fluid): stability envelope measured — drained-vs-undrained CFL governance answered with numbers + subcycle-N degradation curve; removal-step jump consistency (Q mask flips mid-march, no spurious p transient beyond the physical Mandel-Cryer dip); L-floor sweep reproduced on the dynamic path — **DONE, PR #575 (2026-07-14): three predictions flipped, §12 P0 entry governs (rate-form fs1; `-substep` removed; P3 L=0 @ undrained pencil)** |
| **P1** | `LadrunoPorousOverlay` fs1 implicit end-to-end (Q4/T3/H8 regions) | Terzaghi vs analytic ≤ pinned tol at production Δt; **two-leg gate vs monolithic LadrunoUP** (ADR-71 methodology: mutual Δt-convergence, observed order ≥ 1; the exact-equality leg belongs to P2's iterated driver); **E4-in-OpenSees**: real `remove element` crack, both `-onRemoval` policies, curves vs the toy reference; **sequencing contract test** (step-load column p(0⁺) ≈ q — kills the fluid-first p≡0 trap class); L auto twin-checked vs oracle incl. `updateMaterialStage` dirty path; pattern-timing verification (forces bit-constant across a step's Newton iterations); serial DB round-trip; **full adversarial panel** (per [[feedback_adversarial_gate_when]]: core-touch — Domain hook + pattern semantics — and novel framework integration; the split math itself is spike-validated) — **DONE, PR #576 (2026-07-15): all gates green, panel PASS (§11 + §12 P1 entries); xfail-6c tripwire → upstream rho bug found+fixed #577, 6c now a hard gate** |
| **P2** | iterated driver `LadrunoStaggeredAnalyze` (fs-k) | fixed-point gate: iterated staggered ≡ monolithic LadrunoUP same-Δt (E6 leg, target ≤ 1e-6 rel); B4 footing staggered (CB gate under stab); PDMY staged liquefaction column vs the ADR-71 P4 monolithic reference (stage transport → L/stab cache dirty); iteration-count telemetry (mean/max per step) exposed — **DONE, PR #580 (2026-07-18): all gates green with ONE gate transposition (§12 P2 entry item 1 — the ≤1e-6 exact-equality leg is against monolithic BE = the frozen toy `qs_mono` [measured 3.4e-7 p / 1.1e-8 u, both regimes]; vs the real LadrunoUP it is a mutual-Δt-convergence leg [orders 0.99/0.99], because NO stock integrator reproduces BE rates on the fluid row); Opus ×3 panel PASS-WITH-FIXES, all landed; rider: `-pInit` list DB restore CLEAN post-#577 → promoted to hard gate** |
| **P3** | explicit lane | CD + overlay on B2/ZS84 column vs implicit monolithic (two-leg); **incumbent head-to-head** (§3.4): same column under upstream CentralDifference + FourNodeQuadUP — accuracy AND per-step cost/scaling vs the overlay's diagonal-solid + small-SPD-fluid path, plus the S→0 failure demo; measured CFL envelope vs the P0 pins; **overlay-aware Δt_cr advisory gated** (§3.4 item 3: per-element undrained factor √((M+K_f/n)/M) wired into the ADR-65 machinery — naive drained advisory demonstrated 5–8× optimistic, then corrected); `-subcycle auto|$N` sweep gate; pipelined-fluid design note carried (§3.4 item 2 — implementation may land here or P3b); **LEDGER_quirks row**: honest-p LadrunoUP + upstream CentralDifference = Richardson-unstable p (leapfrogged diffusion) — loud test pinning the symptom; energy-balance advisory channel (ADR-69: overlay work terms enter the closure residual — documented, not silently absent) |
| **P3b** *(option)* | `-fluidUpdate explicit` — fully matrix-free lane (§3.4 item 1, Xu-2021 class) | dual-CFL gate: diffusion limit verified slack by orders on realistic k′ sweep; measured undrained-CFL envelope vs the √((M+K_f/n)/M) prediction; equivalence vs implicit-fluid lane at matched Δt (both are O(Δt) splits); mass-scaling composability (SMS against the undrained pencil); MP smoke on the halo-exchange update (the §8 MP-risk dissolution demonstrated, not asserted) |
| **P4** | ecosystem | overlay p-field recorder channels (own response surface — nodal-DOF recorders can't see it; LadrunoRecorder/Monitor topology rows per [[06_quadrature_global_gp_plan]]); user guide (`LadrunoPorousOverlay_guide.md`) incl. division-of-labor table vs LadrunoUP + init recipes; banner line + ledgers → shipped; apeGmsh emitter runway note (companion repo item) |
| **P5** *(reserved)* | meshless/MPM fluid realization (large-deformation upgrade path); MP/partitioned-domain story | own mini-ADR; the spike's E1–E3 refutations bound what any meshless realization must prove first |

Each phase is one PR off `ladruno`, ledgers updated in-PR.

---

## 8. Risks / open questions

- **Pattern-timing assumption — now VERIFIED, restated precisely** ⟨A-5⟩:
  loads are applied once per step at newStep (Domain.cpp:2060/2381), persist
  in nodal unbalance across Newton iterations (formUnbalance reads, never
  re-applies), and the overlay's commit-refreshed values reach the solid at
  the next newStep. P1's bit-constancy test observes `nodeUnbalance` between
  iterations under Static, Newmark, AND central difference. Residual risk is
  now only exotic integrators that call updateDomain(time,dT) mid-step
  (ArcLength family does — documented as unsupported with the overlay at
  P1); fallback remains the contact-engine mechanism (ADR-39
  backing-element-less FE_Elements, incl. the `formEleTangent` routing
  quirk).
- **Static analyses: domain "time" is the LOAD FACTOR** ⟨A-3, CONFIRMED in
  code⟩: LoadControl::newStep calls applyLoadDomain(currentLambda)
  (LoadControl.cpp:130), and Domain::applyLoad sets currentTime = λ
  (Domain.cpp:2060–2065) — an fs1 fluid advance at a static commit would
  integrate the mass balance over Δt = Δλ, silently wrong physics. Repair
  is the `-staticMode hold | steady` surface (§4.1): fluid frozen (forces
  still applied from committed p) or steady-seepage re-solve per static
  commit; transient marching only under transient analyses. The
  gravity-init recipe in the guide sequences: static stage with
  `-staticMode steady` → transient consolidation/shaking.
- **MP absent in v1** — honest but now lane-dependent: the *implicit*-fluid
  overlay needs interface flux exchange across partitions
  (Schur/Neumann-Neumann class) — genuinely research, parked in P5. The
  **explicit-fluid lane (P3b) largely dissolves this**: a local lumped
  update + halo exchange is the same communication pattern the explicit
  solid already uses — P3b's MP smoke demonstrates it rather than asserting
  it. Serial implicit (the 60 GB emit-host class) is served either way.
- **UniformExcitation**: üg enters the solid row via the standard R·üg
  path; the fluid row's ρ_f·üg seepage drive stays omitted (Chan-1988
  class, consistent with ADR-71 P4's `-dynSeepage off` default) — guide
  documents.
- **Iterated-variant degradation** toward the incompressible-impermeable
  limit (rate → 1; measured 3 → 9 iters; Turska–Schrefler Δt/h² floor).
  fs1 stays stable; P2's telemetry makes the cost visible instead of silent.
- **Removal-step consistency**: Q changes discontinuously at removal; the
  toy accepted the one-step jump (physical stress redistribution dominates).
  P0 pins that no *numerical* artifact rides it; the `-onRemoval drain`
  variant re-factors the fluid operator (cheap, region-sized).
- **Recorder story**: overlay p is not a nodal DOF — nothing existing
  records it. P4 owns this; interim: overlay Print + monitor channel.
- **Solid gravity convention**: user applies mixture-ρ self-weight on solid
  elements (unchanged family convention); the overlay adds only seepage
  body terms. The guide repeats ADR-71's quadUP `<b1 b2>` 2×-trap note.
- **One overlay per hydraulically connected water body** (modeling rule,
  guide-normative): separate overlays own separate p fields — there is no
  flow between them by construction. A layered soil profile (the common
  apeGmsh emission case, flagged in the 2026-07-13 review) is therefore ONE
  overlay spanning all layers, with per-layer k and porosity via `-layer`
  blocks (storage n/K_f is per-cell anyway); the L and α-stab moduli are
  already per-element (material initial tangent), so heterogeneity costs
  nothing. Deliberately-separate aquifers/perched zones are the legitimate
  multi-overlay case.

---

## 9. Alternatives considered (and why rejected)

| Option | Verdict | Killer fact |
|---|---|---|
| (a) monolithic meshless-p (any cloud) | ❌ | spike E1–E3: GP cloud structurally spurious (158 modes, radius-independent) + locking; node cloud = non-cure at ×4.8 cost; BCs non-interpolatory |
| (b) staggered as a new Integrator solving masked subsets of one SOE | ❌ | OpenSees has one DOF graph per analysis; subset re-solves need numberer/handler surgery across core classes — maximal blast radius for zero physics gain |
| (c) two-model script coupling (two domains in Tcl/Python) | ❌ | interpreter/domain is a singleton in practice (openseespy hard-singleton); state exchange via files/arrays per step is fragile and unserviceable; no removal hook |
| (d) pressure-overlay as real p-DOF elements solved monolithically | ❌ | that *is* LadrunoUP — re-ties fluid life to element life (loses E4), keeps the unsymmetric coupled tangent, and forecloses the explicit lane |
| (e) **LoadPattern overlay engine + fixed-stress staggered (chosen)** | ✅ | zero element changes; fluid life-cycle decoupled (E4 measured); SPD sub-solves (symmetric solvers return, factor ×3.3); explicit lane = ZPC-1988 for free; every load-bearing number already measured in the spike |

---

## 10. References

**Measured basis (primary):**
`Ladruno_implementation/adr71_meshless_p_spike/` — `meshless_p_toy.py` E1–E6,
`RESULTS.md`, plots (PR #567); ADR-71 §12 log entry 2026-07-13.

**Split theory:**
- Kim J, Tchelepi HA, Juanes R (2011). *Stability and convergence of
  sequential methods for coupled flow and geomechanics: drained and
  undrained splits* + *fixed-stress and fixed-strain splits.* CMAME
  200(23–24):2094–2116 and 200(13–16):1591–1606 — the drained-split
  instability and the fixed-stress cure this ADR measured.
- Mikelić A, Wheeler MF (2013). *Convergence of iterative coupling for
  coupled flow and geomechanics.* Comput Geosci 17:455–461 — contraction
  proof for fixed-stress.
- Turska E, Schrefler BA (1993). *On convergence conditions of partitioned
  solution procedures for consolidation problems.* CMAME 106:51–63 — the
  Δt/h² floor (iterated variant).

**Explicit runway (inherited from ADR-71 §6):**
- Zienkiewicz OC, Paul DK, Chan AHC (1988). IJNME 26(5):1039–1055 —
  implicit-p/explicit-u staggered (the P3 lane).
- Huang M, Zienkiewicz OC (1998). IJNME 43(6):1029–1052 (fractional step);
  Xu et al. (2021) SDEE 141:106452 (fully explicit — rejected v1).
- *Computational Geomechanics* 2nd ed. (2022) §3.2.4, §5.5 ✅ in hand.

**Fork precedents:** [[71_ladruno_up_family_adr]] (kernel, conventions,
benchmark dossier), ADR-39/60 contact engine (Domain hook,
backing-element-less FE fallback), H5DRM pattern (self-managed nodal loads),
ADR-45 P3c (fork-owned solver driver), ADR-65 (Δt machinery),
ADR-69 (energy closure), [[feedback_adversarial_gate_when]].

## 11. Adversarial review log

**2026-07-13 — pre-P0 full sweep of the ADR + implementation plan** (run
inline by MAIN with file:line verification — the Opus 3-critic panel was
infra-blocked; the P1 code panel remains scheduled per §7 and should re-attack
anything marked empirical here). Findings, all repaired in-place (⟨A-n⟩ tags
mark the edits):

- **A-1 (MAJOR, math)**: the "5–8×" undrained-CFL factor was wrong for the
  fork's own soft benchmark soils — √(1+K_f/(n·M_oed)) = ~13× (E=25 MPa) to
  ~21× (E=10 MPa); 5–8× holds only for M_oed ≳ 150 MPa. §3.4 items 1/3
  corrected to the formula + honest ranges.
- **A-2 (MAJOR, theory-vs-measurement)**: the oedometric L default rested on
  two toy geometries; only the classic α²/K_dr carries a general
  unconditional-stability proof (KTJ), and the measured 0.5× cliff shows the
  boundary is near. Default flipped to classic; `-fsL oedometric` is the
  measured fast opt-in; E7 gains the footing-geometry iteration comparison.
- **A-3 (MAJOR, framework — confirmed in code)**: under static analyses the
  domain "time" is the load factor (LoadControl.cpp:130 →
  Domain.cpp:2060–2065), so a static-commit fluid advance would use Δλ as
  seconds. New `-staticMode hold|steady` surface + §8 entry + guide-recipe
  sequence.
- **A-4 (MAJOR, surface)**: plan's pinned command had dropped the ADR's
  `-alpha/-Ks` and BOTH docs lacked `-thick` for 2D regions (dv needs it;
  the overlay cannot read thickness from arbitrary solid elements). Added.
- **A-5 (MAJOR, mechanism precision)**: "re-applied each Newton iteration"
  was mechanically wrong — loads apply ONCE per step at newStep
  (Domain::update(newTime,dT) → applyLoad, Domain.cpp:2381); formUnbalance
  only reads nodal unbalance; the commit hook must sit at the ADR-39 slot
  (Domain.cpp:2233) BEFORE the dT=0 reset (:2247). Seam table rewritten;
  bit-constancy test pinned to `nodeUnbalance`; ArcLength-family
  (mid-step applyLoadDomain callers) documented unsupported at P1.
- Minor: A-6 diffusion-CFL constant h²/(2d·c_v) (2D example 125→~60 s;
  margin argument unchanged); A-7 ZPC-1988 proof does not cover the
  frozen-force variant — prediction hedged; A-8 checkerboard wording (fs1
  not guaranteed clean pre-fixed-point); A-9 Richardson claim CONFIRMED
  rigorous (CD equilibrium at tₙ + leapfrog velocity ⇒ midpoint rule on
  S·ṗ+H·p); A-10 uSnapshot_ NOT serialized — re-derived from committed
  node disps at recvSelf; A-11 test files/Zone pinned; A-12 this section
  created as the log home; A-13 region-overlap guards (two overlays / 
  LadrunoUP elements inside a region → fatal) pinned into WP1.B + battery.
- **Refuted attacks (survived)**: H5DRM mechanism matches the plan
  (addUnbalancedLoad, H5DRMLoadPattern.cpp:897); kernel/provider signatures
  match exactly (LadrunoUPKernel.h:87–203, LadrunoUPShapes.h); Domain::commit
  hook slot real and correctly ordered (Domain.cpp:2213/2233/2247); PATTERN
  33022 free (registry holds 1–7); broker = getNewLoadPattern
  (FEM_ObjectBrokerAllClasses.cpp:2717); fixed-point algebra; sign
  conventions consistent ADR/plan/toy; element removal retains nodes.

**2026-07-14 — P1 code panel (Opus ×3, on the built module)**. Batteries green
first (physics 6/6, framework 8/8 + one xfail); then the three critics.
Verdicts: **math/contract PASS** (no bug; its one hypothesis — removal-forced
partial-window `dpCommitted_` mis-scale — REFUTED by direct measurement: the
unscaled reference holds error below the surrounding splitting error, and the
"corrected" rescale it proposed BLOWS UP → the literal contract rule is right).
**framework-reality** and **robustness** both PASS-WITH-FIXES; all fixes landed
this PR:

- **P-1 (MAJOR, framework — ⟨A-3⟩ trap NOT covered by the shipped guard)**: the
  earlier code assumed a static run yields commit `dt<=0`; FALSE — `LoadControl`
  hands `dT = Δλ > 0` (LoadControl.cpp:127-130 → Domain.cpp:2065-2066), so
  default `SM_MARCH` silently consolidates the fluid on load-factor pseudo-time
  during a gravity stage (reproduced: node p 24.8→49.9 over λ 0.1→1.0, no
  warning). A LoadPattern cannot see the integrator, so auto-detect is
  impossible. FIX: loud one-time advisory on the FIRST real march (per-instance
  `marchNoticeShown_`) naming dt + `-staticMode hold|steady`; the false
  `dt<=0`-detects-static comment/branch rewritten.
- **P-2 (MAJOR, robustness)**: a repeated tag in `-region` built TWO cells for
  one element → silent 2× stiffness/force/storage. FIX: fatal on duplicate
  region tag; also fatal on an element in two `-layer` blocks and on a `-layer`
  element outside `-region`; per-cell `poro > α` fatal (the global check missed
  per-layer `-poro`).
- **P-3 (MINOR)**: setDomain claimed "INERT" but `onDomainCommit` returned −1
  (aborts every commit) — message corrected to state the analysis aborts until
  the region is fixed (fail-loud kept). recvSelf now rejects negative/short
  payload counts (corrupt-datastore guard). Per-instance advisory latch
  replaces the process-`static dtWarned`.
- **1.D BLOCKER RESOLVED — upstream, not overlay** (deformable DB restart
  non-reproducible, test 6c): reproduced overlay-FREE — a body-force solid +
  ONE ordinary `Plain` pattern with a nodal load restarts non-deterministically
  (3/8 fresh procs, rel up to 3e-2; restored operator eigenvalues go NEGATIVE),
  while the same model with NO load pattern is bit-exact 8/8. The overlay's own
  restored state (p/dp/uSnapshot values, per-node d/v/a, SPs) is verified exact
  in diverging runs — it merely AMPLIFIES the upstream FileDatastore/Domain
  transient-restart corruption through its large +Q·p forces. Deformable
  transient DB restart documented UNSUPPORTED pending an upstream-fix decision
  (LEDGER_quirks); frozen-skeleton restart is bit-exact (proves overlay
  serialization incl. `dpCommitted_` is correct). Also fixed en route: a
  same-size FileDatastore clobber (two Vectors on one dbTag) via a second
  payload dbTag (`dbTag2_`), and moved uSnapshot_ re-derivation to lazy
  first-use (restore-time node state is not ordering-safe).
  **→ RESOLVED 2026-07-15 (PR #577)**: root cause = quad/tri elements never
  serialized their ELEMENT-LEVEL `rho` + uninitialized broker ctor → garbage
  mass matrix on every DB/MPI-restored element (the "LoadPattern in the
  stream" correlate was an artifact — the no-pattern control was nearly
  motionless, so the corrupt M went unexercised). The overlay-exoneration
  adjudication above stands verified. Test 6c promoted xfail → HARD GATE per
  its XPASS contract (reconcile PR after #577); deformable transient DB
  restart is SUPPORTED again fork-wide.
- **Refuted at P1 panel**: revertToLastCommit fluid-rollback (committed fluid
  state mutates only inside a successful commit; test 3 bit-identical);
  per-Newton force re-application (applyLoad once/newStep; bit-constancy green
  Static/Newmark/CD); hook ordering (Domain.cpp:2254, before the dT reset);
  near-incompressible ν=0.4999 divergence (L=α²/K_dr scales WITH the split it
  stabilizes); all-cells-dead / k̄=0 / re-added-tag / detJ≤0 (loud SPD abort or
  sentinel-inert, never NaN).

## 12. Implementation log

### P2 — iterated driver `LadrunoStaggeredAnalyze` + stage transport (PR #580, 2026-07-18)

`LadrunoStaggeredAnalyze $n $dt <-tol $t> <-maxIter $k> <-pScale $s> <-verbose>`
(+ `-stats` telemetry query) live in BOTH interpreters; shared core
`LadrunoStaggeredDriver.{h,cpp}` decomposes analyzeStep (analysisStep /
checkDomainChange / newStep / solveCurrentStep) with an inner
{revertToLastCommit + revertToLastStep + re-solve + advanceTrial(false)} loop,
final momentum resolve, all-or-nothing abort. Overlay gained the
external-stepping latch (applyLoad→`+Q·pTrial_`; latched commit =
commitFluid-only sync — the double-advance guard), the free-system residual
(`lastAdvanceRelChange`), and parameter-route moduli transport
(`parameter $p loadPattern $tag E|nu|layerE $i|layerNu $i` →
`moduliDirty_` → lazy `rebuildModuliCaches` re-assembling aS_/aL_ through the
SAME `assembleStorageAndL` path the snapshot uses). Battery
`tests/test_ladruno_overlay_driver.py` 12/12 green; P1 batteries re-run green
(anchor bit-preserved 3.46e-7). **Adjudications/refutations (this entry
governs where it amends §3/§7 wording):**

1. **§7's "iterated ≡ monolithic LadrunoUP same-Δt ≤ 1e-6" is unattainable
   as written and was TRANSPOSED.** The driver's fixed point is single-step
   backward-Euler (§3.1, L cancels — confirmed in-code by the panel), but NO
   Newmark parameterization reproduces BE rates on the fluid row (γ/β = 1
   and γ/2β = 1 are inconsistent), and in-tree `BackwardEuler` is actually
   BDF2 with a trapezoidal bootstrap. The exact-equality leg therefore gates
   against monolithic BE = the frozen toy `qs_mono` on the same stabilized
   operators: **measured 3.44e-7/3.39e-7 (p) and 1.1e-8 (u), compressible /
   near-undrained, at `-tol 1e-10`** — the honest E6 transposition. The real
   LadrunoUP is covered by the P1 mutual-Δt-convergence methodology
   (orders 0.99/0.99, diffs shrinking 6.8e-3 → 1.7e-3).
2. **Driver iteration counts land slightly ABOVE the frozen e76 means, within
   band** — column classic 11.47 mean / 19 max (e76 11.25/18), column oed
   4.25 (e76 3.29), footing classic 5.05 (e76 4.35); no divergence, no
   maxIter hits anywhere. The pins' "the rate-form warm start may only
   IMPROVE on e76" expectation was NOT borne out: the driver's residual is
   the free-system (drained-rows-excluded) norm — STRICTER than the toy's
   full-vector denominator (panel math-6 fix) — and the C++ models are not
   the toy meshes. Bands (≤1.3×) hold with margin.
3. **The unstabilized near-undrained checkerboard STALLS the fixed-point
   iteration** (new measurement): `-stab off` on the B4-like footing plateaus
   the residual at ~1.5e-5 > tol 1e-6, so the all-or-nothing driver refuses
   every step at production tol — itself a loud symptom of the CB mode. The
   CB gate's control runs at tol 1e-4 (above the stall floor); stab-auto vs
   stab-off roughness ratio 5.3× ≥ 3.
4. **PDMY staged column tracks the monolithic ADR-71 P4 reference to
   rel_p 6.1e-3 / rel_u 5.5e-3** at the same Δt (gravity elastic stage →
   `updateMaterialStage` flip → parameter-route `E` re-set scaled by the
   probed constrained-modulus ratio → response march). The parameter-route
   bit-twin gate (update-before-march vs constructed-with-E₂) is exactly
   0.0. MODELING TRAP promoted to LEDGER_quirks: plain quad `b1 b2` is body
   FORCE/VOLUME while `LadrunoUP -body` is an ACCELERATION (×ρ_mix inside
   the element) — the staggered twin's quad needs `b = ρ_mix·accel`, else
   the overlay's hydrostatic `+Q·p` cancels the under-scaled weight to ~zero
   (first battery run measured settlement 1e-9 vs 3.5e-4, rel_u exactly 1.0).
5. **Opus ×3 panel: PASS-WITH-FIXES ×3, all fixes landed:** integrator
   classTag WHITELIST (Newmark/HHT/GeneralizedAlpha; TRBDF2's unrestored
   two-step history silently degrades to trapezoidal per iteration and the
   explicit family inherits a no-op revertToLastStep — loud fatal now),
   free-system residual denominator, NaN-rejecting parameter guards,
   `tol/maxIter/pScale` validation + non-finite-residual abort (a degenerate
   `pScale ≤ 0` could previously NaN-converge silently), catch-up sync
   record row, `layerNu = 0` sentinel-collision rejection. Driver omissions
   documented (not supported): `_RELIABILITY` per-step sensitivities,
   `numSubLevels` subdivision.
6. **2.B implementation deviation (upheld, panel-verified):** the coupling
   term reads TRIAL displacements under the latch (`readRegionDisp(u,
   trial)`) — mid-step `getDisp()` returns the PREVIOUS committed disp, so
   the pinned seam as-frozen would have produced Qᵀ·Δu ≈ 0 and a dead
   fixed point. Unlatched paths still read committed (P1 bit-identical;
   anchor re-measured 3.46e-7 post-driver in-process).
7. **Rider RESOLVED:** `-pInit` list FileDatastore restore is CLEAN
   post-#577 (the P1-era crash was the same upstream element-rho corruption)
   — promoted from observation to a HARD gate (battery (f)); the P1
   framework battery's MEASURED FINDINGS note updated.

Deferred: P3/P3b explicit lanes, P4 recorder channels + guide (the guide
inherits the staggered-twin body-force recipe and the stage-flip moduli
recipe from this entry).

### P1 — fs1 rate-form implicit overlay end-to-end (PR #576, merged 2026-07-15)

`LadrunoPorousOverlay` live: PATTERN 33022 landed in classTags.h; Tcl + python
dispatch, broker, guarded `Domain::commit()` hook; T3/Q4/H8 snapshot over the
ADR-71 kernel; RATE-form advance with the unified reference rule (§12 P0
pin 1); Jacobi-PCG; `-onRemoval keep|drain`; `-subcycle auto` (θ=0.089);
`-staticMode hold|steady` with the loud ⟨A-3⟩ first-march advisory;
`-record`; full send/recvSelf (incl. `dpCommitted_`, payload under its own
`dbTag2_`). Batteries: physics 6/6 (headline anchor: overlay ≡ frozen-toy
`fs1_extrap` to 3.46e-7 under a massless solid; splitting error 1.00× the
frozen `e7_summary` value) + framework 8/8. Opus ×3 adversarial panel PASS
(§11 2026-07-14 entry). Banner + ledgers + manifest row.
**Post-P1 follow-up (PR #577 + reconcile)**: the P1 xfail 6c tripwire led to
the upstream quad/tri element-`rho` serialization bug (garbage restored mass);
fixed in #577, 6c promoted to a hard gate, deformable transient DB restart
supported again. Deferred: P2 iterated driver, P3/P3b explicit lanes, P4
recorder channels + guide.

### P0 — E7 dynamic pins measured (PR #575, merged 2026-07-14)

`adr71_meshless_p_spike/staggered_pins_e7.py` (toy imported as a frozen
library; predictions printed before measurement). Full numbers: the plan's
frozen-constants block + spike `RESULTS.md` §E7. **Three predictions flipped;
the entries below AMEND §3/§4 — where they conflict with the drafted text
above, this log governs:**

1. **fs1 single-pass is NOT "one sweep of the §3.1 iteration."** The
   fixed-POINT L form (RHS `+L·p⁽ᵏ⁻¹⁾`) is O(1)-wrong when run single-pass —
   compliance double-count, effective storage S\*+L+C instead of S\*+C
   (measured: plain 61 % relL2 at every Δt, either L; +resolve classic 41 %;
   E6's "fs1 O(Δt) 0.09 %" was the resolve+oedometric special case — Δu = 0
   lets L impersonate the true Schur compliance, column only). **The v1
   fluid advance is the fixed-stress RATE form**
   `(S* + L + ΔtH)·p₁ = (S* + L)·p₀ + L·Δp_ref − Qᵀ·Δu + Δt·f_seep`,
   with `Δp_ref` = the previous **committed** fluid increment. Measured
   O(Δt) on both L variants and both regimes (late-window order ≥ +1.00;
   L=oed near-deadbeat 5.6e-6). **Unified reference rule (P2-proof):** the
   FIRST fluid advance of a step references the committed Δp (rate form =
   fs1); repeated advances within one step reference the last TRIAL Δp —
   which recovers §3.1 verbatim, so the iterated driver still converges to
   the monolithic BE fixed point with L cancelling. Costs one extra state
   vector (`dpCommitted_`), which must serialize.
2. **`-substep M` is REMOVED from the command surface** (E7.3b refutation:
   error flat 0.2452→0.2534 over M = 1…10 — the coarse-Δt splitting lag
   dominates; fluid sub-resolution buys nothing). `-subcycle N` stays
   (E7.3a: err ~ N^1.2, all N ≤ 50 stable on the L=0 lane; error doubles at
   N·Δt/(h²/c_v) = **θ = 0.089**, pinning `-subcycle auto`).
3. **§3.4 dynamic lane re-pinned.** Fixed-point-L implicit fluid at commit
   is stable to **0.987×** the drained pencil (the ⟨A-7⟩ hedge resolves in
   ZPC's favor) — but carries a Δt-independent O(1) diffusion-rate artifact
   (0.67 relL2 vs exact): stability true, accuracy false. The rate form is
   **CD-unstable at any Δt** (extrapolation = negative damping on
   oscillatory p) — quasi-static/implicit lane only. **P3 default = L = 0
   at Δt ≤ 0.5× the discrete undrained pencil** (consistent + stable;
   the L=0 boundary measured exactly 1.000× the pencil; err 1.5e-3 vs the
   exact reference). The drained-CFL L-lane remains as a documented
   accuracy-degraded opt-in. Explicit fluid (P3b): boundary = **1.32×** the
   pencil on both benchmark soils; the ⟨A-1⟩ material factor
   √(1+K_f/(n·M_oed)) is ~1.85× conservative as a Δt predictor — the Δt_cr
   advisory uses the DISCRETE pencil; the formula is documentation-only.
4. **Removal gates use the fixed-window jump** (window 20·Δt₀ after the
   removal commit; spread 1.47 % across Δt-halving ×4). The single-commit
   jump VANISHES ~Δt^2.6 — no 1/Δt impulse artifact, and no meaningful
   commit-jump invariant to gate on. `-onRemoval drain` kFactor {1,10,100}
   strictly monotone in drainage time ✓.
5. **⟨A-2⟩ confirmed with the footing pair:** classic α²/K_dr — zero
   divergences, zero maxit-hits on column+footing × compressible+
   near-undrained (max 18 iters); oedometric clean too, 3.4× (column) /
   1.6× (footing) fewer iters; 0.5×oed hits maxit (the E6 cliff). Default
   `-fsL classic` frozen; oed opt-in; floor at oed.
