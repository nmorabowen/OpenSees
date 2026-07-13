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
2. **Explicit dynamics is structurally impossible.** The p-rows carry zero
   mass and the pressure equation is parabolic (ADR-71 §6): the existing
   explicit stack (CentralDifferenceLadruno / SMS) cannot advance them, and
   the fork's flagship runs are explicit (18.6 M-element class).

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
   steps at soil coupling strength — measured): auto
   L = α²/(K_dr+4G/3)·M_p from the material initial tangent (3 iters vs 11
   for classic α²/K_dr; 0.5× oedometric diverges — floor guarded).
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
- Auto pin: **L = α²/(K_dr + 4G/3) · M_p** (oedometric variant), K_dr and G
  from the material **initial isotropic elastic tangent** per element —
  exactly the `-stab auto` moduli machinery of ADR-71 §3.3, same
  `updateMaterialStage` dirty-cache treatment. Measured: 3.1–4.6 iters mean
  vs 11.0 for the classic α²/K_dr (Kim–Tchelepi–Juanes), and **0.5× the
  oedometric value diverges** — the optimum sits at a cliff, so the parser
  floors manual `-fsL $scale` at the oedometric value with a loud warning.
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
- **Checkerboard is inherited only at the fixed point**: iterating to
  convergence reproduces the monolithic pair, so the ADR-71 §3.3 α-stab
  (H̃) stays in S\* for equal-order cells in the undrained limit — same
  default, same `-stab auto <α₀>` surface.

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
ZPC-1988 implicit-p/explicit-u operator split; their analysis gives
unconditional stability of the parabolic part with the usual CFL on the
hyperbolic part. **Which wave speed governs the CFL under the split
(drained vs undrained) is MEASURED at P3, not assumed** — the P0 toy
extension emulates the explicit lane first and pins the envelope cheaply.

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
    <-alpha $biotAlpha> <-Ks $Ks> \
    -drained {$nd1 $nd2 ...}                  ;# p-fixed set (≥1 per connected region
                                              ;#   for statics — ADR-71 §3.2 rider)
    <-pInit steady | hydrostatic $gw <$z0> | {$nd $val ...}> \
    <-stab auto <$alpha0> | off | $alpha>     ;# H̃ in S*, ADR-71 §3.3 semantics
    <-fsL auto <$scale>>                      ;# L floor = oedometric; no off
    <-onRemoval keep | drain $kFactor>        ;# fluid life-cycle policy
    <-subcycle $N>                            ;# fluid advance every N commits (explicit lane)
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
| overlay → solid | nodal forces −Q·p (+ f_seep-consistent solid body terms: none — solid gravity stays user-side with mixture ρ, ADR-71 §3.5 convention) | NodalLoads owned by the pattern (H5DRM mechanism) | **constant within a step's Newton loop** (p frozen) — refreshed only by advance() |
| solid → overlay | Δu = uⁿ⁺¹ − uⁿ of region nodes (committed) | direct Node reads | at advance(), post-commit |
| framework → overlay | advance() trigger; removal rescan | Domain::commit hook (ADR-60 trigger precedent; one strictly-additive `// Ladruno` line) — fallback: applyLoad(time) new-step detection if the hook is refused at review | each commit (or every N with -subcycle) |
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
  fully-explicit-both-fields (Xu 2021) → rejected for v1 (lumped-S CFL +
  diffusion limit interacting with ADR-65 machinery; revisit only on
  measured need).
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
| **P0** | toy extension E7 (numpy, no OpenSees build) | exact v1 sequencing pinned: fs1-without-final-resolve error vs monolithic measured over Δt-halving (expect O(Δt); pins the P1 default); explicit-lane emulation (CD solid + implicit fluid): stability envelope measured — drained-vs-undrained CFL governance answered with numbers + subcycle-N degradation curve; removal-step jump consistency (Q mask flips mid-march, no spurious p transient beyond the physical Mandel-Cryer dip); L-floor sweep reproduced on the dynamic path |
| **P1** | `LadrunoPorousOverlay` fs1 implicit end-to-end (Q4/T3/H8 regions) | Terzaghi vs analytic ≤ pinned tol at production Δt; **two-leg gate vs monolithic LadrunoUP** (ADR-71 methodology: mutual Δt-convergence, observed order ≥ 1; the exact-equality leg belongs to P2's iterated driver); **E4-in-OpenSees**: real `remove element` crack, both `-onRemoval` policies, curves vs the toy reference; **sequencing contract test** (step-load column p(0⁺) ≈ q — kills the fluid-first p≡0 trap class); L auto twin-checked vs oracle incl. `updateMaterialStage` dirty path; pattern-timing verification (forces bit-constant across a step's Newton iterations); serial DB round-trip; **full adversarial panel** (per [[feedback_adversarial_gate_when]]: core-touch — Domain hook + pattern semantics — and novel framework integration; the split math itself is spike-validated) |
| **P2** | iterated driver `LadrunoStaggeredAnalyze` (fs-k) | fixed-point gate: iterated staggered ≡ monolithic LadrunoUP same-Δt (E6 leg, target ≤ 1e-6 rel); B4 footing staggered (CB gate under stab); PDMY staged liquefaction column vs the ADR-71 P4 monolithic reference (stage transport → L/stab cache dirty); iteration-count telemetry (mean/max per step) exposed |
| **P3** | explicit lane | CD + overlay on B2/ZS84 column vs implicit monolithic (two-leg); measured CFL envelope vs the P0 pins (criticalTimeStep interplay with ADR-65 machinery documented); `-subcycle` N-sweep gate; energy-balance advisory channel (ADR-69: overlay work terms enter the closure residual — documented, not silently absent) |
| **P4** | ecosystem | overlay p-field recorder channels (own response surface — nodal-DOF recorders can't see it; LadrunoRecorder/Monitor topology rows per [[06_quadrature_global_gp_plan]]); user guide (`LadrunoPorousOverlay_guide.md`) incl. division-of-labor table vs LadrunoUP + init recipes; banner line + ledgers → shipped; apeGmsh emitter runway note (companion repo item) |
| **P5** *(reserved)* | meshless/MPM fluid realization (large-deformation upgrade path); MP/partitioned-domain story | own mini-ADR; the spike's E1–E3 refutations bound what any meshless realization must prove first |

Each phase is one PR off `ladruno`, ledgers updated in-PR.

---

## 8. Risks / open questions

- **Pattern-timing assumption** (the one framework bet): applyLoad is called
  every formUnbalance — the overlay must be a *constant* force source within
  a step. P1 verifies bit-constancy across Newton iterations; if any
  integrator path mutates pattern loads mid-step, the fallback is the
  contact-engine mechanism (backing-element-less FE_Elements contributing
  residual-only forces — ADR-39 precedent, incl. its
  `integrator->formEleTangent` routing quirk).
- **MP absent in v1** — honest and painful: the fork's largest explicit runs
  are MP. A per-partition overlay needs interface flux exchange
  (Schur/Neumann-Neumann class) — genuinely research, parked in P5. Serial
  explicit (the 60 GB emit-host class) is still served.
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

## 11. Implementation log

*(filled as phases land)*
