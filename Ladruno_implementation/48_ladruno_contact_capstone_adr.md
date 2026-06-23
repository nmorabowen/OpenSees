---
title: Ladruno Contact — definitive architecture & unified roadmap (capstone over ADR-39 / ADR-41 / ADR-47)
project: Ladruno
status: draft
priority: high
owner: nmora
tags:
  - implementation
  - contact
  - capstone
  - architecture
  - roadmap
---

# Ladruno Contact — definitive capstone

## What

This is the **capstone** over the three contact ADRs. It is the **single source of truth** for the
fork's contact subsystem: the target end-state, the contracts every formulation must honor, the
authoritative component status, and the **one unified roadmap** that sequences all remaining work.

It **owns no new leaf code** — the per-track ADRs do. It owns the **architecture, the sequencing,
and the cross-cutting decision log**. When status, scope, or sequencing is in question, *this*
document is authoritative; the per-track ADRs hold the detailed design and the phase-by-phase gates.

- **ADR-39** — explicit/robust-first **NTS penalty** `ContactDomain` (the shipped substrate + the
  remaining explicit-stability tiers). Detailed design of record for the NTS lane.
- **ADR-41** — implicit/accuracy-first **mortar + Augmented-Lagrangian (commit-cycle Uzawa)**.
  Detailed design of record for the mortar lane.
- **ADR-47** — the **deferral ledger** (dual/biorthogonal mortar, true-LM saddle-point, self-contact,
  full slide-line smoothing, anisotropic friction, …). *To be created;* this capstone fences its
  scope.

**Committed end-state (dual-lane, decided 2026-06-22):** one `LadrunoContactDomain` engine → broad
phase → **shared header-only kernels** (projection, friction, normal law) → **two formulations**
(`nts`, `mortar`) **× two enforcements** (`penalty`, **commit-cycle ALM**) → **mesh-tying**, all
behind one `-formulation` / `-enforce` selector command. The deferred set is fenced to ADR-47.

## Why

- **Status drift is the concrete failure this prevents.** The 2026-06-22 scoping found ADR-41's
  "ADR-39 shipped P1 only" premise was **7 PRs stale** — ADR-39 had since shipped P1→P3.5. A
  capstone with a single status-of-record table stops three ADRs from disagreeing about reality.
- **A single maintainer needs one backlog, not three.** The remaining work spans ADR-39 (P4–P6),
  ADR-41 (mortar/ALM), and shared-kernel extraction. Sequencing them coherently — with an explicit
  critical path — is the capstone's job.
- **A "definitive" system needs definitive contracts.** Every formulation plugs into the same
  `FE_Element` adapter / stateless-view / Domain-owned-state / null-build-byte-identity rules.
  Writing those once, here, is what lets a new formulation (or ADR-47 upgrade) land safely.

## Where — the target architecture

```
Domain
 └─ LadrunoContactDomain*              engine: owns surfaces, contact defs, broad phase, PATH STATE
     ├─ LadrunoContactSurface          slaveNodes | masterSegments | slaveSegments (facets)
     ├─ LadrunoContactBucketSort        broad phase (spatial grid, 27-neighbour candidates)   [SHIPPED]
     ├─ shared header-only kernels (OpenSees-free, numpy-oracle-tested):
     │    ├─ LadrunoContactProjection   closest-point projection → {ξ̄, gN, n, t1,t2, g[2][2], φ_m}
     │    ├─ LadrunoFrictionKernel      Coulomb+Tresca return map, cone min(μ|tN|+c, τmax), Css/Csl
     │    └─ (normal law)               penalty ⟨−gN⟩ + commit-cycle ALM augmentation
     ├─ per-pair PATH STATE  (Domain-owned, committed/reverted via Domain::commit hooks)
     │    ├─ FrictionState     NTS per-(slave,segment) {gpT, gT0, engaged}                   [SHIPPED]
     │    └─ LadrunoMortarPair mortar per-GP {s_e, λ_N, λ_T, slipFlag, frame}                [PLANNED]
     └─ formulations (narrow phase), selected per contact def:
          ├─ NTS narrow phase   (in LadrunoContactFE::getResidual/getTangent)               [SHIPPED]
          └─ LadrunoMortarSegment  clipped-GP mortar (D, M, weighted gap g̃)                 [PLANNED]

handle()  [LadrunoContactHandler, 33002]  → injects one LadrunoContactFE adapter per pair (runtime tag)
   adapter = STATELESS VIEW; state lives on the Domain; null-build is byte-identical to stock.
enforcement: penalty (both lanes) + commit-cycle ALM (λ updated in LadrunoContactDomain::commit,
   the LadrunoEmbeddedRebar precedent — stock NewtonRaphson, NO custom EquiSolnAlgo for the MVP).
```

### Component inventory & status — SINGLE SOURCE OF TRUTH

| Component | File | Status | Owner ADR |
|---|---|---|---|
| `ContactDomain` engine + lifecycle hooks | `SRC/domain/contact/LadrunoContactDomain.{h,cpp}` | ✅ shipped | 39 (P1) |
| Constraint handler + adapter injection | `SRC/analysis/handler/LadrunoContact{Handler,FE}.{h,cpp}` | ✅ shipped | 39 (P1–P3.5) |
| Surface representation (nodes / segments / **slave segments**) | `SRC/domain/contact/LadrunoContactSurface.{h,cpp}` | ✅ shipped (`SLAVE_SEGMENTS` added C2.0, #373) | 39 + 41 |
| Broad phase (bucket sort) | `SRC/domain/contact/LadrunoContactBucketSort.h` | ✅ shipped | 39 (P2.5) |
| Projection + normal/penalty + gap | `SRC/domain/contact/LadrunoContactKernel.h` | ✅ shipped | 39 (P2b) |
| Coulomb friction return map | `LadrunoContactKernel.h::frictionReturnMap` | ✅ shipped | 39 (P3) |
| Consistent friction tangent (sym + `-consistanttan`) | `LadrunoContactKernel.h::frictionTangentBlock` | ✅ shipped | 39 (P3.5) |
| `-kn auto` penalty sizing | `LadrunoContactKernel.h` / handler | ✅ shipped | 39 (P2b-2b) |
| NTS path state (per-pair friction) | `LadrunoContactDomain::FrictionState` | ✅ shipped | 39 (P3) |
| Shared `LadrunoFrictionKernel.h` (extract + unified `min(μN+c, τmax)` cone + Css/Csl) | `SRC/domain/contact/LadrunoFrictionKernel.h` (relocated from material/nD — contact-specific, include-path) | ✅ shipped (A1, #367) | 48→41 |
| Shared `LadrunoContactProjection.h` (metric `g` + φ_m for mortar GP) | `SRC/domain/contact/LadrunoContactProjection.h` | ✅ shipped (A2, #365) | 48→41 |
| SOFT=1 Courant-stable penalty (explicit) | — | ⏳ pending | 39 (P4) |
| SOFT=2 segment-based penalty (corner/edge) | — | ⏳ pending | 39 (P5) |
| `∂n/∂u` consistent normal tangent + Hertz | — | ⏳ pending | 39 (P2b-2c) |
| Mortar kernel (overlap clip, D/M, g̃) | `SRC/domain/contact/LadrunoMortarKernel.h` | ✅ shipped (C1, #369) | 41 |
| Mortar per-GP state | `SRC/domain/contact/LadrunoMortarPair.{h,cpp}` | 📋 planned | 41 |
| Mortar narrow phase | INLINE in `LadrunoContactFE` MORTAR mode + handler pairing loop (no separate `LadrunoMortarSegment` — mirrors the shipped NTS inline pattern) | ✅ shipped (C2.1, #374) | 41 |
| Commit-cycle ALM + `analyzeAugmented` recipe | `LadrunoContactDomain::MortarNormalState` (λ_N, Uzawa in `commit()`) + `Ladruno_scripts/analyze_augmented.py` + `ladrunoMortarPenetration` query | ✅ shipped (C2.2, #375) | 41 |
| Mesh-tying (`-tie`, zero-gap) | **degenerate-mortar D/M tie (committed route)** — design [[_adr41_c4_design]] (zero-gap limit of contact: active set frozen ON, FULL 3-vector relative-displacement tie, no KKT/cone, ALM `λ_tie` no-clamp; reuses the shipped mortar machinery) | ✅ SHIPPED #381 (+#382 docs, #383 hardening) — `MortarNormalState.{lambdaTie,rtGlobal,isTie}`, `addMortarTieForce` + `epsTie·B̃ᵀB̃⊗I₃` tangent, no-clamp commit Uzawa; `-mortar -tie -epsTie` (⊥ friction); `ladrunoMortarTieResidual` query; oracle `proto_c4_mortar_tie.py` + `test_adr41_mortar_c4_{1,2}` (7); battery 75/75. 3-reviewer adversarial gate PASS | 41 |
| Viscous stabilization (`-visc`) | (normal law + `getDamp`) | 📋 funded option | 41 (Q-VISCOUS) |
| `mu(N,v)` from `frictionModel` wired into contact | — | 🧭 later | 41 (Q-MUDEP) |
| Custom `LadrunoAugmentedNewton` (global Uzawa) | — | 🧭 deferred/trigger-gated | 41 (Q-DRIVER) |
| dual/biorthogonal mortar, true-LM, self-contact, full slide-line smoothing, anisotropic friction | — | 🚫 deferred | 47 |

## How — definitive contracts + the unified roadmap

### Definitive contracts (every formulation MUST honor)

These are the invariants that make the subsystem composable; they are owned here, not per-ADR.

1. **Adapter = stateless view.** `LadrunoContactFE` (and the mortar adapter) compute the narrow
   phase from `Node` trial state in `getResidual`/`getTangent`; **all path-dependent state lives on
   the Domain** (`FrictionState` / `LadrunoMortarPair`) and is committed/reverted via
   `Domain::commit` / `revertToLastCommit`. Adapters are freely destroyed/rebuilt each `handle()`.
   *Mortar-lane requirement (not yet an existing property):* the mortar lane must register its own
   Domain-owned `LadrunoMortarPair` map (keyed per segment-pair GP) **and** extend
   `LadrunoContactDomain::commit()`/`revertToLastCommit()` to iterate it — the shipped hooks
   currently iterate only `theFrictionStates`, so the shared-hook guarantee is a requirement on the
   mortar PR (C2), not an inherited one.
2. **Committed-only multiplier/slip state.** `λ`, elastic slip, and engagement origin are mutated
   **only in `commit`** → a rejected step's `revertToLastCommit` is automatically safe (the
   EmbeddedRebar invariant).
3. **Enforcement = penalty + commit-cycle ALM; no true-LM, no mandatory custom algorithm.** ALM
   augments `λ` once per `Domain::commit` (EmbeddedRebar precedent) on **stock `NewtonRaphson`**;
   within-step augmentation, when a gate proves it needed, is a documented held-load
   `analyzeAugmented` proc (zero-increment re-commits), **not** a bespoke `EquiSolnAlgo`. (Q-DRIVER,
   re-resolved.) Because the sole `LadrunoContactDomain::commit()` chokepoint is `Domain::commit()`
   (which advances `committedTime`, fires recorders, and bumps `commitTag`), the within-step driver
   **MUST NOT re-enter `Domain::commit()`**; it routes through a dedicated augment-only path (or a
   recorder-freeze + `commitTag`-restore protocol) so augmentation passes produce **no spurious
   recorder samples** — this is the D1 gate's "without recorder/load corruption" clause.
4. **`c1` handled once.** The adapter implements the standard `FE_Element` tangent callbacks;
   `getTangent()` delegates to `Integrator::formEleTangent(this)`, so the integrator supplies the
   c-factor as the `fact` argument of `addKtToTang(fact)` (Newmark passes its member `c1`). The
   adapter returns `fact·K_c` — no `getCFactor()` call and no override. Under CDL/explicit,
   `formEleTangent` invokes only `addMtoTang` (a no-op here), so the contact tangent is **identically
   zero** (mass-only LHS).
5. **Conservative-static-superset connectivity over a topology epoch.** Adapter `getID` is a frozen
   superset within an epoch (no `domainChanged()` churn between augmentations / re-solves); active
   set and pairing change only between physical steps.
6. **Null-build byte-identity.** With no contact defined the build is **bitwise-identical to stock**.
   An *active* adapter declares real connectivity, so the numberer permutation differs — the bitwise
   claim holds only for the null case.
7. **Symmetry discipline.** Frictionless / Tresca branches stay **symmetric**; only the Coulomb
   `Csl≠0` branch is **unsymmetric** and opts into `UmfPackGen`/`FullGen` (`-consistanttan`). Default
   tangent is solver-safe symmetric.
8. **Shared kernels are header-only & OpenSees-free**, numpy-oracle-tested build-free (the
   `LadrunoJ2Kernel.h` discipline). One friction return map, two consumers (NTS + mortar).

### The one unified roadmap (sequenced backlog with gates)

Tracks A–E. **Critical path for the new accuracy lane = A → C** (kernels then mortar); **B is
independent** (NTS hardening) and can interleave. The binding kernel edge is **A2 → C1** (the
mortar GP loop needs the metric `g`/`φ_master` the shipped `evalSegment` does not expose); **A1**
(friction extraction) is a soft prerequisite for **C3 only** (mortar may bind the in-place
`LadrunoContactKernel.h` friction per ADR-41 option (a) if A1 slips). The **commit-cycle ALM core
lands with C2**; D1 adds only the *within-step* `analyzeAugmented` held-load refinement.

| Track | Phase | Delivers | Gate | Status |
|---|---|---|---|---|
| **A** shared kernels | A1 | Extract `LadrunoFrictionKernel.h` (in `SRC/domain/contact/`) from `LadrunoContactKernel.h`; add Tresca/`τmax` cap (`min(μN+c, τmax)` cone) + min()-selected `Css`/`Csl`. *Deferred to where consumed:* `-epsT auto` from `γ_crit=0.5%·L` (handler sizing, C2) + the mortar `λ_T`/ALM-aware form (C3). | ADR-39 P3/P3.5 gates green bit-for-bit; oracle FD-checks `Css/Csl` to 1e-6 | ✅ #367 (oracle FD ~1e-12; cohesion/τmax default 0 = bit-for-bit) |
| | A2 | Extract `LadrunoContactProjection.h` (pure fns) returning surface metric `g[2][2]` + `φ_master[4]` for the GP loop | projection vs closed form on tilted + curved facet; bounded-Newton sentinel | ✅ #365 (oracle 7/7, Zone-A green bit-for-bit) |
| **B** NTS lane | B1 | P4 SOFT=1 Courant-stable penalty (explicit `dt_cr` not throttled by `kₙ`) | explicit stability + energy balance | ⏳ |
| | B2 | P5 SOFT=2 segment-based penalty (corner/edge/T-intersection) | corner/edge robustness | ⏳ |
| | B3 | P2b-2c `∂n/∂u` consistent normal tangent + Hertz benchmark | implicit Newton convergence on curved/large-sliding; Hertz `p(r)` | ⏳ |
| **C** mortar lane | C1 | `LadrunoMortarKernel.h`: overlap clip → sub-tri Gauss → `D`, `M`, weighted gap `g̃` (+ linearization stub for C2). **Design/handoff: [[_adr41_c1_design]]** | partition-of-unity `ΣM=ΣD` to 1e-12; **constant-pressure patch test on a non-matched mesh ≤1e-6** | ✅ #369 (oracle 30/30: ΣM=ΣD 1e-12, patch test 9.7e-16, clip-blind 43% biased; C++ standalone bit-for-bit. **Adversarial gate (3 reviewers → PASS/PASS/SALVAGEABLE, no blocker): gates proven structurally robust (reduce to Σφ=1); folded a convexity guard [non-convex facet refused], a back-map convergence flag, and the warped-facet measure-bias doc → C2/ADR-47**) |
| | C2 | Frictionless **commit-cycle ALM** MVP: per-slave-node `λ_N` (Domain-side, Uzawa — NO new DOFs/saddle solver); mortar narrow phase INLINE in a `LadrunoContactFE` MORTAR mode + `SLAVE_SEGMENTS` surface. **Design/handoff: [[_adr41_c2_design]]** (next-up; C1 ships `D,M,g̃`). Phased C2.0 surface → C2.1 penalty adapter → C2.2 Uzawa-on-commit | across-step converged penetration → tol; release→F=0; **Hertz** converges; eqn count constant across augmentations (within-step `maxAug` convergence is gated in D1) | 🔄 in progress — ALM oracle ✅ #372 (18/18); **C2.0** ✅ #373 (`SLAVE_SEGMENTS` + `-mortar` command, inert); **C2.1** ✅ #374 (mortar PENALTY adapter — forces! `MORTAR` LadrunoContactFE mode + handler facet pairing; penetration δ=P/(epsN·A) matched+non-matched, Newton converges, battery 53/53). **C2.2** ✅ #375 (Uzawa-on-commit + per-global-node `λ_N` state — THE headline accuracy win: penetration → an `epsN`-INDEPENDENT tol that penalty cannot reach. `MortarNormalState` on `LadrunoContactDomain` keyed `(contactTag, slaveNodeTag)`, committed-only, augmented once per `Domain::commit()`; held-load `analyze_augmented` proc + `ladrunoMortarPenetration` query. **Crux** [per-node λ but per-facet adapters]: shipped scheme uses LOCAL-gap force + the global multiplier `λ_I` [global gap reserved for the commit Uzawa update + query] — a deviation from the original "lagged-global-gap-in-pressure" recommendation, forced by NewtonRaphson's facet-by-facet residual sweep [see [[LEDGER_quirks]]]. Oracle 28/28 [+T7 crux, +T8 shipped scheme]; `test_adr41_mortar_c2_2` 4/4 [penetration epsN-INDEP, release→F=0 & λ→0, eqn count CONSTANT, ALM on a REAL elastic LadrunoBrick mesh = where the deferred Hertz gate lands]; battery 58/58). C2 COMPLETE. The curved-surface **Hertz `p(r)` refinement study** remains a deeper validation item [D1/validation phase] |
| | C3 | Frictional mortar: adopt shared `LadrunoFrictionKernel` (`λ_T`); Coulomb unsymmetric branch | incline `a=g(sinθ−μcosθ)`; Tresca cap; `Csl` FD-checked; Δt-independent converged answer | 🔄 in progress — **C3.0 design + oracle ✅ #376**: design/handoff [[_adr41_c3_design]] (resolves the architecture: per-GLOBAL-slave-node tangential slip `gpT[3]` + multiplier `λ_T[3]` keyed `(contactTag,slaveNodeTag)` mirroring the shipped C2.2 `λ_N` — NOT the ADR's never-built per-GP `LadrunoMortarPair`; reuse the shipped `LadrunoFrictionKernel` verbatim; phasing C3.1 force → C3.2 tangent → C3.3 optional `λ_T` Uzawa). Oracle `proto_c3_mortar_friction.py` 20/20: return map stick/slip + `K_ss` vs FD both cone sides, per-node D/M assembly self-equilibrium + consistent constant-slip load, T-LOCAL determinism (the C2.2 lesson re-shown for `λ_T`), incline sign gate `a=g(sinθ−μcosθ)`, unified cone (Coulomb/Tresca/cohesion + `Csl`=0 on the τmax branch). **C3.1 mortar friction FORCE ✅ #377** (penalty `λ_T≡0`): `MortarNormalState` gains `gpT/gpTtrial/gT0/engaged` (folded — same key as `λ_N`, committed in `Domain::commit`); `LadrunoContactFE::addMortarFriction` reuses the shipped `frictionReturnMap` per in-contact node with `N_I=−p_normal,I`, scatters `tFric` via D/−M; `-mu/-epsT/-cohesion/-tauMax` command surface. **CRUX FIX:** the tangential slip is built from DISPLACEMENTS not positions (the closest-point projection makes the weighted relative position purely normal — the NTS SEGMENT lesson; see [[LEDGER_quirks]]). Explicit (CDL) tests `test_adr41_mortar_c3_1` 3/3 (Tresca driven `a=(Q−cap·A)/m` N-independent, Coulomb sign-opposes-motion, μ=0 byte-identical); battery 61/61. **C3.2 friction TANGENT ✅ #378** (the SYMMETRIC consistent tangent — implicit frictional Newton converges, singular without it): `addMortarTang` scatters the per-node 3×3 `K_ss` (`frictionTangentBlock`) via `B̃=[D,−M]` — `K[A][B]+=Σ_I(b_IA b_IB/a_I)K_ss_I` — at the SAME slip/N/active-set the residual used. Default symmetric (drop the non-sym Coulomb `Csl` — solver-safe; the full mortar `Csl` couples through the normal-gap operator, deferred with the geometric terms). MAJOR-2 (C3.1 gate) fixed: `gT0/engaged` double-buffered + reverted (now reachable under implicit). Oracle T6 FD-checks the assembled tangent (rel ~1e-10); `test_adr41_mortar_c3_2` 4/4 (static stick + slip[spring-restrained] CONVERGE, symmetric-solver-safe ProfileSPD==FullGeneral, Coulomb static); battery 65/65. **C3.3 ✅ #379 — C3 COMPLETE** (tangential `λ_T` Uzawa + the non-sym mortar `Csl`): `MortarNormalState` gains `lambdaT/lambdaTtrial`; the residual+tangent inject `λ_T` via the OFFSET TRICK `gTeff_eff = gTeff + λ_T/epsT` (reuse the penalty return map), one Uzawa step per commit `λ_T ← −tFric` ⇒ the spurious STICK elastic creep → 0 at FINITE epsT (epsT-INDEPENDENT tangential position — the λ_N headline, tangential side); `-consistanttan` now FUNCTIONAL (the non-sym Coulomb `Csl` scatters through the same `b·b/a` with kn=epsN — `consistent && !initialStiff`). Oracle 30/30 (+T7 consistent Csl FD ~1e-11, +T8 λ_T stick-creep→0 epsT-indep); `test_adr41_mortar_c3_3` 2/2 (held-load augmentation removes the creep epsT-indep; `-consistanttan` Coulomb slip → same root as symmetric); battery 68/68. **C3 COMPLETE** (penalty + ALM friction, symmetric + consistent tangent, explicit + implicit). **C4 ✅ #381 — see the C4 row; the ADR-41 mortar lane is COMPLETE (C1→C4).** |
| | C4 | Mesh-tying (`-tie`, zero-gap = active set frozen ON) — degenerate-mortar D/M tie. **This is the committed route;** ADR-39 P6's `MP_Constraint`/ADR-30-projection tie is a *separate alternative implementation* of the same `-tie` capability (different kernel, different gate), **superseded by C4 unless the explicit-lane momentum-conserving tie is explicitly scoped in** | uniform/linear traction across non-matched interface = single-block stress (exact patch) | ✅ SHIPPED #381 (+#382 docs, +#383 review-hardening). The FULL 3-vec relative-displacement tie `r_I=ΣD u_s−ΣM u_m → 0` (no KKT/clamp/cone): `t_I=λ_tie,I+epsTie·(r_I/a_I)` via D/−M, SPD tangent `epsTie·B̃ᵀB̃⊗I₃`, no-clamp ALM `λ_tie ← λ_tie+epsTie·(rtGlobal/aGlobal)` once per commit ⇒ `‖r‖→0` epsTie-INDEP (== exact Lagrange). **Shared-node pre-req RESOLVED** (LEDGER_quirks MAJOR-1): `r` is a LINEAR accumulation ⇒ order-independent global accumulator `λ_N` uses (NOT the friction last-writer-wins, which stays fenced). `-mortar -tie [-epsTie auto\|val]` ⊥ `-mu/-cohesion/-tauMax`; `ladrunoMortarTieResidual` query; `analyze_augmented(query=…)` hook. NO new kernel/DOF/classTag; -tie-absent byte-identical. Oracle `proto_c4_mortar_tie.py` (PATCH TEST exact, penalty O(1/epsTie), ALM==exact, tangent FD sym+PSD, shared-node order-indep); `test_adr41_mortar_c4_{1,2}` 7/7; battery 68→75/75. **3-reviewer adversarial gate → PASS** (no correctness bug; 1 defensive NIT hardened in #383). **ADR-41 MORTAR LANE COMPLETE (C1→C4).** |
| **D** enforce/robust | D1 | Commit-cycle ALM mechanics + `analyzeAugmented` Tcl/Py proc (held-load re-commit) | within-step `‖g̃‖→augTol` without recorder/load corruption | 📋 |
| | D2 | Viscous stabilization `-visc μ_c` (`p_visc=μ_c·v_rel`) — funded option for pounding/rocking chatter | chatter/snap-through damped; off under any arc-length lane | ✅ **D2 COMPLETE — D2.1 #385 (NTS: RIGID_PLANE + SEGMENT) + D2.2 #TBD (MORTAR contact)** — design [[_adr41_d2_design]]. The velocity-proportional normal damper reuses the penalty operator driven by velocity: `ġ=B·v`, viscous force `f_visc=−μ_c·ġ·B` (in `getResidual` from `getTrialVel`) + the SPD damping tangent `C_visc=μ_c·B Bᵀ` (in `addCtoTang`). Active only in contact (`gap<0` / mortar `p_I<0`); **statics identically inert** (`v≡0` + `StaticIntegrator` never calls `addCtoTang` ⇒ byte-identical to `μ_c=0`). Force-only under mass-only CDL (explicit), force + consistent C tangent under Newmark/HHT (implicit). Command `contact … -visc <μ_c>` / `contactPlane … -visc <μ_c>` (off by default; **REFUSED on `-tie`** — a bond has no contact-chatter regime; works on NTS AND mortar CONTACT). It is a Kelvin–Voigt contact ⇒ exact restitution `e=exp(−ζπ/√(1−ζ²))`, `ζ=μ_c/(2√(kₙ·m))`. **D2.2** added the MORTAR mode: per in-contact node `p_visc_I=μ_c·ḡ̇_I` (`ḡ̇_I=n·(ΣD v_s−ΣM v_m)/a_I`) scattered via D/−M·n like the normal force + the `μ_c·B̃ᵀ diag(W) B̃ ⊗ n⊗n` damping tangent in addCtoTang — the C2 `addMortarTang` block with `epsN→μ_c`. Oracle `proto_d2_viscous.py` (operator FD sym+PSD, static-inert, restitution vs analytic, energy balance 0.005%, active mask; +T6 mortar B̃ operator FD sym+PSD+self-equilibrium). `test_adr41_viscous_d2` 8/8 (explicit-impact Kelvin–Voigt through real CDL, static byte-identical [NTS+mortar], implicit Newmark converges+damps [the C tangent, NTS+mortar], SEGMENT + MORTAR viscous, `-visc`⊥`-tie` refused / `-mortar -visc` accepted); battery 75→83/83. **3 adversarial reviewers total → PASS** (no bug). **D2 COMPLETE (NTS + mortar viscous stabilization).** |
| **E** deferred | — | dual/biorthogonal mortar, true-LM saddle-point + inf-sup, self-contact, full slide-line + nodal-normal smoothing, anisotropic friction, `mu(N,v)` wiring, custom `LadrunoAugmentedNewton` | — | 🚫 ADR-47 |

### Unified command surface (the definitive API)

```tcl
constraints LadrunoContact
contactSurface 1 -kind masterSegments -faces $master
contactSurface 2 -kind slaveSegments  -faces $slave    ;# slaveNodes for NTS-only

# NTS penalty (ADR-39, shipped)
contact 10 -master 1 -slave 2 -formulation nts    -kn auto -kt auto -mu 0.3
# Mortar + commit-cycle ALM (ADR-41)
contact 11 -master 1 -slave 2 -formulation mortar -enforce alm -epsN auto -augTol 1e-8 -maxAug 20 -ngp 2
# Mortar + ALM, frictional Coulomb (unsymmetric)
contact 12 -master 1 -slave 2 -formulation mortar -enforce alm -friction coulomb -mu 0.3 \
        -kt auto -epsT auto -augTol 1e-8 -maxAug 20      ;# pure Coulomb (no Tresca cap); add -tauMax <val> to cap
system UmfPackGen
# Mesh-tying (degenerate mortar, active set frozen ON)
contact 13 -master 1 -slave 2 -formulation mortar -enforce penalty -tie

algorithm Newton                          ;# MVP: λ augments per Domain::commit (EmbeddedRebar pattern)
# within-step augmentation, if a gate needs it:
#   analyzeAugmented $augTol $maxAug       ;# held-load zero-increment re-commits — NOT a custom algorithm
```

## Decision log (resolved architecture questions — pointers to the owning ADR)

| Decision | Resolution | Owner |
|---|---|---|
| Enforcement family | Penalty + **commit-cycle ALM** (Uzawa-over-penalty, zero new DOFs) | 41 Q-AL/Q-DOF |
| Augmentation driver | **Commit-cycle primary** (stock Newton); custom `EquiSolnAlgo` deferred behind a within-step trigger | 41 Q-DRIVER |
| True-LM / saddle-point | **Rejected** (zero-diagonal + active-set re-number); → ADR-47 | 41 Q-DOF |
| Mortar basis | **Standard + overlap clipping** (passes non-matched patch test); dual/biorthogonal → ADR-47 | 41 Q-MORTARLITE |
| Friction cone | **Unified `min(μ|tN|+c, τmax)`** (Coulomb/Tresca/cap one return map); metric-aware `R=g·r` | 41 |
| Tangential penalty sizing | `-epsT auto` from `γ_crit = 0.5%·L_facet` (Abaqus elastic-slip bound) | 41 |
| Friction code ownership | ADR-39 shipped the kernel first → ADR-41 **extracts/shares** it (direction reversed) | 48→41 |
| Solver symmetry | Coulomb branch unsymmetric only; frictionless/Tresca symmetric | 41 |
| Explicit ALM | **Not promised** — full ALM is implicit-only; explicit is single-pass penalty | 41 Q-EXPLICIT |
| Surface-normal smoothing | Faceted per-GP normal in v1; nodal-normal averaging → ADR-47 | 41 Q-NORMAL |

**Open / cross-cutting questions this capstone tracks:** mortar epoch cost vs sliding distance
(Q-EPOCH/Q-GRAN), implicit SOE fill growth (Q-IMPLFILL), MP-constraint composition restriction
(Q-CONSTR), small- vs finite-sliding cost lane (Q-SLIDING). **RESOLVED: D2 viscous stabilization
earns a committed phase** — D2.1 (NTS) shipped #385 as a gated, off-by-default `-visc` option
(the Q-VISCOUS "best ROI of the residue" recommendation, funded); D2.2 (mortar viscous) follows on demand.

## Relationship to ADR-39 / ADR-41 / ADR-47

- **This capstone (48)** owns: target architecture, contracts, the status-of-record table, the
  unified roadmap, the decision log, and the deferred-set fence. Update it **in the same PR** as any
  change that moves a component's status or re-sequences the backlog.
- **ADR-39** owns: the NTS-penalty detailed design and its per-phase gates (P1–P6). Remains the
  implementation record for the explicit-first lane.
- **ADR-41** owns: the mortar/ALM detailed design, the friction-kernel extraction spec, and its
  per-phase gates. Remains the implementation record for the accuracy-first lane.
- **ADR-47** (created — [[47_ladruno_contact_deferrals_adr]]) owns: the deferred-features ledger
  (rejection rationale + re-open trigger per item); each item graduates to its own ADR when re-opened.

## Risks / open questions

- **Capstone drift (meta-risk).** The capstone is only useful if kept current. Mitigation: the
  status-of-record table is the *only* place phase status lives; per-track ADRs link here rather than
  re-asserting global status. Treat a stale capstone row as a build-control defect (same class as a
  stale ledger row — which this scoping already caught once).
- **Single-maintainer load.** Dual-lane is the committed end-state, but the **critical path is A→C**;
  B (NTS hardening) and D2 (viscous) are independent and can be deprioritized without blocking the
  mortar lane. The roadmap is explicitly parallelizable so scope can flex without re-architecting.
- **ADR-47 deferral ledger — created** ([[47_ladruno_contact_deferrals_adr]], 2026-06-23). The
  deferred set now has a documented home with a re-open trigger per item (the C1 patch-test and
  C2 sliding/Hertz gates lean on the non-dual-basis inf-sup and faceted-normal-chatter rationale).
  *Keep it current:* when a deferred item is re-opened it graduates to its own ADR and the ledger
  row is marked promoted.

## Implementation log

*(filled as phases land; each phase updates its row in the status-of-record table above and links
its PR. Move detailed notes to the owning per-track ADR / `Ladruno_internal/`.)*

## Adversarial review log

**Gate 1 — 2026-06-22 (5-lens workflow, 28 agents, refute-by-default verification). Verdict:
SALVAGEABLE → all dispositions folded; now PASS-equivalent.** Five independent lenses
(Architecture/contracts · Roadmap/sequencing · Status-accuracy-vs-tree · Scope/single-maintainer
realism · Cross-ADR consistency) raised **22 findings → 17 confirmed** after each was adversarially
verified against the ADR text *and* the live tree. **0 blockers, 0 architecture-breaking defects** —
the core structure (Domain-owned state, stateless adapters, commit-cycle ALM on stock Newton,
header-only kernels, A→C critical path) held. Two `fix-now` issues, both with clean verbatim fixes:

| ID | Sev | Issue | Fix (folded) |
|---|---|---|---|
| **ARCH-1** | major | Contract #4 misdescribed the `c1` API — named `getCFactor()` (never called; returns 0) + a non-existent `addKtToTang` override | Rewrote to the true `formEleTangent(this) → addKtToTang(fact)` mechanism (verified `LadrunoContactFE.cpp:244-290`, `Newmark.cpp:292`) |
| **SCOPE-3** | major | Mesh-tying left **double-owned** (degenerate-mortar vs ADR-39 P6 MP path) — the exact cross-ADR drift this capstone exists to prevent | Named the degenerate-mortar D/M tie the committed route (C4); ADR-39 P6 reframed as a separate superseded alternative |

`fold` (applied): **ROAD-3** (D1/C2 ordering — commit-cycle core lands with C2, D1 = within-step
refinement), **XADR-2** (`-tauMax 0` example = frictionless; dropped), **SCOPE-4** (ADR-47 stub
due **before C1**), **ARCH-3** (within-step driver must not re-enter `Domain::commit` — recorder/
`commitTag` corruption), **ARCH-4** (mortar state must extend the commit hooks — not inherited),
**XADR-1/XADR-4** (stale friction-direction narrative in ADR-41 marked superseded), **XADR-5**
(ADR-39 Q-AL now names the `LadrunoEmbeddedRebar::commitState` precedent). `reject` (over-literal /
dups / correct-as-written): ROAD-1, ROAD-2, ROAD-4, STATUS-2, XADR-3, XADR-6. `defer` (separate
source-hygiene PR): **STATUS-1** — stale "zero adapter" comment in `LadrunoContactFE.cpp:22`.
