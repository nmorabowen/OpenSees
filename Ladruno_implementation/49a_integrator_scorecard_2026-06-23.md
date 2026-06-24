---
title: Integrator Scorecard — run 2026-06-23 (assessment + crosswalk + gap analysis)
project: Ladruno
status: complete
owner: nmora
tags:
  - integrator
  - study-result
  - gap-analysis
aliases:
  - ADR-49a
  - Integrator Scorecard
updated: 2026-06-23
parent: "[[49_ladruno_integrator_study_workflow_adr]]"
---

# ADR-49a — Integrator scorecard (run of 2026-06-23)

A run of the [[49_ladruno_integrator_study_workflow_adr]] workflow. Phase 2
(assessment) was done by four source-reading agents over `SRC/analysis/
integrator/*.cpp`; Phase 3 (crosswalk) used the `abaqus-theory-procedures` and
`kratos` skills; Phase 4 (gap analysis) is synthesized below.

> **Fidelity:** Abaqus + Kratos legs are skill-backed (high). **LS-DYNA has no
> skill** — its column is knowledge-based and should be source-confirmed before
> any LS-DYNA claim is published externally.

**Census:** ~38 upstream + ~12 fork "true" integrators + ~24-variant
hybrid-sim cluster + 7 eigen/mass helpers. Base hierarchy unchanged:
`Integrator → IncrementalIntegrator → {Static,Transient}Integrator`;
`EigenIntegrator` direct from `Integrator`.

---

## Part 1 — Assessment by family (condensed)

Columns: **Scheme/order · Stability · Maturity · Key limitation · Verdict.**
Full 10-field rubric tables are in the run transcript; condensed here.

### 1.1 Static path-following (`StaticIntegrator`)

| Integrator | Scheme/order | Stability/control | Maturity | Key limitation | Verdict |
|---|---|---|---|---|---|
| `LoadControl` | proportional λ; iter-count step scaling | crude adaptive (Jd/Jlast) | shipped, mature | heuristic step adapt; no error control | **keep** (baseline) |
| `LoadPath` | follows prescribed λ vector | none (kinematic) | legacy/thin | no equilibrium feedback | document / maybe deprecate |
| `StagedLoadControl` | LoadControl + inactive-DOF diag stabilization | hard-coded 1.0 penalty; MPI-aware | shipped, mature | penalty not tunable; no momentum restore | **keep** (staged) |
| `DisplacementControl` | single-DOF disp control, Ramm-scaled | weak to snap-back | shipped, mature; has sensitivity | single-DOF; sensitivity bloats core | **keep**; refactor sensitivity out |
| `DistributedDisplacementControl` | DispControl + MPI plumbing | **stub** — channels allocated, no collectives | shipped but **inert** | MPI send/recv never called | **deprecate or complete** |
| `ArcLength` | spherical quadratic dU²+α²dλ²=s² | robust past snap-through if α≤1 | shipped, mature; heavy DDM | no Ramm adapt; quadratic can have no real root | **keep** (workhorse) |
| `ArcLength1` | arc-length, **linear** corrector | less robust past limit points | shipped, **legacy** | truncation error; near-dup of ArcLength | **consolidate** into ArcLength `-linear` |
| `MinUnbalDispNorm` | Ramm/Chan min-unbalance-norm | norm-based adapt (robust) | shipped, mature | DDM embedded | **keep** |
| `HSConstraint` | hyperspherical scaled arc-length | user ψ_u/ψ_f weights | shipped, research | no default tuning rule; opaque u_ref | keep if tested; document |
| `EQPath` | weighted equal-path (Rezaiee-Pajand) | research method | shipped, **research** | opaque `type`/`m` params, no validation | **mark experimental / validate or drop** |
| `HarmonicSteadyState` | LoadControl + sinusoidal phase | quasi-static (slow assumed) | shipped, niche | no check period ≫ structural | keep if validated |
| `LadrunoArcLength` (33004) | ArcLength + Ramm `-adapt` + viscous `-stabilize` + snapshot/revert | superset of ArcLength | **fork, shipped** | `-stabilize` mutually exclusive w/ quadratic; calib empirical | **keep** (production superset) |
| `LadrunoIndirectControl` (33006) | weighted multi-DOF (CMOD) control | monotone by design | **fork, shipped** | fixed MAXC=64 array; silent fail if ctrl DOF constrained | **keep**; dynamic-alloc fix |

### 1.2 Implicit transient (`TransientIntegrator`)

| Integrator | Scheme/order | Stability | Maturity | Key limitation | Verdict |
|---|---|---|---|---|---|
| `Newmark` | Newmark β,γ; 2nd | cond.; uncond. if γ=½,β≥¼ | mature; **only one with sensitivity** | sensitivity only disp-form | **keep** (core) |
| `Newmark1` | Newmark + Rayleigh at newStep | same | mature | no sensitivity; near-dup | **consolidate** into Newmark |
| `StagedNewmark` | Newmark + inactive-DOF diag fill | inherits | mature | overhead; staged-only | **keep** (staged) |
| `HHT` | HHT-α; 2nd | uncond. α∈[0,⅓]; L-stable α=0 | mature | fixed β(α),γ(α); no sensitivity | **keep** |
| `HHTGeneralized` | indep αI,αF; 2nd | uncond. | mature | overlaps GeneralizedAlpha | keep, but see consolidation |
| `GeneralizedAlpha` | Chung-Hulbert αM,αF; 2nd | uncond. | mature | no sensitivity | **keep** (gold standard) |
| `Collocation` | collocation θ; β poly-fit | stable θ∈[½,1] | legacy | ad-hoc β fit | keep for completeness |
| `WilsonTheta` | Wilson-θ; 2nd | stable θ≥1.37 | legacy | needs θ>1; outdated | keep (historical) |
| `Houbolt` | 3-step Houbolt; 2nd | A-stable (after 2 trap steps) | legacy | brittle step-count logic | keep/low-priority |
| `ParkLMS3` | 3-step Park LMS; 2nd | A-stable | legacy | velocity history, brittle | keep/low-priority |
| `TRBDF2` | composite Trap+BDF2; 2nd | A-stable (Bathe-ref) | mature | accuracy dip first 2 steps | **keep** (energy-conserving) |
| `TRBDF3` | Trap→BDF2→Houbolt; 2nd | A-stable | experimental | 3-cycle brittle; uses Houbolt not Bathe-3 | keep/experimental |
| `BackwardEuler` | BE 1st (2 variants) | A-stable, uncond. | mature | 1st-order accuracy | **keep** (robust fallback) |
| `PFEMIntegrator` | Newmark-like for particle-FEM | cond. (Newmark) | mature (PFEM) | PFEM-domain-specific | **keep** (PFEM) |
| `GimmeMCK` | **utility** — exports M/C/K | n/a | utility | **not a time integrator** | keep; reclassify/doc |

### 1.3 Explicit transient (`TransientIntegrator`)

| Integrator | Scheme/order | Stability | Maturity | Key limitation | Verdict |
|---|---|---|---|---|---|
| `CentralDifference` | leap-frog; 2nd | cond. CFL=2/ω; damping halves Δt_cr | upstream baseline | "assume U_{-1}=U_0" start; Rayleigh on LHS | keep (baseline) |
| `CentralDifferenceAlternative` | leap-frog, M-only | cond. | upstream | discards damping; redundant | **deprecate** |
| `CentralDifferenceNoDamping` | leap-frog, M-only | cond. | upstream | damping-via-residual impractical | **deprecate** (pedagogical) |
| `ExplicitDifference` | leap-frog variant | cond. | upstream, **buggy** | `updateCount>2` malformed (line 274) | **deprecate** |
| `ExplicitDifferenceStatic` | leap-frog + FLAC local damping | cond. | fork-ish (2024) | hard-coded FLAC consts; needs Linear solver | superseded by EB-LNVD? assess |
| `NewmarkExplicit` | explicit Newmark γ | cond. (γ=½→CD) | upstream, mature | — | **keep** |
| `KRAlphaExplicit` | Kolay-Ricles α; 2nd | **unconditionally stable** | upstream, mature | per-step dense eigensolve cost | **keep** (premium stability) |
| `CentralDifferenceLadruno` (33003) | clean robust CD; 2nd | cond.; `-cflAbort`/`-recompute`; constraint projector | **fork, shipped** | none critical | **keep** (recommended baseline) |
| `CentralDifferenceSMS` | CD + selective mass scaling (lumped) | scaled-to-target; rejects `-cflAbort` | **fork, shipped** | MP-slaves unscaled; lumped only | keep (see consolidation) |
| `CentralDifferenceSMSConsistent` | CD + Olovsson consistent SMS (PCG) | scaled; needs Diagonal/MPIDiagonal SOE | **fork, shipped** | strict SOE req; PCG tuning | keep (see consolidation) |
| `ExplicitBathe` (33000) | Noh-Bathe 2-substep; 2nd | cond., **~2× CFL bonus** | **fork, shipped** | 2 solves/step | **keep** (preferred explicit) |
| `ExplicitBatheLNVD` (33002) | Bathe + FLAC local-no-velocity damping | cond. | **fork, shipped** | — | **keep** (explicit quasi-static) |
| `ExplicitBatheSMS` / `…SMSConsistent` | Bathe + SMS (lumped/consistent) | scaled-to-target | fork, shipped | SMS v1 limits | keep (see consolidation) |
| `ExplicitBatheLNVDSMS` / `…SMSConsistent` | Bathe + LNVD + SMS (combinatorial) | scaled | fork, shipped | **cross-layer interaction untested** | **consolidate** (flag explosion) |
| `LadrunoDynamicRelaxation` (33005) | leap-frog + Gershgorin fictitious mass + Cundall kinetic damping | **no physical CFL** (pseudo-dt) | **fork, shipped** | matrix-free; quasi-static only | **keep** (quasi-static) |

### 1.4 Hybrid-simulation cluster (cross-cuts 1.2/1.3) — OpenSees-UNIQUE

~24 variants (Schellenberg / PEER, real-time hybrid testing): `Newmark`/`HHT`/
`Collocation` × {`HSFixedNumIter`, `HSIncrLimit`, `HSIncrReduct`} × {base, `_TP`},
plus `AlphaOS`(`Generalized`)(`_TP`), `HHTExplicit`, `HHTGeneralizedExplicit`,
`KRAlphaExplicit_TP`.

| Sub-family | What it does |
|---|---|
| `_TP` ("table passing") | trapezoidal variant assembling full unbalance (`formUnbalance`/`Put`) for synchronized actuator/controller comms |
| `HSFixedNumIter` | Lagrange-interpolate past disps to hit a fixed iteration count (rate-limited actuators) |
| `HSIncrLimit` | norm-cap accumulated increments (suppress spurious load/unload) |
| `HSIncrReduct` | constant reduction factor on increments (simplest, least smooth) |
| `AlphaOS` | α operator-splitting to decouple numerical vs physical substructures |

**Verdict: score as ONE cluster. No commercial-code analogue** — Abaqus/Kratos/
LS-DYNA do not do real-time hybrid testing / experimental substructuring. This is
an **OpenSees strength**, not a gap. Keep all; no action.

### 1.5 Eigen + mass/utility helpers (enablers, not standalone integrators)

| Item | Role | Standalone/enabler | Maturity |
|---|---|---|---|
| `EigenIntegrator` | K,M assembly for eigenproblem (→ ARPACK/LAPACK) | standalone (modal) | Berkeley core, mature |
| `CriticalTimeStep` | per-element Δt_cr eigenproblem; lumping (RowSum/Diag/HRZ); MPI_MIN | enabler (shared kernel) | fork 2026 |
| `LadrunoFictitiousMass` | Gershgorin/lumped/unity artificial M* | enabler (DR, ArcLength) | fork 2026 |
| `LadrunoMassLumping` | HRZ mass-conserving diagonal lump | enabler | fork 2026 |
| `LadrunoMassScaling` | selective SMS (lumped DT2MS + consistent Olovsson) | enabler | fork 2026 |
| `LadrunoMassScalingEnergy` | registry: integrator publishes M̄ blocks → energy recorder | enabler (glue) | fork 2026 |
| `LadrunoConsistentRefine` | matrix-free PCG for M̃ a = r (serial + MPI) | enabler | fork 2026 |

---

## Part 2 — Cross-code crosswalk (Phase 3)

| OpenSees | Abaqus (skill) | Kratos (skill) | LS-DYNA (knowledge — confirm) | They have / we don't |
|---|---|---|---|---|
| `Newmark` | Newmark (β,γ tied to HHT α) | `newmark` scheme option | implicit Newmark (`*CONTROL_IMPLICIT_DYNAMICS`) | — |
| `HHT`, `HHTGeneralized`, `GeneralizedAlpha` | **HHT-α (default, α≈−0.05)** + generalized | **Bossak (default)** + generalized-α | HHT-α implicit | parity (3 codes converge here) |
| `CentralDifference*`, `*Ladruno`, `ExplicitBathe*` | Abaqus/Explicit central-diff | `explicit_central_differences_scheme` | **core explicit solver** | **bulk viscosity** (b1/b2), **hourglass control**, contact-aware auto Δt |
| `KRAlphaExplicit` | — (no uncond. explicit) | — | — | **we have, they don't** |
| `CentralDifferenceSMS*`, `ExplicitBatheSMS*`, `LadrunoMassScaling` | mass scaling (DT2MS-style) | (explicit solver mass-scaling opts) | DT2MS / selective mass scaling | LS-DYNA mature; **our consistent/Olovsson exceeds DT2MS** |
| `ArcLength`, `LadrunoArcLength`, `MinUnbalDispNorm`, `HSConstraint`, `EQPath` | **modified Riks** | `arc_length` strategy | `*CONTROL_IMPLICIT_AUTO` arc-length | Riks predictor uses normal-plane / per-DOF v_max (we use spherical quadratic) |
| `LadrunoArcLength -stabilize` | `*STATIC, STABILIZE` (energy-fraction) | — | — | Abaqus STABILIZE is a 1st-class procedure |
| `LadrunoIndirectControl` | — (Abaqus uses Riks) | (displacement-control condition) | — | **we have, they don't** (CMOD) |
| `LadrunoDynamicRelaxation` | — (uses STABILIZE) | `formfinding_strategy` (DR-like) | **dynamic relaxation** (`*CONTROL_DYNAMIC_RELAXATION`) | LS-DYNA has 1st-class DR |
| `TRBDF2/3` | — (HHT covers) | — | — | we have composite; **but not modern implicit Bathe β1/β2** |
| `Collocation`,`Wilson`,`Houbolt`,`ParkLMS3` | — (legacy equivalents) | — | — | legacy, low value |
| hybrid-sim cluster (§1.4) | — | — | — | **OpenSees-unique strength** |
| `EigenIntegrator` | Lanczos/subspace (FREQUENCY) | `eigen_value` solver | implicit eigenvalue | parity |
| **— none —** | **half-increment-residual auto-dt** | (time-stepping process) | **`*CONTROL_IMPLICIT_AUTO`** | **adaptive/automatic time stepping — GAP** |
| **— none —** | **bulk viscosity (b1,b2)** | — | **bulk viscosity (std)** | **explicit bulk viscosity — GAP** |

---

## Part 3 — Key cross-cutting findings

1. **Sensitivity (DDM) is almost absent.** Only `Newmark` (disp-form),
   `DisplacementControl`, `ArcLength`, `MinUnbalDispNorm` implement
   `formSensitivityRHS`/`commitSensitivity`. `HHT`, `GeneralizedAlpha`, and the
   entire explicit/fork family have **no sensitivity** → blocks reliability /
   fragility / FORM workflows on the fork's flagship explicit integrators.
2. **Where it exists, sensitivity bloats the core class** (50+ LOC of DDM state
   vectors inline). Refactor to an optional mixin/subclass.
3. **`DistributedDisplacementControl` is an inert stub** — MPI channels
   allocated, no collective calls. Either finish or remove.
4. **Upstream CentralDifference triplet is cruft.** `CentralDifferenceAlternative`
   + `CentralDifferenceNoDamping` are subsumed by `CentralDifference`;
   `ExplicitDifference` has a malformed `updateCount>2` guard (line 274).
5. **ExplicitBathe variant explosion** — 6 classTags from stacked subclassing
   (`EB ← LNVD`, `EB ← SMS`, `← LNVDSMS`, `← Consistent`). Cross-layer
   combinations (esp. LNVD × Consistent-PCG) are **not independently validated**.
   Collapse to one class + orthogonal flags `-lnvd -sms -consistent`.
6. **SMS ⊥ Δt_cr query** (MF-1): all SMS classes hard-reject `-cflAbort`/
   `-recompute` because the per-element eigensolve can't see nodal mass
   augmentation. Acceptable but should report "Δt_cr = pre-scaling estimate".
7. **Parallel-safety dichotomy in consistent SMS** needs audit:
   `CentralDifferenceSMSConsistent` flags "not parallel-safe" yet
   `ExplicitBatheSMSConsistent` claims distributed-PCG support. Reconcile.
8. **Near-duplicates to consolidate:** `ArcLength`/`ArcLength1`,
   `Newmark`/`Newmark1`, `HHTGeneralized`/`GeneralizedAlpha` (overlap).
9. **Research-grade, under-validated:** `EQPath`, `HSConstraint` — opaque
   parameters, no reference tests. Gate as experimental or validate.

---

## Part 4 — Gap analysis → ranked new-integrator candidates

Ranked by (fit-to-fork-goals × demand) ÷ build-cost. Each, if green-lit,
becomes its own design ADR (52+).

### Rank 1 — Automatic / adaptive time stepping (implicit) — **BUILD**
- **Gap:** Both Abaqus (half-increment residual) and LS-DYNA
  (`*CONTROL_IMPLICIT_AUTO`) auto-control Δt on error/convergence. OpenSees has
  **nothing built-in** — users script `analyze` retries in Tcl/Python.
- **Fit:** directly strengthens the robust-solve driver ([[31_ladruno_robust_solve_driver_adr]]);
  cuts the #1 user pain (manual step-cutting).
- **Cost:** medium — wraps existing transient integrators with a residual/error
  estimator + step controller; no new core math.
- **Shape:** a controller layer over `Newmark`/`HHT`, not a new scheme.

### Rank 2 — Explicit bulk viscosity (b1 linear / b2 quadratic) — **BUILD (small)**
- **Gap:** confirmed absent everywhere in OpenSees incl. the fork explicit path;
  Abaqus (b1=0.06, b2=1.2) and LS-DYNA ship it standard to damp the highest
  element mode and smear shocks.
- **Fit:** the fork lives in explicit dynamics + contact (ADR-39/41); bulk
  viscosity is the standard stabilizer there. Pairs with `CriticalTimeStep`.
- **Cost:** low-medium — element-level pressure term in the explicit residual;
  excluded from reported stress. Could be a damping add-on, not a new integrator.

### Rank 3 — Modern implicit composite (Bathe β1/β2, 2012) — **ASSESS THEN MAYBE**
- **Gap:** we have `TRBDF2` (Trap+BDF2) but not the tunable 2-sub-step **implicit**
  Bathe with spectral-radius control — strong for contact/sharp nonlinearity.
- **Fit:** good for contact dynamics; overlaps `TRBDF2`/`HHT`.
- **Cost:** medium. **Action:** first benchmark `TRBDF2` vs Bathe-2012 to confirm
  the gap is real before building.

### Rank 4 — Energy-momentum conserving (Simo-Tarnow / EMC) — **DEFER**
- **Gap:** no code (incl. commercial) ships it natively; research-grade. Value
  for long-duration nonlinear dynamics without damping drift.
- **Fit/cost:** niche / high. Defer unless a project demands it.

### Rank 5 — Explicit↔implicit co-simulation / domain split — **DEFER**
- LS-DYNA's seamless switching; high complexity, low current demand. Defer.

### Not gaps (already covered / fork strengths — do NOT build)
- Generalized-α, HHT, Newmark, central difference, Riks-class arc-length,
  dynamic relaxation, selective/consistent mass scaling, eigen — **covered**.
- Hybrid-sim, consistent-Olovsson SMS, dissipation arc-length, CMOD indirect
  control, FLAC-LNVD explicit quasi-static — **fork strengths beyond the
  commercial baseline**.

---

## Part 5 — Consolidation / cleanup (non-new-integrator work)

Independent of any new build, the study surfaced concrete housekeeping:

- **C1** Merge `ArcLength1` → `ArcLength -linear`; `Newmark1` → `Newmark`.
- **C2** Deprecate `CentralDifferenceAlternative`, `CentralDifferenceNoDamping`;
  fix or remove malformed `ExplicitDifference`.
- **C3** Collapse 6 `ExplicitBathe*` classTags → one class + `-lnvd/-sms/-consistent`.
- **C4** Finish or remove `DistributedDisplacementControl`.
- **C5** Extract DDM sensitivity into a mixin; then **add sensitivity to `HHT`/
  `GeneralizedAlpha`** (prereq for any explicit-integrator sensitivity).
- **C6** Validate or mark-experimental `EQPath`, `HSConstraint`.
- **C7** Audit consistent-SMS parallel safety (finding #7).
- **C8** Reclassify `GimmeMCK` as a utility (not transient integrator) in docs.

---

## Recommendation (one line)

**Build Rank 1 (adaptive implicit Δt) and Rank 2 (explicit bulk viscosity);
assess Rank 3 (implicit Bathe) before committing; defer 4–5. Run the C1–C8
cleanups opportunistically.** New-integrator candidates → promote to ADR-52+.
