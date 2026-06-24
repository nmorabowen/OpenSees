---
title: Integrator Study Workflow — assess all OpenSees integrators, gap-analyse vs Abaqus / Kratos / LS-DYNA
project: Ladruno
status: draft
priority: medium
owner: nmora
tags:
  - implementation
  - integrator
  - study-workflow
  - gap-analysis
  - explicit-dynamics
aliases:
  - Integrator Study Workflow
  - ADR-49
updated: 2026-06-23
---

# ADR-49 — Integrator study workflow

## What

A **repeatable study workflow** for systematically assessing every time/path
integrator in the Ladruno fork (upstream + fork-authored, ~50 classes), scoring
each on a fixed rubric, crosswalking each to its Abaqus / Kratos / LS-DYNA
counterpart, and from the gaps deciding **which new integrators (if any) are
worth building**.

This document is the *method and the census*, not the filled-in results. It
defines:

1. the **inventory/taxonomy** (Section 3 — already populated from the tree),
2. the **per-integrator assessment rubric** (Section 4),
3. the **external crosswalk** procedure and template (Section 5),
4. the **gap-analysis → candidate-ranking** procedure (Section 6),
5. how to **run** it — incrementally by hand or as a one-shot multi-agent
   workflow (Section 7).

**Not in scope here:** actually filling the rubric for all 50 (that is a *run*
of this workflow, Phase 2+); proposing a specific new integrator with a full
design (that is a follow-up ADR seeded by Section 6).

## Why

The fork's active program has pushed hard into **explicit dynamics + selective
mass scaling** (ADRs 04, 05, 35, 36, 37, 38) and **stabilized static
path-following** (ADRs 20, 21, 22). Those were built feature-by-feature without
a single map of *where OpenSees stands integrator-by-integrator versus the
commercial codes we benchmark against*. Before committing to the next
integrator we want one defensible artifact that answers:

- Which integrators are solid, which are stubs/duplicates, which are fork-unique?
- For each, what do Abaqus/Kratos/LS-DYNA do that we don't (order, damping
  control, adaptivity, conservation, cost)?
- Given the fork's goals (nonlinear RC, contact, explicit quasi-statics,
  performance), what is the highest-value *missing* scheme?

Companion docs already cover the fork-authored set in depth — this workflow
*supersets* them with upstream + external comparison:
[[ladruno_integrators_guide]], [[Ladruno_explicit_roadmap]],
[[04_explicit_dynamics_and_energy_balance]], [[40_ladruno_performance_adr]],
[[40a_kratos_crosspollination_amendment]].

## Where

This is a documentation/analysis workflow — it reads source, it does not (yet)
add code.

- New doc: this file + (per run) a results file
  `Ladruno_implementation/49a_integrator_scorecard_<date>.md`.
- Read (assessment source of truth): `SRC/analysis/integrator/*.{h,cpp}` —
  always score from the `.cpp` (`newStep` / `update` / `formTangent` /
  `domainChanged`), never from docs alone.
- Reference skills: `abaqus-theory` → `abaqus-theory-procedures` (Abaqus leg),
  `kratos` (Kratos leg). **No `ls-dyna` skill exists** — the LS-DYNA leg draws
  on the `opensees-performance` HPC library + targeted web research and is
  therefore lower-fidelity/higher-effort (flag every LS-DYNA claim with a
  source).
- Ledgers to reconcile against: [[LEDGER_implementations]] (fork classTags
  33000–33007 + SMS cluster), [[LEDGER_vanilla_files]].

## Section 3 — Inventory & taxonomy (Phase 1, pre-populated)

Base hierarchy: `Integrator` → `IncrementalIntegrator` → {`StaticIntegrator`,
`TransientIntegrator`}; `EigenIntegrator` derives directly from `Integrator`.
Source = `SRC/analysis/integrator/`. **U** = upstream, **F** = fork-authored.

### 3.1 Static path-following (`StaticIntegrator`)

| Integrator | Src | Role / scheme | Notes |
|---|---|---|---|
| `LoadControl` | U | proportional load increment | baseline static |
| `LoadPath` | U | user load-factor path | |
| `StagedLoadControl` | U | staged-construction load control | |
| `DisplacementControl` | U | single-DOF disp control | |
| `DistributedDisplacementControl` | U | parallel disp control | OpenSeesSP/MP |
| `ArcLength` | U | spherical/cylindrical arc-length | limit points |
| `ArcLength1` | U | arc-length variant | duplicate? — confirm in Phase 2 |
| `MinUnbalDispNorm` | U | min unbalanced disp norm | |
| `HSConstraint` | U | hyperspherical constraint arc-length | |
| `EQPath` | U | equilibrium-path control | |
| `HarmonicSteadyState` | U | harmonic steady-state (static base) | |
| `LadrunoArcLength` | F (33004) | stabilized arc-length | ADR-20 |
| `LadrunoIndirectControl` | F (33006) | indirect/displacement-norm control | |
| `LadrunoDissipationArcLength` | F (33007) | dissipation-based arc-length | **scoped, no code** — ADR-22 |

### 3.2 Implicit transient (`TransientIntegrator`)

| Integrator | Src | Family | Notes |
|---|---|---|---|
| `Newmark`, `Newmark1` | U | Newmark-β | `Newmark1` likely legacy — confirm |
| `StagedNewmark` | U | Newmark + staging | |
| `HHT`, `HHT_TP` | U | Hilber-Hughes-Taylor α | `_TP` = "table-passing" hybrid-sim variant |
| `HHTGeneralized`, `HHTGeneralized_TP` | U | generalized HHT | |
| `GeneralizedAlpha` | U | Chung-Hulbert generalized-α | |
| `Collocation` | U | collocation (Hilber) | |
| `WilsonTheta` | U | Wilson-θ | |
| `Houbolt` | U | Houbolt (multistep) | |
| `ParkLMS3` | U | Park 3-step stiffly-stable | |
| `TRBDF2`, `TRBDF3` | U | composite TR-BDF | implicit composite — compare to Bathe |
| `BackwardEuler` | U | backward Euler | L-stable |
| `PFEMIntegrator` | U | particle-FEM fractional step | fluid |
| `GimmeMCK` | U | *utility* — exports M, C, K | not a true integrator; tag separately |

### 3.3 Explicit transient (`TransientIntegrator`)

| Integrator | Src | Scheme | Notes |
|---|---|---|---|
| `CentralDifference` | U | central difference | conditionally stable |
| `CentralDifferenceAlternative` | U | CD variant | |
| `CentralDifferenceNoDamping` | U | CD, no damping | |
| `ExplicitDifference` | U | explicit difference | |
| `ExplicitDifferenceStatic` | U | explicit quasi-static | |
| `NewmarkExplicit` | U | explicit Newmark | |
| `HHTExplicit`, `HHTExplicit_TP` | U | explicit HHT | hybrid-sim |
| `HHTGeneralizedExplicit`, `_TP` | U | explicit generalized HHT | |
| `KRAlphaExplicit`, `_TP` | U | Kolay-Ricles unconditionally-stable explicit | notable |
| `AlphaOS`, `AlphaOS_TP` | U | α operator-splitting | hybrid-sim |
| `AlphaOSGeneralized`, `_TP` | U | generalized α-OS | |
| `CentralDifferenceLadruno` | F (33003) | robust CD | ADR-05 |
| `CentralDifferenceSMS` | F | CD + selective mass scaling | ADR-36 |
| `CentralDifferenceSMSConsistent` | F | CD + consistent SMS | ADR-38 |
| `ExplicitBathe` | F (33000) | Noh-Bathe two-sub-step | ADR-04 |
| `ExplicitBatheLNVD` | F (33002) | Bathe + linear-no-velocity-damping (quasi-static) | |
| `ExplicitBatheSMS`, `ExplicitBatheSMSConsistent` | F | Bathe + SMS | |
| `ExplicitBatheLNVDSMS`, `…SMSConsistent` | F | Bathe LNVD + SMS | |
| `LadrunoDynamicRelaxation` | F (33005) | matrix-free dynamic relaxation | ADR-21 |

### 3.4 Hybrid-simulation cluster (cross-cuts 3.2/3.3)

`*HSFixedNumIter`, `*HSIncrLimit`, `*HSIncrReduct` variants of `Newmark`,
`HHT`, and `Collocation`, plus every `_TP` class. These exist for **real-time
hybrid testing** (Schellenberg et al.) and have **no analogue** in
Abaqus/Kratos/LS-DYNA — flag as a fork/upstream-unique strength, score lightly
(one rubric row for the cluster, not one per variant).

### 3.5 Eigen / shared utilities

| Item | Src | Role |
|---|---|---|
| `EigenIntegrator` | U | modal/eigen RHS assembly |
| `CriticalTimeStep` | F | shared per-element Δt_cr estimator (powers `criticalTimeStep()`) |
| `LadrunoFictitiousMass`, `LadrunoMassLumping`, `LadrunoMassScaling`, `LadrunoMassScalingEnergy`, `LadrunoConsistentRefine` | F | mass-construction helpers consumed by the explicit/SMS integrators — assess as *enablers*, not standalone integrators |

> **Census totals (to verify on each run):** ~38 upstream + ~12 fork "true"
> integrators, plus the HS/`_TP` hybrid-sim cluster and the mass/utility
> helpers. The ~50 figure counts true integrators; helpers and base classes
> are excluded.

## Section 4 — Per-integrator assessment rubric (Phase 2)

Score each integrator (or cluster) from its `.cpp`. One row per integrator in
the results file, these columns:

| Field | What to record | Where to read it |
|---|---|---|
| **Scheme & order** | named method; temporal order of accuracy | header comment + `newStep`/`update` |
| **Stability** | unconditional / conditional; CFL or θ/α-range; L- vs A-stable | algorithm params, `domainChanged` |
| **Numerical damping** | controllable? (ρ∞, α, γ); high-freq dissipation | parameter parsing |
| **Tangent / cost** | needs full K each iter? lumped-mass-only? factorization reuse? matrix-free? | `formTangent`, `formNodTangent` |
| **Convergence robustness** | iterative (needs algorithm) vs one-pass explicit; known stall modes | `update` return paths |
| **Parallel / sensitivity** | OpenSeesSP/MP support; `formSensitivity*` implemented? | `formSensitivityRHS`, MPI guards |
| **Intended regime** | dynamics / quasi-statics / pushover / modal / hybrid-test | docs + usage |
| **Maturity** | shipped / legacy-duplicate / stub / fork-WIP | code + [[LEDGER_implementations]] |
| **Known limitations** | from code comments, ADRs, [[LEDGER_quirks]] | grep `Ladruno`, `// TODO`, `opserr` |
| **Verdict** | keep-as-is / document / deprecate / improve | analyst call |

Scoring guidance: prefer a 1-line evidenced note per field over a number.
Where two classes are near-duplicates (`ArcLength`/`ArcLength1`,
`Newmark`/`Newmark1`), the deliverable is a **consolidation recommendation**.

## Section 5 — External crosswalk (Phase 3)

For each integrator *family* (not each variant), fill one crosswalk row:

| Column | Source | Fidelity |
|---|---|---|
| OpenSees integrator(s) | this tree | — |
| **Abaqus** equivalent | `abaqus-theory-procedures` (HHT-α default implicit, generalized-α, explicit central-difference, Riks arc-length) | high (skill) |
| **Kratos** equivalent | `kratos` skill — `Scheme`/`SolvingStrategy`: `ResidualBasedNewmarkScheme`, `ResidualBasedBossakDisplacementScheme`, `GeneralizedAlpha`, explicit strategies | high (skill) |
| **LS-DYNA** equivalent | opensees-performance HPC refs + web (explicit central-difference w/ hourglass, implicit HHT, *CONTROL_IMPLICIT) | **low — cite every claim** |
| **They have, we don't** | the deltas | — |
| **We have, they don't** | e.g. hybrid-sim `_TP`/HS family, dissipation arc-length | — |

Known anchor mappings to seed the table:
- Newmark/HHT/generalized-α ↔ all three codes have these (Abaqus default = HHT-α;
  Kratos = Bossak/generalized-α; LS-DYNA implicit = HHT).
- Central difference ↔ Abaqus/Explicit & LS-DYNA core explicit solver
  (their advantage: built-in **hourglass control**, **bulk viscosity**,
  **automatic mass scaling**, **contact-aware Δt_cr**).
- Riks arc-length ↔ `ArcLength`/`LadrunoArcLength`.
- Bathe composite (implicit) ↔ we have **explicit** Noh-Bathe only; the
  **implicit** Bathe composite (≈ `TRBDF2`) is a gap to assess.

## Section 6 — Gap analysis → candidate ranking (Phase 4)

From Section 5's "they have, we don't", build a candidate list and rank each by:
(a) fit to fork goals (nonlinear RC, contact, explicit quasi-static,
performance), (b) build cost / reuse of existing machinery, (c) overlap with
something we already ship.

**Seed candidates to evaluate (not commitments):**

1. **Implicit composite / Bathe (β₁/β₂) implicit** — robust for contact &
   sharp nonlinearity; partially covered by `TRBDF2`. *Assess overlap first.*
2. **Automatic / adaptive time stepping** — Δt control on residual or local
   error (Abaqus/LS-DYNA `*CONTROL_IMPLICIT` auto-step). High value for the
   robust-solve driver (ADR-31).
3. **Energy-momentum conserving** (Simo-Tarnow / EMC) — long-duration nonlinear
   dynamics without artificial damping drift.
4. **Bulk viscosity + hourglass-aware explicit Δt** — to close the LS-DYNA/Abaqus
   explicit-quality gap, building on `CriticalTimeStep`.
5. **Co-simulation / explicit-implicit domain split** — niche; likely defer.

Output: a ranked shortlist, each as a one-paragraph stub that, if green-lit,
becomes its own design ADR (52+).

## Section 7 — How to run this workflow

Two modes; pick per appetite.

**A. Incremental (default, no agent spend).** Work family-by-family from
Section 3. For each family: read the `.cpp`s → fill Section-4 rows → fill the
Section-5 crosswalk row (invoke `abaqus-theory`/`kratos` skills) → append to
`49a_integrator_scorecard_<date>.md`. Do Section 6 once 3.1–3.4 are done.

**B. One-shot multi-agent workflow** (when ready to spend tokens). Fan out one
agent per family (3.1, 3.2, 3.3, 3.4+3.5) each returning a structured
Section-4 table; a second stage runs the Section-5 crosswalk per family
(Abaqus+Kratos via skills, LS-DYNA via web); a synthesis agent produces
Section 6. Gate the LS-DYNA leg as "best-effort, cited". Trigger only on
explicit opt-in ("use a workflow" / ultracode).

## Risks / open questions

> [!question]
> LS-DYNA leg has no skill — accept lower fidelity (cited web claims), or drop
> LS-DYNA to a "notes" column and treat Abaqus+Kratos as the authoritative
> baselines?

> [!question]
> Do we score the ~20 hybrid-sim (`_TP`/HS) variants individually, or as one
> cluster row? (Recommended: cluster — they're real-time-hybrid-testing
> machinery, not analysis integrators, and have no commercial-code analogue.)

> [!question]
> Confirm suspected duplicates/legacy (`ArcLength1`, `Newmark1`,
> `CentralDifferenceAlternative`) — consolidate or keep?

- Source-of-truth discipline: score from `.cpp`, reconcile names against
  `SRC/classTags.h` and [[LEDGER_implementations]] before publishing.
- The fork-authored set is already documented in [[ladruno_integrators_guide]];
  this study must *reference* not *duplicate* it.

## Implementation log

- **2026-06-23 — first full run.** Phases 1–4 executed. Phase 2 (assessment of
  all ~50 integrators) by four source-reading agents; Phase 3 crosswalk via
  `abaqus-theory-procedures` + `kratos` skills (LS-DYNA knowledge-based);
  Phase 4 gap analysis + ranking. Result: [[49a_integrator_scorecard_2026-06-23]].
  Headline: **build adaptive-implicit-Δt (Rank 1) and explicit bulk viscosity
  (Rank 2); assess implicit Bathe (Rank 3); defer EMC + co-sim.** Plus 8
  consolidation items (C1–C8: ArcLength1/Newmark1 merges, CentralDifference
  cruft removal, ExplicitBathe flag-collapse, sensitivity mixin, etc.).
