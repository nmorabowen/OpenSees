---
title: "ADR 98 — MPM in OpenSees: what twelve libraries teach, the candidate architectures, and how the best one is chosen"
project: Ladruno
status: scoping — NO code. REVISED 2026-09-25 after adversarial review (§15) — not clean on first pass (1 blocker, 5 majors), repaired here. No provisional winner (A 75 · F 74 · C 72); one prerequisite (§6-0 integrator re-seed contract); M0 re-scoped (§10); decisions requested (§11)
priority: high
owner: nmora
amends: 86_geomechanics_libraries_and_mpm_scoping_report
tags:
  - adr
  - scoping
  - mpm
  - material-point-method
  - large-deformation
  - soil
  - particles
aliases:
  - ADR-98
  - MPM scoping
updated: 2026-09-25
---

# ADR-98 — MPM in OpenSees: twelve libraries, the candidate architectures, and how the best one is chosen

> **Number 98.** 86 is held twice (SANISAND and the geomechanics/MPM report), 87 is
> "depth with width", 88 is the state-transfer draft, and 89 is still named for TIM
> calibration in ADR-86's (now stale) §6.4 text. 89, 98 and 99 were checked free on
> every branch on 2026-09-25; 98 avoids the ADR-86 reference.
>
> **Scope.** ADR-86 D2 decided MPM proceeds as **A + B** — an in-tree lane *and* the
> FE→MPM handoff — and required M0 to re-confirm that decision on corrected facts
> (ADR-86 §5.5). This is the M0 scoping study the owner asked for: *study the other
> MPM libraries, and let the best idea win.* No code. Under ADR-87 every unmeasured
> claim below is a claim, and says so; §10 M0b exists to measure the ones that decide.
>
> **Revision.** The first draft named C the provisional winner at 83/100. A Fable
> adversarial review the same day (§15) found a **blocker** — OpenSees integrators
> cannot be fed grid state through the Nodes (§5 W1) — plus five majors. Every
> finding the author re-checked in-tree held. This revision withdraws C's lead,
> adds a prerequisite (§6-0) and a new finalist (F), and re-scopes M0.
>
> **Dependency status.** ADR-88 (state transfer, Duvaut-Lions) is an **uncommitted
> draft** in the owner's main checkout — it is *not* on `ladruno` as of 2026-09-25.
> Every ADR-88 item cited here (`getHistory/setHistory`, `MATERIAL_STATE`, the
> Duvaut-Lions wrapper) is a dependency on a draft, and the `[[88_…]]` links resolve
> only once it lands. Nothing before M5 depends on it.

---

## 0. Verdict

| Question | Answer |
|---|---|
| **Build MPM natively?** | **Yes, still** — but not for free. The pure kernels (`LadrunoHypoKernel.h`, `LadrunoUPKernel.h`) and the `NDMaterial` library do run unchanged at particles, and a per-step hook exists in every serial driver (`Domain::analysisStep`). What does **not** exist is a way to hand grid state to an integrator: every transient integrator caches its own nodal vectors and re-reads the Nodes only in `domainChanged()` (§5, W1). |
| **What does the field agree on?** | Seven settled ideas (§3) — adopted, not debated. |
| **Prerequisite** | **The integrator re-seed contract (§6-0):** `reseedFromCommitted()` on `TransientIntegrator` plus a `Domain` flag consumed by the drivers. Every candidate that uses OpenSees integrators needs it (A, B, C, and F's later implicit lane). |
| **Winner** | **Not declared.** Revised scores **A 75 · F 74 · C 72** (E 69, B 68, D 57) — three finalists inside the scoring's noise. The first draft's C at 83 rested on two refuted claims (W1; threaded soil updates) and one omitted cost (node-set churn). **M0b decides on measured numbers.** |
| **Finalists** | **A** particle-as-element (the Kratos implicit pattern). **F** a self-stepping explicit engine — OpenSees `NDMaterial` + hypo kernel at particles, no OpenSees integrator, SOE or handler for explicit steps — with implicit MPM via A or C later: the architecture of Berkeley, CB-Geo and Kratos-explicit, added by the review. **C** sparse cell kernels materialized on demand. |
| **Rejected** | **B** dense grid (1–2 orders of magnitude of waste on sparse runout problems); **D** vendor Berkeley's MIT engine (stateless-material contract, its own solve loop, no MSVC path, permanent drift) — kept as a cross-code oracle; **E** co-simulation only (excluded by ADR-86 D2); **C′** handler-emitted adapters (§6, rejected on cost). |
| **The differentiator, at the strength the evidence supports** | No surveyed open-source MPM code ships a cyclic-liquefaction sand model in its own library. The field reaches such models through user DLLs (Anura3D's external-soil-model route) or its own implementations (UBCSAND in CPDI-MPM, Giridharan et al. 2020). **PM4Sand / SANISAND / PDMY at material points was not found in the published record** — and the one group that had PM4Sand at the FEM end chose not to carry it into MPM, citing MPM's stress inaccuracy (§12). |
| **First step** | M0-pre (the §6-0 contract), then a re-scoped M0b of **measurements, not spikes** (§10): A's rebuild cost on the existing fork with a synthetic churn deck, the re-seed cost, and F's explicit engine as the numpy oracle's C++ twin. |

---

## 1. MPM in OpenSees's own vocabulary

MPM is **updated-Lagrangian FEM whose quadrature points are persistent material
points, solved on a background grid that is reset every step.** Stated that way,
the candidates in §6 become comparable, and it is obvious which half is new.

| MPM step | What happens | OpenSees / fork analogue | Exists? |
|---|---|---|---|
| 1. Locate | assign each particle to a grid cell | — | **new** |
| 2. P2G | project mass, momentum (and inertia, implicit) to grid nodes; the grid starts from zero | **none that reaches an integrator** — `ladrunoSetNodeTrial` (ADR-52) sets trial node state and never commits, and integrators ignore node state between steps anyway (§5 W1) | **new — needs §6-0** |
| 3. Grid solve | explicit update, or Newton on nodal increments | `Integrator` + `EquiSolnAlgo` + `LinearSOE` | **yes** |
| 4. Constitutive | strain increment at each particle → stress, tangent | `NDMaterial::setTrialStrain/getStress/getTangent`, with `LadrunoHypoKernel::hypoIncrement` for the objective update | **yes** |
| 5. G2P | update particle x, v, a, F, V | — | **new** |
| 6. Reset | the grid forgets everything | — | **new** |

Steps 3 and 4 are the mechanics that make OpenSees worth anything to MPM, and both
exist. The first draft called the rest "bookkeeping"; the review showed step 2 is
not — handing the P2G state to OpenSees's integrators is itself new infrastructure
(§6-0), unless the explicit update is owned by the MPM engine instead (candidate F).

---

## 2. The field — twelve codes

Licences verified 2026-09-25 from each repository's own licence file, the GitHub
licence API, or Crossref. **"Source"** = the fork may read and reuse code with
attribution; **"papers"** = ideas from the published literature and documentation
only — never open the source while writing fork code (the ADR-86 §2 rule, extended).

| Code | Group | Lang | Licence | Fork may use | Solvers | Bases | Multiphase | Contact | Parallel | Activity |
|---|---|---|---|---|---|---|---|---|---|---|
| **Kratos MPMApplication** | CIMNE | C++/Py | **BSD-4** (advertising clause; no app-level file ⇒ `kratos/license.txt`) | **source** + acknowledgement | explicit USF/USL/MUSL; implicit Bossak; quasi-static; static | linear; PQMPM (explicit only) | mixed u-p element | particle Dirichlet (penalty / Lagrange), Neumann, slip; partitioned FE coupling | OpenMP; MPI = a particle-transfer primitive only | active; checked out locally |
| **geomechanics/mpm** | Berkeley (Soga, Kamrin) | C++17 | **MIT** | **source** | explicit USF/USL/MUSL; **implicit Newmark** (Newton; Eigen / PETSc); semi-implicit Navier–Stokes and two-phase; explicit two-phase | linear, GIMP, B-spline, LME | explicit + semi-implicit two-phase | multi-material frictional | MPI + load balancer + halo particle transfer (dynamic serialization); OpenMP | very active (pushed 2026-09-23) |
| **CB-Geo mpm** | Cambridge / UT Austin | C++17 | **MIT** | **source** | explicit USF/USL | — | — | frictional | MPI + TBB; *replicated* grid + `all_reduce` (2019 paper) | dormant since 2024-11 (continued as Berkeley) |
| **Uintah** | Utah | C++ | **MIT** (`src/LICENSE`) | **source** | explicit; ICE coupling (literature) | GIMP, CPDI (literature) | ICE | multi-material | asynchronous many-task; "hundreds to hundreds of thousands of processors" (README) | active |
| **Karamelo** | Monash | C++ | **GPL-2** (no "or later" found ⇒ GPL-2-only presumption) | papers | ULMPM + **TLMPM** | linear, cubic B-spline, Bernstein (TLMPM) | — | — | MPI; Kokkos GPU branch (2024) | dormant (last push 2025-03-31) |
| **NairnMPM** | Oregon State | C++ | **none found** (no licence file; docs say only "open-source" ⇒ all rights reserved) | papers | explicit; FMPM(k) | splines | — | multimaterial frictional, imperfect interfaces, rigid particles; explicit cracks + cohesive zones | — | active (2026-04) |
| **Anura3D** | community (Deltares origin) | Fortran | **LGPL-3** | papers (ADR-86 §2) | explicit dynamic; quasi-static | linear; MPM-MIXED / UL-FEM | 2-phase, 3-phase, **double-point** | — | — | active; no releases |
| **AMPLE** | Durham (Coombs, Augarde) | MATLAB | **GPL-3** | papers | **implicit quasi-static** (literature) | GIMP, CPDI (literature) | — | — | serial | teaching / research |
| **Taichi MLS-MPM** | MIT CSAIL / graphics | Taichi, C++ | **Apache-2.0** | **source** | explicit | quadratic B-spline (MLS) | — | — | GPU; sparse grids | active |
| **GeoTaichi** | Shi & Guo | Python / Taichi | **GPL-3** | papers | explicit MPM + DEM | — | — | DEM coupling | GPU | active (2026-08) |
| **Matter** | ETH (Blatny) | C++ | **GPL-3** | papers | explicit, finite strain; DP, overstress, critical state, μ(I) (GMD 2025) | — | — | — | — | active (2026-09) |
| **Sparse framework** | arXiv 2605.28525 (2026) | — | code not stated | idea | — | — | — | — | active-node indexing: scan (CPU), hash (GPU); **1–2 orders of magnitude** less time and memory in strongly sparse cases | 2026 |

Constitutive libraries (the differentiator, §0): Kratos MPM ships linear elastic,
hyperelastic, Hencky J2 plasticity, Hencky Mohr–Coulomb (+ strain softening, + u-p),
Hencky–Borja Cam-Clay, Johnson–Cook, Newtonian; Berkeley ships linear elastic,
Hencky hyperelastic, Mohr–Coulomb, Modified Cam-Clay, NorSand, Bingham
(+ viscoplastic), Newtonian. Kratos MPM elements also refuse the GeoMechanics
application's UDSM wrapper at `Check()` (strain-measure mismatch, review §F4).

**Licence consequences, beyond ADR-86 §2.** (i) **Kratos's BSD-4 advertising clause
is GPL-incompatible**: Kratos code may enter the fork (with the acknowledgement
"This product includes Kratos Multi-Physics technology" in advertising material)
but **may not enter apeGmsh**, which is **GPL-3**. (ii) Conversely, apeGmsh *may*
take MIT, Apache, LGPL and GPL-3-compatible material that the fork may not — **but
not GPL-2-only code** (Karamelo, as its licence reads). So "apeGmsh owns geometry,
seeding and post; the fork owns mechanics" (§8) is also the licence-clean split.
(iii) Uintah is MIT and may be read; NairnMPM stays papers-only. (iv) Papers-only
for LGPL code is fork *policy*: the fork's `COPYRIGHT` terms bar incorporating
GPL/LGPL source, while LGPL itself would permit dynamic linking — the policy is
deliberately stricter. An engineer's reading, not legal advice (ADR-86 §2).

---

## 3. Settled — what every serious code agrees on

Adopted as design constraints. Not re-litigated here.

- **S1 — Particles own all history; the grid owns none and is reset every step.**
  Kratos clears nodal mass, momentum, inertia, displacement, velocity and
  acceleration in the scheme's `Predict()` before P2G
  (`MPMApplication/custom_strategies/schemes/mpm_residual_based_bossak_scheme.hpp:237-259`,
  P2G from `:281`); Berkeley re-initialises nodes every step. Note *where* Kratos
  does it: in the scheme — i.e. inside the integrator. That is why Kratos has no
  equivalent of W1.
- **S2 — Solve over active nodes only.** *Implicit / static:* Kratos rebuilds the
  DOF set every step (`python_scripts/mpm_solver.py:172` → `Clear()` →
  `SetDofSetIsInitializedFlag(false)`,
  `kratos/solving_strategies/strategies/residualbased_newton_raphson_strategy.h:707`);
  Berkeley's `reinitialise_matrix()` calls `assign_active_nodes_id()` "at the
  beginning of every time step" (`include/solvers/mpm_implicit.tcc:337-352`).
  *Explicit:* there is **no DOF set or global system at all** — Kratos's explicit
  strategy runs on a `dof_set_dummy` and its `Clear()` only clears the scheme
  (`custom_strategies/strategies/mpm_explicit_strategy.hpp:197-218`); per-step
  explicit cost is the particle search plus a nodal reset. (Corrected after review.)
  The 2026 sparse-framework paper measures 1–2 orders of magnitude on sparse problems.
- **S3 — Implicit MPM = Newton on grid increments, particles frozen at step-start
  positions, G2P at convergence.** Kratos: element `FinalizeSolutionStep` commits
  the law, then `UpdateGaussPoint` moves the particle with Guilkey–Weiss
  trapezoidal kinematics (`custom_elements/mpm_updated_lagrangian.cpp:903, :984`);
  Berkeley: Newmark predictor/corrector (`mpm_implicit.tcc:182-231`).
- **S4 — Linear hat functions are not production grade** (cell-crossing error).
  Every serious code ships an alternative: Kratos PQMPM (partitioned quadrature,
  Wilson et al. 2021 — explicit only); Berkeley GIMP, B-spline, LME; Karamelo
  B-spline, CPDI; NairnMPM splines.
- **S5 — Near-incompressibility needs a remedy.** Kratos: mixed u-p elements and
  `_UP` laws; Berkeley: `particles/anti_locking/particle_bbar`.
- **S6 — Search = neighbour-first, spatial bins as fallback** (Kratos
  `custom_utilities/material_point_search_utility.h:171-205, :286-351`, both
  OpenMP-parallel); structured grids use index arithmetic (Berkeley, Uintah, Taichi).
- **S7 — Body-to-body contact = multi-velocity field with Coulomb friction**
  (Berkeley `map_multimaterial_*` + `contacts/contact_friction`; NairnMPM
  multimaterial mode; Bardenhagen et al. 2000).

---

## 4. Contested — where the best idea has to win

| # | Choice | Positions in the field |
|---|---|---|
| C1 | Particle representation | **particle-as-element** (Kratos) vs **particle store + cell kernel** (Berkeley, CB-Geo, Uintah, Taichi) |
| C2 | Constitutive binding | **stateful law per particle** (Kratos — and OpenSees `NDMaterial`) vs **shared stateless law + per-particle state vector** (Berkeley: `dense_map`) |
| C3 | Transfers | PIC/FLIP blend (most) vs APIC/MLS (Taichi lineage) vs FMPM / XPIC (NairnMPM) |
| C4 | Kinematics | ULMPM (everyone) vs TLMPM (Karamelo — no cell crossing, but no mixing) |
| C5 | Grid storage | dense (CB-Geo, Kratos grid mesh) vs sparse / hashed (Taichi, sparse framework) |
| C6 | Two-phase | single-point (Berkeley, Kratos UP) vs double-point (Anura3D) |
| C7 | Host | inside a general FE framework (Kratos) vs standalone (all others) |
| C8 | Who owns the explicit update | the framework's integrator (Kratos, via its scheme) vs the MPM engine's own loop (Berkeley, CB-Geo, Uintah, Taichi) — *added after review* |

**The fork's answer to C2 is forced** — `NDMaterial` is stateful — which makes
Kratos the most relevant precedent. **C7 and C8 are not forced**, and the review
showed they matter most: Kratos can hand grid state to its integrator because P2G
*lives in* its scheme; OpenSees has no such seam (W1). So C1, C5 and C8 are the
real contest; C3, C4 and C6 are options inside whichever architecture wins.

---

## 5. What the fork already has — and the walls

Verified in-tree 2026-09-25 (by the author, and independently by the review).

| Need | Fork asset | Evidence |
|---|---|---|
| Constitutive update at particles | the whole `nDMaterial` library — PM4Sand, PM4Silt, ManzariDafalias, LadrunoSANISAND (+ IMPL-EX, ADR-92), PDMY, J2, DP, ASDPlastic … | `SRC/material/nD/` |
| Large-strain objective update | **`LadrunoHypoKernel.h`** — "the PURE numerical core of `-geom hypo`… NO OpenSees dependency… each consuming element supplies only its reference shape gradients and nodal displacements"; Hughes–Winget midpoint | `SRC/element/solidTransformation/LadrunoHypoKernel.h:114` `dispGradient`, `:174` `hypoIncrement`, `:280` `pushStress6` |
| Two-phase blocks | **`LadrunoUPKernel.h`** — `addQ/addH/addS/addHtilde/addFseep` over `(N, dN, dv)` arrays, `dofMap` for mixed ndf | `SRC/element/ladrunoUP/LadrunoUPKernel.h:87-203` |
| A per-step hook before the domain-changed check and the predictor | `Domain::analysisStep(dT)` — virtual, returns 0, called at the top of each step by every **serial** driver (Static, DirectIntegration, VariableTimeStep, PFEM, the ADR-73 staggered driver — verified by review); an add/remove inside it is consumed by the same step's `hasDomainChanged()`. In parallel, `PartitionedDomain::analysisStep` calls it once (`:1159`) but `Subdomain`/`ShadowSubdomain` override and delegate to a decomposition-analysis step, so **a fork hook never fires inside a Subdomain**. **It cannot hand state to an integrator (W1).** | `Domain.cpp:2768`; `DirectIntegrationAnalysis.cpp:198 → 206 → 220`; `StaticAnalysis.cpp:136 → 150 → 166`; `PartitionedDomain.cpp:1131-1174`; `Subdomain.cpp:808` |
| An engine attached to the `Domain` | `theContactDomain`: `Domain::commit()` calls its `commit()` and raises `domainChange()` on `needsResort` (ADR-60, #444); also `revertToStart`, `sendSelf/recvSelf` (ADR-78 P2, #739) | `LEDGER_vanilla_files.md:99-101`; `Domain.cpp:2302-2306` |
| Mid-analysis add/remove of Nodes and Elements | `domainChange()` → driver `domainChanged()` rebuilds handler, numberer, SOE and **re-seeds every integrator from committed node state**; PFEM `BackgroundMesh` adds and removes Nodes and Elements at runtime | ADR-88 Context; `BackgroundMesh.cpp:799/813, :1082, :2818` |
| Threaded constitutive updates | **WP-107** OpenMP over `Domain::update` (opt-in `ladrunoThreads`) — **gated by an allowlist that contains no soil model**: only `ElasticIsotropicPlaneStrain2D` returns true; `ManzariDafalias` refuses "unconditionally" (4/4 segfaults in the plastic branch, root cause not located) and `LadrunoSANISAND` delegates to it; PDMY's process-static arrays are unsafe by construction; the allowlist is all-or-nothing per Domain. **Particle constitutive updates do not thread today.** (Corrected after review.) | `ElasticIsotropicPlaneStrain2D.h:96`; `ManzariDafalias.cpp:1094`; `LadrunoSANISAND.cpp:1948-1965`; `107_ladruno_openmp_element_loop.md` §3.3, §5.1 |
| Explicit stack | CentralDifferenceLadruno, `system Diagonal`, `criticalTimeStep`, HRZ / selective mass scaling | ADR-04/05/35–38 |
| Particle-state serialization (handoff, migration) | ADR-88 D1 `getHistory/setHistory` — opaque, versioned, loud −1 (**draft**) | `88_ladruno_state_transfer_adaptivity_adr.md` D1 |
| Regularization for softening | ADR-88 D3 Duvaut-Lions wrapper (**draft**) | ADR-88 D3 |
| Persistent-state + stateless-views pattern | the contact engine: `LadrunoContactDomain` owns entity-keyed state; `LadrunoContactFE` adapters rebuilt per `handle()`; bucket grid — "nodes (particles) persistent and carry all state; elements disposable connectivity" | `60a_pfem_crosspollination_amendment.md:30-47` |
| Seeding and post | apeGmsh (ours, GPL-3), sanctioned as a seeder incl. particle models | `60a_…:134-138` |

| # | Wall | Consequence | Evidence |
|---|---|---|---|
| **W1** | **Every transient integrator caches its own nodal vectors and re-reads the Nodes only in `domainChanged()`** (CentralDifferenceLadruno also on `revertToLastStep`) | nodal state written between steps is **silently discarded** at `newStep()` — P2G and the grid reset cannot be delivered through the Nodes. **Blocker for every candidate that uses OpenSees integrators** until §6-0 lands. Static integrators hold no nodal vectors — quasi-static survives | `Newmark.cpp:206-208` (predict from its own `U, Udot, Udotdot`), `:398-406` (the only re-read), `:497`; `CentralDifferenceLadruno.cpp:360, :518` (`seedFromCommitted` only in `domainChanged` / `revertToLastStep`); same pattern in HHT, GeneralizedAlpha, CentralDifference, DynamicRelaxation (review §F1) |
| W2 | One `ConstraintHandler` per model; `handle()` is the only legal `FE_Element` factory; **the contact handler is Plain-only** | handler-emitted MPM adapters would have to extend `LadrunoContactHandler` (possible — C′, §6) and inherit its Plain-only restriction: an MPM + FE model with `equalDOF` / rigid links is Plain-style or nothing | `LadrunoContactHandler.h:22-29`; Node-less `DOF_Group`s exist (`LagrangeDOF_Group.h:50-52`, review §F7) |
| W3 | OpenSees makes a `DOF_Group` for **every** Domain node | empty or empty-support grid nodes are singular unless excluded (SP ⇒ a `domainChange()` per add/remove) or ghost-regularized | `PlainHandler.cpp:123-132` |
| W4 | **Stock recorders cache `Node*` / `Element*`, and `Domain` never calls `Recorder::domainChanged()`** | retiring a grid node or cell that a stock recorder watches is a use-after-free; only `LadrunoRecorder` re-maps on the domain stamp | `Domain.cpp:2238, 2357, 2430, 3101` (record / restart / flush only); `NodeRecorder.cpp:1159-1187`, `LadrunoRecorder.cpp:398-424` (review §F6) |
| W5 | Removal cleans nothing | `Domain::removeElement/removeNode` leave loads, SP/MP constraints and recorder pointers behind; `remove element`'s cleanup is commented out; an SP on a missing node `exit(-1)`s | `Domain.cpp:1198-1241`; `OpenSeesMiscCommands.cpp:269-275`; `SP_Constraint.cpp:363-374` (review §F6) |
| W6 | `PartitionedDomain::repartition()` is an empty stub; and the subdomain `domainChange()` loop in `PartitionedDomain::analysisStep` is dead code (its `while` requires the flag `false` inside an `if` that requires it `true`) | no dynamic load balancing for MPI MPM; subdomains are never told the domain changed through that path (**upstream bug**) | ADR-88 Context; `PartitionedDomain.cpp:1147-1156` |
| W7 | PFEM `BackgroundMesh`: singleton, no `sendSelf`, SPH-average P2G, PIC G2P on a different basis | precedent that the Domain tolerates add + remove churn — **not** proof that recorders or constraints survive it; not reusable code | ADR-86 §5.1 (corrected) |

---

## 6. The candidates

### 6-0 — The prerequisite: an integrator re-seed contract (added after review)

W1 means MPM cannot hand its grid state to OpenSees by writing Nodes. The minimal
fix (review §F1) is integrator-side:

- `virtual int TransientIntegrator::reseedFromCommitted()` — the DOF loop that
  already exists inside every `domainChanged()`, and as
  `CentralDifferenceLadruno::seedFromCommitted` (`:473-497`);
- a `Domain` flag beside `hasDomainChangedFlag`, raised by the MPM engine's
  `beginStep()` after it writes the P2G velocity / acceleration and zero
  displacement as **committed** node state (`setTrialX` + `commitState`);
- consumed by the drivers right after `analysisStep` (per the review:
  `DirectIntegrationAnalysis.cpp:198-217`,
  `VariableTimeStepDirectIntegrationAnalysis.cpp:87-97`,
  `LadrunoStaggeredDriver.cpp:205-222`).

O(N) per step, no rebuild. Vanilla touches: `TransientIntegrator.h`, Newmark / HHT /
GeneralizedAlpha / CentralDifference, the fork's CentralDifferenceLadruno and
DynamicRelaxation, the drivers → `LEDGER_vanilla_files` rows. Static integrators need
nothing. An explicit-only alternative — injecting the P2G momentum as an equivalent
nodal force `M·v_P2G / Δt` through the cells — works for lumped central difference,
not for Newmark predictors. Either way, "works with any integrator" becomes **"works
with the re-seed-aware integrators"**, and A, B and C all pay for it; F's explicit
lane does not, because F owns its own update. The same hook is the missing primitive
for any future *mid-analysis* state injection (ADR-88 injects only into fresh
models, which is why it never met W1).

### A — Particle-as-element (the Kratos implicit pattern, native)

A new Element `LadrunoMPMPoint` per material point, owning one `NDMaterial` copy,
one integration point, and as its external nodes the support nodes of the cell it
currently occupies (4 / 8 linear; 9 / 27 quadratic B-spline; up to 16 / 64 GIMP).
When a particle changes cell its connectivity is re-pointed and the `Domain` is
marked changed, so the driver runs a full `domainChanged()` — AnalysisModel
`clearAll`, handler, numberer, SOE graph and symbolic factorization. For implicit
this is what Kratos pays with `Clear()` every step. **For explicit it has no
precedent** (corrected after review): Kratos-explicit keeps no DOF set at all (S2),
so A's explicit cost — `FE_Element` creation, numbering and `setSize` at every
crossing — would be worse than any surveyed explicit code. On a rebuild step
`domainChanged()` also re-seeds the integrator, so A escapes W1 exactly when it pays
the rebuild; on non-rebuild steps it needs §6-0 like everyone else.

- **Inherits for free:** per-particle recorders (one GP per element — the
  LadrunoRecorder's shipped `GLOBAL_GP_COORDS` group *is* the particle cloud, and
  ADR-88's proposed `MATERIAL_STATE` group would carry its history),
  `setParameter`, `sendSelf`, both interpreters, any handler.
- **Costs:** an Element and an `FE_Element` per particle; empty grid nodes still
  need exclusion (W3).
- **Killer risk — explicit throughput.** In a flowing body with ~16 particles per
  cell, *some* particle crosses a cell boundary almost every explicit step, so A
  rebuilds almost every step. Whether the rebuild dominates depends on the
  constitutive cost per particle — and since soil models do not thread (§5), that
  cost is serial. **Unmeasured — M0b(i).**

### B — Cell kernels on a fixed dense grid

Grid Nodes and cell Elements created once for the whole box; cells integrate over
whatever particles are inside; empty nodes carry ghost stiffness / mass; nothing is
ever rebuilt — but P2G still goes through §6-0.

- Wins explicit throughput on dense problems.
- **Loses on the target.** Implicit factorization covers the whole box (structural
  zeros cost as much as nonzeros in sparse direct factorization), memory holds an
  OpenSees `Node` for every grid point, and runout — ADR-86 §12.3's use case — is
  exactly the strongly sparse regime where the 2026 paper measures 1–2 orders of
  magnitude of waste. **Rejected.**

### C — Sparse cell kernels, materialized on demand

- **Virtual grid.** A structured Cartesian grid known only by arithmetic (origin,
  spacing, counts). No objects.
- **Materialization.** The engine instantiates real OpenSees `Node`s and
  `LadrunoMPMCell` Elements only for *active* cells (cells whose support touches a
  particle) plus a one-cell halo, and retires a cell only after it has been empty
  for K steps. The *cell* set changes only when the material front moves beyond the
  halo.
- **Particles.** An engine-owned store, `LadrunoMPMDomain` — one per body, **not a
  singleton** — attached to the `Domain` the way `theContactDomain` is. Each particle
  holds x, v, a, V, m, F, phase data and one `NDMaterial` copy, and is *owned* for
  commit/revert by exactly one cell, so the standard `Domain::commit()` /
  `revertToLastCommit()` reach every particle's material through the cells.
- **Cells.** Fixed connectivity (the support patch of the chosen basis).
  `update()` forms each owned particle's strain increment through
  `ladruno_hypo::hypoIncrement` and calls its `NDMaterial` (WP-107 would thread this
  loop only for allowlisted materials — none of the soil models, §5).
  `getResistingForce / getTangentStiff / getMass` integrate over the particles in
  the cell's support.
- **Orchestration through `Domain::analysisStep(dT)` + §6-0.** A small fork hook
  calls the engine's `beginStep(dT)`: G2P from the last *committed* step, re-bucket,
  materialize / retire (marking the domain changed if the set changed), P2G written
  as *committed* node state, reset, and the re-seed flag raised. The driver's own
  `hasDomainChanged()` rebuilds if needed, the re-seed pulls the P2G state into the
  integrator, `newStep()` predicts from it, the algorithm iterates, `commit()`
  commits. **Works with the re-seed-aware integrators** (the first draft's "works
  with any integrator" was refuted, W1), any algorithm, any handler — including
  `constraints LadrunoContact`. Serial only: the hook never fires inside a
  `Subdomain` (§5).
- **Inherits:** every re-seed-aware integrator, every algorithm (KrylovNewton, line
  search, arc-length, IMPL-EX control), every solver (PARDISO, MUMPS, UmfPack,
  Diagonal), the explicit stack (lumped cell mass, `criticalTimeStep`), energy
  balance, both interpreters, the apeGmsh emitter.
- **Costs:** the most machinery of the native options — virtual grid,
  materialization with hysteresis, particle store and ownership, a particle recorder,
  §6-0, and a small vanilla touch (`Domain::analysisStep`) beside the three
  contact-engine hooks in `Domain.cpp`. **Omitted from the first draft (review §F6):**
  (a) the *node* set churns faster than the *cell* set — nodes inside the active
  region empty and fill as particles move, so they need either an engine-managed SP
  lifecycle (a `domainChange()` per add/remove, W3) or ghost mass / stiffness plus a
  mass cutoff; "particle motion inside the active region is free" holds for cells,
  not for the nonsingular node set; (b) **only `LadrunoRecorder` may watch grid
  tags** — stock recorders dangle on retirement (W4), so the engine must refuse them;
  (c) gravity must be a particle body force through P2G — "cells carry no elemental
  loads" is a design *constraint* that W5 makes mandatory, not a fact; (d) never
  reuse a retired tag within a run.
- **Risks:** Domain add/remove churn (W7); under GIMP / B-spline a particle's
  *support* spans several cells while its *ownership* is one — integration and commit
  must not be conflated.

### D — Vendor Berkeley's MIT engine + an `NDMaterial` adapter

Link `geomechanics/mpm` into the fork and write a `Material<Tdim>` adapter that
drives OpenSees `NDMaterial`s.

- **Wins:** the most capability per week — explicit / implicit / semi-implicit,
  two-phase, GIMP / B-spline / LME, B-bar, frictional multi-material contact, MPI
  with load balancing and halo transfer, HDF5 I/O.
- **Loses on fit.** (1) C2 mismatch: Berkeley's law is shared and stateless with
  state in a `dense_map` (`include/materials/material.h:72, :97`); an adapter must
  own one stateful `NDMaterial` per particle and map trial/commit semantics onto a
  per-iteration `compute_stress`; its **HDF5 checkpoint record** caps state at a
  fixed `svars[20]` (`pod_particle.h:55`) — far short of SANISAND's or PDMY's history
  (ADR-88 Context) — so restart / handoff files need a record extension. (Corrected
  after review: MPI migration is *not* capped — `serialize()` packs `nstate_vars`
  dynamically, `particle.tcc:1321-1330`.) (2) MPM runs its **own** loop: no fork
  integrators, algorithms, solvers, recorders, profiler, monitor; no FE and MPM in
  one `Domain`, so ADR-86 §12.4's continuity claim is not delivered. (3) Build:
  C++17 with GNU / Intel / Clang flag branches only, Boost filesystem and Eigen
  required (`CMakeLists.txt:17-68`; CI on CircleCI, no Windows) — no MSVC path,
  while the fork's single build is Windows-first through `build.bat`. (4) Drift: a
  code pushed two days ago, vendored forever.
- **Verdict:** best capability per week, worst fit per year. **Rejected as the
  implementation; retained as a cross-code oracle (§9).**

### E — Co-simulation / handoff only

Excluded by ADR-86 D2 (the owner chose an in-tree lane). Retained only as the
fallback route to an external runout engine.

### F — Self-stepping explicit engine; implicit via A or C later (added by the review)

The architecture of Berkeley, CB-Geo and Kratos-explicit, with OpenSees supplying
only the constitutive half: an engine-owned grid (plain arrays, active-node
indexing), lumped mass with a mass cutoff, USL / MUSL updates and its own particle
loop — calling OpenSees `NDMaterial`s and `ladruno_hypo::hypoIncrement` at
particles, but **no OpenSees integrator, SOE or handler** for explicit steps.
Implicit MPM arrives later via A or C, with §6-0.

- **Wins:** sidesteps W1 for the regime that matters most (runout is explicit) and
  W3–W5 entirely; no per-step rebuild; the engine threads its own particle loop
  (still under the WP-107 allowlist rule — so, today, serial for soils).
- **Costs:** two code paths (the explicit engine and an integrator-based implicit
  lane); explicit FE–MPM coupling becomes staggered co-stepping with
  CentralDifferenceLadruno plus a seam penalty — Kratos's partitioned pattern; the
  explicit lane gives up the fork's explicit stack (mass scaling, `criticalTimeStep`,
  energy balance) unless re-plumbed.

### C′ — Handler-emitted MPM adapters (considered, rejected on cost)

Extend `LadrunoContactHandler` with an MPM emitter and Node-less `DOF_Group`s for
the grid DOFs (precedent: `LagrangeDOF_Group`), avoiding W4 and W5. The first draft's
"MPM *must* use Domain Elements" was a false dichotomy (review §F7). Rejected on
cost, not impossibility: every re-emit is a full `domainChanged()` (ADR-60 D6), §6-0
is still required, and it inherits the contact handler's Plain-only restriction (W2).

---

## 7. Scoring — may the best idea win

Weights are explicit so they can be argued; re-weighting is the owner's call.
Scores 1–5. **Revised after review (§15):** first-draft scores in brackets where
they changed.

| Criterion (weight) | A particle-as-element | B dense cells | C sparse cells | D vendor Berkeley | E co-sim | **F** self-stepping explicit |
|---|---|---|---|---|---|---|
| Solve efficiency — sparse runout, implicit penetration (25) | 3 | 2 | 4 [5] — node-set churn or ghost mass (§6-C) | 4 | 4 | 5 |
| Fit with OpenSees contracts (20) | 4 [5] — §6-0 on non-rebuild steps; no WP-107 for soils | 3 [4] — needs §6-0 | 3 [4] — needs §6-0; no WP-107 for soils | 1 | 2 | 3 |
| FE continuity (15) | 5 | 5 | 5 | 1 | 1 | 3 — staggered co-stepping |
| Cost and risk to a first lane (15) | 4 | 4 | 2 | 3 [2] — MPI record not capped | 5 | 3 — two paths |
| Numerics ceiling (10) | 4 | 4 | 5 | 5 | 5 | 4 |
| Parallel path (10) | 2 | 3 | 2 [3] — ADR-88 D1 is a draft; the hook never fires in a Subdomain | 5 | 5 | 3 |
| Maintenance / drift (5) | 5 | 5 | 5 | 1 | 3 | 5 |
| **Weighted total (/100)** | **75** [79] | **68** [72] | **72** [83] | **57** [54] | **69** | **74** |

The review re-scored A, C and D and added F; it left B at 72. B writes P2G into
nodes exactly as A and C do, so it needs §6-0 too — its fit drops 4 → 3 here
(72 → 68). This is the author's consistency correction, not the reviewer's.

**Reading the table honestly.** Three finalists — **A 75, F 74, C 72** — sit inside
the scoring's noise. The first draft's 11-point lead for C rested on two refuted
claims (free orchestration through the Nodes, W1; threaded soil updates) and one
omitted cost (node-set churn). **No provisional winner is declared.** M0b (§10)
measures what separates them: A's rebuild cost per particle crossing, the §6-0
re-seed cost, and F's explicit throughput.

---

## 8. What we take from each library

"**Code**" = permissive, may be reused with attribution. "**Idea**" = from the
published papers and documentation only.

| From | Take | Path |
|---|---|---|
| Kratos | the framework-embedded design; P2G + reset **inside the scheme's predictor** (the lesson W1 teaches) and G2P in the commit slot; neighbour-first search; PQMPM; mixed u-p; non-conforming particle BCs by penalty (Chandra et al. 2021); partitioned FE–MPM coupling | code (+ acknowledgement), fork only |
| Berkeley | active-node re-initialisation per step; B-bar particles; LME basis; semi-implicit two-phase; frictional multi-material contact; dynamic particle serialization for migration; load-balancer design | code |
| CB-Geo | USF / USL / MUSL as interchangeable scheme classes (`mpm_scheme_usf/usl/musl`) | code |
| Uintah | structured patches + AMR as the route to scale | code |
| Karamelo | TLMPM as an option for solids without mixing | idea |
| NairnMPM | multimaterial contact and imperfect interfaces; FMPM(k) / XPIC to control PIC dissipation | idea |
| Anura3D | double-point two-phase; MPM-MIXED / UL-FEM — FE and MPM coexisting in one mesh | idea |
| AMPLE | implicit quasi-static MPM and the elastic column under self-weight | idea |
| Taichi | APIC / MLS transfers; sparse grids; particle sorting for cache locality | code |
| GeoTaichi | MPM–DEM coupling (much later) | idea |
| Matter | finite-strain granular rheologies incl. μ(I) for runout | idea |
| Sparse framework | active-node indexing — C's virtual grid is its CPU "scan" variant realized as Domain objects; F's grid would be it directly | idea |

### The apeGmsh half

apeGmsh is ours, and ADR-60a already sanctions it as a **seeder** — explicitly
including "an actual OpenSees PFEM fluid model — particle cloud …". For MPM it owns:

- **Seeding** — particles inside Gmsh volumes (n per cell, point-in-volume), or at
  the Gauss points of an FE mesh of the body (Kratos's `GenerateMaterialPointElement`
  pattern: volumes partition exactly), with material tags from physical groups.
- **The grid spec** — origin, spacing, counts, sized from the domain box and a
  declared runout envelope; boundary planes tagged for grid-aligned Dirichlet
  conditions.
- **Bulk transport** — at 10⁵–10⁶ particles per-particle commands are too slow; an
  HDF5 particle file in the LadrunoRecorder GP layout, read by one fork command.
- **The handoff (M5)** — a fresh-session rebuild in the apeGmsh ADR 0103 style:
  particles at the deformed FE Gauss points, carrying the ADR-88 D1 state bundle
  (nearest-GP copy degenerates to an identity copy when particles sit *at* the GPs).
  A fresh model means the first `analyze` runs `domainChanged()`, so the handoff
  itself does not meet W1.
- **Post** — particle clouds from the recorder HDF5 into apeGmsh Results and the
  viewers.

Licence: apeGmsh is GPL-3 — MIT, Apache, LGPL and GPL-3-compatible material may
enter it; Kratos code (BSD-4) and GPL-2-only code (Karamelo) may not (§2).

---

## 9. Verification — ADR-87 D1 / D2 / D5

Oracles, numpy in the `tests/*_reference.py` idiom, no build:

| ID | Oracle | Checks | Phase |
|---|---|---|---|
| V1 | 1D axially vibrating bar — closed form | P2G/G2P, transfer dissipation (PIC vs FLIP), energy, explicit + implicit | M0a |
| V2 | elastic column under self-weight, large deformation — closed-form vertical **Cauchy stress** σ_yy = ρ₀g(l₀ − Y), exact at large deformation and law-independent in 1D (AMPLE documentation case l₀ = 50 m, E = 1 MPa, ν = 0, ρ₀ = 80 kg/m³; Charlton et al. 2017). A *displacement* comparison is law-dependent — the hypo kernel will not coincide with AMPLE's log-strain law beyond small strain | implicit Newton; convergence order per basis (stress) | M0a |
| V3 | method of manufactured solutions for MPM (Kamojjala et al. 2015) | spatial convergence order per basis | M2 |
| V4 | cell-crossing meter — one particle crossing a cell face under uniform strain | the stress jump, linear vs GIMP vs B-spline, **as a number** — run with SANISAND particles too (§12) | M2 |
| V5 | plate with a hole — Howland stress concentration 3.14 | 2D elastic accuracy | M1 |
| V6 | granular column collapse vs Lube et al. (2004/2005) experiments | validation: runout vs aspect ratio | M3 |
| V7 | rigid footing penetration vs ADR-79 §9 (q = 1108–1152 kPa at s/B = 10 %, punching) | FE ↔ MPM on the fork's own measured anchor; ties to TIMs | M3 |
| V8 | Terzaghi consolidation on particles, reusing the LadrunoUP oracle | two-phase | M4 |

- **Cross-code (D5):** Kratos (`pip install KratosMPMApplication`, BSD-4) and
  Berkeley (MIT) run V1, V2, V5, V6 — **Tier-2 screens, never gates** (the ADR-86
  §4.2 cross-code trap applies in full).
- **Mutation gate (D2):** `-DLADRUNO_MUTATE_MPM` with two named mutants, each paired
  with the oracle that catches it: *flip the P2G momentum sign* → V1 red; *drop the
  G2P position update* → **V2** (or a large-amplitude V1) red. A small-amplitude V1
  may stay green on the second mutant, because position enters only through
  shape-function evaluation, at second order in amplitude (corrected after review).

---

## 10. Phases

Model, effort and work-package assignment for every phase: [[_adr98_orchestration_plan]]
(draft, awaiting the owner's approval — ADR-87 requires it decided at kickoff and
recorded). M0 is fully assigned there; M1–M7 are provisional until M0b picks the
architecture.

| Phase | Scope | Gate |
|---|---|---|
| **M0-pre** | **the §6-0 re-seed contract** — `reseedFromCommitted()` on `TransientIntegrator` + Newmark / HHT / GeneralizedAlpha / CentralDifference / CentralDifferenceLadruno / DynamicRelaxation, the `Domain` flag, the drivers; `LEDGER_vanilla_files` rows | byte-identical when the flag is never raised (every existing deck); a deck that writes committed nodal velocity between two steps sees it in the next predictor — the W1 falsifier, red today |
| **M0a** | numpy oracle — 1D/2D MPM, explicit MUSL (the F engine) + implicit Newton, linear / GIMP / quadratic B-spline; V1, V2, V5 references | closed forms reproduced; convergence orders measured |
| **M0b** | **measurements, not spikes** (re-scoped after review — the first draft's "spike of C" needed the hook, the re-seed, materialization, a particle store, an SP lifecycle and a recorder: that is M1): (i) **A's rebuild cost on the existing fork** — a synthetic churn deck adding / removing N elements between `analyze` calls under `system Diagonal` and PARDISO, N = 10³ … 10⁵, zero MPM code; (ii) the **re-seed cost** from M0-pre; (iii) **F's explicit engine** as the M0a oracle's C++ twin — throughput vs particle count with elastic and SANISAND particles | the three numbers written into this ADR; **the winner among A / F / C declared on them** |
| **M0c** | decisions: contact (§11 D-contact), parallel in / out, D2 re-confirmation (ADR-86 §5.5 warning) | owner sign-off |
| M1 | the winner, 2D + 3D, linear basis, explicit + implicit, any `NDMaterial` through the hypo kernel; particle recorder (LadrunoRecorder only on grid tags, W4); both interpreters; apeGmsh seeding v0; manifest + mutation gate + guide stub | V1, V2, V5 green; mutation red; both interpreters |
| M2 | GIMP + quadratic B-spline | V3 orders; V4 crossing numbers |
| M3 | locking — B-bar / F-bar particle option (Berkeley precedent); undrained | the ADR-81 discipline: plateau or no ship; V6, V7 |
| M4 | two-phase single-point u-p particles through `LadrunoUPKernel.h` | V8 |
| M5 | FE → MPM handoff — apeGmsh fresh-session rebuild + ADR-88 D1/D2 state at particles | identity round trip (the ADR-88 INV-1 analogue): FE state → particles → one zero-load MPM step, no drift |
| M6 | contact — multi-velocity-field MPM–MPM with friction; FE–MPM per M0c | per M0c |
| M7 | parallel, if IN | per M0c |

---

## 11. Decisions requested

- **D-arch.** No provisional winner. Accept **A, F and C** as finalists and the
  re-scoped M0b as the decider; accept **§6-0 as a prerequisite regardless of which
  finalist wins** (A, B and C all need it; F's implicit lane will).
- **D-contact — reframing ADR-86 §5.3, corrected after review.** FE–FE stays with
  `LadrunoContactDomain`. MPM–MPM uses the field's settled answer (S7,
  multi-velocity field with friction) — **which is a second frictional contact
  implementation in the tree**, the duplication ADR-86 §5.5 obligation 1 objected
  to; this ADR *accepts* the duplication, it does not dissolve it. FE–MPM: (i) **rigid
  or prescribed-motion bodies as grid Dirichlet / friction conditions** for v1 — the
  first draft's "shared grid, automatic no-slip" is incoherent for a *deformable* FE
  body, because grid nodes are reset every step while FE nodes carry total
  displacement; (ii) **penalty coupling of FE seam nodes to the particles / grid**
  for deformable bodies — Kratos's partitioned pattern
  (`python_scripts/apply_mpm_coupling_interface_dirichlet_condition_process.py`,
  CoSimulationApplication `mpm_dirichlet_wrapper.py`); (iii) `LadrunoContactDomain`
  with particles as slaves rejected — particles are not `Node`s.
- **D-parallel.** OUT for M1–M6: serial. WP-107 threads nothing soil-shaped today
  (§5), so "serial + threads" means serial for the models that matter until the
  ManzariDafalias threading defect is found and fixed. IN later needs particle
  migration (ADR-88 D1, still a draft), a partition mechanism (`repartition()` is a
  stub), and a hook that fires inside Subdomains (it does not; W6's dead loop too).
- **D-RITSS.** RITSS (ADR-88 + apeGmsh ADR 0103) and MPM are complementary, not
  rivals: RITSS where a remesh is cheap and history transfer is tolerable
  (quasi-static footing penetration); MPM where remeshing is impractical (flow,
  runout, fragmentation). They share ADR-88 D1 (state) and D3 (regularization).
  MPM's M5 depends on ADR-88 S1–S2 having landed; M0–M4 do not.
- **D-skill.** ADR-86 D3 approved a `material-point-method` skill. §2–§8 of this ADR
  are its seed; write it alongside M0, as ADR-86 recommended.

---

## 12. Risks

- **The value proposition itself.** Sordo, Rathje & Kumar had PM4Sand at the FEM
  end (OpenSees) and deliberately used total-stress Mohr–Coulomb with strain
  softening for the MPM runout, citing MPM's stress inaccuracy (arXiv 2412.08040, per
  the review). If MPM stress noise makes cyclic sand models unusable at particles,
  the fork's differentiator is theoretical. V4 must be run with SANISAND particles,
  not only elastic ones, and V6 / V7 must use the soil models — early, before M3's
  investment is committed.
- **The scoring is a judgment.** Its weights are explicit so they can be argued;
  the bake-off exists because A, F and C are within its noise.
- **Unmeasured claims (ADR-87):** rebuild cost, re-seed cost, memory per particle,
  explicit throughput — all measured at M0b, none asserted.
- **WP-107's ManzariDafalias defect** (4/4 segfaults in the plastic branch, root
  cause not located) is a prerequisite for any threaded SANISAND-family MPM, in any
  architecture.
- **Vanilla footprint grows.** §6-0 touches `TransientIntegrator.h`, four stock
  integrators and the drivers; the hook touches `Domain.cpp/.h` — each a
  `LEDGER_vanilla_files` row and a `// Ladruno` marker.
- **Kratos BSD-4:** any Kratos code in the fork carries the advertising
  acknowledgement; prefer re-derivation to keep the fork's licence surface simple.
- **Small-mass nodes:** the classic explicit MPM instability at nodes touched by a
  particle's support edge; the fork's generic explicit integrator has no mass cutoff,
  so the cell / engine needs one (F owns this directly).
- **ADR-86 §5.5 obligation 2, restated:** native MPM is justified by *fork materials
  at particles* (at the strength §0 states), *fork algorithms and solvers*, and *FE
  continuity in one `Domain`* — never by out-running Berkeley on throughput. The
  claim is perishable; re-test it at every gate.
- **Review status.** Reviewed once (§15), not clean, repaired here. The revision
  changed the plan materially (a prerequisite, a new finalist, a re-scoped M0), so
  it should itself be re-reviewed before M0b is scheduled.

---

## 13. References

**Libraries (licences as verified in §2)**
- Kratos MPMApplication — <https://github.com/KratosMultiphysics/Kratos/tree/master/applications/MPMApplication>; local checkout `C:\Users\nmb\Documents\Github\Kratos`
- geomechanics/mpm (Berkeley) — <https://github.com/geomechanics/mpm>
- CB-Geo mpm — <https://github.com/cb-geo/mpm>
- Uintah — <https://github.com/Uintah/Uintah>
- Karamelo — <https://github.com/adevaucorbeil/karamelo>
- NairnMPM — <https://github.com/nairnj/nairn-mpm-fea>; <https://osupdocs.forestry.oregonstate.edu/index.php/NairnMPM>
- Anura3D — <https://github.com/Anura3D/Anura3D_OpenSource>
- AMPLE — <https://github.com/WMCoombs/AMPLE>
- Taichi — <https://github.com/taichi-dev/taichi>
- GeoTaichi — <https://github.com/Yihao-Shi/GeoTaichi>
- Matter — <https://github.com/larsblatny/matter>
- Sparse framework — arXiv:2605.28525

**Method**
- Sulsky, Chen, Schreyer (1994), CMAME 118:179–196 — the original MPM.
- Bardenhagen, Brackbill, Sulsky (2000), CMAME 187:529–541 — multi-velocity-field contact.
- Guilkey, Weiss (2003), IJNME 57:1323–1338 — implicit MPM.
- Sulsky, Kaul (2004), CMAME 193:1137–1170 — implicit dynamics in MPM.
- Bardenhagen, Kober (2004), CMES 5(6):477–496 — GIMP.
- Steffen, Kirby, Berzins (2008), IJNME 76:922–948 — quadrature errors, B-splines.
- Sadeghirad, Brannon, Burghardt (2011), IJNME 86:1435–1456 — CPDI.
- Nairn (2013), CMES 92(3):271–299 — imperfect interfaces via multimaterial MPM.
- Jiang et al. (2015), ACM TOG 34(4) — APIC.
- Kamojjala, Brannon, Sadeghirad, Guilkey (2015), Eng. Comput. 31(2):193–213 — verification tests incl. MMS.
- Charlton, Coombs, Augarde (2017), Comput. Struct. 190:108–125 — iGIMP.
- Hammerquist, Nairn (2017), CMAME 318:724–738 — XPIC.
- Coombs et al. (2018), CMAME 333:1–21 — volumetric locking in MPM.
- Hu et al. (2018), ACM TOG 37(4) — MLS-MPM.
- Fern, Rohe, Soga, Alonso (eds.) (2019), *The Material Point Method for Geotechnical Engineering*, CRC.
- Iaconeta, Larese, Rossi, Oñate (2019), Comput. Mech. 63:1243–1260 — stabilized mixed implicit MPM (Kratos).
- Kumar et al. (2019), arXiv:1909.13380 — CB-Geo scalability (shared-memory 3.28× on 16 threads; distributed tests "required").
- Coombs, Augarde (2020), Adv. Eng. Softw. 139:102748 — AMPLE.
- Giridharan, Gowda, Stolle, Moormann (2020), *Soils and Foundations* 60(4):989–1000, doi:10.1016/j.sandf.2020.06.001 — UBCSAND and hypoplasticity in CPDI two-phase MPM.
- de Vaucorbeil, Nguyen, Hutchinson (2020), CMAME 360:112783 — TLMPM.
- de Vaucorbeil et al. (2020), Adv. Appl. Mech. 53:185–398 — MPM after 25 years.
- Chandra et al. (2021), Acta Geotech. 16:2315–2335 — non-conforming Dirichlet BCs (Kratos).
- de Vaucorbeil et al. (2021), Comput. Part. Mech. 8:767–789 — Karamelo.
- Solowski et al. (2021), Adv. Appl. Mech. 54:113–204 — MPM overview and challenges.
- Wilson, Wüchner, Fernando (2021), IJNME 122:1513–1537 — PQMPM.
- Tjung, Soga, arXiv:2111.13584 — liquefaction-induced dam failure in CB-Geo (Mohr–Coulomb, residual strengths).
- Sordo, Rathje, Kumar, arXiv:2412.08040 — hybrid FEM (PM4Sand) → MPM (Mohr–Coulomb with softening) for tailings-dam failure.
- Blatny et al. (2025), GMD 18:9149–9166 — Matter.
- Lube et al. (2004), JFM 508:175–199; (2005), Phys. Rev. E 72:041301 — granular column collapse.
- Hu, Randolph (1998), IJNAMG 22:327–350 — RITSS.

**Fork**
- [[86_geomechanics_libraries_and_mpm_scoping_report]], [[88_ladruno_state_transfer_adaptivity_adr]] (draft, uncommitted),
  [[87_ladruno_depth_with_width_adr]], [[60a_pfem_crosspollination_amendment]],
  [[79_ladruno_geom_hypo_adr]], [[71_ladruno_up_family_adr]], [[107_ladruno_openmp_element_loop]],
  [[75b_ladruno_threaded_assembly_adr]], [[39_ladruno_contact_domain_adr]]

---

## 14. Implementation log

*(empty — no code authored under this ADR)*

---

## 15. Adversarial review log

### 2026-09-25 — Fable, post-draft gate. **NOT CLEAN — one blocker, five majors.**

Pointed at eleven claims, instructed to verify against source and never to open
GPL/LGPL source. Findings and disposition:

| # | Finding | Verdict | Disposition |
|---|---|---|---|
| F1 | **Every transient integrator caches its nodal vectors and re-reads the Nodes only in `domainChanged()`** — P2G written into Nodes is discarded at `newStep()`; "works with any integrator" is false | **BLOCKER — REFUTED** | W1 added (§5); §6-0 re-seed contract made a prerequisite; §6-C orchestration rewritten; M0-pre added |
| F2 | Kratos-explicit keeps no DOF set; A's "Kratos pays this too" precedent does not exist for explicit | MAJOR — WEAKENED | S2 and §6-A corrected |
| F3 | FE–MPM "shared grid, automatic no-slip" is incoherent for a deformable FE body; S7 is a second contact engine | MAJOR — REFUTED | §11 D-contact rewritten |
| F4 | Differentiator overstated: UBCSAND has run at particles; Anura3D's external-soil-model route admits user DLLs; Sordo et al. chose not to carry PM4Sand into MPM | MAJOR — WEAKENED | §0 restated; value-proposition risk added (§12) |
| F5 | "Rides WP-107" is false — the allowlist holds no soil model | MAJOR — REFUTED | §5, §6, §11 D-parallel corrected |
| F6 | Materialize / retire omitted node-set churn, recorder pointer lifetime and removal cleanup | MAJOR — WEAKENED | W3–W5 added; §6-C costs rewritten |
| F7 | The handler wall is a false dichotomy; the real wall is the contact handler being Plain-only | MINOR | W2 rephrased; C′ added, rejected on cost |
| F8 | Uintah is MIT; Berkeley's `svars[20]` caps the HDF5 checkpoint record, not MPI migration; Karamelo GPL-2-only and dormant; NairnMPM has no licence; PQMPM is explicit-only | MINOR | §2, §3, §6-D corrected |
| F9 | The scoring flips inside the noise; F (self-stepping explicit) was never scored | MAJOR | §7 re-scored (A 75, F 74, C 72); F added (§6) |
| F10 | V2's closed form is the Cauchy stress, not a displacement; the G2P mutant escapes a small-amplitude V1 | MINOR | §9 corrected |
| F11 | GPL-2-only is incompatible with GPL-3 apeGmsh; the hook is overridden in parallel; M0b was M1 in disguise | MINOR | §2, §5, §10 corrected |

**Reviewer's verdict:** "Not safe to act on as written" — C no longer beats A; only
"inside the noise, M0b decides" survives. **The one thing to change before M0:**
replace "orchestration through `Domain::analysisStep` works with any integrator" with
an explicit integrator-side contract (§6-0).

**Author's verification before acceptance** — not taken on the reviewer's word:
F1 (`Newmark.cpp:206-208` predicts from its own `U, Udot, Udotdot`; `:398-406`
re-reads nodes only in `domainChanged`; `CentralDifferenceLadruno.cpp:360, :518`);
F2 (`mpm_explicit_strategy.hpp:197-218`: `dof_set_dummy`, scheme-only `Clear()`);
F5 (`ManzariDafalias.cpp:1094` refuses unconditionally; `LadrunoSANISAND.cpp:1948-1965`
delegates; only `ElasticIsotropicPlaneStrain2D.h:96` returns true); F6 / W4
(`Domain.cpp` calls recorders' `record / restart / flush` only — never
`domainChanged`); F8 (Uintah `src/LICENSE` is MIT); F11 (`PartitionedDomain.cpp:1131`,
`Subdomain.cpp:808` override `analysisStep`; the dead subdomain loop at
`PartitionedDomain.cpp:1147-1156`, recorded as W6). **All confirmed.**

**One correction to the review:** it cited the UBCSAND-in-MPM paper as *Computers
and Geotechnics* 2020. Crossref: Giridharan, Gowda, Stolle, Moormann, *Soils and
Foundations* 60(4):989–1000 (2020), doi:10.1016/j.sandf.2020.06.001.

**The failure mode, again.** F1 was already known to this fork: ADR-88 D2 records
that "Newmark reads committed u, v, a in `Newmark::domainChanged`"
(`Newmark.cpp:390-418`) — the author had read it in the same session and did not
carry it into ADR-98. Same class as ADR-86's blocker: a conclusion drawn without
re-reading one's own neighbouring documents.

**Owed ledger rows** (land with this ADR's PR): `LEDGER_quirks` — W1 (integrators
discard nodal state written between steps), W4 (stock recorders cache pointers;
`Domain` never calls `Recorder::domainChanged()`), W6 (the dead subdomain
`domainChange()` loop, an upstream bug).
