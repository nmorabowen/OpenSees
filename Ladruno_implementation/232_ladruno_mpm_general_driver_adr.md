---
title: "WP-232 — A general MPM driver for OpenSees: architecture, contracts, apeGmsh integration, and the first slice"
project: Ladruno
status: draft — NO code. Supersedes the unmerged ADR-98 draft. Four adversarial reviews folded (§14.1–§14.4)
priority: high
owner: nmora
wp: "WP-232 (provisional — confirmed when the MPM program issue claims it)"
amends: 86_geomechanics_libraries_and_mpm_scoping_report
tags:
  - adr
  - mpm
  - material-point-method
  - large-deformation
  - apegmsh
  - driver
aliases:
  - WP-232
  - MPM general driver
updated: 2026-11-13
---

# WP-232 — A general MPM driver for OpenSees

> **Lineage.** This document supersedes the ADR-98 MPM scoping draft, which was never merged
> (rescued as commit `a1c0b4a14` on `origin/rescue/adr98-mpm-scoping`); the number 98 has since
> been taken on `ladruno` by the membrane load-distribution ADR. It keeps that draft's library
> survey, licence map and settled ideas, and replaces its architecture contest, which two
> adversarial reviews dismantled (§14.1, §14.2). A third review rewrote this document's apeGmsh
> half and its footing gate (§14.3).
>
> **Number.** WP-232 is the next number above every WP in use on 2026-11-13 (highest: 231, from
> branches, PR titles and ledger fragments). It is provisional until the MPM `program` issue
> claims it in its title, as #1024 did with "(WP-211…)". That is practice, not rule:
> `Ladruno_internal/ISSUES.md` describes program issues but says nothing about reserving numbers.
>
> **Goal, set by the owner.** A **general MPM driver** for OpenSees — any material the driver admits,
> any problem class (runout, flows, impact, penetration, quasi-static large deformation). The
> TIMs shallow footing is one test in the suite, not the purpose.
>
> **Amends ADR-86 D2.** D2 reads "A **and** B": the MPM lane built in-tree, *and* driven from the
> FEM initiation lane through the sequential handoff (`86_…report.md:500-501`). This document keeps
> A and moves B — the FE → MPM handoff — from the v1 driver to phase 9 (§9), because the
> state-transfer primitive it needs has not landed (§1). The owner approved this on 2026-11-13
> (§10, D7), and ADR-86 D2 carries the amendment note.
>
> **Discipline.** The three documents before this one each shipped a claim refuted by something the
> fork already said about itself (an implementation log, an integrator's cached state, a
> regularization report, an ADR's own limits section). This document's first draft shipped five
> more (§14.3): two contradicted apeGmsh's own current ADRs, one inverted ADR-79's own sentence, one
> misread the precedent's header comment, and one cited an apeGmsh ADR number that means something
> else. The rule is therefore stricter: every reused artifact is cited next to its own scope or
> limits, **quoted, not recalled** (§4); fork claims are checked on `origin/ladruno` @ `aa6229bf9`,
> apeGmsh claims on apeGmsh `origin/main` @ `cf2e06593`; every ADR number cited was opened on the
> branch it is cited from. Under ADR-87, an unmeasured claim is a claim, and says so.

---

## 0. Decision summary

| Question | Answer |
|---|---|
| **Architecture** | MPM is a **driver**: an `mpm …` command family that owns the step loop — one row of `LadrunoCommandTable.h`, feeding all three interpreter engines. The precedent is the fork's `LadrunoStaggeredAnalyze`, a command that decomposes `analyzeStep` (§5.1). Two solve modes. |
| **Explicit mode** | **F — a native engine**: engine-owned grid arrays with active-node indexing, lumped mass with a cutoff, MUSL, PIC/FLIP. No OpenSees integrator, SOE or handler, so no per-step rebuild and no integrator-state problem (W1). The engine object is **hosted as a `LoadPattern` subclass**, after ADR-73's `LadrunoPorousOverlay`: `wipe` deletes it, the `parameter … pattern` target reaches it, and no `FE_Element` is ever built for it (W15). It runs G2P itself before calling `Domain::commit`, so it needs **no Domain hook** (§5.2). |
| **Implicit and quasi-static mode** | **A — particle-as-element**, connectivity re-pointed and `Domain::domainChange()` raised **unconditionally every step** — Kratos and Berkeley both rebuild every step. Transient integrators are restricted by class tag to Newmark, HHT and Generalized-α, whose own `domainChanged()` re-reads committed disp, vel and accel — **no new integrator hook**. Quasi-static calls `StaticAnalysis::analyze(1, false)`, because `StaticAnalysis` has no `analyzeStep` (W16, §5.3). |
| **Deferred / rejected** | C (sparse cell kernels) deferred as an optimization, only if M0b's numbers demand it. Rejected: B (dense grid), C′ (handler adapters), D (vendor Berkeley), E (co-simulation only), a zero-node Domain element as the explicit host (W15), and a Domain hook for mode F's commit (§5.5). |
| **Particle and material contract** | Each particle carries its own F, rotation and strain feed; a **new incremental entry** to the hypo kernel (its numerics reused, its entry point not); an **admitted-material table** with loud refusal for the rest — v1 admits `ElasticIsotropic`, `LadrunoJ2` and the post-`31322a47a` `DruckerPrager`. The engine zeroes `ops_TheActiveElement` around particle updates, and v1 offers no regularization at particles (§6). |
| **apeGmsh** | Rewritten against apeGmsh's current architecture (its ADRs 0112, 0116, 0121, 0122, 0131). **Particle results are a new group inside `.ladruno`** — an additive change: `FORMAT_MINOR_VERSION` 1 → 2, `FORMAT_VERSION` stays 1, so today's apeGmsh keeps opening the files — read by apeGmsh's `.ladruno` readers in Python and TypeScript and drawn as a three.js points layer in apeGmshViewer — not a new results file, and not the frozen Qt viewers. The particle set is a field-by-field contract tied to apeGmsh's emit (charter P4: the bridge owns tags). All of it under a new apeGmsh ADR that the maintainer accepts. Bridge primitives only after the fork slice merges **and** a fork release is pinned in apeGmsh's `FORK_PIN` (§7). |
| **Seeding** | **Lattice first**: n particles per grid cell per axis, kept where they fall inside a body — positive volumes, uniform particles per cell, GIMP half-lengths by construction. Gauss-point seeding later, with positive-weight rules only: apeGmsh's tetrahedron rule has a negative weight (§7.4). |
| **Instant Meshes** | Studied at the owner's request: a BSD-3-style field-aligned **surface** remesher with a batch mode and point-cloud input. Useful later for MPM → FE remeshing; not for seeding or the grid. Adopt only if it beats Gmsh's own discrete-surface remeshing on the same surfaces, passes a closure check, and comes through a binding whose licence is verified (§7.5). |
| **First slice** | Explicit, 2D plane strain and 3D, serial; linear + GIMP; MUSL + PIC/FLIP; three admitted materials; one stage path (`DruckerPrager` `materialState`); gravity ramp, grid planes, frictional wall, quiet K0 start, local damping for quasi-static checks, a critical-Δt check; particle output into `.ladruno` — **with the apeGmsh ADR, lattice seeding, the readers and the points layer in the same slice** (§9). |
| **Differentiator, at honest strength** | No surveyed open-source MPM code ships a cyclic-liquefaction sand model in its own library; PM4Sand / SANISAND / PDMY at material points was not found in the published record. NorSand at particles already exists (Berkeley ships it), so `LadrunoNORSAND` serves as a matched-material cross-code oracle, not a differentiator. |

---

## 1. Goal and non-goals

**Goal.** A general MPM capability in the fork: one driver, explicit and implicit, any admitted
material, general loads and boundaries, several bodies, optional two-phase particles, and a
first-class front and back end in apeGmsh. FE and MPM in one model is a later phase (§9).

**Non-goals for v1** (each is a later phase, §9):

- **Regularization for softening soils.** No soil model in the fork consumes a length (W9), and the
  only soil length-scale work in existence (the TIMs SANISAND band scaling, unpushed) fails several
  of its own pre-registered gates. ADR-90's report found that rate regularization gives no band
  width in quasi-statics. Softening soils are open research, not a v1 feature.
- **An explicit FE–MPM seam.** v1 couples MPM only to rigid or prescribed-motion bodies expressed as
  grid conditions, and explicit mode **refuses a Domain with free FE DOFs**. A deformable FE–MPM
  seam comes later, implicit first (§9).
- **Parallel.** Serial in v1; the engine threads its own particle loop only for materials on the
  WP-107 allowlist, which today contains no soil model (W14).
- **FE → MPM handoff.** It needs a material state-transfer primitive (`getHistory` / `setHistory`).
  That primitive exists only as a draft: an uncommitted ADR in the owner's main checkout
  (`88_ladruno_state_transfer_adaptivity_adr.md`, 2026-09-02 — its number is taken on `ladruno` by
  the merged H5DRM higher-order-elements work, #778 / #800), and phases M1b and M6 of ADR-86 §6.3.
  Nothing has landed.
- **Arc-length on a reset grid.** Unprecedented in either reference code; quasi-static uses load
  control and rigid-body displacement control.
- **Fluids.** A dam-break validation needs a fluid law the fork does not obviously have; it is
  listed only if one is admitted.
- **Staged live runs from apeGmsh.** apeGmsh's live path runs flat (unstaged) models only (its ADR
  0121, ruling 2), so a staged K0 → plastic MPM run is deck-only.
- **Database restart of particle state.** Refused while particles exist; it arrives with the restart
  test in phase 1.

---

## 2. The field — twelve codes, and the licence map

Licences verified from each repository's own licence file, the GitHub licence API, or Crossref.
**Source** = may be read and reused with attribution. **Papers** = ideas from publications and
documentation only; never open the source while writing code.

| Code | Licence | May use | What it teaches this driver |
|---|---|---|---|
| Kratos MPMApplication (CIMNE) | BSD-4 (advertising clause) | source, **fork only** | MPM inside a general FE framework; P2G + reset inside the scheme's predictor, G2P in `FinalizeSolutionStep`; implicit rebuilds the DOF set every step; explicit keeps **no** DOF set; PQMPM (explicit only); mixed u-p; partitioned FE–MPM coupling by penalty / Lagrange |
| geomechanics/mpm (Berkeley) | MIT | source | explicit USF/USL/MUSL, implicit Newmark, semi-implicit two-phase; GIMP, B-spline, LME; B-bar particles; frictional multi-material contact; MPI with load balancing; **ships NorSand** |
| CB-Geo mpm | MIT | source | USF / USL / MUSL as interchangeable scheme classes; dormant since 2024-11 |
| Uintah | MIT | source | structured patches + AMR as the route to scale |
| Taichi MLS-MPM | Apache-2.0 | source | APIC / MLS transfers, sparse grids, particle sorting |
| Karamelo | GPL-2 (no "or later") | papers | TLMPM |
| NairnMPM | none found (all rights reserved) | papers | multimaterial contact, imperfect interfaces, FMPM / XPIC |
| Anura3D | LGPL-3 | papers | double-point two-phase; MPM-MIXED / UL-FEM |
| AMPLE | GPL-3 | papers | implicit quasi-static MPM; the self-weight column |
| GeoTaichi | GPL-3 | papers | MPM–DEM coupling |
| Matter | GPL-3 | papers | finite-strain granular rheologies, μ(I) |
| Sparse framework (arXiv 2605.28525) | not stated | idea | active-node indexing; 1–2 orders of magnitude on sparse problems |

**Tools on the apeGmsh side.** apeGmsh itself is GPL-3. Instant Meshes (Jakob, Tarini, Panozzo,
Sorkine-Hornung) and QuadriFlow (Huang et al.) carry BSD-3-style licence texts in their
`LICENSE.txt` files (the GitHub licence API reports NOASSERTION for both). The Python binding
`pynanoinstantmeshes` has **no licence file**: the GitHub API returns none, and PyPI shows
`license: None` with only the classifier "BSD License". By this section's own standard its licence
is unverified.

**Consequences.** Kratos code may enter the fork with the acknowledgement its BSD-4 requires, but
**never apeGmsh** (the advertising clause is GPL-incompatible). GPL-2-only code (Karamelo, as its
licence reads) is incompatible with GPL-3 and enters neither. BSD-3 and MIT material may enter both.
Papers-only for GPL / LGPL code in the fork is policy: the fork's `COPYRIGHT` file never mentions
the GPL, so barring GPL / LGPL source is an inference from its University of California
non-commercial terms, not a sentence it contains. This is an engineer's reading, not legal advice.

---

## 3. Settled — what every serious code agrees on

Adopted as constraints, not re-litigated.

- **S1** Particles own all history; the grid owns none and is reset every step.
- **S2** Solve over active nodes only — implicit codes rebuild the system every step; explicit codes
  keep no global system at all.
- **S3** Implicit MPM is Newton on grid increments with particles frozen at step-start positions,
  G2P at convergence (Guilkey–Weiss kinematics); quasi-static is the implicit loop minus inertia
  (Kratos `MPMQuasiStaticSolver`, Berkeley's `quasi_static_` flag).
- **S4** Linear hat functions are not production grade (cell crossing); GIMP and/or B-spline is the
  minimum.
- **S5** Near-incompressibility needs a remedy (Kratos mixed u-p; Berkeley B-bar particles).
- **S6** Locate particles neighbour-first, bins as fallback, or by index arithmetic on a structured grid.
- **S7** Body-to-body contact is a multi-velocity field with Coulomb friction.

---

## 4. The fork — assets and walls, each with its own scope

Verified on `origin/ladruno` @ `aa6229bf9` (2026-11-13). Every reused artifact is listed with the
limit that bounds its reuse, in its own words where it has them.

| Asset | Reused for | Its own scope / limit |
|---|---|---|
| `LadrunoHypoKernel.h` — `polar3`, Hughes–Winget increment, `pushStress6`, `pushTangent6` | the objective large-strain update at particles | **Entry point needs adding:** `hypoIncrement(nNodes, dNdX, un, u1, …)` builds F from reference gradients and committed / trial **nodal** displacements (`:43`, `:166-176`); a reset grid has neither. The third review reads the existing algebra as G = (F₁ − Fₙ)·F_mid⁻¹ (`:166-200`), so an entry on (Fₙ, F₁) is a refactor of it; the oracle decides. The `-geom hypo` lane is **3D-only** (`OPS_LadrunoUP.cpp:511-514`) and refuses B-bar (ADR-79). Its gated material evidence is PDMY and Drucker–Prager (ADR-79 §8–9) |
| `LadrunoUPKernel.h` — `addQ/addH/addS/addHtilde/addFseep` over `(N, dN, dv)` | two-phase particles | equal-order stabilized u-p inflates the first undrained settlement **1.71×** in TIMs' C0-up runs (stabilization off: 0.94, but an early checkerboard of 1.35 p₀) |
| `Newmark` / `HHT` / `GeneralizedAlpha` `domainChanged()` | re-seeding the implicit integrator after a per-step rebuild | each re-reads committed disp, vel and accel per DOF group (`Newmark.cpp:398-414`, `HHT.cpp:329-345`, `GeneralizedAlpha.cpp:337-353`) — but only when the domain stamp moved (`DirectIntegrationAnalysis.cpp:207-217`). `LadrunoHHT` (33013) and `LadrunoGeneralizedAlpha` (33014) are subclasses with their own tags and inherit it (no override in their files). **Ill-defined** for CentralDifferenceLadruno (a re-seed forces a starter step, `firstStep = true` at `:362`, `:523`), the SMS integrators, LadrunoDynamicRelaxation and the multi-step methods (TRBDF2, Houbolt, …) |
| `LadrunoStaggeredAnalyze` (`SRC/domain/pattern/ladrunoPorousOverlay/LadrunoStaggeredDriver.{h,cpp}`) | the driver pattern: a command that composes the user's analysis | its header's own scope: it "**decomposes** `analyzeStep`" and "omits analyzeStep's `_RELIABILITY` sensitivity pass and the numSubLevels dt subdivision" (`.h:44-48`). It is **transient-only** — an active static analysis is a hard error (`OpenSeesCommands.cpp:4083-4092`). Its Newmark / HHT / Generalized-α whitelist exists because it repeats the same step through `revertToLastStep` + `newStep` (`.cpp:152-172`), detected by `getClassTag()`; it says nothing about re-seeding. One row of `LadrunoCommandTable.h` (`:101`) |
| `LadrunoPorousOverlay` (ADR-73) | the host for the explicit engine | a `LoadPattern` subclass. `Domain::applyLoad` calls every pattern's `applyLoad` (`Domain.cpp:2148`); `wipe` deletes every pattern (`Domain::clearAll`, `:1114`); the `parameter … loadPattern <tag>` target reaches a pattern's `setParameter` in Python (`OpenSeesParameterCommands.cpp:156-176`) and in classic Tcl through its generic branch (`TclParameterCommands.cpp:235`, `:248`) — classic Tcl's `pattern` keyword means only `pattern <tag> lambda` (`:134-135`), so the spelling is `loadPattern`. `Domain::commit` drives the overlay by its class tag (`:2386-2389`), a dispatch the MPM host does not need. `revertToLastCommit` and `revertToStart` reach a pattern only as a re-`applyLoad` (`:2458`, `:2503`) |
| WP-99 commit-refusal counter | detecting a particle material's refusal at commit | a process-global atomic (`LadrunoMaterialStatus.h:162-165`), read and cleared only inside `Domain::commit`, after the element loop and before recorders (`Domain.cpp:2343-2345`); nothing resets it at the start of a commit. A refusal is seen after sibling points have committed, so it cannot be undone: `Domain::commit` calls the state "INCONSISTENT" and aborts the analysis (W17). Only `LadrunoConcrete3D`, `LadrunoHS` and `LadrunoSANISAND` declare one |
| `Domain::commit` order | where G2P must sit | nodes → elements → refusal check → contact → overlay → committed time → recorders (`Domain.cpp:2294-2423`); element `commitState()` return codes are dropped |
| The `nDMaterial` library | particle constitutive updates | admitted by table, not wholesale (§6); stage switches reach materials only through Domain elements (W12) |
| `LadrunoRecorder` and the `.ladruno` format | particle output (Contract 2) | the fork owns the format versions: `FORMAT_VERSION = 1` (`LadrunoRecorder.cpp:989`) and `FORMAT_MINOR_VERSION = 1` for additive changes (`:990-994`); `INFO/PROVENANCE` (`:1033`). apeGmsh treats a major bump as breaking and opens only `FORMAT_VERSION` 1 (its `schema/_versions.py:73-90`). A new model stage opens when the topology fingerprint of tags **and object addresses** changes, on a user `domainChange` command, or on a static↔transient flip (`:664-694`) — not on a C++ `Domain::domainChange()` that leaves the fingerprint alone. Mode-A re-pointing keeps the same objects, so element output would keep stale connectivity; per-step node churn would write one model stage per step. Particle output is therefore its own group |
| `LadrunoHS` (WP-210, from 2026-10-24) | a later particle material | young: fix rounds ran to 2026-11-09 (#1069) |
| `LadrunoNORSAND` (WP-144, draft PR #896) | the matched-material cross-code track with Berkeley | not merged; Berkeley's NorSand formulation must be aligned on an element test first |
| TIMs C0-up (`tims-workbench`, 2026-11-08) | Terzaghi and Mandel closed forms for two-phase checks | LadrunoUP Q4 passes Mandel (peak 1.0754 vs 1.0757) |

| # | Wall | Consequence for this driver | Evidence |
|---|---|---|---|
| W1 | Transient integrators cache U, U̇, Ü and re-read the Nodes only in `domainChanged()` | node writes between steps are discarded — so explicit mode does not use an OpenSees integrator, and implicit mode rebuilds every step | `Newmark.cpp:206-208`, `:398-414`; `CentralDifferenceLadruno.cpp:360`, `:518` |
| W2 | One `ConstraintHandler` per model; the contact handler is Plain-only | handler-emitted MPM adapters (C′) inherit the Plain-only restriction | `LadrunoContactHandler.h:22-29` |
| W3 | A `DOF_Group` for every Domain node | empty-support grid nodes are singular unless excluded or ghost-regularized | `PlainHandler.cpp:123-132` |
| W4 | Stock recorders cache `Node*` / `Element*`; `Domain` never calls `Recorder::domainChanged()`. `LadrunoRecorder` opens a new stage only on a topology-fingerprint change, a user `domainChange` command or a static↔transient flip | stock recorders must not watch grid entities; no recorder follows an element whose nodes change in place; the engine raises `Domain::domainChange()` from C++, never through the user command, or `recorder ladruno` opens a stage per step | `Domain.cpp` calls `record / restart / flush` only; `LadrunoRecorder.cpp:664-694` |
| W5 | `removeElement` / `removeNode` clean nothing (loads, SP/MP, recorder pointers) | gravity is a particle body force through P2G by design; never reuse a retired tag | `Domain.cpp:1198-1241`; `OpenSeesMiscCommands.cpp:269-275` |
| W6 | `PartitionedDomain::repartition()` is a stub, and its subdomain `domainChange()` loop is dead code | no parallel MPM without new infrastructure; an upstream bug | `PartitionedDomain.cpp:1147-1156` |
| W7 | PFEM `BackgroundMesh` | precedent that the Domain tolerates churn — not reusable code, and not proof that recorders or constraints survive it | ADR-86 §5.1 |
| W8 | The hypo kernel's entry is reference-configuration and 3D-only | a new incremental entry and a plane-strain embedding are new code, each with an oracle | see the asset row above |
| W9 | Materials get an element length only by reading the thread-local global `ops_TheActiveElement`, which the `Element` constructor, `Domain::addElement` and the element-update loops all write | once any element has existed the global is non-null, and after `wipe` it can dangle; consumers are concrete, J2-damage and steel (`LadrunoJ2.cpp:353-354` dereferences it when non-null); **no soil model** reads it. A host that is not an element must set it to 0 on every worker thread around particle updates | `Element.cpp:61`, `:131`; `Domain.cpp:493`, `:2534`, `:2776-2819`; quirk `LEGACY-ops-theactiveelement-is-a-mutable-global` |
| W10 | `DirectIntegrationAnalysis::analyze` sub-steps internally when `numSubLevels != 0` | wrapping `analyze()` gives several commits with no grid reset between them; the driver calls `analyzeStep(dt)` itself | `DirectIntegrationAnalysis.cpp:176-177`, `:280`, `:293` |
| W11 | Recorders run inside `Domain::commit`, last; element `commitState()` return codes are dropped | G2P must precede recorders, and a G2P inside `commitState()` cannot fail the step | `Domain.cpp:2294-2423` |
| W12 | `updateMaterialStage` and `parameter … material` reach materials only via `Domain::getElements()`; `MaterialStageParameter` stops at the first element that answers | mode-A particles are elements and are reached; mode-F particles only through the host's `parameter … pattern` target; `updateMaterialStage` stays blind to mode F — and none of the v1 materials answers it (`DruckerPrager.cpp:1063-1064` returns −1; its switch is `materialState`, `:1045`) | `MaterialStageParameter.cpp:60-81`; `MatParameter.cpp` (element walk) |
| W13 | LadrunoRecorder's `GLOBAL_GP_COORDS` is a mesh-bound layout | a particle dataset is a new output, on both sides | `LadrunoRecorder.cpp` (review §14.2 F13) |
| W14 | The WP-107 thread allowlist contains no soil model | `ManzariDafalias` refuses unconditionally (`:1235`); `LadrunoSANISAND` delegates (`:2758`); the allowlisted nD materials are `ElasticIsotropic` (3D, PlaneStrain2D) and `ASDConcrete3D` | allowlist on `origin/ladruno` |
| W15 | `FE_Element` exits the process on an element with no DOFs, and `PlainHandler` builds an `FE_Element` for every Domain element, active or not | no zero-node element may sit in a Domain: the first analysis set up there kills the process. The explicit host is a `LoadPattern`, never an element | `FE_Element.cpp:126-129`; `PlainHandler.cpp:291-327` (other handlers unverified) |
| W16 | `StaticAnalysis` exposes only `analyze(numSteps, flush)` — no `analyzeStep`, `getModel` or `checkDomainChange` — and has no sub-levels | quasi-static mode calls `analyze(1, false)` and flushes recorders itself | `StaticAnalysis.h:66`; `StaticAnalysis.cpp:128-236` |
| W17 | A material's commit refusal cannot be undone: `Domain::commit` sees it after sibling points committed and aborts the analysis ("INCONSISTENT … aborted") | a refusal at particles aborts the run in both modes — no retry; no v1 material declares one (only `LadrunoConcrete3D`, `LadrunoHS`, `LadrunoSANISAND` do) | `Domain.cpp:2343-2355`; `LadrunoMaterialStatus.h:162-165` |

---

## 5. Architecture — the MPM driver

### 5.1 A command family, not an analysis type

`mpm grid …`, `mpm body …`, `mpm particles -file …`, `mpm analyze …` — one `LADRUNO_COMMAND` row in
`SRC/interpreter/LadrunoCommandTable.h`, which registers it in the three engines (classic Tcl, DL
Tcl, Python); `ci/check_ladruno_commands.py` fails a command registered anywhere else. Two of the
three run in shipped binaries — classic Tcl (`OpenSees.exe`) and Python (`opensees.pyd`); the DL Tcl
engine is compiled but never instantiated (`TclInterpreter.cpp` is commented out,
`SRC/tcl/CMakeLists.txt:94`). Classic Tcl holds its analysis in different engine state from the
`cmds` singleton (the table's header comment), so mode A's `mpm analyze` needs two bridges to the
user's analysis, as `LadrunoStaggeredAnalyze` has (phase 1); mode F has no analysis and needs none.
v1 refuses `mpm` in OpenSeesSP and OpenSeesMP.

- **Mode A, transient:** composes the user's `Transient` analysis through `DirectIntegrationAnalysis`'s
  public `analyzeStep`, `checkDomainChange` and `getModel` (`DirectIntegrationAnalysis.h:71`, `:87`,
  `:92`). `VariableTransient` derives from `DirectIntegrationAnalysis`; it is detected by
  `dynamic_cast` and refused.
- **Mode A, quasi-static:** `StaticAnalysis::analyze(1, false)`, then the driver flushes recorders
  itself (W16). Static integrators hold no nodal vectors, so W1 does not arise.
- **Mode F:** creates no Analysis object and drives its own loop (§5.2).

Outside `mpm analyze` the host pattern's `applyLoad` is inert, so a later FE `analyze` in the same
Domain is unaffected by it.

### 5.2 Explicit mode — the native engine (F)

Engine-owned grid arrays with active-node indexing; lumped mass with a mass cutoff; MUSL stress
update; PIC/FLIP blend.

**The host.** The engine object is a `LoadPattern` subclass (provisional name `LadrunoMPMPattern`;
a class tag from the pattern band, ledgered), after ADR-73's `LadrunoPorousOverlay`. Without any
Domain edit, that gives:

- `wipe` deletes it with the other patterns (`Domain.cpp:1114`), and `remove loadPattern` removes it;
- the `parameter … loadPattern <tag>` target reaches its `setParameter` in both shipped engines,
  and it forwards to its particles' materials (W12);
- its own `TimeSeries` drives gravity ramps and prescribed grid motion.

What it does not give: `reset` / `revertToStart` reach a pattern only as a re-`applyLoad` (`:2503`),
so particles return to their initial state only through `mpm reset` — tested, and documented as
such. Two alternatives were rejected: a zero-node Domain element, proposed by the third review
(W15: it kills the process at the first analysis set up in that Domain), and a hook in
`Domain::commit` beside the contact engine, from this document's first draft (unnecessary — the
engine orders G2P before `Domain::commit` itself — and an unadmitted vanilla edit).

**Each step, in order.** Everything a particle carries is trial until step 7.

1. **Stability.** Δt must not exceed the engine's critical-Δt estimate — cell size over the largest
   particle wave speed (from each admitted material's current tangent) plus particle speed — times a
   safety factor; otherwise the step is refused. `mpm criticalTimeStep` reports the estimate.
2. **Loads and time.** `Domain::applyLoad(t + Δt)` sets the domain time and `ops_Dt`
   (`Domain.cpp:2170`) and calls the host's `applyLoad`, which reads its own time-series factor.
   Domain load patterns and `fix` / `sp` never reach the engine's grid (they act on Domain nodes and
   elements only); grid conditions are the engine's own (`mpm fix …`, planes, the frictional wall).
3. **Grid reset and P2G:** mass, momentum, external and internal force over active nodes only; nodes
   below the mass cutoff are dropped from the update.
4. **Grid momentum update** and grid conditions.
5. **Particle update (MUSL):** particle velocity from the updated grid (PIC/FLIP); a second
   projection of particle momentum; grid velocities; ΔF and rotation through the new kernel entry
   (§6), with shape functions evaluated at the step-start positions; `ops_TheActiveElement = 0` on
   every worker thread (W9); the constitutive update. The particle loop threads only for allowlisted
   materials (W14) and must give the same result at any thread count.
6. **Commit particle materials.** A declared refusal cannot be undone once sibling points have
   committed (W17): if `ladrunoPendingCommitRefusals() > 0`, the engine aborts the run, as
   `Domain::commit` does — before any particle moves. No v1 material declares one.
7. **G2P kinematics:** position, velocity, F and volume take their trial values. A failed step
   leaves them untouched.
8. **`Domain::commit()`:** commits FE nodes and elements (none free in v1), sets the committed time
   and commit tag, and runs recorders — after step 7, so positions and stresses come from the same
   step (W11).
9. **Profiler rows:** the `OPS_PROFILE_SCOPE("step")` bracket and the `OPS_PROFILE_STEP` row are
   replicated (they live in `analyzeStep`, `DirectIntegrationAnalysis.cpp:194`, `:273`), so the
   profiler viewers see MPM steps.

At the end of `mpm analyze`: `Domain::flushRecorders()` (it lives in `analyze()`, not in `commit`).

**v1 refuses:** a Domain with free FE DOFs, a non-admitted material, database save or restore
while particles exist, and `mpm` in OpenSeesSP / OpenSeesMP.

### 5.3 Implicit and quasi-static mode — particle-as-element (A)

Each particle is an element (`LadrunoMPMPoint`, provisional name) with one integration point, one
`NDMaterial` copy, and as nodes the support of the cell it occupies. Every step:

1. **Re-point connectivity** without remove / re-add: the element refetches its `Node*` itself, and
   on `domainChanged()` the handler rebuilds every `FE_Element` from `getExternalNodes()`.
2. **Raise `Domain::domainChange()` unconditionally.** `analyzeStep` runs `domainChanged()` only
   when the stamp moved (`DirectIntegrationAnalysis.cpp:207-217`), so on a step in which no particle
   changes cell the projected grid state would otherwise be ignored — W1 itself.
3. Let the analysis rebuild, which re-seeds the integrator through its own `domainChanged()`.

**Integrators.** Newmark, HHT and Generalized-α, detected by `getClassTag()` as the staggered driver
does. The justification is their own: each `domainChanged()` re-reads committed disp, vel and accel
(§4). The precedent's reason — a clean retry through `revertToLastStep` + `newStep` — also holds.
`LadrunoHHT` and `LadrunoGeneralizedAlpha` inherit the same re-seed but stay refused until a test
admits them; everything else is refused loudly.

**Element obligations.**
- G2P sits in `commitState()` (the Kratos `FinalizeSolutionStep` pattern), which runs before
  recorders and only on an accepted step. Because `Domain::commit` drops its return code (W11), it
  must not fail. A step that fails before commit (Newton divergence) is cut and retried by the
  driver with no G2P to undo; a material refusal at commit aborts the analysis (W17).
- Raise `Domain::domainChange()` from C++, never through the user `domainChange` command, which
  would open a `recorder ladruno` stage every step (W4).
- The element does not opt into the WP-192 threaded commit path (`ladrunoThreadSafeForm`).
- Particle output is the particle group (§7.3), never the element / GP layout (W4, W13).

**Quasi-static.** `StaticAnalysis::analyze(1, false)` (W16) under load control or displacement
control on a rigid body or FE node, never on a grid node; arc-length is out of scope. Empty-support
nodes are excluded or ghost-regularized — the choice is measured at M0b.

### 5.4 The step order, summarized

1. Mode A transient: the driver calls `analyzeStep(dt)` directly (no sub-levels, W10) and owns any
   step-size cutting. Mode A quasi-static: `StaticAnalysis::analyze(1, false)` (it has no
   sub-levels). Mode F: its own loop (§5.2).
2. P2G and the grid reset happen before the step; mode A raises `Domain::domainChange()`
   unconditionally (§5.3).
3. G2P happens before recorders: in each particle element's `commitState()` (mode A), or before the
   engine calls `Domain::commit()` (mode F). No Domain hook in either mode.
4. G2P runs once per accepted step. Mode A retries only steps that failed before commit, so there
   is nothing to undo; a refusal at commit aborts the run in both modes (W17).

### 5.5 What was not chosen

| Candidate | Status | Reason |
|---|---|---|
| C — sparse cell kernels, materialized on demand | **deferred** | its advantage is explicit throughput, which F now owns; revisit only if M0b shows the per-step rebuild of A dominates implicit runs |
| B — dense fixed grid | rejected | full-box factorization; 1–2 orders of magnitude of waste on sparse problems |
| C′ — handler-emitted adapters | rejected on cost | a full `domainChanged()` per re-emit and the Plain-only restriction (W2) |
| D — vendor Berkeley's engine | rejected as the implementation | stateless-material contract, its own solve loop, no MSVC path, permanent drift; kept as a cross-code oracle |
| E — co-simulation only | excluded | ADR-86 D2 chose an in-tree lane |
| A zero-node Domain element as the explicit host (§14.3 M6) | rejected | W15: the first analysis set up in that Domain exits the process |
| A `Domain::commit` hook for mode F (first draft) | rejected | unnecessary — the engine orders G2P before `Domain::commit` — and an unadmitted vanilla edit |

---

## 6. The particle kinematics and material contract

**Particle state:** x, V, m, v, a (implicit), F, R, the accumulated strain feed, σ, GIMP
half-lengths, one owned `NDMaterial` copy, a stable particle id, body id, material tag, flags.

**Kinematics.** A new incremental entry to the hypo kernel takes the particle-carried Fₙ and the
step's ΔF, and returns the strain increment and rotation in the kernel's unrotated frame. It reuses
`polar3`, the Hughes–Winget midpoint and `pushStress6`; it gets its own numpy oracle before any C++.
The plane-strain embedding is new work (the shipped lane is 3D-only).

**Admitted materials.** An explicit table; everything not in it is refused loudly. Cells marked
"232c" are filled and tested in that WP.

| Material | Type strings | Hypo use | Stage path | Process statics | WP-107 threads | Length consumer | v1 |
|---|---|---|---|---|---|---|---|
| `ElasticIsotropic` | ThreeDimensional, PlaneStrain | yes | none | 232c | allowlisted (3D, PlaneStrain2D) | no | **admitted** |
| `LadrunoJ2` | 232c | yes | none | 232c | not allowlisted → serial | only under `-autoRegularization` (`LadrunoJ2.cpp:352-356`), which the driver refuses at particles; with the global at 0 it falls back to its reference length | **admitted** |
| `DruckerPrager` (UW, post-`31322a47a`) | ThreeDimensional, PlaneStrain (`DruckerPrager3D`, `DruckerPragerPlaneStrain`) | gated in the `-geom hypo` lane (ADR-79 §8–9) | `materialState` (`DruckerPrager.cpp:1045`, response 5); `updateMaterialStage` returns −1 (`:1063-1064`) | 232c | not allowlisted → serial | no | **admitted** |
| `J2Plasticity` | — | — | none | yes (`J2Plasticity.h:137-139`) | not allowlisted | no | not in v1 (`LadrunoJ2` covers J2) |
| concrete and steel length consumers | — | — | — | — | — | yes, through the global | not in v1 |
| SANISAND, HS, NorSand, PDMY, PM4Sand | — | — | — | — | not allowlisted | no | later phase, with the SANISAND cell-crossing meter first (§11) |

**Parameters and stages.** Mode A: particles are elements, so the stock element walks
(`updateMaterialStage`, `parameter … material`) reach them (W12). Mode F:
`parameter <t> loadPattern <mpmTag> …` reaches the host, which forwards to its particles' materials;
`updateMaterialStage` does not reach mode F, and the manual says so. v1's one stage path is
`DruckerPrager` `materialState` through `parameter` + `updateParameter`. `MaterialStageParameter`
stops at the first element that answers, which relies on process-static stage flags; the table's
"process statics" column records which materials that bites.

**Length.** v1 offers no regularization at particles. The engine zeroes `ops_TheActiveElement`
(W9), and the only admitted consumer (`LadrunoJ2`) then uses its reference length. `NDMaterial` has
no length API; a length channel is later work, and before any length is declared for MPM, the band
width is measured as a function of basis and particles per cell — cell size is not assumed to be the
right length.

---

## 7. apeGmsh integration

Checked on apeGmsh `origin/main` @ `cf2e06593` (2026-11-13). The first draft of this section was
written against a pre-October apeGmsh and was refuted by its current ADRs (§14.3 B1, B2, M7–M10);
this is the rewrite.

### 7.1 Division of labour, and the apeGmsh rules that bind it

apeGmsh owns geometry, seeding, transport, results reading and viewing; the fork owns the mechanics
and the results format. This is the licence-clean split (§2), and it matches ADR-60a's standing
decision that apeGmsh is a **seeder, not a runtime remesher**. The apeGmsh rules this work must
follow:

- **ADR 0112 D1 — "The files are the model."** Anything apeGmsh produces for a run belongs to its
  emitted files, not to a side channel.
- **ADR 0122 (Accepted 2026-10-15) and ADR 0131 (Accepted 2026-11-03).** `.ladruno` is the canonical
  results format; apeGmsh reads one results format; no new native results file, and no new writer
  for native `results.h5`. ADR 0122 supersedes ADR 0113 for the results zone.
- **ADR 0112 D8 and ADR 0116.** The Qt / VTK viewers are frozen to bug fixes and in sunset; new
  interactive 3D views are three.js inside apeGmshViewer (TypeScript / Electron).
- **ADR 0121 D4 (proposed; direction ruled 2026-10-11) — the freeze is a lock.** A registry grows
  only through an ADR whose status is Accepted, cited by the new census entry; `scripts/land_pr.py`
  refuses a diff that sets any ADR's status to Accepted, so an agent cannot accept its own exception.
- **Charter P4 — the bridge owns tag allocation;** users never type tags.
- **Live runs are flat only** (ADR 0121, ruling 2): staged runs are deck-only.
- **Fork consumption.** A fork-only command is consumed only after its fork PR has merged into
  `ladruno`, and apeGmsh's `live-fork` lane tests the tagged fork release named in `FORK_PIN` (today
  `ladruno-v20261103`; its four keys are bumped together). Fork-only primitives carry the
  `ladruno_fork` marker and a row in `docs/concepts/backend-capabilities.md`; capability probes
  return `None` when inconclusive.
- **apeGmsh's own program board** (#1203) carries its slices; the fork's board (#937) carries the
  fork's.

### 7.2 Contract 1 — the particle set (apeGmsh → fork)

At 10⁵–10⁶ particles, per-particle commands are not an option: the set travels as HDF5, read by
`mpm particles -file {…}` (the path brace-quoted in emitted Tcl, per apeGmsh's emitted-deck rules).

| Field | Shape / type | Rule |
|---|---|---|
| format version, provenance | attributes | versioned; provenance names the apeGmsh commit and the session |
| units | attributes | declared, never implied |
| `ndm`, plane-strain thickness | attributes | thickness is required when `ndm` = 2 |
| grid | origin, cell size, extents | the engine's grid and the viewer's; the seeding is checked against it |
| particle id | int64 [N] | stable and unique, never reused within a run |
| position | float64 [N, ndm] | |
| volume | float64 [N] | positive; zero or negative is refused |
| mass or density | float64 [N] | exactly one of the two, declared |
| GIMP half-lengths | float64 [N, ndm] | from the seeding cell; GIMP is the v1 basis and volume alone does not define the particle domain |
| initial velocity | float64 [N, ndm] | zero by default |
| initial stress | float64 [N, 6] (4 in plane strain) | Voigt order and sign convention stated |
| body id | int32 [N] | |
| material | int32 [N] | a **tag resolved at emit** from the body's material label (charter P4) — the emitter writes this column, or the whole file |
| flags | uint32 [N] | |

**Placement.** Either a new FEMData stream inside `model.h5` — then apeGmsh's `compose-streams`
rule and its zone, corpus and app-target bumps apply — or a sidecar recorded as an explicit
exception to ADR 0112 D1. The apeGmsh ADR decides (§7.4); this contract does not. Datasets are
chunked and compressed.

### 7.3 Contract 2 — particle results (fork → apeGmsh)

A **new particle group inside `.ladruno`**, written by `LadrunoRecorder` or by an MPM writer that
shares its file. It is an **additive** change: `FORMAT_MINOR_VERSION` goes from 1 to 2 and
`FORMAT_VERSION` stays 1 (`LadrunoRecorder.cpp:989-994`). A major bump would be wrong: apeGmsh
treats one as a breaking layout change and its readers refuse anything outside `FORMAT_VERSION` 1
(`schema/_versions.py:73-90`, `_ladruno.py:633-637`, `ladruno.ts:56`), so every new file would stop
opening. With the minor bump, today's readers keep opening the files and ignore the new group until
232d teaches them to read it. The change also needs a fork format-document PR and the fork's own
format tooling updated (`Ladruno_scripts/ladruno_recorder_tests/ladruno_format.py`, which requires
`FORMAT_VERSION == 1`, and `make_synthetic.py`). apeGmsh reads the group in
`src/apeGmsh/results/readers/_ladruno.py` and in `apeGmshViewer/src/reader/` (the `ladruno*.ts`
readers, which already read `FORMAT_MINOR_VERSION`).

- **Per output:** time and step; the particle id list (particles may be removed or added; ids are
  never reused); positions (float64); velocities; stress; named state variables; volumes or GIMP
  half-lengths if they evolve. Grid metadata once per stage.
- **Cadence: never every step.** 10⁶ particles × about 20 doubles is about 160 MB per output, and
  explicit runs take 10⁴–10⁶ steps. Output every N steps or every Δt_out, with field selection;
  chunked, compressed datasets.
- **Mode A** writes this group too, not the element / GP layout (W4, W13).

### 7.4 The apeGmsh work packages

1. **An apeGmsh ADR**, numbered from apeGmsh's decisions directory on `origin/main` at the time of
   writing (the highest today is 0133), for the particle-set stream and a particle-cloud result
   ontology. **The maintainer accepts it**; an agent cannot (ADR 0121 D4). It gates 232a.
2. **Seeding — lattice first.** n particles per grid cell per axis, kept where the point lies inside
   a body (Gmsh's inside test, or the body's own mesh): positive volumes, uniform particles per cell,
   GIMP half-lengths of cell / (2n) by construction. Gauss-point seeding comes later and needs new
   rules. apeGmsh's `fem/_quadrature.py` states its own scope — "Minimal Gauss quadrature rules",
   sufficient for HRZ mass lumping, where "the absolute scale of the quadrature weights cancels" —
   and its tetrahedron rule is Keast's 11-point rule with a **negative centroid weight**
   (−0.013155…), which would make a negative-volume, negative-mass particle in every tet. So
   Gauss-point seeding needs positive-weight rules (triangle 1 / 3 / 4 points, tetrahedron 1 / 4;
   pyramids refused — `_shape_functions.py` has none, nor wedge15 / 18), a per-family "supported for
   seeding" table with loud refusal, a signed-Jacobian check (`compute_jacobian_dets` returns
   `abs(det)`, `_shape_functions.py:903`), and a particles-per-cell uniformity report against the
   declared grid.
3. **Readers:** the `.ladruno` particle group in Python (`_ladruno.py`) and TypeScript
   (`apeGmshViewer/src/reader/`).
4. **Viewer:** a three.js `THREE.Points` layer in apeGmshViewer with per-output position buffers,
   coloured by stress or state, after the node layer (`apeGmshViewer/src/renderer/nodes.ts`,
   `makeNodePoints`), with apeGmsh's visual check. Not the Qt viewers, which are frozen; ADR 0087's
   "dots" are Qt chrome glyphs whose own row says they are never used for any other point cloud.
   This is TypeScript feature work, not pattern-following.
5. **Bridge primitives** for `mpm` — only after 232c merges **and** a fork release containing `mpm`
   is cut and `FORK_PIN` is bumped; the `ladruno_fork` marker; the capability-map row; probes that
   return `None` when inconclusive; staged runs deck-only.
6. **Later:** the FE → MPM handoff (§1) and MPM → FE remeshing (§7.5).

### 7.5 Instant Meshes — studied at the owner's request (2026-11-13)

**What it is.** A field-aligned remesher (Jakob, Tarini, Panozzo, Sorkine-Hornung, *Instant
field-aligned meshes*, ACM TOG 34(6), 2015) that turns a surface into an isotropic, quad- or
triangle-dominant surface mesh, fast and robustly. The original, `wjakob/instant-meshes`
(BSD-3-style `LICENSE.txt`; last pushed 2022-01-03), has an interactive GUI **and a batch mode**
(`-o, --output`, `src/main.cpp:150`), and it accepts **point clouds** as input (`main.cpp:148`, with
`--knn`). The Python binding `pynanoinstantmeshes` takes a triangle mesh only, in float32 — a
precision problem for engineering coordinates unless they are shifted to a local origin. Its single
release (0.0.3, 2025-01-03) ships wheels for Windows, Linux and macOS arm64 (cp312 is abi3), and its
licence is unverified (§2). QuadriFlow (Huang et al. 2018) is a more robust quadrangulation built on
the same fields.

**Where it would help:** MPM → FE, turning a deformed body back into an FE model — boundary
particles extracted from the volumetric cloud (the real missing step: the cloud is volumetric, so
its surface must be found first) → a surface from them, or Instant Meshes' own point-cloud mode →
a remeshed surface → a Gmsh volume mesh → an FE model. Also the deformed-boundary remeshing a
RITSS loop needs (Hu & Randolph 1998).

**Where it would not:** seeding (a volume problem) and the MPM grid (structured arithmetic).

**Gate before adopting:** compare it with Gmsh's own discrete-surface remeshing
(`classifySurfaces` / `createGeometry`) on the same deformed surfaces; check that the output is
closed and manifold before any volume mesh (Instant Meshes guarantees neither); adopt it only if it
wins on robustness or quality, and only through a binding with a verified licence or the original's
batch executable.

**Not to be confused with** TencentARC's *InstantMesh* (2024, Apache-2.0), a generative model that
produces a 3D mesh from a single image — irrelevant to engineering geometry.

---

## 8. Verification and tests — a general suite

Slice column: **S1** = slice 1; **P*n*** = phase *n* of §9.

| Kind | Test | What it checks | Slice |
|---|---|---|---|
| Verification | 1D vibrating bar | transfers; PIC vs FLIP dissipation | S1 |
| | elastic column under self-weight — closed-form vertical Cauchy stress σ = ρ₀g(l₀ − Y) (Charlton et al. 2017; AMPLE) | large deformation; explicit with local damping in S1, implicit quasi-static in P1; convergence per basis | S1, P1 |
| | method of manufactured solutions (Wallstedt & Guilkey 2008; Kamojjala et al. 2015) | spatial convergence order per basis | P3 |
| | cell-crossing meter | the stress jump, as a number | S1 |
| | no-cell-crossing re-seed: a mode-A step in which no particle changes cell equals the re-seeded result | the unconditional `domainChange()` | P1 |
| | plate with a hole — state the hole-diameter / width ratio (Howland's factor depends on it; 3.0 is the infinite-plate limit) | 2D accuracy | P1 |
| | rigid-body translation and rotation of a stressed body | objectivity of the new kernel entry; zero spurious stress | S1 |
| | two-body elastic impact (Sulsky et al. 1994) | momentum and energy | S1 |
| | explicit energy ledger | conservation, FLIP vs PIC | S1 |
| | block on a frictional incline | the frictional-wall condition | S1 |
| | large-deflection cantilever against FE | the implicit mode | P1 |
| | quiet start under K0 | staged initial stress | S1 |
| | small-mass stability | the mass cutoff | S1 |
| | critical Δt: stable just below the estimate; refused above it | the CFL check | S1 |
| | failed-step retry (mode A): a step that fails before commit is cut and retried, and G2P runs once | retry without a double G2P | P1 |
| | refusal at commit aborts the run cleanly, particles unmoved — with the first refusing material | W17 | P5 |
| | recorder synchrony | positions and stresses from the same step | S1 |
| | lifecycle: `wipe`, `reset` with and without `mpm reset`, `remove loadPattern` | the host's lifecycle | S1 |
| | refusals: each non-admitted material, a non-whitelisted integrator, free FE DOFs, `VariableTransient`, database save with particles, `mpm` in OpenSeesSP / MP | loud refusal, never a wrong answer | S1, P1 |
| | engine parity: classic Tcl and Python give the same particle output (DL Tcl is registered but runs in no shipped binary) | registration; the `loadPattern` stage path | S1 |
| | thread-count determinism (allowlisted materials) | the threaded particle loop | S1 |
| | contract round trips: a fixture file each way, read by the other side | both contracts | S1 |
| | particle-id stability across outputs, with removal | Contract 2 | S1 |
| | seeding guard: zero or negative volume refused | Contract 1 | S1 |
| | restart of particle state (database) | serialization | P1 |
| | Terzaghi + Mandel on particles (+ TIMs' first-settlement diagnostic) | two-phase | P6 |
| Validation | granular column collapse (Lube et al. 2004, 2005) | runout against experiment | S1 |
| | dam break (Martin & Moyce 1952) | only if a fluid law is admitted | — |
| Applications | shallow footing — gated on the **non-associated** mild-cone (φ_ps = 27.47°) plane-strain Prandtl–Reissner control, exact q₀·N_q (138.9 kPa). The FE reference is a monotone sequence, not a converged value: H8 B-bar with the post-`31322a47a` `DruckerPrager` gives 1.0020 / 0.9517 / 0.9260 of exact at 4 / 8 / 16 elements across B, the finest rung softening after its peak (ADR-79 `:560-562`, `:588-590`). The MPM gate is a tolerance band around exact, fixed at that phase's kickoff from a resolution study (cells across B, particles per cell) — never one FE number. ADR-95's gate (#722: 1.0842 / 0.9938 / 0.9513 at h₀ = 1.0 / 0.5 / 0.25) is a different gate and is not this one. The associated leg is not a gate: it never plateaus (ADR-79 `:564-566`) | the TIMs geometry | P1 |
| | tailings flow liquefaction → runout (NorSand) | the flow-slide case | P5 |
| | cross-code against Berkeley and Kratos; NorSand matched-material track once WP-144 lands | Tier-2 screens, never gates | — |

**Mutation gate (ADR-87 D2)**, through the fork's existing mechanism: a new family `MPM` in the
`LADRUNO_MUTATE_FAMILY` list (`CMakeLists.txt`, a vanilla-ledgered line), in the options of
`.github/workflows/ladruno_mutation.yml` and in `Ladruno_scripts/mutation_gate.py`; the mutant is
built with `Ladruno_scripts/mutation_build.bat` (the orchestrator launches it) or the CI lane. The
engine defines what each mode sabotages:

- `ZERO` — the P2G internal force; the vibrating bar catches it;
- `SCALE` — the G2P position increment; the self-weight column or a large-amplitude bar catches it
  (a small-amplitude bar can stay green).

Ordering hazards — G2P after `Domain::commit`, and a conditional `domainChange()` in mode A — are
not expressible as these modes. Their regression tests (recorder synchrony; the no-cell-crossing
re-seed) guard them. The `ops_TheActiveElement = 0` guard stays a rule without a test in v1: no
admitted material observes it (`LadrunoJ2` reads the global only under `-autoRegularization`,
which the driver refuses); its test arrives with the first admitted length consumer.

---

## 9. The first slice and the phase order

**Slice 1 — the fork and apeGmsh together, contract first.**

- **Gate before contract-bound code (232a):** the apeGmsh ADR accepted by its maintainer; both
  contracts as field-by-field schemas with both repos' checklists ticked; a fixture file each way.
  The engine's numerics start earlier, on 232b's frozen interface; only the particle-set reader and
  the `.ladruno` writer wait for 232a.
- **Fork:** explicit mode only, hosted as a `LoadPattern`; 2D plane strain and 3D; serial (threads
  only for allowlisted materials); linear + GIMP; MUSL with PIC/FLIP; `ElasticIsotropic`,
  `LadrunoJ2` and `DruckerPrager` through the new kernel entry; the `materialState` stage path
  through `parameter … pattern`; gravity ramp through the host's time series, grid planes, a
  frictional wall, quiet K0 start; local non-viscous damping (the FLAC form) for quasi-static
  checks; the critical-Δt estimate and refusal; `mpm` as one `LadrunoCommandTable.h` row, refused in
  OpenSeesSP / MP; particle output as a `.ladruno` group with the `FORMAT_MINOR_VERSION` bump and
  its format document.
- **apeGmsh:** lattice seeding writing the particle set; the `.ladruno` particle readers (Python,
  TypeScript); the three.js points layer. Bridge primitives follow once 232c has merged and a fork
  release is pinned.
- **Tests:** the S1 rows of §8 and the mutation gate.

**Then, in order:**

1. implicit dynamic and quasi-static through A (unconditional rebuild, integrator whitelist) — the
   cantilever, the non-associated mild-cone footing, retry with kinematics restored, the
   no-cell-crossing re-seed test, recorder synchrony, restart;
2. **M0b measurements:** the per-step `domainChanged()` cost against one Newton iteration at
   10³–10⁵ particle-elements (PARDISO and UmfPack, with the empty-node treatment); explicit
   throughput across materials (elastic → Drucker–Prager → SANISAND, builds pinned); memory per
   particle-element — deciding whether C is ever needed;
3. B-spline basis and the manufactured-solution orders;
4. anti-locking in rate form;
5. soil models (SANISAND, HS, NorSand), with the SANISAND cell-crossing meter first;
6. two-phase particles (Terzaghi + Mandel; stabilized equal-order versus a semi-implicit route);
7. multi-body frictional contact;
8. the FE–MPM seam, implicit first;
9. the FE → MPM handoff (needs the state-transfer primitive, §1);
10. MPM → FE remeshing (Instant Meshes against Gmsh, §7.5);
11. parallel.

Model and effort per work package: [[_wp232_orchestration_plan]]. The operating skill — boot,
dispatch briefs, worker rules, fork checklist — is `.claude/skills/ladruno-mpm/` (`/ladruno-mpm <wp>`).

---

## 10. Decisions requested

- **D1** Baseline F (explicit, hosted as a `LoadPattern`) + A (implicit and quasi-static); C deferred.
- **D2** v1 scope: no explicit FE–MPM seam (free FE DOFs refused), no soil regularization, serial.
- **D3** Contract-first co-design under a new apeGmsh ADR that its maintainer accepts: particle
  results as a `.ladruno` group, the particle set tied to the emit, lattice seeding first.
- **D4** Open the MPM `program` issue claiming WP-232 onward, with a row on the fork board (#937) and
  a link on apeGmsh's (#1203).
- **D5** Instant Meshes as a later apeGmsh evaluation, gated against Gmsh's own remeshing and on a
  verified binding licence.
- **D6** The orchestrator session (see the plan).
- **D7** Amend ADR-86 D2: the FE → MPM handoff (D2's "and B") moves from the v1 driver to phase 9.
  **Approved by the owner, 2026-11-13.**

---

## 11. Risks

- **The value proposition.** Sordo, Rathje & Kumar had PM4Sand at the FEM end and deliberately used
  Mohr–Coulomb with softening for the MPM runout, citing MPM's stress inaccuracy; and in FE, the
  TIMs SANISAND footings are 2.4–4× too strong against Siddiquee. Cyclic sand models at particles
  must earn their place: the SANISAND cell-crossing meter runs first in the soil phase.
- **New numerics.** The incremental kernel entry and the plane-strain embedding each need an oracle
  before C++.
- **The host.** A `LoadPattern` host has one precedent (ADR-73). `reset` reaches it only as a
  re-`applyLoad`, so `mpm reset` exists and is tested; the stage path must be spelled
  `loadPattern` (classic Tcl's `pattern` keyword means only `pattern <tag> lambda`).
- **Stage routing.** `updateMaterialStage` does not reach mode F; soil models staged in mode F use
  `parameter … pattern`.
- **Empty-support nodes** in mode A; decided by measurement.
- **Sequencing across repos.** A fork release and a `FORK_PIN` bump sit between 232c and 232e, and
  the apeGmsh ADR needs its maintainer's acceptance before 232a can close.
- **Output volume.** Particle output at full cadence is untenable (§7.3).
- **Numbering.** WP-232 is provisional. ADR numbers have now collided three times in this lineage:
  98 (the membrane ADR), 88 (the merged H5DRM work against the uncommitted state-transfer draft), and
  apeGmsh 0103 (the SANISAND deck contract, which the state-transfer draft names as its sibling).
- **The recurring failure mode.** Nine claims across this lineage were refuted by facts already on
  record — four before this document, five in its first draft. The mitigation is §4's rule, now
  extended to apeGmsh `origin/main` and to every ADR number cited.

---

## 12. References

**Libraries:** Kratos MPMApplication · geomechanics/mpm · CB-Geo mpm · Uintah · Karamelo · NairnMPM
· Anura3D · AMPLE · Taichi · GeoTaichi · Matter (URLs and verification in the predecessor draft's §2,
`origin/rescue/adr98-mpm-scoping`) · Instant Meshes `github.com/wjakob/instant-meshes` ·
QuadriFlow `github.com/hjwdzh/QuadriFlow` · `pypi.org/project/pynanoinstantmeshes`.

**Method:** Sulsky, Chen & Schreyer (1994) CMAME 118 · Bardenhagen, Brackbill & Sulsky (2000) CMAME
187 · Guilkey & Weiss (2003) IJNME 57 · Bardenhagen & Kober (2004) CMES 5(6) · Steffen, Kirby &
Berzins (2008) IJNME 76 · Wallstedt & Guilkey (2008) J. Comput. Phys. 227 · Sadeghirad, Brannon &
Burghardt (2011) IJNME 86 · Kamojjala et al. (2015) Eng. Comput. 31 · Jakob, Tarini, Panozzo &
Sorkine-Hornung (2015) ACM TOG 34(6) · Charlton, Coombs & Augarde (2017) Comput. Struct. 190 ·
Huang et al. (2018) Comput. Graph. Forum 37(5) · Coombs & Augarde (2020) Adv. Eng. Softw. 139 ·
Giridharan, Gowda, Stolle & Moormann (2020) Soils Found. 60(4):989–1000 · Wilson, Wüchner & Fernando
(2021) IJNME 122 · Tjung & Soga, arXiv:2111.13584 · Sordo, Rathje & Kumar, arXiv:2412.08040 · Lube
et al. (2004) JFM 508; (2005) PRE 72 · Martin & Moyce (1952) Phil. Trans. R. Soc. A 244 · Hu &
Randolph (1998) IJNAMG 22.

**Fork:** [[86_geomechanics_libraries_and_mpm_scoping_report]] · [[87_ladruno_depth_with_width_adr]]
· [[79_ladruno_geom_hypo_adr]] · [[90_ladruno_regularization_tims_report]] ·
[[95_prandtl_reissner_quadratic_root_cause]] · [[71_ladruno_up_family_adr]] ·
[[73_ladruno_porous_overlay_adr]] · [[107_ladruno_openmp_element_loop]] ·
[[52_ladruno_integrator_strengthening_adr]] · [[60a_pfem_crosspollination_amendment]] ·
[[39_ladruno_contact_domain_adr]] · the uncommitted state-transfer draft
`88_ladruno_state_transfer_adaptivity_adr.md` (owner's main checkout, 2026-09-02; not the merged
ADR-88).

**apeGmsh** (`origin/main` @ `cf2e06593`): ADRs 0112 (the files are the model; the read-only app),
0116 (render technologies), 0121 (charter v2), 0122 (`.ladruno` the canonical results format), 0131
(all in `.ladruno`); `architecture/charter.md` P4; `FORK_PIN`; `docs/concepts/backend-capabilities.md`.

---

## 13. Implementation log

*No code is authored under WP-232 itself; the code WPs log here as they land.*

- **2026-11-13** — the owner's go. Orchestrated from the authoring session (Opus 5.5) at the
  owner's choice; fork merges delegated to the orchestrator; apeGmsh PRs and merges go through
  apeGmsh's T coordinator session. 232a and 232b launched in parallel.
- **2026-11-13** — the owner approved D7: ADR-86 D2's FE → MPM handoff moves to phase 9 (amendment
  note added to ADR-86 D2).

---

## 14. Review logs

### 14.1 First review — 2026-09-25, Fable, of the ADR-98 draft. Not clean: one blocker, five majors.

The blocker was W1: OpenSees integrators discard nodal state written between steps, so the draft's
"orchestrate through the Nodes, works with any integrator" was false. The majors: Kratos-explicit has
no DOF set; a "shared grid" FE–MPM coupling is incoherent for deformable bodies; the differentiator
was overstated; WP-107 threads no soil model; the materialize/retire design omitted node churn and
recorder lifetimes. All verified by the author before acceptance; full table in the predecessor's §15.

### 14.2 Second review — 2026-11-13, Fable, of the general-driver brief. Not clean: two blockers, eight majors.

| # | Finding | Disposition here |
|---|---|---|
| F1 | **Blocker:** the hypo kernel's entry needs reference gradients and nodal displacements; a reset grid has neither; the lane is 3D-only | W8; §6 new incremental entry + oracle; plane strain is new work |
| F2 | **Blocker:** "regularization via `getCharacteristicLength`" — the channel is a global, no soil model reads it, the TIMs band scaling fails its own gates | W9; §1 non-goal; §6 no regularization at particles in v1 (W9's wording corrected in §14.3, M5) |
| F3 | the precedent is the fork's staggered driver, not `PFEMAnalysis`; registration is per interpreter | §5.1 command family (its scope and registration corrected in §14.3, M1 and m1) |
| F4 | the step order was wrong in four places (domain change, sub-stepping and retry, recorder timing, explicit commit/loads/Δt) | §5.2, §5.4 (completed in §14.3, M3 and M4) |
| F5 | the new integrator hook is unnecessary for Newmark/HHT/GA and ill-defined for the rest | §5.3 whitelist; hook dropped |
| F6 | "any NDMaterial" breaks on stage routing, process statics and the allowlist | §6 admitted-material table; W12, W14 (mechanism and stage path in §14.3, M6) |
| F7 | explicit FE + MPM coexistence unspecified | §1 non-goal for v1; free FE DOFs refused (§14.3, M4) |
| F8 | the 1108–1152 kPa footing anchor is a mesh-specific B-bar number on an ill-posed problem | §8 footing gated on the **non-associated** mild-cone control — the first fold mislabelled it "associated" (corrected in §14.3, B3) |
| F9 | the test suite was not a general-driver suite | §8 (gaps closed in §14.3, m8 and m9) |
| F10 | the A/F/C scores no longer meant anything | §0 baseline F + A; §9 M0b redefined |
| F11 | solve modes under-specified | §5.3 (quasi-static corrected in §14.3, M2) |
| F12 | brief §4 verified, with corrections (HS is young; the band work failed gates) | §4 asset rows |
| F13 | the GP layout is not a moving-particle format | W13; §7.3 (the apeGmsh side rewritten in §14.3, B1 and B2) |

**Author's verification before acceptance:** F1 (`LadrunoHypoKernel.h:43, :166-176`;
`OPS_LadrunoUP.cpp:512`), F2 (`Element.cpp:49-61`; consumer list; zero soil hits), F3
(`LadrunoStaggeredDriver.{h,cpp}`), F4 (`analyzeSubLevel` at
`DirectIntegrationAnalysis.cpp:177/280/293`; recorder order in `Domain::commit`), F5
(`firstStep = true` at `CentralDifferenceLadruno.cpp:362/523`), F6 (allowlist on `origin/ladruno`),
F8 (ADR-79 `:541`, `:581`, `:595-598`). All confirmed — and F8's label was still inverted in the
fold, which the third review caught.

### 14.3 Third review — 2026-11-13, Fable, of this document's first draft. Not clean: three blockers, ten majors, twelve minors.

Reviewed on fork `origin/ladruno` @ `aa6229bf9` and apeGmsh `origin/main` @ `7630ee636`; the
author re-verified on apeGmsh `cf2e06593`. The review's verdict: the fork-side architecture largely
survives; the apeGmsh half was written against a pre-October apeGmsh and had to be redone; the
footing gate had to be relabelled; opening the draft PR and the program issue is fine.

| # | Finding | Author's check | Disposition |
|---|---|---|---|
| B1 | **Blocker:** Contract 2 as a separate file contradicts apeGmsh ADRs 0122 / 0131 — `.ladruno` is the one results format | confirmed: both Accepted; `LadrunoRecorder.cpp:989` owns `FORMAT_VERSION` | §7.3: a particle group inside `.ladruno` |
| B2 | **Blocker:** the viewer targeted the frozen Qt viewers through ADR 0087, a chrome design system | confirmed: ADR 0112 D8, ADR 0116; `nodes.ts` `makeNodePoints` | §7.4(4): a three.js points layer under the new apeGmsh ADR; plan 232d re-scoped |
| B3 | **Blocker:** the footing gate was labelled "associated"; the numbers are the non-associated control, drifting 7.4 % | confirmed: ADR-79 `:560-566`, `:588-590`; ADR-95 `:41` | §8 relabelled; tolerance band at kickoff; build named; ADR-95's numbers separated; §14.2 F8 corrected |
| M1 | the staggered driver whitelists for clean retry, not re-seeding; it is transient-only | confirmed: `LadrunoStaggeredDriver.h:44-48`, `.cpp:152-172`; `OpenSeesCommands.cpp:4083-4092` | §4 quotes its scope; §5.3 justifies the whitelist on its own evidence |
| M2 | `StaticAnalysis` has no `analyzeStep` | confirmed: `StaticAnalysis.h:66`; no sub-levels in `StaticAnalysis.cpp` | W16; §5.1 `analyze(1, false)` + own flush; `VariableTransient` refused |
| M3 | "if they change the domain" re-creates W1 | confirmed: `DirectIntegrationAnalysis.cpp:207-217` | §5.3 unconditional `domainChange()`; test and mutant |
| M4 | the explicit step list was incomplete and partly wrong | confirmed: `Domain.cpp:2148`, `:2170`, `:2343-2345`; `DirectIntegrationAnalysis.cpp:184-186`, `:194`, `:273` | §5.2 ordered list; free FE DOFs refused; critical Δt; no Domain hook |
| M5 | W9's "no active element in explicit mode" is false | confirmed: `Element.cpp:131`; `Domain.cpp:493`; `LadrunoJ2.cpp:353-354`; `J2Plasticity.h:137-139` | W9 reworded; the engine zeroes the global per thread; v1 J2 is `LadrunoJ2` |
| M6 | "registers as a parameter target" had no mechanism; the v1 stage path was misnamed | confirmed: the element walks; `DruckerPrager.cpp:1045`, `:1063-1064` | **fix partly declined:** the proposed zero-node Domain element would exit the process at the first analysis (W15, `FE_Element.cpp:126-129`, found while verifying); the host is a `LoadPattern` instead; stage path `materialState` via `parameter … pattern` |
| M7 | `_quadrature.py` is not a seeding library; Keast's centroid weight is negative | confirmed: −0.013155… at the centroid; `_shape_functions.py:903` | §7.4(2): lattice first; positive-weight rules for Gauss-point seeding |
| M8 | the contracts were not specified enough to design before code | accepted | §7.2 field table; §7.3 cadence and precision |
| M9 | apeGmsh's protocol and sequencing were under-specified | confirmed: `FORK_PIN` `ladruno-v20261103`; ADR 0121 D4 — **proposed, direction ruled 2026-10-11**, not ratified as the review implied | §7.1 rules; the plan adds the release, the pin bump, the apeGmsh ADR gate and #1203 |
| M10 | apeGmsh ADR 0103 is the SANISAND deck contract, not a transfer bundle | confirmed | removed; the citation came from the uncommitted state-transfer draft, whose sibling pointer is stale |
| m1 | registration is one `LadrunoCommandTable.h` row feeding three engines | confirmed: `:101`; `tests/test_ladruno_command_registry.py` | §5.1 |
| m2 | caveats for mode-A re-pointing | confirmed: `Domain.cpp:2294-2423` | §5.3 element obligations |
| m3 | `LadrunoRecorder` keys on a topology fingerprint, not the stamp | confirmed: `LadrunoRecorder.cpp:664-694` | W4; the recorder quirk and plan rule 4 corrected |
| m4 | the integrator quirk partly duplicates WP-215's | confirmed | cross-referenced; WP-215's dangling "ADR 98" pointer fixed |
| m5 | the length quirk overlaps two LEGACY entries | confirmed | turned into a follow-up on `LEGACY-ops-theactiveelement-is-a-mutable-global` |
| m6 | sub-stepping only when `numSubLevels != 0` | confirmed: `:176` | W10 and the quirk |
| m7 | Instant Meshes: batch mode, point clouds, binding licence, float32, closure | confirmed: `main.cpp:148-150`; licence API null for the binding | §2, §7.5 |
| m8 | the explicit self-weight column needs damping | accepted | slice 1 adds local damping; the implicit column is P1 |
| m9 | §8 lacked general-driver tests | accepted | §8 rows added |
| m10 | WP-232 moves ADR-86 D2's "and B" without saying so | confirmed: `86_…report.md:500-501` | the lineage note; D7 |
| m11 | the `Agent` tool takes an effort parameter | confirmed in this session's tool schema | plan §0 corrected |
| m12 | nits: number reservation, the `COPYRIGHT` inference, the goal sentence | accepted | front matter, §2, §1 |

**Found by the author while verifying:** W15 (above); and a third number collision — on `ladruno`,
ADR-88 is the H5DRM higher-order-elements work (#778 / #800), while the state-transfer draft the
first version cited as "ADR-88" (and "ADR-88 D5") is uncommitted in the owner's main checkout. It is
now cited by file name. The review did not re-verify arXiv 2605.28525 or Sordo et al.; they come
from the predecessor's survey.

### 14.4 Fourth review — 2026-11-13, Fable, of #1074 @ `bc8b8f11c` (ADR, plan, skill). Not clean: two blockers, six majors, nine minors.

Read on fork `origin/ladruno` @ `0103c1017` and apeGmsh `origin/main` @ `963d8ac9d`. The verdict:
the architecture and the declined M6 fix survive (W15 confirmed; the `LoadPattern` host sound for
`wipe`, `applyLoad` and reset); #1074 can merge once findings 1–8 are fixed in the text. Every
finding was checked by the author before acceptance.

| # | Finding | Author's check | Disposition |
|---|---|---|---|
| 1 | **Blocker:** the "`-DLADRUNO_MUTATE_MPM=<n>`" mutation gate does not exist — the mechanism is a closed family list with ZERO / SCALE / IDENT modes; two of the five mutants could not go red | confirmed: `CMakeLists.txt:243-268`; `ladruno_mutation.yml` options; `Ladruno_scripts/mutation_gate.py`, `mutation_build.bat` | §8: a new `MPM` family through the existing mechanism; ZERO and SCALE defined; ordering hazards guarded by regression tests |
| 2 | **Blocker:** "`FORMAT_VERSION` 1 → 2" would make every new `.ladruno` unreadable by current apeGmsh | confirmed: `LadrunoRecorder.cpp:989-994` has `FORMAT_MINOR_VERSION`; apeGmsh `_versions.py:73-90`, `_ladruno.py:633-637`, `ladruno.ts:56`; fork `ladruno_format.py` | §7.3, §0, §9: a minor bump (1 → 2), the major stays 1; the running 232a worker corrected before it froze the contract |
| 3 | classic Tcl's `parameter … pattern` means only `pattern <tag> lambda` | confirmed: `TclParameterCommands.cpp:134-135`; the generic branch matches `loadPattern` (`:235`, `:248`) | spelled `loadPattern` everywhere |
| 4 | "three interpreter engines": the DL Tcl engine runs in no shipped binary | confirmed: `SRC/tcl/CMakeLists.txt:94` (`TclInterpreter.cpp` commented out); `build.bat:156` | §5.1, §8: registered in three, run in two |
| 5 | "Zone-A green" is unreachable as briefed: CI runs `pytest -m zone_a` on Ubuntu without h5py | confirmed: `ladruno.yml:287`; its pip line installs no h5py | fork checklist: the `zone_a` marker, h5py added to the CI install, the GCC / Linux build, and the other gates (`check_duplicate_defs.py`, ruff F811, `stamp_headers.py --check`) |
| 6 | "revert the particles' trial states and fail the step" is impossible: a refusal is seen after siblings committed, and `Domain::commit` aborts | confirmed: `Domain.cpp:2343-2355`; only `LadrunoConcrete3D`, `LadrunoHS` and `LadrunoSANISAND` declare refusals | W17; §5.2 step 6, §5.3, §5.4: a refusal aborts the run; mode A retries only pre-commit failures; the refusal test moves to P5 |
| 7 | the stale-length test and its mutant cannot go red: `LadrunoJ2` reads the global only under `-autoRegularization`, which the driver refuses | confirmed: `LadrunoJ2.cpp:352-356`; `dmgRegOn` private (`LadrunoJ2.h:128`) | test and mutant dropped; the guard stays a rule until the first admitted length consumer |
| 8 | the merge delegation contradicts AGENTS.md ("The owner merges; agents do not") and the plan's own rule 2 | confirmed | AGENTS.md carries the scoped exception (a program issue records the owner's delegation); the skill map and plan rule 2 aligned |
| m | the D4 verdict and handoff pack missing from the warrant; 232c's start vs §9's gate; the plan's front matter still naming a Fable orchestrator; "never edits `SRC/`" vs "takes it on the third round"; 232b's Voigt conventions; `mpm` in OpenSeesSP / MP and the object broker; "both analysis bridges" in slice 1; W4's other stage triggers; `LadrunoDamage.h` listed as a length reader | all confirmed (`87_…adr.md:146-153`; `LadrunoRecorder.cpp:664-694`; no `ops_TheActiveElement` in `LadrunoDamage.h`) | each fixed in the ADR, plan, skill or the LEGACY follow-up; 232b sent the convention (`tests/hypo_reference.py:15`) while running |
