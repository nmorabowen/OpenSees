---
title: Upstream PR campaign — jaabell/OpenSees (ladruño)
project: Ladruno
tags:
  - upstream
  - planning
  - provenance
---

# Upstream PR campaign — jaabell/OpenSees `ladruño`

Plan for porting the Ladruno fork's work to **`github.com/jaabell/OpenSees`,
branch `ladruño`**, as a sequence of consolidated, squash-merged pull requests
authored by the Ladruno team (Nicolas Mora Bowen, Patricio Palacios, José Abell).

## 0. Ground truth (measured 2026-07-22)

- `jaabell/ladruño` = jaabell/master + 45 upstream-sync commits — effectively
  **current upstream OpenSees as of 2026-07-13** (contains upstream PRs up to
  #1765: CreepMaterial removal, elastic-beam shear option, UmfPack numeric
  reuse #1762, VTKHDF region filtering, H5DRM load-const, …).
- Our `ladruno` diverged from that line **1018 commits ago** (merge-base
  `e1237189a`). Net footprint vs merge-base: **1203 files, +373k/−26k lines**;
  of that, **241 new `SRC/` files** and **142 modified vanilla `SRC/` files**.
  The rest (817 files) is fork-internal (`Ladruno_*`, `tests/`, `.claude/`,
  CLAUDE.md) and **never goes upstream**.
- Consequence: **our git history is not portable.** No cherry-picking. Every
  package is built as a *fresh branch off `jaabell/ladruño`*: new files copied
  and adapted, vanilla edits re-applied hunk-by-hunk (they are all findable via
  `grep -rn "// Ladruno" SRC/` + `LEDGER_vanilla_files.md`), then squashed into
  one commit / squash-merged PR.

## 1. Hard rules for every upstream PR

1. **Authorship — team only, no AI traces.**
   - Commit author: Nicolas Mora Bowen. Squash commit / PR-branch commits carry
     `Co-authored-by: Patricio Palacios <…>` and
     `Co-authored-by: Jose A. Abell <…>` trailers. *(Fill real emails once —
     they must match GitHub accounts to get avatar credit.)*
   - **No** `Co-Authored-By: Claude` trailers, **no** "Generated with Claude
     Code" lines, anywhere in the port branches — the squash commit inherits
     trailers from branch commits, so the branch must be clean from the first
     commit, not scrubbed at the end.
   - File headers: re-stamp **without "Guppi"** (199 SRC files currently say
     `Nicolas Mora Bowen · Patricio Palacios · José Abell · Guppi`). Add an
     `--upstream` mode to `Ladruno_scripts/stamp_headers.py` that emits a sober
     header: standard OpenSees/PEER copyright block + a short
     `// Developed by: N. Mora Bowen, P. Palacios, J.A. Abell (Ladruño project)`
     credit + the literature references for that class. Recommend dropping the
     ASCII-art banner block for upstream files (jaabell's call; default = drop).
   - Nothing from `.claude/`, `CLAUDE.md`, `Ladruno_implementation/`,
     `Ladruno_internal/`, `Ladruno_scripts/` (except ported test assets),
     `banner_*` ships in a package.
2. **Documentation with references is part of "done".** Each ported class gets:
   full citations in the header (author, title, journal, year, DOI where
   known); a theory summary + verification evidence in the PR body; inline
   comments where the code implements a specific equation ("Eq. (12) of
   Noh & Bathe (2013)"-style). Rewrite fork-internal `ADR-##` comment
   references (127 files) into either the paper citation or plain prose — the
   ADR docs won't exist upstream.
3. **Squash merge** every PR → one consolidated commit per package on `ladruño`.
4. **Each package must build and pass its ported tests against the `ladruño`
   base**, not against our fork. Keep a dedicated integration worktree checked
   out on jaabell-based branches; port the relevant Zone-A pytest subset into
   the package (scrubbed of fork-internal paths).
5. **Ledger discipline continues**: when a package merges upstream, mark the
   corresponding rows in `LEDGER_implementations.md` / `LEDGER_vanilla_files.md`
   with the upstream PR number.

## 2. Handling the vanilla-side changes

Three distinct classes (from `LEDGER_vanilla_files.md`), handled differently:

- **A. Pure bug fixes** → their own small, early PRs (Wave 0 below). Highest
  trust-building value; reviewable in minutes.
- **B. Registration/wiring** (classTags, broker cases, interpreter dispatch,
  CMake) → travels **inside the feature package it wires**. Never separate.
- **C. Behavioral hooks features depend on** (LinearSOE virtuals, Domain
  contact hooks, Element base virtuals, beam `localAxes` responses…) → travel
  **with the first package that needs them**, called out explicitly in the PR
  body as "infrastructure this feature requires".

Reconcile-before-port list (upstream moved under us):
- UmfPack numeric-factorization reuse: upstream #1762 (gaaraujo) vs our ADR-40
  numeric-persist + strategy AUTO change — diff the two, keep upstream's,
  port only what's genuinely missing.
- CreepMaterial is gone upstream; elastic beams gained shear terms; check every
  vanilla file we touch for drift before re-applying hunks.

Flag loudly in PR bodies anything **not byte-identical by default**: H5DRM
z-flip removal + hold-final (jaabell-validated), DirectIntegrationAnalysis
error-return honoring, quad/tri rho serialization (wire-format +1 slot),
GmshRecorder hex20 output format.

## 3. The packages, in order of importance

### Wave 0 — vanilla bug fixes (small PRs, first; ~1 week of porting)
| # | Package | Content |
|---|---|---|
| 0.0 | **TenNodeTetrahedron 6× stiffness fix** | `shp3d` double-applied 1/6 Jacobian factor → all stiffness/reactions 6× too soft (jaabell's own element; found unledgered by the 2026-07-22 audit). Single-hunk, patch-test-verified — the ideal first PR |
| 0.1 | Crash/UB fixes | H5DRM `stuff[12]` identity init (3 parsers); FE_Element subtype-ctor scratch guard; PythonStream `"%s"` format; GmshRecorder hex20 type-17 + permutation; DistributedSuperLU `stat`→`superlu_stat` MSVC collision |
| 0.2 | Serialization fixes | quad/Tri31 element-rho send/recv (wire-format change — flag); missing broker entries audit |
| 0.3 | Error-contract fixes | DirectIntegrationAnalysis + TransientDomainDecomposition return honoring; `Domain::clearAll()` EQ leak; Mumps `-opt` parse guard |
| 0.4 | Registration gaps | Lysmer loader/triangle interpreter registration; H5DRM openseespy dispatch; InitStrainNDMaterial dimension-general; ASDPlasticMaterial3D setResponse labels |
| 0.5 | H5DRM behavior | z-flip removal + hold-final-displacement + tend overrun (jaabell's own patches, validated in the DRM study) |
| 0.6 | Byte-identical perf fixes | ADR-74 harvest, output-identical + suite-gated: MPIDiagonalSOE `setSize` O(N²)-quicksort → `std::sort` (161×); TransformationConstraintHandler hash-membership `handle()` (42×); TransformationDOF_Group redundant SP sweep removal (strip the fork profiler brackets when porting) |

### Wave 1 — small additive features, zero/near-zero vanilla footprint
| # | Package | Content | Deps |
|---|---|---|---|
| 1.1 | Plane-strain σ_zz | `NDMaterial::getStressZZ` + `stressesPlaneStrain` responses | — |
| 1.2 | Beam localAxes | response id 30 on the 10 beam classes | — |
| 1.3 | DDM integrators | LadrunoHHT + LadrunoGeneralizedAlpha (header promotions only) | — |
| 1.4 | Robust statics | LadrunoArcLength (Ramm + STABILIZE), DynamicRelaxation, IndirectControl, StabilizedUnbalance test | — |
| 1.5 | ASDPlastic geomaterials | Hoek–Brown (rock) + StiffSoil shear/cap components for jaabell's ASDPlasticMaterial3D framework (10 new YF/PF/EL/hardening headers + regenerated registries, +2.9k lines, `test_HoekBrown.cpp`). His framework — coordinate directly; found unledgered by the audit | — |

### Wave 2 — core method infrastructure (the fork's flagship)
| # | Package | Content | Deps |
|---|---|---|---|
| 2.1 | Explicit dynamics I | CentralDifferenceLadruno, ExplicitBathe (collapsed class, **no deprecated alias tags**), HRZ mass-conserving lumping (ADR-35, Hinton–Rock–Zienkiewicz 1976), CriticalTimeStep + queryable dt_cr, bulk viscosity | — |
| 2.2 | Mass scaling (explicit II) | The full SMS stack: lumped selective mass scaling (ADR-36, DT2MS-style; `CentralDifferenceSMS` + `-sms` on ExplicitBathe) + **consistent Olovsson SMS** (ADR-38; `…SMSConsistent`/`-sms -consistent`, centroidal M̄ + matrix-free PCG) + the 4 `LinearSOE.h` base virtuals + MPIDiagonalSOE distributed-PCG overrides (V5) + KE_ms energy channel + the classic-Tcl SMS registration. Refs: Olovsson–Simonsson–Unosson 2005; ADR-37 validation battery = the ported tests. **Caveat: the V5 parallel PCG is locally validated but not CI-gated — state it in the PR body or ship the MPI part later** | 2.1 |
| 2.3 | Projection handler | LadrunoProjectionHandler + projector (ADR-30) | 2.1 |
| 2.4 | Finite-strain infra | FiniteStrainNDMaterial base, LogStrainNDMaterial, LogStrain2D, InitDefGrad, StagedStrain | — (can go in Wave 1) |
| 2.5 | Energy + recorders | EnergyBalanceRecorder + kernel/channels, LadrunoRecorder (HDF5), MonitorRecorder | 2.1 (energy), 1.2 (localAxes); HDF5 ≥1.12 |

### Wave 3 — elements & materials catalog
| # | Package | Content | Deps |
|---|---|---|---|
| 3.1 | Solids | LadrunoBrick (+EAS/hourglass), SolidTransformation layer, Brick20, SolidShell, plane family (Quad/CST/LST/CSTPair) | 2.4 |
| 3.2 | J2/steel | LadrunoJ2 (+kernel/hardening/Lemaitre/IMPL-EX), J2Finite, UniaxialJ2, RebarBuckling, BondSlip, CohesiveHinge(+Biaxial) | 2.4 |
| 3.3 | Concrete | LadrunoRCConcrete, RCFiniteStrain, LadrunoConcrete3D (CDPM2-grade) | 3.1 (lch handshake), 2.4 |
| 3.4 | Beams | LadrunoIMKBeam 2d/3d, LadrunoDispBeamColumn 2d/3d (regularization + hinge) | — |
| 3.5 | Coupling/embedded | RBE2/RBE3 couplings, EmbeddedRebar, EmbeddedNode (**minus** experimental UR/UP/D9/corot modes) + Element base virtuals they need | 2.1 (dt_cr virtual) |
| 3.6 | Bézier elements | BezierTri6, BezierTet10 (+corot/F-bar) | 2.4 |
| 3.7 | Modal/eigen | complexEigen + damping assembler, modalResponseHistory, RSA -combine, frequency/random response | — |

### Wave 4 — hold / defer / coordinate with José first
- **Contact + LadrunoTie** — biggest single subsystem, but **serial-only**
  (handler send/recv stubs) and ADR-60 re-emit still draft. Ship only when
  José wants it as-is (serial) or after parallelization.
- **FEAST eigen** — hard MKL dependency; needs an upstream-acceptable
  build-gating story (`LADRUNO_MKL_FEAST` opt-in exists).
- **Porous media (LadrunoUP + PorousOverlay)** — newest, deepest fork coupling;
  let Waves 2–3 land first since it rides on them.
- **LadrunoRigidBody** — explicit-only, serial-only v1.
- **OpenSeesPyMP + build patches** — jaabell has his own OpenSeesPyMP branch;
  reconcile with him rather than PR blind.
- **Profiler, ADR-74 numberer/perf track** — fork-workflow tooling / in-flight;
  numberer isn't even ledgered yet. Later.
- **Splash banner** — fork branding; off by default in packages.

## 4. Per-package workflow (repeat for each)

1. Branch `up/<package>` off `jaabell/ladruño` in the integration worktree.
2. Copy new files; re-apply vanilla hunks onto the drifted base
   (`grep -rn "// Ladruno"` + ledger rows = the hunk list).
3. Scrub: re-stamp headers (no Guppi, references in), rewrite ADR-## comments,
   drop banner/ledger edits, no fork-internal paths.
4. Port the package's test subset; build (`OpenSees`, `openseespy`, MP where
   relevant) and run tests **on the jaabell base**.
5. Single clean commit (team trailers) → PR to `jaabell:ladruño` with the
   documented body (theory, references, verification, vanilla-touch list,
   byte-identical statement) → **squash merge**.
6. Mark ledgers with the upstream PR number.

## 5. Completeness audit (2026-07-22) — every artifact → a package

Method: `git diff e1237189a origin/ladruno` (241 new SRC files + 142 modified
vanilla SRC files) cross-checked against both ledgers, `banner_features.txt`,
the `classTags.h` 33000+ band, and the ADR index. **Result: 100% bucketed.**

> [!warning] Port-list process note
> The `// Ladruno` marker grep is NOT a complete port list — the audit found
> **19 modified vanilla files with no marker** (TenNodeTet 6× fix,
> DistributedSuperLU MSVC fix, the ASDPlastic extension set, profiler-macro
> hooks, `TransientIntegrator.h getCriticalTimeStep` virtual, build wiring).
> The authoritative port list for every package is
> `git diff --name-only e1237189a origin/ladruno`, classified by this table.

### New-file directories → packages
| Artifact | Package |
|---|---|
| `analysis/analysis/LadrunoComplexEigen|DampingAssembler|ModalResponse|ModalCombination` | 3.7 |
| `analysis/handler/LadrunoProjection*` | 2.3 |
| `analysis/handler/LadrunoContact*`, `domain/contact/*`, `domain/constraints/LadrunoTie*` | W4 contact |
| `analysis/integrator/CentralDifferenceLadruno|CriticalTimeStep|LadrunoMassLumping|LadrunoEnergyChannels` (+ ExplicitBathe rework) | 2.1 |
| `analysis/integrator/CentralDifferenceSMS*|LadrunoMassScaling*|LadrunoConsistentRefine` | 2.2 |
| `analysis/integrator/LadrunoArcLength|DynamicRelaxation|IndirectControl|FictitiousMass`, `convergenceTest/LadrunoStabilizedUnbalance` | 1.4 |
| `analysis/integrator/LadrunoHHT|LadrunoGeneralizedAlpha` | 1.3 |
| `analysis/numberer/LadrunoParallelNumberer` | W4 perf (byte-identical harvest → 0.6) |
| `domain/pattern/ladrunoPorousOverlay/*`, `element/ladrunoUP/*`, `recorder/Ladruno_OverlayResults*` | W4 porous |
| `element/bezierTriangle|bezierTetrahedron` | 3.6 |
| `element/ladrunoBrick|solidTransformation|ladrunoSolidShell|ladrunoPlane` | 3.1 |
| `element/ladrunoDispBeamColumn|ladrunoIMKBeam` | 3.4 |
| `element/ladrunoEmbeddedRebar|ladrunoEmbeddedNode|ladrunoDistributingCoupling|ladrunoKinematicCoupling` | 3.5 |
| `element/ladrunoRigidBody` | W4 rigid body |
| `interpreter/PythonMPIModule.cpp` | W4 OpenSeesPyMP |
| `material/nD/ASDPlasticMaterial3D/{HoekBrown,StiffSoil}*` | **1.5 (audit find)** |
| `material/nD/FiniteStrainND*|LogStrain*|InitDefGrad*|StagedStrain*` | 2.4 |
| `material/nD/LadrunoJ2*|LadrunoHardening|LadrunoDamage`, `material/uniaxial/LadrunoUniaxialJ2|RebarBuckling|BondSlip|CohesiveHinge`, `nD/LadrunoCohesiveHingeBiaxial` | 3.2 |
| `material/nD/LadrunoRCConcrete|RCFiniteStrain|RCKernel|LadrunoConcrete3D*` | 3.3 |
| `recorder/EnergyBalance*|LadrunoRecorder|Ladruno_*|LadrunoMonitor*` | 2.5 |
| `system_of_eqn/eigenSOE/Feast*|LadrunoBlockZ*|LadrunoDistBlockZ*|LadrunoFeastInnerSolve` | W4 FEAST |
| `utility/profiler/*` | W4 profiler |

### Modified vanilla files (all 142) → packages
0.0 TenNodeTetrahedron · 0.1 H5DRM parsers ×3 / FE_Element / PythonStream /
GmshRecorder / DistributedSuperLU · 0.2 fourNodeQuad×5 + Tri31 rho ·
0.3 DirectIntegration/TransientDD analyses, Domain clearAll, OpenSeesMiscCommands ·
0.4 TclModelBuilder + Lysmer, OpenSeesPatternCommands, InitStrainNDMaterial,
ASDPlasticMaterial3D.h · 0.6 MPIDiagonalSOE(sort part) /
TransformationConstraintHandler / TransformationDOF_Group · 1.1 NDMaterial +
4 plane-strain materials + quad/tri responses · 1.2 the 10 beam classes ·
1.3 HHT.h / GeneralizedAlpha.h · 1.5 ASDP registries ×9 · 2.1 ExplicitBathe.{h,cpp},
TransientIntegrator.{h,cpp} (getCriticalTimeStep) · 2.2 LinearSOE.h,
MPIDiagonalSOE(hooks), DiagonalSOE.h · 2.5 Node energy/localAxes consumers,
Lysmer/ASDAbsorbing publishers, recorder CMake · 3.1 VTK/VTKHDF/PVD vtktypes ·
3.7 ResponseSpectrumAnalysis, NodeRecorder, Node complex modes ·
W4: Mumps* (FEAST), Umfpack* (reconcile #1762), Matrix/Vector/ID/MovableObject +
algorithm/analysis profiler hooks, ParallelNumberer, classTags/broker/interpreter/
CMake wiring rows (travel with their features), tclMain/PythonModule banner (excluded).

### ADRs with no SRC artifacts (studies — nothing to port)
ADR-42 buckling, 49/49a integrator study, 50 AEM scoping, 51 element removal,
54 FDEM-lite, 55 contact runtime discovery, 59 gradient concrete, 65 dt strategies,
67/68 perf studies, modal_gap_study, apeGmsh scoping docs.

### Fork work explicitly OUT of upstream scope (complete list, so nothing is "forgotten")
Inno Setup installer + `Ladruno_scripts/` build/banner/stamp tooling ·
`Ladruno_tools/` (profiler viewer FastAPI+React, AnalysisLog driver) ·
`tests/` harness (281 files — mined per-package for ported tests, not shipped wholesale) ·
`Ladruno_implementation`/`Ladruno_internal` docs (30+ user guides = source material
for PR documentation) · splash banner · robust-solve Python driver (ADR-31) ·
apeGmsh integration contract · external skill repos.

## 6. Status board — THE campaign memory (update in every session that touches it)

This document is the single source of truth for the upstream campaign. Any
session (human or agent) that ports, opens, merges, or re-scopes a package
**must** update this table and add a dated entry to §7.

| Package | Content (short) | Branch | Upstream PR | Status |
|---|---|---|---|---|
| 0.0 | TenNodeTet 6× stiffness fix | `up/00-tennodetet-shp3d` | [jaabell#29](https://github.com/jaabell/OpenSees/pull/29) | **PR open** (2026-07-22) |
| 0.1 | Crash/UB fixes (H5DRM init, FE_Element, PythonStream, Gmsh hex20, SuperLU MSVC) | — | — | not started |
| 0.2 | quad/tri rho serialization | — | — | not started |
| 0.3 | Error-contract fixes | — | — | not started |
| 0.4 | Registration gaps | — | — | not started |
| 0.5 | H5DRM behavior (z-flip, hold-final) | — | — | not started |
| 0.6 | Byte-identical perf fixes (ADR-74 harvest) | — | — | not started |
| 1.1 | Plane-strain σ_zz | — | — | not started |
| 1.2 | Beam localAxes responses | — | — | not started |
| 1.3 | DDM HHT/GeneralizedAlpha | — | — | not started |
| 1.4 | Robust statics (ArcLength/DR/IndirectControl) | — | — | not started |
| 1.5 | ASDPlastic Hoek–Brown + StiffSoil | — | — | not started (coordinate w/ José) |
| 2.1 | Explicit dynamics I (CD-Ladruno, ExplicitBathe, HRZ, dt_cr) | — | — | not started |
| 2.2 | Mass scaling (SMS lumped + consistent Olovsson + LinearSOE virtuals) | — | — | not started |
| 2.3 | Projection handler (ADR-30) | — | — | not started |
| 2.4 | Finite-strain material infra (LogStrain, staged wrappers) | — | — | not started |
| 2.5 | Energy + recorders (EnergyBalance, .ladruno HDF5, Monitor) | — | — | not started |
| 3.1 | Solids (Brick/EAS, SolidTransformation, Brick20, SolidShell, plane family) | — | — | not started |
| 3.2 | J2/steel materials | — | — | not started |
| 3.3 | Concrete (RCConcrete, Concrete3D) | — | — | not started |
| 3.4 | Beams (IMK, regularized DispBeamColumn) | — | — | not started |
| 3.5 | Coupling/embedded (RBE2/RBE3, embedded rebar/node) | — | — | not started |
| 3.6 | Bézier elements | — | — | not started |
| 3.7 | Modal/eigen (complexEigen, modal family) | — | — | not started |
| W4 | contact+tie / FEAST / porous / rigid body / OpenSeesPyMP / profiler / numberer | — | — | deferred (decide w/ José) |

## 7. Decision & session log (append-only, newest first)

- **2026-07-22 — package 0.0 shipped as [jaabell#29](https://github.com/jaabell/OpenSees/pull/29).**
  First upstream PR of the campaign. Port mechanics validated end-to-end:
  worktree `up-00-tennodetet` off `jaabell/ladruño`, single clean commit,
  author Nicolas Mora Bowen `<nmorabowen@gmail.com>`, trailers
  `Co-authored-by: Patricio Palacios <ppalacios92@gmail.com>` +
  `Co-authored-by: Jose A. Abell <jaabell@uandes.cl>` (emails harvested from
  git history — reuse these for all packages), zero AI traces. Physics note
  for the PR narrative: K and M scaled down TOGETHER, so eigenfrequencies were
  ~unchanged by the bug — it only shows in absolute response (displacements 6×
  too large, mass/weight 6× under-counted); don't claim frequency shifts.

- **2026-07-22 — mass-scaling placement confirmed.** The whole mass-scaling
  stack upstreams as package 2.2 (lumped ADR-36 + consistent-Olovsson ADR-38 +
  the LinearSOE PCG virtuals), with HRZ lumping riding 2.1 and the KE_ms
  recorder channel in 2.5. Known soft spot: V5 distributed PCG not CI-gated.
  The porous-overlay SMS composability stays with W4 porous.
- **2026-07-22 — completeness audit done.** All 241 new + 142 modified SRC
  files bucketed (§5). Found + ledgered: TenNodeTet 6× fix (→0.0),
  DistributedSuperLU MSVC fix (→0.1), ASDP Hoek–Brown/StiffSoil (→1.5).
  Port lists come from the merge-base diff, NOT the `// Ladruno` grep
  (19 unmarked files).
- **2026-07-22 — campaign plan created.** Target `jaabell/ladruño`
  (= upstream master 2026-07-13). Fresh-branch ports, squash merges, team-only
  authorship (no AI traces, no Guppi in headers), documentation-with-references
  required, Waves 0→4.
