# WP-120 — Code-health survey: duplication and dead code in fork-authored sources

Revision 1. Not yet adversarially reviewed. One Opus subagent built the §4 reuse inventory (read-only); every
claim from it that this doc repeats was re-checked by running the command (see "Verification of delegated claims").

Status: **survey complete; PR #855 marked ready 2026-09-26. Read-only — no production code changed.** The
survey's follow-ups are under way as their own WPs: WP-122 (#857, merged) for R1; WP-123 (#858, refactor
candidate 1); WP-124 (#859, candidate 2); WP-126 (#861, open question 1). The owner merges.

Scoped 2026-09-25. Branch `wp/120-code-health-survey`, cut from `ladruno` @ `bc5c33453`. This is
Phase 0 ("measure before building") of the `agent-surface` method, applied to code health. Tooling:
[`wp120_code_health/`](wp120_code_health/README.md); raw outputs in `wp120_code_health/results/`.

## Problem

WP-115…119 (merged 2026-09-24) found that bugs multiplied through copies. The Rayleigh-snapshot defect lived
in parallel copies of LadrunoQuad/CST/LST/CSTPair (#562), and again in BezierTet10/Tri6, LadrunoIMKBeam/2d and
LadrunoDispBeamColumn2d/3d (hazards converted in WP-115/116). The ground-motion double subtraction lived in
three near-identical vanilla beams (#854). Every fix had to be applied N times. The founding entry is
`LEDGER_quirks` "`getResistingForceIncInertia` MUST snapshot the shared static `resid`".

Before building any refactor or gate, this WP measures: how much fork code is duplicated, where, whether the
clusters have a history of fixes replicated across copies, how much dead code exists, and which shared kernels
already exist to reuse.

## Shape

1. **Scope.** Fork-authored = files carrying `LADRUNO-HEADER-START` (237 files; the set the quirk lint scans).
   Cross-checked against git: `SRC` files on `ladruno` that are neither on `upstream/master` nor at the merge-base
   (`inventory.py`).
   *Accept:* the drift between the stamp, `stamp_headers.py` GLOBS and git is listed file by file.
2. **Duplication.** Neither PMD CPD (no Java) nor jscpd is installed, so `clones.py` is a dependency-free
   re-implementation of CPD's method. It reuses the quirk lint's `clean()`/`functions()` rather than writing a
   second C++ scanner. It runs a rolling hash over 100-token windows (CPD's default), extends matches to maximal
   length, and verifies every pair token by token. Two modes: `exact` (type-1) and `norm` (identifiers and
   literals normalised, type-2). Scopes: fork, vanilla (baseline), and cross (fork code copied from vanilla).
   *Accept:* on the tree just before the #562 fix (`4a975edee`) it must group the four plane elements, with
   `getResistingForceIncInertia` among the cloned functions. **Passed** (a 7-file family: Quad, CST, LST,
   CSTPair, Brick, BezierTri6, BezierTet10; the function is listed).
3. **Bug history.** `history.py` takes each clone family and finds (a) first-parent PRs on `ladruno` that touch
   ≥ 2 members and whose title reads as a fix, and (b) `LEDGER_quirks` sections naming ≥ 2 member classes. Every
   candidate PR was then **hand-classified** from its diff's hunk functions: *replicated defect* (the same defect
   fixed in ≥ 2 copies) or *co-edit* (feature or docs).
   *Accept:* a family counts as a refactor candidate only if it has ≥ 1 hand-verified replicated defect.
4. **Dead code.** `deadcode.py`: (a) unbuilt files; (b) functions never referenced; (c) always-false guards and
   never-defined `#ifdef`s; commented-out statements. Fork and vanilla are counted separately. (d) Coverage:
   none is available (below).
5. **Reuse.** An inventory of shared kernels and seams, their consumers, and re-implementations of what a kernel
   already provides.

## Results

### R1 — Scope: the stamp misses 31 fork files

| Set | Files | Note |
|---|---:|---|
| Fork-added per git (not on `upstream/master`, not at merge-base `e1237189a`) | 266 | |
| Stamped (`LADRUNO-HEADER-START`) = lint scope = this survey's scope | 237 | 235 fork-added + `ExplicitBathe.{h,cpp}`, which upstream has (Abell, 2026-01-15) and the fork stamped after heavy edits |
| **Fork-added but unstamped (invisible to the quirk lint)** | **31 (8,042 lines)** | list in `wp120_code_health/unstamped_fork_files.txt` |
| Stamped but missing from `stamp_headers.py` GLOBS | 10 | `LadrunoAutoPenaltyReduce.*`, `LadrunoContactAbort.*`, `LadrunoSolverQuery.h`, `LadrunoCohesiveHinge{,Biaxial}.*`, `LadrunoParallelBuild.cpp`. A re-stamp would not maintain them |
| GLOBS entries matching no file | 5 | `ExplicitBathe{SMS,SMSConsistent,LNVD,LNVDSMS,LNVDSMSConsistent}.*`, deleted by #419 |

The 31 unstamped files include shared seams: `SRC/element/LadrunoMassCache.h` (5 element includers),
`LadrunoResponseTokens.h` (17), `LadrunoDamage.h`, `LadrunoThreads.{h,cpp}`, `Ladruno_mutation.h`,
`CriticalTimeStep.{h,cpp}`, `LadrunoHHT`, `LadrunoGeneralizedAlpha`, `LadrunoLoadControl`,
`LadrunoParallelNumberer`, `PythonMPIModule.cpp`, `DRMHigherOrderNode.h`, and 13 ASDPlasticMaterial3D kit headers
(Hoek–Brown, StiffSoil, MohrCoulombTensionCutoff). This is the WP-116 failure mode (a fork class the lint never
saw), ×31. Including them changes the duplication totals by < 0.2 points (R2); the gap matters for the lint, not
for this survey.

### R2 — Duplication headline

| Metric (W = 100 tokens) | Fork, 237 stamped | Fork + 31 unstamped | Vanilla, 3,333 files |
|---|---:|---:|---:|
| Code lines (non-blank after cleaning) | 77,053 | 80,867 | 867,749 |
| Exact clones (type-1), lines | **11,159 (14.5 %)** | 11,803 (14.6 %) | 331,946 (38.3 %) |
| Renamed clones (type-2), lines | **20,388 (26.5 %)** | 21,510 (26.6 %) | 538,626 (62.1 %) |
| …of which cross-file only (exact / type-2) | 7,519 (9.8 %) / 13,775 (17.9 %) | | |
| Fork lines copied from vanilla (exact) | 5,567 (7.2 %) | 6,045 (7.5 %) | — |
| Clone pairs (exact) | 553 | 576 | 42,517 |

Window sensitivity (exact, fork): W = 50 → 26.2 %; W = 100 → 14.5 %; W = 200 → 5.0 %. The family ranking is
the same at every window. Fork code is duplicated at about 40 % of vanilla's rate.

**Trend** (same tool, stamped scope, `git archive` of the first-parent commit before each date):

| Tree | Files | Exact | Type-2 |
|---|---:|---:|---:|
| 2026-06-01 `6fe4bb6ba` | 0 | — (stamp not yet introduced) | — |
| 2026-07-01 `6a5089382` | 152 | 14.8 % | 27.7 % |
| 2026-07-12 `4a975edee` (pre-#562) | 189 | 14.6 % | — |
| 2026-08-01 `d3c18d2b3` | 217 | 13.1 % | 24.5 % |
| 2026-09-01 `a551f4ed5` | 229 | 12.8 % | 24.1 % |
| 2026-09-25 `bc5c33453` | 237 | 14.5 % | 26.5 % |
| …excluding LadrunoDispBeamColumn2d/3d | 233 | 12.6 % | |

The September jump is a **scope change, not new copying**: WP-116 stamped LadrunoDispBeamColumn (2,003 duplicated
lines of 4,405). Without it the rate is flat at about 12.6–12.8 %. So any ratchet must be keyed per file (R8).

### R3 — Clone clusters × replicated-fix history

"Replicated defects" are hand-verified from each PR's diff (which hunk functions it changed in which member files).

| # | Family (copies) | Shared lines, exact / type-2 (sum of pairs) | Replicated defects (hand-verified PRs) | Quirk sections naming ≥ 2 members | Verdict |
|---|---|---:|---|---:|---|
| 1 | **Continuum element shells** (8): LadrunoQuad, CST, LST, CSTPair, LadrunoBrick, Brick20, BezierTri6, BezierTet10 | 1,002 / 2,365 | **8 PRs:** #562 Rayleigh P-clobber in `getResistingForceIncInertia` ×4 · #228 `getInitialStiff` not returning the cached `*Ki` ×2 (Quad, CST) · #670 mass cache not invalidated in `recvSelf` ×8 (+ SolidShell) · #683 EAS inner-Newton warning spam ×2 (Brick, Quad `formEAStrue`) · #588 EAS degeneracy guard blind to axis collapse ×2 (`buildEAStrue`) · #709 unsymmetric plastic tangent symmetrised ×2 (Bezier) · #224 `setParameter` not forwarded to GP materials ×2 (Bezier) · #852 ground-motion sign ×2 (Bezier) | 22 | **Refactor candidate (1st by evidence)** |
| 2 | **Coupling / embedded elements** (4): LadrunoDistributingCoupling, KinematicCoupling, EmbeddedNode, EmbeddedRebar (+ their 4 `OPS_` parsers) | 44–129 per pair / 594 (+ parsers 455) | **Implicit-transient crash (Rayleigh ignored but no `getDamp` override) fixed twice in 2 days:** #219 (RBE3, 2026-06-08) → #220 (EmbeddedNode + EmbeddedRebar, 2026-06-09); #221 (RBE2) was born with the override pre-applied. The 4 bodies are identical today (`C0->Zero(); return *C0;`) | 7 (incl. "overrides `setRayleighDampingFactors` … makes 11 `Element` methods dead", all 4 named) | **Refactor candidate (cheapest)** |
| 3 | **Explicit integrators** (5): CentralDifferenceLadruno, CD-SMS, CD-SMSConsistent, ExplicitBathe, LadrunoDynamicRelaxation | 206 / 359 | #394 W1-E3a `-cflAbort/-recompute` under SMS ×4 (parser + `domainChanged`) · #468 Rayleigh Δt sizing ignoring `betaKinit/betaKcomm` (CDL, Bathe) and a NaN breaker blind to NaN, because `pNorm(0)` skips NaN (CDL, DR) · #472 `recvSelf` on a live scaled object leaking the injected ΔM ×2 (CD-SMS, Bathe) · #475 SMS guards (mostly absorbed **once** by the `LadrunoMassScaling.h` kernel — evidence the kernel pattern works) | 9 | **Candidate, lower priority**: CD-SMS/SMSC already subclass CDL; #419 already collapsed ExplicitBathe 6 → 1 |
| 4 | LadrunoIMKBeam / IMKBeam2d (2) | 122 / 393 | #853 `betaKc` frozen (`commitState` skipped `Element::commitState`) ×2; WP-115 Rayleigh hazard ×2 (hazard, not a defect) | 0 | Watch. Hinge math is already shared (`LadrunoIMKHinge.h`); the element shells are not |
| 5 | LadrunoDispBeamColumn2d / 3d (2) | 771 / 1,369 (max pair) | none in fork history (WP-116 conversion ×2 = hazard). 48–56 % copied from vanilla `DispBeamColumn2d/3d`; **no upstream change to those since the copies (2026-06-16/17)** | 0 | Not a candidate. Track upstream drift instead |
| 6 | Materials (6): LadrunoJ2, J2Finite, UniaxialJ2, RCConcrete, RCFiniteStrain, Concrete3D | 642 / 951 (RC ↔ RCFinite alone 337: `setResponse/getResponse/Print/sendSelf`) | none | 3 (feature-level) | Not a candidate |
| 7 | Ladruno recorder ← vanilla `MPCORecorder.cpp` | 1,536 copied | none | 0 | Not a candidate: byte-faithful by design (1e-12 parity gate against the frozen recorder) |
| 8 | Plane `OPS_` parsers (4); `SolidTransformation{Linear,Hypo,Finite}`; `LogStrain`/`InitDefGrad`/`StagedStrain` wrappers; SANISAND wrappers | 236; 144; 163; 39 | 0 / feature co-edits only | 0–4 | Not candidates |
| 9 | Intra-file: `LadrunoContactFE.cpp` (608 lines: 2D/3D mortar + friction variants), `Ladruno_NodeResults.cpp` (549), `LadrunoModalResponse.cpp` (380), `LadrunoDispBeamColumn3d.cpp` (376) | — | no replicated fix identified (ContactFE's 3 fix commits not checked hunk-by-hunk for 2D↔3D pairs) | — | Not candidates |

What family 1's copies share is the **`Element`-contract shell**, not the physics. The plane cores already
share `LadrunoFiniteStrain2DKernel.h`. The replicated fixes all sit in six places:
`getResistingForceIncInertia` (f − Q + M·a + Rayleigh, snapshot order), `addInertiaLoadToUnbalance` (sign and
mass source), the cached `getInitialStiff`, mass-cache invalidation on `recvSelf`, `setParameter` forwarding to
GP materials, and the `setResponse` token ladder.

### R4 — Dead code

| Category | Fork (237 files, 120,794 raw lines) | Vanilla (3,333 files, 1,360,671 raw lines) |
|---|---|---|
| (a) `.cpp` in no CMake file | 1: `profiler_selftest.cpp` (a standalone harness) | 118 |
| (a) `.h` included by nothing and in no CMake file | 0 | 23 |
| (b) functions never referenced (hand-checked) | **37** names, 51 definitions, 157 lines; plus 12 more called only from `tests/*.cpp` harnesses (test API, legitimate) | not surveyed |
| (c) `const bool X = false` guard flags | 0 | 6 flags guarding 111 uses: `ManzariDafalias::debugFlag` 53 live (61 raw, incl. commented-out), `PM4Sand`/`PM4Silt` 19 each, `SAniSandMS` 15, `cms` 3, `verbose` 2 |
| (c) `#define X 0` flags in use | 0 | 12 flags, 37 uses |
| (c) `if (0)` / `if (false)` / `#if 0` lines | 0 / 0 | 5 / 55 |
| (c) `#ifdef X` with X defined nowhere | 3 (deliberate opt-ins: `LADRUNO_J2FINITE_CHANNELB_NUMERIC` ×2, `LADRUNO_USE_SWMR` ×1) — never compiled by any build | 95 |
| Commented-out statements (`// f(x);`, `// a = b;`) | 34 (0.28 / KLOC); 16 sit in the two DispBeamColumn copies (11 of the same lines exist in vanilla `DispBeamColumn2d/3d`) | 6,398 (4.70 / KLOC) |
| (d) code no test executes | **not measurable:** no Zone-A coverage artifact exists (`ladruno.yml` has no gcov/lcov step) | — |

The 37 never-referenced functions are mostly accessors never called: 9 `FeastEigenSOE::get*`, 4
`LadrunoContactDomain::getNum*` + `buildAdapterCount`, `numFramesWritten`, `lastSolveSeconds`. Five are worth a look (`deadcode.py` output lists
all 37):

- **`requiresPartitionReduction()`** (`Ladruno_ResultIO.h:82`; 8 definitions across the recorder sources). The header comment
  (`Ladruno_NodeResults.h:186`) says the reaction sources need "a per-step partition reduction … before sink
  accumulation". **Nothing calls it.** Either the reduction is done unconditionally elsewhere, or reactions at
  partition-boundary nodes are not summed under OpenSeesMP. Not verified; see Open questions.
- **`frameTimeVarying()`** (×4 in the SolidTransformation seam) and **`getJ()`** (×5). A seam contract with no
  consumer.
- **`linearizePair`** (`LadrunoMortarKernel.h:511`) and **`lateralResidual`** (`LadrunoConcrete3DKernel.h:625`).
  Unreached kernel code.

On dead code the fork is clean: no dead guards, no `#if 0`, almost no commented-out code. The dead-guard
pattern (`debugFlag`) is a vanilla phenomenon that the fork inherits by subclassing (LadrunoSANISAND derives from
ManzariDafalias).

### R5 — Reuse: existing shared kernels

| Kernel / seam | Provides | Fork consumers |
|---|---|---|
| `ladrunoPlane/LadrunoFiniteStrain2DKernel.h` | 2D F, F-bar, finite internal force/tangent (`ladruno_fs2d::`) | Quad, CST, LST, CSTPair |
| `solidTransformation/SolidTransformation.h` (+Linear/Corot/Finite/Hypo) | geometry seam for solids | BezierTet10, LadrunoBrick, LadrunoUP |
| `solidTransformation/LadrunoHypoKernel.h` | `ladruno_hypo::` det3, inv3, polar3, objective increments | LadrunoBrick, LadrunoUP |
| `ladrunoIMKBeam/LadrunoIMKHinge.h` | series-hinge solve | IMKBeam, IMKBeam2d |
| `ladrunoEmbeddedRebar/LadrunoEmbeddedKernel.*` | gap/penalty coupling assembly, Δt_cr | all 4 coupling/embedded elements |
| `ladrunoUP/LadrunoUPKernel.h` + `LadrunoUPShapes.h` | Q/H/S blocks; T3/Q4/H8/BT6/BTET10 shapes + GP tables | LadrunoUP, LadrunoPorousOverlay |
| `ladrunoBrick/LadrunoHex20Shape.h` | hex20 shape/B/consistentMass | Brick20 only |
| `domain/contact/Ladruno{ContactProjection,ContactKernel,Contact2DKernel,EdgeKernel,MortarKernel,FrictionKernel,ContactBucketSort}.h` | NTS/mortar/friction/broad-phase | ContactFE, ContactHandler, LadrunoTie, ContactDomain (FrictionKernel: ContactFE only) |
| `material/nD/LadrunoJ2Kernel.h`, `LadrunoRCKernel.h`, `LadrunoConcrete3DKernel.h`, `LogStrainKernel.h`, `LadrunoHardening.h` | return maps, backbones, Hencky/log-strain, Voce/linear hardening | J2+J2Finite; RC+RCFinite; Concrete3D only; J2Finite+RCFinite+LogStrain; J2Kernel+UniaxialJ2 |
| `material/LadrunoMaterialStatus.h` | `LADRUNO_MATERIAL_REFUSED`, refusal counter | LadrunoBrick, LadrunoSANISAND (+ vanilla Domain, ASDPlastic) |
| `analysis/integrator/LadrunoMassLumping.h` | HRZ lumping with guards | Brick20 only (+ vanilla CriticalTimeStep) |
| `analysis/integrator/LadrunoMassScaling.h`, `LadrunoConsistentRefine.h`, `LadrunoFictitiousMass.h`, `LadrunoEnergyChannels.h` | SMS, consistent refine, Gershgorin diagonal, energy registry | CD-SMS/SMSC/Bathe; SMSC/Bathe; ArcLength/DR; Bathe + EnergyBalance |
| `element/LadrunoMassCache.h` (**unstamped**) | per-instance mass cache | BezierTet10/Tri6, LST, Quad, SolidShell |
| `element/LadrunoResponseTokens.h` (**unstamped**) | response-token table | 17 includers |

**Re-implementations of what a kernel already provides** (verified lines):

| Duplicated | Where | Kernel that has it | Drop-in? |
|---|---|---|---|
| 3×3 inverse/det, row-major `[9]` | `LadrunoBrick.cpp:1479 invert3x3`, which **already includes `LadrunoHypoKernel.h`**; `BezierTet10.cpp:883`; `LogStrainKernel.h:83`; `LadrunoUPShapes.h:98`; `SolidTransformationCorot.cpp:57`; `InitDefGradNDMaterial.cpp:62` | `ladruno_hypo::inv3` (`LadrunoHypoKernel.h:96`) | Brick/Tet10/LogStrain/UP yes; Corot and InitDefGrad differ on singular input |
| 3×3 inverse, `[3][3]` | `LadrunoHex20Shape.h:205` ≡ `LadrunoSolidShell.cpp:403`; inline in Brick, Tet10, Concrete3DKernel | none | — |
| Symmetric 3×3 Jacobi eigen | `LadrunoRCKernel.h:312`, `LadrunoConcrete3DKernel.h:640`, `LogStrainKernel.h:96`, `LadrunoDistributingCoupling.cpp:51` (4 copies) | none shared | no (sorting and layout differ) |
| Polar decomposition | `SolidTransformationCorot.cpp:138`; `LadrunoHypoKernel.h:129` (its comment at `:126` names the Corot algorithm as the source) | HypoKernel | return type differs |
| `dot3`/`cross3`/`norm3` | `LadrunoContactProjection.h`, `LadrunoEdgeKernel.h`, `LadrunoFrictionKernel.h` | ContactProjection | yes |
| Gauss tables and shape functions (Q4, H8, T6, Tet10) | LadrunoQuad (×2 in its two constructors), LST, Brick, SolidShell (×4), Bezier*, Ladruno_ElementResults.h | `LadrunoUPShapes.h` (used only by UP) | no (node order and layout differ) |
| HRZ lumping | `LadrunoLST.cpp:344-374` inline | `LadrunoMassLumping.h::hrzLump` | no (kernel needs full M + direction table) |
| `rowsum\|diagonal\|hrz` lumping-flag parser | CDL, CD-SMS, CD-SMSC, ExplicitBathe (×3) | none | — |
| zero `getDamp` | the 4 coupling/embedded elements, identical | `LadrunoEmbeddedKernel` (all 4 already include it) | would be |
| `addInertiaLoadToUnbalance` | 14 fork elements, 2 variants (diagonal `M(i,i)`: Quad/CST/LST, where mass is lumped so it is exact; full `M·ra`: Tri6/SolidShell/Brick20) | none | — |

### Verification of delegated claims

The subagent's claims that this doc repeats were re-run: `LadrunoBrick.cpp:68` includes `LadrunoHypoKernel.h`
and `:1479` defines its own `invert3x3`; `LadrunoMassLumping.h`'s only fork includer is `LadrunoBrick20.cpp`;
`LadrunoFrictionKernel.h`'s only `.cpp` includer is `LadrunoContactFE.cpp`; the four `getDamp` bodies are
identical (`rg -A4 '::getDamp\('`); two separate `jacobi3` definitions (LogStrainKernel, DistributingCoupling)
plus the RC and Concrete3D eigen routines; and every `file:line` in the R5 re-implementation table (`sed -n`
on each). Two of its claims were corrected on re-run: the NaN-breaker fix in #468 touched CDL and DR, not
Bathe; and `getJ` has 5 definitions, not 3. Claims not re-run are not repeated here.

## Recommendations (for the owner — nothing here is started)

**(i) Refactor WPs, only where a replicated defect is proven.** Each must be proven bit-identical, like the
WP-115 conversions (fingerprint the recorded outputs before and after; keep operation order), and must add a
test that fails on a one-line mutation of the shared code.

1. **Coupling/embedded "undamped constraint element" helper (smallest, do first).** Move the Rayleigh-ignoring
   trio (`setRayleighDampingFactors` / `getDamp` / `getRayleighDampingForces`) into `LadrunoEmbeddedKernel`,
   which all four already include. Evidence: #219 → #220 → #221 in 2 days. About 40 lines; low risk.
2. **Continuum element-shell helpers (largest payoff).** Free functions (not a new base class — the eight
   elements have different bases and DOF layouts) for the dynamic residual (`f − Q + M·a + Rayleigh` with
   snapshot and a fixed operation order), ground-motion `Q` accumulation, cached initial stiffness, mass-cache
   invalidation, and `setParameter` forwarding. Evidence: 8 replicated-fix PRs across 8 elements. Stage it:
   plane four first (they already share a kernel), then Bezier, then Brick/Brick20.
3. **Explicit integrator guards (medium).** One helper for damped Δt sizing (all three βK slots), the NaN/Inf
   breaker, the lumping-flag parser, and SMS inject/restore on `recvSelf`. Share it between CDL and ExplicitBathe
   (independent `TransientIntegrator`s). Evidence: #394, #468, #472; #475 shows the pattern working once it lives
   in a kernel. Do **not** redo #419's class collapse: CD-SMS/SMSC already subclass CDL.

Not recommended: DispBeamColumn2d/3d, the material family, the recorder↔MPCO copy, parsers. None has replicated-bug
history. IMK beams: fold in with (2) only if a second replicated defect appears.

**(ii) A CI ratchet: yes, but per-file and advisory-first.** The totals move with scope (R2: WP-116's stamp
alone added 2,003 duplicated lines), so a global-percentage ratchet would fire on hygiene work. Proposal:
commit `baseline.json` (per fork file: exact cross-file duplicated lines, never-referenced function count). In
the `static-gates` job, **fail only if a file already in the baseline increases**. New files are reported, not
failed. A PR that lowers a count updates the baseline. Runtime is about 12 s (fork scope). Never an absolute gate.
Add the stamp-drift check from `inventory.py` too (fork-added-but-unstamped count must not increase): it is the
one number here that directly limits the quirk lint.

**(iii) "Reuse before writing" — proposed addition to AGENTS.md and the `ladruno-new-element` guide:**

> Before writing a helper, check these (grep the function name):
> 3×3 det/inv/polar → `ladruno_hypo::` in `solidTransformation/LadrunoHypoKernel.h`;
> 2-D finite strain / F-bar → `ladrunoPlane/LadrunoFiniteStrain2DKernel.h`;
> shapes + Gauss tables (T3/Q4/H8/T6/Tet10) → `ladrunoUP/LadrunoUPShapes.h`;
> HRZ lumping → `analysis/integrator/LadrunoMassLumping.h`; per-instance mass cache → `element/LadrunoMassCache.h`;
> Voce/linear hardening → `material/nD/LadrunoHardening.h`; Hencky/log strain, `jacobi3` → `LogStrainKernel.h`;
> material refusal → `material/LadrunoMaterialStatus.h`; `dot3/cross3/norm3` → `LadrunoContactProjection.h`;
> SMS / consistent refine → `LadrunoMassScaling.h` / `LadrunoConsistentRefine.h`.
> A new element's `getResistingForceIncInertia` / `addInertiaLoadToUnbalance`: copy the current plane-family
> shape (post-#562, post-#852), never an older sibling's.

## Rejected approaches

- **Install jscpd via `npx` or PMD.** Neither is installed; running an unvetted downloaded package in the dev
  session was avoided. The CPD algorithm is small enough to re-implement, and reusing the lint's C++ cleaner keeps
  one scanner in the repo.
- **Count "fix PRs touching ≥ 2 members" as evidence directly.** The regex over PR titles has false positives
  (#787 "docs+test", #650 feature). Every counted replicated defect was classified by hand from the hunk functions.
- **Use per-file line counts to call a function dead.** The first heuristic ("≤ 2 occurrences") flagged 203
  functions, mostly helpers defined once and called once. It was replaced by per-occurrence classification
  (definition header / prototype / use), then a whole-repo grep that moved 12 to "test-only API".
- **A global duplication-percentage ratchet.** Refuted by the trend table: scope changes move it by ~2 points.
- **Collapse the CentralDifference classes like #419.** CD-SMS/SMSC already inherit CDL; the replicated fixes
  were in parsers and sizing/guard code, which a shared helper covers without a class-tag change.
- **Refactor DispBeamColumn2d/3d (the largest single pair, 771 lines).** No replicated-bug history. Size alone
  is not evidence.

## Open questions

1. **`requiresPartitionReduction()` has no caller.** Are reaction forces at partition-boundary nodes summed in
   the Ladruno recorder under OpenSeesMP? A 2-partition run with a reaction recorder on a shared node would
   answer it. This may be a live bug; it deserves its own WP rather than a note.
   **→ It was a live bug: WP-126 (#861).** Reproduced 2026-09-25. apeGmsh's stitch returned (0, 10) where the
   serial reaction is (20, 30). The fix is in progress (a `PARTITION_REDUCTION` attribute plus an apeGmsh PR).
2. ~~**Stamp the 31 fork files?**~~ **Done: WP-122 (#857, merged 2026-09-26).** All 31 files are stamped (the
   change is comment-only, proven by fingerprint); GLOBS is +41/−5. The lint found 0 issues, because no rule
   has anything to read in these files yet. `stamp_headers.py --check` is now a CI gate, including the
   upstream-manifest and dead-entry rules.
3. **Coverage.** A gcov/lcov build of Zone-A (Ubuntu) would answer 3(d) and show which of the 37 unreferenced
   functions and the 3 `#ifdef` opt-ins are untested as well as uncalled. Cost: one more CI build (~2× build time).
4. **ContactFE 2D/3D intra-file clones (608 lines).** Not checked hunk-by-hunk for replicated 2D↔3D fixes.
5. **Phase 6 re-measure.** Re-run `clones.py`, `history.py` and `inventory.py` after ~10 WPs; compare against
   R2/R3.
