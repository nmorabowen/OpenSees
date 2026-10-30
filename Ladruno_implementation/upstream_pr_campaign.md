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
   - Commit author: Nicolas Mora Bowen `<nmorabowen@gmail.com>`. Squash commit /
     PR-branch commits carry — **canonical emails, confirmed by Nicolas 2026-07-22**:
     `Co-authored-by: Patricio Palacios <pxpalacios@miuandes.cl>`
     `Co-authored-by: Jose A. Abell <jaabell@miuandes.cl>`
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
   - Nothing from `.claude/`, `CLAUDE.md`, `AGENTS.md`, `ci/check_quirk_patterns.py`,
     `ci/check_viewer_ledger.py`, `Ladruno_implementation/`,
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
6. **Upstream writing rules (José's request, 2026-10-06).** PR bodies, commit
   messages and code comments follow the office tone: impersonal, declarative,
   numbers first (defect, reproduction, value before / after, behaviour change,
   test). No internal references anywhere: no `// Ladruno` markers, ADR-##,
   WP-###, F##, fork PR numbers, LEDGER/quirk names, TIMs. No "we"/"you", no
   em-dash pauses, no bold lead-ins. Model commit: `c6be72ca1` on `ladruño`.
7. **How José lands our PRs.** He does not merge the PR: he cherry-picks each
   commit (authorship kept) onto a `fix/<slug>` branch cut from current upstream,
   merges that into `ladruño`, and closes our PR with a comment. So one commit
   per independently droppable fix, and no fixups after review: amend instead.
8. **Building his base on Windows needs the build fixes** (package 0.13). Until
   he merges them, the integration worktree carries them as uncommitted local
   edits; they never ride inside a bug-fix PR.


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
| 0.0 | **TenNodeTetrahedron: 6× stiffness fix + `getResponse` heap overrun** | Two independent defects in jaabell's own element, both single-hunk and both still live upstream. **(a)** `shp3d` double-applied 1/6 Jacobian factor → all stiffness/reactions 6× too soft (found unledgered by the 2026-07-22 audit); patch-test-verified. **(b)** `getResponse`'s `static Vector stresses(6)` is written 4 GP × 6 = **24 doubles** by both the `stresses` and `strains` branches — a 144-byte heap overrun per recorder step (unchecked `Vector::operator()` outside `_G3DEBUG`); fix to `6*NumGaussPoints`, which is what `setResponse` already advertises. Copy-paste slip from `FourNodeTetrahedron` (1 GP). ⚠ **(a) and (b) canNOT ship together** — (a) is ALREADY upstream as [jaabell#29](https://github.com/jaabell/OpenSees/pull/29) (open since 2026-07-22, scoped to the `shp3d` hunk alone). (b) was found later (fork PR [#692](https://github.com/nmorabowen/OpenSees/pull/692), 2026-08-04) and needs its OWN package **0.0b** — do not widen #29's scope under review |
| 0.1 | Portability & latent-crash fixes | **SCOPED DOWN to 3 pure non-behavioral fixes** (shipped): FE_Element subtype-ctor scratch guard; PythonStream `"%s"` format; DistributedSuperLU `stat`→`superlu_stat` MSVC collision. *H5DRM `stuff[12]` init moved to 0.5 (H5DRM is build-gated + its work is behavioral); GmshRecorder hex20 moved to a recorder-fixes PR (needs the type-17 permutation table)* |
| 0.2 | Serialization fixes | quad/Tri31 element-rho send/recv (wire-format change — flag); missing broker entries audit; **+ GmshRecorder hex20 type-17 + mid-edge permutation (moved here from 0.1)** |
| 0.3 | Error-contract fixes | DirectIntegrationAnalysis + TransientDomainDecomposition return honoring; `Domain::clearAll()` EQ leak; Mumps `-opt` parse guard |
| 0.4 | Registration gaps | Lysmer loader/triangle interpreter registration; H5DRM openseespy dispatch; InitStrainNDMaterial dimension-general; ASDPlasticMaterial3D setResponse labels |
| 0.5 | H5DRM (all changes together) | `stuff[12]` identity init (3 parsers, pure fix) **+** z-flip removal + hold-final-displacement + tend overrun (jaabell's own patches, validated in the DRM study). All H5DRM edits ship as one PR since they touch the same build-gated subsystem |
| 0.6 | Byte-identical perf fixes | ADR-74 harvest, output-identical + suite-gated: MPIDiagonalSOE `setSize` O(N²)-quicksort → `std::sort` (161×); TransformationConstraintHandler hash-membership `handle()` (42×); TransformationDOF_Group redundant SP sweep removal (strip the fork profiler brackets when porting) |
| 0.7 | `CorotCrdTransf3d` static-`T` cross-element aliasing | `CorotCrdTransf3d::T` is **`static`** (`.h:135`), shared by every instance, and `DispBeamColumn3d::getInitialStiff` never calls `crdTransf->update()` first — so `getInitialGlobalStiffMatrix` triple-products through whatever `T` the LAST element to update left behind (ADR-76 §3 audit finding). Fix = make `T` a member (or recompute locally in the initial-stiff path). ⚠ Scope the behavioral claim honestly in the PR: it changes `getInitialStiff()` answers on multi-element corot3d models (wrong→right), which shows up anywhere `-initial`/`betaK0` consumes it. ⚠ Do NOT bundle a freeze of `Ki` at the reference configuration — that was tried on the fork and REVERTED (`betaK0*getInitialStiff()` enters the residual; a frozen `K0` changes converged answers and diverges under `-initial` past ~2-8° chord rotation — see `_adr76_session_handoff.md`). The aliasing fix and the "which configuration" question are separable; ship only the aliasing fix |
| 0.8 | **Dense SOEs `exit(-1)` on a model with ZERO free equations** | Six SOEs (`FullGenLinSOE`, `BandGenLinSOE`, `BandSPDLinSOE`, `ProfileSPDLinSOE`, `SProfileSPDLinSOE`, `DiagonalSOE`) build their `vectX`/`vectB` wrappers only under `if (size != oldSize)`, so a model whose every DOF is fixed or `sp`-constrained (`size == oldSize == 0`) leaves them null and the first `getX()`/`getB()` takes the `FATAL … exit(-1)` branch — the **whole process** dies with no traceback. `UmfPack` handles the same model fine. Six one-line guards (`|| vectX == 0`) plus a `size > 0` guard on the two ProfileSPD variants' `iDiagLoc[size-1]` read (UB when an existing SOE is resized down to zero). No fork concepts, no behaviour change for any nonzero-equation model (the new branch is reachable only at the first `setSize` with `size == 0`). Fork PR (this one) with a 13-case zone_a gate that verifiably fails 6 pre-fix (all the zero-equation-from-the-start cases; the 3 -> 0 shrink half never reproduced the crash and is a guard on the `iDiagLoc` UB only) — the test ports cleanly. See [[LEDGER_quirks]] "A fully-prescribed model (zero free DOFs) ... `FullGenLinSOE::getX - vectX == 0`" |

| 0.9 | **SOE accessors kill the process when the SOE was never sized (`printA`/`printB`)** | The twin of 0.8, reached by the other door: 0.8 is `setSize()` running but skipping its wrapper block at zero free equations; this is `setSize()` **never running at all**, which any failed `domainChanged()` produces. **22 SOE classes, three distinct defect shapes.** (a) `exit(-1)` in `getA`/`getX`/`getB` on a null wrapper -- whole-process death, and from Python a silent one (clean `exit()` => no traceback, no exception, `faulthandler` mute, `opserr` redirected). `getA` exists in only two classes (`FullGenLinSOE`, `DiagonalSOE`); `getX`/`getB` in 17. Fix: warn + return an empty result -- `getA` returns 0, which `OPS_printA` and `LinearSOE::saveSparseA` **already** branch on, and `getX`/`getB` return a size-0 `Vector`, which `OPS_printB` already branches on. (b) two null dereferences in `SymSparseLinSOE` -- `zeroA()` (reached from `formTangent()`, i.e. BEFORE any accessor, so fixing accessors alone was not enough) and the DESTRUCTOR's check-after-use `while (1) { if (blkPtr->next == blkPtr) { if (blkPtr != NULL)`. (c) **five parallel `getB()` overrides with no null check at all** (`MumpsParallelSOE` + four `Distributed*`), invisible to a `grep exit(` survey because they carry no FATAL text -- they just crash (measured `0xC0000005` under `mpiexec -n 2`). Their guard must be rank-uniform or the collective hangs; it is, because `setSize()` is only reachable through the collective `domainChanged()`. Also folded in: `MPIDiagonalSOE::getpartofA`'s wrong-method FATAL text, and `DiagonalSOE`'s sized ctor testing `size > 0` before `size = N` (so the overload allocated nothing -- dead code, no callers). No fork concepts anywhere; provably inert for any SOE that was sized. Gate: `tests/test_printa_unsized_soe.py` (zone_a, subprocess-isolated, parametrized over every serial `system` x printA/printB, plus `mpiexec -n 2` rows). **Port 0.8 and 0.9 together** -- same family, adjacent lines, and 0.9's gate supersedes 0.8's. See [[LEDGER_quirks]] "`printA` / `printB` KILL the interpreter" |
### Wave 1 — small additive features, zero/near-zero vanilla footprint
| # | Package | Content | Deps |
|---|---|---|---|
| 0.14 | LAPACK `return -info+1` reports a zero first pivot as success (BandGen/FullGen/BandSPD); `algorithm` keeps the old one on a null factory | `up/10-lapack-singular-and-algorithm-null` | [jaabell#40](https://github.com/jaabell/OpenSees/pull/40) | **taken** 2026-10-13 (José applied it to `ladruño`, PR closed; opened 2026-10-06). Test 4 fail on base → 7 pass |
| 0.15 | Repeated openseespy `eigen` without analysis; `printA -sparse -ret` | `up/15-interpreter-small-fixes` | [jaabell#41](https://github.com/jaabell/OpenSees/pull/41) | **taken** 2026-10-13 (José applied it to `ladruño`, PR closed; opened 2026-10-06). `OPS_GetStringFromAll` Tcl fix DROPPED: no vanilla caller reaches the buffer path, not reproducible on his base |
| 0.16 | TenNodeTet stray shape-function print; `update()` ignores `setTrialStrain` failure | `up/16-tet10-print-and-update` | [jaabell#42](https://github.com/jaabell/OpenSees/pull/42) | **taken** 2026-10-13 (José applied it to `ladruño`, PR closed; opened 2026-10-06). Broker cases for Tet10/Brick20 DROPPED: their `sendSelf`/`recvSelf` are broken (Tet10 sends 4 of 10 node tags; Brick20 blank ctor leaves `materialPointers` null), so the case would turn a clean restore error into a crash. Fix send/recv first |
| 0.17 | Arpack: unconverged modes returned as computed (0.0) with success; ArpackSolver work arrays not resized when n grows (access violation); `getNCV` min(2nev, nev+8) skips copies of repeated eigenvalues | `up/17-arpack-unconverged-and-resize` | [jaabell#43](https://github.com/jaabell/OpenSees/pull/43) | **taken** 2026-10-13 (José applied it to `ladruño`, PR closed; opened 2026-10-06), 2 commits. Tests: base 6 fail; commit 1 → 4 pass / 2 fail; both → 6 pass. `getNCV` now nev+8 (no memory growth for nev ≥ 8). **Still live on the fork** — port back to `ladruno` as its own WP. Also seen: `-symmBandLapack` returns empty on the test chain (base and fixed), not investigated |
| 0.18 | `system Pardiso` (MKL): harden the never-compiled upstream Gen pair, build + register (Python + Tcl), `-matrixType`, `-krylov`, `-stats`, `-deterministic`/`-cbwr`; Sym pair left unbuilt | `up/18-pardiso-solver` | [jaabell#44](https://github.com/jaabell/OpenSees/pull/44) | **taken** 2026-10-29 (`feat/pardiso-solver`; opened 2026-10-06), 4 commits. Tag `SOLVER_TAGS_PARDISOGenLinSolver` = 34. CMake: Windows on with MKL; `MKL_PARDISO_LINUX` / `_THREADED` default OFF. `tests/test_pardiso.py` 15 pass (skip on base); Tcl agrees with UmfPack to 1e-15 |
| 0.19 | UW DruckerPrager two-surface return map (cutoff never assembled; tangent divides by returned norm) + DruckerPragerPlaneStrain `getInitialTangent` returns mCep | `up/19-drucker-prager-return-map` | [jaabell#45](https://github.com/jaabell/OpenSees/pull/45) | **taken** 2026-10-13 (José applied it to `ladruño`, PR closed; opened 2026-10-06), 2 commits. Test 6 fail on base (hydrostatic case hangs) → 6 pass; `tests/` 155 passed. Evidence: TenNodeTet Prandtl ratio 0.391 → 1.170. Follow-up: `DruckerPragerThermal.cpp` has the same unreachable `Jact(i) == 2` arm |
| 0.20 | `LoadControl -tangentPredictor` (stock integrator option, neutral name; `-extrapolate` dropped: failed its gate) + FE_Element/DOF_Group hooks | `up/20-loadcontrol-tangent-predictor` | [jaabell#46](https://github.com/jaabell/OpenSees/pull/46) | **taken** 2026-10-13 (José applied it to `ladruño`, PR closed; opened 2026-10-06), 3 commits. Test 6 fail on base → 7 pass. Acceptance model (vanilla J2Plasticity): cutbacks 23 → 0, iterations 224 → 12 |
| 0.21 | RBE2/RBE3 coupling elements as `KinematicCoupling` / `DistributingCoupling` (neutral names, upstream tags) | `up/21-kinematic-distributing-coupling` | [jaabell#47](https://github.com/jaabell/OpenSees/pull/47) | **taken** 2026-10-13 (José applied it to `ladruño`, PR closed; opened 2026-10-06), 3 commits. Tags 275/276; 32 tests pass (skip on base). Dropped: bipenalty/dt_cr/mass-penalty, per-iteration AL; tie force renamed `couplingForce`. Classic Tcl dispatch placed ahead of the element chain (MSVC C1061 nesting limit) |
| 0.22 | Parallel numberer: stock `ParallelNumberer` improved in place where byte-identical (Lagrange ref<0 fusion), neutral new class only if the ordering differs; DOF_Numberer MP index (#598) | `up/22-parallel-numberer` | [jaabell#48](https://github.com/jaabell/OpenSees/pull/48) | **taken** 2026-10-13 (José applied it to `ladruño`, PR closed; opened 2026-10-06), 7 commits; stock `ParallelNumberer` improved in place (numbering byte-identical, no new class/tag). Lagrange fusion: base 2.75e-3 vs exact 9.09e-4 → fixed. Full 4-target build, `tests/` 201 passed |
| 0.23 | ElasticBeam2d, ElasticTimoshenkoBeam2d/3d subtract the ground-motion inertia load twice | `up/23-elastic-beam-ground-motion-load` | [jaabell#49](https://github.com/jaabell/OpenSees/pull/49) | **taken** 2026-10-29 (José cherry-picked it onto a `fix/` branch, merged to `ladruño`, PR closed; opened 2026-10-15). 2 commits. Base 12 fail / 4 pass (ElasticBeam3d control) → 16 pass; `tests/` 253 |
| 0.24 | Constraint FEs (Penalty/Lagrange/Auto) contribute mass in `eigen` / `modalDamping` | `up/24-constraint-fe-eigen-mass` | [jaabell#57](https://github.com/jaabell/OpenSees/pull/57) | **taken** 2026-10-29 (José cherry-picked it onto a `fix/` branch, merged to `ladruño`, PR closed; opened 2026-10-15). 2 commits. Base 10 fail → 10 pass; `tests/` 247. Test model enlarged to 12 oscillators so every system exceeds ARPACK ncv = nev+8 (n = 5 aborted with info -9999 once M was correct) |
| 0.25 | AutoConstraintHandler `applyLoad` misses `updateElement`; `sp -subtractInit` no-op in Python; PathSeries drops `-useLast` | `up/25-sp-load-and-path-series` | [jaabell#50](https://github.com/jaabell/OpenSees/pull/50) | **taken** 2026-10-29 (José cherry-picked it onto a `fix/` branch, merged to `ladruño`, PR closed; opened 2026-10-15). 3 commits. Base 6 fail / 5 pass → 11 pass; `tests/` 248 |
| 0.26 | Refused `system`/`numberer`/`constraints`/`test`/`integrator` pass silently (Python) or crash (Tcl) | `up/26-refused-analysis-objects` | [jaabell#56](https://github.com/jaabell/OpenSees/pull/56) | **taken** 2026-10-29 (José cherry-picked it onto a `fix/` branch, merged to `ladruño`, PR closed; opened 2026-10-15). 2 commits (Python; Tcl droppable). Base 17 fail (Tcl: 0xC0000409, 0xC0000005) → 17 pass; `tests/` 254 |
| 0.27 | ZeroLengthND symmetrizes the tangent; ContactMaterial2D/3D initial tangent aliases the current one | `up/27-zerolength-nd-and-contact-tangents` | [jaabell#51](https://github.com/jaabell/OpenSees/pull/51) | **taken** 2026-10-29 (José cherry-picked it onto a `fix/` branch, merged to `ladruño`, PR closed; opened 2026-10-15). 2 commits; the ContactMaterial commit has no test (no command reaches `getInitialTangent`). Base 1 fail / 1 pass → 2 pass; `tests/` 239 |
| 0.28 | ASDAbsorbingBoundary3D integrates with signed det J; degenerate bricks accepted | `up/28-asd-absorbing-detj` | [jaabell#52](https://github.com/jaabell/OpenSees/pull/52) | **taken** 2026-10-29 (José cherry-picked it onto a `fix/` branch, merged to `ladruño`, PR closed; opened 2026-10-15). 2 commits. Base 9 fail / 1 pass → 10 pass; `tests/` 247 |
| 0.29 | dt = 0 accepted by TRBDF2/TRBDF3/Houbolt/BackwardEuler (NaN state); no Hall tangent in BackwardEuler/Newmark1/Collocation | `up/29-transient-integrator-guards` | [jaabell#53](https://github.com/jaabell/OpenSees/pull/53) | **taken** 2026-10-29 (José cherry-picked it onto a `fix/` branch, merged to `ladruño`, PR closed; opened 2026-10-15). 2 commits. Base 7 fail → 7 pass; `tests/` 244 |
| 0.30 | Element response 444444 depends on evaluation order; Rayleigh self-heal calls the virtual setter | `up/30-element-damping-and-inertia-response` | [jaabell#54](https://github.com/jaabell/OpenSees/pull/54) | **taken** 2026-10-29 (José cherry-picked it onto a `fix/` branch, merged to `ladruño`, PR closed; opened 2026-10-15). 2 commits. Test passes on the MSVC base (order is compiler-dependent: minimal model gives -1.2 vs 3.0 under GCC 15, 3.0 under MSVC); GCC base run not made here. Self-heal commit defensive, no reachable path, no test |
| 0.31 | Tcl `OPS_GetStringFromAll` never fills the buffer; Python `Py_AtExit` registered per import (crash at exit) | `up/31-interpreter-string-and-reimport` | [jaabell#55](https://github.com/jaabell/OpenSees/pull/55) | **taken** 2026-10-29 (José cherry-picked it onto a `fix/` branch, merged to `ladruño`, PR closed; opened 2026-10-15). 2 commits. Base 2 fail / 1 pass → 3 pass; `tests/` 240. Reverses 0.15's drop: KinematicCoupling (0.21) now reaches the buffer path |
| 0.32 | LoadPath / ArcLength / ArcLength1 / HSConstraint discard `updateDomain()` failures | `up/32-loadpath-arclength-update-domain` | [jaabell#58](https://github.com/jaabell/OpenSees/pull/58) | **taken** 2026-10-29 (José cherry-picked it onto a `fix/` branch, merged to `ladruño`, PR closed; opened 2026-10-15). 2 commits. Base 3 fail / 3 pass (elastic controls) → 6 pass; `tests/` 243. LoadPath has no interpreter command (review only) |
| 0.33 | PressureDependMultiYield / 02 / 03 uncapped substep count (huge increment hangs) | `up/33-pdmy-substep-cap` | [jaabell#60](https://github.com/jaabell/OpenSees/pull/60) | **taken** 2026-10-29 (José cherry-picked it onto a `fix/` branch, merged to `ladruño`, PR closed; opened 2026-10-15). 3 commits; cap 1e5 substeps, setTrialStrain returns -1. Base: 3 cases still running at 20 s, 3 small-increment controls pass → 6 pass; `tests/` 243 |
| 0.34 | ProfileSPD / BandSPD / SProfileSPD silently drop the unsymmetric part of element matrices | `up/34-spd-asymmetry-warning` | [jaabell#59](https://github.com/jaabell/OpenSees/pull/59) | **taken** 2026-10-29 (José cherry-picked it onto a `fix/` branch, merged to `ladruño`, PR closed; opened 2026-10-15). 3 commits; warning only, results unchanged. Base 3 fail / 4 pass → 7 pass; `tests/` 244 |
| 1.1 | Plane-strain σ_zz | `NDMaterial::getStressZZ` + `stressesPlaneStrain` responses | — |
| 1.2 | Beam localAxes | response id 30 on the 10 beam classes | — |
| 1.3 | DDM integrators | LadrunoHHT + LadrunoGeneralizedAlpha (header promotions only) | — |
| 1.4 | Robust statics | LadrunoArcLength (Ramm + STABILIZE), DynamicRelaxation, IndirectControl, StabilizedUnbalance test | — |
| 1.5 | ASDPlastic geomaterials | Hoek–Brown (rock) + StiffSoil shear/cap components for jaabell's ASDPlasticMaterial3D framework (10 new YF/PF/EL/hardening headers + regenerated registries, +2.9k lines, `test_HoekBrown.cpp`). His framework — coordinate directly; found unledgered by the audit. `HoekBrown_YF.h` is now realigned with `jaabell/ASDP` `60d9b9b23` (composite tension yield surface ported wp/94d, ADR-94 H10a/M4); `HoekBrown_PF.h`/`_Utils.h`/`_ParameterTypes.h` and the StiffSoil family remain unreconciled with his tree | — |

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
> **19 modified vanilla files with no marker** (TenNodeTet 6× fix — *marker
> since added, see 0.0*; DistributedSuperLU MSVC fix, the ASDPlastic extension
> set, profiler-macro hooks, `TransientIntegrator.h getCriticalTimeStep`
> virtual, build wiring).
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
| 0.0 | TenNodeTet 6× stiffness fix | `up/00-tennodetet-shp3d` | [jaabell#29](https://github.com/jaabell/OpenSees/pull/29) | **superseded** 2026-10-05: José fixed it his own way on `fix/tet10-integration` (removes `/6.0` from Jdet); we are credited as co-authors upstream. Merging ours on top would have made the element 6× too stiff |
| 0.0b | TenNodeTet `getResponse` heap overrun | — | — | **done by José** (`110bb4dcc`, with `setParameter` all-GP and `-doInitDisp` fixes) |
| 0.1 | Portability & crash (FE_Element, PythonStream, SuperLU MSVC) | `up/01-portability-crash-fixes` | [jaabell#30](https://github.com/jaabell/OpenSees/pull/30) | **taken** 2026-10-05 (FE_Element + PythonStream cherry-picked onto `fix/fe-element-pythonstream`; SuperLU `stat` landed upstream as `73200b004`) |
| 0.2 | quad/tri rho serialization | `up/02-quad-tri-rho-serialization` | [jaabell#31](https://github.com/jaabell/OpenSees/pull/31) | **taken** 2026-10-05 (`fix/quad-tri-rho-serialization`) |
| 0.2b | GmshRecorder hex20 | `up/03-gmsh-hex20` | [jaabell#32](https://github.com/jaabell/OpenSees/pull/32) | **taken** 2026-10-05 (`fix/gmsh-hex20`) |
| 0.3a | Domain::clearAll EQ_Constraint leak | `up/04-domain-clearall-eq-leak` | [jaabell#33](https://github.com/jaabell/OpenSees/pull/33) | **taken** 2026-10-05 (`fix/domain-clearall-eq`) |
| 0.3b | Error-return honoring (DirectIntegration/TransientDD) | — | — | HELD for José (policy change) |
| 0.3c | Mumps `-opt` parse guard | — | — | HELD (low value) |
| 0.4 | Registration gaps (Lysmer, InitStrain dim-general, ASDP setResponse) | — | — | HELD for José |
| 0.5 | H5DRM | — | — | **partly done by José** (cfactor, hold-final, tend dataspace, 6-DOF skip: `feat/h5drm-cfactor-hold-final`). Still live on `ladruño`: `stuff[12]` uninitialized in `TclPatternCommand.cpp:539` + runtime parser 3-arg ctor |
| 0.6 | Byte-identical perf fixes | — | — | HELD (profiler strip + suite pass first) |
| 0.7 | CorotCrdTransf3d static `T` | — | — | **not fixed on the fork either** (`.h:135` still `static` on both); nothing to port until fixed here |
| 0.8 + 0.9 (+ #759) | SOE accessors: zero-equation wrappers, unsized `exit(-1)` ×17 classes, SymSparse null derefs, 5 parallel `getB`, FullGen `Bsize` sizing | `up/05-soe-unsized-and-zero-equation` | [jaabell#35](https://github.com/jaabell/OpenSees/pull/35) | **taken** 2026-10-13 (José applied it to `ladruño`, PR closed; opened 2026-10-06). New test 16 fail / 15 pass on base → 31 pass |
| 0.10 | GeneralizedAlpha `update()` discards `Ualphadotdot` | `up/06-generalizedalpha-alpham-inertia` | [jaabell#36](https://github.com/jaabell/OpenSees/pull/36) | **taken** 2026-10-13 (José applied it to `ladruño`, PR closed; opened 2026-10-06). SDOF order test: 3 of 4 fail on base → 4 pass. Results change for every alphaM ≠ 1 |
| 0.11 | Analysis-object pools skip slot `[MAX_NUM_DOF]` (11 loops, 4 files) | `up/07-analysis-pool-max-dof-slot` | [jaabell#37](https://github.com/jaabell/OpenSees/pull/37) | **taken** 2026-10-13 (José applied it to `ladruño`, PR closed; opened 2026-10-06). No portable test (heap-state dependent); body carries the fork's crash table |
| 0.12 | Newmark file-scope `static bool converged` | `up/08-newmark-static-state` | [jaabell#38](https://github.com/jaabell/OpenSees/pull/38) | **taken** 2026-10-13 (José applied it to `ladruño`, PR closed; opened 2026-10-06) |
| 0.13 | Windows MSVC + ifx + MUMPS build fixes (8 commits: `/bigobj`, version define and globbed includes C/C++-only, MPI 8.3 paths, MUMPS `.lib`, LP64 ScaLAPACK, `MUMPS_INCLUDE_DIR`, per-exe Tcl domain sources) | `up/09-windows-msvc-ifx-build-fixes` | [jaabell#39](https://github.com/jaabell/OpenSees/pull/39) | **taken** 2026-10-30 (`fix/windows-msvc-build`; opened 2026-10-06). Windows re-check of `efd1be55a` posted 2026-10-15 at José's request: base fails all four targets (every Fortran source crashes ifx 0xC0000005; SP/MP miss `dmumps`), with the PR all four build, `tests/` 219 passed / 2 skipped; `/bigobj` no longer needed on that base. OpenSees/SP/MP/Py build on his base; `tests/` 112 passed; SP and MP on 2 ranks with Mumps match serial UmfPack to 1e-15. Proven needed: his base does not compile here without the first four. Old patches 1 and 4 already upstream; OpenSeesPyMP target excluded (feature) |
| 0.35 | UmfPack `setSize` appends the new pattern after a graph change | `up/35-umfpack-setsize-pattern` | [jaabell#62](https://github.com/jaabell/OpenSees/pull/62) | **PR open** 2026-10-30. 2 commits. Base 3 fail → 3 pass |
| 0.36 | ASDConcrete1D/3D hardening-law registry survives `wipe` | `up/36-asdconcrete-wipe-storage` | [jaabell#63](https://github.com/jaabell/OpenSees/pull/63) | **PR open** 2026-10-30. 2 commits. Base 4 fail → 4 pass |
| 0.37 | GradientInelasticBeamColumn section locations with L = 0; DispBeamColumnNL3d 12-DOF localForce | `up/37-gradient-inelastic-and-dispnl3d-response` | [jaabell#64](https://github.com/jaabell/OpenSees/pull/64) | **PR open** 2026-10-30. 2 commits. Base 2 fail → 2 pass. The NL3d commit has no test: the parser requires ndm 2, so the 3D element cannot be created. `responseID == 9` claim dropped, reachable on his base |
| 0.38 | `responseSpectrumAnalysis -scale` ignored | `up/38-rsa-scale` | [jaabell#65](https://github.com/jaabell/OpenSees/pull/65) | **PR open** 2026-10-30. 1 commit. Base 5 fail / 1 pass → 6 pass |
| 0.39 | `modalProperties` drops the mass of MP-constrained DOFs under Transformation | `up/39-modalprops-mp-retained-mass` | [jaabell#66](https://github.com/jaabell/OpenSees/pull/66) | **PR open** 2026-10-30. 1 commit. Base 1 fail / 1 pass → 2 pass |
| 0.40 | TransformationFE `dofData` overrun above 64 nodes | `up/40-transformationfe-dofdata-size` | [jaabell#67](https://github.com/jaabell/OpenSees/pull/67) | **PR open** 2026-10-30. 1 commit. Base 0xC0000374 → pass; allocator dependent |
| 0.41 | ElasticTimoshenkoBeam3d ignores its geomTransf | `up/41-timoshenko3d-geomtransf` | [jaabell#68](https://github.com/jaabell/OpenSees/pull/68) | **PR open** 2026-10-30. 1 commit, results-changing. Base 10 fail → 10 pass. Calls `update()` before forming Ki, which the fork does not (nmorabowen/OpenSees issue filed) |
| 0.42 | ASDShellQ4/T3 corotational tangent has no geometric stiffness (`-pDelta` not sent) | `up/42-asdshell-corotational-tangent` | [jaabell#69](https://github.com/jaabell/OpenSees/pull/69) | **PR open** 2026-10-30. 2 commits. Base 6 fail / 2 pass → 8 pass; Newton 33 → 4, period ratio 1.000 → 1.405 |
| 0.43 | `system Pardiso` does not compile against oneMKL 2026.1 | `up/43-pardiso-onemkl-2026` | [jaabell#61](https://github.com/jaabell/OpenSees/pull/61) | **PR open** 2026-10-30. 1 commit. Needed to build his base here |
| 1.1 | Plane-strain σ_zz | `up/44-plane-strain-stress-zz` | [jaabell#70](https://github.com/jaabell/OpenSees/pull/70) | **PR open** 2026-10-30. 2 commits. Base 7 fail → 7 pass |
| 1.2 | Beam localAxes responses | `up/45-beam-local-axes` | [jaabell#71](https://github.com/jaabell/OpenSees/pull/71) | **PR open** 2026-10-30. 1 commit. Base 30 fail → 30 pass |
| 1.3 | DDM HHT/GeneralizedAlpha | `up/46-ladruno-ddm-integrators` | [jaabell#72](https://github.com/jaabell/OpenSees/pull/72) | **PR open** 2026-10-30. 2 commits, tags 33013/33014. Base 12 fail → 12 pass. Re-forms the current tangent before the sensitivity solves; the fork does not (ModifiedNewton/`-initial` gradients off by about 1e-3) |
| 1.4 | Robust statics (ArcLength/DR/IndirectControl) | `up/47-ladruno-robust-statics` | [jaabell#73](https://github.com/jaabell/OpenSees/pull/73) | **PR open** 2026-10-30. 4 commits, tags 33000, 33004–33006. Base 75 fail / 5 pass / 3 skip → 83 pass. classTags.h one-hunk conflict with #72 (keep both) |
| 1.5 | ASDPlastic Hoek–Brown + StiffSoil | — | — | **done by José** (`feat/asdp-hoekbrown`, `feat/asdp-stiffsoil`, 2026-10-01). Our later ASDP review fixes are a separate package (see §8) |
| 2.1 | Explicit dynamics I | — | — | not started. ⚠ José merged his own `ExplicitBathe -lnvd` and `ExplicitDifferenceStatic` rework (2026-10-02): reconcile before porting |
| 2.2 | Mass scaling | — | — | not started |
| 2.3 | Projection handler | — | — | not started |
| 2.4 | Finite-strain material infra | `up/48-finite-strain-infra` | [jaabell#74](https://github.com/jaabell/OpenSees/pull/74) | **PR open** 2026-10-30. 5 commits, nD tags 33010/33013/33014/33016. Base 25 fail / 4 pass → 29 pass; no stock element drives the finite-strain path, so command-level and kernel tests only |
| 2.5 | Energy + recorders | — | proposal emailed 2026-10-06 | **proposed** (LadrunoRecorder). ⚠ José merged his own `EnergyBalanceRecorder` (tag 26, 482 lines, no shared kernel) on 2026-10-01; the email proposes reconciling both into `EnergyBalanceKernel.h`. Pending his answer on kernel, command name (`ladruno` vs neutral) and one or two PRs |
| 3.1–3.6 | Element / material catalog | — | — | not started |
| 3.7 | Modal/eigen (complexEigen, modalResponseHistory, frequency/random response, RSA `-combine`, per-mode RSA) | `up/49-ladruno-modal` | [jaabell#75](https://github.com/jaabell/OpenSees/pull/75) | **PR open** 2026-10-30. 7 commits. Base 76 fail / 33 pass → 109 pass. `wipe` clears eigenvalues; per-mode output lacks the Ritz keys |
| W4 | contact+tie / FEAST / porous / rigid body / OpenSeesPyMP / profiler / numberer | — | — | deferred (decide w/ José) |

## 7. Decision & session log (append-only, newest first)

- **2026-10-30 — round 5 sent: jaabell#70–#75 (packages 1.1, 1.2, 1.3, 1.4, 2.4, 3.7).** José asked for faster sends. Method: one porting agent per package in its own worktree on `5deb418e5`, a lint over added lines, commit messages and PR bodies (`Ladruno_scripts/upstream_lint.py`: internal references, trailers, fork paths, em dashes, we/you), one stacked build of all branches per round (OpenSees, OpenSeesSP, OpenSeesMP, OpenSeesPy) plus a per-package run of the new tests on the unmodified base, and a web check of every citation. Stacked `tests/`: 716 passed, 6 skipped; base 439 passed, 13 skipped. **Owner decision: classes keep their Ladruno names upstream** (text rules still apply). Port fixes not yet in the fork: LadrunoHHT/GeneralizedAlpha tangent re-form before DDM; optional-number parsing in LadrunoArcLength/DynamicRelaxation/IndirectControl (a failed read no longer eats the next flag); DR non-finite guard; InitDefGrad refuses singular F0; LogStrain reports inner `setTrialStrain` failure, and revert resets the staged trial.

- **2026-10-30 — round 4 sent: jaabell#61–#69.** José took every earlier PR by 2026-10-30 (#39, #44, #49–#60), none left open. Fork fixes dated after the 2026-10-06 inventory were re-audited on `5deb418e5`; nine packages (0.35–0.43) were ported and verified together: new tests fail on the base, pass with the commits, `tests/` 483 passed, 5 skipped. Not built branch by branch. Still held: ManzariDafalias stage-static fix (#1005, owner HOLD on the family), H5DRM vertical flip and chunk cache (WP-186/187, coordinate first), ASDShell `-pDelta` (feature). José opens his own upstream OpenSees PRs from these (#1856 merged by Michael Scott on 2026-10-29).

- **2026-10-15 — batch D sent: jaabell#58–#60.** The three small coordinate-list fixes that need no design decision (static-integrator update failures, PDMY substep cap, SPD asymmetry warning), verified on `8ddefbbec` like round 3. A coordination email (names, class tags, base-class hooks, recorders, contact, elements, PR size, open 2026-10-06 questions) was drafted for the owner to send.

- **2026-10-15 — naming decision (owner): element technology and recorders keep the `Ladruno` name.** Wave 3 elements/materials (LadrunoBrick, LadrunoSolidShell, LadrunoQuad/CST/LST, Bezier*, LadrunoJ2, LadrunoConcrete3D, ...) and the 2.5 recorders (`recorder ladruno`, `.ladruno` HDF5 format, monitor) ship upstream under their Ladruno names. This overrides the neutral-name default used for the couplings (0.21) and the numberer (0.22) for these two areas only; contact, ties and other kinematic constraints keep the neutral-name default unless the owner says otherwise. Rule 6 still applies to comments/commit text (no ADR/WP refs); the class and command names are the exception.

- **2026-10-15 — round 3 sent: jaabell#49–#57, #39 re-checked.** José took #35–#38, #40–#43 and #45–#48 on 2026-10-13 (applied by hand, PRs closed); #39 (build fixes) and #44 (Pardiso) stay open. At his request #39 was re-verified on Windows against `efd1be55a` (base fails all four targets, PR builds all four, `tests/` 219 passed). A re-audit of the §8 backlog against `8ddefbbec` found every item still live except the Matrix::Solve alias (fork-only threading regression); nine packages (0.23–0.31) were ported on fresh `up/` branches by parallel agents and each verified centrally on his base plus the #39 build fixes: new tests fail on the base, pass with the commits, full `tests/` green. TenNodeTetrahedron: nothing left to send (José took or superseded all of it); the reverse sync into the fork is nmorabowen/OpenSees#991.

- **2026-10-06 — coordination email sent to José (ASDP + explicit dynamics).** Issues and discussions are disabled on jaabell/OpenSees, so it went by email. ASDP: A1–A5 fixes (special_return tangent, strict_convergence, DP dilatant apex, Numerical_Algorithmic on the committed map, StiffSoil NaN), A6 MC tension cutoff, A7 Closest_Point/Algorithmic (23/46; footing 7/10→10/10, 7/12→12/12, 6/10→9/10, 6/12→11/12); open findings HoekBrown_PF::g Tresca collapse (26 %) and static shared state. Explicit: build on his ExplicitBathe -lnvd (acceleration form is better than ours); offered E1 -sms, E2 -consistent (LinearSOE virtuals, MPI PCG not CI-gated), E3 CD + HRZ, E4 criticalTimeStep, E5 energy channels (waits on recorder). Six questions pending. Round paused for his answers.

- **2026-10-06 — round 2: four bug-fix PRs open, build fixes building, recorder proposed.**
  José closed #29–#33 on 2026-10-05: all taken by cherry-pick (authorship kept)
  except #29, superseded by his own Tet10 fix. He asked that PR text and code
  follow the office tone with no internal references (rule 6). Since 2026-10-01
  he also ported on his own: Tet10 heap/setParameter/initDisp, H5DRM
  cfactor/hold-final, ASDP Hoek–Brown + StiffSoil, a separate
  EnergyBalanceRecorder, ExplicitBathe `-lnvd`, ExplicitDifferenceStatic, Tcl
  `-lumped` for (MPI)Diagonal, optional MUMPS in CMake. A fresh inventory of
  fork vanilla fixes still live on `ladruño` found 32 candidates (§8).
  Opened jaabell#35–#38 (0.8+0.9, 0.10, 0.11, 0.12) and #39 (0.13); round 3 block 1 opened as #40–#42 (0.14–0.16), each failing on his base and passing with the fixes, `tests/` 128 passed, each verified on his base:
  unmodified base built with the local Windows fixes, new tests fail there,
  `tests/` 147 passed / 2 skipped with all four merged. His base needed `/bigobj`
  (his own ASDP registry overflows MSVC's section limit) and three ifx fixes to
  compile at all, which became package 0.13. The recorder (2.5) was proposed
  by email to José, signed "Ladruno Guppi Team".

- **2026-07-22 — Wave 0 clean bug-fixes shipped, remainder HELD for José.**
  Six upstream PRs now open: jaabell#29 (TenNodeTet 6×), #30 (portability trio),
  #31 (quad/tri rho), #32 (GmshRecorder hex20), #33 (Domain::clearAll EQ leak).
  These are the unambiguous, verifiable, mostly-byte-identical bug fixes.
  **Deliberately HELD** (per the user's "we'll wait for José" — and on merit):
  0.3b error-return honoring (a *behavioral policy change* — aborts setups that
  stock ran on; motivated by fork handlers that don't exist upstream);
  0.3c Mumps `-opt` (extraction risk from the FEAST/numberer-forked
  OpenSeesCommands.cpp, low value); 0.4 registration gaps (additive features,
  InitStrain dim-general is really an enhancement); 0.5 H5DRM (`_H5DRM`
  build-gated so unverifiable here, and the z-flip/hold-final are José's own
  patches); 0.6 perf harvest (needs profiler-strip on parallel code + a build
  pass first). Port mechanics note: repo has **mixed CRLF/LF** — the rho port
  script had to be newline-aware per file; github.com DNS kept dropping, so all
  push/PR calls run in a retry loop.
- **2026-07-22 — package 0.1 shipped as [jaabell#30](https://github.com/jaabell/OpenSees/pull/30)**,
  and re-scoped. Kept as 3 pure, non-behavioral, always-relevant fixes
  (FE_Element dead guard, PythonStream `%s`, DistributedSuperLU MSVC `stat`
  rename), one commit each so José can drop any individually. **Pulled OUT of
  0.1:** H5DRM `stuff[12]` init → folded into 0.5 (H5DRM is `_H5DRM` build-gated
  and its other changes are behavioral — one subsystem, one PR); GmshRecorder
  hex20 → 0.2 (needs the type-17 mid-edge permutation table, not a one-liner).
  Rationale: keep the first "bundle" PR trivially reviewable and free of
  build-gated code the reviewer may not compile.
- **2026-07-22 — miuandes.cl correction.** Co-author emails are
  `pxpalacios@miuandes.cl` + `jaabell@miuandes.cl` (NOT `@uandes.cl`).
  jaabell#29 amended again. Supersedes the prior `@uandes.cl` entry.
- **2026-07-22 — canonical co-author emails set by Nicolas** (supersedes the
  git-history-harvested ones): Patricio `pxpalacios@uandes.cl` (NOT
  ppalacios92@gmail.com), José `jaabell@uandes.cl`, Nicolas
  `nmorabowen@gmail.com`. jaabell#29's commit amended + force-pushed accordingly.

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

## 8. Round 3 backlog (inventory of 2026-10-06, all live on `ladruño`)

> [!note] 2026-10-15: the class-A groups `up/11`–`up/14`, `up/17` (minus DruckerPrager, sent as 0.19) and the `up/15` leftovers were sent as 0.23–0.31 (jaabell#49–#57); `up/18` stays on the owner HOLD. The coordinate-with-José list below is unchanged and still live on `8ddefbbec`.

> [!warning] HOLD (owner, 2026-10-06): SANISAND and its ManzariDafalias-family
> seams are NOT to be sent upstream while the SANISAND work (re-seat, tension
> cutoff, regularization, footing campaign) is in progress. This covers
> LadrunoSANISAND and, until the owner says otherwise, the ManzariDafalias
> packages below (`up/18-manzari-fspm-platerebar` crash fixes and the
> results-changing set), since they touch the same family.


Pure fixes (class A), proposed grouping:
- `up/10-lapack-singular-and-algorithm-null`: LAPACK `return -info+1` reports a singular matrix as success (BandGen/FullGen/BandSPD); `OPS_Algorithm` returns 0 on a null factory (#642).
- `up/11-elastic-beam-ground-motion-double-inertia`: ElasticBeam2d / ElasticTimoshenkoBeam2d/3d subtract the ground-motion load twice (#854). Results change.
- `up/12-element-scratch-and-eval-order`: Element.cpp `setRayleighDampingFactors` self-heal qualifiers (#676; upstream reachability unclear) + response 444444 evaluation order (#859).
- `up/13-transient-integrator-guards`: dt = 0 guards in TRBDF2/TRBDF3/Houbolt/BackwardEuler; `HALL_TANGENT` branch in BackwardEuler/Newmark1/Collocation (#650). New tests needed.
- `up/14-sp-constraint-and-path-series`: AutoConstraintHandler `applyLoad` misses `updateElement` (#697, results change); `OPS_SP -subtractInit` inverted (#675); PathSeries `-useLast` dropped (both routes).
- `up/15-interpreter-small-fixes`: repeated `eigen` (#609), `printA -sparse -ret` (#761), `OPS_GetStringFromAll` under Tcl (#840); Python re-import `Py_AtExit` (#712) separately.
- `up/16-tet10-and-quadratic-broker`: Tet10 stray `std::cout`, `update()` discards `setTrialStrain`, broker cases for Tet10 and 20-node brick. Rebase on José's Tet10 code.
- `up/17-nd-initial-tangent-and-zerolength-symmetrize`: ZeroLengthND lower-triangle mirror, DruckerPragerPlaneStrain `getInitialTangent` returns mCep, ContactMaterial2D/3D tangent aliasing (#720). Results change.
- `up/18-manzari-fspm-platerebar`: ManzariDafalias parser overrun, ForwardEuler shadowed `r`, uninitialized `nG,nK` in MaxStrainInc/MaxEnergyInc (#901, #914); FluidSolidPorous `getCopy` on the UW family; PlateRebar `recvSelf` angle.

Coordinate with José first (class B/C): H5DRM `stuff[12]` (0.5 remainder); PDMY substep cap (#874, uses a fork return code → -1); AutoConstraintHandler MPI `KAVG` (#733); ManzariDafalias results-changing set (void-ratio interpolant, ME clamp, `ToCovariant` 2×, `mElastFlag`); LoadPath/ArcLength `updateDomain()` return (#792); DruckerPrager two-surface return-map repair (#803); ProfileSPD lower-triangle discard; ASDPlasticMaterial3D review fixes (his framework).
