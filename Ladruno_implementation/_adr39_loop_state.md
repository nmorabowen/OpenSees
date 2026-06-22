# ADR 39 ContactDomain — implementation loop state

> Persistent driver state for the phased ContactDomain build. Survives context
> resets. One source of truth for "where are we, what's next." Mirrors the
> ADR-30 `_adr30_loop_state.md` pattern. ADR = `39_ladruno_contact_domain_adr.md`.

## ⇒ HANDOFF POINT (2026-06-21): P1 shipped in PR #345; resume at P2a.
See `contact_p2_handoff.md`. P1 = ContactDomain skeleton + handler + lifecycle +
parser (zero-force), 7/7 green. P2 designed + gated (split P2a/P2b). NEXT = code P2a
(rigid plane, -kn val) — study PenaltySP_FE/TransformationFE for the custom-FE
connectivity. Verify #345 merged before stacking.

## Scheme

Phased loop, explicit-first ship at P3. Each phase:
**design → (adversarial design gate at critical junctions) → implement → build →
test → adversarial code gate → commit/PR → advance.**
Critical junctions requiring an adversarial Workflow gate: **P1** (architecture
wiring — handler injection + conservative-static connectivity, the BLOCKER-risk),
**P2** (narrow-phase mechanics — projection/penalty/tangent), **P3** (friction
return map). P0 = validation only (light verification, no heavy gate).

## Environment (verified 2026-06-21)

- Run python: `C:/Users/nmora/AppData/Local/Python/pythoncore-3.12-64/python.exe`
  (numpy 2.4.6, pytest). Bootstrap: `os.add_dll_directory(DIST); sys.path.insert(0,DIST)`
  where `DIST = C:/Users/nmora/Github/OpenSees_Compile/OpenSees/dist/bin`.
- Built binary = commit 8c8edb21 (may lag the rebased worktree HEAD; fine for
  no-SRC P0). `OpenSees.exe` + `opensees.pyd` both present. `ZeroLengthContactASDimplex`
  available.
- Build cmd (P1+): `Ladruno_scripts\build.bat OpenSees OpenSeesSP OpenSeesMP`
  (Windows, heavy ~30min, can fail). Use the conan cmake full-path if ZLIB-not-found
  (see [[project_build_env_cmake43_conan_zlib]]).
- Tests: `tests/`, `from _testbed import ops`, pytest markers `zone_a`/`zone_b`.

## Phase ledger

| Phase | State | Gate | PR | Notes |
|---|---|---|---|---|
| P0 falsify/baseline (no SRC) | **DONE ✓** | light verify ✓ | local | both protos pass; 2 design rules extracted |
| P1a FE+handler+empty-conn zero, bitwise | **DONE ✓** (5eb3b810) | design gate ✓ + code gate ✓ | local | rebuilt + 3/3 bitwise green w/ fixes; committed |
| P1b ContactDomain+surface+lifecycle hooks | **DONE ✓** (344a9c86) | light review ✓ (test-covered + clean compile; Domain non-copyable → no double-free) | local | 7/7 green; build exit 0; committed |
| P2 NTS penalty narrow phase (frictionless) | DESIGN GATE ✓ → SPLIT P2a/P2b | full gate ✓ (SALVAGEABLE; folded) | — | design revised; P2a next (rigid plane, -kn val) |
| P2a rigid+inclined plane, -kn val | **GREEN ✓** (3/3 + P1 regression 10/10) | code gate PENDING (fold into P2b gate) | guppi/contact-p2a | mechanics exact; committing |
| P2b faceted+deformable+Hertz+auto-kn+∂n/∂u | NOT STARTED | full gate already covers; code gate | — | the deformable mechanics |
| P2 NTS penalty frictionless | NOT STARTED | **adversarial gate** | — | rigid-plane rung first |
| P2.5 bucket sort drop-in | NOT STARTED | verify==brute force | — | — |
| P3 IMPL-EX Coulomb — SHIP | NOT STARTED | **adversarial gate** | — | v1 ship |
| P3.5 implicit frictional Newton | DEFERRED | — | — | post-v1 |

## Loop log

### Iteration 1 — 2026-06-21 — P0 start
- Env recon complete; runtime verified (opensees.pyd + numpy + ZLC-ASDimplex).
- Loop infra created. Starting P0 prototypes under `contact_prototypes/`.

### Iteration 1 (cont.) — P0 DONE ✓
- `proto_bucket_sort.py` — ALL PASS. Caught the LS-DYNA §26.8.1 deficiency LIVE:
  centroid-only bucketing (type-4) missed 610/2000 (30%) true-near pairs. Fixed
  with **segment-span registration (type-13, eq 26.44)** → 0 misses. Runaway-node
  bbox-clip guard verified (real surface collapses to 1 bucket without it,
  recovers to 100 with). 33× fewer checks than brute force.
- `proto_baseline_asdimplex.py` — ALL PASS. Penalty penetration δ=P/Kn EXACT
  (0.00% err). Explicit impact STABLE, restitution **e=1.000** (elastic,
  energy-conserving), max-penetration bounded = 2·gap. Confirms force-injection /
  robust-by-construction thesis empirically.
- **C++ design rules extracted (→ feed P1/P2.5):**
  1. Bucket sort = segment-span registration (corner bucket-pointer bbox, eq 26.44),
     NOT centroid-only. Centroid-only is the documented-deficient type-4.
  2. `LMAX = max(median_seg_diagonal, r_search) · safety_margin` — cell ≥ search
     radius for the ±1-neighbour superset guarantee. r_search = gap_tol + max slip/epoch.
  3. Runaway-node bbox percentile-clip guard (one fly-off node else collapses buckets).
  4. Penalty penetration = P/Kn is exact; explicit penalty contact is energy-
     conserving (e=1.0) — Q4-a force injection validated.
- Gate decision: P0 is no-SRC validation → light self-review only; the heavy
  adversarial Workflow gate is reserved for P1 (wiring) / P2 (mechanics) / P3 (friction).

### Iteration 2 — P1 design + adversarial design gate
- Read real interfaces: ConstraintHandler (handle() is the only FE factory, called
  in domainChanged after clearAll), FE_Element bare ctor (tag,numDOF,ndof) for
  adapters, **ConstraintHandler::update() is NOT called in the transient step loop**
  (loop = newStep→solveCurrentStep→commit) → per-step contact trigger can't be there.
- KEY DESIGN DECISION: **self-contained adapter** — LadrunoContactFE computes its
  own narrow phase inside getResidual()/getTangent() (reads Node::getTrialDisp like
  a normal Element) → NO per-step global trigger, NO vanilla integrator edit. Broad
  phase only at epoch re-emit.
- Wrote `_adr39_p1_design.md` with 5 open questions (Q-P1-1 trigger/trial-disp,
  Q-P1-2 handler delegation, Q-P1-3 zero-force-NOT-byte-identical-graph, Q-P1-4
  pair-state-commit, Q-P1-5 conservative connectivity).
- **Launched adversarial DESIGN gate Workflow `wv35sge9t`** (4 source-grounded
  reviewers: trigger · handler/numbering · state/lifecycle · completeness → synth).
  Running in background. NEXT on completion: fold the verdict + fix list into the
  P1 design, THEN write C++.

### Iteration 3 — P1 design gate verdict (SALVAGEABLE-WITH-CHANGES) FOLDED
Gate caught 3 real bugs the zero-force test would NOT expose (but corrupt state at P3):
- **BLOCKER-1:** commit-hook story refuted by all 4 — ConstraintHandler never called
  at commit; its only per-step entry applyLoad() is pre-solve on REJECTABLE trial
  state. FIX: state on Domain-owned LadrunoContactDomain, committed via `// Ladruno`
  hook in `Domain::commit()`; adapters = STATELESS VIEWS re-bound by pair-key in handle().
- **BLOCKER-2:** missing `Domain::revertToLastCommit()` hook (failed implicit steps
  revert → friction state leaks into retry). Added.
- **MAJOR-3:** implicit getTangent must return `c1·K_c` (Newmark addKtToTang(c1)), not raw K_c.
- **MAJOR-4:** "explicit never calls getTangent" imprecise — getTangent IS called,
  delegates to mass-only formEleTangent under CDL. Also myEle==0 ⇒ base helpers
  early-return/exit(-1) ⇒ adapter MUST own its Vector/Matrix buffers + override both.
- **Q-P1-2:** REPLICATE PlainHandler (not compose) — FE-tag collision + MP -4 elim.
- **Q-P1-3:** EMPTY connectivity (size-0 getID) for inactive P1 adapter → genuinely
  BITWISE identical (connectivity perturbs the numberer graph even at zero force).
- **MINOR-10:** getExplicitCriticalTimeStep DEAD on adapters (CriticalTimeStep scans
  Domain elements only) → removed from P1, route via ContactDomain at P3.
- classTags: HANDLER_TAG_LadrunoContactHandler 33002, ELE_TAG_LadrunoContactFE 33015.
- **P1 SPLIT: P1a** (FE + handler + `constraints LadrunoContact`, ≈1 vanilla edit,
  bitwise gate) → **P1b** (ContactDomain on Domain + surface + commit/revert hooks
  + parser). Folded into `_adr39_p1_design.md` + ADR table + hybrid-table precision.
- NEXT: write P1a C++ (read LadrunoProjectionHandler::handle + PlainHandler::handle
  as the replicate template), build (background), test bitwise under CDL+Newmark.

### Iteration 4 — P1a C++ WRITTEN + committed (e6943c3f), build running
Files (all in SRC/analysis/handler/ for P1a; LadrunoContactFE moves to SRC/domain/contact at P2):
- `LadrunoContactFE.{h,cpp}` — FE_Element adapter, bare ctor (0,0)=empty conn,
  owns Vector resid(0)+Matrix tang(0,0), overrides getResidual/getTangent (zero)
  + addMtoTang no-op. (myEle==0 ⇒ base helpers exit(-1), so owning buffers mandatory.)
- `LadrunoContactHandler.{h,cpp}` — ConstraintHandler, replicates LadrunoProjectionHandler
  DOF/FE loop, tracks numFe, injects ONE LadrunoContactFE(numFe++). MP present⇒warn
  (P1b delegates). OPS_LadrunoContactHandler factory + sendSelf/recvSelf→0 stubs.
- Wiring: classTags HANDLER_TAG_LadrunoContactHandler 33002; `constraints LadrunoContact`
  branch in OpenSeesCommands.cpp (mirrors LadrunoProjection, 1 vanilla edit); CMake
  (handler dir); stamp_headers globs + stamped (4 files).
- `tests/test_adr39_contact_p1.py` — bitwise gate (truss chain, CDL+Newmark, LadrunoContact
  ==Plain exactly) + rebuild smoke (remove element mid-run).
- BUILD: `build.bat OpenSeesPy` (worktree = COLD build, ~30-45min, runs conan). bg bcc7upkiq.
  build.bat roots at script location → BUILD_DIR=worktree/build/Release (no warm cache here).
- NEXT on build done: run pytest test_adr39_contact_p1.py; if green → P1a code gate
  (adversarial) → ledger row + banner? (handler, not banner-worthy yet) → commit → P1b.
- NOTE: only the OpenSeesPy target built (test uses opensees.pyd). Full OpenSees/SP/MP
  rebuild + ledger LEDGER_implementations row deferred to the P1 PR.

### Iteration 5 — BUILD ✓ (exit 0) + P1a TEST ✓ 3/3 + code gate ✓
- Cold worktree build succeeded (bcc7upkiq, exit 0); worktree now has its OWN
  dist/bin/opensees.pyd (must run tests against THIS dist, not the compile-root one).
- `test_adr39_contact_p1.py` **3/3 PASS**: bitwise-identical LadrunoContact==Plain
  under BOTH CentralDifferenceLadruno (explicit) AND Newmark (implicit) + rebuild smoke.
  The hybrid FE_Element-adapter injection is validated in compiled code.
- **P1a code gate** (1 adversarial agent, source-grounded) → NO BLOCKER; core memory/
  ownership/graph-neutrality SOUND. 4 fixes applied (rebuild pending):
  1. (MAJOR) broker: `FEM_ObjectBrokerAllClasses::getNewConstraintHandler` case
     HANDLER_TAG_LadrunoContactHandler + include (parallel/db restore).
  2. (MAJOR) Tcl: `commands.cpp specifyConstraintHandler` `LadrunoContact` branch
     (was openseespy-only → Tcl users got a hard error).
  3. (REQUIRED) ledgers: LEDGER_implementations row + 5 LEDGER_vanilla_files rows.
  4. (MINOR) LadrunoContactFE: override getM_Force/getC_Force/getTangForce/getK_Force/
     getKi_Force → size-0 (base warns on myEle==0 under modal damping / non-diag mass).
- NEXT: rebuild OpenSeesPy (incremental, warm cache → fast), re-run test_adr39
  (still green), commit P1a fixes, then START P1b (ContactDomain on Domain + surface
  + Domain::commit/revertToLastCommit hooks + contactSurface/contact parser).
  (Build/commit briefly blocked on command classifier outage — retrying.)

### Iteration 6 — build/commit blocked (persistent classifier outage); P1b design de-risked
- Could not rebuild/commit the P1a fixes (Bash + PowerShell both gated on the
  command-safety classifier, temporarily unavailable). All work below is files-only.
- **P1b design RESOLVED to build-ready** (`_adr39_p1b_design.md`):
  - Q-P1b-1 (lifetime): delete/reset `theContactDomain` in `Domain::clearAll()`
    (`Domain.cpp:1041`, the wipe path) — mirrors the ADR-30 `theEQs` leak-fix at :1054.
    domainChanged (AnalysisModel::clearAll) correctly leaves it alone.
  - Q-P1b-3 (cross-lib): LadrunoContactDomain + LadrunoContactSurface → `SRC/domain/contact/`
    (OPS_Domain); handler already links OPS_Domain (uses Domain.h) ⇒ no cycle. FE stays
    in handler dir (OPS_Analysis) for P1b.
  - commit/revert hook insertion points: `Domain.cpp:2185` (commit) + `:2211` (revert).
- HELD: did NOT write P1b vanilla Domain.{h,cpp} edits — would stack unbuildable/
  unverifiable code on the un-rebuilt P1a fixes. Resume order when classifier back:
  (1) rebuild OpenSeesPy, (2) re-run test_adr39 (confirm green), (3) commit P1a fixes
  + ledgers + P1b design, (4) write+build+test P1b.
- WORKING-TREE UNCOMMITTED (durable via this doc + git): P1a fixes
  [FEM_ObjectBrokerAllClasses.cpp, tcl/commands.cpp, LadrunoContactFE.{h,cpp}],
  ledgers [LEDGER_implementations.md, LEDGER_vanilla_files.md], designs
  [_adr39_p1b_design.md], this loop-state.

### Iteration 7 — classifier back; P1a fixes committed (5eb3b810); P1b WRITTEN + building
- Rebuilt OpenSeesPy w/ P1a fixes (bzxri5l4l, exit 0); test_adr39 3/3 green; COMMITTED 5eb3b810.
- P1b code written (per build-ready design):
  - NEW `SRC/domain/contact/`: LadrunoContactDomain.{h,cpp} (surfaces + contact defs +
    commit/revert counters; buildAdapterCount) + LadrunoContactSurface.{h,cpp} + CMakeLists
    (OPS_Domain) + add_subdirectory(contact).
  - VANILLA Domain.{h,cpp}: theContactDomain ptr (declared LAST → no -Wreorder; init all 4
    ctors; delete in dtor + clearAll-wipe mirroring ADR-30 theEQs fix; setter/getter; commit()
    + revertToLastCommit() hooks drive the engine).
  - Handler: inject buildAdapterCount() adapters (0 if no engine = pure Plain).
  - PARSER: contactSurface / contact / ladrunoContactInfo — OPS_ fns (OpenSeesOutputCommands.cpp)
    + decls (OpenSeesCommands.h) + Py + Tcl wrappers (dual-wired per the P1a code-gate lesson).
  - TEST: +4 P1b cases (defined-contact bitwise+commit-fires CDL/Newmark, wipe-clears, null=pure-Plain).
  - Stamped 4 contact files. Build bo5f2u6cy live.
- NEXT on build done: full test_adr39 (P1a 3 + P1b 4 = 7) → light review → ledger P1b vanilla
  rows (Domain.{h,cpp}, interpreter ×4) → commit P1b → P2 (NTS penalty narrow phase).

### Iteration 8 — P2 full design gate (wp3cr60mf) → SALVAGEABLE, folded
4-lens mechanics gate caught a SELF-CONTRADICTION + 2 silent-wrong BLOCKERs before any C++:
- **Q-P2-tan (decisive):** main-term-only K_c CANNOT pass the design's own FD-on-rotated
  gate — ∂n/∂u≠0 for a FLAT segment under rigid rotation (n rotates w/ nodes; what
  vanishes for flat is CURVATURE, not ∂n/∂u). FIX: include the ∂n/∂u block (honest
  Laursen/Wriggers NTS tangent, SYMMETRIC for frictionless), drop only O(gₙ·κ). Formula
  fix: main term = **kₙ BᵀB** not kₙ Bᵀ(n⊗n)B (B already carries n).
- **BLOCKER-1 (silent-wrong):** outward normal must be DERIVED from master element
  centroid, not trusted from node winding (flip → contact silently passes through).
- **BLOCKER-2:** bounded projection Newton (cap 10 + detK guard + non-converge sentinel);
  SimpleContact3D's unbounded `while` would hang the Domain scan.
- MAJOR: concave-corner tie-break (max-penetration in-bounds, else edge/vertex avg normal);
  gₙ0 default = no-offset + ABORT-on-init-penetration (stress-free = opt-in); -kn auto
  DEMOTED to P2b (K/A/V has no Element API — cache at setDomain from material GP + nodal V);
  implicit active-set freeze-per-step (or NewtonLineSearch); getID immutable+asserted;
  implicit RE-PROJECTS per Newton iter, explicit freezes per step; first-order-exact is
  interior-only.
- REFERENCE (not 1e-12 abs): rigid-plane analytic g=P/kₙ rel-1e-8 [P2a]; hand-placed
  ZeroLengthContactASDimplex pair cross-check rel-1e-6 [P2b oracle]; Hertz convergence.
- **SPLIT MANDATORY: P2a** (rigid+inclined plane, -kn val, connectivity={slave}, B=nᵀ,
  no projection kernel) **→ P2b** (faceted+deformable+Hertz+auto-kn+∂n/∂u tangent FD).
- Kernel: header-only OpenSees-free `LadrunoContactKernel.h` (mirrors LadrunoJ2Kernel);
  FE stays in OPS_Analysis (Analysis→Domain). ContactMaterial3D has NO penalty tangent
  to lift (author fresh).
All folded into `_adr39_p2_design.md` (committed bae2456c → revised). NEXT: code P2a.
  P2a needs the custom-FE_Element connectivity pattern (study PenaltySP_FE/TransformationFE
  — how a constraint-handler FE_Element subtype sets myDOF_Groups + myID + returns tang/resid).

### Iteration 9 — P1 MERGED #345; P2a coded on fresh branch guppi/contact-p2a
- P1 merged to ladruno (squash cfd35f458); CI caught a missing manifest.yaml row
  (fixed: testbed/manifest.yaml LadrunoContactHandler row) → green → merged. LESSON:
  add the manifest.yaml row in the SAME commit as any new classTag (not just the LEDGER table).
- Branched fresh off ladruno (guppi/contact-p2a) per the stacked-PR pitfall.
- **Custom-FE_Element connectivity pattern NAILED (PenaltySP_FE):** `FE_Element(tag, 1, ndof)`
  + `myDOF_Groups(0) = node->getDOF_GroupPtr()->getTag()`; base `setID()` fills `myID` with
  the node's first `numDOF` equation numbers (so coupling the leading ndm translational DOFs
  is automatic). SIGN: `getTangent` returns **+K** (PenaltySP returns +alpha, NOT c1-scaled);
  `getResidual` = force driving toward the constraint. c1 RESOLVED: residual determines the
  answer (exact, c1-independent); bare K_c tangent → exact in statics (c1=1), ignored in
  explicit, right-answer-but-maybe-more-Newton-iters in implicit-dynamic. So bare K_c is SAFE.
- **P2a code (rigid analytical plane, inline math — kernel header deferred to P2b):**
  - LadrunoContactDomain: `addRigidPlane(tag, slaveSurf, p0, n, kn)` (n normalized) + RigidPlane
    struct + getNumRigidPlanes/getRigidPlane.
  - LadrunoContactFE: 2nd ctor `(tag, Node* slave, ndm, p0, n, kn)` → FE_Element(tag,1,ndm),
    binds slave DOF_Group; rigidPlaneGap()=n·(X+u−p0); getResidual=−kn·g·n (g<0); getTangent=
    kn·n⊗n (active). EMPTY ctor kept for P1b path.
  - Handler: per rigid plane → per slave node → bound adapter (ndm from node getCrds().Size()).
  - Parser: `contactPlane tag slaveSurf nx ny nz px py pz kn` (OPS_ + Py + Tcl, dual-wired).
  - Test test_adr39_contact_p2a.py: axis-aligned pen=P/kn (static, rel 1e-6) + inclined (3-4-5
    off-axis n⊗n) + release→F≈0 + explicit impact stable/e≈1. Node STARTS just-penetrated
    (−1e-8·n) so static contact is active from step 1 (else singular — no out-of-contact stiffness).
  - Build bge579pp8 live.
- NEXT on build: run test_adr39_contact_p2a.py; if green → P2a code gate → ledger P2a vanilla
  rows (interpreter contactPlane) + manifest? (no new classTag, P2a rides 33002) → commit → PR
  (or stack) → P2b (faceted+deformable+Hertz+auto-kn+∂n/∂u tangent + ZeroLengthContactASDimplex oracle).

### Iteration 10 — P2a debug: build int-arg fix + STATIC PASS + explicit mass-pollution bug FIXED
- Build fix 1: `OPS_GetDoubleInput(&nd, ...)` needed int* count; had `double nd` → C2664. Fixed.
- **STATIC penetration=P/kn test PASSES** (1e-6) → adapter/penalty/connectivity/parser all correct.
- **EXPLICIT pass-through BUG (important, generalizes to any backing-element-less FE_Element):**
  node coasted THROUGH the plane at constant v. ROOT CAUSE: my getTangent returned K_c (contact
  stiffness) DIRECTLY; the `Linear` algorithm re-forms the tangent EVERY step, and CDL assembles
  getTangent into the explicit MASS (Diagonal) → M_xx = m + kn ≈ 1e6 → a_x ≈ 0 → node coasts.
  (getResidual/contact force were FINE; the huge wrong mass made them ineffective.)
  **FIX:** getTangent must route through `theIntegrator->formEleTangent(this)` (NOT return K
  directly), + override the virtual `zeroTangent`/`addKtToTang`/`addCtoTang`/`addMtoTang` to feed
  MY tang buffer. Then the INTEGRATOR decides: CDL formEleTangent = addMtoTang only (no-op) → tang=0
  → NO mass pollution; Newmark = addKtToTang(c1) → c1·K_c; statics = addKtToTang(1) → K_c. This
  ALSO resolves the c1-scaling correctly (better than the earlier "bare K_c" rationalization).
  LESSON for P2b/handoff: a custom FE_Element (myEle==0) MUST go through formEleTangent; returning
  K from getTangent corrupts the explicit mass matrix.

### Iteration 11 — explicit physics FIXED (e=1.0, pen=0.002 exact); 3rd bug = teardown null-deref
- After the formEleTangent fix: explicit impact CORRECT — penetration 0.002 (=v·√(m/kn) exact),
  restitution e=1.0. Static still exact. MECHANICS ALL CORRECT.
- Crash on a SECOND contact analysis in one process (first OK, 2nd ops.wipe crashes); NOT
  static-specific — ANY two contact analyses; impact-alone passes.
- **3rd bug (teardown null-deref, generalizes to any contact-ONLY model):** impact model has
  node+mass+contact adapter but NO element-backed FE_Element → the base FE_Element shared static
  scratch (theMatrices/theVectors), only allocated by the element-backed ctor, stays NULL → base
  ~FE_Element loops `theVectors[i]` when last FE destroyed (numFEs→0) → null deref → segfault on
  the 2nd wipe. P1 never hit it (Truss elements present). FIX: `ensureSharedScratch()` in both
  LadrunoContactFE ctors — guard-allocate if null. LESSON: a custom FE_Element in a possibly-
  element-free model must ensure the base shared scratch exists.
- Rebuild b92sm7ub5 live. NEXT: full test_adr39_contact_p2a.py (expect static + 2 impacts green).
- ensureSharedScratch guard FAILED to compile: theMatrices/theVectors are PRIVATE in FE_Element
  (not protected) → subclass can't touch them. REVERTED. Instead fixed at the TEST level: each
  P2a model carries a NEGLIGIBLE z-direction anchor truss (tiny EA, no test has z-motion → zero
  physics effect) which gives the model an element-backed FE_Element so the base shared scratch
  allocates. (A real contact model always has structural elements, so this is realistic, not a
  hack.) Rebuild bmtn93h21 live. NEXT: run the battery.
- KNOWN-LIMITATION for handoff: a true element-FREE contact model (single mass on a plane, no
  elements at all) still crashes on teardown (base FE_Element private shared-scratch null-deref).
  Options if needed: (a) make FE_Element theMatrices/theVectors protected (1-line vanilla edit,
  ledger), or (b) document that contact models need ≥1 structural element. Deferred — not on the v1 path.

### Iteration 12 — P2a ADVERSARIAL CODE GATE (wpqtt6tut) → SALVAGEABLE→PASS, fixes folded
Full multi-agent gate (4 source-grounded reviewers → adversarial verify each → synth;
25 agents, 20 findings, all verified against REAL upstream source). Verdict: SALVAGEABLE.
- **Hardest claim VERIFIED**: getTangent→formEleTangent(this) routing is correct — no
  double-count, no re-entrancy hazard, CDL truly mass-only, Newmark truly c1·K_c. Sign
  restoring (never attracts), residual/tangent active-set consistent within a Newton seq,
  parser n/p0/kn threading correct, 1e-6 tol conservative (not masking).
- **1 BLOCKER (B1) — element-free teardown null-deref** [FIXED at source]: a contact-ONLY
  model (mass + contactPlane, zero Domain Elements) segfaulted on teardown / 2nd analysis.
  ROOT CAUSE (real upstream bug): `FE_Element(tag,numDOF_Group,ndof)` subtype ctor did
  `numFEs++` BEFORE the `if(numFEs==0)` scratch-alloc guard → guard dead → theMatrices/
  theVectors never allocated when the only FEs are subtype adapters → `~FE_Element` derefs
  null. FIX = move `numFEs++` BELOW the guard (mirrors the element-backed ctor; zero
  regression — real models have an element-backed FE first). Better than the gate's
  "make members protected" workaround. Un-anchored `test_p2a_element_free_teardown` added
  (runs the model TWICE → the exact repro); the anchor-truss was MASKING this from CI.
- **In-scope folds** (all applied): kₙ≤0 guard + surface-Kind!=SLAVE_NODES guard in
  `addRigidPlane`; ndf<ndm skip-guard in the handler (silent equation-0 mis-assembly + OOB
  trial-disp read); `addKiToTang` override mirroring addKtToTang (else contact stiffness
  vanishes under Newton -initial / ModifiedNewton / HALL_TANGENT — residual still correct
  so never silently-wrong, just degraded). `test_p2a_inclined_penetration_static` added.
- **Ledgers**: LEDGER_vanilla_files += FE_Element.cpp row; LEDGER_implementations P2a row
  rewritten (was P1 zero-force) + PR #346.
- **Deferred to P2b (documented scope)**: no implicit active-set freeze (P2a = monotone-
  penetrating de-risk rung — anti-chatter is a P2b/implicit-ship deliverable); operator-
  splitting integrators (AlphaOS) unsupported (getK_Force is a zero stub); null-DOF_Group
  bind (unreachable in current handle() order). Static release→F=0 NOT added as a test —
  out-of-contact penalty model is rank-deficient/singular; release is covered by the
  explicit rebound tests (active set on→off, e≈1).
- **NEXT**: rebuild OpenSeesPy on guppi/contact-p2a → run test_adr39_contact_p2a.py (expect
  static + static-inclined + 2 impacts + element-free-teardown green) + p1 regression →
  stamp_headers (no new files) → commit fixes to #346 → P2b.

### Iteration 13 — gate fixes BUILT + 12/12 green; #346 AUTO-MERGED mid-work → recovered via #350
- Incremental build (silly-raman worktree) exit 0; **full battery 12/12 GREEN** (5 P2a incl.
  element-free-teardown + static-inclined + 7 P1 regression). B1 fix CONFIRMED: the un-anchored
  element-free model runs TWICE in one process with no segfault.
- Committed fixes to guppi/contact-p2a (7b3fa11cd) and pushed — but **PR #346 had AUTO-MERGED
  (squash, f8e614b44) while the code gate was running**, so 7b3fa11cd STRANDED on the merged
  branch (classic [[feedback_stranded_commits_after_automerge]]: d9ec2010d NOT an ancestor of
  ladruno; squash ⇒ branch orphaned). LESSON RE-CONFIRMED: re-check `gh pr view <n> --json state`
  IMMEDIATELY before any follow-up push; on this fork the window between "build" and "push" is
  enough for an auto-merge.
- RECOVERY: fresh branch `guppi/contact-p2a-gatefix` off origin/ladruno + cherry-pick 7b3fa11cd
  (clean, 9 files; the only post-#346 commit #349 just renamed an ADR doc — no overlap) → PR
  **#350** (base ladruno). P2a (rigid plane) is now FULLY SHIPPED once #350 merges (#346 base +
  #350 gate fixes). silly-raman's guppi/contact-p2a branch is now orphaned/ignore it.
- **NEXT = P2b** (the deformable mechanics, where the design-gate BLOCKERs live): faceted-master
  projection + 2 deformable LadrunoBrick blocks + Hertz + `-kn auto` + SOFT floor + ∂n/∂u tangent
  + bounded projection Newton + normal-from-element-centroid. Oracle = hand-placed
  ZeroLengthContactASDimplex pair (rel 1e-6). Verify #350 merged before stacking. See
  contact_p2_handoff.md "P2b" section.

### Iteration 14 — P2a MERGED (#350); P2b STARTED: plan + numpy oracle (7/7) BEFORE C++
- #350 squash-merged to ladruno (P2a rigid-plane contact + gate fixes now fully shipped).
  Branched guppi/contact-p2b off the merged ladruno (P2a verified present).
- **P2b SPLIT into sub-rungs** (`_adr39_p2b_design.md`): **P2b-1** = projection kernel +
  single slave vs single FIXED master segment (real bilinear/linear projection + derived
  outward normal + faceted B + kₙBᵀB(+∂n/∂u) tangent; `-kn $val`; ASDimplex node-pair oracle).
  **P2b-2** = deformable master + `-kn auto` + SOFT floor + Hertz + FD-on-rotated tangent gate.
- Precedent recon (Explore): SimpleContact3D (bilinear projection/tangent algebra + the
  unbounded-`while` :600-635 to FIX with a bounded Newton); ZeroLengthContactASDimplex ctor
  `(tag, masterNd, slaveNd, Kn, Kt, mu, ndm, itype, xN,yN,zN)` = the node-pair ORACLE;
  ContactMaterial3D tangent is LAGRANGE (author penalty tangent fresh); LadrunoBrick exposes
  `materialPointers[gp]->getInitialTangent()` + `getCharacteristicLength()=cbrt(V)` for -kn auto.
- **P1b parser ALREADY has the surface/contact plumbing**: `contactSurface $tag -master $nps
  $n1..` stores nodesPerSeg; `contact $tag $masterSurf $slaveSurf $kn $kt $mu` parses kn. So
  P2b-1 needs mainly: kernel + surface coord/normal cache + FE segment mode + handler + tests.
- **NUMPY ORACLE `contact_prototypes/proto_p2b_nts.py` — 7/7 PASS** (oracle-first, the fork
  discipline; no build): interior pen=exact, self-equilibrium |ΣF|=0/|ΣM|=1e-12, **winding-flip→
  identical force** (BLOCKER-1), oblique-30° normal, out-of-bounds→0, tangent slave-block FD==
  +kₙn⊗n (1.5e-11) + FULL tangent SYMMETRIC (4e-11, frictionless ✓) + ∂n/∂u rel-weight 3.6%
  (non-negligible — main-term-only fails a tight FD gate, as the design gate predicted), tri-3+quad-4.
  KEY: **OpenSees tangent convention K_assembled = −∂r/∂q** (verified vs P2a: r_s=+tₙn,
  addKtToTang=+kₙn⊗n) — kernel main term = +kₙBᵀB. Math LOCKED; C++ transcription de-risked.
- **NEXT (P2b-1 C++)**: transcribe to header-only `SRC/domain/contact/LadrunoContactKernel.h`
  (pure fns, mirror LadrunoJ2Kernel) → LadrunoContactSurface setDomain coord/normal cache →
  LadrunoContactFE SEGMENT ctor (project per-iter implicit/per-step explicit via kernel) →
  handler builds segment adapters → `tests/test_adr39_contact_p2b.py` (ASDimplex oracle rel-1e-6,
  self-eq, winding-flip, oblique, oob) → build OpenSeesPy → code gate → PR (base ladruno). The
  full ∂n/∂u analytic block + FD-on-rotated gate lands in P2b-2 (deformable/rotated master);
  P2b-1's fixed master makes the slave-block main term sufficient for its solve.

### Iteration 15 — P2b-1 C++ WRITTEN (faceted NTS, fixed master), building
- Transcribed the validated oracle to **NEW header-only `SRC/domain/contact/LadrunoContactKernel.h`**
  (OpenSees-free pure fns: shape tri-3/quad-4, BOUNDED projection Newton w/ degenerate-segment
  reject, DIRECTION-oriented winding-immune normal, gap, penalty). Stamped + stamp_headers glob.
- `LadrunoContactFE` SEGMENT mode (3rd ctor): conn = {slave}∪{seg nodes} (ndof=3*(1+nps)), per-eval
  projection via kernel, resid=Bᵀtn, tangent=kn·BᵀB (routed thru formEleTangent like P2a;
  ∂n/∂u block deferred to P2b-2 — fixed master ⇒ slave-block kn·n⊗n exact). setID maps each
  ndf==3 node's DOF_Group → 3 myID slots (verified FE_Element::setID concatenates full group IDs).
- `LadrunoContactHandler`: one adapter per (slave node, master SEGMENT) pair — gate-sanctioned
  BRUTE-FORCE pairing (bucket sort = P2.5); guards kind + ndf==3 + degenerate; computes orientDir.
- `LadrunoContactDomain`: Contact struct made public + `getContact(i)` + optional `outward[3]`.
- Parser: `contact ... -outward ox oy oz` (orientation dir toward allowed half-space; robust to
  just-penetrated starts; peek/un-read via OPS_ResetCurrentInputArg(-1) keeps kn kt mu positional).
- **ORIENTATION = explicit DIRECTION** (not slave-position): auto-orient-toward-slave-pos FAILS for
  a just-penetrated start (slave is on the forbidden side by −1e-8 → normal flips). -outward fixes it
  robustly + is winding-immune. Handler auto-derives orientDir = slave_ref − seg_centroid when no
  -outward (for clearly-separated starts).
- **P1b REGRESSION SURVIVES UNTOUCHED:** its `contact` defines a COLLINEAR master (nodes 1,2,3 on
  x-axis) → degenerate → kernel projection sentinel returns -1 → zero force; numberer Plain is
  connectivity-independent ⇒ still bitwise. (Made `contact` load-bearing without breaking P1b.)
- TEST `tests/test_adr39_contact_p2b.py`: static pen=P/kn (quad+tri, 1e-6), oblique-30°, winding-flip
  immunity, OOB pass-through, explicit e≈1, **ZeroLengthContactASDimplex cross-check** (rel 1e-6).
- BUILD: cold worktree (peaceful-bassi) hit the MUMPS-download gotcha (no local archive + curl reset)
  → copied mumps_src.tar.gz from compile-root → rebuilding (bakbc6n5c). LESSON: a fresh worktree needs
  `mumps-archive/mumps_src.tar.gz` copied in before build.bat (offline-safe).
- **NEXT on build**: run test_adr39_contact_p2b + p2a + p1 (expect all green) → P2b-1 code gate →
  ledger (kernel new-file row + P2b in the contact row) → commit → PR (base ladruno) → P2b-2
  (deformable master + -kn auto + Hertz + FD-on-rotated ∂n/∂u tangent gate).

### Iteration 16 — P2b-1 BUILT + 19/19 green; code gate (wxh6f8ki0) PASS, fixes folded
- Cold worktree build exit 0 (after copying mumps_src.tar.gz). First run: 16/19 — caught a
  TEST bug (static tests freed the WRONG DOF: `fix(1,0,1,1)` frees x, but normal/load are +z →
  singular matrix; winding-flip + ASDimplex passed VACUOUSLY comparing two frozen 1e-8 values).
  FIX = `fix(1,1,1,0)` (free z) + oblique redesigned to free-z-only (pen=P/(kn·nz²)). **19/19 green**
  (7 P2b + 5 P2a + 7 P1). C++ was correct; only the test DOF-freeing was wrong.
- Committed P2b-1 (5eafac6ae) → **adversarial code gate (4 reviewers→verify→synth, 17 agents):
  verdict PASS** (no BLOCKER, no MAJOR). Kernel faithful to oracle; tangent-drop sound for FIXED
  master; P1b-regression-preservation PROVEN; no vacuous tests remain (re-audited).
- **Folded (gate-recommended, in-scope):** H1 = addContact reject kn<0 (kn==0 still OK for P1b
  zero-force); PARSE-2/H1 = handler warn+skip kn<=0 SEGMENT contact; H2 = kernel normalOriented
  FAIL-SAFE on ambiguous orientation (|n·refDir|≈0 → false, refuse rather than guess a sign) +
  NEW `test_p2b_auto_orientation_no_outward` covering the previously-untested auto-orient path;
  KMF-3 = traction() internal Macaulay clamp (gap>=0 → 0, no adhesive foot-gun).
- **Deferred to P2b-2 (gate-scoped):** TANG-1 = free-master ∂n/∂u consistent-tangent block (the
  P2b-2 deliverable; optional P2b-1 diagnostic = warn on a free segment DOF). NITs deferred:
  PARSE-1 (peek-idiom hardening — works on all tested paths), PARSE-3 (unknown-option swallow).
- Rebuild b2bcxvto3 (warm) live. **NEXT on build**: re-run battery (expect 20/20 w/ the auto-orient
  test) → commit fold → PR to ladruno (verify #350 merged; this branch is stacked on it) → P2b-2.

### Iteration 17 — P2b-1 MERGED #354; P2b-2a deformable-master VALIDATED (no new C++)
- #354 (P2b-1) merged to ladruno (squash). Branched guppi/contact-p2b2 off it.
- **KEY EMPIRICAL FINDING (probe_p2b2_deformable.py):** the MERGED P2b-1 NTS code ALREADY
  handles a DEFORMABLE master — the FE/handler read master trial-disp + assemble the master
  tangent blocks, so a free master "just works" with the main-term kn·BᵀB tangent. Single
  fixed-bottom LadrunoBrick, top face = master, slave pressed down: static Newton CONVERGES,
  contact gap = P/kn EXACT, brick compression = PL/EA EXACT (series). So the deferred ∂n/∂u
  block (gate TANG-1) is a CONVERGENCE REFINEMENT (large-rotation/curved), NOT a correctness gap.
- **`tests/test_adr39_contact_p2b2.py` 2/2 green** (test the MERGED code, zero new C++):
  slave-on-deformable-brick (gap=P/kn + comp=PL/EA) + block-on-block (deformable-vs-deformable
  interface, both bricks compress, per-node interface pen=(P/4)/kn). Block-on-block needed the
  top block to START just-penetrated (else the unsupported block falls → diverge).
- **P2b-2 REORGANIZED by the finding:**
  - **P2b-2a = deformable-master validation — DONE** (this iteration; tests + probe; no C++).
  - **P2b-2b = `-kn auto`** (NEXT real feature): cache K/A/V at setDomain — K from a master
    `LadrunoBrick materialPointers[gp]->getInitialTangent()(0,0)`, V via getCharacteristicLength()
    =cbrt(V) or nodal coords, retain the master element ptr; kn=f_si·K·A²/V (26.14a) + SOFT floor
    (26.15). Needs a surface→master-element link (the surface stores node tags only today).
  - **P2b-2c = ∂n/∂u consistent-tangent block + Hertz** (the hard math; FD-on-rotated gate;
    convergence refinement for large-rotation/curved — oracle already has the FD ground truth in
    proto_p2b_nts.py, ∂n/∂u≈3.6%). Lower priority since explicit-ship never forms the tangent.
- **NEXT**: commit P2b-2a (tests + probe + plan) → PR to ladruno (validation-only, fast) → P2b-2b -kn auto.

### Iteration 18 — P2b-2b `-kn auto` design DECIDED (user gate); oracle next
- Resumed at P2b-2b (`-kn auto`) — verified #355 (P2b-2a) merged + is ladruno HEAD.
- **DATA-MODEL DECISION (user-gated, both recs accepted):**
  - **kn source = GENERIC element stiffness** (NOT LadrunoBrick-specific). `LadrunoBrick::
    materialPointers[]` is PRIVATE → the handoff's `materialPointers[gp]->getInitialTangent()`
    sketch is unreachable without a new accessor/dynamic_cast. Instead use base-`Element`
    virtuals only: auto-detect the owning solid element per segment by face-node-subset match
    (handler has the Domain at handle()), pull `getInitialStiff()` (24×24 for a hex), and
    reduce **kn = f_si · mean_over_face_nodes( nᵀ K_block_node n )** where K_block_node is the
    3×3 diagonal block at that node's DOFs and n is the segment normal. Works for ANY solid
    element (LadrunoBrick/SSPbrick/stdBrick); NO element-type coupling; NO vanilla edit.
  - **Dimensionally ≡ LS-DYNA 26.14a** `f·K·A²/V`: for a 3D solid of size L, modulus E,
    K_diag ~ E·L and A²/V ~ L⁴/L³ = L ⇒ E·A²/V ~ E·L ~ K_diag (constant absorbed in f_si).
    Documented as the OpenSees adaptation; f_si default 0.10 (LS-DYNA SLSFAC).
  - **SOFT Courant floor (26.15) DEFERRED to P2b-2c** (needs nodal mass + Δt_c; mainly
    explicit anti-chatter — keep this rung tight + gate-able).
- Resolution point = `handle()` (has Domain → element lookup). Parser parses literal `auto`
  in the kn slot → `knAuto` flag on Contact; handler resolves per pair. Reject auto for the
  rigid plane (P2a — no deformable master). `-kn $val` path unchanged.
- **NEXT**: oracle `contact_prototypes/proto_p2b2b_autokn.py` (reduction arithmetic on a
  synthetic SPD block + scaling laws kn∝E, kn∝L, g=P/kn small) → C++ → build → test → gate → PR.
- ORACLE `proto_p2b2b_autokn.py` **6/6 PASS** (assembles a real trilinear-hex K in numpy =
  `LadrunoBrick -formulation std` getInitialStiff; validates reduction arithmetic T1, kn∝E
  T2, kn∝L T3, penetration g/L<10% T4, oblique-normal positivity T5, mesh-objectivity T6).
- **C++ WRITTEN** (cold worktree build live, b11fgz5d7; copied mumps_src.tar.gz first):
  - `LadrunoContactDomain.{h,cpp}`: Contact struct += `bool knAuto`; `addContact(...,knAuto=false)`.
  - PARSER `OpenSeesOutputCommands.cpp OPS_LadrunoContact`: kn slot accepts literal `auto`
    (peek/classify; numeric `kn kt mu` path unchanged; `auto` may carry optional kt mu).
  - HANDLER `LadrunoContactHandler.cpp`: file-static `ladrunoResolveAutoKn(dom,segNodes,nps,
    orientDir)` — ref normal via kernel `normalOriented` at seg center, auto-detect owning
    solid (ElementIter: getNumDOF==3·nNodes + getExternalNodes ⊇ seg tags), reduce
    kn=f_si·mean(nᵀ getInitialStiff()_block n), f_si=0.10. Segment loop: `!ct.knAuto && kn<=0`
    skip guard + per-pair resolve (knUse; skip+warn on <=0). NO new vanilla file, NO ele coupling.
  - TEST `tests/test_adr39_contact_p2b2b.py` (4): converges+usable, kn∝E, ABSOLUTE P/pen==oracle
    (imports proto_p2b2b_autokn — single source of truth), no-owning-solid→skip(spring-only).
  - LEDGERS updated (impl row P2b-2b notes+files+#355 fix; vanilla P2b-2b parser-keyword row;
    `// Ladruno ADR-39 P2b-2b` grep-marker in the parser). No stamp_headers change (no new SRC file).
- **NEXT on build done**: run test_adr39_contact_p2b2b + p2b2 + p2b + p2a + p1 (expect all green) →
  light code gate (heuristic, no novel math/core edit) → commit → PR base ladruno (verify #355 merged).
- **BUILD ✓ (b11fgz5d7, exit 0)** cold worktree; re-wired the site-packages `.pth` to THIS
  worktree's dist via `Ladruno_scripts/wire_venv_pth.py <wt>/dist/bin` (each fresh worktree has
  its own dist → must re-wire before pytest; `import opensees` else fails). **BATTERY 26/26 GREEN**
  (7 P1 + 5 P2a + 8 P2b-1 + 2 P2b-2a + 4 P2b-2b). One test self-correction: the "reasonable"
  test mis-scaled the load (P=100 on E=2e4 = 0.5% pressure → 21% penetration; the heuristic is
  CORRECT — penalty contact at f_si=0.10 ≈ LS-DYNA's soft default, penetration ∝ pressure/modulus)
  → fixed to load a representative 0.1%-modulus pressure (P=p·E·L²) + threshold pen/L<0.10. The
  ABSOLUTE oracle-match test (P/pen == imported proto_p2b2b auto_kn, rel<2%) passed unchanged.
- **CODE GATE (1 focused source-grounded reviewer): NO BLOCKER, NO MAJOR.** Core verified sound:
  PSD diagonal block ⇒ nᵀKn ≥ 0 (negative kn impossible), nᵀKn sign-insensitive (ref-vs-deformed
  normal sign disagreement is harmless), generic [3·loc,3·loc+3) DOF-block extraction valid for any
  3-DOF/node solid, getNumDOF()==3·nn filter CANNOT match a non-solid in a 3D model (truss has 2
  nodes → face-node-subset fails), getInitialStiff() static-buffer use alias-safe (read before any
  re-entrant call), shared ElementIter re-reset safe (outer ele loop fully drained first), knAuto
  default-arg init clean, skip-guard + per-pair resolve-failure `continue` correct. **3 MINOR:**
  (1) FOLDED — parser silently swallowed stray trailing tokens (also a deferred P2b-1 NIT) → now
  errors on an unexpected token (`// Ladruno ADR-39 P2b-2b` gate MINOR-1). (2) DEFERRED to P2b-2c —
  no floor on a tiny-positive kn (exactly what the SOFT Courant floor addresses). (3) DOC-ONLY —
  first-match owning element is order-dependent for an INTERIOR/shared face (correct for the
  intended EXTERIOR contact-master face; single owner). Gate-fix rebuild bo8d7gkui live.
- **NEXT on gate-fix build**: re-run battery (expect 26/26) → commit P2b-2b → PR base ladruno
  (re-check `gh pr view 355 --json state`==MERGED immediately before push; fork auto-merges fast).

- **BUILD GOTCHA (new, important):** after `wire_venv_pth.py` points the site-packages `.pth`
  at a VALID dist, EVERY `python` startup eagerly `import opensees` (the openseespy alias in the
  boot) → prints the LADRUNO banner to stdout. CMake's `find_package(Python)` probe captures the
  banner instead of a version string → "Could NOT find Python: Found unsuitable version <ASCII art>"
  + a CMake syntax error → configure FAILS. FIX = `set LADRUNO_OPENSEES_QUIET=1` (user env) BEFORE
  building whenever the .pth is wired (cold build worked only because its .pth pointed at a dist
  where `import opensees` failed silently). Candidate for BUILD_GOTCHAS.md / LEDGER_quirks.

### Iteration 18 RESULT — P2b-2b `-kn auto` SHIPPED-pending-PR (battery 26/26, gate PASS)

### Iteration 19 — P2.5 bucket-sort broad phase (user-chosen next rung); coded, building
- User picked P2.5 (bucket sort) over P3/P2b-2c. Replaces the handler's brute-force
  O(nSlave·nSeg) pairing with the §26.11 spatial bucket sort. Algorithm already validated
  oracle-first in `proto_bucket_sort.py` (P0a 4/4: superset contract 0-miss, jittered mesh,
  runaway guard recovers, 33× pruning) — re-confirmed PASS this iteration.
- Design `_adr39_p25_design.md`. **C++ WRITTEN** (building on top of the unmerged P2b-2b
  commits in THIS worktree; will cherry-pick onto a fresh ladruno branch once #357 merges):
  - NEW header-only `SRC/domain/contact/LadrunoContactBucketSort.h` (OpenSees-free Grid class:
    median-diag cell, bbox + runaway percentile clip, cap NX·NY·NZ≤min(nSeg,5000), segment-SPAN
    registration type-13, 27-neighbour candidates() with sparse-set stamp de-dup). Stamped + glob.
  - Handler segment loop: build Grid from master seg REFERENCE coords (missing-node backfill =
    superset-preserving), per slave `grid.candidates()` → loop candidates not all nSeg. Narrow
    phase UNCHANGED ⇒ result identical for any kept pair.
  - `-cell <frac>` knob (Contact.cellFrac, parser, default 1.0): a HUGE value ⇒ 1 bucket ⇒ every
    segment candidate = BRUTE FORCE (the equivalence-gate mechanism). Also a real distorted-mesh knob.
  - TEST `tests/test_adr39_contact_p25_bucketsort.py`: **verify==brute force** (5×5=25-seg grid,
    default cell vs `-cell 1e12`, bitwise-identical slave disps — the gate) + interior-segment
    correctness (4×4) + single-segment regression.
- **SCOPING (documented in the header):** grid built at handle() from REFERENCE coords ⇒
  bitwise == brute force for STATIC/SMALL-MOTION (the validated regime + gate meshes). LARGE
  sliding needs epoch re-emit (re-handle on membership change; LS-DYNA re-sorts every 10–15
  cycles) — the follow-on rung. Build ba8eujt42 live (LADRUNO_OPENSEES_QUIET=1 set User-env).
- **NEXT on build**: run test_p25 + full battery (expect 26+3 green) → code gate → commit P2.5
  as ONE separable commit → after #357 merges, fresh branch off ladruno + cherry-pick → PR.
- #357 (P2b-2b) status when last checked: OPEN, CLEAN/green (Zone-A passed), not yet auto-merged.
- **BUILD ✓ (ba8eujt42); BATTERY 29/29 GREEN** (3 P2.5 + 26 prior). The verify==brute-force
  gate passed: default bucket grid gives BITWISE-identical slave disps to `-cell 1e12` (brute
  force) on a 5×5=25-seg grid ⇒ broad phase drops no near pair.
- **CODE GATE (1 focused reviewer, 250k-probe superset re-derivation): NO BLOCKER/MAJOR.**
  Superset guarantee PROVEN (penetration ⇒ slave in segment corner-AABB ⇒ in a span-registered
  cell; ±1 neighbour = extra margin; 0 misses across skewed/ramped/mixed meshes). clampP edge
  cases sound, linId in-range, cap-shrink-before-registration consistent, degenerate-seg floor
  ok, missing-node backfill cannot leak HUGE_VAL, `-cell 1e12`⇒1 bucket truly brute force, out[]
  no overflow (dedup), Grid↔handler seg-index aligned. **2 MINOR FOLDED:** (1) de-dup counter
  `stampTick_`/`stamp_` widened int→long long (UB-after-2³¹ hardening, unreachable but clean);
  (2) design-doc reconciled — code uses cell=lmaxFrac·median_diag (no r_search term); harmless
  for the static/small-motion scope (gate-proven). Rebuild bm9q08wdk live; re-test → commit.
- **COMMIT PLAN:** P2.5 commits stack on the unmerged P2b-2b (#357) in this worktree. Commit P2.5
  as ONE separable commit; once #357 merges, fresh branch off ladruno + cherry-pick → PR base ladruno.

### Iteration 19 RESULT — P2.5 MERGED #358 (Zone-A confirmation in flight)
- #357 (P2b-2b) was green/CLEAN but the fork auto-merge stalled ~45 min → user approved
  `gh pr merge --squash --auto`. On THIS repo (no required status checks) `--auto` merges
  IMMEDIATELY when mergeable. #357 merged (ladruno c49de8773).
- Branched fresh `guppi/contact-p25-bucketsort` off origin/ladruno (now incl. P2b-2b) +
  cherry-picked the P2.5 commit (clean, 11 files) → PR **#358** → enabled `--auto` → MERGED.
- **LESSON (`--auto` on this repo = merge-NOW):** with no required status checks, `gh pr merge
  --auto` merges the instant the PR is mergeable — it does NOT wait for Zone-A. #358 merged while
  Zone-A (Ubuntu) was still PENDING ⇒ P2.5 landed in ladruno before Linux verification. For #357
  I'd waited for Zone-A first; for #358 I enabled --auto early. FUTURE: wait for Zone-A green
  THEN merge (or accept the merge-before-Linux risk + fix-forward). Monitoring Zone-A on #358
  (buinif8wt); code is portable header-only STL + local battery 29/29, so low risk, but this
  fork has had Win→Linux drift (dsygv_ incident) — confirm green, fix-forward if red.

### Iteration 19 FINAL — P2.5 Zone-A GREEN (9m7s); fully shipped, no Win→Linux drift.

### Iteration 20 — P3 IMPL-EX Coulomb friction (the v1 SHIP) STARTED: oracle done
- User: "lets do it" → P3 (the v1 ship, critical junction). ladruno confirmed sound (P2.5 Zone-A green).
- Read the reference `ZeroLengthContactASDimplex::updateInternal` (damage-formulated return map +
  its quarantined bugs) — Ladruno authors the CLEAN Coulomb radial return instead (ADR grounding).
- **ORACLE `proto_p3_friction.py` 6/6 PASS** (friction math LOCKED): stick (‖tT‖<μN, no slip),
  slip caps at μN + correct flow dir, dissipation≥0, **non-symmetric** consistent tangent vs FD
  (∂tT/∂gT=(μN·kt/‖tT*‖)(I−n̂n̂) + the ∂tT/∂gN=−μ·kn·n̂ pressure-coupling column with NO symmetric
  partner = rigorous non-symmetry), IMPL-EX→implicit in steady slip (6.8e-13), sliding-block-on-
  incline a=g(sinθ−μcosθ) exact. Oracle caught 2 of my bugs (dropped-kt in ktan; FD ground truth)
  + surfaced 2 design facts: (a) tangent non-symmetry is real (∂tT/∂gN coupling), (b) IMPL-EX
  OVERSHOOTS one step at stick→slip onset (dlam_old still pre-slip 0) then re-syncs — documented.
- Wrote `_adr39_p3_design.md` + ran the **adversarial DESIGN gate (2 source-grounded reviewers:
  state/lifecycle + mechanics/signs) → SALVAGEABLE**, all 9 findings FOLDED into the design:
  - **BLOCKER (sign):** `addB` sums the residual (no flip); `+tT` on the slave ACCELERATES it →
    `a=g(sinθ+μcosθ)` energy injection. FIX: friction OPPOSES `n̂*`, kernel returns already-negated
    applied force; **incline test MUST run through real CDL+addB** (not the oracle's hand-subtraction
    — that's how the sign bug would ship green). The single most important guardrail.
  - **MAJOR (engagement ref):** global-reference `gT` corrupts late-engaging contact (pre-contact
    drift → spurious stick); store `gT0` at first activation.
  - **MAJOR (DROP IMPL-EX from v1):** explicit discards the tangent so IMPL-EX buys nothing but costs
    the onset-overshoot impulse; use the DIRECT return map (exact, no overshoot, removes dlam state) —
    also dissolves the firstStep-double-getResidual BLOCKER. IMPL-EX kept in header for P3.5 only.
  - **MAJOR:** dead-slot GC each handle() (ADR-30 theEQs leak class); mu≤0 short-circuit before any
    slot touch (byte-identical P2b + dodges 0/0); explicit has no auto-retry → retry wrapper must revert.
  - **MINOR:** segIndex=GLOBAL ordinal (not bucket candidate ci); cross-contact shared-slave warn;
    adapter lazy-refetches engine via Domain* (wipe deletes engine→dangling); current-N for v1
    (committed-N cap=P3.5); near-zero guard physically scaled; segIndex-flip slip-loss + kt≤kn note.
  - Verdict: core architecture (engine-keyed state, commit-hook, clean return map) SOUND; design HARDENED.
- **NEXT = code P3** (kernel direct Coulomb return map returning the NEGATED applied force +
  engagement gT0 + engine per-pair state {gpT,gT0,engaged} + GC + mu≤0 short-circuit + adapter
  lazy engine-refetch) → build → test (INCLINE through real CDL+addB = the gate; stick/slip/energy/
  frictionless-regression/commit-revert/wipe-reanalyze) → code gate → PR base ladruno. Branch fresh
  off ladruno (it has P2.5; HEAD after #358 merge). Oracle `proto_p3_friction.py` 6/6 is the math ref.

### Iteration 21 — P3 friction C++ WRITTEN (explicit FORCE-only ship), cold worktree building
- New worktree `vigilant-solomon-8c93bd` (branch guppi/vigilant-solomon-8c93bd, off ladruno
  HEAD fe83cbe40 = #359 P3 design). Cold build (no dist/build/mumps) — copied
  mumps_src.tar.gz in + set LADRUNO_OPENSEES_QUIET=1 (User env) before build.bat.
- **P3 C++ transcribed from the gate-hardened design + oracle (proto_p3_friction.py 6/6):**
  - **Kernel** `LadrunoContactKernel.h`: `tangentPart()` + `frictionReturnMap()` (clean
    Coulomb DIRECT return map; returns NEGATED applied force `tFric=−tT` — SIGN in ONE place;
    `gpTtrial` PURE fn of committed ⇒ firstStep-double-eval idempotent; near-zero guard is
    physically scaled = slip needs ‖tT*‖>cap>0 so n̂ always normalizable). IMPL-EX stays in
    the oracle only (explicit discards tangent).
  - **Engine** `LadrunoContactDomain.{h,cpp}`: `FrictionState{gpT,gpTtrial,gT0,engaged}` in a
    `std::map<PairKey(contactTag,slaveTag,segIndex),...>`; `getOrCreateFrictionState` (lazy);
    real `commit()`(gpT=gpTtrial)/`revertToLastCommit()`(gpTtrial=gpT) REPLACING the P1b
    counters (counters kept for the info command); `frictionGCBegin/Mark/End` (live key-set
    GC each handle, ADR-30 theEQs leak class). `getNumFrictionStates()`.
  - **Adapter** `LadrunoContactFE.{h,cpp}`: SEGMENT ctor += `kt,mu,Domain*,contactTag,segIndex`;
    `segmentActive` += optional `gTvec` (tangential rel-pos = (xs−x̄)_⊥n at the SAME projection);
    getResidual friction block — `mu>0 && theDomain` guard (mu≤0 SHORT-CIRCUITS ⇒ byte-identical
    P2b), lazy `getLadrunoContactDomain()` re-fetch (wipe deletes engine), capture `gT0` at first
    activation, `gTeff=gTvec−gT0`, return map w/ N=tn, write gpTtrial, MIRROR normal block
    (slave+=tFric, master_i+=−N_i tFric). Tangent UNCHANGED (friction tangent=P3.5).
  - **Handler** `LadrunoContactHandler.cpp`: thread `ct.kt,ct.mu,theDomain,ct.tag,seg` into the
    SEGMENT ctor; `frictionGCBegin` + per-frictional-pair `frictionGCMark` + `frictionGCEnd`
    bracketing the contact loop; cross-contact shared-slave multiset warning.
  - No NEW vanilla file (Domain commit/revert hooks already wired P1b); no new classTag
    (rides HANDLER 33002); no parser change (kt/mu already parsed P1b). No new authored file
    (stamp_headers --check: all 146 current).
  - **TEST** `tests/test_adr39_contact_p3_friction.py` (7): incline a=g(sinθ−μcosθ) through
    REAL CDL+addB (THE sign gate; +μ=0 leg = g·sinθ), slip-caps-μN (a=(Q−μN)/m, f=μN, two
    drives), stick (Q<μN held), dissipation (v_fric/v_free≈(Q−μN)/Q), frictionless regression
    (mu=0 vs 0.5 pen=P/kn identical), commit-fires (numCommits==steps), wipe-reanalyze (GC+lifetime).
  - Ledger: LEDGER_implementations contact row += P3 paragraph + P3 test file + #358/#359 PR fixes.
- BUILD bqw9r1yfp live (cold, configure passed = QUIET flag worked, ~640/1871 at last check).
- **CRITICAL SELF-CAUGHT BUG (pre-test, oracle-style numpy check):** my first gTvec used
  POSITIONS (x_s − x̄)_tan — but the closest-point projection makes (x_s − x̄) ∥ n, so that
  is IDENTICALLY ZERO ⇒ zero slip ⇒ NO friction. FIX = displacement-based slip
  `d = u_s − Σ N_i u_i`, tangential part (master DISPLACEMENTS, not positions). Verified
  numerically (pos→0 always, disp→δ) before rebuild.
- **BUILD ✓** cold (bqw9r1yfp exit 0, QUIET flag → CMake Python probe clean) + incremental
  rebuild after the slip fix. Re-wired .pth to this worktree dist.
- **BATTERY 36/36 GREEN** (7 P3 + 7 P1 + 5 P2a + 8 P2b-1 + 2 P2b-2a + 4 P2b-2b + 3 P2.5).
  P3 7/7: incline a=g(sinθ−μcosθ) (sign gate, +μ=0 leg g·sinθ), slip-caps-μN (2 drives),
  stick, dissipation (v ratio 0.5), frictionless regression (mu=0 vs 0.5 pen identical),
  commit-fires (numCommits==50), wipe-reanalyze. NO regression in the 29 prior.
- **P3 CODE GATE (3 parallel source-grounded adversarial reviewers): 2× PASS + 1 SALVAGEABLE.**
  - **sign/mechanics → PASS:** traced the SIGN end-to-end through REAL `addB`/CDL (could not
    prove it wrong) — kernel `tFric=−tT` opposes motion ⇒ a=g(sinθ−μcosθ); return map
    oracle-EXACT; cone N=tn safe; slip branch `norm>cap>0` ⇒ no 0/0; displacement-based slip
    correct (position-based ≡0 by closest-point ⊥), convective migration a documented limit.
  - **state/lifecycle → PASS** (all MINOR/NIT): firstStep idempotent (`gpTtrial` pure set of
    committed); commit/revert hooks fire; GC bracket correct + `erase(it++)` safe; mu≤0
    byte-identity; key = global `seg` (not `ci`), collision-free; lazy refetch wipe-safe; all
    3 ctors init in declaration order.
  - **tests → SALVAGEABLE** (the 2 load-bearing tests are genuinely strong; folded the gaps):
    + **MAJOR revert-path** — was untested (design MAJOR-7); ADDED `test_p3_revert_path`
      (failed implicit step → numReverts increments; verified the hook fires, probe confirmed).
    + **MAJOR deformable-master friction** — the `resid_master_i += −N_i·tFric` block was never
      exercised (all masters fixed); ADDED `test_p3_deformable_master_drag` (slave slides on a
      deformable LadrunoBrick top face → friction DRAGS the brick top +x; sign-decisive).
    + **MINOR late-engagement gT0** — the #1 design fix was untested (all tests pre-penetrate);
      ADDED `test_p3_late_engagement_gt0` (slave drifts +x while separated, then penetrates →
      a_x≈0 at engagement vs the −kt·Δ/m spurious-stick bug).
    + dissipation reframed as an independent ENERGY-BALANCE check (drive work = KE + μN·x);
      stick bound TIGHTENED to 1.2·(2Q/kt); frictionless-regression clarified (a statically-
      ACTIVE friction leg is INFEASIBLE in v1 — free tangential DOF is singular w/o the friction
      tangent, which is WHY the ship is explicit; active friction covered by slide/stick/incline).
- **BATTERY 39/39 GREEN** (10 P3 + 29 prior). No regression.

### Iteration 21 RESULT — P3 friction MERGED #360 (squash 0139fd33a, Zone-A green 8m43s)
- P3 (v1 explicit friction ship) merged to ladruno. **ADR-39 v1 contact ship COMPLETE:**
  P1 → P2a → P2b-1/2a/2b → P2.5 → P3 friction. classTag/manifest + Zone-A all green; no
  Win→Linux drift. User authorized pr+merge.

### Iteration 22 — P3.5 implicit frictional Newton (the consistent NON-SYMMETRIC tangent) STARTED
- Branched guppi/contact-p35-implicit-friction off origin/ladruno (has P3).
- **WHY:** P3 ships friction FORCE-only; under IMPLICIT (static/Newmark) a free tangential DOF
  is SINGULAR without the friction tangent (observed: FullGeneral U(0,0)=0 in P3 testing). P3.5
  assembles the consistent friction tangent so implicit Newton converges. The per-traction
  tangent blocks are already oracle-validated (proto_p3_friction.py friction_tangent, 6/6):
  STICK ∂tT/∂gT=kt·P_t, ∂tT/∂gN=0 (SYMMETRIC); SLIP ∂tT/∂gT=(μN·kt/‖tT*‖)(P_t−n̂⊗n̂),
  ∂tT/∂gN=−μ·kn·n̂ (the pressure-coupling column ⇒ NON-SYMMETRIC).
- **Assembled tangent (the new 3D part to pin + code):** K_fric = Gᵀ[D_TT·P_t + d_TN⊗n]G where
  G=[I|−N_i I] (rel-disp operator), P_t=(I−n⊗n) tangent projector, D_TT=∂tT/∂gT, d_TN=∂tT/∂gN.
  Lives in `addKtToTang` (implicit path only; explicit CDL's addMtoTang skips it ⇒ P3 explicit
  byte-identical). NON-SYMMETRIC ⇒ needs FullGeneral/UmfPack (document; cannot detect solver in FE).
- **Open design Qs:** (1) current-N (non-sym, rigorous) vs committed-N (∂tT/∂gN=0 ⇒ SYMMETRIC,
  robust, design's smoothness option) — pick/knob; (2) active-set per-Newton-iter re-projection vs
  freeze-per-step (chatter); (3) IMPL-EX symmetric-secant variant = P3.5b or defer.
- ORACLE `proto_p35_implicit_tangent.py` **4/4 PASS**: slip K_ss==FD (1.5e-9) + NON-symmetric
  (asym/kt=4.0); stick==kt·P_t symmetric; committed-N (drop d_TN)==symmetric; full assembly
  self-equilibrium (K·u_rigid≈0). 3D assembled tangent LOCKED.
- Wrote `_adr39_p35_design.md` + ran **DESIGN GATE (2 source-grounded reviewers): tangent-mechanics
  PASS + solver/active-set SALVAGEABLE.** Folded:
  - **GATE-Q2 BLOCKER (decisive):** non-symmetric default = silent-wrong-answer foot-gun — symmetric
    SOEs (ProfileSPDLinSOE.cpp:327, BandSPDLinSOE.cpp:269, SymSparse) read ONLY the upper triangle +
    silently drop the lower; the FE can't see the SOE type. **FLIPPED the default to the SYMMETRIC
    tangent** (drop the d_TN⊗n column ⇒ correct on ANY solver, superlinear); non-sym consistent
    tangent (quadratic) = OPT-IN `-consistanttan` + a parser WARNING re FullGeneral/UmfPack.
  - **GATE-Q3 MAJOR (deferred to P3.5b):** per-step active-set freeze; symmetric default removes the
    d_TN chatter driver; document `NewtonLineSearch`. **Q4 resolved:** residual ALWAYS current-N
    (correct equilibrium); symmetric default = approximate tangent for the exact residual (same root,
    +1 iter). **Q5 addKiToTang stick-only = PASS** (SPD contraction). Explicit byte-identity = PASS.
  - tangent-mechanics reviewer FD-validated the FULL ndof assembly to 1.5e-9; sign/scatter/c1 all exact.
- **NEXT = code P3.5:** symmetric friction tangent in `addKtToTang` (Gᵀ D_TT P_t G; stick kt·P_t /
  slip (μN kt/‖tT*‖)(P_t−n̂⊗n̂), current N, drop d_TN) + `-consistanttan` opt-in adds d_TN⊗n + parser
  warning; `addKiToTang` stick-only; tangent reads committed gpT NOT gpTtrial. → build → test
  (static stick converges = THE gate (singular in P3); static slip; superlinear default vs quadratic
  `-consistanttan` on FullGeneral; Newmark dynamic; explicit byte-identity; FD-at-slip tangent==∂resid)
  → code gate → PR base ladruno.

### Iteration 22 RESULT — P3.5 implicit friction tangent CODED + 6/6 green + code gate PASS
- C++: kernel `frictionTangentBlock` (3×3 K_ss; stick kt·P_t / slip (μN kt/‖tT*‖)(P_t−n̂⊗n̂);
  +d_TN⊗n iff consistent) + FE `addFrictionTang` (Gᵀ K_ss G scatter, w=[1,−N_i]) wired into
  `addKtToTang` (reads COMMITTED gpT) + `addKiToTang` (stick-only kt·GᵀP_tG); `consistentTan`
  member (default false=symmetric); Contact.consistentTan + addContact arg; handler threads it;
  parser `contact … -consistanttan` + the FullGeneral/UmfPack WARNING. No new file/classTag.
- Incremental build exit 0. **P3.5 6/6 green**: static stick converges (THE gate — singular in
  P3) + symmetric-solver-safe (ProfileSPD) + `-consistanttan` ≤-iters (FullGeneral) + static slip
  w/ anchor (friction=μN) + Newmark dynamic (a=(Q−μN)/m) + explicit byte-identity. **Full battery
  45/45** (6 P3.5 + 39 prior), no regression.
- **CODE GATE (1 source-grounded reviewer): PASS, NO findings.** Transcription cross-checks against
  the oracle EXACTLY (0.0) for stick/slip/G-scatter/symmetric-vs-nonsym/committed-N≡consistent=false;
  predicate identity (tangent stick switch == residual's); committed gpT read; addKiToTang stick-only
  SPD; explicit byte-identity (no friction-tangent path under CDL); all init/threading correct. Only
  non-bug: design-doc naming drift `-symtan`→`-consistanttan` (fixed w/ a superseding note).
- **NEXT**: commit P3.5 → PR base ladruno (wait Zone-A green before merge).

## Deferred / parked
- **P3.5b**: per-step active-set FREEZE / chatter detector (the design-gate Q3 MAJOR, deferred —
  symmetric default removes the d_TN chatter driver; NewtonLineSearch is the v1 mitigation).
- P2b-2c ∂n/∂u normal tangent + Hertz + SOFT Courant floor; P4 SOFT, P5 segment-based, P6 tied,
  AL upgrade (Q-AL), MPI — all per ADR.
