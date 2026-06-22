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

## Deferred / parked
- P4 SOFT, P5 segment-based, P6 tied, AL upgrade (Q-AL), MPI — all post-v1 per ADR.
