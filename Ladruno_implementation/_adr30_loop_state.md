# ADR-30 implementation loop — progress ledger

**This file is the resumable source of truth for the self-paced `/loop` building
`LadrunoProjectionHandler` (ADR 30). Read it FIRST every iteration; update it LAST.**
Context gets summarized across a multi-day run — trust this file, not memory.

Spec: `30_ladruno_explicit_constraint_projection_adr.md`. Host integrator:
`SRC/analysis/integrator/CentralDifferenceLadruno.cpp` (classTag 33003).

## Loop policy (locked with user 2026-06-19)
- **PR cadence:** approval gate per phase. Loop drives a phase autonomously
  (implement → build → test → adversarial panel → fix → bookkeeping) then STOPS
  and surfaces the PR for the user to review/merge. Do NOT auto-merge.
- **Scope:** full v1 = P0 → P1 → P2 → P3. P4 stays deferred per ADR §7.
- **Adversarial gates:** multi-lens `Workflow` panel per gate (skeptics prompted to
  REFUTE + verify-the-finding follow-ups). CONFIRMED correctness bug ⇒ BLOCK
  (loop back to implement). Refuted findings logged + dismissed (no thrash).
- Each phase = fresh branch off `ladruno`. Before any follow-up push check
  `gh pr view <n> --json state` == OPEN (stranded-commit trap).

## Build / test recipe
- Build (P1+ only; P0 needs no SRC change): `cmd /c build.bat OpenSeesPy` via the
  PowerShell tool. `classTags.h` + analysis-core edits ⇒ near-full recompile (slow).
  Seed `mumps-install` + `mumps-archive` from the main clone first if missing.
- Test interpreter: `C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe`
  (numpy 2.4.6, pytest 9.0.3). Bootstrap: `os.add_dll_directory(DIST)` +
  `PYTHONPATH=DIST`, `DIST = <clone>\dist\bin`.
- P0 runs against the existing built pyd in the MAIN clone
  (`C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin`, build 8c8edb2) — it
  only exercises already-shipped pieces (CDL 33003 + Transformation + Diagonal/BandGen).

## Seams (verified against code 2026-06-19)
- New fork files → `SRC/analysis/handler/`: `LadrunoProjectionHandler.{h,cpp}`,
  `LadrunoConstraintProjector.h`, `LadrunoConstraintGraph.{h,cpp}`,
  `LadrunoProjectionConsumer.h` (abstract 1-method interface).
- CDL hook (fork-owned, NO vanilla): starter projection between
  `*Aprev = theLinSOE->getX()` (CentralDifferenceLadruno.cpp:446) and the
  `v_{-1/2}` seed (:447); main-path projection in `update()` before `Vfull` (~:507);
  IC check + mass-cache invalidation in `domainChanged()`.
- Vanilla touches (LEDGER_vanilla rows, `// Ladruno` comment at each edit):
  `classTags.h` HANDLER_TAG 33001; `OPS_ConstraintHandler()` branch at
  `SRC/interpreter/OpenSeesCommands.cpp:1604`; Tcl `SRC/tcl/commands.cpp`;
  `SRC/interpreter/PythonAnalysisBuilder.cpp`; `FEM_ObjectBroker` handler entry.

## Phase status
| Phase | State | PR | Tests | Gate |
|---|---|---|---|---|
| P0 falsify & baseline (no SRC) | ✅ MERGED | #300 (671b5236) | 3 passed (hardened) | Gate-0 SOUND |
| P1 core (equalDOF) | ✅ MERGED | #305 (668ad8e9) | 11 P1 + 3 P0; ZoneA 870 | Gate-A + Gate-B SOUND |
| P2 general C (rigidLink/diaphragm) | ✅ MERGED | #307 (9990e206) | 6 P2 + general-review fixes | Gate-C + general review SOUND |
| P3 queries + EQ | ✅ MERGED | #312 (2760cf37) + #313 handoff | T8/T9/EQ/guard/missing-node; ZoneA 938 | Gate-D SOUND + EQ-leak upstream bug fixed |
| P4a tie-force recorder | ✅ MERGED | #317 (61341837) | p4 records+h5py readback ±F·m₂/M; ZoneA 942 | focused + general review BOTH SAFE-TO-MERGE |

## Adversarial findings log
_(append per gate: finding · lens · REAL/REFUTED · resolution)_

### Gate-0 (P0)
- 2026-06-19: `tests/test_adr30_projection_p0.py` 3 passed against built pyd (8c8edb2).
  - **T6 CONFIRMED**: CD+Transformation+Diagonal grossly disagrees with both
    CD+Transformation+FullGeneral and implicit Newmark+Full on the coupling-induced
    uy; Diagonal uy suppressed <10% of correct. Diagonal silently drops the
    condensed `T^T M T` (uy,rz) off-diagonal → wrong physics. ADR D2 premise holds.
  - **Massless (§2.5) observed behavior**: zero-mass free DOF under CD+Diagonal →
    "starter: SOE solve failed" → analyze returns -2 → state stays 0 (NOT silent
    NaN garbage). Still not the correct finite response + cryptic message →
    confirms handler must give a NAMED handle()-time error. (Refine wording in P1.)
  - **Baseline**: made informational (timing noise-dominated at n=400; the real
    win is O(n) vs O(n^2) dense memory at 50k DOF — infeasible to allocate dense).
  - **Gate-0 panel (wf_853e5055-d32, 11 agents): verdict SOUND-WITH-FIXES on all 4
    lenses.** Core T6 falsification CONFIRMED SOUND — multiple agents independently
    hand-built the condensed 3x3 system + closed-form modal solution, matched the
    FullGeneral runs to <1%; D2 premise holds, not refuted. 6 MAJOR confirmed-blocking
    findings, all in the massless test (not T6). RESOLUTIONS applied:
    - massless test fully rewritten → `test_massless_dof_is_not_policeable_by_the_soe_layer`:
      asserts on analyze RETURN CODE (not the old permissive (0,0)); positive-control
      (mass=1) isolates the cause; demonstrates the EMPIRICAL finding — Diagonal aborts
      opaquely (rc=-2 aii=0) while FullGeneral/BandGeneral swallow the singular factor
      and return rc=0 with garbage (SILENTLY WRONG). This STRENGTHENS §2.5 → handle()-time
      scan is mandatory; ADR §2.5 amended with the finding.
    - T6 hardened: added handler-independent analytic eigen-anchor (`_analytic_uy`,
      OpenSees-free modal solution, asserted <2e-2); fixed false "rz tracks closely"
      comment (Diagonal detunes rz, now asserted); added rigidLink-vs-rigidDiaphragm
      stand-in note (ADR §6 amended).
    - baseline narrative corrected (timing is environment-noisy: came out 0.8x AND the
      panel saw 5-17x → genuinely non-deterministic; real arg is O(n) vs O(n^2) memory).
    - LEDGER_quirks: T6 + massless-SOE quirk row appended.
  - Minor findings (non-blocking, noted): IC kick is constraint-non-compliant (shared
    by all 3 runs, doesn't bias the comparison; analytic anchor sidesteps it); the two
    OpenSees refs share the Transformation handler (now anchored by the OpenSees-free leg).
  - **P0 verdict: COMPLETE, ready for user approval. Doc/test-only, no SRC change.**

## Running notes / decisions
- 2026-06-19: loop started. P0 test = `tests/test_adr30_projection_p0.py`.

### Gate-A (P1 pre-code math, wf_09671e08, 15 agents)
- **Verdict: algebra core SOUND** (C1 M-orthogonal/idempotent, C2 momentum, C5 equalDOF→mass-avg
  all rigorously verified; operator correctly general-C → P2 not blocked). All blockers in seam/hook.
- 10 confirmed fixes folded into `_adr30_p1_design.md` "Gate-A resolutions" (R1–R10). Load-bearing:
  - R1 IC check uses MP offset form `(u_c−Uc0)=C(u_r−Ur0)` (verified accessors exist MP_Constraint.h:82-83).
  - R2 mass source = Diagonal-SOE diagonal ONLY (Node::getMass omits element mass); read BETWEEN
    formTangent and solve (verified DiagonalDirectSolver.cpp:124 overwrites A→1/aii on factor).
  - R3 `*Aprev=Aproj` at CDL line 533 is the load-bearing manifold write.
  - R4 SP-fixed retained → zero-overwrite orphaned slaves (not empty projector).
  - R5 consistent-mass guard (refuse cMass on tied DOFs).
  - R6 IC check covers velocity, every domainChanged, null-projector→warn.
- Code facts verified: DiagonalDirectSolver in-place reciprocal; MP Uc0/Ur0 accessors; DiagonalSOE getA.
- **NEXT: implement handler+projector+CDL hook, build, run T1/T2/T4/T5, then Gate-B (code panel).**

### P1 build + test (2026-06-20)
- Full OpenSeesPy build GREEN (worktree dist/bin; mumps seeded from main clone).
- P1 battery 9/9 + P0 3/3 = 12 passed. Dynamics regression 43 passed/1 xfail.
- 2 BUGS caught by tests before Gate-B (both fixed + rebuilt):
  1. doneNumberingDOF override skipped base FE_Element::setID() → added `this->ConstraintHandler::doneNumberingDOF()`.
  2. LATENT UPSTREAM: DirectIntegrationAnalysis::domainChanged() ignored handle()/doneNumberingDOF()
     return codes → handler conflict diagnostics printed but analysis ran on. Fixed driver to honor
     `<0`=error contract (vanilla edit, ledgered). T5 conflict battery caught it. (analyze checks
     domainChanged()<0 at :209 → propagates.)
- Vanilla edits this phase: classTags.h, DiagonalSOE.h (getDiagonalA), OpenSeesCommands.cpp,
  tcl/commands.cpp, FEM_ObjectBrokerAllClasses.cpp, analysis/handler/CMakeLists.txt,
  analysis/analysis/DirectIntegrationAnalysis.cpp. All ledgered (LEDGER_vanilla).
- Full Zone-A regression: **870 passed, 1 xfail, 0 regressions** (382s) — the universal
  DirectIntegrationAnalysis driver edit + CDL hook are safe across the whole suite.
- Gate-B code panel: wf_aa96881d-4e1 (running).

### Gate-B (P1 code panel, wf_aa96881d-4e1, 10 agents)
- **Verdict: projector math + implementation SOUND** (verified EXACT: tie_err=0, momentum 4e-16
  on unequal-mass multi-slave star; memory/ownership clean; mass-read-before-factor timing correct;
  no false aborts of other handlers; consistentMassGuard no false positives). 2 BLOCKERS, both the
  same error-contract family as the handle()/doneNumberingDOF fix:
  1. **IC-abort swallowed** — `DirectIntegrationAnalysis::domainChanged()` line 461 left
     `theIntegrator->domainChanged()` UNCHECKED, so CDL's enforceIC -1 (the -projectICs gate) was
     dropped → off-manifold ICs ran on silently. FIX: check `<0 → return -1` (mirrors handle path).
     [A Gate-B verify agent applied this edit to the source itself — audited via git diff, correct +
     clean, kept. PROCESS NOTE: Gate-B used agentType general-purpose (Edit access); future gates
     should use read-only Explore to avoid agents mutating the tree.]
  2. **Zero-group rebuild UAF** — doneNumberingDOF freed the old projector but skipped the consumer
     push under `numGroups()>0`, leaving the integrator with a dangling pointer (MPs removed mid-run).
     FIX: re-bind the consumer UNCONDITIONALLY — push projector if groups, push 0 if zero (CDL guards
     all derefs on `theProjector != 0`); clearAll comment corrected.
- Tests added: T3 (IC-violation aborts; -projectICs snaps+runs), zero-group-rebuild UAF guard
  (uses remove('mp',2)). Rebuild blj5yulya running.
- Minor/nit (deferred, non-blocking): near-singular probe only catches exact zero pivot; project()
  solve-failure `continue` is quiet; dead delta-size code; parallel send/recv scope (v1 single-proc).

### P2 (general-C transport) — 2026-06-20
- **KEY: P1 was already general-C** (L=[I;C], handler builds L from full Ccr) → transport needs
  NO new projector code. P2 = validation + boundary-pinning.
- Probes proved it: 2D rigidLink-beam under LadrunoProjection+Diagonal MATCHES Transformation+
  FullGeneral (rel≈0, <5e-7) — FIXES the P0/T6 mass-drop; transport gap u_c-(C u_r)=8.6e-16.
  3D rigidDiaphragm matches dense ref + stays rigid.
- §2.5 boundary CONCRETE: every TIED DOF keeps its own equation → needs nonzero lumped mass.
  rigidLink-beam / 3D rigidDiaphragm tie the slave PERPENDICULAR ROTATION (rz) → slave rz needs
  rotational mass (Transformation eliminates it, so it was free there). Massless rotational tie →
  refused (buildMass per-row m<=0). SP-fixing a diaphragm-controlled DOF → SP-on-slave refused.
  Frozen small-rotation Ccr = SAME limit Transformation has (not a regression). P4 = SOE-coop
  elimination to relax massless-tied-slave.
- Tests tests/test_adr30_projection_p2.py 5/5 (T6fix, T7, T3, massless-rot-tie refused, SP-on-tied
  refused) — pass on the SHIPPED P1 build (no rebuild).
- Bookkeeping: banner_features.txt + patch_banner (LadrunoProjection line, P2 = first user-visible);
  ADR §7 P2 marked DONE; impl ledger row updated.
- Gate-C: wf_dda26141-35a — READ-ONLY Explore agents (no tree mutation, fixing the Gate-B process gap).
- OPEN: near-singular LtML rank check only catches EXACT zero pivot (Gate-B/Gate-C minor) — decide if
  P2 hardens it (needs rebuild) or defers. NOTE: P2 added no SRC yet, so banner/ADR rebuild pending
  only if banner must show; tests don't need it.

### Gate-C (P2 general-C, wf_dda26141-35a, 13 READ-ONLY Explore agents)
- **Verdict: general-C transport mathematically SOUND** (3 of 4 lenses + verify stage confirm by
  code-read + probes; L built faithfully from multi-column Ccr, no sign/index error; momentum +
  manifold to machine eps; FIXES P0/T6 — uy RMS 0.020 vs diagonal-wrong 0.0).
- The lone "critical SP-fixed-master truncation" claim (transport lens) was REFUTED: dropping a
  HOMOGENEOUS-SP master column is CORRECT (that DOF's accel ≡ 0, so the term is genuinely zero);
  the edge-cases lens independently verified partial-column survival is right.
- 4 "confirmed_blocking" entries are NOT correctness bugs (per the verify stage's own fixes):
  2 = validation-strategy notes ("PASS / not a bug"), 2 = doc/UX gaps. RESOLUTIONS:
  - Added **T6b**: OpenSees-INDEPENDENT analytic anchor (forced modal response incl. slave rz mass)
    → proj ≈ analytic, closes the triangle (addresses "reference not independent").
  - Improved the massless-tie error message (rotational-tie guidance: add ~0.01-0.1% rot mass).
  - LEDGER_quirks: tied-DOF-needs-mass + frozen-Ccr small-rotation (SAME limit as Transformation,
    not a regression) + near-singular LtML (ADR O1, deferred).
- DEFERRED (non-blocking, documented): runtime frozen-Ccr staleness guard (P2b/P4); condition-number
  rank check (ADR O1). Both are the SAME limits Transformation carries.
- Rebuild bkom5fq9z (projector message + banner regen). Then run P2 6/6 + regression → PR.

### General pre-merge review (P2, wf_0980ef20-3b6, 11 READ-ONLY Explore agents)
- **Verdict: SAFE-TO-MERGE** — projection pipeline fundamentally correct, no bugs in handler/
  projector/hook; house-rules compliant (ledgers, header stamps, classTag 33001 no collision,
  banner); vanilla DirectIntegrationAnalysis edit verified backward-compatible (870 ZoneA, 0 regr).
- 2 REAL consistency gaps found (same error-contract family as the main-path fix, missed elsewhere)
  — both FIXED in this branch:
  1. DirectIntegrationAnalysis::setIntegrator()/setAlgorithm() ignored domainChanged() returns
     (mid-session swap path) → added the `<0`⇒return -1 guard to both setters.
  2. TransientDomainDecompositionAnalysis::domainChanged() (parallel driver) dropped the
     doneNumberingDOF() return (overwritten by setSize) → added the check. Defense-in-depth
     (v1 projection handler is partition-interior only).
  Both ledgered (LEDGER_vanilla). The other 2 "confirmed_blocking" were verify-stage CONFIRMATIONS
  of correct claims (not bugs). Minor/deferred: condition-number rank check (ADR O1), frozen-Ccr
  staleness guard (P4) — both documented, both = same limits Transformation carries.
- Rebuild bi1sves9v; then test + transient regression → push to PR #307 → merge.

### P3 (queries + EQ) — 2026-06-20, branch guppi/adr30-p3-queries
- Tie-force query (option 1 + recorder): projector caches f_tie = M(a_raw-a_proj) per eqn in
  project(); handler getTieForce(node,dof) maps node/dof→eqn; command `ladrunoProjectionTieForce
  node dof` wired (decl OpenSeesCommands.h, impl OpenSeesOutputCommands.cpp via OPS_GetHandler+
  dynamic_cast, Py wrapper+register PythonWrapper.cpp, Tcl wrapper+register TclWrapper.cpp).
- EQ_Constraint support: buildGroups now also iterates theDomain->getEQs() (single cDOF tied to a
  coeff VECTOR of retained DOFs = multi-master general-C row the projector already handles). Two
  loops mirrored (vertex/edge build + slave re-walk). delta = Uc0 - sum c_k Ur0_k.
- -verbose: per-group retained/constrained/fixed summary at doneNumberingDOF.
- Tests tests/test_adr30_projection_p3.py: T8 (tie force = ±F*m2/M, equal/opposite), T9 (energy
  closure, drift<1e-3), EQ (equationConstraint enforced + tie-force query on EQ), GUARD (query
  refuses non-projection handler). Build bt3av58si running.
- RECORDER finding: LadrunoRecorder is node-based (iterates domain nodes) and CANNOT reach the
  handler-owned projector → native recording needs either the reaction slot (option 2, user leaned
  away) or a new nodal buffer + recorder source. DECISION: ship the lean query in P3; surface the
  recorder-source as a scoped follow-up at the P3 gate (don't force an awkward coupling).
- Vanilla touches added (ledger debt to record): OpenSeesCommands.h, OpenSeesOutputCommands.cpp,
  PythonWrapper.cpp, TclWrapper.cpp (the query command). EQ support is in the fork-owned handler.
- NEXT: build → run P3 tests → ledger/bookkeeping → Gate-D → PR.

### Gate-D (P3 code panel, wf_becac856-3c3, 5 READ-ONLY Explore agents)
- **Verdict: 3 of 4 lenses PASS (no bugs)** — tie-force cache (sign/timing/sizing correct, a_raw
  captured before overwrite, fixedEqn intentionally uncached), command plumbing (dynamic_cast null-
  safe, dof 1→0-based, registered Py+Tcl), EQ support (mirrored loops, delta=Uc0-Σc_k Ur0_k, chain/
  double diagnostics fire for EQ + mixed groups), tests rigorous, vanilla fully ledgered.
- 1 BLOCKER (confirmed): ADR §4 partition-crossing refusal was UNIMPLEMENTED — a cross-partition /
  missing node tag would be silently treated as SP-fixed (mis-assembly). FIXED: added a node-
  existence guard in buildGroups (every vtx node must be in this domain, else named refusal). Test
  test_constraint_to_missing_node_refused added. Rebuild b8x5icdc3.
- ~20 minor findings all CONFIRM correctness. No other action.
- **v1 (P0-P3) functionally COMPLETE.** Deferred to P4 (documented): native LadrunoRecorder tie-force
  source, prescribed-motion overwrite, MP-chain composition, ExplicitBathe/LNVD adoption, near-singular
  condition-number rank check, frozen-Ccr runtime staleness guard.

### P3 post-Gate-D: regression caught a REAL upstream bug (EQ leak)
- Gate-D partition guard added → full-suite run FAILED (9) but each test PASSED ALONE: classic
  test-ordering leaked-global-state signature. Root cause: `Domain::clearAll()` (called by `wipe`)
  never cleared `theEQs` — EQ_Constraint was a later upstream add, omitted from clearAll. The EQ test
  left a stale equationConstraint; the NEXT test's P3 handler (now reads getEQs()) picked it up →
  mis-assembled. INVISIBLE pre-P3 (no handler read EQs). FIX: `theEQs->clearAll();` in Domain.cpp
  (one line, upstreamable; ledgered LEDGER_vanilla + LEDGER_quirks).
- After fix: previously-failing order p3,p0,p1,p2 = 25 passed; dynamics regr 30+1xfail; full Zone-A
  running (b9zs7381a) — Domain::clearAll is universal so full regression is the gate before PR.
- **v1 (P0-P3) COMPLETE pending full-ZoneA + PR.** Vanilla P3 touches: OpenSeesCommands.h,
  OpenSeesOutputCommands.cpp, PythonWrapper.cpp, TclWrapper.cpp (query cmd) + Domain.cpp (EQ fix);
  fork handler/projector/CDL. All ledgered.

---
## ✅ v1 COMPLETE (2026-06-20)
P0 #300 · P1 #305 · P2 #307 · P3 #312 — ALL MERGED to ladruno. LadrunoProjectionHandler
(HANDLER 33001) ships explicit constraint projection for equalDOF + rigidLink + rigidDiaphragm
+ equationConstraint, with the tie-force query, IC manifold guard, full diagnostic battery, and
the deferred-P4 backlog documented. Handoff: `projection_handler_handoff.md`. Loop ends.

### P4a — native tie-force recorder (2026-06-20, branch guppi/adr30-p4-recorder)
- User chose "dedicated nodal buffer + source" (not reactions). Implemented:
  - Node.{h,cpp} (VANILLA): lazily-alloc projTieForce slot (in-class `=0` init → safe across all
    ctors w/o touching init lists; freed in ~Node), get/setProjectionTieForce. Ledgered.
  - CentralDifferenceLadruno::commit(): scatter projector tieForce → nodes (via DOF_Group getID)
    before commitDomain (recorders read after commit); only when projector active+massReady.
  - Recorder (fork): Ladruno_Types.h enum ConstraintTieForce; Ladruno_NodeResults
    ConstraintTieForceSource (vec3 force, TFx/y/z); LadrunoRecorder keyword `constraintTieForce`/
    `tieForce` + switch case. Writes CONSTRAINT_TIE_FORCE/DATA to .ladruno HDF5.
- Test tests/test_adr30_projection_p4.py: records + h5py readback DATA==±F*m2/M equal/opposite.
  Test bug found+fixed (h5py walk matched sibling STEP metadata → read /DATA only).
- ADR-30 battery+recorder 38 passed. Full Zone-A (buw5d7w6e) + focused P4 review (wf_ce3743b4) running.
- Correct by construction: recorder reads Node tie force = integrator scatter = projector cache =
  same source as the P3 query (validated vs analytics in T8).

### P4a general review (wf_454fdba8, 4 read-only lenses) + MERGED #317 (61341837)
- **Verdict: all 4 lenses SAFE-TO-MERGE; 0 serious, 0 confirmed-blocking; 19 findings all nits**
  (mostly correctness confirmations): Node core no blast radius (reactions/sensitivity/parallel/DB/
  copy unaffected, not serialized, no leak/double-free); scatter guard airtight (plain CDL + SMS pay
  nothing); recorder robust across envelope/f32/region/partition + enum backward-compatible; test
  binds (catches sign-flip + zero-scatter). No fixes needed.
- Minor follow-ups (optional, non-blocking): add a plain-CDL-no-projection zeros recording test;
  h5py-in-CI nit (test importorskips → could silently skip on a box without h5py).
- **P4a MERGED → ladruno (squash 61341837).** Native tie-force recording (`recorder ladruno ...
  constraintTieForce`) ships. ADR-30: v1 (P0-P3) + P4a recorder all on ladruno.
- Remaining P4 (deferred, independent): prescribed-motion overwrite, MP-chain composition,
  ExplicitBathe/LNVD adoption, near-singular LtML condition gate, frozen-Ccr runtime guard.

### P4b — prescribed-motion overwrite (2026-06-20, branch guppi/adr30-p4b-prescribed-motion)
- User picked this as the next phase. Design doc `_adr30_p4b_design.md`. Splits into
  **Tier 1** (prescribed DOF NOT in a constraint group — common multi-support/base-excitation;
  handler-only) and **Tier 2** (prescribed DOF IS an MP master — the literal "overwrite a
  before projecting", needs projector known-RHS surgery).
- **PRE-CODE design gate wf_5afb0c9a-8cb (5 refute lenses, read-only Explore, + verify):**
  Tier-1 mechanism SOUND on every lens (no approach refutation; the 22 "confirmed-blocking"
  are all "not implemented yet" = the to-do list). Confirmed: implement handler applyLoad()
  = `node->setTrialDisp(sp->getValue()+sp->getInitialValue(), dof)` per recorded non-homog SP,
  called from AnalysisModel::updateDomain→handler->applyLoad (runs AFTER domain applyLoad, so
  ImposedMotionSP's node vel/accel survive); keep prescribed DOFs at eqn=-1; integrator
  UNCHANGED. Critical safety check: a prescribed DOF used as MP/EQ **master** is silently
  dropped today (eqn<0 → fixedEqns → slave accel forced 0) → MUST refuse with named error.
  **Tier 2 DEFERRED to P4c** (gate confirmed Q1 real: starter project(a0) runs BEFORE that
  step's applyLoad, so the node-read of a_p fails at the starter; needs read-from-SP/GM-at-t0
  + free/prescribed column split + C_p·a_p RHS). Refuted: buildMass-wrong + IC-check-misfire
  (prescribed DOFs are eqn<0, not in the u/v vectors checkIC reads).
- **SCOPE P4b = Tier 1 + prescribed-master refusal.** Plan: handle() records non-homog SPs +
  refuses prescribed slave (overconstraint) and prescribed master (Tier-2 deferred); applyLoad()
  override sets prescribed disp; clearAll() clears records. Tests test_adr30_projection_p4b.py:
  base-excitation vs Transformation+BandGen, constant non-homog SP offset, prescribed-master
  refusal, prescribed-slave refusal.
- **IMPLEMENTED (handler-only, ZERO vanilla touch).** LadrunoProjectionHandler.{h,cpp}:
  handle() records non-homog SPs as PrescribedDOF{Node*,dof,SP*} + prescribedKey set (was
  warn-and-drop); applyLoad() override = setTrialDisp(getValue()+getInitialValue(),dof) per
  prescribed DOF; doneNumberingDOF refuses a prescribed MASTER (else-if prescribedKey when
  e<0, before the silent fixedEqns drop); SP-on-slave message now distinguishes PRESCRIBED
  SLAVE; clearAll() clears the records. Integrator UNCHANGED.
- **BUILD green** (full first build of the worktree, exit 0). **Tests: P4b 4/4** (TP1 base-
  excitation rel<1e-6 vs Transformation+FullGeneral + response nontrivial, TP2 settlement,
  TP3 prescribed-slave refused, TP4 prescribed-master refused — both named errors fire).
  **ADR-30 battery 30/30** (P0-P4 + P4b, no projection regression). **Full Zone-A: 949
  passed, 1 xfailed, 2 FAILED** — both pre-existing g++-standalone-kernel-compile env
  failures (test_hrz_standalone_kernel AssertionError; test_cpp_kernel_matches_oracle_dump
  WinError-2 CreateProcess) UNRELATED to P4b (diff = LadrunoProjectionHandler.{cpp,h} +88/-10
  only; CI Ubuntu runs g++ properly). Ledgers updated (impl row + quirks; no vanilla row).
- **Code-review gate wf_34b51896-209 (4 read-only lenses + verify, 17 agents) DONE.** 7
  "confirmed" findings triaged; 6 REJECTED on the mechanics, 1 ADOPTED:
  - REJECTED "homogeneous-fix drift" (blocking): an eqn=−1 DOF is never integrated (not in the
    eqn vector), setResponse skips it → stays at committed 0; unchanged from shipped P0-P4a +
    all Zone-A fixed-DOF models. The verify agent rubber-stamped a wrong premise.
  - REJECTED "MP-added-after-setup stale slave gate" + "missing sendSelf of prescribed records"
    + "Tier-2 IC" + "pre-analyze query shows un-imposed disp" (the proposed doneNumberingDOF
    fix is itself buggy — a load-pattern SP's valueC isn't set yet → would impose 0).
  - ADOPTED (findings 3+5, robustness/convention): PrescribedDOF stored a raw Node*; the rest
    of the handler stores node TAGS (vtxNode) and resolves via getNode. Switched PrescribedDOF
    to {int nodeTag; int dof; SP_Constraint* sp} + resolve `theDomain->getNode(nodeTag)` fresh
    in applyLoad (skip if null) — removes the theoretical dangling-Node* UAF, matches convention.
  - Rebuild bof4yzu7j (incremental). Re-run P4b + ADR-30 battery, then PR.
