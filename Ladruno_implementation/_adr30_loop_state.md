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
| P1 core (equalDOF) | BUILT, 12/12 green; Gate-B running | guppi/adr30-p1-core | T1,T2,T4,T5 + P0 all pass | Gate-A done; Gate-B + full-ZoneA running |
| P2 general C (rigidLink/diaphragm) | not started | — | T3,T6,T7 | Gate-C |
| P3 queries + EQ | not started | — | T8,T9 | Gate-D |

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
