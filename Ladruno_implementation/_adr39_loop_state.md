# ADR 39 ContactDomain — implementation loop state

> Persistent driver state for the phased ContactDomain build. Survives context
> resets. One source of truth for "where are we, what's next." Mirrors the
> ADR-30 `_adr30_loop_state.md` pattern. ADR = `39_ladruno_contact_domain_adr.md`.

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
| P1a FE+handler+empty-conn zero, bitwise | DESIGN DONE → CODING | design gate ✓ (SALVAGEABLE) | — | gate folded; replicate-handler + empty-conn + Domain commit/revert hooks |
| P1b ContactDomain+surface+lifecycle hooks | NOT STARTED | — | — | after P1a green |
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

## Deferred / parked
- P4 SOFT, P5 segment-based, P6 tied, AL upgrade (Q-AL), MPI — all post-v1 per ADR.
