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
| P1 skeleton+handler+brute-force, zero-force | NOT STARTED | **adversarial design gate FIRST** | — | the BLOCKER-risk wiring |
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

### Iteration 2 — P1 design + adversarial design gate (NEXT)
- Write `_adr39_p1_design.md` (FE_Element adapter + custom ConstraintHandler
  injection + conservative-static connectivity + data structures).
- Launch the P1 adversarial DESIGN gate (Workflow) BEFORE writing C++ — mirrors
  ADR-30's pre-code Gate-A that caught 10 seam fixes. This is THE critical junction
  (the B1/B2 BLOCKER-risk wiring).

## Deferred / parked
- P4 SOFT, P5 segment-based, P6 tied, AL upgrade (Q-AL), MPI — all post-v1 per ADR.
