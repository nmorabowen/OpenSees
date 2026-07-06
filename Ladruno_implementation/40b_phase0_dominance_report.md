---
title: Phase-0 dominance report — measured per-lane performance breakdown
project: Ladruno
status: measured (first pass, 2026-07-06)
priority: high
owner: nmora
amends: 40_ladruno_performance_adr
tags:
  - performance
  - measurement
  - phase-0
  - profiler
  - benchmark
---

# ADR-40b — Phase-0 dominance report (first measured pass)

> **The Phase-0 deliverable ADR-40 gates everything on.** Five lanes measured 2026-07-06 with
> the SHIPPED profiler (`profiler start -deep -perStep` → `report` → h5) on the existing
> `dist\bin\opensees.pyd` (2026-06-18) — **zero rebuild**. Five parallel measurement agents,
> one per lane; models + parser preserved in `Ladruno_files/testbed/perf/phase0/`.
> Run recipe: `phase0/PHASE0_RECIPE.md` (python 3.12 `-S`, manual pyd wiring, threads pinned).
> Caveats: single machine, single run per lane (no median-of-7 yet), moderate model sizes
> (0.7k–11.5k DOF), `cpu_ms_total` is a Win32 stub (all splits are relative wall_ns).
> **These numbers gate; they are not baselines.** The rank-3 bench extension turns them into
> locked baselines later.

## Headline: dominance is per-lane — three different winners

| Lane | Model (script) | Dominant bucket | Element state-det¹ | linearSolve | Attribution closure |
|---|---|---|---|---|---|
| **A fiber frame** | 510 forceBeamColumn, 136-fiber RC sections, 810 DOF, pushover 93 steps (`laneA_fiber_frame.py`) | **`update` 64.5%** = force-based element interior iteration | 1.8% *attributed* (the real cost hides in `update`) | 2.9% | step-level clean; `update` opaque |
| **B 3D solid J2** | 3375 LadrunoBrick + LadrunoJ2, 11.5k DOF, 15 steps, 2.9 it/step (`laneB_model.py`) | **UmfPack linearSolve 66.4%** | 30.1% (tangent 26.2 / residual 3.9) | 66.4% | 99.3% |
| **C plate-fiber shell** | 1600 ASDShellQ4 + LayeredShell(8× LadrunoJ2 PlateFiber), 9.8k DOF, 14 steps (`laneC_shell_model.py`) | **element 55.9%** — but it's shell/section machinery, not J2 | 55.9% (tangent 38.9 / residual 17.0) | 30.7% | 99.6% |
| **D explicit CDL** | 1000 LadrunoBrick -lumped + J2, 3.6k DOF, 2500 steps @0.8dt_cr (`laneD_model.py`) | **element 48.9%** (residual 33.3 / tangent 15.6) | 48.9% | ~0% (Diagonal) | **99.9%** |
| **E IMK frame** | 420 LadrunoIMKBeam2d (Bilin hinges), 660 DOF, cyclic 614 steps, 18× ductility (`laneE_imk_frame.py`) | **`update` 35.4%** (the hinge Newton) + solve 30.6% | 6.9% *attributed* (real ≈40%+ incl. update share) | 30.6% | `update`+`newStep` opaque ≈25% |

¹ *element state determination as attributed by the existing `elem.tangent`/`elem.residual` per-classTag scopes.*

Per-element costs worth keeping: LadrunoBrick implicit tangent **40.6 µs/ele** vs residual 5.4 (7.5×);
explicit residual 7.76 vs tangent 3.45; ASDShellQ4 tangent 47.8 vs residual 19.6 (2.4×);
forceBeamColumn/IMK getTangent/getResisting ≈0.1–0.5 µs/ele (trivially cheap — their cost is in `update`).

## Finding 1 — the #1 instrumentation gap is `Element::update()`, not the gaps ADR-40 listed

On the two frame lanes the dominant cost is **invisible to per-classTag attribution**: the
force-based interior iteration (lane A: 64.5%) and the IMK hinge Newton (lane E: bulk of 35.4%)
both run inside `element->update()` under the integrator's domain-update, which the profiler books
as one opaque `update` leaf with **no `elem_by_type` bucket**. A classTag bench reads lane E's
element cost as 2.8% when it is really a top-two consumer (~40%). This is exactly the ADR-40a §3.2
pathology, now measured. **Phase-0 scope work, priority order (revised):**
1. `OPS_PROFILE_FE_ELEM_SCOPE` on the domain-update element loop (`elem.update` bucket) — unlocks
   attribution on BOTH frame lanes (the fork's primary lane).
2. factor-vs-triangular-solve split in the SOE layer (gates rank 8/10; lane B makes it urgent).
3. ModifiedNewton scopes — still zero, still needed for the rank-5/8 reuse question.
- **REFUTED:** ADR-40's "the explicit CentralDifferenceLadruno path is unscoped" — lane D closes
  to 99.9% with full deep buckets. No explicit-path scope work is needed.

## Finding 2 — per-lane gate verdicts (ADR-40 ranks / ADR-68 targets)

| Item | Verdict | Measured basis |
|---|---|---|
| **rank 2** UMFPACK STRATEGY_AUTO + `-strategy`/`-pivotTol` | **PROMOTED — do first.** | Lane B is 66.4% linearSolve on the exact UmfPack path rank 2 touches; already regression-scoped, no gate needed |
| **rank 8/10** factorization reuse | Elevated interest, **still blocked on the factor-vs-solve scope** | Lane B 66.4% — but full Newton refactorizes by design; the lever pays on ModifiedNewton, and we can't split factor from solve yet |
| **rank 7** OpenMP element loop (explicit v1) | **Gate PASSED** (>40% ele) | Lane D: 48.9% element; residual-led. Blocked on the de-static audit (T3 + FE_Element/DOF_Group pools) |
| **ADR-68 T1** J2 `tang_flag` lazy tangent | **NOT authorized** — no lane shows a J2-tangent-bound profile | Lane B: J2-bearing residual is 3.9% of step; lane D's tangent waste is integrator-level (see Finding 3), not material |
| **ADR-68 T2** plane-stress closed-form return | **DEMOTED — ceiling ≈2% of wall** | Lane C elastic-baseline diff: condensation ≈ +2.64 µs/ele on residual ≈ 2.2% of step; iteration-count inflation (2→6.6 forms/step), not the material return, owns the plasticity cost. The G-EQUIV oracle ceremony is not worth ~1–2% |
| **ADR-68 T3** brick geometry cache (+ de-static) | **Conditionally authorized** — needs a micro-drill | Brick tangent is the element cost on lanes B (26.2%) and D; unknown how much of 40.6 µs/ele is geometry recompute vs B·D·Bᵀ. De-static proceeds regardless (thread-safety, C16 carve-out) |
| **ADR-68 T6** IMK hinge Newton | Lane is real (top-two consumer) — **blocked on the `elem.update` scope** | Lane E: cost is structurally unmeasurable per-classTag today; instrument, re-measure, then warm-start |
| **ADR-68 T4/T5** EAS inner Newton / LogStrain tax | Unmeasured this pass | Lanes used `std` formulation + small strain; the Lemaitre `profiler_dump.md` (bbar-vs-eas) is prior art for T4's drill-down |
| Lane A optimization (forceBeamColumn interior) | **Out of fork budget for now** — vanilla element, vanilla file | Instrument first (`elem.update` bucket); any optimization there is upstream-facing work per the change-budget policy |

## Finding 3 — two NEW explicit-integrator items (feed to ADR-67)

Lane D surfaced two costs no ADR lists, together ≈40% of explicit step time:
1. **CDL forms element tangents every step — 15.6% of wall** (2.5M `getTangent` evals) — while the
   Diagonal solve consumes no stiffness. If this is the constant lumped-mass path being reassembled
   per step, caching it is a **G-BYTE** skip worth ~15% of explicit wall. Needs a source drill-down
   (what does CDL's formTangent actually rebuild per step, and is it provably constant under
   fixed-M assumptions + no domainChanged?).
2. **`newStep` is 23.9%** — heavier than the leap-frog update kinematics (21.9%). Drill down what
   the CDL starter/mass path does per step before naming a fix.
Both are integrator-side (ADR-67's axis) and were only findable because the explicit path turned
out to be fully scoped.

## What Phase 0 says the program should do next (priority order)

1. **rank 2** (UMFPACK AUTO) — promoted, regression-free, lane B says it's the biggest single lever measured.
2. **`elem.update` profiler scope** — the one instrument that unblocks the fork's frame lanes (A, E) + T6.
3. **Lane-D drill-down** (CDL per-step tangent rebuild + newStep) — potential ~15%+ G-BYTE explicit win.
4. **T3 micro-drill** (geometry share of the brick tangent µs) → then the cache, then rank 7 v1.
5. **factor-vs-solve SOE scope** + ModifiedNewton scopes — before any rank-8 work.
6. **Drop/park:** T1 (no authorizing lane), T2 (ceiling ~2%), T4/T5 (await their benches).

## Addendum — `elem.update` scope built + frame lanes re-measured (2026-07-06, same day)

The Finding-1 instrument now exists: `Domain::update(void)`'s element loop carries a deep
`elem.update` per-classTag bucket (`Domain.cpp`, ledgered; raw-`Element*` macro variant
`OPS_PROFILE_ELEMPTR_SCOPE` in `ProfilerMacros.h`; inert without `profiler start -deep`).
Lanes A and E re-run on the instrumented build:

- **Lane A (fiber frame) — hypothesis CONFIRMED and quantified.** The formerly opaque `update`
  phase is now 92% attributed: `ForceBeamColumn2d` (classTag 73) `update()` costs **57.6 µs/ele**
  (~200× its `getTangent` at 0.28 µs/ele), and summing its three update sites
  (solveCurrentStep/update 12.9 s + newStep 3.1 s + DisplacementControl direct 2.0 s) books
  **~82% of step time to classTag-73 state determination**. The fiber lane is force-based-
  interior-bound, per the classTag record. Optimization of that interior loop is **vanilla,
  upstream-facing work** (change-budget policy) — but it is now measurable, so any future
  upstream contribution carries a number.
- **Lane E (IMK) — hypothesis REFUTED; ADR-68 T6 DEMOTED.** The hinge Newton is **cheap**:
  `LadrunoIMKBeam2d` `elem.update` = **0.62 µs/ele**, only ~14% of the `update` phase
  (486 of 3428 ms); total classTag-33004 work across ALL buckets ≈ 12% of step. The 35.4%
  "update" band is **integrator/domain update machinery OUTSIDE the element loop** (~1.8 ms per
  update call at 660 DOF — DOF_Group/node trial-state propagation, model-update glue), which the
  original pass mis-read as hidden element cost. **T6 (hinge warm-start) payoff ceiling < 1% —
  dropped.** The real lane-E cost structure: linearSolve ~30%, update machinery ~30%,
  newStep/step glue ~25% — a small-model *per-step fixed overhead* profile, not a kernel profile.
  New drill-down candidate: what the integrator update spends 1.8 ms/call on at 660 DOF.

Lesson worth keeping: the first pass correctly located the opaque band but **guessed its owner
wrong on lane E** — attribute, then rank; never optimize an unattributed band.

## Addendum — MUMPS parallel lane (measured 2026-07-06, same day)

Strong scaling of the implicit parallel path (`openseesmp`, hand-partitioned z-slabs, shared
interface planes): fixed 5488-element LadrunoBrick + LadrunoJ2 block, 18 900 DOF, 8 Newton steps,
`system Mumps` (auto ordering ICNTL7=7, matType 0, **uniform `-ICNTL14 2000`** — see the quirks
ledger), `numberer ParallelRCM`. 16-logical-core box. Script:
`Ladruno_files/testbed/perf/phase0/mumps_scaling.py`.

| np | wall (s) | speedup | efficiency | rank-0 linearSolve share | rank-0 formTangent share |
|----|----------|---------|------------|--------------------------|--------------------------|
| 1  | 65.3 | 1.00× | 100% | 57% | 35% |
| 2  | 38.1 | 1.71× | 86%  | 63% | 31% |
| 4  | 28.2 | 2.32× | 58%  | 74% | 21% |
| 8  | 33.0 | 1.98× | 25%  | 89% | 9%  |

Tip displacement agrees across np to ~1e-15 relative (partition correct; MUMPS ordering roundoff
only); all steps converged at every np.

**Verdicts:**
- **The implicit parallel path is solve-bound and gets MORE so with rank count** (57%→89%):
  element assembly strong-scales near-ideally (23.1 s → 2.9 s ≈ 8×), the distributed
  factorization does not — its wall *grows* past np=4 (20.7 → 29.2 s) as communication/fill
  overhead overtakes the shrinking per-rank arithmetic on a 19k-DOF matrix.
- **np=4 is the saturation point at this problem size** (2.32×, 58%); np=8 regresses. Bigger
  matrices will push the knee outward, but the shape confirms ADR-40's framing: parallel speed
  on the implicit path is a *solver* question, not an element-loop question — consistent with
  the serial lane-B result and with rank 11's demand-gating (no ParMETIS work is justified by a
  solver-bound profile at this scale).
- ADR-30 relevance: numbering/partition quality affects the factorization the profiler bills as
  `linearSolve`; any parallel-numberer claim must be re-measured against this baseline, at a
  production problem size, before ranking.

## Implementation log
- **2026-07-06 — first measured pass.** 5 lanes, 5 parallel agents, shipped binary, recipe +
  scripts preserved under `Ladruno_files/testbed/perf/phase0/`. h5 artifacts in session scratch
  (regenerate via the scripts — deterministic models). Numbers above are single-run wall_ns splits;
  re-run median-of-7 via the rank-3 bench extension before locking baselines.
