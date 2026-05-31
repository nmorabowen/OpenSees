---
title: Stack Profiler (performance + memory)
project: Ladruno
status: draft
priority: medium
owner: nmora
tags:
  - implementation
  - profiler
  - solver
  - analysis
---

# Stack Profiler (performance + memory)

## What

A consistent, low-overhead profiler that instruments the OpenSees solution stack to capture:

1. **Performance along the stack** — per-phase wall + CPU time for one analysis step, broken into the same taxonomy for **both implicit (Newton) and explicit (Linear/central-difference)** paths, plus per-element-assembly aggregates.
2. **Memory** — live/peak/cumulative allocation of the core numeric objects (`Matrix`/`Vector`/`ID`), a live-object census of domain components (leak detection across analysis boundaries), and an optional Windows CRT checkpoint bracket for true leak reports in debug builds.
3. **A machine-parseable output** (JSON) — an aggregate rollup tree (flame-graph friendly) plus an optional per-step time series (time-history friendly), emitted through the existing `OPS_Stream` file path and controlled by a new `profiler(...)` interpreter command.

**Not in scope (follow-up plan):** the visualization layer. The output schema is designed *for* visualization (hierarchical rollup + time series) but the viewer/parser is a separate deliverable.

## Why

OpenSees has **no usable profiler today**:

- The `Timer` class (`SRC/utility/Timer.h/.cpp`) is a **no-op stub on Windows** (`#ifdef _WIN32 // fill in later` — `getReal()`/`getCPU()` return `0.0`). On this Intel-oneAPI/Windows build there is effectively no wall-clock instrumentation.
- The only per-phase timing that exists (`getTotalTimeCPU`/`getSolveTimeCPU`/`getAccelTimeCPU` on `EquiSolnAlgo`, exposed as `totalCPU`/`solveCPU`/`accelCPU` Python commands) is implemented in **only two algorithms** — `ModifiedNewton` and `AcceleratedNewton` — and rides on the broken `Timer`. `NewtonRaphson`, `Linear`, and every explicit integrator return `0.0`.
- There is **zero** memory instrumentation anywhere in `SRC` (no allocation counters, no CRT hooks, no arena).

For the Ladruno explicit-dynamics work (CentralDifference/ExplicitBatheLNVD, energy balance) and large nonlinear implicit runs, we need a *consistent* way to answer "where is time going" and "is memory growing" across both solver families, with output we can parse and plot.

## Where

Architecture map (verified against source — file:line):

### Shared driver (implicit + explicit both flow through here)
- `SRC/analysis/analysis/DirectIntegrationAnalysis.cpp:168` `analyze()` → `:191` `analyzeStep()`: `newStep` (:213) → `solveCurrentStep` (:221) → `commit` (:248).
- `SRC/analysis/analysis/StaticAnalysis.cpp:127` — same shape, inline loop.

### Implicit inner loop (Newton)
- `SRC/analysis/algorithm/equiSolnAlgo/NewtonRaphson.cpp:127` `solveCurrentStep()` — the 6 per-iteration call sites are the cleanest, backend-agnostic instrumentation plane:
  - `formUnbalance` (:143 initial, :198 per-iter)
  - `formTangent` (:181)
  - `theSOE->solve()` linear solve (:187)
  - `update()` (:193) → `Newmark::update` → `theModel->updateDomain()` (element state recompute)
  - `theTest->test()` (:204)

### Explicit inner loop (Linear algorithm)
- `SRC/analysis/algorithm/equiSolnAlgo/Linear.cpp:91` `solveCurrentStep()` — `formTangent` (:109, lumped mass, `factorOnce`), `formUnbalance` (:119, **dominant cost**), `solve` (:125, trivial M⁻¹), `update` (:133, predictor/corrector).
- `ExplicitBatheLNVD` runs **two** residual+solve sub-steps per analyze step (`SRC/analysis/integrator/ExplicitBatheLNVD.cpp:424-426`, override `formUnbalance` at :464) — the profiler must count sub-steps, not assume one solve per step.

### Shared assembly (one edit covers all integrators)
- `SRC/analysis/integrator/IncrementalIntegrator.cpp` — `formTangent` (:85, static), `formUnbalance` (:130), `formElementResidual` (:267, element loop + `addB`), `formNodalUnbalance` (:246).
- `SRC/analysis/integrator/TransientIntegrator.cpp` — `formTangent` (:68, adds DOF_Group mass/damping loop), `formUnbalance` (:122). **All stock explicit integrators reuse these base bodies** — they do not override them (only `ExplicitBatheLNVD` does, and it delegates to the base first).
- Per-element seam: `SRC/analysis/fe_ele/FE_Element.cpp` `getTangent` (:287), `getResidual` (:311) → into element `getTangentStiff`/`getResistingForce(IncInertia)`. **Instrument as an aggregate around the loop, not per call** (element counts are large).

### Memory seams
- `SRC/matrix/Matrix.cpp:104` / `Vector.cpp:69` / `ID.cpp` — the `new (nothrow) double[]`/`int[]` sites every numeric buffer flows through (live allocator; `malloc`/`free` paths are commented out). **Best single seam for allocation counting** (~6 lines each).
- `SRC/tagged/TaggedObject.h:45` ctor/dtor — live-object census (Element/Node/Material via DomainComponent). Diff across analyses → leak detection.
- Caveat: `FE_Element` `theMatrices[]`/`theVectors[]` (`:103`), `DOF_Group` statics, and `Truss` `trussM*/trussV*` statics are **deliberate never-freed scratch singletons** — must be whitelisted so CRT leak dumps / census diffs don't false-positive.

### Control + output (mirror existing `start`/`stop` timer commands)
- Register `profiler`: `SRC/interpreter/PythonWrapper.cpp:3166` (next to `start`/`stop`), `SRC/interpreter/TclWrapper.cpp:1830`. Body `OPS_profiler()` in `SRC/interpreter/OpenSeesCommands.cpp:3196` (next to `OPS_startTimer`), declared in `OpenSeesCommands.h:440`.
- Active-analysis handles via the `cmds` singleton (`OpenSeesCommands.cpp:153`): `getStaticAnalysis()`/`getTransientAnalysis()`, `OPS_GetDomain()` (:1213), `getSimulationInformation()` (run header).
- Output via `OPS_Stream`/`FileStream` (`SRC/handler/`). No JSON stream exists today — the registry serializes JSON itself. Run-level header (start time, units, files) can be pulled from the live `SimulationInformation`.

### New code
- `SRC/utility/profiler/PerfClock.{h,cpp}` — portable wall+CPU clock (fixes the Win32 stub).
- `SRC/utility/profiler/Profiler.{h,cpp}` — registry (node tree), `ScopedTimer`, memory counters, JSON dump.
- `SRC/utility/profiler/ProfilerMacros.h` — `OPS_PROFILE_SCOPE(...)` etc., no-op when compiled out.

## How

### Three layers

**Layer 1 — Portable timing core.** A `PerfClock` using `std::chrono::steady_clock` for wall time and `GetProcessTimes` (Windows) / `getrusage` (POSIX) for CPU time — replaces the broken `Timer` Win32 branch. A process-global `Profiler` singleton holds a **tree of named `ProfileNode`s**; each node accumulates `{callCount, wallNs(total/min/max), cpuNs, allocBytesDelta}`. A thread-local "current node" pointer makes nesting automatic. RAII:

```cpp
{ OPS_PROFILE_SCOPE("formTangent");   // pushes child, starts clock
  ... }                               // dtor stops + accumulates
```

When the build flag is off, `OPS_PROFILE_SCOPE` expands to nothing → **zero overhead**.

**Layer 2 — Consistent phase taxonomy (same tree for implicit & explicit).**

```
step                       (one analyze increment)
├─ newStep
├─ solveCurrentStep
│  ├─ formTangent          implicit: K assembly │ explicit: lumped-mass build (once)
│  │  └─ elem.tangent      aggregate over FE_Element loop
│  ├─ linearSolve          implicit: factor+solve │ explicit: trivial M⁻¹
│  ├─ update               → updateDomain (element state recompute)
│  ├─ formUnbalance        residual; explicit's dominant cost
│  │  ├─ elem.residual     aggregate over FE_Element loop
│  │  └─ nodalUnbalance
│  └─ convTest             implicit only (Linear/explicit: absent)
└─ commit                  → commitDomain (material commitState)
```

Instrumentation strategy: **call-site scopes** in `NewtonRaphson::solveCurrentStep` and `Linear::solveCurrentStep` (backend-agnostic, covers solve/test which have no shared base body) **+ base-body scopes** in `IncrementalIntegrator`/`TransientIntegrator` `formTangent`/`formElementResidual` (one edit covers every integrator; aggregate element loops). `analyze`/`analyzeStep`/`commit` scopes in the analysis driver. This combination gives an identical phase tree for both solver families — the consistency is the whole point.

**Layer 3 — Memory.**
- Atomic counters in `Matrix`/`Vector`/`ID` ctor+dtor: `liveCount`, `liveBytes`, `peakBytes`, `cumulativeAllocs`. Captures transient churn (the per-element `tmp(numDOF)` vectors in `FE_Element`/`DOF_Group`) and steady growth.
- `TaggedObject` ctor/dtor census: live Element/Node/Material counts; `profiler` can snapshot and diff across analysis boundaries → "after teardown, N components still live" leak signal.
- Optional `_CrtMemCheckpoint`/`_CrtMemDifference` bracket (debug builds) for true allocation-site leak reports, with the deliberate static singletons whitelisted.

### Public API (Python; Tcl mirrors)

```python
profiler('start')                              # begin accumulating
profiler('stop')
profiler('reset')
profiler('config', '-perStep', True)           # also keep per-step time series
profiler('config', '-warmupSteps', 5)          # drop first K steps from rollup (P1#8)
profiler('report', 'runs.h5', '-run', 'caseA') # append this run as /runs/caseA in the HDF5 dataset
profiler('memory')                             # -> live/peak bytes snapshot (also queryable mid-run)
# comparison/diff is done offline in the Python backend, not in-engine:
#   ProfilerResults('runs.h5').diff('caseA', 'caseB')
```
`report` appends a new run group to the `.h5` (never mutates existing runs), so repeated analyses accumulate into one dataset.

### Output schema (logical view)

The schema below is the **logical** structure. At rest it is materialized as **HDF5** (see "Storage & visualization architecture"): each node = an HDF5 group, metrics = group attributes, the per-step `series` = chunked datasets, one run = one group under `/runs/<id>`. The JSON shown here is the **wire format** the Python backend serves to the React frontend, not the on-disk format.

```jsonc
{
  "meta": { "title": "...", "units": "...", "start": "...", "nSteps": 1200 },
  "rollup": {                       // hierarchical -> flame graph
    "name": "step", "calls": 1200, "wall_ms": 8421.0, "cpu_ms": 8390.0,
    "children": [
      { "name": "solveCurrentStep", "calls": 1200, "wall_ms": 8100.0, "children": [
        { "name": "formTangent",  "calls": 3600, "wall_ms": 3950.0,
          "children": [ { "name": "elem.tangent", "calls": 0, "wall_ms": 3700.0 } ] },
        { "name": "linearSolve",  "calls": 3600, "wall_ms": 2600.0 },
        { "name": "formUnbalance","calls": 7200, "wall_ms": 1200.0 },
        { "name": "convTest",     "calls": 7200, "wall_ms": 60.0 }
      ]}
    ]
  },
  "series": [                       // optional -> time-history plots
    { "step": 1, "t": 0.005, "dt": 0.005, "iters": 3,
      "wall_ms": { "formTangent": 3.2, "linearSolve": 2.1, "formUnbalance": 1.0 },
      "mem": { "live_bytes": 10485760, "peak_bytes": 12582912 } }
  ],
  "memory": { "matrix_live": 1422, "vector_live": 880, "peak_bytes": 12582912,
              "components_live": { "Element": 5000, "Node": 5151 } }
}
```

### Testing
- Smoke: a small cantilever, run 10 steps implicit + 10 explicit, `profiler('report')`, assert JSON parses and the phase tree is non-empty for both.
- Overhead: same model with profiler off (compile flag) vs on; assert <~2-3% phase-level overhead.
- Memory: run an analysis, tear it down, re-run; assert `components_live` returns to baseline (no leak) and `matrix_live` does not monotonically grow across steps.

## Decisions (locked 2026-05-30)

1. **Gating — single binary, all runtime-gated (revised 2026-05-30).** The whole profiler — coarse phase timing, deep per-element buckets, *and* memory counters — is compiled into the **one shipped binary** and toggled by **runtime** flags (`enabled` / `deep` / `mem`), so the user controls all three layers per run with no special build. `OPS_PROFILE_SCOPE` gates on `enabled()`; `OPS_PROFILE_SCOPE_DEEP` on `enabled() && deep()`; `OPS_PROFILE_COUNT_ALLOC/FREE` on `mem()`. Cost when off is one predicted-not-taken branch (the `mem()` branch is ~1 ns against a heap `new[]` ~100 ns — lost in the noise). *Why not compile-gate deep/memory:* "firing" (the increment) and "checking" (the branch) are different costs; a runtime branch is negligible, and a single binary removes the two-build / "which build do I have?" complexity the user explicitly wanted to avoid. The build switch flips to an **opt-out**: `-DOPS_PROFILE_DISABLE` compiles every `OPS_PROFILE_*` macro to `((void)0)` — literally zero, not even a branch — for a benchmark/purist build only. `profiler('start')` = coarse; `profiler('start','-deep')` adds per-element; `'-memory'` adds alloc counters.
2. **Storage — HDF5 (revised 2026-05-30).** The engine writes an HDF5 dataset (the engine already links HDF5 for `.mpco`), **one group per run** under `/runs/<id>`, the node tree as a group hierarchy (metrics = attributes), and the per-step `series` as chunked datasets. *Not* JSON-at-rest: JSON is only the wire format the Python backend serves to the frontend. HDF5 wins because (a) it handles the large per-step `series` (10⁵–10⁶ rows) that JSON bloats on, (b) the node tree maps 1:1 onto HDF5 groups so the stable path keys (P0#6) become literal HDF5 paths — diffing two runs = walk both group trees by path, (c) it keeps the toolchain consistent with `.mpco`/STKO_to_python.
3. **Visualization — Python backend + React frontend (decided 2026-05-30).** Three tiers: C++ → `profile.h5`; a **Python backend** (`h5py` + a small FastAPI/Flask app) reads HDF5, applies normalizers/diff, serves JSON; a **React** frontend renders flame/icicle (rollup), time-series (series), run-comparison/diff, and dt_cr/leak badges. The frontend never reads HDF5 directly — that is the backend's job. Lives outside the OpenSees source tree (e.g. `Ladruno_tools/profiler_viewer/`); only the HDF5 writer is in-engine.
4. **Memory v1 — in-process counters.** `Matrix`/`Vector`/`ID` allocation counters + `TaggedObject` live-component census. Cross-platform, cheap. Windows CRT bracket and Intel Inspector are explicitly **deferred** (debug-only add-on / external workflow).

### Storage & visualization architecture

```
C++ engine ──writes──► profile.h5  (HDF5 at rest; one group per run)
                          │
                  Python backend (h5py + FastAPI)  ──normalize / diff──► JSON over HTTP
                          │
                  React frontend  ──flame graph · time-series · run comparison · badges
```

HDF5 layout:
```
profile.h5
  /runs/<run_id>                         one group per run (the "dataset" of many runs)
    attrs: model, engine_sha, integrator, algorithm, solver, units, timestamp,
           dt_min, dt_max, nSteps, nDOF, nElem, nNode, nnz, threads,
           dt_cr, oversample_ratio, wall_ms_total, cpu_ms_total   ← normalizers (P0#2,#5)
    /rollup/...                          node tree as nested groups; attrs = {calls, wall_ns,
                                         wall_ns_min/max, cpu_ns}; path = stable key (P0#6)
       .../formTangent/elem_by_type      dataset [classTag, count, wall_ns, wall_ns_per_elem, fb_coupled]  (P0#1)
       .../linearSolve/{factor,triSolve} (P1#7)
    /series/                             chunked datasets length nSteps: step, t, dt, iters,
                                         wall_ms[nSteps×nPhase], mem_live_bytes, mem_peak_bytes
    /memory/                             attrs + components_live dataset [classTag, count]
```
Backend endpoints (sketch): `GET /runs` (manifest from group names+attrs), `GET /runs/{id}/rollup`, `/series`, `GET /diff?base=&cand=` (walk both rollups by HDF5 path → per-node Δ%/Δshare — this *is* P0#6).

**Crash-durability (heed [[project_mpco_exit_crash]]):** profiling runs are exactly the ones that may crash, and `.mpco` is known to be left unflushed on a kernel `exit(-1)`. v1 = buffer in memory, write the `.h5` at `profiler('report')`/`stop` (simplest); fast-follow = periodic chunked flush of `/series` + a `try/finally` flush hook so a mid-run crash keeps partial data.

### Overhead budget

- **Phase-level (always-on):** ~2 clock reads (`steady_clock`→`QueryPerformanceCounter`, ~20 ns) + thread-local push/pop ≈ 50–100 ns/scope; ~20 scopes/step → ~2 µs/step against a ms-scale step ⇒ **<0.5%, usually <0.1%**. Leave on by default.
- **Deep per-element (`OPS_PROFILE` + runtime toggle):** one scope per element per iteration; overhead is **inversely proportional to element cost** — ~5–10% on cheap trusses (~1 µs), **<1% on force-based fiber beams (10–100 µs)**, i.e. cheap exactly when you'd profile. Opt-in.
- **Compiled out:** literally zero (macros expand to nothing).
- **Real contention risk = the memory atomics** (global `std::atomic` in `Matrix.cpp:104`/`Vector.cpp:69` fire per scratch alloc under threaded assembly, inflating `elem.*`) → per-thread counters summed at report (P1#9).
- **Windows CPU granularity ~15 ms** (`GetProcessTimes`) → per-step `series` records **wall only**; CPU meaningful only at the run rollup (wall−cpu parallelism signal). Ship the self-measured `profiler.overhead` node (P2#14) so the budget is data, not a promise.

## Phased build

**Engine (C++):**
- **P1 — Timing core (no engine edits).** `PerfClock`, `Profiler` registry, `ScopedTimer`, macros, **HDF5 writer** (groups=nodes, attrs=metrics, chunked `series`). Self-test against a synthetic node tree. New CMake handling for the `OPS_PROFILE` flag + HDF5 link (touch [[../Ladruno_internal/01_compilation_journal]]).
- **P2 — Phase seams.** Scopes in `DirectIntegrationAnalysis`/`StaticAnalysis` (`analyze`/`analyzeStep`/`commit`), `NewtonRaphson::solveCurrentStep` (6 sites), `Linear::solveCurrentStep` (4 sites). Yields the full phase tree for implicit + explicit. **+P0#3 DONE (PR #52):** per-step series — `OPS_PROFILE_STEP(getCommitTag(), getCurrentTime(), dt, theAlgorithm->getNumIterations())` after the commit in both `DirectIntegrationAnalysis::analyzeStep` (dt = the integrator `dT`) and `StaticAnalysis::analyze` (dt = 0; pseudo-time = load factor). Armed by `profiler('start','-perStep')` (the command/`config().perStep`, writer `writeSeries`, and `Profiler::recordStep` scalar+mem path already existed since P1/P5 — this wires the missing call site). `getCommitTag()` is monotonic across the whole run (incremented in `Domain::commit()`). Gated on `enabled()` so the getters are only called during a run; macro short-circuits, `-DOPS_PROFILE_DISABLE`→`((void)0)`. Activates the frontend time-series (iters/step + dt/step). Smoke `Ladruno_tools/profiler_viewer/series_smoke_{model,check}.py` (transient SDOF, Newmark, `analyze(20,0.01)`): 20 rows, monotonic step 1..20, dt=0.01, Newton iters=2/step. **Per-phase `wall_ms[nSteps][nPhase]` series still deferred** (recordStep appends a zero-width row; needs per-step phase-wall capture — harder, lower value). **+P0#5** dt/dt_cr/oversample into the run header (still pending).
- **P3 — Assembly aggregates (deep). DONE (per-element-class), PR #49.** `OPS_PROFILE_SCOPE_DEEP_NAMED` named deep scopes `elem.tangent` (in `IncrementalIntegrator::formTangent` + `TransientIntegrator::formTangent`'s FE loop) and `elem.residual` (in `IncrementalIntegrator::formElementResidual`, reused by the transient path). Per element, `OPS_PROFILE_FE_ELEM_SCOPE(tmr, elePtr)` times the `getTangent`/`getResidual` call and folds the wall into the **+P0#1** per-element-class `elem_by_type` bucket keyed by `getElement()->getClassTag()`, with the `fb_coupled` tag from `isForceBasedClassTag()` (classTag whitelist of the force-based beam-column families; heuristic, kept off the vanilla Element API). New RAII `ElemScope` + `ScopedTimerOpt::engaged()` accessor: the macro short-circuits on `tmr.engaged()` so `getElement()` (a non-inline call) is touched ONLY when the deep gate is live — an unprofiled run pays nothing; `-DOPS_PROFILE_DISABLE` compiles it all to `((void)0)`. `count` is the cumulative element-evaluation count (N_elem × n_calls), so `wall_ns_per_elem` = per-evaluation kernel cost. Smoke `Ladruno_tools/profiler_viewer/deep_smoke_{model,check}.py` (20-truss static, `-deep`): `elem.tangent` Truss×400 @ 44.5 ns/elem, `elem.residual` Truss×600 @ 45.2 ns/elem, `fb_coupled` False. **Deferred:** **+P0#4** ExplicitBatheLNVD sub-step scope placement (`:424-425`) — explicit-integrator-specific, separate follow-up.
- **P4 — Memory counters. COMPLETE (byte counters + census).** Per-thread (not global-atomic, P1#9) counters in `Matrix`/`Vector`/`ID` ctor+dtor; `TaggedObject` census; static-scratch whitelist.
  - **Byte counters DONE (Matrix + Vector), PR #46; ID, PR #47.** `AllocType`-tagged `countAlloc/countFree` (`OPS_PROFILE_COUNT_ALLOC/FREE(bytes, type)`) wired at every per-object `new[]`/`delete[]` in `Matrix.cpp`/`Vector.cpp`/`ID.cpp` — alloc unconditional after `new[]`, free under the `fromFree==0` ownership guard before `delete[]`, so external-data views and moves stay balanced; the static `matrixWork`/`intWork` scratch is whitelisted (not counted). `ID` uses `arraySize` as the allocated capacity so its frees are exact; bare-`if` deletes were brace-converted and `ID::operator=` captures the old `arraySize` in a temp before its overwrite. `buildMemorySnapshot` lowers the per-thread `liveBytesByType[]` into the snapshot's `matrix_live`/`vector_live`/`id_live` + aggregate `peak_bytes`.
  - **TaggedObject census DONE, PR #48.** `components_live` `[{classTag,count}]` filled from a per-thread `std::map<int,int64_t> componentLive` bumped by `OPS_PROFILE_CENSUS_BORN/DIED(classTag)` (gated on `mem()`). **Seam is `MovableObject` (NOT `TaggedObject`):** `TaggedObject` stores only the object `theTag` and `getClassTag()` is virtual (unsafe in a base ctor/dtor); `MovableObject::classTag` is a plain member valid through the dtor, and Element/Node/Material all derive from it. `buildMemorySnapshot` merges the per-thread maps in ascending-classTag order and emits one `ComponentLive` row per net-positive count. Counts ALL MovableObjects by raw classTag — the integer space is shared across families (e.g. `NOD_TAG_Node == MAT_TAG_ElasticMaterial == 1`) so a row is a superset; the viewer filters by classTag band. The census reflects components **born while `mem()` is armed**, so arm `-memory` before building the domain for a full live census. Smoke `Ladruno_tools/profiler_viewer/mem_smoke_{model,check}.py` (20-truss static, armed before build): `matrix_live=160, vector_live=1456, id_live=892, peak_bytes=2844`; `components_live` = `{1:72, 2:4, 6:1, 12:20}` (classTag 12 = `ELE_TAG_Truss` × 20; classTag 1 = Node/ElasticMaterial superset), all balanced.
  - **Known approximation:** `Vector::resize` shrink-keeps its buffer while lowering `sz`, so a later free under-counts the shrunk delta (rare; `liveCount` stays exact).
- **P5 — Command + output.** `OPS_profiler()` body + Python/Tcl registration; `start`/`stop`/`reset`/`config`/`report`/`memory` subcommands; HDF5 append (one group per run); **+P0#2** size normalizers from the live SOE/census; pull run header from `SimulationInformation`; crash-flush hook. Fold/retire the legacy `getTotalTimeCPU` algorithm counters.
- **P6 — Engine tests.** Cantilever smoke (implicit + explicit), overhead budget check, leak/no-growth check, HDF5 round-trip.

**Tooling (outside the source tree, `Ladruno_tools/profiler_viewer/`):**
- **P7 — Python backend. DONE, PR #50.** `ProfilerResults` loader (`h5py`) + normalizers + `diff(base, cand)` (walk rollup groups by path) — the headless analysis API (Jupyter-usable) — **plus the FastAPI HTTP layer `profiler_api.py`**: `GET /health`, `/runs` (manifest), `/runs/{id}/header|rollup|series|memory`, `/diff?base=&cand=`, and the auto Swagger `/docs`. Owns no analysis logic — every endpoint forwards to a `ProfilerResults` method, so the wire format is exactly what the loader documents. Opens the HDF5 **per request** (open→read→close) so a run appended by a later `report` shows up without a restart and there is no long-lived file lock; permissive CORS (local dev tool); 404 on unknown run / missing series, 503 on a missing dataset. `create_app(path)` factory + `$PROFILER_H5`-driven module `app` (`uvicorn profiler_api:app`) + a `--file/--host/--port` CLI. `test_api.py` exercises every endpoint via Starlette `TestClient` (no live server); a live `uvicorn` boot was smoke-tested against `profile_sample.h5` (`/diff` reproduces the caseA→caseB formTangent −67%). Deps added to `requirements.txt` (`fastapi`, `uvicorn`).
- **P8 — React frontend. DONE, PR #51.** Vite + React + TypeScript viewer in `Ladruno_tools/profiler_viewer/frontend/`, a pure client of the P7 API (never reads HDF5). Panels: run picker (view + vs-baseline selection); **hand-rolled SVG icicle/flame** of the rollup (children laid proportional to the parent's `wall_ms` so self-time shows as a trailing gap; click a frame → `NodeDetails` with the `elem_by_type` per-element-class table + force-based/disp-based chip, P0#1); time-series (per-phase wall / iters / dt via a generic downsampling SVG `LineChart`); memory panel (per-type live bytes + peak + the TaggedObject census as sorted bars); run-to-run diff table (per-node Δ along the stable path, green=faster, P0#6); header badges (scheme / nDOF / nElem / nSteps / threads / `dt vs dt_cr` + oversample, with a "could raise dt" hint when oversample>1.5). Typed `api.ts` mirrors the wire format; permissive-CORS direct fetch to the backend (`VITE_API_BASE` override). No charting dependency. Verified: `npm run build` (tsc + vite) clean; live dev-server boot against `profile_sample.h5` exercised every panel through the backend (icicle drill-in shows classTag 74 force-based dominating; diff colors caseA→caseB; 3 series charts; census sorted). TS strictness honored (`verbatimModuleSyntax` → `import type`; `erasableSyntaxOnly` → string-union tabs, no enums).
- **P9 — Advisory + polish.** P2 deltas (advisory ruleset, exclusive/self time, `setup/numbering` node, accuracy co-signal).

## FEM-workflow panel review (2026-05-30)

A four-lens expert panel (element-kernel · linear-algebra · time-integration · methodology) stress-tested whether this design lets us *assess a FEM workflow AND propose/validate improvements*, not just time it.

**Verdict: partial-but-fixable.** The skeleton (phase taxonomy, call-site seams, implicit/explicit unification, assembly-vs-solve fork, explicit trivial-M⁻¹ negative result, leak yes/no) is correct and endorsed by all four lenses. But the deliberate "aggregate one scope around the whole element loop" decision discards the **element-type → material → integration-order** axis on which every formulation fix is made, and the absence of a run-to-run diff + size normalizers means improvements can be observed but not *proven*. Six P0 additions close the gap.

**Essential gaps (each breaks a specific, common prescription):**
- **G1 — no per-element-TYPE attribution.** A mixed model's loop average can point at the wrong element class (false confidence). Blocks "switch this element's formulation / integration order."
- **G2 — no run-to-run diff / stable node-path keys.** Absolute wall-clock with no diff contract → improvement claims unfalsifiable. Blocks *proving* a fix.
- **G3 — no size normalizers** (`theSOE->getNumEqn()`, element/node count, nnz/profile). "ms/step" non-portable, "elem.tangent ms" uninterpretable, no scaling study. Blocks ProfileSPD→sparse / renumber decisions.
- **G4 — `linearSolve` is one opaque timer.** No factor-vs-triangular-solve split (there is **no** `factor()` seam at `LinearSOE` — must cut inside each concrete `solve()`). Blocks full-vs-modified-Newton / factor-reuse quantification.
- **G5 — per-step iteration count not wired** (exists at `NewtonRaphson.cpp:205`, schema `series.iters` exists). Blocks the integrator-fix vs solver-fix fork (total solves = iters×steps).
- **G6 — no dt / dt_cr / oversample for explicit.** Engine already computes dt_cr (`damped_minimum_critical_timestep`, `ExplicitBatheLNVD.cpp:320`) but never surfaces it. Blocks the highest-value explicit lever (mass-scale / raise dt).
- **G7 — ExplicitBatheLNVD's 2nd residual mis-placed.** It executes *inside* `update()` (`ExplicitBatheLNVD.cpp:424-425`), so it lands as an opaque child of `update`, not a sibling `formUnbalance` — analyst mis-concludes the kinematics are expensive. A scope-*placement* fix, not just a count.

**Adjudicated tensions:** (a) per-type bucketing overhead objection does **not** survive the ratio — `getClassTag()`+lookup is ns vs a μs `∫BᵀDB`; do it (DEEP-gated). (b) The real contention risk is the **memory atomics** in `Matrix.cpp:104`/`Vector.cpp:69` (Heisenberg inflation of the very `elem.*` numbers) → move to per-thread counters summed at report time. (c) For force-based elements (`forceBeamColumn`) the tangent/residual split is an artifact of coupled state determination → tag buckets force-based vs displacement-based. (d) Timing in an accuracy vacuum is a trap for explicit → co-report a cheap energy-balance/unbalance signal so a "fast but unstable-dt" run isn't green-lit.

### Prioritized delta (merge into Decisions / Phased build / schema)

**P0 — must-have before this is a profiler you can write a PR from:**
1. **Per-element-class buckets** at the loop seam (`IncrementalIntegrator.cpp:267`, `TransientIntegrator::formTangent`), keyed by `getClassTag()`. Schema: `elem_by_type:[{classTag,count,wall_ms,wall_ms_per_elem,fb_coupled}]` under `elem.tangent`/`elem.residual`. *Amends P3.* (G1, force-based tag)
2. **Size normalizers** in `meta`: `nDOF=theSOE->getNumEqn()` (free), `nElem`/`nNode` (promote census), optional `profileSize`/`nnz` once per `domainChanged`; emit derived `wall_ms_per_step` / `_per_dof` / per-elem-per-call. *Amends P5.* (G3)
3. **Per-step iteration count, non-optional** — wire `numIterations` (`NewtonRaphson.cpp:205`) into the step node; count Bathe sub-steps. *Amends P5.* (G5)
4. **ExplicitBatheLNVD sub-step scope placement** — wrap the internal `formUnbalance()`/`solve()` at `ExplicitBatheLNVD.cpp:424-425` in the *same* named scopes (siblings, not children of `update`) + `subStep` counter. *Amends P3.* (G7)
5. **dt / dt_cr / oversample** into the always-on rollup header `{dt_min,dt_max,nSteps,dt_cr_estimate,oversample_ratio}`, routing the already-computed dt_cr. *Amends P5.* (G6)
6. **Run-to-run diff contract** — stable hierarchical path keys, deterministic child ordering, offline `profiler('diff', base, cand)` → per-node Δabs/Δ%/Δshare; declare `min` (or trimmed-mean) the comparison statistic. *New, amends Decisions + P5.* (G2)

**P1 — high-value:** 7. split `linearSolve`→`solve.factor`+`solve.triSolve` inside each concrete solver (G4). 8. steady-state windowing (`config -warmupSteps K`) to defeat the `FE_Element.cpp:123` first-touch outlier + turbo ramp. 9. per-thread memory counters (replace hot global atomics). 10. per-class material/state-determination sub-cost on the `update`/`commit` path. 11. normalized intensity (`wall_ns/numIP`, `/nDOF²`) for the integration-order decision. 12. thread count in `meta` + tag thread-summed nodes.

**P2 — later:** 13. `setup/numbering` node sibling to `step` (RCM + `CriticalTimeStep` CFL sweep at `ExplicitBatheLNVD.cpp:280`, which live *outside* `solveCurrentStep`). 14. exclusive/self time per node + measured `profiler.overhead` node + `wall−cpu` field. 15. compute-vs-scatter split (`getResidual` vs `addB` at `IncrementalIntegrator.cpp:277`). 16. stability/accuracy co-signal per step (unbalance norm, KE-ratio, energy-balance). 17. thin advisory ruleset over the rollup (never hides the raw tree).

*The three that change the tool's character: P0#1 (per-type → prescriptive on the kernel), P0#6 (diff → can prove a fix), P0#5 (dt_cr → prescriptive on explicit). The other P0s are near-free and remove specific misreads.*

## Risks / open questions

- Threading: explicit/parallel solvers and MKL threading — counters must be `std::atomic`; the node-tree current-pointer is thread-local (per-rank under MP). Cross-MPI-rank aggregation is a later concern.
- Reuse, don't duplicate: P5 must fold/replace the existing `getTotalTimeCPU` algorithm counters rather than leaving two parallel timing systems.
- Static-scratch false positives in any leak diff (FE_Element/DOF_Group/Truss statics) must be whitelisted.

## Implementation log

- **2026-05-30 — ✅ END-TO-END SMOKE TEST PASSED (built OpenSeesPy + ran a real profiled analysis).** After the Fortran fix, the full build reached the final link and failed on 3 BLAS/LAPACK unresolved symbols (`DSYGV` in `CriticalTimeStep.cpp`, `BLAS_DGEMV_X/2_X`) — because the configure ran without the **MKL** env (`Could NOT find BLAS/LAPACK`). Fix: activate the **full oneAPI env** (compiler **+ MKL** `mkl\latest\env\vars.bat`) so `find_package(MKL)` resolves (MKL 2025.1.0, static, intel_lp64) → **`OpenSeesPy.dll` linked, 35 MB**. Smoke test (Python 3.14, the build's linked Python): copied `OpenSeesPy.dll`→`opensees.pyd`, `os.add_dll_directory(oneAPI\compiler\2025.1\bin)` for `libiomp5md.dll` (Python 3.8+ ignores PATH for ext-module DLLs), built a 2-node truss + `Static`/`Newton`/`LoadControl`, ran `profiler('start'); analyze(10); profiler('stop'); profiler('report','smoke.h5','-run','smoke')`. Result: `analyze()=0`, ux=0.001, `smoke.h5` written. Read back by Python `ProfilerResults`: **measured phase tree exact** — `step`(10) ▸ newStep(10)/commit(10)/solveCurrentStep(10) ▸ formTangent(20)/linearSolve(20)/update(20)/convTest(20)/**formUnbalance(30 = 3×10, proving the initial+2/iter seam placement)**; nDOF=1. THE WHOLE PIPELINE WORKS: P1 core → P2 seams → P5 command → HDF5 → Python backend, on a real OpenSees run. Deferred (refinement): run-header integrator/algorithm name strings (empty), per-step series.
- **2026-05-30 — Diagnosed + fixed the `ifx` Fortran build crash (upstream CMake bug, NOT the profiler).** The full-build blocker (`ifx: error #10273 ... 0xC0000005` on `.f` files) was root-caused by bisection: `CMakeLists.txt:213` did `add_definitions(-DOPENSEES_VERSION="${GIT_VERSION}")`, which passes the quoted 40-char git-SHA define to **all** languages including Fortran. `ifx`'s `/fpp` **buffer-overruns** parsing that define from a **response file** (`CMAKE_NINJA_FORCE_RESPONSE_FILE=ON` + long worktree paths force the rsp path; isolated bisection flipped the crash to `0xC0000409` STATUS_STACK_BUFFER_OVERRUN — smoking gun). Proof: removing only `OPENSEES_VERSION` from the rsp compiled the file cleanly. **Fix:** `add_compile_definitions($<$<COMPILE_LANGUAGE:C,CXX>:OPENSEES_VERSION="${GIT_VERSION}">)` — scope to C/C++ (Fortran never uses it). Verified: `DoddRestrepo.f.obj` (was crashing) now builds (exit 0, 13110 B). NOTE: Ninja writes `.rsp` files at **build** time, not configure — must rebuild to regenerate. This is a general fork build fix (belongs in the CMake-gotchas ledger), surfaced by the profiler build but independent of it.
- **2026-05-30 — Full in-tree build verification (configure + compile).** Drove the real Conan/Intel toolchain against the worktree: `conan install` (deps cached) + `cmake` configure → **exit 0**, build files written, and the profiler STATUS lines confirmed in full context (`HDF5 writer enabled (HDF5 1.14.0)`, `present, runtime-gated`). Then compiled the three C++ targets carrying all profiler code — **`OPS_Utilities`** (P1: PerfClock/Profiler/ProfilerHDF5Writer), **`OPS_Analysis`** (P2: NewtonRaphson/Linear/DirectIntegrationAnalysis/StaticAnalysis), **`OPS_InterpPyCmds`** (P5: OpenSeesCommands) — `143/143`, **zero errors**; plus the two wrapper objects directly (`PythonWrapper.cpp` in OpenSeesPy, `TclWrapper.cpp` in OPS_InterpTcl). **Every P1+P2+P5 .cpp compiles in-tree with cl/MSVC + the real include set + HDF5 + the actual CMake wiring.** ⚠️ The full **OpenSeesPy link / live smoke test is blocked by a pre-existing Intel `ifx` Fortran crash** (`error #10273 ... terminated by 0xc0000005`) on unrelated `.f` files (DoddRestrepo.f, STEELDR.f, drain/*.f) — an oneAPI/`ifx` environment issue, NOT the profiler (all profiler code is C++; the build reached 458/2481 with C++ objects compiling fine before hitting Fortran). Worktree build dir retained at `build/Release` for future incremental builds. Findings for the user's build flow: `makeWIN.bat` is stale (MUMPS at `OpenSees/mumps-build` not `..\..\mumps\build`; existing cache was Debug-typed); current `ifx` cannot compile OpenSees Fortran.
- **2026-05-30 — P5 command landed (drivable end-to-end).** `OPS_profiler()` in `OpenSeesCommands.cpp` (+decl in `.h`) + Python (`Py_ops_profiler` shim + `addCommand("profiler", ...)`) and Tcl (`Tcl_ops_profiler` + `addCommand(interp,"profiler", ...)`) registration, mirroring the existing `start`/`stop` timer commands. Subcommands: `profiler('start' [,'-deep'] [,'-memory'] [,'-perStep'])`, `'stop'`, `'reset'`, `'report', <file> [,'-run', <id>]`, `'memory'` (→ peak bytes). Uses the lazy process-global `theProfiler()` — same instance the P2 hot-path scopes accumulate into (no `OpenSeesCommands` ctor/ownership change needed). `report` assembles `mergedRollup()`/`buildMeta()`/`series()`/`buildMemorySnapshot()`, wires the **nDOF normalizer** (`cmds->getSOE()->getNumEqn()`, P0#2), and appends one run group via `ProfilerHDF5Writer`. Verified: the profiler-side API call sequence type-checks against the real headers under g++ (`-I SRC/utility`); `<profiler/Profiler.h>` resolves in the interpreter target like `<Timer.h>` does. **End-to-end path now exists** (P1 core + P2 seams + P5 command): `ops.profiler('start'); ops.analyze(n); ops.profiler('stop'); ops.profiler('report','run.h5','-run','caseA')` → readable by the Python `ProfilerResults`. Full compile/link gate = next `makeWIN` build. (Other run-header fields — model/integrator/algorithm/solver names, nElem/nnz, dt_cr — left defaulted; refinement.)
- **2026-05-30 — P2 phase seams landed (implicit + explicit).** Instrumented 4 engine hot-path files with `OPS_PROFILE_SCOPE` (include `<profiler/ProfilerMacros.h>`, resolves via SRC/utility on the path — same as `<Timer.h>`, no CMake change): `NewtonRaphson::solveCurrentStep` (formUnbalance/formTangent/linearSolve/update/formUnbalance/convTest — 6), `Linear::solveCurrentStep` (formTangent/formUnbalance/linearSolve/update — 4), `DirectIntegrationAnalysis::analyzeStep` (step/newStep/solveCurrentStep/commit — 4), `StaticAnalysis::analyze` (step-per-iteration/newStep/solveCurrentStep/commit — 4). Yields the identical phase tree for both solver families. Critical invariant honored: **one scope per `{ }` block** (the RAII timer pops the thread-local cursor before the next sibling), so phases are siblings not nested; the outer `step` is a function-level (Direct) / per-iteration-block (Static) parent. Verified: diffs reviewed, brace balance intact, `<profiler/ProfilerMacros.h>` resolves + macros expand under g++ in both default and `-DOPS_PROFILE_DISABLE` modes. Edits were behavior-preserving (only added include + wrapping blocks; touched lines had trailing whitespace normalized — cosmetic). **Per-step series (`recordStep`: iters/dt — P0#3) deliberately deferred** to keep P2 to the rollup phase tree; the tree is fully usable once the P5 command toggles `enabled` + reports. Full compile/link gate = next `makeWIN` build (same as P1).
- **2026-05-30 — P1 in-tree build verification.** Staged the module into the main tree's configured Ninja/Conan/cl build and reconfigured: `add_subdirectory(profiler)` processed cleanly, the `HDF5_FOUND` guard passed (HDF5 1.14.0), and the profiler STATUS lines printed — **CMake wiring configures**. (The broader reconfigure then died on an *unrelated* `find_package(ZLIB)` "NOT USING CONAN" issue — a bare `cmake .` doesn't re-apply the Conan toolchain; reproducing the full `makeWIN.bat` configure, `conan install` + `-DCMAKE_TOOLCHAIN_FILE`, is needed for a whole-program build, not anything in the profiler.) Then compiled all three TUs to objects with the **production `cl`/MSVC** compiler (the writer against the real Conan HDF5 include, and `Profiler.cpp` under `-DOPS_PROFILE_DISABLE`) — all clean. Combined with the writer link+run (self-test) and Python parity, P1 is proven at configure + compile + link + parity; only the full whole-program link awaits a normal `makeWIN` build. Main tree restored clean (staging removed, `SRC/utility/CMakeLists.txt` reverted via git).
- **2026-05-30 — P1 landed + verified (agent sweep).** Engine module `SRC/utility/profiler/`: `PerfClock.{h,cpp}` (portable wall=steady_clock/QPC + cpu=GetProcessTimes/getrusage — replaces the Win32 no-op `Timer` stub), `Profiler.{h,cpp}` (ProfileNode tree, thread-local cursor + per-thread merge, `ScopedTimer` RAII, runtime gate), `ProfilerMacros.h` (single-binary runtime gating: `OPS_PROFILE_SCOPE`→`enabled()`, `_DEEP`→`enabled()&&deep()`, `COUNT_ALLOC`→`mem()`; opt-out `-DOPS_PROFILE_DISABLE` compiles all to nothing — revised from the earlier compile-gated `OPS_PROFILE`, verified g++ both modes), `ProfilerHDF5Writer.{h,cpp}` (sole `#include <hdf5.h>` TU; mirrors the Python `write_run` contract), `CMakeLists.txt` (core unconditional on `OPS_Utilities`, writer gated by `HDF5_FOUND`, `option(OPS_PROFILE OFF)`), `profiler_selftest.cpp`. **Verified end-to-end:** g++ syntax-checks the HDF5-free core in both gating modes; the writer compiled with MSVC + Conan static HDF5 (`cl /MT`), and the Python `ProfilerResults` reads the C++-written `.h5` identically to `make_sample` caseA (rollup, `linearSolve` split, `elem_by_type` fb_coupled flags, diff=0.0 on identical trees); `test_contract.py` still green. Sweep's adversarial review found defects, all fixed: **C1** (illegal chunk on the `nPhase==0` default-P1 series → contiguous fallback), the masked **empty-`phases` `H5Awrite(NULL)`** rejection in HDF5 1.14 (caught by an added nPhase==0 regression run), **H1** (cross-thread merge race → documented precondition + warn-if-enabled), **M1** (`ScopedTimer` ctor `noexcept` vs throwing `child()` → degrade-to-null), **M2** (explicit move-delete). Parity review = OK (field-by-field table). Build artifacts gitignored. **No engine hot-path edits yet — that is P2.**
- **2026-05-30 — Schema-first prototype (P7 head start, pre-engine).** Pinned the HDF5 contract + Python backend before any C++. `Ladruno_tools/profiler_viewer/`: `profiler_schema.py` (reference `write_run()` = the contract the C++ writer mirrors), `profiler_results.py` (`ProfilerResults`: manifest / normalized rollup / series / memory / `diff` by stable HDF5 path), `make_sample.py` (3-run fixture), `test_contract.py` (passing: round-trip, run immutability, normalizers, `elem_by_type`, diff). Validated end-to-end: full-Newton→modified-Newton diff shows formTangent −67% / solve.factor −67% / step −48.7%; per-class breakdown shows 200 force-based beams (15.8 ms/elem) dominate 20k trusses (9.9 µs/elem) — i.e. G1+G2+P0#1/#6 demonstrated in Python. README lists the C++ writer's obligations. Engine phases P1–P6 unchanged; the on-disk layout is now locked.
