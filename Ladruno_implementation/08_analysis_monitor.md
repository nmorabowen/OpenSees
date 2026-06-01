---
title: Analysis Monitor (live streaming response + profiler telemetry)
project: Ladruno
status: draft
priority: medium
owner: nmora
tags:
  - implementation
  - recorder
  - profiler
  - tooling
---

# Analysis Monitor (live streaming)

## What

A **live monitoring device**: pick points / regions / response quantities before a
run, call `analyze(N)` once, and watch them plotted **live as the analysis steps**,
without waiting for the run to finish. Two feeds share one mechanism:

1. **Analysis monitor** (the star) — response at *selected* points/regions:
   nodal disp/vel/accel/reaction per DOF, element force/stress/strain, drift, base
   shear. Scalar **time-histories** only in v1.
2. **Profiler telemetry** — the per-step series the profiler already computes
   (phase wall/CPU, memory live/peak, iters, dt, residual norm), streamed live
   instead of only flushed at `profiler report`.

**Decided scope (from scoping session 2026-05-31):**
- **Shape B** — the engine pushes from inside the C++ step loop, so a single
  uninterrupted `analyze(N)` streams live (no Python `analyze(1)` babysitting loop).
- **Local-first, network outside C++** — the engine only ever does local SWMR-HDF5
  append. Watching a remote/headless/MPI run from a laptop (P1.5) is done by running
  the FastAPI sidecar *next to* the engine and crossing the network with the SSE/WS
  hop only — no socket code in the engine, ever.
- **Scalar time-histories only.** Spatial / deformed-shape / contour live view is
  explicitly a **P2 follow-up** (per-node field every frame is much heavier).

**Not in scope (follow-up):** spatial/deformed live view (P2); network transport
(P1.5); writing the picker selection *back* into a running analysis (selection is
fixed at recorder creation in Shape B).

## Why

The profiler answers "where did time go" only **after** a run (`profiler report`
→ HDF5 → FastAPI → React, all post-mortem). For long nonlinear implicit runs and
especially explicit dynamics (10⁴–10⁶ steps), you want to *watch the physics and the
cost as they happen* — catch a diverging drift, a runaway residual, a memory leak, or
a wrong base-shear trend at step 2,000 instead of discovering it an hour later. The
"pick points/regions" framing is also the **backpressure control**: monitoring only
chosen channels is what makes per-step streaming tractable at explicit step rates.

## Where

Architecture (Shape B, reusing the existing profiler viewer stack):

```
  C++ engine                    sidecar (existing FastAPI)        existing React
 ┌───────────┐  append-only    ┌────────────────────────┐  SSE/  ┌──────────────┐
 │ Monitor   │ ───────────────▶│ tails live sink,        │ ─WS──▶ │ live plots,  │
 │ recorder  │  SWMR-HDF5  or  │ serves /live/.../stream │        │ channel      │
 │ + sink    │  NDJSON tail    │ (decimate, backpressure)│        │ picker       │
 └───────────┘                 └────────────────────────┘        └──────────────┘
```

### Reuse map — almost everything already exists
- **Recorder selection parsing** — `NodeRecorder`/`ElementRecorder`
  (`SRC/recorder/NodeRecorder.h`, `ElementRecorder.h`) already parse node tags, DOFs,
  element tags, response strings, and `-region`. The "pick points/regions" picker is a
  solved parsing problem; the monitor recorder reuses it.
- **Per-step push seam** — `Recorder::record(commitTag, timeStamp)`
  (`SRC/recorder/Recorder.h:50`) is already called at every `commit()` by the domain.
  This is the emission point; no driver edit needed for the response feed.
- **Profiler per-step seam** — `OPS_PROFILE_STEP` already fires in
  `SRC/analysis/analysis/StaticAnalysis.cpp:231` and
  `DirectIntegrationAnalysis.cpp:269`. The telemetry feed is a tee off this existing
  call site into the same live sink.
- **Original live socket precedent** — `SRC/handler/TCP_Stream.{h,cpp}` (OpenSees's old
  MATLAB/Java GUI seam: an `OPS_Stream` that opens a socket and writes `Vector`s).
  Reference for the *network* transport (P1.5), NOT v1 — raw blob protocol, no framing,
  blocks on the reader (bad backpressure).
- **Viewer** — `Ladruno_tools/profiler_viewer/profiler_api.py` (FastAPI, per-request
  HDF5 reader) + `frontend/` (Vite/React, hand-rolled SVG time-series, no charting dep).
  Live mode = another data source + an SSE subscription, not a new app.

### New code
- `SRC/recorder/LadrunoMonitorRecorder.{h,cpp}` — recorder that reuses Node/Element
  selection, and on each `record()` pushes selected scalars to a live sink.
- `SRC/recorder/LadrunoLiveSink.{h,cpp}` — append-only **local SWMR-HDF5** sink
  (chunked, unlimited-`maxshape` frame dataset; `H5Fstart_swmr_write` after layout is
  defined). Single impl — the engine never speaks network (see resolved decision #3),
  so the "transport seam" is just the file-path contract between engine and sidecar.
- Tee from `OPS_PROFILE_STEP` (in `Profiler`) into the same sink when a live run is armed.
- Tooling: `profiler_api.py` gains `/live/{run}/stream` (SSE); `frontend/` gains a live
  mode + channel picker. Pure `Ladruno_tools/` — no SRC/ change.

### Modify (minimal, `// Ladruno` hooks)
- `SRC/utility/profiler/Profiler.{h,cpp}` — optional live-sink handle + tee at the
  existing per-step record point. No new call sites.
- Command registration for `monitor ...` (mirror `recorder`/`profiler` wiring in
  `SRC/interpreter/OpenSeesCommands.cpp`).
- Class tag: **allocate from `LEDGER_implementations.md`** (Ladruno band ≥33000) to
  avoid collision — do NOT hardcode before reserving the row.

## How

### Public API (sketch)
```python
# Shape B: arm once, then a single uninterrupted analyze(N) streams live.
ops.monitor('start', '-sink', 'live.h5',          # SWMR-HDF5
            '-every', 10,                          # step decimation: 1 frame / 10 steps
            '-hz', 30,                             # wall-clock throttle: <=30 frames/s
            '-node', 5, 7, '-dof', 1, 2, '-resp', 'disp',
            '-ele', 101, 102, '-resp', 'force',
            '-region', 3, '-resp', 'drift',
            '-profiler')                           # also tee profiler telemetry
ops.analyze(100000)                                # streams the whole time; never blocks on a reader
ops.monitor('stop')
```
Then `python Ladruno_tools/profiler_viewer/profiler_api.py` (already running) +
the React viewer in **Live** mode tails `live.h5` and plots the picked channels.

### Data flow
1. `monitor start` builds a `LadrunoMonitorRecorder`, adds it to the domain (so
   `record()` fires each commit), and opens a `LadrunoLiveSink`.
2. Each `record()`: gather selected scalars (reusing Node/Element response code),
   apply `-every` decimation, append one frame `{step, time, [values...]}` to the sink.
   **Engine never waits on a reader** — SWMR/NDJSON append is fire-and-forget, so a
   slow/absent viewer can never stall `analyze()` (the key failure mode to avoid).
3. Profiler telemetry tees into the same sink at the existing `OPS_PROFILE_STEP`.
4. FastAPI tails the sink (HDF5-SWMR read needs `HDF5_USE_FILE_LOCKING=FALSE` — see
   `LEDGER_quirks.md`), applies wall-clock throttle (~30 Hz cap), pushes SSE frames.
5. React subscribes via `EventSource`; existing SVG time-series panels append points.

### One file, two lifetimes
The live sink **is** a valid post-run profiler/results HDF5 — when the run ends it's
readable by the *existing* post-run API and viewer unchanged. No second format, no
second viewer; "live" is just "tail instead of read-all."

### Three knobs that are mandatory, not optional
1. **Selection** — only picked points/regions (the core ask; bounds the data).
2. **Decimation** — `-every K` and/or a wall-clock throttle in the sidecar.
3. **Backpressure** — viewer can't keep up → drop frames in the sidecar, never in the
   engine. The append-only/file-tail design enforces this structurally.

### Testing
- Smoke: a 2-DOF transient, `monitor start -node ... -sink live.h5`, `analyze(100)`,
  then h5py-read `live.h5` (`HDF5_USE_FILE_LOCKING=FALSE`) and assert ~`100/every`
  frames with monotonic step + plausible disp.
- Overhead microbench (acceptance gate): monitored vs unmonitored identical transient,
  added wall time within the profiler's <0.5% phase-level budget.
- Parity: monitored channel values must equal a parallel `NodeRecorder`/`ElementRecorder`
  on the same selection (the monitor must not *change* results, only observe them).
- Backpressure: kill the reader mid-run → `analyze()` completes at full speed, sink
  still well-formed.
- Frontend: verify via Claude_Preview MCP (`preview_start` + `preview_eval` DOM asserts;
  screenshots flaky — use eval), per the profiler-viewer P8 pattern.

## Resolved decisions (scoping session 2026-05-31)

> [!check] Sink format: **SWMR-HDF5**.
> One file, both lifetimes — the live sink is byte-valid as a post-run results file,
> so the existing viewer/API opens it unchanged when the run ends. No NDJSON, no
> converter. Cost: HDF5 SWMR (`H5Fstart_swmr_write` after defining datasets; chunked,
> unlimited-`maxshape` frame dataset) + the `HDF5_USE_FILE_LOCKING=FALSE` reader
> gotcha (already in `LEDGER_quirks.md`). SWMR constraint to respect: **no new
> datasets/attrs may be created after `start_swmr_write`** — so the monitor must
> define the full channel layout from the `-node`/`-ele`/`-resp` selection *before*
> arming SWMR (selection is fixed at `monitor start`, which we already require).

> [!check] Overhead: **throttle in the engine, not just the sidecar** (my call on #2).
> The append path is gated in the recorder itself by two cheap guards before any HDF5
> touch: `-every K` step decimation **and** a wall-clock min-interval (default ~33 ms /
> 30 Hz) using the profiler's existing `PerfClock` — so even a 10⁶-step explicit run
> appends at most ~30 frames/s regardless of step rate. The sidecar throttle stays as a
> second line of defense for the *reader*. Acceptance gate: a microbench (monitored vs
> unmonitored identical transient) must keep added wall time within the profiler's
> <0.5% phase-level budget; if SWMR-append alone blows it, fall back to batching frames
> (append every M frames in one `H5Dwrite`). Add this microbench to the test suite.

> [!check] Network (P1.5): **kept outside C++.** No TCP/WebSocket sink in the engine.
> The engine only ever does local SWMR-append. For remote/headless/MPI runs, the
> FastAPI sidecar runs *next to* the engine on the remote box, tails the local `.h5`,
> and only the SSE/WS hop crosses the network to the laptop. `LadrunoLiveSink` therefore
> needs **one** impl (local SWMR-HDF5); the "transport seam" is really just the
> file-path contract between engine and sidecar. Simpler than originally scoped.

- Class tag must be reserved in `LEDGER_implementations.md` before coding (≥33000 band).
- Ledgers: new recorder → `LEDGER_implementations.md` row; any `// Ladruno` hook in
  `Profiler.cpp`/`OpenSeesCommands.cpp` → `LEDGER_vanilla_files.md` row + in-source
  comment; SWMR/locking lessons → `LEDGER_quirks.md`. Banner line for the shipped row.

## Implementation log

- 2026-05-31 — Plan drafted from scoping session. Decisions: Shape B (engine push),
  local-first sink behind a transport seam, scalar time-histories only for v1.
- 2026-05-31 — **v1 BUILT + TESTED (green).** Shipped the Shape-B nodal monitor:
  `RECORDER_TAGS_LadrunoMonitorRecorder=33002`; `LadrunoMonitorSink.{h,cpp}`
  (SWMR-HDF5, chunked unlimited `FRAMES`/`STEP`/`TIME` + self-describing `COLUMNS`,
  `H5Fstart_swmr_write` after layout, per-frame `H5Dflush`, fire-and-forget append);
  `LadrunoMonitorRecorder.{h,cpp}` (node/dof + disp|vel|accel|reaction, `-every K`
  + `-hz H` steady_clock gates, lazy open on first commit). Registered via the
  existing `recorder` command as type **`Monitor`** (reuses addRecorder /
  `remove recorder $tag` = stop / tag-as-result) — minimal vanilla footprint
  (classTags.h, OpenSeesOutputCommands.cpp recordersMap, recorder/CMakeLists.txt).
  Built clean on Windows/oneAPI (OpenSeesPy). `Ladruno_scripts/test_monitor_smoke.py`
  green: SMOKE (200 frames, columns) · PARITY 5.55e-17 vs full-precision NodeRecorder ·
  SWMR reopen+refresh in `swmr=True` · `-every 10` → 20 frames stride-10 · overhead
  ~70 us/emitted frame (3× H5Dflush dominated; gates bound it).
  Quirk found: `getCommitTag()` is a GLOBAL monotonic counter — does NOT reset on
  `wipe()` (the monitor `STEP` axis keeps climbing across analyses in one session).
- **NEXT (P1):** FastAPI `/live/{run}/stream` SSE tailing the SWMR file + React live
  mode (true cross-process concurrent tail exercised there). Then element/region
  channels, the profiler-telemetry tee, and the parallel path. Overhead-tuning option
  if needed: batch frames (append every M, one `H5Dwrite`) instead of per-frame flush.
- 2026-05-31 — Open questions resolved: (1) sink = SWMR-HDF5 (one file, both
  lifetimes); (2) throttle in the engine via `-every` + `-hz` wall-clock gate on
  `PerfClock`, with an overhead microbench as acceptance gate; (3) network kept
  outside C++ — single local sink impl, sidecar-next-to-engine + SSE for remote.
  Plan is now decision-complete and ready to implement.
