# Profiler frontend (P8)

React + TypeScript + Vite viewer for the OpenSees/Ladruno stack profiler. It is a
pure client of the P7 FastAPI backend (`../profiler_api.py`) — it never opens
HDF5 itself; every panel renders the JSON the backend serves.

## Panels

- **Run picker** — manifest table; pick a run to *view* and (optionally) a run to
  compare *vs* for the diff.
- **Rollup** — hand-rolled SVG icicle/flame graph of the phase tree. Click a frame
  to inspect its metrics and, in deep runs, its per-element-class breakdown
  (`elem_by_type`, with the force-based / disp-based tag).
- **Series** — per-step time history: per-phase wall time, iterations/step, dt/step. With a *vs* run selected, that run's curves overlay dashed/faded for comparison.
- **Memory** — per-type live byte counters + peak, and the TaggedObject
  live-component census (classTag → count). With a *vs* run selected, a leak
  badge + Δ column flag classTags with more live objects than the baseline.
- **Diff** — per-node wall-time delta along the stable rollup path (green = the
  viewed run is faster than the baseline).
- **Badges** — scheme / size / `dt vs dt_cr` / oversample, with a hint when the
  step is finer than critical (a cheap explicit speed-up).

No charting dependency: the icicle and the line charts are plain SVG.

## Run it

```bash
npm install

# 1) start the backend (sibling dir) against any profile.h5:
#      python ../profiler_api.py --file ../profile_sample.h5
#    (generate a demo file first with:  python ../make_sample.py)
# 2) start the dev server:
npm run dev          # http://localhost:5173
```

The backend defaults to `http://127.0.0.1:8000` and sets permissive CORS, so the
dev server talks to it directly. Point the UI at a different backend with
`VITE_API_BASE`:

```bash
VITE_API_BASE=http://host:9000 npm run dev
```

## Build / check

```bash
npm run build        # tsc -b (typecheck) + vite production build -> dist/
npm run lint
```

## Layout

```
src/
  api.ts                 typed fetch client + wire-format types (mirror profiler_results.py)
  format.ts              ms / bytes / pct / sci formatters + frame color
  App.tsx                orchestration: run selection, data loading, tabs
  components/
    RunPicker.tsx        manifest table (view + vs selection)
    Badges.tsx           header badges (dt_cr / oversample levers)
    Icicle.tsx           SVG flame/icicle of the rollup tree
    NodeDetails.tsx      selected-node metrics + elem_by_type table
    LineChart.tsx        generic responsive SVG multi-line chart (downsamples long series)
    TimeSeries.tsx       per-phase wall / iters / dt charts
    MemoryPanel.tsx      live byte counters + census bars
    DiffTable.tsx        per-node run-to-run delta
```
