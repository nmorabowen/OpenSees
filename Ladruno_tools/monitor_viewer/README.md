# Ladruno Monitor Viewer

Viewers for the **Analysis monitor** stream — the SWMR-HDF5 sink written live by
`recorder Monitor` (`ladruno::MonitorSink`) while an analysis runs. Two options,
both pure tooling (no OpenSees rebuild):

1. **`monitor_view.py`** — a zero-web-stack matplotlib plotter (static or live tail).
2. **`monitor_server.py` + `monitor_page.html`** — a live browser dashboard that
   polls for new frames and draws them as auto-updating line charts.

Both sit on **`monitor_reader.py`** (`MonitorReader`), which opens the sink in
SWMR read mode and refreshes on each read, so it tails a file another process is
still writing.

## Sink format (read-only contract)

```
/            attrs: FORMAT="ladruno-monitor", FORMAT_VERSION, GENERATOR, units?
COLUMNS      [nCols]      vlen strings   self-describing channels, e.g. node2.disp.dof1
STEP         [T]          int32          commitTag per frame
TIME         [T]          float64        pseudo-time per frame
FRAMES       [T, nCols]   float64        one row per emitted frame
```

## Produce a sink

```python
ops.recorder("Monitor", "-node", 2, 5, "-dof", 1, "-resp", "disp",
             "-sink", "run_live.h5", "-every", 1)   # or "-hz", 40
# ... analyze ...  (frames stream to run_live.h5 as the run commits steps)
```

## CLI plotter

```bash
python monitor_view.py run_live.h5                 # all channels vs TIME (window)
python monitor_view.py run_live.h5 --save out.png  # headless -> PNG
python monitor_view.py run_live.h5 --x step        # x-axis = STEP index
python monitor_view.py run_live.h5 --watch         # live tail; redraw as frames arrive
python monitor_view.py run_live.h5 --channels node2.disp.dof1,node5.disp.dof1
```

## Live web dashboard

```bash
python monitor_server.py run_live.h5 --port 8030   # opens http://127.0.0.1:8030
# or:  MONITOR_H5=run_live.h5 uvicorn monitor_server:app --port 8030
```

Open the URL while the analysis is running — the page shows a **LIVE** pill that
pulses while frames arrive (and flips to *idle* when the run is at rest), a
multi-channel line chart (toggle x-axis time/step, Follow, Pause), and a per-
channel latest / min / max table. It works the same on an at-rest sink file.

Endpoints (thin layer over `MonitorReader`, file opened per request):

| route | returns |
|---|---|
| `GET /` | the dashboard page |
| `GET /health` | `{status, file, nframes}` |
| `GET /api/meta` | columns + counts + units |
| `GET /api/frames?since=N` | incremental frames `[N:n)` — the poll endpoint |

## Requirements

`h5py` + `numpy` for the reader/CLI; `fastapi` + `uvicorn` for the server;
`matplotlib` for the CLI plots. Set `HDF5_USE_FILE_LOCKING=FALSE` (the reader
does this on import) so a reader can open a file the writer still holds.

## Test

```bash
python test_monitor_view.py   # binary-free: synthetic sink -> reader + API endpoints
```

## Not this viewer

Perf/profiling (phase tree, flame graph, memory) is a **different** stream — see
`../profiler_viewer/`. This viewer is only for the nodal-results monitor sink.
