---
name: ladruno-viewer-tools
description: >
  Checklist for changing the Ladruno VIEWER TOOLS: Ladruno_tools/profiler_viewer/ (ProfilerResults
  loader, FastAPI profiler_api, React/TypeScript/Vite frontend, launch.py launcher, profiler_monitor,
  smoke checks) and Ladruno_tools/monitor_viewer/ (MonitorReader, monitor_server + monitor_page.html
  dashboard, monitor_view CLI). Use before editing any of them, adding a panel, endpoint or tool, or
  verifying a viewer change in a browser. Covers the data to test against, the checks to run (none
  run in CI), visual verification, and the ledger rule the viewer ledger gate enforces.
---

# Viewer tools — checklist

Read this before changing anything under `Ladruno_tools/`. Items name the
`Ladruno_implementation/LEDGER_quirks.md` heading to grep for; read that entry when the item
applies. **[lint]** items are enforced by `python ci/check_viewer_ledger.py`. Out of scope: using
the viewers (their READMEs) and the C++ engine side (the profiler and `recorder Monitor` sources).

## Contracts other code reads

- [ ] `profiler_schema.py` is the on-disk contract the C++ `ProfilerHDF5Writer` mirrors. A schema
      change is an engine change (build + ledgers), not a viewer-only PR.
- [ ] `frontend/src/api.ts` mirrors the dicts `profiler_results.py` returns. Change both together.
- [ ] The monitor sink layout (`COLUMNS`/`STEP`/`TIME`/`FRAMES`) is written by the engine's
      `LadrunoMonitorSink` (`08_analysis_monitor.md`); `monitor_reader.py` only reads it.
- [ ] apeGmsh imports `profiler_results.ProfilerResults` and runs `launch.py <h5>` from this
      directory (`apeGmsh/src/apeGmsh/profiler.py`, `open()` / `show_web()`). Renaming the module
      or class, or changing `launch.py`'s positional argument, needs a matching apeGmsh change.

## The data you check against

- [ ] An engine-shaped file, not only `make_sample.py`'s: the engine's rollup top is an UNTIMED
      `root`, and panels that look fine on the sample were empty on real runs twice. Quirks:
      "writes a TIMED top node". Binary-free engine shape: `test_contract.py` run "C",
      `monitor_smoke.py`. Real file: a `*_smoke_model.py` run on the built pyd.
- [ ] `STEP` is a monotonic id across every analysis in the process, not "step N of this run".
      Quirks: "`getCommitTag()` is a GLOBAL monotonic counter".
- [ ] A `profile.h5` from several runs in one process holds SUMS unless the script called
      `profiler reset`. Quirks: "The Profiler is a process-global singleton and `ops.wipe()` does
      NOT reset it".

## Checks — no C++ build needed, and none of them run in CI (you are the gate)

```
cd Ladruno_tools/profiler_viewer
python test_contract.py            # h5py + numpy
python monitor_smoke.py            # h5py + numpy; fake ops, engine-shaped rollup
<venv-python> test_api.py          # needs fastapi/httpx (requirements.txt, or launch.py's .viewer_venv)
cd frontend && npm run build && npm run lint   # tsc -b is the type gate
cd ../../monitor_viewer && <venv-python> test_monitor_view.py
```

- [ ] Real engine file: the BUILD Python (3.12, `os.add_dll_directory`; `BUILD_GOTCHAS.md` §0/§4)
      runs `<x>_smoke_model.py <dist\bin> <out>`, the venv runs `<x>_smoke_check.py <out>`.
- [ ] TypeScript: `verbatimModuleSyntax` needs `import type`; `erasableSyntaxOnly` bans enums (use
      string unions). `tsc` reports both; do not relax `tsconfig.app.json` to get green.

## Visual verification — for any change a user can see

- [ ] Launch what you will look at, on a free port:
      profiler `python launch.py <file.h5> --port <p> --no-browser` — it serves the EXISTING
      `frontend/dist`, so after a frontend edit add `--rebuild` or you are looking at the old
      bundle. Dev loop: `python profiler_api.py --file <h5> --port 8000` plus `npm run dev -- --port
      5173 --strictPort` (`api.ts` finds the backend only when the page is on 5173).
      Monitor: `<venv-python> monitor_server.py <sink.h5> --port <p>`; for the live path, write the
      sink from a slowed analysis loop.
- [ ] Reference vs candidate: the same file, viewport and selection (view run, vs run, tab) on
      `ladruno` (a second worktree) and on your branch. Look for empty or single-bar panels,
      clipped or overlapping labels and legends, NaN / Infinity / absurd percentages, and stale
      data after switching runs.
- [ ] Drive every control you touched from code (tab, run / vs pickers, Follow, Pause, x-axis
      toggle) and assert the DOM changed. #487's Follow button rendered and did nothing; its own
      self-review caught it (commit `0ec96549a`).
- [ ] Make the DOM the evidence, not a screenshot. Quirks: "Browser-pane screenshots of the viewers
      timed out in two sessions".
- [ ] Stop what you launched, by PID: `launch.py` starts a uvicorn child and `npm run dev` a node
      process. On Windows `taskkill /PID <pid> /T /F`; then check the port is free.

## Ledgers and docs (same PR)

- [ ] **[lint]** A new file under `Ladruno_tools/` needs an edit to its tool's row in
      `LEDGER_implementations.md`. #35, #53, #485 and #487 skipped it; the last two after the
      lesson was written down (`121_viewer_agent_surface.md`).
- [ ] The gate cannot see a modification-only change. When behaviour changes, update the row and
      the plan-doc log anyway (`06_profiler.md`, `08_analysis_monitor.md`); #55 needed a follow-up
      doc PR (#56) for this.

Found a new trap? Add it to `LEDGER_quirks.md`, then one line here pointing to it. If it names a
greppable pattern, make it a check instead.
