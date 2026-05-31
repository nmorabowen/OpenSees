# Profiler viewer / backend (`Ladruno_tools/profiler_viewer`)

Python backend + (later) React frontend for the OpenSees/Ladruno stack profiler.
Design doc: [`../../Ladruno_implementation/06_profiler.md`](../../Ladruno_implementation/06_profiler.md).

This is the **schema-first prototype**: it pins the on-disk HDF5 contract and the
analysis/diff API *before* any C++ is written, so the in-engine writer (phase P1) and
this backend (P7) agree by construction. It runs today with only `h5py` + `numpy`.

## Files

| File | Role |
|------|------|
| `profiler_schema.py` | **The contract.** HDF5 layout, dtypes, attr names, and the reference `write_run()`. The C++ writer must mirror this exactly. |
| `profiler_results.py` | Backend loader `ProfilerResults`: manifest, normalized rollup, series, memory, and `diff(base, cand)`. Returns the JSON wire format the React frontend consumes. |
| `make_sample.py` | Synthesizes `profile_sample.h5` (3 runs) and prints manifest + diff. Frontend dev fixture before the C++ writer exists. |
| `test_contract.py` | Plain-assert guard: round-trip, run immutability, normalizers, `elem_by_type`, diff. |

## Quick start

```bash
python make_sample.py     # writes profile_sample.h5, prints manifest + a full→modified-Newton diff
python test_contract.py   # asserts the contract holds
```

```python
from profiler_results import ProfilerResults
with ProfilerResults("profile_sample.h5") as pr:
    pr.manifest()              # picker rows  (GET /runs)
    pr.rollup("caseA")         # flame graph  (GET /runs/{id}/rollup)
    pr.series("caseA")         # time history (GET /runs/{id}/series)
    pr.diff("caseA", "caseB")  # prove a fix  (GET /diff?base=&cand=)
```

## Data flow

```
C++ engine ──write_run()──► profile.h5 ──ProfilerResults──► JSON ──► React frontend
            (HDF5 at rest, one group per run)        (this package)   (P8, not yet built)
```

## What the C++ writer (P1) must honor

- **One immutable group per run** under `/runs/<id>`; `report` *appends*, never mutates.
- **Node path is the diff key.** Emit the rollup as nested groups so
  `/rollup/step/solveCurrentStep/formTangent` is byte-identical across runs (P0#6).
  Keep child-group creation order deterministic.
- **Times in nanoseconds** (`wall_ns`, `cpu_ns`, min/max) as node attrs; the backend
  converts to ms and derives `share` / per-step / per-DOF (P0#2). Do not pre-normalize.
- **Run-header normalizers are mandatory**: `nDOF` (`theSOE->getNumEqn()`), `nElem`,
  `nNode`, `nnz`, `dt_cr`, `oversample_ratio`, `threads` (P0#2/#5). Without them the
  rollup numbers are not portable across runs/machines.
- **`series` is wall-only per phase** (Windows CPU clock is ~15 ms coarse); chunk the
  datasets and flush periodically so a mid-run crash keeps partial data
  (`project_mpco_exit_crash` lesson — `.mpco` is left unflushed on `exit(-1)`).
- **`elem_by_type`** (deep mode) attaches to `formTangent`/`formUnbalance` with the
  `fb_coupled` flag so force-based vs displacement-based elements aren't misread (P0#1, T4).

## Not yet built

- FastAPI wrapper exposing the `ProfilerResults` methods as HTTP endpoints (P7).
- React frontend: flame/icicle, time-series, run-comparison, dt_cr/leak badges (P8).
- Steady-state windowing / `min`-statistic comparison surfaced in the diff (P1#8).
