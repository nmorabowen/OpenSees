# Phase-0 lane measurement — proven recipe (do not deviate)

## Run environment (PROVEN working 2026-07-06, smoke test passed)
- Python: `C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe` with `-S` flag (MANDATORY — the boot .pth preloads a stale pyd otherwise)
- Wire the pyd manually at the top of every model script:
```python
import sys, os
DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
sys.path.insert(0, DIST)
os.add_dll_directory(DIST)
os.environ["PMI_RANK"] = "0"
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import opensees as ops
assert ops.__file__.lower().startswith(DIST.lower()), f"WRONG PYD: {ops.__file__}"
```
- Profiler driving (exact syntax, proven):
```python
ops.profiler("start", "-deep", "-perStep")
ops.analyze(...)          # the timed work ONLY — build the model before start
ops.profiler("stop")
ops.profiler("report", OUT_H5_ABS_PATH)
```
- Parse the h5 with `C:\Users\nmora\AppData\Local\Programs\Python\Python311\python.exe`
  (has h5py) using the shared parser:
  `<scratchpad>\parse_prof.py <file.h5>` — prints phase rollup (ms/calls) + per-classTag
  `elem_by_type` buckets (count, wall ms, us/ele).

## Known caveats (do not re-discover)
- `cpu_ms_total` is a Win32 stub — unreliable. Trust wall_ns and RELATIVE splits.
- `nElem`/`nNode` read 0 on this binary (pre-patch) — expected, ignore.
- The binary is `dist\bin\opensees.pyd` (2026-06-18). Do NOT rebuild anything.
- Element/material Tcl names: if an arg signature fails, grep the OPS_ parser in
  `C:\Users\nmora\Github\OpenSees_Compile\OpenSees\.claude\worktrees\competent-engelbart-ab6503\SRC\`
  (worktree source ≈ binary within days; parser signatures are stable).
- Wrap `ops.analyze` in a loop with failure check; if Newton fails mid-run, reduce step size.
- Keep each profiled run < ~3 min wall. Time the whole analyze block with time.perf_counter too.

## Reporting format (per lane)
1. Model card: elements (type/count), DOF, material, steps, algorithm/solver/integrator.
2. Phase table: wall ms + calls for every scope printed by parse_prof.py.
3. elem_by_type rows.
4. Derived split: % of step time in {element state determination (elem.tangent+elem.residual),
   linearSolve, update, convTest, commit, everything-else}. State clearly what the profiler
   could NOT attribute (gap = step total - sum of children).
5. One-line dominance verdict for the lane.
