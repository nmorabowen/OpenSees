"""
profiler_schema.py  —  canonical on-disk contract for the OpenSees/Ladruno profiler.

This module IS the contract. The in-engine C++ HDF5 writer (phase P1) must produce
exactly this layout, and the Python backend (P7) reads it via profiler_results.py.
`write_run()` below is the reference writer — keep it and the C++ writer in lockstep.

HDF5 layout
-----------
profile.h5
  / (root)
      attrs: schema_version : str

  /runs/<run_id>                         one group per analysis (the "dataset" of runs)
      attrs (run header / normalizers):
        model, engine_sha, integrator, algorithm, solver, units, timestamp : str
        nSteps, nDOF, nElem, nNode, nnz, threads                           : int
        dt_min, dt_max, dt_cr, oversample_ratio                            : float
        wall_ms_total, cpu_ms_total                                        : float

    /rollup/<root>/...                   profile node tree as nested groups
        per-node attrs: calls, wall_ns, wall_ns_min, wall_ns_max, cpu_ns
        a node may carry one optional child dataset:
          elem_by_type : compound [classTag, count, wall_ns, wall_ns_per_elem, fb_coupled]
                         (deep mode, P0#1; attached to formTangent / formUnbalance)
        node *path* (e.g. /rollup/step/solveCurrentStep/formTangent) is the STABLE
        DIFF KEY (P0#6): comparing two runs == walking both trees by identical path.

    /series/                             per-step time history (chunked); WALL ONLY for
                                         per-phase time (Windows CPU clock is ~15 ms coarse)
        step           : int64  [nSteps]
        t              : float64[nSteps]
        dt             : float64[nSteps]
        iters          : int32  [nSteps]      (P0#3)
        mem_live_bytes : int64  [nSteps]
        mem_peak_bytes : int64  [nSteps]
        wall_ms        : float32[nSteps, nPhase]   attr 'phases' = phase name list

    /memory/                             end-of-run memory snapshot
        attrs: matrix_live, vector_live, id_live, peak_bytes : int
        components_live : compound [classTag, count]   (TaggedObject census)
"""

from __future__ import annotations
import numpy as np

SCHEMA_VERSION = "0.1.0"

# Compound dtypes (shared by writer and reader so the binary layout is single-sourced).
ELEM_BY_TYPE_DT = np.dtype([
    ("classTag",          "i4"),
    ("count",             "i8"),
    ("wall_ns",           "i8"),
    ("wall_ns_per_elem",  "f8"),
    ("fb_coupled",        "u1"),   # 1 = force-based (tangent/residual coupled), 0 = displacement-based
])

COMPONENTS_LIVE_DT = np.dtype([
    ("classTag", "i4"),
    ("count",    "i8"),
])

# Run-header attribute keys, grouped by type so the C++ writer and reader agree on names.
RUN_ATTRS_STR   = ("model", "engine_sha", "integrator", "algorithm", "solver", "units", "timestamp")
RUN_ATTRS_INT   = ("nSteps", "nDOF", "nElem", "nNode", "nnz", "threads")
RUN_ATTRS_FLOAT = ("dt_min", "dt_max", "dt_cr", "oversample_ratio", "wall_ms_total", "cpu_ms_total")

NODE_ATTRS_INT = ("calls", "wall_ns", "wall_ns_min", "wall_ns_max", "cpu_ns")

SERIES_SCALAR = {  # name -> dtype
    "step": "i8", "t": "f8", "dt": "f8", "iters": "i4",
    "mem_live_bytes": "i8", "mem_peak_bytes": "i8",
}


# --------------------------------------------------------------------------- writer
def write_run(h5file, run_id: str, *, meta: dict, rollup: dict,
              series: dict | None = None, memory: dict | None = None) -> None:
    """Append one run as /runs/<run_id>. Never mutates existing runs (P0 multi-run model).

    rollup node dict: {name, calls, wall_ns, wall_ns_min, wall_ns_max, cpu_ns,
                       children:[...], elem_by_type:[(classTag,count,wall_ns,per_elem,fb)] }
    """
    h5file.attrs.setdefault("schema_version", SCHEMA_VERSION)
    runs = h5file.require_group("runs")
    if run_id in runs:
        raise ValueError(f"run '{run_id}' already exists; runs are immutable — pick a new id")
    g = runs.create_group(run_id)

    for k in RUN_ATTRS_STR:
        g.attrs[k] = str(meta.get(k, ""))
    for k in RUN_ATTRS_INT:
        g.attrs[k] = int(meta.get(k, 0))
    for k in RUN_ATTRS_FLOAT:
        g.attrs[k] = float(meta.get(k, 0.0))

    _write_node(g.create_group("rollup"), rollup)

    if series:
        sg = g.create_group("series")
        n = len(series["step"])
        for name, dt in SERIES_SCALAR.items():
            sg.create_dataset(name, data=np.asarray(series[name], dtype=dt),
                              chunks=(min(n, 4096),) if n else None)
        wall = np.asarray(series["wall_ms"], dtype="f4")          # [nSteps, nPhase]
        ds = sg.create_dataset("wall_ms", data=wall,
                               chunks=(min(n, 4096), wall.shape[1]) if n else None)
        ds.attrs["phases"] = list(series["phases"])

    mg = g.create_group("memory")
    memory = memory or {}
    for k in ("matrix_live", "vector_live", "id_live", "peak_bytes"):
        mg.attrs[k] = int(memory.get(k, 0))
    comps = np.asarray(memory.get("components_live", []), dtype=COMPONENTS_LIVE_DT)
    mg.create_dataset("components_live", data=comps)


def _write_node(parent, node: dict) -> None:
    ng = parent.create_group(node["name"])
    for k in NODE_ATTRS_INT:
        ng.attrs[k] = int(node.get(k, 0))
    if node.get("elem_by_type"):
        arr = np.asarray([tuple(r) for r in node["elem_by_type"]], dtype=ELEM_BY_TYPE_DT)
        ng.create_dataset("elem_by_type", data=arr)
    for child in node.get("children", []):
        _write_node(ng, child)
