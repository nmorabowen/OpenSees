"""
profiler_results.py  —  Python backend loader for profiler HDF5 datasets.

This is the headless analysis API (usable directly from Jupyter) AND the data layer the
FastAPI endpoints (P7) wrap. It reads the contract defined in profiler_schema.py, applies
the normalizers the FEM panel required (P0#2: per-step, per-DOF), and implements the
run-to-run diff (P0#6) by walking the rollup group trees along their stable HDF5 paths.

The dicts returned here ARE the JSON wire format served to the React frontend (P8); the
frontend never opens HDF5 itself.

    from profiler_results import ProfilerResults
    pr = ProfilerResults("profile.h5")
    pr.manifest()              # -> [{run, model, nDOF, wall_ms_total, ...}]
    pr.rollup("caseA")         # -> nested dict (flame graph)
    pr.series("caseA")         # -> {phases, step, t, dt, iters, wall_ms, mem_*}
    pr.diff("caseA", "caseB")  # -> [{path, base_ms, cand_ms, d_abs, d_pct, d_share}]
"""

from __future__ import annotations
import h5py
import numpy as np
from profiler_schema import RUN_ATTRS_STR, RUN_ATTRS_INT, RUN_ATTRS_FLOAT, NODE_ATTRS_INT

_NS_PER_MS = 1.0e6


class ProfilerResults:
    def __init__(self, path: str):
        self.path = path
        self._f = h5py.File(path, "r")
        if "runs" not in self._f:
            raise ValueError(f"{path}: no /runs group — not a profiler dataset")

    # -- lifecycle ----------------------------------------------------------
    def close(self):
        self._f.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    @property
    def schema_version(self) -> str:
        return self._f.attrs.get("schema_version", "?")

    # -- runs ---------------------------------------------------------------
    def runs(self) -> list[str]:
        return list(self._f["runs"].keys())

    def header(self, run: str) -> dict:
        g = self._f["runs"][run]
        h = {"run": run}
        for k in RUN_ATTRS_STR:
            h[k] = g.attrs.get(k, "")
        for k in RUN_ATTRS_INT:
            h[k] = int(g.attrs.get(k, 0))
        for k in RUN_ATTRS_FLOAT:
            h[k] = float(g.attrs.get(k, 0.0))
        return h

    def manifest(self) -> list[dict]:
        """Backend GET /runs — one summary row per run for the picker."""
        keep = ("run", "model", "integrator", "algorithm", "solver",
                "nDOF", "nElem", "nSteps", "threads",
                "dt_min", "dt_cr", "oversample_ratio", "wall_ms_total", "timestamp")
        return [{k: self.header(r)[k] for k in keep} for r in self.runs()]

    # -- rollup -------------------------------------------------------------
    def rollup(self, run: str) -> dict:
        """Backend GET /runs/{run}/rollup — normalized nested tree (flame graph)."""
        hdr = self.header(run)
        rg = self._f["runs"][run]["rollup"]
        roots = list(rg.keys())
        root = self._read_node(rg[roots[0]], hdr, root_wall_ns=None)
        return root

    def _read_node(self, ng, hdr: dict, root_wall_ns: int | None) -> dict:
        wall_ns = int(ng.attrs.get("wall_ns", 0))
        if root_wall_ns is None:
            root_wall_ns = wall_ns or 1
        nSteps = max(hdr["nSteps"], 1)
        nDOF = max(hdr["nDOF"], 1)
        node = {
            "name": ng.name.rsplit("/", 1)[-1],
            "path": ng.name,                       # stable diff key
            "calls": int(ng.attrs.get("calls", 0)),
            "wall_ms": wall_ns / _NS_PER_MS,
            "wall_ms_min": int(ng.attrs.get("wall_ns_min", 0)) / _NS_PER_MS,
            "wall_ms_max": int(ng.attrs.get("wall_ns_max", 0)) / _NS_PER_MS,
            "cpu_ms": int(ng.attrs.get("cpu_ns", 0)) / _NS_PER_MS,
            "share": wall_ns / root_wall_ns,                     # fraction of root
            "wall_ms_per_step": wall_ns / _NS_PER_MS / nSteps,   # P0#2 normalizer
            "wall_ms_per_dof": wall_ns / _NS_PER_MS / nDOF,      # P0#2 normalizer
        }
        if "elem_by_type" in ng:                                 # P0#1 per-class breakdown
            arr = ng["elem_by_type"][()]
            node["elem_by_type"] = [
                {"classTag": int(r["classTag"]), "count": int(r["count"]),
                 "wall_ms": int(r["wall_ns"]) / _NS_PER_MS,
                 "wall_ns_per_elem": float(r["wall_ns_per_elem"]),
                 "fb_coupled": bool(r["fb_coupled"])}
                for r in arr
            ]
        node["children"] = [self._read_node(ng[c], hdr, root_wall_ns)
                            for c in ng if isinstance(ng[c], h5py.Group)]
        return node

    # -- series -------------------------------------------------------------
    def series(self, run: str) -> dict | None:
        """Backend GET /runs/{run}/series — per-step arrays (time-history plots)."""
        g = self._f["runs"][run]
        if "series" not in g:
            return None
        sg = g["series"]
        out = {k: sg[k][()].tolist() for k in sg if k != "wall_ms"}
        wall = sg["wall_ms"]
        out["phases"] = list(wall.attrs.get("phases", []))
        out["wall_ms"] = wall[()].tolist()         # [nSteps][nPhase]
        return out

    def memory(self, run: str) -> dict:
        mg = self._f["runs"][run]["memory"]
        out = {k: int(mg.attrs.get(k, 0))
               for k in ("matrix_live", "vector_live", "id_live", "peak_bytes")}
        comps = mg["components_live"][()]
        out["components_live"] = [{"classTag": int(r["classTag"]), "count": int(r["count"])}
                                  for r in comps]
        return out

    # -- diff (P0#6) --------------------------------------------------------
    def diff(self, base: str, cand: str) -> list[dict]:
        """Backend GET /diff?base=&cand= — per-node delta by stable path.

        Walks both rollup trees, unions their node paths (relative to each root), and
        reports absolute / percent / share deltas in wall time. This is what turns the
        profiler from 'observe one run' into 'prove a change helped'.
        """
        b = _flatten(self.rollup(base))
        c = _flatten(self.rollup(cand))
        paths = sorted(set(b) | set(c), key=lambda p: (p.count("/"), p))
        rows = []
        for p in paths:
            bn, cn = b.get(p), c.get(p)
            b_ms = bn["wall_ms"] if bn else 0.0
            c_ms = cn["wall_ms"] if cn else 0.0
            d_abs = c_ms - b_ms
            d_pct = (d_abs / b_ms * 100.0) if b_ms else float("inf")
            rows.append({
                "path": p,
                "name": (cn or bn)["name"],
                "base_ms": round(b_ms, 4),
                "cand_ms": round(c_ms, 4),
                "d_abs_ms": round(d_abs, 4),
                "d_pct": (round(d_pct, 2) if d_pct != float("inf") else None),
                "d_share": round((cn["share"] if cn else 0) - (bn["share"] if bn else 0), 4),
                "status": "added" if not bn else "removed" if not cn else "changed",
            })
        return rows


def _flatten(node: dict, rel: str = "", out: dict | None = None) -> dict:
    """Map relative-path -> node, so two runs align even if absolute /runs/<id> differs."""
    out = {} if out is None else out
    rel = f"{rel}/{node['name']}" if rel else node["name"]
    out[rel] = node
    for ch in node.get("children", []):
        _flatten(ch, rel, out)
    return out
