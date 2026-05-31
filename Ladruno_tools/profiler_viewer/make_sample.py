"""
make_sample.py  —  synthesize a profile.h5 and validate the contract end-to-end.

Generates three runs into one dataset:
  caseA  : implicit, full Newton            (baseline)
  caseB  : implicit, modified Newton        (candidate: tangent reused -> formTangent cheaper)
  expl   : explicit central-difference      (same taxonomy: residual-dominated, trivial solve)

Then reads them back through ProfilerResults and prints manifest + diff(caseA, caseB).
This is also the fixture the React frontend (P8) develops against before the C++ writer exists.

Run:  python make_sample.py
"""

from __future__ import annotations
import numpy as np
import h5py
from profiler_schema import write_run
from profiler_results import ProfilerResults

MS = 1.0e6  # ns per ms
rng = np.random.default_rng(0)


def node(name, wall_ms, calls, children=None, elem_by_type=None, cpu_frac=0.99):
    wall_ns = int(wall_ms * MS)
    return {
        "name": name, "calls": calls,
        "wall_ns": wall_ns,
        "wall_ns_min": int(wall_ns / max(calls, 1) * 0.7),
        "wall_ns_max": int(wall_ns / max(calls, 1) * 1.6),
        "cpu_ns": int(wall_ns * cpu_frac),
        "children": children or [],
        "elem_by_type": elem_by_type or [],
    }


def elem_rows(rows):
    # rows: [(classTag, count, wall_ms, fb_coupled)]
    out = []
    for tag, cnt, wms, fb in rows:
        wns = int(wms * MS)
        out.append((tag, cnt, wns, wns / max(cnt, 1), fb))
    return out


# Per-class assembly breakdown (P0#1): 200 force-based beams dominate despite 20k cheap trusses.
def tangent_elems(scale):
    return elem_rows([
        (74,   200,   scale * 0.80, 1),   # forceBeamColumn — fb_coupled, expensive per-elem
        (53,   400,   scale * 0.15, 0),   # ShellMITC4
        (12, 20000,   scale * 0.05, 0),   # truss — many, cheap
    ])


def residual_elems(scale):
    return elem_rows([
        (74,   200,   scale * 0.78, 1),
        (53,   400,   scale * 0.16, 0),
        (12, 20000,   scale * 0.06, 0),
    ])


def implicit_rollup(form_tangent_ms, factor_ms):
    """One step-rollup tree. form_tangent_ms/factor_ms vary full vs modified Newton."""
    return node("step", form_tangent_ms + factor_ms + 2_326.0, 1200, children=[
        node("newStep", 50.0, 1200),
        node("solveCurrentStep", form_tangent_ms + factor_ms + 2_006.0, 1200, children=[
            node("formTangent", form_tangent_ms, 3600, elem_by_type=tangent_elems(form_tangent_ms)),
            node("linearSolve", factor_ms + 350.0, 3600, children=[
                node("solve.factor",  factor_ms, 3600),     # P1#7 split
                node("solve.triSolve", 350.0,   3600),
            ]),
            node("update", 1_000.0, 3600),
            node("formUnbalance", 1_200.0, 7200, elem_by_type=residual_elems(1_200.0)),
            node("convTest", 60.0, 7200),
        ]),
        node("commit", 270.0, 1200),
    ])


def explicit_rollup():
    # residual-dominated, trivial M^-1 solve, no convergence test
    return node("step", 9_050.0, 20000, children=[
        node("newStep", 400.0, 20000),
        node("solveCurrentStep", 8_400.0, 20000, children=[
            node("formTangent", 30.0, 1, elem_by_type=tangent_elems(30.0)),  # factorOnce: lumped mass
            node("formUnbalance", 7_600.0, 20000, elem_by_type=residual_elems(7_600.0)),
            node("linearSolve", 120.0, 20000),                               # trivial diagonal
            node("update", 650.0, 20000),
        ]),
        node("commit", 250.0, 20000),
    ])


def make_series(nSteps, dt, phase_means):
    phases = list(phase_means.keys())
    wall = np.column_stack([
        np.clip(rng.normal(m, m * 0.15, nSteps), 0, None) for m in phase_means.values()
    ]).astype("f4")
    return {
        "phases": phases,
        "step": np.arange(1, nSteps + 1, dtype="i8"),
        "t": np.arange(1, nSteps + 1) * dt,
        "dt": np.full(nSteps, dt),
        "iters": rng.integers(2, 5, nSteps).astype("i4"),
        "mem_live_bytes": (10_000_000 + rng.integers(0, 200_000, nSteps)).astype("i8"),  # flat = no leak
        "mem_peak_bytes": np.full(nSteps, 12_582_912, dtype="i8"),
        "wall_ms": wall,
    }


def base_meta(**kw):
    m = dict(model="frame_5story", engine_sha="deadbeef", units="N,mm,s",
             timestamp="2026-05-30T12:00:00", nDOF=50000, nElem=20600, nNode=8400,
             nnz=2_400_000, threads=8, cpu_ms_total=0.0)
    m.update(kw)
    return m


def main():
    path = "profile_sample.h5"
    with h5py.File(path, "w") as f:
        write_run(f, "caseA",
                  meta=base_meta(integrator="Newmark", algorithm="Newton", solver="Mumps",
                                 nSteps=1200, dt_min=0.005, dt_max=0.005, dt_cr=0.0,
                                 oversample_ratio=0.0, wall_ms_total=8421.0),
                  rollup=implicit_rollup(form_tangent_ms=3_950.0, factor_ms=2_250.0),
                  series=make_series(1200, 0.005,
                                     {"formTangent": 3.29, "linearSolve": 2.17,
                                      "update": 0.83, "formUnbalance": 1.0, "convTest": 0.05}),
                  memory=dict(matrix_live=1422, vector_live=880, id_live=300, peak_bytes=12_582_912,
                              components_live=[(74, 200), (53, 400), (12, 20000), (1, 8400)]))

        write_run(f, "caseB",
                  meta=base_meta(integrator="Newmark", algorithm="ModifiedNewton", solver="Mumps",
                                 nSteps=1200, dt_min=0.005, dt_max=0.005, dt_cr=0.0,
                                 oversample_ratio=0.0, wall_ms_total=5910.0),
                  rollup=implicit_rollup(form_tangent_ms=1_300.0, factor_ms=750.0),
                  series=make_series(1200, 0.005,
                                     {"formTangent": 1.08, "linearSolve": 0.63,
                                      "update": 0.83, "formUnbalance": 1.0, "convTest": 0.05}),
                  memory=dict(matrix_live=1390, vector_live=860, id_live=300, peak_bytes=12_582_912,
                              components_live=[(74, 200), (53, 400), (12, 20000), (1, 8400)]))

        write_run(f, "expl",
                  meta=base_meta(integrator="CentralDifference", algorithm="Linear", solver="Diagonal",
                                 nSteps=20000, dt_min=2e-6, dt_max=2e-6, dt_cr=6e-6,
                                 oversample_ratio=3.0, wall_ms_total=9050.0),
                  rollup=explicit_rollup(),
                  series=make_series(20000, 2e-6,
                                     {"formTangent": 0.0, "linearSolve": 0.006,
                                      "update": 0.033, "formUnbalance": 0.38, "convTest": 0.0}),
                  memory=dict(matrix_live=900, vector_live=700, id_live=300, peak_bytes=9_400_000,
                              components_live=[(74, 200), (53, 400), (12, 20000), (1, 8400)]))

    print(f"wrote {path}\n")
    pr = ProfilerResults(path)
    print(f"schema_version: {pr.schema_version}")
    print("\n=== manifest ===")
    for row in pr.manifest():
        print(f"  {row['run']:7s} {row['integrator']:>16s}/{row['algorithm']:<14s} "
              f"nDOF={row['nDOF']:>6d} nSteps={row['nSteps']:>5d} "
              f"oversample={row['oversample_ratio']:>4.1f} wall={row['wall_ms_total']:>8.1f} ms")

    print("\n=== diff caseA (full Newton) -> caseB (modified Newton) ===")
    print(f"  {'node':<42s} {'base_ms':>10s} {'cand_ms':>10s} {'d_pct':>8s} {'dshare':>8s}")
    for r in pr.diff("caseA", "caseB"):
        pct = "n/a" if r["d_pct"] is None else f"{r['d_pct']:+.1f}%"
        print(f"  {r['path']:<42s} {r['base_ms']:>10.1f} {r['cand_ms']:>10.1f} "
              f"{pct:>8s} {r['d_share']:>+8.3f}")

    print("\n=== caseA formTangent per-element-class (P0#1) ===")
    ft = next(c for c in pr.rollup("caseA")["children"]
              if c["name"] == "solveCurrentStep")
    ft = next(c for c in ft["children"] if c["name"] == "formTangent")
    for e in ft["elem_by_type"]:
        flag = "force-based" if e["fb_coupled"] else "disp-based"
        print(f"  classTag {e['classTag']:>4d}  n={e['count']:>6d}  {e['wall_ms']:>8.1f} ms  "
              f"{e['wall_ns_per_elem']:>10.0f} ns/elem  [{flag}]")
    pr.close()


if __name__ == "__main__":
    main()
