"""ADR-75 P1c — does the desktop-PARDISO win GROW with model size?

The production regime for this fork is **huge solid nonlinear models**, but every
ADR-75 gate so far was set on an 11.5k-DOF Lane-B model. A 3D sparse-direct
factorization scales ~O(N^1.5-2) in flops while assembly is ~O(N), so the solve
fraction — and therefore the value of a better/threaded solver — should RISE
with N. This script tests that directly instead of assuming it.

Sweeps a cubic LadrunoBrick + LadrunoJ2 mesh over n = 15, 20, 25, 30
(~11.5k / 26.5k / 51k / 87k DOF) and times UmfPack vs threaded PARDISO on each.

Metric: wall seconds for NSTEPS load steps (build excluded), median of REPEATS.
Correctness: PARDISO must match UmfPack's tip displacement to 1e-9 relative at
every size, else the timing is meaningless.

Run:
  OPS_PYD=<dist\\bin> MKL_NUM_THREADS=4 python3.12 -S laneB_scaling.py
Env knobs: OPS_SCALE_SIZES (default "15,20,25,30"), OPS_SCALE_STEPS (3),
           OPS_SCALE_REPEATS (3).
"""
from __future__ import annotations
import json
import os
import statistics
import sys
import time
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", os.environ["MKL_NUM_THREADS"])
os.environ.setdefault("PMI_RANK", "0")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

DIST = os.environ.get("OPS_PYD", r"C:\Users\nmb\Documents\Github\OpenSees"
                                 r"\.claude\worktrees\mumps-opensees-study-f833bf\dist\bin")
sys.path.insert(0, DIST)
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)
import opensees as ops  # noqa: E402

HERE = Path(__file__).resolve().parent
SIZES = [int(s) for s in os.environ.get("OPS_SCALE_SIZES", "15,20,25,30").split(",")]
NSTEPS = int(os.environ.get("OPS_SCALE_STEPS", "3"))
REPEATS = int(os.environ.get("OPS_SCALE_REPEATS", "3"))
THREADS = os.environ["MKL_NUM_THREADS"]

E, nu = 200000.0, 0.3
K = E / (3.0 * (1.0 - 2.0 * nu))
G = E / (2.0 * (1.0 + nu))
s0, Qinf, b, Hiso = 250.0, 0.0, 0.0, 2000.0
L = 100.0


def build(n, system_args):
    nx = ny = nz = n
    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", s0, Qinf, b, Hiso)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(nid(i, j, k), i * L, j * L, k * L)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.fix(nid(i, j, 0), 1, 1, 1)
    e = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ops.element("LadrunoBrick", e,
                            nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k),
                            nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j+1, k+1), nid(i, j+1, k+1),
                            1)
                e += 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.load(nid(i, j, nz), 4.0e5, 0.0, -4.0e5)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system(*system_args)
    ops.test("NormDispIncr", 1.0e-7, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / 15.0)   # same increment as the P1 bench
    ops.analysis("Static")
    nnode = (nx + 1) * (ny + 1) * (nz + 1)
    ndof = (nnode - (nx + 1) * (ny + 1)) * 3
    return nid(nx // 2, ny // 2, nz), e - 1, ndof


def run_once(n, system_args):
    tip, nele, ndof = build(n, system_args)
    t0 = time.perf_counter()
    for s in range(NSTEPS):
        if ops.analyze(1) != 0:
            raise RuntimeError(f"analyze failed at step {s+1}")
    wall = time.perf_counter() - t0
    return wall, ops.nodeDisp(tip, 1), nele, ndof


def main():
    print(f"Lane-B SCALING | sizes={SIZES} steps={NSTEPS} repeats={REPEATS} "
          f"MKL_NUM_THREADS={THREADS}")
    print(f"{'n':>3} {'nele':>7} {'ndof':>8} {'UmfPack':>9} {'PARDISO':>9} "
          f"{'speedup':>8} {'relerr':>9}")
    out = {"threads": THREADS, "steps": NSTEPS, "repeats": REPEATS, "rows": []}
    for n in SIZES:
        try:
            res = {}
            meta = None
            _pairs = [("UmfPack", ["UmfPack"]), ("Pardiso", ["Pardiso"])]
            _only = [s.strip() for s in os.environ.get("OPS_SCALE_SOLVERS", "").split(",") if s.strip()]
            if _only:
                _pairs = [(n_, a_) for (n_, a_) in _pairs if n_ in _only]
            for name, args in _pairs:
                samples, ux = [], None
                for _ in range(REPEATS):
                    w, ux, nele, ndof = run_once(n, args)
                    samples.append(w)
                    meta = (nele, ndof)
                res[name] = (statistics.median(samples), ux, sorted(samples))
            nele, ndof = meta
            # Either solver may be absent when OPS_SCALE_SOLVERS filters (e.g. a
            # PARDISO-only probe at a size where UmfPack runs out of memory).
            um, uux, us = res.get("UmfPack", (None, None, None))
            pm, pux, ps = res.get("Pardiso", (None, None, None))
            spd = (um / pm) if (um and pm) else None
            err = (abs(pux - uux) / abs(uux)) if (uux and pux) else None
            print(f"{n:>3} {nele:>7} {ndof:>8} "
                  f"{('%9.3f' % um) if um else '        -'} "
                  f"{('%9.3f' % pm) if pm else '        -'} "
                  f"{('%7.2fx' % spd) if spd else '       -'} "
                  f"{('%9.1e' % err) if err is not None else '        -'}")
            out["rows"].append({"n": n, "nele": nele, "ndof": ndof,
                                "umfpack_s": um, "pardiso_s": pm,
                                "speedup": spd, "rel_err": err,
                                "umfpack_samples": us, "pardiso_samples": ps})
        except Exception as ex:  # OOM or solver failure at the top sizes
            print(f"{n:>3} FAILED: {str(ex)[:90]}")
            out["rows"].append({"n": n, "failed": str(ex)[:200]})
            break
    (HERE / "laneB_scaling.json").write_text(json.dumps(out, indent=2))
    print("wrote laneB_scaling.json")


if __name__ == "__main__":
    main()
