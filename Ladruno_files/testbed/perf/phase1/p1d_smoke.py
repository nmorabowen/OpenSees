"""ADR-75 P1d smoke: does symmetric PARDISO give the RIGHT answer, and how fast?

Deliberately SMALL (nx=ny=nz=8, ~1.9k DOF) and per-solver timed with a hard
progress print, because the full Lane-B bench went dark for 25+ minutes and gave
no way to tell which of the four solvers was responsible. Correctness first,
timing second — a wrong answer makes a timing meaningless.

Run:
  OPS_PYD=<worktree>\\dist\\bin python3.12 -S p1d_smoke.py
"""
from __future__ import annotations
import os
import sys
import time

sys.stdout.reconfigure(line_buffering=True)
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("PMI_RANK", "0")

DIST = os.environ.get("OPS_PYD", r"C:\Users\nmb\Documents\Github\OpenSees\dist\bin")
sys.path.insert(0, DIST)
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)
import opensees as ops  # noqa: E402
assert ops.__file__.lower().startswith(DIST.lower()), f"WRONG PYD: {ops.__file__}"
print("pyd OK:", ops.__file__, flush=True)

E, nu = 200000.0, 0.3
K = E / (3.0 * (1.0 - 2.0 * nu))
G = E / (2.0 * (1.0 + nu))
s0, Qinf, b, Hiso = 250.0, 0.0, 0.0, 2000.0
NX = int(os.environ.get("P1D_N", "8"))
nx = ny = nz = NX
L = 100.0
NSTEPS = 5


def _nid(i, j, k):
    return 1 + i + (nx + 1) * (j + (ny + 1) * k)


def build(system_args):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", s0, Qinf, b, Hiso)
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(_nid(i, j, k), i * L, j * L, k * L)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.fix(_nid(i, j, 0), 1, 1, 1)
    eTag = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ops.element(
                    "LadrunoBrick", eTag,
                    _nid(i, j, k), _nid(i + 1, j, k), _nid(i + 1, j + 1, k), _nid(i, j + 1, k),
                    _nid(i, j, k + 1), _nid(i + 1, j, k + 1), _nid(i + 1, j + 1, k + 1), _nid(i, j + 1, k + 1),
                    1)
                eTag += 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.load(_nid(i, j, nz), 4.0e5, 0.0, -4.0e5)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system(*system_args)
    ops.test("NormDispIncr", 1e-7, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / NSTEPS)
    ops.analysis("Static")
    return _nid(nx // 2, ny // 2, nz)


CASES = [
    ("UmfPack",    ["UmfPack"]),
    ("Pardiso",    ["Pardiso"]),
    ("PardisoSym", ["Pardiso", "-matrixType", 2]),
    ("PardisoSPD", ["Pardiso", "-matrixType", 1]),
]

base = None
print(f"model {nx}^3, {NSTEPS} steps, MKL_NUM_THREADS={os.environ['MKL_NUM_THREADS']}")
for name, args in CASES:
    print(f"--- {name} ...", flush=True)
    try:
        tip = build(args)
        t0 = time.perf_counter()
        for s in range(NSTEPS):
            if ops.analyze(1) != 0:
                raise RuntimeError(f"analyze failed at step {s + 1}")
        wall = time.perf_counter() - t0
        ux = ops.nodeDisp(tip, 1)
        if base is None:
            base = ux
        err = abs(ux - base) / abs(base) if base else float("nan")
        print(f"    {name:11s} {wall:8.3f} s   ux={ux:.12e}  rel_err={err:.3e}", flush=True)
    except Exception as e:  # noqa: BLE001
        print(f"    {name:11s} FAILED: {e}", flush=True)
print("SMOKE-DONE")
