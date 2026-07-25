"""ADR-75 P1d — does symmetric PARDISO actually save MEMORY, not just time?

P1c established that on this fork's regime (huge 3D solid nonlinear models)
MEMORY, not time, is the binding constraint: UmfPack ran OUT OF MEMORY at 86,490
DOF while PARDISO solved it. So the symmetric path's memory claim has to be
measured, not inferred from "half the entries".

P2b's MUMPS BLR result is the cautionary precedent: the stored factors shrank
21.8% while PEAK memory moved only 4.6%, because peak is dominated by the active
frontal/working space rather than the stored factors. `-stats` reports BOTH here
(iparm[14]/[15]/[16] + nnz in factors) so the same distinction is visible.

One factorization per configuration is enough — the counters are determined by
the sparsity pattern, not by the Newton history.

Run:
  OPS_PYD=<worktree>\\dist\\bin python3.12 -S p1d_memory.py
Optional: P1D_N (default 15 => the Lane-B model, ~11.5k DOF).
"""
from __future__ import annotations
import os
import sys

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
print("pyd OK:", ops.__file__)

E, nu = 200000.0, 0.3
K = E / (3.0 * (1.0 - 2.0 * nu))
G = E / (2.0 * (1.0 + nu))
s0, Qinf, b, Hiso = 250.0, 0.0, 0.0, 2000.0
NX = int(os.environ.get("P1D_N", "15"))
nx = ny = nz = NX
L = 100.0


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
    ops.integrator("LoadControl", 1.0 / 15)
    ops.analysis("Static")


CASES = [
    ("unsymmetric  (mtype 11)", ["Pardiso", "-stats"]),
    ("symmetric    (mtype -2)", ["Pardiso", "-matrixType", 2, "-stats"]),
    ("SPD          (mtype  2)", ["Pardiso", "-matrixType", 1, "-stats"]),
]

print(f"model {nx}^3 LadrunoBrick+LadrunoJ2, ONE load step, "
      f"MKL_NUM_THREADS={os.environ['MKL_NUM_THREADS']}\n")
for label, args in CASES:
    print(f"===== {label} =====")
    build(args)
    rc = ops.analyze(1)
    print(f"  analyze -> {rc}\n")
print("MEMORY-DONE")
