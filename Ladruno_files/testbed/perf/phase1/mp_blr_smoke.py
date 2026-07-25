"""ADR-75 P2 — MUMPS Block Low-Rank (-BLR) smoke test under MPI.

Verifies that `system Mumps -BLR <eps>` (ICNTL(35)=1 + CNTL(7)=eps) actually
runs on a partitioned model and lands on the SAME answer as the stock full-rank
factorization, to within the BLR dropping tolerance. BLR is an APPROXIMATE
factorization, so unlike every other ADR-75 gate this one is NOT expected to be
bit-identical -- the check is "agrees to ~eps", and a bit-identical result at a
loose eps would actually mean BLR never engaged.

Launch (proven template, from phase0/mumps_scaling.py):
  dist/openseesmp/mpiexec.exe -n 2 python -S mp_blr_smoke.py dist/openseesmp
"""
from __future__ import annotations
import os
import sys
import time

MODDIR = sys.argv[1]

_IMPI = r"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\bin"
_LIBFABRIC = r"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\opt\mpi\libfabric\bin"
_MKL = r"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin"
_ICOMP = r"C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin"
_DISTBIN = os.path.join(os.path.dirname(MODDIR.rstrip("\\/")), "bin")
for d in (MODDIR, _DISTBIN, _IMPI, _LIBFABRIC, _MKL, _ICOMP):
    if d and os.path.isdir(d):
        os.add_dll_directory(d)
os.environ["PATH"] = os.pathsep.join(
    [p for p in (_IMPI, _LIBFABRIC, MODDIR, os.environ.get("PATH", "")) if p])
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
sys.path.insert(0, MODDIR)

import openseesmp as ops  # noqa: E402

pid = ops.getPID()
npr = ops.getNP()
ops.start()

NX = int(os.environ.get("BLR_NX", "8")); NY = int(os.environ.get("BLR_NY", "8")); NZ = int(os.environ.get("BLR_NZ", "16"))
H = 100.0
E, NU = 200.0e3, 0.3
K = E / (3.0 * (1.0 - 2.0 * NU))
G = E / (2.0 * (1.0 + NU))
NXP, NYP = NX + 1, NY + 1


def nid(i, j, k):
    return 1 + i + NXP * j + NXP * NYP * k


def run(system_args, label):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", 250.0, 0.0, 0.0, 2000.0)
    k0 = pid * NZ // npr
    k1 = (pid + 1) * NZ // npr
    mine = set()
    for k in range(k0, k1 + 1):
        for j in range(NYP):
            for i in range(NXP):
                ops.node(nid(i, j, k), i * H, j * H, k * H)
                mine.add(nid(i, j, k))
    for k in range(k0, k1):
        for j in range(NY):
            for i in range(NX):
                e = 1 + i + NX * j + NX * NY * k
                ops.element("LadrunoBrick", e,
                            nid(i, j, k), nid(i+1, j, k), nid(i+1, j+1, k), nid(i, j+1, k),
                            nid(i, j, k+1), nid(i+1, j, k+1), nid(i+1, j+1, k+1), nid(i, j+1, k+1),
                            1)
    if k0 == 0:
        for j in range(NYP):
            for i in range(NXP):
                ops.fix(nid(i, j, 0), 1, 1, 1)
    tip = nid(NX // 2, NY // 2, NZ)
    if k1 == NZ:
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for j in range(NYP):
            for i in range(NXP):
                ops.load(nid(i, j, NZ), 4.0e5, 0.0, -4.0e5)
    ops.constraints("Plain")
    try:
        ops.numberer("LadrunoParallelRCM")
    except Exception:
        ops.numberer("ParallelRCM")
    ops.system(*system_args)
    ops.test("NormDispIncr", 1.0e-7, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / 4.0)
    ops.analysis("Static")
    t0 = time.perf_counter()
    ok = 0
    for _ in range(2):
        ok = ops.analyze(1)
        if ok != 0:
            break
    wall = time.perf_counter() - t0
    ux = ops.nodeDisp(tip, 1) if tip in mine else None
    if pid == 0:
        print(f"{label:28s} ok={ok} wall={wall:7.3f}s", flush=True)
    return ok, ux, wall


base = ["Mumps", "-ICNTL14", 200]
r_full = run(base, "Mumps (full-rank)")
r_blr = run(base + ["-BLR", 1.0e-9], "Mumps -BLR 1e-9")
r_blr_loose = run(base + ["-BLR", 1.0e-4], "Mumps -BLR 1e-4 (loose)")

for label, r in (("full", r_full), ("blr1e-9", r_blr), ("blr1e-4", r_blr_loose)):
    if r[1] is not None:
        print(f"TIP pid={pid} {label:8s} ux={r[1]:.17e}", flush=True)

if r_full[1] is not None:
    u0 = r_full[1]
    for label, r in (("blr1e-9", r_blr), ("blr1e-4", r_blr_loose)):
        if r[1] is not None and u0:
            print(f"AGREEMENT {label}: rel_diff={abs(r[1]-u0)/abs(u0):.3e}", flush=True)

if pid == 0:
    print(f"RESULT np={npr} all_ok={r_full[0]==0 and r_blr[0]==0 and r_blr_loose[0]==0}",
          flush=True)
ops.wipe()
