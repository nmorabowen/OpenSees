"""ADR-40 rank-3 MUMPS parallel strong-scaling measurement (OpenSeesMP).

Implicit-static 3D solid, FIXED global mesh (nx=14, ny=14, nz=28 -> 5488
LadrunoBrick elements, ~20k DOF), vary ranks 1/2/4/8. Partition = z-slabs by pid.

Launch (per the proven template):
  dist/openseesmp/mpiexec.exe -n N  python -S  mumps_scaling.py  dist/openseesmp
Rank 0 prints RESULT np=<N> wall=<s> steps=<done>/8 and the global top-center
tip displacement (17 sig figs). Every rank prints its partition-local sum|u|.
"""
from __future__ import annotations
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
MODDIR = sys.argv[1]  # dir holding openseesmp.pyd

# --- DLL wiring (verbatim from sms_bar_mp.py) ------------------------------
_IMPI = r"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\bin"
_LIBFABRIC = r"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\opt\mpi\libfabric\bin"
_MKL = r"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin"
_ICOMP = r"C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin"
_DISTBIN = os.path.join(os.path.dirname(MODDIR.rstrip("\\/")), "bin")  # sibling dist\bin
for d in (MODDIR, _DISTBIN, _IMPI, _LIBFABRIC, _MKL, _ICOMP):
    if d and os.path.isdir(d):
        os.add_dll_directory(d)
os.environ["PATH"] = os.pathsep.join(
    [p for p in (_IMPI, _LIBFABRIC, MODDIR, os.environ.get("PATH", "")) if p])
sys.path.insert(0, MODDIR)

import openseesmp as ops  # noqa: E402

pid = ops.getPID()
npr = ops.getNP()
ops.start()

# --- geometry / material parameters ----------------------------------------
NX, NY, NZ = 14, 14, 28
H = 100.0                       # element edge (mm)
E, NU = 200.0e3, 0.3            # 200 GPa (N/mm^2), nu=0.3
K = E / (3.0 * (1.0 - 2.0 * NU))   # 166666.6667
G = E / (2.0 * (1.0 + NU))         # 76923.0769
MATTAG = 1
FX, FZ = 4.0e5, -4.0e5          # per top-face node

# node numbering: nid(i,j,k) = 1 + i + (NX+1)*j + (NX+1)*(NY+1)*k
NXP, NYP = NX + 1, NY + 1


def nid(i, j, k):
    return 1 + i + NXP * j + NXP * NYP * k


def xyz(i, j, k):
    return (i * H, j * H, k * H)


# --- z-slab ownership: rank owns elements with k in [k0, k1) ---------------
k0 = pid * NZ // npr
k1 = (pid + 1) * NZ // npr

ops.model("basic", "-ndm", 3, "-ndf", 3)
ops.nDMaterial("LadrunoJ2", MATTAG, K, G, "-iso", "voce", 250.0, 0.0, 0.0, 2000.0)

# define every node touched by this rank's owned element layers.
# element layer k touches node planes k and k+1 -> node planes [k0, k1] inclusive.
my_nodes = set()
for k in range(k0, k1 + 1):
    for j in range(NYP):
        for i in range(NXP):
            n = nid(i, j, k)
            ops.node(n, *xyz(i, j, k))
            my_nodes.add(n)

# elements owned by this rank (element tag encodes the ordinal)
for k in range(k0, k1):
    for j in range(NY):
        for i in range(NX):
            etag = 1 + i + NX * j + NX * NY * k
            # standard hex node order: bottom face CCW then top face CCW
            n1 = nid(i, j, k)
            n2 = nid(i + 1, j, k)
            n3 = nid(i + 1, j + 1, k)
            n4 = nid(i, j + 1, k)
            n5 = nid(i, j, k + 1)
            n6 = nid(i + 1, j, k + 1)
            n7 = nid(i + 1, j + 1, k + 1)
            n8 = nid(i, j + 1, k + 1)
            ops.element("LadrunoBrick", etag, n1, n2, n3, n4, n5, n6, n7, n8, MATTAG)

# --- base plane k=0 fully fixed (only ranks that own those nodes) ----------
if k0 == 0:
    for j in range(NYP):
        for i in range(NXP):
            ops.fix(nid(i, j, 0), 1, 1, 1)

# --- top-face k=NZ nodes loaded (only rank owning that plane) --------------
top_center = nid(NX // 2, NY // 2, NZ)   # global tip node
if k1 == NZ:
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(NYP):
        for i in range(NXP):
            ops.load(nid(i, j, NZ), FX, 0.0, FZ)

# --- analysis --------------------------------------------------------------
ops.constraints("Plain")
numberer_used = "ParallelRCM"
try:
    ops.numberer("ParallelRCM")
except Exception:
    ops.numberer("ParallelPlain")
    numberer_used = "ParallelPlain(fallback)"

# -ICNTL14 = MUMPS workspace relaxation %. Default 20 gives Error -9 (work array
# too small) on the partitioned (np>=4) factorization of this ~19k-DOF solid;
# 200 gives ample headroom and is applied uniformly at every np for a fair wall.
_icntl14 = int(os.environ.get("MUMPS_ICNTL14", "200"))
ops.system("Mumps", "-ICNTL14", _icntl14)
ops.test("NormDispIncr", 1.0e-7, 25)
ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0 / 8.0)
ops.analysis("Static")

# optional profiler on rank 0 only
prof_on = False
if pid == 0:
    try:
        ops.profiler("start")
        prof_on = True
    except Exception:
        prof_on = False

done = 0
t0 = time.perf_counter()
for step in range(8):
    ok = ops.analyze(1)
    if ok != 0:
        break
    done += 1
wall = time.perf_counter() - t0

if pid == 0 and prof_on:
    try:
        ops.profiler("stop")
        ops.profiler("report", os.path.join(HERE, f"mumps_prof_np{npr}_r0.h5"))
    except Exception as e:
        print(f"PROFILER-REPORT-FAILED: {e}")

# --- fingerprints ----------------------------------------------------------
local_sum = 0.0
for n in my_nodes:
    for dof in (1, 2, 3):
        local_sum += abs(ops.nodeDisp(n, dof))
print(f"FINGERPRINT pid={pid} np={npr} sum|u|={local_sum:.10e} nnodes={len(my_nodes)}")

# tip displacement from whichever rank owns the global top-center node
if top_center in my_nodes:
    ux = ops.nodeDisp(top_center, 1)
    uy = ops.nodeDisp(top_center, 2)
    uz = ops.nodeDisp(top_center, 3)
    print(f"TIPDISP np={npr} node={top_center} "
          f"ux={ux:.17e} uy={uy:.17e} uz={uz:.17e}")

if pid == 0:
    print(f"RESULT np={npr} wall={wall:.4f} steps={done}/8 "
          f"numberer={numberer_used} prof={prof_on}")

ops.wipe()
