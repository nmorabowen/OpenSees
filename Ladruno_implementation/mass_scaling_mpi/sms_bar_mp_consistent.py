"""Parallel CONSISTENT (Olovsson) SMS shared-node validation (OpenSeesMP).

Same 1-D fixed-free bar as sms_bar_mp.py (20 truss elements, fine 4-element zone
9..12 straddling the central node 11), but driven by CentralDifferenceSMSConsistent.
The consistent scaling injects per-element centroidal blocks M_bar_e = beta[diag(m)
- m m^T/M_e]; under a 2-rank split the elements either side of node 11 contribute
OFF-DIAGONAL coupling that crosses the partition. The step solves M_tilde a = r by
a matrix-free Jacobi-PCG; in parallel the PCG must (1) sum shared-DOF matvec entries
across ranks and (2) form GLOBAL inner products -- otherwise the boundary coupling is
dropped and np=2 diverges from the np=1 whole-model reference.

Validation: run np=1 (whole bar; the distributed PCG degenerates to the serial PCG
-- no neighbours, identity reduce) vs np=2 (split at node 11). Matching tip-disp
histories prove the distributed consistent PCG reproduces the reference solve.

Run:  mpiexec -n 1 python sms_bar_mp_consistent.py <dist\openseesmp>
      mpiexec -n 2 python sms_bar_mp_consistent.py <dist\openseesmp>
Writes tipc_np<NP>.out (time  tipDisp) next to this script.
"""
from __future__ import annotations
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MODDIR = sys.argv[1]  # dir holding openseesmp.pyd

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
sys.path.insert(0, MODDIR)

import openseesmp as ops  # noqa: E402

pid = ops.getPID()
npr = ops.getNP()
ops.start()

E, A, rho = 1.0e4, 1.0, 1.0
N = 20
hB, hF = 1.0, 0.1
dtTar, dtRun = 0.004, 0.0032
nSteps, Fapp = 150, 1.0
fLo, fHi = N // 2 - 1, N // 2 + 2
tipNode = N + 1

x = {1: 0.0}
for e in range(1, N + 1):
    h = hF if (fLo <= e <= fHi) else hB
    x[e + 1] = x[e] + h

ops.model("basic", "-ndm", 1, "-ndf", 1)
ops.uniaxialMaterial("Elastic", 1, E)

e0 = pid * N // npr + 1
e1 = (pid + 1) * N // npr
nLo, nHi = e0, e1 + 1

for n in range(nLo, nHi + 1):
    ops.node(n, x[n])
for e in range(e0, e1 + 1):
    ops.element("Truss", e, e, e + 1, A, 1, "-rho", rho)

if nLo <= 1 <= nHi:
    ops.fix(1, 1)

if nLo <= tipNode <= nHi:
    ops.timeSeries("Constant", 1, "-factor", 1.0)
    ops.pattern("Plain", 1, 1)
    ops.load(tipNode, Fapp)
    ops.recorder("Node", "-file", os.path.join(HERE, f"tipc_np{npr}.out"),
                 "-time", "-node", tipNode, "-dof", 1, "disp")

ops.constraints("Transformation")
ops.numberer("ParallelPlain")
ops.algorithm("Linear")
ops.system("MPIDiagonal")
ops.integrator("CentralDifferenceSMSConsistent", dtTar)
ops.analysis("Transient")
ops.analyze(nSteps, dtRun)

print(f"RANK {pid}/{npr} DONE tipNode={tipNode}")
ops.wipe()
