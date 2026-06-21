"""Parallel lumped-SMS shared-node validation (OpenSeesMP / openseespy).

Mirrors sms_bar_mp.tcl: a 1-D fixed-free truss bar of 20 elements with a fine
4-element zone (elements 9..12) straddling the central node 11. Under a 2-rank
manual partition the bar splits at node 11 (SHARED). The two fine elements either
side of it -- element 10 (rank 0) and element 11 (rank 1) -- are both below
dtTarget, so CentralDifferenceSMS injects fictitious mass into node 11 from BOTH
ranks. Correct behaviour requires the MPIDiagonal solver to SUM those per-rank
contributions across the partition.

Run:  mpiexec -n 1 python sms_bar_mp.py <dist\bin>
      mpiexec -n 2 python sms_bar_mp.py <dist\bin>
Writes tip_np<NP>.out (time  tipDisp) next to this script.
"""
from __future__ import annotations
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MODDIR = sys.argv[1]  # dir holding openseesmp.pyd

# openseesmp.pyd links Intel MPI / libfabric / MKL / ifx runtime; register every
# oneAPI bin that exists (DLLs resolve from add_dll_directory + a process PATH for
# libfabric). Matches Ladruno_scripts/ladruno_recorder_tests/mp_parallel_model.py.
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

# ---- parameters (identical to the Tcl version) ---------------------------
E, A, rho = 1.0e4, 1.0, 1.0
N = 20
hB, hF = 1.0, 0.1
dtTar, dtRun = 0.004, 0.0032
nSteps, Fapp = 150, 1.0
fLo, fHi = N // 2 - 1, N // 2 + 2          # fine elements 9..12
tipNode = N + 1

# global node x-coordinates (same on every rank)
x = {1: 0.0}
for e in range(1, N + 1):
    h = hF if (fLo <= e <= fHi) else hB
    x[e + 1] = x[e] + h

ops.model("basic", "-ndm", 1, "-ndf", 1)
ops.uniaxialMaterial("Elastic", 1, E)

# ownership: rank r owns elements [r*N/np+1 .. (r+1)*N/np]; touches nodes e0..e1+1
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
    ops.recorder("Node", "-file", os.path.join(HERE, f"tip_np{npr}.out"),
                 "-time", "-node", tipNode, "-dof", 1, "disp")

ops.constraints("Transformation")
ops.numberer("ParallelPlain")
ops.algorithm("Linear")
ops.system("MPIDiagonal")
ops.integrator("CentralDifferenceSMS", dtTar)
ops.analysis("Transient")
ops.analyze(nSteps, dtRun)

print(f"RANK {pid}/{npr} DONE tipNode={tipNode}")
ops.wipe()  # flush the recorder
