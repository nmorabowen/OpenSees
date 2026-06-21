"""Serial CONSISTENT-SMS reference (plain opensees.pyd, system Diagonal).

The SAME 1-D bar as sms_bar_mp_consistent.py, run serially through the
Zone-A-validated serial path: DiagonalSOE + the serial Ladruno::consistentPCG.
Its tip-disp history (tipc_serial.out) is the independent gold reference for the
MPI distributed PCG (tipc_np2.out) -- if they match, the distributed
consistentParPCG reproduces the trusted serial consistentPCG, not just itself.

Run:  python sms_bar_serial_consistent.py <dist\bin>
"""
from __future__ import annotations
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MODDIR = sys.argv[1]  # dir holding opensees.pyd (dist\bin)

_MKL = r"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin"
_ICOMP = r"C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin"
for d in (MODDIR, _MKL, _ICOMP):
    if d and os.path.isdir(d):
        os.add_dll_directory(d)
sys.path.insert(0, MODDIR)

import opensees as ops  # noqa: E402

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

ops.wipe()
ops.model("basic", "-ndm", 1, "-ndf", 1)
ops.uniaxialMaterial("Elastic", 1, E)
for n in range(1, N + 2):
    ops.node(n, x[n])
for e in range(1, N + 1):
    ops.element("Truss", e, e, e + 1, A, 1, "-rho", rho)
ops.fix(1, 1)
ops.timeSeries("Constant", 1, "-factor", 1.0)
ops.pattern("Plain", 1, 1)
ops.load(tipNode, Fapp)
ops.recorder("Node", "-file", os.path.join(HERE, "tipc_serial.out"),
             "-time", "-node", tipNode, "-dof", 1, "disp")
ops.constraints("Transformation")
ops.numberer("Plain")
ops.algorithm("Linear")
ops.system("Diagonal")
ops.integrator("CentralDifferenceSMSConsistent", dtTar)
ops.analysis("Transient")
ops.analyze(nSteps, dtRun)
ops.wipe()
print("SERIAL consistent reference DONE")
