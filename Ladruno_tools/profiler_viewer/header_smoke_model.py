"""Run-header smoke -- model runner (BUILD python, no h5py needed).

Runs an explicit transient SDOF with CentralDifferenceLadruno (which overrides
getCriticalTimeStep) so the report command populates the run header: the
algorithm / integrator / solver name strings AND the explicit-dynamics dt_cr +
oversample lever (P0#5). The checker (venv python) reads the header back.

    python header_smoke_model.py <dir_with_opensees_pyd> <out_dir>
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP

os.add_dll_directory(TMP)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402

out_h5 = os.path.join(OUT, "header_smoke.h5")
if os.path.exists(out_h5):
    os.remove(out_h5)

DT = 0.005  # < dt_cr -> oversample > 1

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
N = 8  # chain of N trusses -> multi-DOF eigenproblem for dt_cr
for i in range(N + 1):
    ops.node(i + 1, float(i), 0.0)
ops.fix(1, 1, 1)
for i in range(2, N + 2):
    ops.fix(i, 0, 1)            # free in x, fixed in y
ops.uniaxialMaterial("Elastic", 1, 100.0)
for i in range(1, N + 1):
    ops.element("Truss", i, i, i + 1, 1.0, 1)
for i in range(2, N + 2):
    ops.mass(i, 1.0, 0.0)

ops.timeSeries("Constant", 1)
ops.pattern("Plain", 1, 1)
ops.load(N + 1, 1.0, 0.0)
ops.constraints("Plain")
ops.numberer("RCM")
ops.system("BandGeneral")
ops.test("NormDispIncr", 1e-8, 10, 0)
ops.algorithm("Linear")
ops.integrator("CentralDifferenceLadruno")   # explicit; overrides getCriticalTimeStep
ops.analysis("Transient")

ops.profiler("start", "-perStep")
ops.analyze(30, DT)
ops.profiler("stop")
ops.profiler("report", out_h5, "-run", "explicit")
ops.wipe()  # flush + close

print("DT", DT)
print("NEW:", out_h5, os.path.exists(out_h5))
print("HEADER_SMOKE_MODEL OK")
