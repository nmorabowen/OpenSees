"""P0#3 per-step series smoke -- model runner (BUILD python, no h5py needed).

Runs a small transient SDOF analysis with the per-step series armed
(`profiler('start','-perStep')`) so the analysis driver's OPS_PROFILE_STEP hook
fires once per committed step, recording step/t/dt/iters. The checker (venv
python) reads the /series group back and asserts it has one row per step with a
monotonic step index, the expected dt, and a positive Newton iteration count.

    python series_smoke_model.py <dir_with_opensees_pyd> <out_dir>
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP

os.add_dll_directory(TMP)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402

out_h5 = os.path.join(OUT, "series_smoke.h5")
if os.path.exists(out_h5):
    os.remove(out_h5)

N = 20
DT = 0.01

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.node(1, 0.0, 0.0)
ops.node(2, 1.0, 0.0)
ops.fix(1, 1, 1)
ops.fix(2, 0, 1)             # node 2: free in x only -> 1 DOF
ops.uniaxialMaterial("Elastic", 1, 100.0)
ops.element("Truss", 2, 1, 2, 1.0, 1)
ops.mass(2, 1.0, 0.0)        # x mass -> SDOF oscillator (omega = sqrt(k/m) = 10 rad/s)

ops.timeSeries("Constant", 1)
ops.pattern("Plain", 1, 1)
ops.load(2, 1.0, 0.0)
ops.constraints("Plain")
ops.numberer("RCM")
ops.system("BandGeneral")
ops.test("NormDispIncr", 1e-10, 10, 0)
ops.algorithm("Newton")
ops.integrator("Newmark", 0.5, 0.25)
ops.analysis("Transient")

ops.profiler("start", "-perStep")
ops.analyze(N, DT)
ops.profiler("stop")
ops.profiler("report", out_h5, "-run", "transient")
ops.wipe()  # flush + close

print("N_STEPS", N, "DT", DT)
print("NEW:", out_h5, os.path.exists(out_h5))
print("SERIES_SMOKE_MODEL OK")
