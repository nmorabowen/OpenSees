"""P3 deep per-element-type smoke -- model runner (BUILD python, no h5py needed).

Runs a small static truss analysis with the profiler DEEP layer armed
(`profiler('start','-deep')`) so the per-element timing hooks fire inside
IncrementalIntegrator::formTangent / formElementResidual. Each element's
getTangent/getResidual is timed and folded into a per-classTag bucket
(elem_by_type) under the `elem.tangent` / `elem.residual` rollup nodes. The
checker (venv python) reads it back and asserts those buckets exist and carry
the Truss classTag with the expected element count.

    python deep_smoke_model.py <dir_with_opensees_pyd> <out_dir>
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP

os.add_dll_directory(TMP)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402

out_h5 = os.path.join(OUT, "deep_smoke.h5")
if os.path.exists(out_h5):
    os.remove(out_h5)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)

N = 20  # chain of N trusses -> N per-type tangent/residual samples
for i in range(N + 1):
    ops.node(i + 1, float(i), 0.0)
ops.fix(1, 1, 1)
for i in range(2, N + 2):
    ops.fix(i, 0, 1)
ops.uniaxialMaterial("Elastic", 1, 1000.0)
for i in range(1, N + 1):
    ops.element("Truss", i, i, i + 1, 1.0, 1)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(N + 1, 1.0, 0.0)
ops.system("BandGeneral")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 0.1)
ops.test("NormDispIncr", 1e-8, 10, 0)
ops.algorithm("Newton")
ops.analysis("Static")

# arm the DEEP layer BEFORE the analysis so the per-element hooks fire
ops.profiler("start", "-deep")
ops.analyze(10)
ops.profiler("stop")
ops.profiler("report", out_h5, "-run", "deepcase")
ops.wipe()  # flush + close

print("N_TRUSS", N)
print("NEW:", out_h5, os.path.exists(out_h5))
print("DEEP_SMOKE_MODEL OK")
