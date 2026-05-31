"""P4 memory-counter smoke -- model runner (BUILD python, no h5py needed).

Runs a small static truss analysis with the profiler memory counters armed
(`profiler('start','-memory')`) so the Matrix/Vector alloc hooks fire during
formTangent/solve, then writes a profile run via the in-engine HDF5 writer. The
checker (venv python) reads it back and asserts matrix_live/vector_live/peak are
populated and self-consistent.

    python mem_smoke_model.py <dir_with_opensees_pyd> <out_dir>
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP

os.add_dll_directory(TMP)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402

out_h5 = os.path.join(OUT, "mem_smoke.h5")
if os.path.exists(out_h5):
    os.remove(out_h5)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)

# Arm the memory counters BEFORE building the domain. The Matrix/Vector/ID byte
# counters fire on allocations during analyze() either way, but the TaggedObject
# census (MovableObject ctor/dtor) only counts components BORN while mem() is on,
# so the Node/Truss/ElasticMaterial objects must be constructed under the armed
# window for components_live to be non-empty at report time.
ops.profiler("start", "-memory")

N = 20  # chain of N trusses -> a real K/residual to allocate
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

ops.analyze(10)
peak = ops.profiler("memory")        # -> peak bytes (int)
ops.profiler("stop")
ops.profiler("report", out_h5, "-run", "memcase")
ops.wipe()  # flush + close

print("PEAK_BYTES", peak)
print("NEW:", out_h5, os.path.exists(out_h5))
print("MEM_SMOKE_MODEL OK")
