"""Multi-stage regression — model runner.

Drives ONE process through TWO model stages with BOTH recorders:
    recorder mpco        -> ms_ref.mpco       (frozen value oracle)
    recorder mpcoLadruno -> ms_test.ladruno   (the new recorder)

A topology change (a new node + truss element added between the two analysis
runs) bumps Domain::hasDomainChanged(), so the second analyze() lands in a NEW
MODEL_STAGE. This exercises the fix for the C1/C2 multi-stage bug: record() must
re-detect the stamp, re-write the model, and rebuild every source/sink (otherwise
only the first MODEL_STAGE is written and cached Element*/Response* dangle).

Both recorders see the identical model, so the per-stage nodal values must match
to 1e-12 (checked by multistage_check.py), AND both files must contain the SAME
set of MODEL_STAGE groups.

Run with the BUILD python (no boot .pth), same recipe as parity_model.py:
    python multistage_model.py <dir_with_opensees_pyd> <out_dir>
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]                      # dir holding the fresh opensees.pyd
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP
DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"

os.add_dll_directory(DIST)             # MKL etc. (HDF5 is static-linked)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402

ref = os.path.join(OUT, "ms_ref.mpco")
new = os.path.join(OUT, "ms_test.ladruno")
for p in (ref, new):
    if os.path.exists(p):
        os.remove(p)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)

# --- Stage 1 topology: a small 2-D truss frame ---------------------------
ops.node(1, 0.0, 0.0)
ops.node(2, 2.0, 0.0)
ops.node(3, 4.0, 0.0)
ops.node(4, 2.0, 2.0)
ops.fix(1, 1, 1)
ops.fix(3, 1, 1)
ops.uniaxialMaterial("Elastic", 1, 2.0e5)
ops.element("truss", 1, 1, 2, 1.0, 1)
ops.element("truss", 2, 2, 3, 1.0, 1)
ops.element("truss", 3, 1, 4, 1.0, 1)
ops.element("truss", 4, 2, 4, 1.0, 1)
ops.element("truss", 5, 3, 4, 1.0, 1)

# BOTH recorders, identical requests, record every step
NODE_REQ = ["-N", "displacement", "reactionForce"]
ops.recorder("mpco", ref, *NODE_REQ, "-T", "dt", 0.0)
ops.recorder("mpcoLadruno", new, *NODE_REQ, "-T", "dt", 0.0)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(4, 0.0, -1.0e3)
ops.system("BandSPD")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 0.2)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(3)                         # STAGE 1: 3 steps

# --- Topology change -> a new MODEL_STAGE ---------------------------------
# Add node 5 (fixed) above node 4 and a truss tying 4->5. This changes the
# domain geometry (Domain::hasDomainChanged() bumps), keeps the system non-
# singular, and adds node 5 to the recorded set on the rebuild.
ops.node(5, 2.0, 4.0)
ops.fix(5, 1, 1)
ops.element("truss", 6, 4, 5, 1.0, 1)

ops.analyze(3)                         # STAGE 2: 3 more steps (new stamp)
ops.wipe()                             # flush + close both recorders

print("REF:", ref, os.path.exists(ref))
print("NEW:", new, os.path.exists(new))
print("MULTISTAGE_MODEL OK")
