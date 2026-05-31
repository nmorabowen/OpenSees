"""Parity gate — model runner (L1 of the test plan).

Runs ONE canonical model with BOTH recorders in a single process:
    recorder mpco        -> ref.mpco       (the frozen value oracle)
    recorder ladruno -> test.ladruno   (the new recorder)

so the emitted files can be value-diffed by parity_check.py to 1e-12.

Run with the BUILD python (no boot .pth), pointing sys.path at a copy of the
fresh OpenSeesPy.dll -> opensees.pyd (same two-step recipe as smoke_recorder.py):
    python parity_model.py <dir_with_opensees_pyd> <out_dir>
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

ref = os.path.join(OUT, "ref.mpco")
new = os.path.join(OUT, "test.ladruno")
for p in (ref, new):
    if os.path.exists(p):
        os.remove(p)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
# a small but non-trivial 2-D truss frame (stable, multi-node, real reactions)
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

# BOTH recorders, identical requests
NODE_REQ = ["-N", "displacement", "reactionForce"]
ops.recorder("mpco", ref, *NODE_REQ, "-T", "dt", 0.0)
ops.recorder("ladruno", new, *NODE_REQ, "-T", "dt", 0.0)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(4, 0.0, -1.0e3)
ops.system("BandSPD")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 0.2)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(5)
ops.wipe()  # flush + close both recorders

print("REF:", ref, os.path.exists(ref))
print("NEW:", new, os.path.exists(new))
print("PARITY_MODEL OK")
