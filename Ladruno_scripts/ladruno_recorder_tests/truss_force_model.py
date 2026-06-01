"""Element-force parity on TRUSS + ZEROLENGTH — model runner (L1, element channel).

New element-type coverage (truss/zeroLength end forces) for the ladruno recorder.
Runs ONE model with BOTH recorders requesting element FORCE responses, so the
emitted files diff to 1e-12 via the generic element_parity_check.py.

    recorder mpco        -> tref.mpco      (frozen value oracle)
    recorder ladruno     -> ttest.ladruno  (the new recorder)

Requests `-E axialForce basicForces force`:
  - Truss::setResponse  "axialForce" (scalar) + "basicForces"
  - ZeroLength::setResponse "force" (the spring's material force vector)
Each recorder forwards the SAME token list to every element's setResponse; an
element that does not implement a token simply yields no bucket for it (mpco and
ladruno behave identically), so the intersection still diffs cleanly.

Run with the BUILD python, TMP = dir holding opensees.pyd + the MKL DLLs:
    python truss_force_model.py <dir_with_opensees_pyd_and_mkl> <out_dir>
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP

os.add_dll_directory(TMP)          # opensees.pyd + MKL co-located (dist/bin)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402

ref = os.path.join(OUT, "tref.mpco")
new = os.path.join(OUT, "ttest.ladruno")
for p in (ref, new):
    if os.path.exists(p):
        os.remove(p)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)

# A stable 2-D truss frame (nodes 1..4) + a grounded zeroLength spring at node 5.
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

# zeroLength spring: node 5 grounded to node 6 (coincident), X+Y springs.
ops.node(5, 0.0, 4.0)
ops.node(6, 0.0, 4.0)
ops.fix(6, 1, 1)
ops.uniaxialMaterial("Elastic", 2, 5.0e3)
ops.element("zeroLength", 20, 6, 5, "-mat", 2, 2, "-dir", 1, 2)

# BOTH recorders, identical element request.
EREQ = ["-E", "axialForce", "basicForces", "force"]
ops.recorder("mpco", ref, *EREQ, "-N", "displacement", "-T", "dt", 0.0)
ops.recorder("ladruno", new, *EREQ, "-N", "displacement", "-T", "dt", 0.0)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(4, 0.0, -1.0e3)
ops.load(5, 1.0e2, -2.0e2)
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
print("TRUSS_FORCE_MODEL OK")
