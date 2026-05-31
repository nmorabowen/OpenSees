"""BezierTri6 recorder gate -- model runner.

Exercises the new self-describing-basis path: a 6-node quadratic Bezier
triangle (ELE_TAG_BezierTri6) is a CUSTOM-rule element whose geometry the
recorder cannot guess from a class-tag table. It must instead consume the
element's "basisInfo" / "integrationPoints" / "integrationWeights" probes.

Two BezierTri6 triangles tile a 2x2 square (straight-sided, v1), sharing the
diagonal edge and its mid-edge node (node 6). Plane stress, ElasticIsotropic.

Runs BOTH recorders so the element-stress channel can be value-diffed:
    recorder mpco        -> bezier_ref.mpco      (frozen value oracle)
    recorder ladruno -> bezier_test.ladruno  (the new recorder)

Run with the BUILD python (no boot .pth), same recipe as element_parity_model:
    python bezier_model.py <dir_with_opensees_pyd> [out_dir]
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP

os.add_dll_directory(TMP)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402

ref = os.path.join(OUT, "bezier_ref.mpco")
new = os.path.join(OUT, "bezier_test.ladruno")
for p in (ref, new):
    if os.path.exists(p):
        os.remove(p)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)

# Two BezierTri6 triangles over a 2x2 square, split on the (2,0)-(0,2) diagonal.
#   corners:  1=(0,0) 2=(2,0) 3=(0,2) 4=(2,2)
#   mid-edge: 5=(1,0)[1-2]  6=(1,1)[2-3, SHARED]  7=(0,1)[3-1]
#             8=(2,1)[2-4]  9=(1,2)[4-3]
ops.node(1, 0.0, 0.0)
ops.node(2, 2.0, 0.0)
ops.node(3, 0.0, 2.0)
ops.node(4, 2.0, 2.0)
ops.node(5, 1.0, 0.0)
ops.node(6, 1.0, 1.0)
ops.node(7, 0.0, 1.0)
ops.node(8, 2.0, 1.0)
ops.node(9, 1.0, 2.0)

# Clamp the x=0 edge (nodes 1,3,7) fully -> a stable cantilever.
ops.fix(1, 1, 1)
ops.fix(3, 1, 1)
ops.fix(7, 1, 1)

ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)

# element BezierTri6 tag c1 c2 c3 mid12 mid23 mid31 thick type matTag
# Triangle A: corners (1,2,3), mids (5,6,7)
ops.element("BezierTri6", 1, 1, 2, 3, 5, 6, 7, 1.0, "PlaneStress", 1)
# Triangle B: corners (2,4,3), mids (8,9,6)  -> mid(3-2)=6 shared
ops.element("BezierTri6", 2, 2, 4, 3, 8, 9, 6, 1.0, "PlaneStress", 1)

# BOTH recorders. Element stress (3 GP x 3 comps = 9 cols) + nodal disp.
ops.recorder("mpco", ref, "-N", "displacement", "-E", "stress", "-T", "dt", 0.0)
ops.recorder("ladruno", new, "-N", "displacement", "-E", "stress", "-T", "dt", 0.0)

# Pull the x=2 edge (nodes 2,4,8) in +x.
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(2, 1.0e2, 0.0)
ops.load(4, 1.0e2, 0.0)
ops.load(8, 1.0e2, 0.0)

ops.system("BandSPD")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 0.25)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(4)
ops.wipe()  # flush + close both recorders

print("REF:", ref, os.path.exists(ref))
print("NEW:", new, os.path.exists(new))
print("BEZIER_MODEL OK")
