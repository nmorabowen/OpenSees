"""LOCAL_AXES + KIND gate -- model runner.

Two ElasticBeam3d members oriented along the global diagonal (1,1,0) so their
local frame is clearly NON-identity -- the whole point of MODEL/LOCAL_AXES is to
kill the silent identity-quaternion fallback, so the test geometry must make a
wrong identity frame fail loudly. With `-kind transient` we also exercise the
KIND stage attribute.

    recorder mpcoLadruno -> localaxes.ladruno

Run with the BUILD python (no boot .pth):
    python localaxes_model.py <dir_with_opensees_pyd> [out_dir]
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP

os.add_dll_directory(TMP)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402

new = os.path.join(OUT, "localaxes.ladruno")
if os.path.exists(new):
    os.remove(new)

ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 6)
ops.node(1, 0.0, 0.0, 0.0)
ops.node(2, 1.0, 1.0, 0.0)   # local x = (1,1,0)/sqrt(2)
ops.node(3, 2.0, 2.0, 0.0)
ops.fix(1, 1, 1, 1, 1, 1, 1)

# vecxz = (0,0,1): local z lies in the x-(0,0,1) plane -> frame is non-identity.
ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
# elasticBeamColumn tag n1 n2  A E G J Iy Iz  transfTag
ops.element("elasticBeamColumn", 1, 1, 2, 1.0, 1000.0, 400.0, 1.0, 1.0, 1.0, 1)
ops.element("elasticBeamColumn", 2, 2, 3, 1.0, 1000.0, 400.0, 1.0, 1.0, 1.0, 1)

ops.recorder("mpcoLadruno", new, "-N", "displacement", "-kind", "transient",
             "-T", "dt", 0.0)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(3, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0)
ops.system("BandSPD")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 1.0)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(1)
ops.wipe()

print("NEW:", new, os.path.exists(new))
print("LOCALAXES_MODEL OK")
