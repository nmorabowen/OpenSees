"""LOCAL_AXES + KIND gate -- 2D model runner.

Same idea as localaxes_model.py but for the planar (ndm=2,ndf=3) beam families
that now answer the Ladruno "localAxes" setResponse (id 30):

    ElasticBeam2d   (elasticBeamColumn)
    DispBeamColumn2d (dispBeamColumn)
    ForceBeamColumn2d (forceBeamColumn)

Each member runs along local x = (1,1)/sqrt(2) so a wrong identity frame fails
loudly. 2D CrdTransf::getLocalAxes returns full 3D direction cosines
(vx=(cos,sin,0), vy=(-sin,cos,0), vz=(0,0,1)) so the same quaternion check works;
the checker pads 2D node coordinates to 3D before comparing.

    recorder mpcoLadruno -> localaxes_2d.ladruno

Run with the BUILD python (no boot .pth):
    python localaxes_model_2d.py <dir_with_opensees_pyd> [out_dir]
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP

os.add_dll_directory(TMP)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402

new = os.path.join(OUT, "localaxes_2d.ladruno")
if os.path.exists(new):
    os.remove(new)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 3)

#   ElasticBeam2d : nodes 1-2
#   DispBeamColumn2d : nodes 3-4
#   ForceBeamColumn2d : nodes 5-6
ops.node(1, 0.0, 0.0)
ops.node(2, 1.0, 1.0)
ops.node(3, 0.0, 2.0)
ops.node(4, 1.0, 3.0)
ops.node(5, 0.0, 4.0)
ops.node(6, 1.0, 5.0)
for n in (1, 3, 5):
    ops.fix(n, 1, 1, 1)

ops.geomTransf("Linear", 1)

# ElasticBeam2d
# elasticBeamColumn tag n1 n2  A E Iz  transfTag
ops.element("elasticBeamColumn", 1, 1, 2, 1.0, 1000.0, 1.0, 1)

# section Elastic tag E A Iz
ops.section("Elastic", 1, 1000.0, 1.0, 1.0)
ops.beamIntegration("Legendre", 1, 1, 3)
ops.element("dispBeamColumn", 2, 3, 4, 1, 1)
ops.element("forceBeamColumn", 3, 5, 6, 1, 1)

ops.recorder("mpcoLadruno", new, "-N", "displacement", "-kind", "transient",
             "-T", "dt", 0.0)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(2, 1.0, 0.0, 0.0)
ops.load(4, 1.0, 0.0, 0.0)
ops.load(6, 1.0, 0.0, 0.0)
ops.system("BandGeneral")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 1.0)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(1)
ops.wipe()

print("NEW:", new, os.path.exists(new))
print("LOCALAXES_MODEL_2D OK")
