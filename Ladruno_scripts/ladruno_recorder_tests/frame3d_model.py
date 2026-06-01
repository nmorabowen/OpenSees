"""Element-force parity on a 3-D FRAME — model runner (L1, element + nodal).

New element-type coverage (3-D beam-column end forces, 6-DOF nodes) for the
ladruno recorder. Runs ONE model with BOTH recorders, diffed to 1e-12 by the
generic element_parity_check.py (element channel) and parity_check.py (nodal).

    recorder mpco        -> fref.mpco      (frozen value oracle)
    recorder ladruno     -> ftest.ladruno  (the new recorder)

Mixes an elasticBeamColumn (3d) and a forceBeamColumn (3d, fiber section) so both
the closed-form and the integration-point force paths are exercised. Requests
`-E localForce` (12 local end-force components) + `-N displacement rotation`.

Run with the BUILD python, TMP = dir holding opensees.pyd + the MKL DLLs:
    python frame3d_model.py <dir_with_opensees_pyd_and_mkl> <out_dir>
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP

os.add_dll_directory(TMP)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402

ref = os.path.join(OUT, "fref.mpco")
new = os.path.join(OUT, "ftest.ladruno")
for p in (ref, new):
    if os.path.exists(p):
        os.remove(p)

ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 6)

# L-shaped cantilever: 1 -> 2 along X (elastic), 2 -> 3 up Z (force/fiber).
ops.node(1, 0.0, 0.0, 0.0)
ops.node(2, 2.0, 0.0, 0.0)
ops.node(3, 2.0, 0.0, 2.0)
ops.fix(1, 1, 1, 1, 1, 1, 1)

# elasticBeamColumn 3d on element 10 (vecxz not parallel to local x = global X)
E, G = 2.0e5, 8.0e4
A, Jx, Iy, Iz = 1.0e2, 1.0e3, 5.0e2, 5.0e2
ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)   # for member along X
ops.element("elasticBeamColumn", 10, 1, 2, A, E, G, Jx, Iy, Iz, 1)

# forceBeamColumn 3d on element 11 with a small fiber section
ops.uniaxialMaterial("Elastic", 1, E)
# -GJ supplies torsional stiffness so the 3-D fiber section is not singular
ops.section("Fiber", 1, "-GJ", G * Jx)
ops.patch("rect", 1, 2, 2, -0.5, -0.5, 0.5, 0.5)
ops.geomTransf("Linear", 2, 1.0, 0.0, 0.0)   # for member along Z
ops.beamIntegration("Lobatto", 1, 1, 3)
ops.element("forceBeamColumn", 11, 2, 3, 2, 1)

EREQ = ["-E", "localForce"]
ops.recorder("mpco", ref, *EREQ, "-N", "displacement", "rotation", "-T", "dt", 0.0)
ops.recorder("ladruno", new, *EREQ, "-N", "displacement", "rotation", "-T", "dt", 0.0)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(3, 1.0e2, 5.0e1, -2.0e2, 0.0, 0.0, 0.0)
ops.system("BandGeneral")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 0.2)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(5)
ops.wipe()

print("REF:", ref, os.path.exists(ref))
print("NEW:", new, os.path.exists(new))
print("FRAME3D_MODEL OK")
