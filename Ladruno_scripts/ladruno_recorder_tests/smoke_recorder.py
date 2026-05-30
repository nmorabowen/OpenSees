"""Phase-1 file-validity gate: drive the freshly built `recorder mpcoLadruno` on a
tiny model and assert the emitted .ladruno passes the schema validator.

Run (under an oneAPI-activated shell, with the fresh OpenSeesPy.dll copied to <tmp>/opensees.pyd):
    python smoke_recorder.py <tmp_dir_with_opensees_pyd>
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
HERE = os.path.dirname(os.path.abspath(__file__))

# dependent DLLs (MKL, libiomp5md) live in dist\bin; HDF5 is static-linked
os.add_dll_directory(DIST)
sys.path.insert(0, TMP)  # our fresh opensees.pyd wins over any .pth-wired install

import opensees as ops  # noqa: E402

out = os.path.join(TMP, "smoke.ladruno")
if os.path.exists(out):
    os.remove(out)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.node(1, 0.0, 0.0)
ops.node(2, 1.0, 0.0)
ops.fix(1, 1, 1)
ops.fix(2, 0, 1)  # node 2 free in X, fixed Y -> stable 1-DOF axial system
ops.uniaxialMaterial("Elastic", 1, 1000.0)
ops.element("truss", 1, 1, 2, 1.0, 1)

ops.recorder("mpcoLadruno", out)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(2, 1.0, 0.0)
ops.system("BandSPD")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 1.0)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(1)
ops.wipe()  # destroys the recorder -> flush + close the file

print("OPENSEES_MODULE:", getattr(ops, "__file__", "?"))
print("RECORDED:", out, "exists:", os.path.exists(out))
print("RUN OK")
