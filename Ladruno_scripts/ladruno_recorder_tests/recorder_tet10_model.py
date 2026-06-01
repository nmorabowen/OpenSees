"""Recorder round-trip — BezierTet10 model runner.

Drives `recorder ladruno` on a BezierTet10 model and writes a .ladruno
file. Verifies the inverted-dispatch wiring: the recorder should self-declare
the element basis from its `basisInfo`/`integrationPoints`/`integrationWeights`
responses (FAMILY=bernstein, topology=tet, NUM_CTRL=10, NUM_GP=4, exact
GP_PARAM/GP_WEIGHT) with no recorder class-specific code.

Run with the BUILD python (no boot .pth), sys.path at the fresh opensees.pyd:
    python recorder_tet10_model.py <dir_with_opensees_pyd> <out_dir>
Check the emitted file with recorder_tet10_check.py (venv python + h5py).
"""
import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP
DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\.claude\worktrees\keen-mestorf-50bccd\dist\bin"

os.add_dll_directory(DIST)
sys.path.insert(0, TMP)
import opensees as ops  # noqa: E402

out = os.path.join(OUT, "test_tet10.ladruno")
if os.path.exists(out):
    os.remove(out)

# straight-sided reference tet, node order = TenNodeTetrahedron (D2′)
V = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 1, 0), 4: (0, 0, 1)}
EDGES = [(1, 2), (2, 3), (1, 3), (1, 4), (3, 4), (2, 4)]
NODES = dict(V)
for i, (a, b) in enumerate(EDGES):
    NODES[5 + i] = tuple(0.5 * (V[a][k] + V[b][k]) for k in range(3))

ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)
ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
for n, (x, y, z) in NODES.items():
    ops.node(n, float(x), float(y), float(z))
ops.element("BezierTet10", 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1)

# fix the three base vertices fully; load the apex
for n in (1, 2, 3):
    ops.fix(n, 1, 1, 1)

ops.recorder("ladruno", out, "-N", "displacement", "-E", "stresses",
             "-T", "dt", 0.0)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(4, 0.0, 0.0, -1.0)
ops.system("FullGeneral")
ops.numberer("Plain")
ops.constraints("Plain")
ops.integrator("LoadControl", 0.5)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(2)
ops.wipe()  # flush + close recorder

print("NEW:", out, os.path.exists(out))
print("RECORDER_TET10_MODEL OK")
