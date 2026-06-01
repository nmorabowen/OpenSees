"""Element-result parity on SHELLS (ShellMITC4) — model runner (L1).

New element-type coverage (shell stress resultants, 6-DOF nodes) for the ladruno
recorder. Runs ONE model with BOTH recorders, diffed to 1e-12 by the generic
element_parity_check.py (element channel) + parity_check.py (nodal).

    recorder mpco        -> shref.mpco       (frozen value oracle)
    recorder ladruno     -> shtest.ladruno   (the new recorder)

Model: a 2x2 grid of ShellMITC4 over an ElasticMembranePlateSection, cantilevered
along one edge and tip-loaded out of plane. Requests `-E forces` (the 8 shell
stress resultants per integration point) + `-N displacement rotation`.

Run with the BUILD python, TMP = dir holding opensees.pyd + the MKL DLLs:
    python shell_model.py <dir_with_opensees_pyd_and_mkl> <out_dir>
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP

os.add_dll_directory(TMP)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402

ref = os.path.join(OUT, "shref.mpco")
new = os.path.join(OUT, "shtest.ladruno")
for p in (ref, new):
    if os.path.exists(p):
        os.remove(p)

ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 6)

# 3x3 node grid -> 2x2 ShellMITC4 in the z=0 plane.
nid = 0
for j in range(3):
    for i in range(3):
        nid += 1
        ops.node(nid, float(i), float(j), 0.0)
# clamp the j=0 edge (nodes 1,2,3)
for n in (1, 2, 3):
    ops.fix(n, 1, 1, 1, 1, 1, 1)

# ElasticMembranePlateSection: secTag E nu h rho
ops.section("ElasticMembranePlateSection", 1, 2.0e5, 0.25, 0.1, 0.0)

# 2x2 quad of shells, CCW
def n(i, j):
    return j * 3 + i + 1

etag = 100
for j in range(2):
    for i in range(2):
        etag += 1
        ops.element("ShellMITC4", etag, n(i, j), n(i + 1, j), n(i + 1, j + 1), n(i, j + 1), 1)

EREQ = ["-E", "forces"]
ops.recorder("mpco", ref, *EREQ, "-N", "displacement", "rotation", "-T", "dt", 0.0)
ops.recorder("ladruno", new, *EREQ, "-N", "displacement", "rotation", "-T", "dt", 0.0)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
# out-of-plane tip load on the free edge (nodes 7,8,9)
for nn in (7, 8, 9):
    ops.load(nn, 0.0, 0.0, -1.0e1, 0.0, 0.0, 0.0)
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
print("SHELL_MODEL OK")
