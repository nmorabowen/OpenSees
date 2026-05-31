"""Element-result parity — model runner (L1, element channel).

Runs ONE canonical model with BOTH recorders, requesting ELEMENT results:
    recorder mpco        -> eref.mpco       (frozen value oracle)
    recorder ladruno -> etest.ladruno   (the new recorder)

Element results are the ported-verbatim-but-value-unverified channel: the frozen
recorder and ladruno share the same element engine and packing order, so the
per-element result row Vector must be identical -> the matrices diff to 1e-12.

Model: a 2-element FourNodeQuad plane-stress strip (4 Gauss pts x 3 stress comps
= 12 columns), plus the truss frame so nodal results are also present. Requests
`-E stress` (FourNodeQuad::setResponse "stress", FourNodeQuad.cpp:1370).

Run with the BUILD python (no boot .pth), same recipe as parity_model.py:
    python element_parity_model.py <dir_with_opensees_pyd> <out_dir>
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP
DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"

os.add_dll_directory(DIST)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402

ref = os.path.join(OUT, "eref.mpco")
new = os.path.join(OUT, "etest.ladruno")
for p in (ref, new):
    if os.path.exists(p):
        os.remove(p)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)

# 2x1 strip of FourNodeQuad (nodes 1..6), plane stress
#   4---5---6
#   |   |   |
#   1---2---3
ops.node(1, 0.0, 0.0)
ops.node(2, 1.0, 0.0)
ops.node(3, 2.0, 0.0)
ops.node(4, 0.0, 1.0)
ops.node(5, 1.0, 1.0)
ops.node(6, 2.0, 1.0)
ops.fix(1, 1, 1)
ops.fix(4, 1, 1)

ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)
# element quad tag n1 n2 n3 n4 thick type matTag
ops.element("quad", 10, 1, 2, 5, 4, 1.0, "PlaneStress", 1)
ops.element("quad", 11, 2, 3, 6, 5, 1.0, "PlaneStress", 1)

# BOTH recorders, identical element request. Record every step.
ops.recorder("mpco", ref, "-E", "stress", "-T", "dt", 0.0)
ops.recorder("ladruno", new, "-E", "stress", "-T", "dt", 0.0)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(3, 1.0e2, 0.0)
ops.load(6, 1.0e2, 0.0)
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

# =========================================================================
# Block B — fiber-section beam, the deeper OutputDescriptor path
# (element -> section[Gauss] -> fiber -> material -> stress). This is the
# case the schema was explicitly designed around (COLUMN_MAP multiplicity /
# per-block GAUSS_ID). Needs ndf=3, so it runs as a separate model.
#
# KNOWN (separate, non-fatal): with a fiber section, ladruno::writeSections
# emits ~3 HDF5-DIAG diagnostics to stderr while writing the SECTION_<n> fiber
# datasets (LadrunoRecorder.cpp ~792-807). The frozen recorder does not. It
# does NOT corrupt data — SECTION_ASSIGNMENTS (ID/ASSIGNMENT/FIBER_DATA/
# FIBER_MATERIALS) and all element values still match the oracle to 1e-12 (this
# gate passes). Tracked as a writeSections follow-up bug.
# =========================================================================
bref = os.path.join(OUT, "eref_beam.mpco")
bnew = os.path.join(OUT, "etest_beam.ladruno")
for p in (bref, bnew):
    if os.path.exists(p):
        os.remove(p)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 3)
ops.node(1, 0.0, 0.0)
ops.node(2, 1.0, 0.0)
ops.node(3, 2.0, 0.0)
ops.fix(1, 1, 1, 1)

ops.uniaxialMaterial("Elastic", 1, 2.0e5)
# rectangular fiber section: 2x3 = 6 fibers (rich enough for multiplicity)
ops.section("Fiber", 1)
ops.patch("rect", 1, 2, 3, -0.5, -0.25, 0.5, 0.25)
ops.geomTransf("Linear", 1)
ops.beamIntegration("Lobatto", 1, 1, 3)  # 3 integration points
ops.element("dispBeamColumn", 10, 1, 2, 1, 1)
ops.element("dispBeamColumn", 11, 2, 3, 1, 1)

ops.recorder("mpco", bref, "-E", "section.fiber.stress", "-T", "dt", 0.0)
ops.recorder("ladruno", bnew, "-E", "section.fiber.stress", "-T", "dt", 0.0)

ops.timeSeries("Linear", 2)
ops.pattern("Plain", 2, 2)
ops.load(3, 0.0, -1.0e1, 0.0)
ops.system("BandSPD")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 0.25)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(4)
ops.wipe()

print("BEAM REF:", bref, os.path.exists(bref))
print("BEAM NEW:", bnew, os.path.exists(bnew))
print("ELEMENT_PARITY_MODEL OK")
