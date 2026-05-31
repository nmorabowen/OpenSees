"""Element-result ENVELOPE gate -- model runner.

Exercises the per-column COLUMN_MAP on element ENVELOPES (the gauss/section/fiber
column structure, beyond the result-level COMPONENTS attr -- which is empty for
element results). Each model runs ONE transient with TWO ladruno recorders on
the SAME model so they see identical commitTags:

    recorder ladruno <stem>_ts.ladruno  -E <req>            (time series)
    recorder ladruno <stem>_env.ladruno -E <req> -envelope  (envelope)

A sine load drives sign-changing element response so MIN<0<MAX and ARG_STEP is
non-trivial. The checker reduces the time series to componentwise extremes and
matches the ENVELOPES datasets, and verifies the envelope COLUMN_MAP round-trips
the time-series COLUMN_MAP.

Two blocks:
  A. 2-element FourNodeQuad plane-stress strip, `-E stress` (4 GP x 3 comps =
     12 cols) -- the simple gauss-only COLUMN_MAP.
  B. fiber-section dispBeamColumn, `-E section.fiber.stress` -- the rich
     element->section[Gauss]->fiber->material path (COLUMN_MAP multiplicity /
     per-block GAUSS_ID), the case the schema was designed around.

Run with the BUILD python (no boot .pth):
    python element_envelope_model.py <dir_with_opensees_pyd> [out_dir]
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP

os.add_dll_directory(TMP)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402

NSTEP = 40
DT = 0.02


def _paths(stem):
    ts = os.path.join(OUT, f"{stem}_ts.ladruno")
    env = os.path.join(OUT, f"{stem}_env.ladruno")
    for p in (ts, env):
        if os.path.exists(p):
            os.remove(p)
    return ts, env


# =========================================================================
# Block A -- FourNodeQuad plane-stress strip, -E stress
# =========================================================================
qa_ts, qa_env = _paths("ee_quad")

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
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
for n in (2, 3, 5, 6):
    ops.mass(n, 1.0, 1.0)

ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)
ops.element("quad", 10, 1, 2, 5, 4, 1.0, "PlaneStress", 1)
ops.element("quad", 11, 2, 3, 6, 5, 1.0, "PlaneStress", 1)

ops.recorder("ladruno", qa_ts, "-E", "stress", "-kind", "transient",
             "-T", "dt", 0.0)
ops.recorder("ladruno", qa_env, "-E", "stress", "-kind", "transient",
             "-T", "dt", 0.0, "-envelope")

ops.timeSeries("Trig", 1, 0.0, 100.0, 0.5)  # sine -> sign-changing stress
ops.pattern("Plain", 1, 1)
ops.load(3, 1.0e2, 0.0)
ops.load(6, 1.0e2, 0.0)
ops.constraints("Plain")
ops.numberer("RCM")
ops.system("BandGeneral")
ops.test("NormDispIncr", 1e-8, 10, 0)
ops.algorithm("Linear")
ops.integrator("Newmark", 0.5, 0.25)
ops.analysis("Transient")
for _ in range(NSTEP):
    ops.analyze(1, DT)
ops.wipe()

print("QUAD TS :", qa_ts, os.path.exists(qa_ts))
print("QUAD ENV:", qa_env, os.path.exists(qa_env))

# =========================================================================
# Block B -- fiber-section dispBeamColumn, -E section.fiber.stress
# =========================================================================
ba_ts, ba_env = _paths("ee_beam")

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 3)
ops.node(1, 0.0, 0.0)
ops.node(2, 1.0, 0.0)
ops.node(3, 2.0, 0.0)
ops.fix(1, 1, 1, 1)
for n in (2, 3):
    ops.mass(n, 1.0, 1.0, 0.1)

ops.uniaxialMaterial("Elastic", 1, 2.0e5)
ops.section("Fiber", 1)
ops.patch("rect", 1, 2, 3, -0.5, -0.25, 0.5, 0.25)  # 2x3 = 6 fibers
ops.geomTransf("Linear", 1)
ops.beamIntegration("Lobatto", 1, 1, 3)  # 3 integration points
ops.element("dispBeamColumn", 10, 1, 2, 1, 1)
ops.element("dispBeamColumn", 11, 2, 3, 1, 1)

ops.recorder("ladruno", ba_ts, "-E", "section.fiber.stress",
             "-kind", "transient", "-T", "dt", 0.0)
ops.recorder("ladruno", ba_env, "-E", "section.fiber.stress",
             "-kind", "transient", "-T", "dt", 0.0, "-envelope")

ops.timeSeries("Trig", 2, 0.0, 100.0, 0.5)
ops.pattern("Plain", 2, 2)
ops.load(3, 0.0, -1.0e1, 0.0)
ops.constraints("Plain")
ops.numberer("RCM")
ops.system("BandGeneral")
ops.test("NormDispIncr", 1e-8, 10, 0)
ops.algorithm("Linear")
ops.integrator("Newmark", 0.5, 0.25)
ops.analysis("Transient")
for _ in range(NSTEP):
    ops.analyze(1, DT)
ops.wipe()

print("BEAM TS :", ba_ts, os.path.exists(ba_ts))
print("BEAM ENV:", ba_env, os.path.exists(ba_env))
print("ELEMENT_ENVELOPE_MODEL OK")
