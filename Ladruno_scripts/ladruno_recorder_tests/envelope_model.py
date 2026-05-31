"""Envelope (ADR D7) gate -- model runner.

ONE transient with TWO recorders on the SAME model: a time series (env_ts.ladruno)
and an `-envelope` recorder (env_env.ladruno). They must see identical commitTags
(running them in separate analyses would offset ARG_STEP). A sine load drives the
free nodes so displacement changes sign over the run -> MIN<0<MAX and a non-trivial
ARG_STEP. The checker reduces the time series to componentwise extremes and compares
to the ENVELOPES datasets (the envelope must equal the reduction it summarizes).

Run with the BUILD python (no boot .pth):
    python envelope_model.py <dir_with_opensees_pyd> [out_dir]
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

ts_out = os.path.join(OUT, "env_ts.ladruno")
env_out = os.path.join(OUT, "env_env.ladruno")
for p in (ts_out, env_out):
    if os.path.exists(p):
        os.remove(p)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.node(1, 0.0, 0.0)
ops.node(2, 100.0, 0.0)
ops.node(3, 200.0, 0.0)
ops.fix(1, 1, 1)
ops.fix(2, 0, 1)
ops.fix(3, 0, 1)
ops.mass(2, 1.0, 1.0)
ops.mass(3, 1.0, 1.0)
ops.uniaxialMaterial("Elastic", 1, 1000.0)
ops.element("Truss", 1, 1, 2, 1.0, 1)
ops.element("Truss", 2, 2, 3, 1.0, 1)

# both recorders on the same model -> identical commitTags
ops.recorder("ladruno", ts_out, "-N", "displacement", "-kind", "transient",
             "-T", "dt", 0.0)
ops.recorder("ladruno", env_out, "-N", "displacement", "-kind", "transient",
             "-T", "dt", 0.0, "-envelope")

ops.timeSeries("Trig", 1, 0.0, 100.0, 0.5)  # sine -> sign-changing displacement
ops.pattern("Plain", 1, 1)
ops.load(3, 5.0, 0.0)
ops.constraints("Plain")
ops.numberer("Plain")
ops.system("BandGeneral")
ops.test("NormDispIncr", 1e-8, 10, 0)
ops.algorithm("Linear")
ops.integrator("Newmark", 0.5, 0.25)
ops.analysis("Transient")
for _ in range(NSTEP):
    ops.analyze(1, DT)
ops.wipe()  # flush + close both recorders

print("NEW:", ts_out, os.path.exists(ts_out))
print("NEW:", env_out, os.path.exists(env_out))
print("ENVELOPE_MODEL OK")
