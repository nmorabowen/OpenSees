"""Per-stage KIND + stage-boundary regression -- openseespy runner (apeGmsh #2034).

The openseespy twin of stage_kind_model.tcl for the cases whose code path
differs between the interpreters: `ops.domainChange()` goes through
OPS_domainChange (interpreter/OpenSeesMiscCommands.cpp), not the classic Tcl
command, and the running analysis is read through OpenSeesCommands' accessors.

    python stage_kind_model.py <dir_with_opensees_pyd> <out.ladruno> <case>

Cases (same expectations as the Tcl deck, see stage_kind_check.py):
    dc_static_transient, flip_static_transient, dc_static_static, internal_restamp
Prints "STAGE_KIND_MODEL OK <case>" last on success.
"""

from __future__ import annotations

import os
import sys

PYD_DIR = sys.argv[1]
OUT = sys.argv[2]
CASE = sys.argv[3]

os.add_dll_directory(PYD_DIR)          # MKL etc. sit beside opensees.pyd in dist\bin
sys.path.insert(0, PYD_DIR)

import opensees as ops  # noqa: E402

if os.path.exists(OUT):
    os.remove(OUT)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 3)
ops.node(1, 0.0, 0.0)
ops.node(2, 0.0, 3.0)
ops.fix(1, 1, 1, 1)
ops.mass(2, 10.0, 10.0, 0.0)
ops.geomTransf("Linear", 1)
ops.element("elasticBeamColumn", 1, 1, 2, 0.1, 2.0e8, 1.0e-3, 1)
ops.recorder("ladruno", OUT, "-N", "displacement", "-T", "dt", 0.0)


def chain(kind):
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("BandGen")
    ops.test("NormDispIncr", 1.0e-10, 10)
    ops.algorithm("Newton")
    if kind == "static":
        ops.integrator("LoadControl", 0.25)
        ops.analysis("Static")
    else:
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")


def gravity():
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -100.0, 0.0)


def lateral():
    ops.timeSeries("Sine", 2, 0.0, 10.0, 0.5)
    ops.pattern("Plain", 2, 2)
    ops.load(2, 100.0, 0.0, 0.0)


def ok(rc):
    if rc != 0:
        raise SystemExit(f"analyze failed ({rc}) in case {CASE}")


if CASE in ("dc_static_transient", "flip_static_transient"):
    gravity(); chain("static"); ok(ops.analyze(4))
    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()
    if CASE == "dc_static_transient":
        ops.domainChange()
    lateral(); chain("transient"); ok(ops.analyze(5, 0.01))
elif CASE == "dc_static_static":
    gravity(); chain("static"); ok(ops.analyze(4))
    ops.loadConst("-time", 0.0)
    ops.timeSeries("Linear", 3)
    ops.pattern("Plain", 3, 3)
    ops.load(2, 50.0, 0.0, 0.0)
    ops.domainChange()
    ok(ops.analyze(3))
elif CASE == "internal_restamp":
    gravity(); chain("static"); ok(ops.analyze(2))
    ops.pattern("Plain", 5, 1)
    ops.sp(2, 3, 0.0)
    ok(ops.analyze(2))
    ops.pattern("Plain", 6, 1)
    ops.eleLoad("-ele", 1, "-type", "-beamUniform", -1.0)
    ok(ops.analyze(2))
else:
    raise SystemExit(f"unknown case {CASE}")

ops.remove("recorders")
print(f"STAGE_KIND_MODEL OK {CASE}", flush=True)
