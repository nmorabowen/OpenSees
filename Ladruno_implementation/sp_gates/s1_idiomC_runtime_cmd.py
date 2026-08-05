"""ADR-80 S1 — IDIOM C: adaptive ladder driven through `ladrunoLoadControl setDeltaLambda`.

Idiom A (re-issuing `integrator LadrunoLoadControl ...` every step) was measured
DEAD: frac=1 was byte-identical to frac=0 and to stock, because each re-issue
constructs a new integrator object and destroys the predictor state. This script
runs the SAME adaptive ladder but resizes the step through the runtime command,
so the object -- and its prediction -- survives across steps and cutbacks.

Python rather than Tcl because `ladrunoLoadControl` is registered in the
interpreter wrappers (PythonWrapper/TclWrapper) and NOT in SRC/tcl/commands.cpp,
so OpenSees.exe does not have it. That gap is pre-existing -- `ladrunoArcLength`
is missing from the .exe for the same reason.

  PYTHONPATH=<worktree>\\dist\\bin python3.12 s1_idiomC_runtime_cmd.py
"""
import json
import sys

import opensees as ops

L, AX, E = 100.0, 10.0, 200000.0
FY, HISO, DELTA = 379.5, 2000.0, 0.15


def build(N, mat):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    h = L / N
    n = 0
    for k in range(N + 1):
        z = k * h
        for x, y in ((0.0, 0.0), (AX, 0.0), (AX, AX), (0.0, AX)):
            n += 1
            ops.node(n, x, y, z)
    if mat == "j2":
        ops.nDMaterial("LadrunoJ2", 1, E / 3.0, E / 2.0,
                       "-iso", "voce", FY, 0.0, 0.0, HISO)
    else:
        ops.nDMaterial("ElasticIsotropic", 1, E, 0.0)
    for k in range(N):
        b = 4 * k
        ops.element("stdBrick", k + 1, b + 1, b + 2, b + 3, b + 4,
                    b + 5, b + 6, b + 7, b + 8, 1)
    for k in range(N + 1):
        b = 4 * k
        for j in range(1, 5):
            if k == 0:
                ops.fix(b + j, 1, 1, 1)
            else:
                ops.fix(b + j, 1, 1, 0)
    return 4 * N


def setup(N, mat, maxIter):
    top = build(N, mat)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(1, 5):
        ops.sp(top + j, 3, DELTA)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-10, maxIter, 0)
    ops.algorithm("KrylovNewton")
    ops.analysis("Static")
    return top


def march(N, mat, mode, frac, maxIter=5):
    """mode: 'reissue' (new object each step) | 'cmd' (setDeltaLambda)."""
    top = setup(N, mat, maxIter)
    dl, lam = 0.1, 0.0
    dlmin, dlmax = 1e-4, 0.2
    inc = cuts = iters = good = guard = 0
    fired = 0
    if mode == "cmd":
        ops.integrator("LadrunoLoadControl", dl, 1, 1e-8, 1.0, "-extrapolate", frac)
    while lam < 1.0 - 1e-12 and guard < 4000:
        guard += 1
        if lam + dl > 1.0:
            dl = 1.0 - lam
        if mode == "cmd":
            ops.ladrunoLoadControl("setDeltaLambda", dl)
            if ops.ladrunoLoadControl("armed") > 0.5:
                fired += 1
        else:
            ops.integrator("LadrunoLoadControl", dl, 1, dl, dl, "-extrapolate", frac)
        rc = ops.analyze(1)
        iters += ops.testIter()
        if rc == 0:
            lam += dl
            inc += 1
            good += 1
            if good >= 2:
                dl = min(dl * 2.0, dlmax)
                good = 0
        else:
            cuts += 1
            good = 0
            dl /= 2.0
            if dl < dlmin:
                return inc, cuts, iters, -1.0, fired
    return inc, cuts, iters, ops.nodeDisp(top + 1, 3), fired


rows = []
print("mode      frac  mat      inc  cutbacks  iters  armed  u_top")
for mode in ("reissue", "cmd"):
    for frac in (0.0, 1.0):
        for mat in ("elastic", "j2"):
            inc, cuts, iters, utop, fired = march(20, mat, mode, frac)
            print(f"{mode:9s} {frac:<5} {mat:8s} {inc:3d}  {cuts:8d}  {iters:5d}  {fired:5d}  {utop:.6f}")
            rows.append(dict(mode=mode, frac=frac, mat=mat, increments=inc,
                             cutbacks=cuts, iters=iters, armed_steps=fired,
                             u_top=utop))

with open("s1_idiomC_runtime_cmd_2026-08-04.json", "w") as fh:
    json.dump({"gate": "ADR-80 S1 idiom C — runtime setDeltaLambda keeps the predictor alive",
               "rows": rows}, fh, indent=2)
print("WROTE s1_idiomC_runtime_cmd_2026-08-04.json")
