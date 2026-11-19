"""T10: single-phase MPM answers from one binary, as JSON, for a bit-for-bit comparison.
Usage: python t10_single.py <dist\\bin> <out.json>"""
import json
import os
import sys

B = sys.argv[1]
os.add_dll_directory(B)
sys.path.insert(0, B)
import opensees as ops  # noqa: E402


def lst(v):
    return [float(x) for x in v] if isinstance(v, (list, tuple)) else [float(v)]


def deck(mode, basis, mat, al="off"):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    if mat == "elastic":
        ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.25, 2.0)
    else:
        ops.nDMaterial("DruckerPrager", 1, 5555.56, 4166.67, 5.0, 0.2, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 0.5, "-cells", 4, 8, "-basis", basis)
    ops.mpm("mode", mode)
    ops.mpm("body", 1, "box", 0.0, 0.0, 2.0, 3.0, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    ops.mpm("fix", "xmin", "x")
    ops.mpm("fix", "xmax", "x")
    ops.mpm("fix", "ymin", "y")
    ops.mpm("gravity", 0.0, -9.81)
    if al != "off":
        ops.mpm("antiLocking", "fbar")
    if mode == "explicit":
        ops.mpm("damping", "local", 0.3)
    else:
        ops.constraints("Transformation")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("NormDispIncr", 1e-10, 25)
        ops.algorithm("Newton")
        if mode == "implicit":
            ops.integrator("Newmark", 0.5, 0.25)
            ops.analysis("Transient")
        else:
            ops.integrator("LoadControl", 0.1)
            ops.analysis("Static")


out = {"build": ops.ladrunoBuild()}
for mode, dt, n in (("explicit", 2e-4, 200), ("implicit", 1e-3, 30), ("quasiStatic", 0.1, 10)):
    for basis in ("gimp", "bspline2"):
        for mat in ("elastic", "dp"):
            for al in ("off", "fbar"):
                deck(mode, basis, mat, al)
                ops.mpm("analyze", n, dt)
                out["%s_%s_%s_%s" % (mode, basis, mat, al)] = {f: lst(ops.mpm("particleData", f))
                                                              for f in ("x", "v", "sig")}
json.dump(out, open(sys.argv[2], "w"))
print(out["build"], len(out) - 1)
