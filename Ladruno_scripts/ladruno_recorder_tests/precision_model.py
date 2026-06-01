"""Precision-mode model runner: run the SAME model TWICE — once with the default
lossless f64 DATA, once with the opt-in `-precision f32` lossy DATA — so the two
.ladruno files differ ONLY in on-disk result precision. precision_check.py then
asserts (a) f32 stays within a float32 error bound of f64, (b) the f32 file is
meaningfully smaller, (c) STORED_PRECISION / on-disk dtype are tagged correctly.

  python precision_model.py <pyd_dir> <out_dir>
"""
import os
import sys

TMP, OUT = sys.argv[1], sys.argv[2]
DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\.claude\worktrees\nifty-pasteur-c1ad8f\dist\bin"
os.add_dll_directory(DIST)
sys.path.insert(0, TMP)
import opensees as ops  # noqa: E402

f64 = os.path.join(OUT, "prec_f64.ladruno")
f32 = os.path.join(OUT, "prec_f32.ladruno")
for p in (f64, f32):
    if os.path.exists(p):
        os.remove(p)

NX, NY = 20, 6
L, H = 6.0, 2.0
dx, dy = L / NX, H / NY


def nid(i, j):
    return 1 + i + j * (NX + 1)


def make_mesh():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for j in range(NY + 1):
        for i in range(NX + 1):
            ops.node(nid(i, j), i * dx, j * dy)
    for j in range(NY + 1):
        ops.fix(nid(0, j), 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e7, 0.25, 2.0)
    etag = 1
    for j in range(NY):
        for i in range(NX):
            ops.element("quad", etag, nid(i, j), nid(i + 1, j),
                        nid(i + 1, j + 1), nid(i, j + 1), 0.1, "PlaneStress", 1)
            etag += 1


def run_load():
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(nid(NX, NY), 0.0, -40.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    ops.analyze(120, 1.0e-3)
    ops.wipe()  # flush + close the recorder


# run 1: default lossless f64
make_mesh()
ops.recorder("ladruno", f64, "-N", "displacement", "-E", "stress",
             "-precision", "f64", "-T", "dt", 0.0)
run_load()

# run 2: opt-in lossy f32
make_mesh()
ops.recorder("ladruno", f32, "-N", "displacement", "-E", "stress",
             "-precision", "f32", "-T", "dt", 0.0)
run_load()

print("F64:", f64, os.path.exists(f64), os.path.getsize(f64))
print("F32:", f32, os.path.exists(f32), os.path.getsize(f32))
print("PRECISION_MODEL OK")
