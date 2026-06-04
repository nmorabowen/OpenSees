"""Shell per-LAYER material-stress emission — model runner (L1).

Regression for the fix that re-enabled the per-layer read through the natural
``material.fiber.<resp>`` verb on a layered shell. Before the fix the
``do_all_materials`` recorder path never expanded the fiber index, so
``-E material.fiber.stress`` silently produced NO element bucket; only the
``section.fiber.*`` path (which the recorder internally swaps section->material
for shells) worked.

This model writes the SAME layered ShellMITC4 result under BOTH verbs:

    recorder ladruno  -E material.fiber.stress   -> shlayer_material.ladruno
    recorder ladruno  -E section.fiber.stress    -> shlayer_section.ladruno

shell_layer_check.py then asserts the material verb now emits a real per-layer
bucket byte-identical to the section verb.

Run with the BUILD python, TMP = dir holding opensees.pyd + the MKL DLLs:
    python shell_layer_model.py <dir_with_opensees_pyd_and_mkl> <out_dir>
"""
from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP

os.add_dll_directory(TMP)
sys.path.insert(0, TMP)
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

import opensees as ops  # noqa: E402

f_mat = os.path.join(OUT, "shlayer_material.ladruno")
f_sec = os.path.join(OUT, "shlayer_section.ladruno")
for p in (f_mat, f_sec):
    if os.path.exists(p):
        os.remove(p)

ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 6)
nid = 0
for j in range(3):
    for i in range(3):
        nid += 1
        ops.node(nid, float(i), float(j), 0.0)
for n in (1, 2, 3):
    ops.fix(n, 1, 1, 1, 1, 1, 1)

# layered shell: 3 plate-fiber layers over an elastic 3D spine
ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.25)
ops.nDMaterial("PlateFiber", 2, 1)
NLAYERS = 3
layers = []
for _ in range(NLAYERS):
    layers += [2, 0.1 / NLAYERS]
ops.section("LayeredShell", 1, NLAYERS, *layers)


def n(i, j):
    return j * 3 + i + 1


etag = 100
for j in range(2):
    for i in range(2):
        etag += 1
        ops.element("ShellMITC4", etag, n(i, j), n(i + 1, j), n(i + 1, j + 1), n(i, j + 1), 1)

# NB: the -E verb is a single dot-joined token (the recorder splits on '.'),
# NOT space-separated args.
ops.recorder("ladruno", f_mat, "-E", "material.fiber.stress",
             "-N", "displacement", "-T", "dt", 0.0)
ops.recorder("ladruno", f_sec, "-E", "section.fiber.stress",
             "-N", "displacement", "-T", "dt", 0.0)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
for nn in (7, 8, 9):
    ops.load(nn, 0.0, 0.0, -1.0e1, 0.0, 0.0, 0.0)
ops.system("BandGeneral")
ops.numberer("RCM")
ops.constraints("Plain")
ops.integrator("LoadControl", 0.5)
ops.algorithm("Linear")
ops.analysis("Static")
ops.analyze(2)
ops.wipe()

print("MAT:", f_mat, os.path.exists(f_mat))
print("SEC:", f_sec, os.path.exists(f_sec))
print("SHELL_LAYER_MODEL OK")
