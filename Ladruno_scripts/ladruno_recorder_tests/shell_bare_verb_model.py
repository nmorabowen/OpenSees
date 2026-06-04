"""Bare element-verb robustness — model runner (L1).

Regression for the crash where a malformed element-result request of just the
bare keyword (`-E section` / `-E material`, no sub-verb) segfaulted: the recorder
forwarded a zero-length arg list to the section's setResponse, which dereferenced
argv[0]. The fix guards the non-fiber setResponse call sites so such a request is
simply skipped (no bucket) instead of crashing.

This model attaches BOTH bare-verb recorders to a ShellMITC4 model and runs the
analysis. If the guard is missing the model runner crashes (gate fails); if
present, analyze() completes and the files carry no element bucket.

    recorder ladruno  -E section    -> bare_section.ladruno
    recorder ladruno  -E material   -> bare_material.ladruno

Run with the BUILD python, TMP = dir holding opensees.pyd + the MKL DLLs:
    python shell_bare_verb_model.py <dir_with_opensees_pyd_and_mkl> <out_dir>
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

f_sec = os.path.join(OUT, "bare_section.ladruno")
f_mat = os.path.join(OUT, "bare_material.ladruno")
for p in (f_sec, f_mat):
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

# layered shell (so the bare verb has fibers it *could* have expanded, had a
# sub-verb been given) — the bare keyword must still not crash.
ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.25)
ops.nDMaterial("PlateFiber", 2, 1)
layers = []
for _ in range(3):
    layers += [2, 0.1 / 3]
ops.section("LayeredShell", 1, 3, *layers)


def n(i, j):
    return j * 3 + i + 1


etag = 100
for j in range(2):
    for i in range(2):
        etag += 1
        ops.element("ShellMITC4", etag, n(i, j), n(i + 1, j), n(i + 1, j + 1), n(i, j + 1), 1)

# the malformed bare keywords (no sub-verb) + a valid nodal request so the file
# is otherwise well-formed.
ops.recorder("ladruno", f_sec, "-E", "section", "-N", "displacement", "-T", "dt", 0.0)
ops.recorder("ladruno", f_mat, "-E", "material", "-N", "displacement", "-T", "dt", 0.0)

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

print("SEC:", f_sec, os.path.exists(f_sec))
print("MAT:", f_mat, os.path.exists(f_mat))
print("SHELL_BARE_VERB_MODEL OK")
