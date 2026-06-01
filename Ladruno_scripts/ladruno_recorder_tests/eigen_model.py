"""Eigen / modesOfVibration coverage — model runner (L2 validate path).

Exercises the recorder's MODAL output path (KIND='eigen', ON_NODES/
MODESOFVIBRATION) which the per-step parity gate deliberately skips. A small
portal frame with lumped mass, an eigen solve for 3 modes, then a forced
ops.record() so the recorder captures the eigenvectors.

    recorder ladruno -> eig_test.ladruno
    recorder mpco    -> eig_ref.mpco   (kept for manual side-by-side; not value-diffed)

Run with the BUILD python, TMP = dir holding opensees.pyd + the MKL DLLs:
    python eigen_model.py <dir_with_opensees_pyd_and_mkl> <out_dir>
"""

from __future__ import annotations

import os
import sys

TMP = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP

os.add_dll_directory(TMP)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402

ref = os.path.join(OUT, "eig_ref.mpco")
new = os.path.join(OUT, "eig_test.ladruno")
for p in (ref, new):
    if os.path.exists(p):
        os.remove(p)

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 3)

# single-bay single-story portal frame
ops.node(1, 0.0, 0.0)
ops.node(2, 4.0, 0.0)
# small rotational mass (3rd dof) keeps the generalized mass matrix non-singular
ops.node(3, 0.0, 3.0, "-mass", 2.0, 2.0, 1.0e-4)
ops.node(4, 4.0, 3.0, "-mass", 2.0, 2.0, 1.0e-4)
ops.fix(1, 1, 1, 1)
ops.fix(2, 1, 1, 1)

E = 2.0e5
A, I = 1.0e2, 1.0e3
ops.geomTransf("Linear", 1)
ops.element("elasticBeamColumn", 1, 1, 3, A, E, I, 1)   # left column
ops.element("elasticBeamColumn", 2, 2, 4, A, E, I, 1)   # right column
ops.element("elasticBeamColumn", 3, 3, 4, A, E, I, 1)   # beam

ops.recorder("mpco", ref, "-N", "modesOfVibration", "-T", "dt", 0.0)
ops.recorder("ladruno", new, "-N", "modesOfVibration", "-T", "dt", 0.0)

ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("BandGeneral")
vals = ops.eigen(3)
print("EIGENVALUES:", vals)
ops.record()   # force the recorders to capture the eigenvectors
ops.wipe()

print("REF:", ref, os.path.exists(ref))
print("NEW:", new, os.path.exists(new))
print("EIGEN_MODEL OK")
