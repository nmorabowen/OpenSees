"""ADR-30 Phase 4 (recorder) — tie force is recordable via LadrunoRecorder.

The lean query `ladrunoProjectionTieForce` (P3) is now joined by NATIVE recording: the
projector caches f = M(a_raw - a_proj); CentralDifferenceLadruno scatters it onto the
nodes each commit (Node::setProjectionTieForce); and a new LadrunoRecorder node source
`constraintTieForce` writes it to the .ladruno HDF5 (group CONSTRAINT_TIE_FORCE). The
recorded value comes from the SAME source as the P3 query (validated against analytics in
T8), so this test confirms the record CHAIN end-to-end and the written values.
"""
import glob
import os

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"
PROJ = "LadrunoProjection"


def test_P4_tie_force_recordable_via_ladruno_recorder(tmp_path):
    h5py = pytest.importorskip("h5py")

    m1, m2, F = 2.0, 3.0, 10.0
    f_exact = F * m2 / (m1 + m2)          # 6.0 — the internal tie force (±)
    out = str(tmp_path / "tieforce.ladruno")

    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.mass(1, m1)
    ops.node(2, 1.0); ops.mass(2, m2)
    ops.equalDOF(1, 2, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, F)
    ops.constraints(PROJ)
    ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.recorder("ladruno", out, "-N", "constraintTieForce")
    ops.analysis("Transient")
    for _ in range(10):
        assert ops.analyze(1, 1e-3) == 0
    ops.wipe()                            # destroy recorders -> flush/close the HDF5

    files = glob.glob(out) or glob.glob(out.replace(".ladruno", "*.ladruno"))
    assert files, "no .ladruno output produced"
    assert os.path.getsize(files[0]) > 0

    # readback: the field lives at .../CONSTRAINT_TIE_FORCE/DATA (shape [nstep, nnode,
    # ncomp]); the sibling ID/STEP/TIME datasets are metadata — read DATA only.
    found = [None]
    with h5py.File(files[0], "r") as h:
        def visit(name, obj):
            up = name.upper()
            if (isinstance(obj, h5py.Dataset) and "CONSTRAINT_TIE_FORCE" in up
                    and up.rstrip("/").endswith("DATA")):
                found[0] = np.asarray(obj[...], dtype=float)
        h.visititems(visit)

    assert found[0] is not None, "CONSTRAINT_TIE_FORCE/DATA not found in the .ladruno output"
    data = found[0]
    assert data.size, "tie-force DATA empty"
    # both tied nodes carry the internal constraint force ±F*m2/M, equal-and-opposite
    assert abs(np.max(np.abs(data)) - f_exact) < 1e-6, (
        f"recorded tie force max {np.max(np.abs(data))} != expected {f_exact}")
    # per step the two nodes are equal-and-opposite (sum over nodes ~ 0)
    assert np.max(np.abs(data.sum(axis=tuple(range(1, data.ndim))))) < 1e-9
