"""ASDShellQ4 / ASDShellT3 ``localAxes`` response and ``MODEL/LOCAL_AXES`` (WP-242) — Zone-A.

The Ladruno recorder writes ``MODEL_STAGE[k]/MODEL/LOCAL_AXES/<tag>-<Class>/{ID, FRAME}``
for every element answering ``"localAxes"`` (9 packed direction cosines e1|e2|e3,
LadrunoRecorder::writeModelLocalAxes). The shells now answer it with their MATERIAL
frame, the one their sections (and so ``section.fiber.crackPattern``) work in:

  * ``-local x y z``: e1 = the vector projected on the shell plane, e3 = the normal,
    e2 = e3 x e1 (ASDShellQ4/T3::setDomain, the ``m_local_x`` branch);
  * no ``-local``: Q4 e1 = (P2+P3)/2 - (P1+P4)/2, T3 e1 = P2 - P1.

The frame comes from the undeformed coordinates, so a ``-corotational`` shell reports
the same frame before and after it moves.

What it proves: the frames equal the hand-computed e1/e2/e3 (a skewed quad catches a
sign error in the rotation by ``m_angle``), a vertical wall with ``-local 1 0 0`` has
|e2_z| > 0.99, the recorder writes the same frames as quaternions, and the existing
force/stress/strain responses keep their sizes.
"""
import glob
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

S2 = 1.0 / math.sqrt(2.0)

# tag: (class, nodes, -local or None, expected e1, e2, e3)
CASES = {
    # horizontal plate, no -local: e1 = global X
    1: ("Q4", [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)], None,
        (1, 0, 0), (0, 1, 0), (0, 0, 1)),
    # skewed plate, no -local: e1 along the mid-edge line (2, 0.5, 0) != the
    # element's Jacobian x-axis, so m_angle != 0 and its sign is exercised
    2: ("Q4", [(0, 0, 0), (2, 0.5, 0), (2, 1.5, 0), (0, 1, 0)], None,
        (2 / math.hypot(2, 0.5), 0.5 / math.hypot(2, 0.5), 0),
        (-0.5 / math.hypot(2, 0.5), 2 / math.hypot(2, 0.5), 0), (0, 0, 1)),
    # wall in the XZ plane, -local 1 0 0: e2 is vertical
    3: ("Q4", [(0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)], (1, 0, 0),
        (1, 0, 0), (0, 0, 1), (0, -1, 0)),
    # the same wall at y=1 with the nodes listed from the top: normal flips, e2 = -Z
    4: ("Q4", [(0, 1, 1), (1, 1, 1), (1, 1, 0), (0, 1, 0)], (1, 0, 0),
        (1, 0, 0), (0, 0, -1), (0, 1, 0)),
    # wall in the YZ plane, -local 0 -1 0
    5: ("Q4", [(2, 0, 0), (2, 1, 0), (2, 1, 1), (2, 0, 1)], (0, -1, 0),
        (0, -1, 0), (0, 0, -1), (1, 0, 0)),
    # horizontal plate, -local off the plane: projected to (1, 1, 0)/sqrt2
    6: ("Q4", [(3, 0, 0), (4, 0, 0), (4, 1, 0), (3, 1, 0)], (1, 1, 5),
        (S2, S2, 0), (-S2, S2, 0), (0, 0, 1)),
    # triangle, no -local: e1 = P2 - P1
    7: ("T3", [(0, 2, 0), (0, 3, 0), (-1, 2, 0)], None,
        (0, 1, 0), (-1, 0, 0), (0, 0, 1)),
    # triangle wall in the XZ plane, -local 1 0 0
    8: ("T3", [(5, 0, 0), (6, 0, 0), (5, 0, 1)], (1, 0, 0),
        (1, 0, 0), (0, 0, 1), (0, -1, 0)),
}
COROT = 9  # copy of case 3 with -corotational, loaded and moved


def _build():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.section("ElasticMembranePlateSection", 1, 3.0e4, 0.2, 0.2, 0.0)
    ntag = 0
    for etag, (kind, xyz, local, *_e) in list(CASES.items()) + [
            (COROT, ("Q4", [(0, 5, 0), (1, 5, 0), (1, 5, 1), (0, 5, 1)], (1, 0, 0)))]:
        tags = []
        for p in xyz:
            ntag += 1
            ops.node(ntag, *map(float, p))
            tags.append(ntag)
        args = ["ASDShellQ4" if kind == "Q4" else "ASDShellT3", etag, *tags, 1]
        if etag == COROT:
            args.append("-corotational")
        if local is not None:
            args += ["-local", *map(float, local)]
        ops.element(*args)
    return ntag


def _assert_frame(la, e1, e2, e3, label):
    la = np.asarray(la, dtype=float)
    assert la.shape == (9,), f"{label}: localAxes must be 9 cosines, got {la.shape}"
    for got, exp, name in ((la[0:3], e1, "e1"), (la[3:6], e2, "e2"), (la[6:9], e3, "e3")):
        assert np.allclose(got, exp, atol=1e-12), f"{label}: {name}={got}, expected {exp}"


def test_asdshell_localaxes_response():
    _build()
    for etag, (kind, _xyz, local, e1, e2, e3) in CASES.items():
        _assert_frame(ops.eleResponse(etag, "localAxes"), e1, e2, e3, f"ele {etag} {kind}")
    # vertical walls with -local 1 0 0: crack directions need e2 vertical
    for etag in (3, 4, 8):
        assert abs(ops.eleResponse(etag, "localAxes")[5]) > 0.99


def test_asdshell_existing_responses_unchanged():
    _build()
    for etag, (kind, *_r) in CASES.items():
        nn = 4 if kind == "Q4" else 3
        assert len(ops.eleResponse(etag, "forces")) == 6 * nn
        assert len(ops.eleResponse(etag, "stresses")) == 8 * nn
        assert len(ops.eleResponse(etag, "strains")) == 8 * nn


def test_asdshell_corotational_frame_is_the_reference_one():
    last = _build()
    before = np.asarray(ops.eleResponse(COROT, "localAxes"))
    _assert_frame(before, *CASES[3][3:], "corotational before")
    # clamp the bottom edge, push the top edge out of plane and rotate the wall
    for n in range(1, last - 1):  # the other shells, and the wall's bottom edge
        ops.fix(n, 1, 1, 1, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (last - 1, last):
        ops.load(n, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0)
    ops.constraints("Plain"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-10, 50); ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1); ops.analysis("Static")
    assert ops.analyze(10) == 0
    assert abs(ops.nodeDisp(last, 2)) > 1e-2  # it actually moved
    after = np.asarray(ops.eleResponse(COROT, "localAxes"))
    assert np.allclose(before, after, atol=1e-14)


def _quat_to_rows(q):
    # schema v1 sec.5: FRAME is a unit quaternion (w, x, y, z) taking global ->
    # element-local, i.e. the matrix whose ROWS are e1, e2, e3
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


def test_asdshell_local_axes_in_recorder(tmp_path):
    h5py = pytest.importorskip("h5py")
    out = str(tmp_path / "shells.ladruno")
    nnodes = _build()
    for n in range(1, nnodes + 1):
        ops.fix(n, 1, 1, 1, 1, 1, 1)
    ops.recorder("ladruno", out, "-E", "stresses")
    ops.constraints("Plain"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-10, 10); ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0); ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.wipe()  # close the HDF5

    files = glob.glob(out) or glob.glob(out.replace(".ladruno", "*.ladruno"))
    assert files, "no .ladruno output produced"
    expected = dict(CASES)
    expected[COROT] = CASES[3]
    seen = {}
    with h5py.File(files[0], "r") as h:
        root = h["MODEL_STAGE[1]/MODEL/LOCAL_AXES"]
        groups = {k: root[k] for k in root}
        assert any(k.endswith("-ASDShellQ4") for k in groups), list(groups)
        assert any(k.endswith("-ASDShellT3") for k in groups), list(groups)
        for name, g in groups.items():
            ids = np.asarray(g["ID"][...]).ravel()
            fr = g["FRAME"][...]
            if fr.dtype.names:  # compound {w, x, y, z}
                fr = np.stack([fr[c] for c in fr.dtype.names], axis=-1)
            fr = np.asarray(fr, dtype=float).reshape(len(ids), 4)
            for eid, q in zip(ids, fr):
                seen[int(eid)] = (name, _quat_to_rows(q))
    assert set(seen) == set(expected), sorted(seen)
    for eid, (name, R) in seen.items():
        _k, _x, _l, e1, e2, e3 = expected[eid]
        for row, exp, lab in ((0, e1, "e1"), (1, e2, "e2"), (2, e3, "e3")):
            assert np.allclose(R[row], exp, atol=1e-9), \
                f"{name} ele {eid}: FRAME {lab}={R[row]}, expected {exp}"
