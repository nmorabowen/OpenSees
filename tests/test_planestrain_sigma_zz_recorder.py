"""stressesPlaneStrain records into `.ladruno` — end-to-end (PR #525 follow-up).

The Ladruno recorder's `-E` verb is generic: any element-response token is
probed via `setResponse` and self-describes through the OutputDescriptor XML,
so the new 4-component `"stressesPlaneStrain"` response (PR #525) needs ZERO
recorder edits. This battery pins that contract so a recorder-side regression
(a grammar tightening, a bucket-drop workaround misfiring) can't silently take
the sigma_zz channel away from apeGmsh:

  * one bucket per element class under RESULTS/ON_ELEMENTS/stressesPlaneStrain,
    NUM_COLUMNS = 4*nGP, for all five plane elements in ONE model,
  * COLUMN_MAP/COMP_NAMES carries "sigma11,sigma22,sigma12,sigma33" per GP,
  * recorded last-step DATA == the live eleResponse (same material query),
  * recorded szz == nu*(sxx+syy) (elastic plane strain),
  * a plane-STRESS material round-trips its NaN sigma_zz through HDF5 intact
    (the "not exposed" contract survives the file).
"""
import glob
import math
import os

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, NU, DELTA = 30_000.0, 0.3, 1.0e-4

# (classTag-prefixed bucket key fragment, element tag, nGP)
BUCKETS = [
    ("FourNodeQuad", 1, 4),
    ("Tri31", 2, 1),
    ("LadrunoQuad", 3, 4),
    ("LadrunoCST", 4, 1),
    ("BezierTri6", 5, 3),
]


def _patch(x0, shape):
    """Standalone unit patch at x-offset x0: all x fixed, bottom y fixed.
    Returns (node_tags, top_node_tags)."""
    crds = {
        "quad": [(0, 0), (1, 0), (1, 1), (0, 1)],
        "tri": [(0, 0), (1, 0), (0, 1)],
        "tri6": [(0, 0), (1, 0), (0, 1), (0.5, 0), (0.5, 0.5), (0, 0.5)],
    }[shape]
    tags, tops = [], []
    for (x, y) in crds:
        t = ops.getNodeTags()[-1] + 1 if ops.getNodeTags() else 1
        ops.node(t, x0 + x, float(y))
        ops.fix(t, 1, 1 if y == 0 else 0)
        if y == 1:
            tops.append(t)
        tags.append(t)
    return tags, tops


def _drive(tops, nsteps=5):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in tops:
        ops.sp(t, 2, -DELTA)
    # Penalty: fully-prescribed rig — Transformation would condense to a
    # zero-size SOE and fatally exit (see LEDGER_quirks)
    ops.constraints("Penalty", 1.0e15, 1.0e15)
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    assert ops.analyze(nsteps) == 0


def _find_buckets(h5file, h5py, result="stressesPlaneStrain"):
    """{class-name-fragment: group} for every bucket under the result."""
    found = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Group) and f"ON_ELEMENTS/{result}/" in name \
                and name.count("/") == 4:
            found[name.rsplit("/", 1)[-1]] = obj

    h5file.visititems(visit)
    return found


def test_ladruno_recorder_stresses_plane_strain(tmp_path):
    h5py = pytest.importorskip("h5py")
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    out = str(tmp_path / "szz.ladruno")

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)

    all_tops = []
    t, tops = _patch(0.0, "quad")
    ops.element("quad", 1, *t, 1.0, "PlaneStrain", 1); all_tops += tops
    t, tops = _patch(3.0, "tri")
    ops.element("Tri31", 2, *t, 1.0, "PlaneStrain", 1); all_tops += tops
    t, tops = _patch(6.0, "quad")
    ops.element("LadrunoQuad", 3, *t, 1, "-type", "PlaneStrain", "-thick", 1.0); all_tops += tops
    t, tops = _patch(9.0, "tri")
    ops.element("LadrunoCST", 4, *t, 1, "-type", "PlaneStrain", "-thick", 1.0); all_tops += tops
    t, tops = _patch(12.0, "tri6")
    ops.element("BezierTri6", 5, *t, 1.0, "PlaneStrain", 1); all_tops += tops

    ops.recorder("ladruno", out, "-E", "stressesPlaneStrain")
    _drive(all_tops)

    # live values BEFORE the wipe closes the recorder
    live = {name: list(ops.eleResponse(tag, "stressesPlaneStrain"))
            for (name, tag, _) in BUCKETS}
    ops.wipe()  # destroy recorders -> flush/close the HDF5

    files = glob.glob(out) or glob.glob(out + "*")
    assert files, "no .ladruno output produced"

    with h5py.File(files[0], "r") as h:
        buckets = _find_buckets(h, h5py)
        for (name, tag, ngp) in BUCKETS:
            keys = [k for k in buckets if name in k]
            assert keys, f"no {name} bucket in {sorted(buckets)}"
            g = buckets[keys[0]]

            # shape: [nsteps, nEle, 4*nGP]
            data = np.asarray(g["DATA"][...], dtype=float)
            assert data.shape[-1] == 4 * ngp, (name, data.shape)
            assert int(g.attrs["NUM_COLUMNS"][0]) == 4 * ngp

            # per-GP component naming carries sigma33
            comp = g["COLUMN_MAP"].attrs["COMP_NAMES"]
            comp = comp[0].decode() if isinstance(comp, np.ndarray) else str(comp)
            rows = [r for r in comp.splitlines() if r.strip()]
            assert len(rows) == ngp, (name, comp)
            assert all(r == "sigma11,sigma22,sigma12,sigma33" for r in rows), (name, comp)

            # recorded last step == live eleResponse (same material query)
            rec = data[-1].reshape(-1)
            assert rec == pytest.approx(live[name], rel=1e-12, abs=1e-14), name

            # and the physics: szz == nu*(sxx+syy) at every GP
            for gp in range(ngp):
                sxx, syy, _, szz = rec[4 * gp: 4 * gp + 4]
                assert szz == pytest.approx(NU * (sxx + syy), rel=1e-8)
                assert abs(szz) > 0.0


def test_ladruno_recorder_nan_survives_hdf5(tmp_path):
    """Plane-stress material: the NaN sigma_zz must round-trip the file as NaN
    (never silently become 0)."""
    h5py = pytest.importorskip("h5py")
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    out = str(tmp_path / "szz_nan.ladruno")

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    t, tops = _patch(0.0, "quad")
    ops.element("quad", 1, *t, 1.0, "PlaneStress", 1)
    ops.recorder("ladruno", out, "-E", "stressesPlaneStrain")
    _drive(tops)
    ops.wipe()

    files = glob.glob(out) or glob.glob(out + "*")
    assert files
    with h5py.File(files[0], "r") as h:
        buckets = _find_buckets(h, h5py)
        keys = [k for k in buckets if "FourNodeQuad" in k]
        assert keys, f"no FourNodeQuad bucket in {sorted(buckets)}"
        rec = np.asarray(buckets[keys[0]]["DATA"][...], dtype=float)[-1].reshape(-1)
        for gp in range(4):
            sxx, syy, sxy, szz = rec[4 * gp: 4 * gp + 4]
            assert math.isnan(szz), f"gp{gp}: NaN sigma_zz lost in file (got {szz})"
            assert all(math.isfinite(v) for v in (sxx, syy, sxy))
