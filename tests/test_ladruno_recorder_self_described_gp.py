"""WP-227 (#1067) — `recorder ladruno` QUADRATURE/GP_PARAM for the Ladruno solids
WP-226 left out.

After WP-226 (#1066), five Ladruno solids still had no arm in
getGeometryAndIntRuleByClassTag: their MODEL/ELEMENTS group was TOPOLOGY custom
with no QUADRATURE, so readers could not place their Gauss points. This file
pins each fix against the element's OWN Gauss loop (count AND order):

* LadrunoLST — the 3-pt interior rule in lst_init_rule order
  (2/3,1/6),(1/6,2/3),(1/6,1/6) as the newly tabulated Triangle_GaussLegendre_2,
  + a write-side quadratic Lagrange T6 basis for GLOBAL_GP_COORDS.
* LadrunoSolidShell — 2x2 in-plane x nz (gauss|lobatto) per element, answered by
  the element ("integrationPoints"/"integrationWeights"), layer-major.
* LadrunoUP — five topologies under one class tag; the element self-describes
  ("basisInfo" + GP probes) and the writer derives the geometry PER BUCKET, so a
  T3 + Q4 mesh lands in two rectangular buckets. Its H8 order (bottom CCW, then
  top) is NOT Brick's i,j,k order.
* Bernstein elements now get GLOBAL_GP_COORDS from a Bernstein basis (curved
  BezierTri6), never from the Lagrange T6 that shares the geometry enum.

Each element's per-GP responses now carry GaussPoint tags, so DATA is bucketed
GAUSS_ID 0..nGP-1 and is checked against eleResponse.

Out of scope (documented, not fixed): LadrunoCSTPair (its GPs are the two
sub-triangle centroids, not fixed points of the quad parent domain on a
distorted quad) and LadrunoCPML (no per-GP responses at all).
"""
import glob
import os

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

A = 1.0 / np.sqrt(3.0)
HEX8_NODES = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                       [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], dtype=float)
SKEW3 = {1: (0.0, 0.0, 0.0), 2: (1.2, 0.1, 0.0), 3: (1.1, 1.3, 0.05), 4: (-0.1, 0.9, 0.0),
         5: (0.05, 0.0, 0.3), 6: (1.0, 0.1, 0.35), 7: (1.2, 1.1, 0.32), 8: (0.0, 1.0, 0.28)}


# ----------------------------------------------------------------- helpers
def _h5py():
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    return pytest.importorskip("h5py")


def _file(out):
    files = glob.glob(out) or glob.glob(out + "*")
    assert files, "no .ladruno output produced"
    return files[0]


def _stage(h):
    return h[[s for s in h if s.startswith("MODEL_STAGE")][0]]


def _groups(h, cls):
    """{name: group} of MODEL/ELEMENTS/<tag>-<cls>[r:c]."""
    els = _stage(h)["MODEL/ELEMENTS"]
    return {k: els[k] for k in els if f"-{cls}[" in k}


def _bucket(h, resp, name):
    """The results bucket of model group <tag>-<cls>[r:c]: <tag>-<cls>[r:c:header]."""
    on = _stage(h)[f"RESULTS/ON_ELEMENTS/{resp}"]
    hits = [k for k in on if k.startswith(name[:-1] + ":")]
    assert len(hits) == 1, (name, sorted(on))
    return on[hits[0]]


def _s(v):
    v = np.atleast_1d(v).reshape(-1)[0]
    return v.decode() if isinstance(v, bytes) else str(v)


def _i(v):
    return int(np.atleast_1d(v).reshape(-1)[0])


def _check_data(h, resp, name, live_by_tag, ngp):
    """DATA[-1] of bucket <name> == eleResponse per element, GAUSS_ID 0..ngp-1."""
    b = _bucket(h, resp, name)
    gid = np.asarray(b["COLUMN_MAP/GAUSS_ID"]).reshape(-1)
    ncol = gid.size
    assert ncol % ngp == 0
    want_gid = np.repeat(np.arange(ngp), ncol // ngp)
    assert gid.tolist() == want_gid.tolist()
    ids = np.asarray(b["ID"][...]).reshape(-1)
    data = np.asarray(b["DATA"][...], dtype=float)[-1]
    for r, tag in enumerate(ids):
        assert data[r] == pytest.approx(live_by_tag[int(tag)], rel=1e-10, abs=1e-12)


def _static(system="FullGeneral"):
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system(system)
    ops.test("NormDispIncr", 1e-10, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


def _t6_lagrange(x, e):
    z = 1.0 - x - e
    return np.array([x * (2 * x - 1), e * (2 * e - 1), z * (2 * z - 1),
                     4 * x * e, 4 * e * z, 4 * z * x])


def _t6_bernstein(x1, x2):
    x3 = 1.0 - x1 - x2
    return np.array([x3 * x3, x1 * x1, x2 * x2, 2 * x1 * x3, 2 * x1 * x2, 2 * x2 * x3])


# ------------------------------------------------------------------ LadrunoLST
LST_GP = np.array([[2 / 3, 1 / 6], [1 / 6, 2 / 3], [1 / 6, 1 / 6]])
# corners 1-3, mid-sides 4:(1-2) 5:(2-3) 6:(3-1); node 4 and 5 pushed OFF the
# chords so the element is curved and x(xi) is genuinely quadratic
LST_NODES = {1: (0.0, 0.0), 2: (2.0, 0.1), 3: (0.3, 1.8),
             4: (1.0, -0.15), 5: (1.3, 1.1), 6: (0.12, 0.9)}


@pytest.mark.parametrize("geom", ["linear", "finite"])
def test_lst_writes_tri_gl2_in_element_order(tmp_path, geom):
    h5py = _h5py()
    out = str(tmp_path / f"lst_{geom}.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    if geom == "finite":   # -geom finite needs a setTrialF material (LogStrain2D)
        ops.nDMaterial("ElasticIsotropic", 9, 3.0e4, 0.2)
        ops.nDMaterial("LogStrain2D", 1, 9)
    else:
        ops.nDMaterial("ElasticIsotropic", 1, 3.0e4, 0.2)
    for n, xy in LST_NODES.items():
        ops.node(n, *xy)
    ops.element("LadrunoLST", 1, *LST_NODES.keys(), 1, "-type", "PlaneStrain",
                "-geom", geom)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 1.0, -0.4)
    ops.load(5, -0.3, 0.7)
    ops.recorder("ladruno", out, "-E", "stress")
    _static()
    live = np.asarray(ops.eleResponse(1, "stress"), dtype=float)
    assert live.size == 9
    ops.wipe()

    with h5py.File(_file(out), "r") as h:
        groups = _groups(h, "LadrunoLST")
        assert len(groups) == 1, sorted(groups)
        name, g = next(iter(groups.items()))
        assert _s(g.attrs["TOPOLOGY"]) == "tri"
        assert _s(g.attrs["FAMILY"]) == "lagrange"
        assert _i(g.attrs["NUM_GP"]) == 3 and _i(g.attrs["NDIR"]) == 2
        gp = np.asarray(g["QUADRATURE/GP_PARAM"][...], dtype=float)
        assert np.allclose(gp, LST_GP, atol=1e-15)
        w = np.asarray(g["QUADRATURE/GP_WEIGHT"][...], dtype=float).reshape(-1)
        assert np.allclose(w, 1.0 / 6.0) and w.size == 3
        X = np.asarray([LST_NODES[n] for n in LST_NODES], dtype=float)
        ggp = np.asarray(g["GLOBAL_GP_COORDS"][...], dtype=float).reshape(3, 2)
        for k in range(3):
            assert np.allclose(ggp[k], _t6_lagrange(*gp[k]) @ X, atol=1e-12)
        _check_data(h, "stress", name, {1: live}, 3)


# ------------------------------------------------------------ LadrunoSolidShell
G2 = 0.5773502691896257
IP4 = np.array([[-G2, -G2], [G2, -G2], [G2, G2], [-G2, G2]])
ZRULES = {("gauss", 2): ([-G2, G2], [1.0, 1.0]),
          ("lobatto", 3): ([-1.0, 0.0, 1.0], [1 / 3, 4 / 3, 1 / 3])}


def _solidshell_rule(quadz, nz):
    zp, zw = ZRULES[(quadz, nz)]
    pts = np.array([[IP4[q, 0], IP4[q, 1], zp[l]] for l in range(nz) for q in range(4)])
    wts = np.array([zw[l] for l in range(nz) for q in range(4)])
    return pts, wts


def test_solidshell_custom_rule_per_element(tmp_path):
    """Two shells with different (quadz, nz) in one model -> two buckets, each
    GP_PARAM in the element's layer-major loop order, trilinear GLOBAL_GP_COORDS
    on a skewed brick, and DATA == eleResponse per GP."""
    h5py = _h5py()
    out = str(tmp_path / "solidshell.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 3.0e4, 0.2)
    for n, xyz in SKEW3.items():
        ops.node(n, *xyz)
    for n, xyz in SKEW3.items():          # a second, shifted copy
        ops.node(10 + n, xyz[0] + 2.0, xyz[1], xyz[2])
    ops.element("LadrunoSolidShell", 1, *SKEW3.keys(), 1, "-nz", 2, "-quadz", "gauss")
    ops.element("LadrunoSolidShell", 2, *[10 + n for n in SKEW3], 1, "-nz", 3, "-quadz", "lobatto")
    for n in (1, 2, 3, 4, 11, 12, 13, 14):
        ops.fix(n, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for base in (0, 10):
        ops.load(base + 7, 2.0, -1.0, -5.0)
        ops.load(base + 6, 0.0, 3.0, 0.0)
    ops.recorder("ladruno", out, "-E", "stress")
    _static()
    live = {t: np.asarray(ops.eleResponse(t, "stress"), dtype=float) for t in (1, 2)}
    assert live[1].size == 48 and live[2].size == 72
    ops.wipe()

    with h5py.File(_file(out), "r") as h:
        groups = _groups(h, "LadrunoSolidShell")
        assert len(groups) == 2, sorted(groups)
        seen = set()
        for name, g in groups.items():
            assert _s(g.attrs["TOPOLOGY"]) == "hex"
            ngp = _i(g.attrs["NUM_GP"])
            conn = np.asarray(g["CONNECTIVITY"][...])
            tag = int(conn[0, 0])
            quadz, nz = ("gauss", 2) if tag == 1 else ("lobatto", 3)
            assert ngp == 4 * nz
            pts, wts = _solidshell_rule(quadz, nz)
            gp = np.asarray(g["QUADRATURE/GP_PARAM"][...], dtype=float)
            assert gp.shape == (ngp, 3) and np.allclose(gp, pts, atol=1e-14)
            w = np.asarray(g["QUADRATURE/GP_WEIGHT"][...], dtype=float).reshape(-1)
            assert np.allclose(w, wts, atol=1e-14) and w.sum() == pytest.approx(8.0)
            X = np.asarray([SKEW3[n] for n in SKEW3], dtype=float)
            if tag == 2:
                X = X + np.array([2.0, 0.0, 0.0])
            ggp = np.asarray(g["GLOBAL_GP_COORDS"][...], dtype=float).reshape(ngp, 3)
            for L in range(ngp):
                N = np.prod(1.0 + HEX8_NODES * gp[L], axis=1) / 8.0
                assert np.allclose(ggp[L], N @ X, atol=1e-12)
            _check_data(h, "stress", name, live, ngp)
            seen.add(tag)
        assert seen == {1, 2}


# ------------------------------------------------------------------ LadrunoUP
UP_OPTS = ["-Kf", 2.2e6, "-poro", 0.4, "-rhoF", 1.0, "-perm", 1e-4, 1e-4]
UP_OPTS3 = ["-Kf", 2.2e6, "-poro", 0.4, "-rhoF", 1.0, "-perm", 1e-4, 1e-4, 1e-4]
# 2D: a distorted Q4 (1,2,3,4) + a T3 (2,5,3) sharing the 2-3 edge
UP2_NODES = {1: (0.0, 0.0), 2: (1.0, 0.1), 3: (1.1, 1.0), 4: (-0.1, 0.9), 5: (2.0, 0.6)}
Q4_GP = np.array([[-A, -A], [A, -A], [A, A], [-A, A]])
# UP's H8 order: bottom layer CCW, then top (LadrunoUPShapes.h H8::gpCoord)
UP_H8_GP = np.array([[-A, -A, -A], [A, -A, -A], [A, A, -A], [-A, A, -A],
                     [-A, -A, A], [A, -A, A], [A, A, A], [-A, A, A]])


def test_up_mixed_t3_q4_buckets_self_described(tmp_path):
    h5py = _h5py()
    out = str(tmp_path / "up2d.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 3.0e4, 0.2)
    for n, xy in UP2_NODES.items():
        ops.node(n, *xy)
    ops.element("LadrunoUP", 1, 1, 2, 3, 4, 1, *UP_OPTS)
    ops.element("LadrunoUP", 2, 2, 5, 3, 1, *UP_OPTS)
    ops.fix(1, 1, 1, 1)
    ops.fix(4, 1, 1, 1)
    for n in (2, 3, 5):
        ops.fix(n, 0, 0, 1)                # drained: p = 0 everywhere
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 1.0, -0.5, 0.0)
    ops.load(3, -0.2, 0.4, 0.0)
    ops.recorder("ladruno", out, "-E", "stress")
    _static()
    live = {t: np.asarray(ops.eleResponse(t, "stress"), dtype=float) for t in (1, 2)}
    assert live[1].size == 12 and live[2].size == 3
    ops.wipe()

    with h5py.File(_file(out), "r") as h:
        groups = _groups(h, "LadrunoUP")
        assert len(groups) == 2, sorted(groups)
        for name, g in groups.items():
            assert "CONNECTIVITY_LAYOUT" not in g.attrs     # one node count per bucket
            conn = np.asarray(g["CONNECTIVITY"][...])
            tag = int(conn[0, 0])
            gp = np.asarray(g["QUADRATURE/GP_PARAM"][...], dtype=float)
            ggp = np.asarray(g["GLOBAL_GP_COORDS"][...], dtype=float).reshape(gp.shape[0], 2)
            assert _s(g.attrs["FAMILY"]) == "lagrange"
            if tag == 1:
                assert conn.shape == (1, 5) and _s(g.attrs["TOPOLOGY"]) == "quad"
                assert _i(g.attrs["NUM_CTRL"]) == 4
                assert np.allclose(gp, Q4_GP, atol=1e-15)
                X = np.asarray([UP2_NODES[n] for n in (1, 2, 3, 4)])
                for k in range(4):
                    s, t = gp[k]
                    N = 0.25 * np.array([(1 - s) * (1 - t), (1 + s) * (1 - t),
                                         (1 + s) * (1 + t), (1 - s) * (1 + t)])
                    assert np.allclose(ggp[k], N @ X, atol=1e-12)
                _check_data(h, "stress", name, live, 4)
            else:
                assert conn.shape == (1, 4) and _s(g.attrs["TOPOLOGY"]) == "tri"
                assert _i(g.attrs["NUM_CTRL"]) == 3
                assert np.allclose(gp, [[1 / 3, 1 / 3]], atol=1e-15)
                X = np.asarray([UP2_NODES[n] for n in (2, 5, 3)])
                assert np.allclose(ggp[0], X.mean(axis=0), atol=1e-12)
                _check_data(h, "stress", name, live, 1)


def test_up_h8_order_is_not_brick_order(tmp_path):
    h5py = _h5py()
    out = str(tmp_path / "up3d.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 4)
    ops.nDMaterial("ElasticIsotropic", 1, 3.0e4, 0.2)
    for n, xyz in SKEW3.items():
        ops.node(n, *xyz)
    ops.element("LadrunoUP", 1, *SKEW3.keys(), 1, *UP_OPTS3)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1, 1)
    for n in (5, 6, 7, 8):
        ops.fix(n, 0, 0, 0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(7, 2.0, -1.0, -5.0, 0.0)
    ops.load(6, 0.0, 3.0, 0.0, 0.0)
    ops.recorder("ladruno", out, "-E", "stress")
    _static()
    live = np.asarray(ops.eleResponse(1, "stress"), dtype=float)
    assert live.size == 48
    ops.wipe()

    with h5py.File(_file(out), "r") as h:
        groups = _groups(h, "LadrunoUP")
        assert len(groups) == 1
        name, g = next(iter(groups.items()))
        assert _s(g.attrs["TOPOLOGY"]) == "hex"
        gp = np.asarray(g["QUADRATURE/GP_PARAM"][...], dtype=float)
        assert np.allclose(gp, UP_H8_GP, atol=1e-15)
        X = np.asarray([SKEW3[n] for n in SKEW3], dtype=float)
        ggp = np.asarray(g["GLOBAL_GP_COORDS"][...], dtype=float).reshape(8, 3)
        for L in range(8):
            N = np.prod(1.0 + HEX8_NODES * gp[L], axis=1) / 8.0
            assert np.allclose(ggp[L], N @ X, atol=1e-12)
        _check_data(h, "stress", name, {1: live}, 8)


# -------------------------------------------------- Bernstein GLOBAL_GP_COORDS
def test_curved_beziertri6_global_gp_uses_bernstein_basis(tmp_path):
    """Triangle_6N now has a Lagrange basis (for LST); a BezierTri6 (FAMILY
    bernstein) must still map its GPs with the Bernstein basis. On a curved
    element the two differ, so this pins the family switch."""
    h5py = _h5py()
    out = str(tmp_path / "bt6.ladruno")
    nodes = {1: (0.0, 0.0), 2: (2.0, 0.0), 3: (0.0, 2.0),
             5: (1.0, -0.3), 6: (1.2, 1.2), 7: (-0.2, 1.0)}   # c1 c2 c3 m12 m23 m31
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25)
    for n, xy in nodes.items():
        ops.node(n, *xy)
    ops.element("BezierTri6", 1, 1, 2, 3, 5, 6, 7, 1.0, "PlaneStress", 1)
    ops.fix(1, 1, 1)
    ops.fix(3, 1, 1)
    ops.fix(7, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0, 0.0)
    ops.recorder("ladruno", out, "-N", "displacement")
    _static()
    ops.wipe()
    with h5py.File(_file(out), "r") as h:
        g = next(iter(_groups(h, "BezierTri6").values()))
        assert _s(g.attrs["FAMILY"]) == "bernstein"
        gp = np.asarray(g["QUADRATURE/GP_PARAM"][...], dtype=float)
        X = np.asarray([nodes[n] for n in (1, 2, 3, 5, 6, 7)], dtype=float)
        ggp = np.asarray(g["GLOBAL_GP_COORDS"][...], dtype=float).reshape(gp.shape[0], 2)
        for k in range(gp.shape[0]):
            bern = _t6_bernstein(*gp[k]) @ X
            assert np.allclose(ggp[k], bern, atol=1e-12)
            assert not np.allclose(bern, _t6_lagrange(*gp[k]) @ X, atol=1e-6)
