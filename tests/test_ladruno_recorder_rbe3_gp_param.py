"""WP-226 — `recorder ladruno` MODEL/ELEMENTS gaps found on the YOO Guayaquil mat.

1. A class whose elements carry DIFFERENT node counts (LadrunoDistributingCoupling,
   one footprint per RBE3) used to log "elements with different number of nodes
   exist within the same class tag" once per element and leave every odd-sized
   element out of MODEL/ELEMENTS. Such a class is now written with the schema
   §3.4 CSR layout (ID + CONNECTIVITY_OFFSETS + CONNECTIVITY_NODES, attr
   CONNECTIVITY_LAYOUT="csr", NUM_CTRL=-1), and the log stays clean.
2. LadrunoBrick had no arm in getGeometryAndIntRuleByClassTag, so its group had
   TOPOLOGY custom and no QUADRATURE/GP_PARAM: readers could not place the 8 GPs
   of "stresses". It now shares Brick's Hexahedron_GaussLegendre_2 arm (8 GPs,
   i,j,k lexicographic, zeta fastest) for every -formulation. LadrunoQuad
   (Quad GL2, CCW) and LadrunoCST (Tri GL1, centroid) had the same gap.
"""
import glob
import os

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

A = 1.0 / np.sqrt(3.0)
# Hexahedron_GaussLegendre_2 / Brick / LadrunoBrick: L = 4i + 2j + k over sg = (-A, +A)
HEX8_GP = np.array([[-A, -A, -A], [-A, -A, A], [-A, A, -A], [-A, A, A],
                    [A, -A, -A], [A, -A, A], [A, A, -A], [A, A, A]])
# Brick node natural coordinates (shp3d order)
HEX8_NODES = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                       [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], dtype=float)
# a skewed hex so GLOBAL_GP_COORDS is not a trivial affine map
SKEW = {1: (0.0, 0.0, 0.0), 2: (1.2, 0.1, 0.0), 3: (1.1, 1.3, 0.05), 4: (-0.1, 0.9, 0.0),
        5: (0.05, 0.0, 1.0), 6: (1.0, 0.1, 1.2), 7: (1.2, 1.1, 1.1), 8: (0.0, 1.0, 0.9)}


def _h5py():
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    return pytest.importorskip("h5py")


def _model_groups(h):
    """{group name: h5py.Group} of every MODEL_STAGE[*]/MODEL/ELEMENTS/<name>."""
    out = {}
    for stage in h:
        if not stage.startswith("MODEL_STAGE"):
            continue
        els = h[stage].get("MODEL/ELEMENTS")
        for name in (els if els is not None else ()):
            out[name] = els[name]
    return out


def _file(out):
    files = glob.glob(out) or glob.glob(out + "*")
    assert files, "no .ladruno output produced"
    return files[0]


def _hex(formulation, tag=1, mat=1):
    for n, xyz in SKEW.items():
        ops.node(n, *xyz)
    args = ["LadrunoBrick", tag, *SKEW.keys(), mat]
    if formulation is not None:
        args += ["-formulation", formulation]
    ops.element(*args)


def _static(system="FullGeneral"):
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system(system)
    ops.test("NormDispIncr", 1e-10, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


# ------------------------------------------------------------------ fix 1: RBE3
def test_rbe3_of_different_node_counts_written_csr_and_log_clean(tmp_path, capfd):
    h5py = _h5py()
    out = str(tmp_path / "rbe3.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 3.0e4, 0.2)
    _hex("eas")
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    ops.node(100, 0.6, 0.6, 1.4, "-ndf", 6)     # top-face reference: 4 independents
    ops.node(101, 0.6, 0.0, 0.5, "-ndf", 6)     # side reference: 3 independents
    ops.element("LadrunoDistributingCoupling", 11, 100, 4, 5, 6, 7, 8, "-k", 1.0e8)
    ops.element("LadrunoDistributingCoupling", 12, 101, 3, 1, 2, 6, "-k", 1.0e8)
    want = {11: [100, 5, 6, 7, 8], 12: [101, 1, 2, 6]}
    for t, nodes in want.items():
        assert [int(x) for x in ops.eleNodes(t)] == nodes
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(100, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0)
    ops.load(101, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0)
    capfd.readouterr()
    ops.recorder("ladruno", out, "-N", "displacement", "-E", "stresses")
    _static()
    ops.wipe()
    log = capfd.readouterr()
    text = log.out + log.err
    assert "different number of nodes" not in text, text

    with h5py.File(_file(out), "r") as h:
        groups = _model_groups(h)
        rbe3 = [g for k, g in groups.items() if "LadrunoDistributingCoupling" in k]
        assert len(rbe3) == 1, sorted(groups)
        g = rbe3[0]
        lay = g.attrs["CONNECTIVITY_LAYOUT"]
        assert (lay.decode() if isinstance(lay, bytes) else str(lay)) == "csr"
        assert "CONNECTIVITY" not in g          # pre-§3.4 readers skip, never misread
        assert int(np.atleast_1d(g.attrs["NUM_CTRL"])[0]) == -1
        ids = np.asarray(g["ID"][...]).reshape(-1)
        off = np.asarray(g["CONNECTIVITY_OFFSETS"][...]).reshape(-1)
        nod = np.asarray(g["CONNECTIVITY_NODES"][...]).reshape(-1)
        assert off.shape == (ids.size + 1,) and off[0] == 0 and off[-1] == nod.size
        got = {int(ids[e]): [int(x) for x in nod[off[e]:off[e + 1]]]
               for e in range(ids.size)}
        assert got == want
        assert "GLOBAL_GP_COORDS" not in g
        # the rectangular brick group is untouched by the CSR switch
        brick = [g for k, g in groups.items() if "LadrunoBrick" in k]
        assert len(brick) == 1 and "CONNECTIVITY" in brick[0]
        assert np.asarray(brick[0]["CONNECTIVITY"]).shape == (1, 9)


def test_rbe3_of_one_node_count_keeps_rectangular_connectivity(tmp_path):
    """The CSR layout is data-driven: a class whose elements all share one node
    count is written exactly as before (rectangular CONNECTIVITY)."""
    h5py = _h5py()
    out = str(tmp_path / "rbe3_same.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 3.0e4, 0.2)
    _hex("std")
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    ops.node(100, 0.6, 0.6, 1.4, "-ndf", 6)
    ops.node(101, 0.6, 0.0, 0.5, "-ndf", 6)
    ops.element("LadrunoDistributingCoupling", 11, 100, 3, 5, 6, 7, "-k", 1.0e8)
    ops.element("LadrunoDistributingCoupling", 12, 101, 3, 1, 2, 6, "-k", 1.0e8)
    ops.recorder("ladruno", out, "-N", "displacement")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(100, 0.0, 0.0, -10.0, 0.0, 0.0, 0.0)
    _static()
    ops.wipe()
    with h5py.File(_file(out), "r") as h:
        g = [g for k, g in _model_groups(h).items() if "LadrunoDistributingCoupling" in k][0]
        assert "CONNECTIVITY_LAYOUT" not in g.attrs
        conn = np.asarray(g["CONNECTIVITY"][...])
        assert conn.tolist() == [[11, 100, 5, 6, 7], [12, 101, 1, 2, 6]]


# ------------------------------------------------------- fix 2: LadrunoBrick GP_PARAM
@pytest.mark.parametrize("formulation", ["eas", "std", "bbar"])
def test_ladrunobrick_writes_hex_gl2_gp_param(tmp_path, formulation):
    h5py = _h5py()
    out = str(tmp_path / f"brick_{formulation}.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 3.0e4, 0.2)
    _hex(formulation)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(7, 2.0, -1.0, -5.0)               # unsymmetric: GP stresses all differ
    ops.load(6, 0.0, 3.0, 0.0)
    ops.recorder("ladruno", out, "-E", "stresses")
    _static()
    live = np.asarray(ops.eleResponse(1, "stresses"), dtype=float)
    assert live.size == 48
    ops.wipe()

    with h5py.File(_file(out), "r") as h:
        groups = _model_groups(h)
        keys = [k for k in groups if "-LadrunoBrick[" in k]
        assert len(keys) == 1, sorted(groups)
        g = groups[keys[0]]
        topo = g.attrs["TOPOLOGY"]
        assert (topo.decode() if isinstance(topo, bytes) else str(topo)) == "hex"
        assert int(np.atleast_1d(g.attrs["NUM_GP"])[0]) == 8
        assert int(np.atleast_1d(g.attrs["NDIR"])[0]) == 3
        gp = np.asarray(g["QUADRATURE/GP_PARAM"][...], dtype=float)
        assert gp.shape == (8, 3)
        assert np.allclose(gp, HEX8_GP, atol=1e-15)
        w = np.asarray(g["QUADRATURE/GP_WEIGHT"][...], dtype=float).reshape(-1)
        assert np.allclose(w, 1.0) and w.size == 8

        # GLOBAL_GP_COORDS == trilinear x(xi_L) on the skewed hex
        X = np.asarray([SKEW[n] for n in SKEW], dtype=float)
        ggp = np.asarray(g["GLOBAL_GP_COORDS"][...], dtype=float).reshape(8, 3)
        for L in range(8):
            N = np.prod(1.0 + HEX8_NODES * gp[L], axis=1) / 8.0
            assert np.allclose(ggp[L], N @ X, atol=1e-12)

        # and the stresses bucket pairs GAUSS_ID 0..7 with those stations
        on_el = h[[s for s in h if s.startswith("MODEL_STAGE")][0]]["RESULTS/ON_ELEMENTS/stresses"]
        b = [on_el[k] for k in on_el if "LadrunoBrick" in k][0]
        assert list(np.asarray(b["COLUMN_MAP/GAUSS_ID"]).reshape(-1)) == list(range(8))
        data = np.asarray(b["DATA"][...], dtype=float)
        assert data[-1].reshape(-1) == pytest.approx(live, rel=1e-10, abs=1e-12)


@pytest.mark.parametrize("cls,nodes,gp", [
    ("LadrunoQuad", {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.1, 1.0), 4: (0.0, 0.9)},
     np.array([[-A, -A], [A, -A], [A, A], [-A, A]])),
    ("LadrunoCST", {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (0.2, 1.0)},
     np.array([[1.0 / 3.0, 1.0 / 3.0]])),
])
def test_ladruno_plane_elements_write_gp_param(tmp_path, cls, nodes, gp):
    h5py = _h5py()
    out = str(tmp_path / f"{cls}.ladruno")
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 3.0e4, 0.2)
    for n, xy in nodes.items():
        ops.node(n, *xy)
    ops.element(cls, 1, *nodes.keys(), 1, "-type", "PlaneStrain")
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 1.0, 0.0)
    ops.recorder("ladruno", out, "-N", "displacement")
    _static()
    ops.wipe()
    with h5py.File(_file(out), "r") as h:
        groups = _model_groups(h)
        keys = [k for k in groups if f"-{cls}[" in k]
        assert len(keys) == 1, sorted(groups)
        got = np.asarray(groups[keys[0]]["QUADRATURE/GP_PARAM"][...], dtype=float)
        assert got.shape == gp.shape and np.allclose(got, gp, atol=1e-15)
