"""WP-234: `recorder ladruno ... -R <eleRegion> -RN <nodeRegion>`.

A MeshRegion cannot hold an explicit node list and an element list together:
MeshRegion::setNodes and ::setElements each delete both lists, and `-ele`
rebuilds the nodes from the elements (OpenSees only warns "both elements & nodes
set, ONLY set using elements"). So one `-R` region could never say "these
elements plus these extra nodes" -- embedded rebar nodes whose minted truss
elements must stay out of the file, for example (apeGmsh example 02). `-RN`
gives the node channel its own region: MODEL/NODES and the -N rows hold the -R
region's nodes plus the -RN nodes, while -R alone still restricts MODEL/ELEMENTS
and the -E rows.
"""
import os

import numpy as np
import pytest

from _testbed import ops

h5py = pytest.importorskip("h5py")
pytestmark = [pytest.mark.zone_a]

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

# 1 -- 2 -- 3 -- 4 along x; elements 1 (1-2), 2 (2-3), 3 (3-4).
# Region 1 = element 1 (its nodes 1, 2). Region 2 = node 4 only.
MODEL = """
ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.uniaxialMaterial("Elastic", 1, 1.0e6)
ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
for n in (2, 3, 4):
    ops.node(n, float(n - 1), 0.0); ops.fix(n, 0, 1)
ops.element("Truss", 1, 1, 2, 1.0, 1)
ops.element("Truss", 2, 2, 3, 1.0, 1)
ops.element("Truss", 3, 3, 4, 1.0, 1)
ops.region(1, "-ele", 1)
ops.region(2, "-node", 4)
"""

STATIC = """
ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(4, 1.0, 0.0)
ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
ops.test("NormDispIncr", 1.0e-12, 10); ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0); ops.analysis("Static")
"""


def _stage(f):
    return f[[k for k in f if k.startswith("MODEL_STAGE")][0]]


def _record(out, *args):
    exec(MODEL, {"ops": ops})
    ok = ops.recorder("ladruno", out, *args)
    exec(STATIC, {"ops": ops})
    ops.analyze(2)
    ops.wipe()
    return ok


def _elem_ids(stage):
    ids = []
    for obj in stage["MODEL/ELEMENTS"].values():
        if isinstance(obj, h5py.Dataset):            # rows: (elemTag, nodes...)
            ids.extend(int(r[0]) for r in np.asarray(obj))
        elif "ID" in obj:                            # CSR group
            ids.extend(int(v) for v in np.asarray(obj["ID"]).ravel())
        else:
            ids.extend(int(r[0]) for r in np.asarray(obj["CONNECTIVITY"]))
    return sorted(ids)


@pytest.mark.parametrize("order", ["R_first", "RN_first"])
def test_rn_adds_nodes_without_adding_elements(tmp_path, order):
    out = str(tmp_path / "rn.ladruno")
    flags = ["-R", 1, "-RN", 2] if order == "R_first" else ["-RN", 2, "-R", 1]
    _record(out, "-N", "displacement", "-E", "force", *flags)
    with h5py.File(out, "r") as f:
        g = _stage(f)
        assert sorted(np.asarray(g["MODEL/NODES/ID"]).ravel().tolist()) == [1, 2, 4]
        disp = g["RESULTS/ON_NODES/DISPLACEMENT"]
        assert sorted(np.asarray(disp["ID"]).ravel().tolist()) == [1, 2, 4]
        data = np.asarray(disp["DATA"])
        assert data.shape == (2, 3, 2)                          # T x nodes x dof
        # node 4 is the loaded free end: u = 3 P L / EA, P = 2 after two steps
        row = np.asarray(disp["ID"]).ravel().tolist().index(4)
        assert data[-1, row, 0] == pytest.approx(6.0e-6, rel=1e-9)
        assert _elem_ids(g) == [1]                              # 2 and 3 stay out
        for bucket in g["RESULTS/ON_ELEMENTS/force"].values():
            assert np.asarray(bucket["ID"]).ravel().tolist() == [1]
        sets = g["MODEL/SETS"]
        assert sorted(sets) == ["SET_1", "SET_2"]
        assert np.asarray(sets["SET_2/NODES"]).ravel().tolist() == [4]
        # a -node region has no elements; the SETS writer skips an empty dataset
        assert ("ELEMENTS" not in sets["SET_2"]
                or np.asarray(sets["SET_2/ELEMENTS"]).size == 0)


def test_plain_R_is_unchanged(tmp_path):
    out = str(tmp_path / "r.ladruno")
    _record(out, "-N", "displacement", "-E", "force", "-R", 1)
    with h5py.File(out, "r") as f:
        g = _stage(f)
        assert sorted(np.asarray(g["MODEL/NODES/ID"]).ravel().tolist()) == [1, 2]
        assert _elem_ids(g) == [1]
        assert sorted(g["MODEL/SETS"]) == ["SET_1"]


def test_rn_without_R_is_refused_by_name(tmp_path, capfd):
    out = str(tmp_path / "rn_alone.ladruno")
    exec(MODEL, {"ops": ops})
    try:
        ops.recorder("ladruno", out, "-N", "displacement", "-RN", 2)
    except Exception:
        pass                         # openseespy raises on a refused command
    ops.wipe()
    log = "".join(capfd.readouterr())
    assert "-RN <nodeRegionTag> needs -R" in log, log


def test_rn_unknown_region_is_refused(tmp_path, capfd):
    out = str(tmp_path / "rn_bad.ladruno")
    exec(MODEL, {"ops": ops})
    try:
        ops.recorder("ladruno", out, "-N", "displacement", "-R", 1, "-RN", 99)
    except Exception:
        pass
    ops.wipe()
    assert "-RN region 99 is null" in "".join(capfd.readouterr())


def test_capability_is_listed():
    """apeGmsh gates its two-region emit on this name (absence = older build)."""
    import json
    caps = json.loads(ops.ladrunoSchema())["capabilities"]
    assert "recorderLadrunoNodeRegion" in caps, caps
