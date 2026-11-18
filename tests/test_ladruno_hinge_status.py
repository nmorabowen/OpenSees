"""`hingeStatus`: one element response for hinge presence and damage (WP-235).

`LadrunoForceBeamColumn`, `LadrunoDispBeamColumn` and `LadrunoIMKBeam`, 2-D and
3-D, all answer the same token with the same 8 columns:

    kindZ_i, damageZ_i, kindZ_j, damageZ_j, kindY_i, damageY_i, kindY_j, damageY_j

kind: 0 absent, 1 uniaxial law, 2 biaxial law (one law, both slots of that end).
damage: the law's own `damage` response, NaN where it is not defined.

The oracle for every damage value is the per-material token that already exists
(`hingeI damage`, `hingeBiaxialJ damage`, `hinge 2 damage`, ...), and the two
laws in each model differ, so a swapped slot or end cannot pass.
"""
import glob
import math
import os

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

COLUMNS = [
    "kindZ_i", "damageZ_i", "kindZ_j", "damageZ_j",
    "kindY_i", "damageY_i", "kindY_j", "damageY_j",
]
SLOTS = ["Zi", "Zj", "Yi", "Yj"]

E, A, Iz, Iy, G, Jx, L = 2.0e8, 0.01, 1.0e-4, 5.0e-5, 8.0e7, 1.0e-4, 3.0
THETA = 0.02
NSTEPS = 4

# material tags
COH_A, COH_B, STEEL, FATIGUE, BI_A, BI_B = 10, 11, 12, 13, 20, 21

FAMILIES = ["LadrunoForceBeamColumn", "LadrunoDispBeamColumn", "LadrunoIMKBeam"]
ALL = [(family, ndm) for family in FAMILIES for ndm in (2, 3)]
ALL_IDS = [f"{family}-{ndm}d" for family, ndm in ALL]
ALL_3D = [(family, 3) for family in FAMILIES]
ALL_3D_IDS = [f"{family}-3d" for family, _ in ALL_3D]


def _is_imk(family):
    return family == "LadrunoIMKBeam"


def _hinge_args(family, ndm, uniaxial, biaxial):
    """Element options for {slot: matTag} and {'i'|'j': ndMatTag}."""
    args = []
    if _is_imk(family):
        if ndm == 2:
            # The 2-D IMK beam takes one law and the ends that carry it.
            tags = set(uniaxial.values())
            assert len(tags) <= 1, "2-D LadrunoIMKBeam has one -matZ"
            if uniaxial:
                ends = {"Zi": "i", "Zj": "j"}
                where = "both" if len(uniaxial) == 2 else ends[next(iter(uniaxial))]
                args += ["-hinge", where, "-matZ", tags.pop()]
            return args
        flag = {"Zi": "-matZi", "Zj": "-matZj", "Yi": "-matYi", "Yj": "-matYj"}
    else:
        flag = {"Zi": "-hingeI", "Zj": "-hingeJ", "Yi": "-hingeYI", "Yj": "-hingeYJ"}
    for slot, tag in uniaxial.items():
        args += [flag[slot], tag]
    for end, tag in biaxial.items():
        args += [{"i": "-hingeBiaxialI", "j": "-hingeBiaxialJ"}[end], tag]
    return args


def _build(family, ndm, uniaxial=None, biaxial=None):
    uniaxial = uniaxial or {}
    biaxial = biaxial or {}
    ops.wipe()
    if ndm == 2:
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.geomTransf("Linear", 1)
        ops.section("Elastic", 1, E, A, Iz)
    else:
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.section("Elastic", 1, E, A, Iz, Iy, G, Jx)
    ops.beamIntegration("Lobatto", 1, 1, 5)

    ops.uniaxialMaterial("LadrunoCohesiveHinge", COH_A, 100.0, 50.0, "-linear")
    ops.uniaxialMaterial("LadrunoCohesiveHinge", COH_B, 150.0, 50.0, "-linear")
    ops.uniaxialMaterial("Steel01", STEEL, 1.0e6, 1.0e8, 0.02)   # stays elastic
    ops.uniaxialMaterial("Fatigue", FATIGUE, STEEL)
    ops.nDMaterial("LadrunoCohesiveHingeBiaxial", BI_A, 100.0, 50.0, 80.0, 50.0, "-linear")
    ops.nDMaterial("LadrunoCohesiveHingeBiaxial", BI_B, 150.0, 50.0, 120.0, 50.0, "-linear")

    hinges = _hinge_args(family, ndm, uniaxial, biaxial)
    if _is_imk(family):
        props = [A, E, Iz] if ndm == 2 else [A, E, G, Jx, Iy, Iz]
        ops.element(family, 1, 1, 2, *props, 1, *hinges)
    else:
        ops.element(family, 1, 1, 2, 1, 1, *hinges)


def _analyze(ndm, steps=NSTEPS, each_step=None):
    """Fixed-fixed beam, far end rotated about every bending axis."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    if ndm == 2:
        prescribed = [0.0, 0.0, THETA]
    else:
        prescribed = [0.0, 0.0, 0.0, 0.0, THETA, THETA]
    for dof, value in enumerate(prescribed, start=1):
        ops.sp(2, dof, value)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0 / steps)
    ops.test("NormDispIncr", 1.0e-10, 100)
    ops.algorithm("Newton")
    ops.analysis("Static")
    for _ in range(steps):
        assert ops.analyze(1) == 0
        if each_step is not None:
            each_step()


def _status():
    values = ops.eleResponse(1, "hingeStatus")
    assert len(values) == len(COLUMNS)
    return {
        slot: (values[2 * n], values[2 * n + 1]) for n, slot in enumerate(SLOTS)
    }


def _uniaxial_token(family, slot):
    if _is_imk(family):
        return ["hinge", str(SLOTS.index(slot))]
    return [{"Zi": "hingeI", "Zj": "hingeJ", "Yi": "hingeYI", "Yj": "hingeYJ"}[slot]]


def _material_damage(family, slot, biaxial_end=None):
    if biaxial_end is not None:
        token = ["hingeBiaxialI" if biaxial_end == "i" else "hingeBiaxialJ"]
    else:
        token = _uniaxial_token(family, slot)
    return ops.eleResponse(1, *token, "damage")[0]


def _assert_absent(status, slots):
    for slot in slots:
        kind, damage = status[slot]
        assert kind == 0, slot
        assert math.isnan(damage), slot


@pytest.mark.parametrize("family,ndm", ALL, ids=ALL_IDS)
def test_no_hinge_reads_absent_and_not_zero(family, ndm):
    _build(family, ndm)
    _assert_absent(_status(), SLOTS)
    _analyze(ndm)
    _assert_absent(_status(), SLOTS)
    ops.wipe()


@pytest.mark.parametrize("family,ndm", ALL, ids=ALL_IDS)
def test_unloaded_hinge_differs_from_absent(family, ndm):
    _build(family, ndm, uniaxial={"Zj": COH_A})
    status = _status()
    assert status["Zj"] == (1.0, 0.0)
    _assert_absent(status, ["Zi", "Yi", "Yj"])
    ops.wipe()


@pytest.mark.parametrize("family,ndm", ALL, ids=ALL_IDS)
def test_uniaxial_damage_matches_material_token(family, ndm):
    if ndm == 3:
        hinges = {"Zi": COH_A, "Zj": COH_B, "Yi": COH_B, "Yj": COH_A}
    elif _is_imk(family):
        hinges = {"Zi": COH_A, "Zj": COH_A}   # one law on the 2-D IMK beam
    else:
        hinges = {"Zi": COH_A, "Zj": COH_B}
    _build(family, ndm, uniaxial=hinges)
    _analyze(ndm)

    status = _status()
    damages = []
    for slot in hinges:
        kind, damage = status[slot]
        assert kind == 1, slot
        assert damage == _material_damage(family, slot), slot
        damages.append(damage)
    assert max(damages) > 0.0, "the load must damage at least one hinge"
    _assert_absent(status, [slot for slot in SLOTS if slot not in hinges])
    ops.wipe()


@pytest.mark.parametrize("family,ndm", ALL, ids=ALL_IDS)
def test_material_without_damage_reads_nan(family, ndm):
    _build(family, ndm, uniaxial={"Zj": STEEL})
    _analyze(ndm)
    status = _status()
    kind, damage = status["Zj"]
    assert kind == 1
    assert math.isnan(damage), "Steel01 defines no damage; 0 would be a fake value"
    _assert_absent(status, ["Zi", "Yi", "Yj"])
    ops.wipe()


@pytest.mark.parametrize("family,ndm", ALL, ids=ALL_IDS)
def test_any_material_defining_damage_is_reported(family, ndm):
    """Damage is not special-cased to the cohesive laws."""
    _build(family, ndm, uniaxial={"Zj": FATIGUE})
    _analyze(ndm)
    kind, damage = _status()["Zj"]
    assert kind == 1
    assert not math.isnan(damage)
    assert damage == _material_damage(family, "Zj")
    ops.wipe()


@pytest.mark.parametrize("family,ndm", ALL_3D, ids=ALL_3D_IDS)
def test_biaxial_law_fills_both_slots_of_its_end(family, ndm):
    _build(family, ndm, biaxial={"i": BI_A, "j": BI_B})
    assert _status()["Zi"] == (2.0, 0.0)
    _analyze(ndm)

    status = _status()
    for end, z_slot, y_slot in (("i", "Zi", "Yi"), ("j", "Zj", "Yj")):
        expected = _material_damage(family, z_slot, biaxial_end=end)
        assert status[z_slot] == (2.0, expected)
        assert status[y_slot] == (2.0, expected)
    assert max(status["Zi"][1], status["Zj"][1]) > 0.0
    ops.wipe()


@pytest.mark.parametrize("family,ndm", ALL_3D, ids=ALL_3D_IDS)
def test_biaxial_at_one_end_leaves_the_other_absent(family, ndm):
    _build(family, ndm, biaxial={"j": BI_A})
    _analyze(ndm)
    status = _status()
    assert status["Zj"][0] == 2.0 and status["Yj"][0] == 2.0
    _assert_absent(status, ["Zi", "Yi"])
    ops.wipe()


def _read_recorded(h5py, path):
    with h5py.File(path, "r") as handle:
        stage = next(key for key in handle if key.startswith("MODEL_STAGE"))
        groups = handle[f"{stage}/RESULTS/ON_ELEMENTS/hingeStatus"]
        assert len(groups) == 1
        group = groups[next(iter(groups))]
        names = group["COLUMN_MAP"].attrs["COMP_NAMES"]
        if isinstance(names, np.ndarray):
            names = names.flat[0]
        if isinstance(names, bytes):
            names = names.decode()
        return str(names), np.array(group["DATA"][...])


RECORDED = [
    ("none", {}, {}),
    ("cohesive_and_steel", {"Zi": STEEL, "Zj": COH_A}, {}),
    ("biaxial", {}, {"i": BI_A, "j": BI_B}),
]


def _recorded_cases():
    for family, ndm in ALL:
        for name, uniaxial, biaxial in RECORDED:
            if ndm == 2 and biaxial:
                continue
            if ndm == 2 and _is_imk(family) and len(set(uniaxial.values())) > 1:
                uniaxial = {"Zj": COH_A}   # one law on the 2-D IMK beam
            yield pytest.param(
                family, ndm, uniaxial, biaxial, id=f"{family}-{ndm}d-{name}"
            )


@pytest.mark.parametrize("family,ndm,uniaxial,biaxial", list(_recorded_cases()))
def test_ladruno_recorder_columns_and_values(tmp_path, family, ndm, uniaxial, biaxial):
    h5py = pytest.importorskip("h5py")
    path = str(tmp_path / "hinge_status.ladruno")

    _build(family, ndm, uniaxial=uniaxial, biaxial=biaxial)
    ops.recorder("ladruno", path, "-E", "hingeStatus", "-flush", 0, "-compress", 0)
    live = []
    _analyze(ndm, each_step=lambda: live.append(ops.eleResponse(1, "hingeStatus")))
    ops.wipe()

    files = glob.glob(str(tmp_path / "*"))
    assert files, "the recorder wrote no file"
    names, data = _read_recorded(h5py, files[0])
    assert names == ",".join(COLUMNS)

    assert data.shape == (NSTEPS, 1, len(COLUMNS))
    np.testing.assert_array_equal(data[:, 0, :], np.array(live))  # NaN == NaN here
