"""WP-232a: the MPM contracts -- fixtures, validators, malformed variants.

Pure h5py + numpy, no OpenSees build (run by the static-gates job of ladruno.yml, which
installs h5py: the import is hard on purpose, so a missing h5py fails the gate instead of
skipping it green). The contracts are Ladruno_implementation/232a_mpm_contracts.md; the
validators are tests/mpm_contracts.py (Contract 1) and
Ladruno_scripts/ladruno_recorder_tests/ladruno_format.py (Contract 2, through
mpm_contracts.validate_particle_results). Every rule code the page names has a refusal
case here (test_every_rule_code_has_a_refusal_case).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mpm_contract_fixtures as fx  # noqa: E402
import mpm_contracts as mc  # noqa: E402

lf = mc.lf
DOC = Path(__file__).resolve().parents[1] / "Ladruno_implementation" / "232a_mpm_contracts.md"
STAGE = "MODEL_STAGE[1]"
OP = f"{STAGE}/RESULTS/ON_PARTICLES"
MP = f"{STAGE}/MODEL/PARTICLES"


# --------------------------------------------------------------------------- helpers


def _content(path) -> dict:
    """Every object's dtype / shape / values / filters and every attribute."""
    out = {}

    def norm(v):
        v = np.asarray(v)
        if v.dtype.kind == "O":
            v = np.array([x.decode() if isinstance(x, bytes) else x for x in v.reshape(-1)])
        return v

    with h5py.File(path, "r") as f:
        def visit(name, obj):
            out[name + "@attrs"] = {k: norm(obj.attrs[k]) for k in obj.attrs}
            if isinstance(obj, h5py.Dataset):
                out[name] = (str(obj.dtype), obj.shape, norm(obj[...]), lf.filters_of(obj))
        out["/@attrs"] = {k: norm(f.attrs[k]) for k in f.attrs}
        f.visititems(visit)
    return out


def _same(a, b) -> bool:
    if isinstance(a, dict):
        return a.keys() == b.keys() and all(_same(a[k], b[k]) for k in a)
    if isinstance(a, tuple):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    if isinstance(a, np.ndarray):
        if a.dtype.kind == "f" and b.dtype.kind == "f":
            return a.shape == b.shape and np.array_equal(a, b, equal_nan=True)
        return a.shape == b.shape and bool(np.all(a == b))
    return a == b


def _codes(problems) -> set[str]:
    return {p.split(":", 1)[0] for p in problems}


# --------------------------------------------------------------------------- the fixtures


def test_committed_fixtures_are_conformant():
    warnings: list[str] = []
    assert mc.validate_particle_set(str(fx.PARTICLE_SET), warnings) == []
    assert warnings == []                      # lattice seeding: VOLUME == GIMP domain
    assert mc.validate_particle_results(str(fx.PARTICLE_RESULTS)) == []


def test_committed_fixtures_match_the_generator(tmp_path):
    """The committed files are what tests/mpm_contract_fixtures.py writes today."""
    fx.write_contract1(tmp_path / "c1.h5")
    fx.write_contract2(tmp_path / "c2.ladruno")
    fx.write_mixed(tmp_path / "mixed.ladruno")
    assert _same(_content(fx.PARTICLE_SET), _content(tmp_path / "c1.h5"))
    assert _same(_content(fx.PARTICLE_RESULTS), _content(tmp_path / "c2.ladruno"))
    assert _same(_content(fx.PARTICLE_MIXED), _content(tmp_path / "mixed.ladruno"))


def test_mixed_fixture_fe_stage_then_particle_stage():
    """232a §4: an FE stage today's readers open, then a particle stage they ignore."""
    assert mc.validate_particle_results(str(fx.PARTICLE_MIXED)) == []
    assert lf.round_trip_oracle(str(fx.PARTICLE_MIXED)) == []
    with lf.LadrunoReader(str(fx.PARTICLE_MIXED)) as r:
        assert int(r.info()["FORMAT_VERSION"]) == 1
        assert int(r.info()["FORMAT_MINOR_VERSION"]) == 2
        assert sorted(r.stages(), key=lf._stage_id) == ["MODEL_STAGE[0]", "MODEL_STAGE[1]"]
        assert r.particle_stages() == ["MODEL_STAGE[1]"]
        fe = r.stage_attrs("MODEL_STAGE[0]")
        assert fe["KIND"] == "static" and "EMPTY_PARTITION" not in fe
        assert r.nodes("MODEL_STAGE[0]")[0].shape[0] > 0
        assert len(r.element_groups("MODEL_STAGE[0]")) > 0
        ps = r.stage_attrs("MODEL_STAGE[1]")
        assert ps["KIND"] == "transient" and int(ps["EMPTY_PARTITION"]) == 1
        ids, _ = r.particle_output("MODEL_STAGE[1]", "POSITION", 0)
        assert ids.tolist() == list(range(12))


@pytest.mark.parametrize("mode,kind", sorted(lf.MODE_KIND.items()))
def test_particle_stage_kind_is_pinned_from_mode(tmp_path, mode, kind):
    """232a §2.3: transient for explicit and implicit-dynamic, static for quasi-static."""
    p = tmp_path / "c2.ladruno"
    fx.write_contract2(p, mode=mode)
    with lf.LadrunoReader(str(p)) as r:
        assert r.stage_attrs(STAGE)["KIND"] == kind
    assert mc.validate_particle_results(str(p)) == []


def test_contract2_is_an_additive_minor_not_a_major():
    """FORMAT_VERSION stays 1; the minor is 2 (232a §2.1)."""
    with lf.LadrunoReader(str(fx.PARTICLE_RESULTS)) as r:
        info = r.info()
    assert int(info["FORMAT_VERSION"]) == lf.FORMAT_VERSION == 1
    assert int(info["FORMAT_MINOR_VERSION"]) == lf.FORMAT_MINOR_VERSION_PARTICLES == 2


def test_contract2_stage_is_the_recorders_empty_stage():
    """232a §2.3: EMPTY_PARTITION = 1, NODES/ID int32 [0 x 1], empty FE groups."""
    with h5py.File(fx.PARTICLE_RESULTS, "r") as f:
        st = f[STAGE]
        assert int(st.attrs["EMPTY_PARTITION"][0]) == 1
        assert int(st.attrs["PARTICLE_EPOCH"][0]) == 0
        assert st["MODEL/NODES/ID"].shape == (0, 1)
        assert st["MODEL/NODES/ID"].dtype == np.dtype("<i4")
        assert len(st["MODEL/ELEMENTS"]) == 0
        for g in lf.EMPTY_RESULT_GROUPS:
            assert len(st[f"RESULTS/{g}"]) == 0


def test_contract2_fixture_uses_the_fork_attribute_encodings():
    """1-element int32 arrays; fixed-length NULLTERM strings of exactly len bytes;
    empty strings absent (232a §0, Ladruno_Hdf5.h writei / writes)."""
    with h5py.File(fx.PARTICLE_RESULTS, "r") as f:
        for name in ("FORMAT_VERSION", "FORMAT_MINOR_VERSION"):
            v = f["INFO"].attrs[name]
            assert v.shape == (1,) and v.dtype == np.dtype("<i4")
        kind = f[STAGE].attrs["KIND"]
        assert kind.shape == (1,) and kind.dtype.kind == "S" and kind[0] == b"transient"
        t = f[STAGE].attrs.get_id("KIND").get_type()
        assert t.get_strpad() == h5py.h5t.STR_NULLTERM and t.get_size() == len(b"transient")
        fg = f[f"{OP}/DEFORMATION_GRADIENT"]
        assert "DIMENSION" not in fg.attrs and "DESCRIPTION" not in fg.attrs
        assert f[f"{MP}/BODIES/NAME"].dtype.kind == "S"


def test_round_trip_contract1_to_contract2():
    """232a §3: what the fork carries from the particle set into the results file."""
    with h5py.File(fx.PARTICLE_SET, "r") as c1, lf.LadrunoReader(str(fx.PARTICLE_RESULTS)) as r:
        stage = r.particle_stages()[0]
        mp = r.f[f"{stage}/MODEL/PARTICLES"]
        assert lf._attr(mp, "SOURCE_SET_ID") == c1.attrs["SET_ID"]
        assert int(r.info()["SPATIAL_DIM"]) == int(c1.attrs["NDM"])
        for k in ("UNIT_LENGTH", "UNIT_FORCE", "UNIT_TIME", "UNIT_SYSTEM"):
            assert lf._attr(mp, k) == c1.attrs[k]
        assert lf._attr(mp, "PLANE_STRAIN_THICKNESS") == c1.attrs["PLANE_STRAIN_THICKNESS"]
        for k in ("ORIGIN", "CELL_SIZE", "NUM_CELLS", "BASIS"):
            assert np.all(np.asarray(lf._attr(mp["GRID"], k)) == np.asarray(c1["GRID"].attrs[k]))
        assert mp["BODIES/ID"][...].tolist() == c1["BODIES/ID"][...].tolist()
        assert lf.str_dataset(mp["BODIES/NAME"]) == lf.str_dataset(c1["BODIES/NAME"])
        P = c1["PARTICLES"]
        for k in ("ID", "BODY_ID", "MATERIAL_TAG"):
            assert mp[k][...].tolist() == P[k][...].tolist()
        assert np.allclose(mp["MASS"][...], P["DENSITY"][...] * P["VOLUME"][...])
        ids, x0 = r.particle_output(stage, "POSITION", 0)
        assert ids.tolist() == P["ID"][...].tolist()
        assert np.array_equal(x0, P["POSITION"][...])
        _, s0 = r.particle_output(stage, "STRESS", 0)
        assert s0.shape[1] == 3 and np.array_equal(s0, P["STRESS"][...])   # [11, 22, 12]
        _, v0 = r.particle_output(stage, "VOLUME", 0)
        assert np.array_equal(v0[:, 0], P["VOLUME"][...])
        assert "STRESS33" not in r.f[f"{stage}/RESULTS/ON_PARTICLES"]   # no v1 material exposes it


def test_contract2_removal_and_state_nan():
    with lf.LadrunoReader(str(fx.PARTICLE_RESULTS)) as r:
        stage = r.particle_stages()[0]
        ids1, _ = r.particle_output(stage, "POSITION", 1)
        ids2, _ = r.particle_output(stage, "POSITION", 2)
        _, eps = r.particle_output(stage, "STATE_eqPlasticStrain", 1)
        body = r.f[f"{stage}/MODEL/PARTICLES/BODY_ID"][...]
    assert 11 in ids1.tolist() and 11 not in ids2.tolist()
    assert np.all(np.isnan(eps[body[ids1] == 1])) and np.all(np.isfinite(eps[body[ids1] == 0]))


def test_contract1_accepts_the_cpp_attribute_encodings(tmp_path):
    """A reader accepts 1-element arrays and fixed-length strings too (232a §0)."""
    p = tmp_path / "c1.h5"
    fx.write_contract1(p)
    with h5py.File(p, "r+") as f:
        for k in list(f.attrs):
            v = f.attrs[k]
            del f.attrs[k]
            if isinstance(v, str):
                b = v.encode()
                f.attrs[k] = np.array([b], dtype=f"S{len(b)}")
            else:
                f.attrs[k] = np.atleast_1d(v)
    assert mc.validate_particle_set(str(p)) == []


# --------------------------------------------------------------------------- valid variants


def test_contract1_mass_source_and_optional_fields(tmp_path):
    m = fx.deep_copy_model()
    p = m["particles"]
    p["MASS"] = p.pop("DENSITY") * p["VOLUME"]
    m["attrs"]["MASS_SOURCE"] = "mass"
    p["VELOCITY"] = np.zeros_like(p["POSITION"])
    p["FLAGS"] = np.zeros(len(p["ID"]), dtype="<u4")
    del p["STRESS33"]
    fx.write_contract1(tmp_path / "c1.h5", m)
    assert mc.validate_particle_set(str(tmp_path / "c1.h5")) == []


def _three_d(m):
    a, p = m["attrs"], m["particles"]
    a["NDM"] = np.int32(3)
    a["VOIGT_ORDER"] = mc.VOIGT[3]
    del a["PLANE_STRAIN_THICKNESS"]
    m["grid"] = {"ORIGIN": np.zeros(3), "CELL_SIZE": np.full(3, 0.5),
                 "NUM_CELLS": np.array([4, 4, 4], dtype=np.int32), "BASIS": "gimp"}
    n = len(p["ID"])
    p["POSITION"] = np.c_[p["POSITION"], np.full(n, 0.125)]
    p["GIMP_HALF_LENGTH"] = np.full((n, 3), 0.125)
    p["VOLUME"] = np.full(n, 0.25 ** 3)
    p["STRESS"] = np.zeros((n, 6))
    del p["STRESS33"]


def test_contract1_three_dimensional_set(tmp_path):
    m = fx.deep_copy_model()
    _three_d(m)
    fx.write_contract1(tmp_path / "c1.h5", m)
    assert mc.validate_particle_set(str(tmp_path / "c1.h5")) == []


def test_contract2_accepts_an_added_particle(tmp_path):
    """Addition: a new id above every used one, appended to the registry (232a §2.7)."""
    m = fx.deep_copy_model()
    outs = fx.contract2_outputs(m)
    p = m["particles"]
    for k in ("ID", "BODY_ID", "MATERIAL_TAG", "VOLUME", "DENSITY", "STRESS33"):
        p[k] = np.r_[p[k], p[k][-1:] + (1 if k == "ID" else 0)]
    for k in ("POSITION", "GIMP_HALF_LENGTH", "STRESS"):
        p[k] = np.vstack([p[k], p[k][-1:]])
    last = outs[-1]
    last["ids"] = np.r_[last["ids"], 12]
    last["fields"] = {k: np.vstack([v, v[-1:]]) for k, v in last["fields"].items()}
    fx.write_contract2(tmp_path / "c2.ladruno", m, outs)
    assert mc.validate_particle_results(str(tmp_path / "c2.ladruno")) == []


def test_contract2_time_may_restart_and_source_is_optional(tmp_path):
    """TIME is only finite (`loadConst -time 0`, `setTime`); SOURCE_SET_ID is absent when
    no contract-1 file was loaded (232a §2.4, §2.5)."""
    outs = fx.contract2_outputs()
    outs[2]["time"] = 0.0
    fx.write_contract2(tmp_path / "c2.ladruno", outputs=outs, source=False)
    assert mc.validate_particle_results(str(tmp_path / "c2.ladruno")) == []


def test_contract2_quasi_static_mode_runs_under_a_static_stage(tmp_path):
    fx.write_contract2(tmp_path / "c2.ladruno", mode="quasi-static")
    assert mc.validate_particle_results(str(tmp_path / "c2.ladruno")) == []


def _second_stage(f, epoch, ids_back):
    """Copy MODEL_STAGE[1] to MODEL_STAGE[2] with PARTICLE_EPOCH = epoch; its outputs carry
    `ids_back` again (particle 11 was removed in stage 1)."""
    f.copy(STAGE, "MODEL_STAGE[2]")
    st = f["MODEL_STAGE[2]"]
    st.attrs["PARTICLE_EPOCH"] = np.array([epoch], dtype="<i4")
    op = st["RESULTS/ON_PARTICLES"]
    off = op["OFFSETS"][...]
    keep = slice(int(off[0]), int(off[1]))           # output 0 of stage 1: all 12 ids
    n = keep.stop - keep.start
    for name in list(op):
        if isinstance(op[name], h5py.Group):
            op[name]["DATA"].resize((n,) + op[name]["DATA"].shape[1:])
    op["ID"].resize((n,))
    op["ID"][...] = ids_back
    for k, v in (("STEP", [30]), ("TIME", [3e-3]), ("OFFSETS", [0, n])):
        op[k].resize((len(v),))
        op[k][...] = v


def test_mpm_reset_starts_a_new_epoch(tmp_path):
    """`mpm reset` brings removed ids back: legal in a new epoch (232a §2.7)."""
    p = tmp_path / "c2.ladruno"
    fx.write_contract2(p)
    with h5py.File(p, "r+") as f:
        _second_stage(f, 1, np.arange(12))
    assert mc.validate_particle_results(str(p)) == []


def test_contract2_shared_file_case_in_make_synthetic(tmp_path):
    lf_dir = Path(lf.__file__).parent
    sys.path.insert(0, str(lf_dir))
    import make_synthetic
    p = str(tmp_path / "shared.ladruno")
    make_synthetic.build(p, particles=True)
    assert lf.validate(p) == []


# --------------------------------------------------------------------------- Contract 1 refusals


def _c1(tmp_path, mutate=None, h5_mutate=None):
    m = fx.deep_copy_model()
    if mutate:
        mutate(m)
    p = tmp_path / "bad.h5"
    fx.write_contract1(p, m)
    if h5_mutate:
        with h5py.File(p, "r+") as f:
            h5_mutate(f)
    return _codes(mc.validate_particle_set(str(p)))


def _set(path, value):
    def f(m):
        d = m
        for k in path[:-1]:
            d = d[k]
        d[path[-1]] = value
    return f


def _pop(path):
    def f(m):
        d = m
        for k in path[:-1]:
            d = d[k]
        d.pop(path[-1])
    return f


def _particle(name, fn):
    def f(m):
        m["particles"][name] = fn(m["particles"][name].copy())
    return f


def _first(value, col=None):
    def f(a):
        if col is None:
            a[0] = value
        else:
            a[0, col] = value
        return a
    return f


def _repeat_id(a):
    a[1] = a[0]
    return a


def _add_flags(m):
    m["particles"]["FLAGS"] = np.r_[1, np.zeros(len(m["particles"]["ID"]) - 1)].astype("<u4")


def _both_mass(m):
    p = m["particles"]
    p["MASS"] = p["DENSITY"] * p["VOLUME"]


def _scaleoffset(f):
    """A filter outside shuffle/deflate/fletcher32 (lzf would do too, but h5py's lzf
    crashes this Windows python3.12 with a stack overflow)."""
    v = f["PARTICLES/VOLUME"][...]
    del f["PARTICLES/VOLUME"]
    f["PARTICLES"].create_dataset("VOLUME", data=v, chunks=(len(v),), scaleoffset=6)


def _stress33_in_3d(m):
    _three_d(m)
    m["particles"]["STRESS33"] = np.zeros(len(m["particles"]["ID"]))


def _ndm_int64(f):
    f.attrs["NDM"] = np.int64(2)


C1_CASES = {
    "missing VOLUME": (dict(mutate=_pop(["particles", "VOLUME"])), "C1-REQ"),
    "missing GRID": (dict(h5_mutate=lambda f: f.__delitem__("GRID")), "C1-REQ"),
    "wrong FORMAT": (dict(mutate=_set(["attrs", "FORMAT"], "ladruno-particles")), "C1-FMT"),
    "format version 2": (dict(mutate=_set(["attrs", "FORMAT_VERSION"], np.int32(2))), "C1-VER"),
    "empty GENERATOR": (dict(h5_mutate=lambda f: f.attrs.__setitem__("GENERATOR", "")),
                        "C1-GEN"),
    "bad set id": (dict(mutate=_set(["attrs", "SET_ID"], "not-a-uuid")), "C1-SETID"),
    "NDM 1": (dict(mutate=_set(["attrs", "NDM"], np.int32(1))), "C1-NDM"),
    "NDM stored int64": (dict(h5_mutate=_ndm_int64), "C1-DTYPE"),
    "empty UNIT_FORCE": (dict(h5_mutate=lambda f: f.attrs.__setitem__("UNIT_FORCE", "")),
                         "C1-UNIT"),
    "NUM_PARTICLES wrong": (dict(mutate=_set(["attrs", "NUM_PARTICLES"], np.int64(13))),
                            "C1-N"),
    "2D without thickness": (dict(mutate=_pop(["attrs", "PLANE_STRAIN_THICKNESS"])), "C1-THK"),
    "mass and density": (dict(mutate=_both_mass), "C1-MASS"),
    "compression-positive sign": (dict(mutate=_set(["attrs", "STRESS_SIGN"],
                                                   "compression-positive")), "C1-SIGN"),
    "old 2D Voigt order": (dict(mutate=_set(["attrs", "VOIGT_ORDER"], "11,22,33,12")),
                           "C1-VOIGT"),
    "zero cell size": (dict(mutate=_set(["grid", "CELL_SIZE"], np.array([0.5, 0.0]))),
                       "C1-GRID"),
    "unknown basis": (dict(mutate=_set(["grid", "BASIS"], "bspline")), "C1-BASIS"),
    "duplicate body name": (dict(mutate=_set(["bodies"], [(0, "soil"), (1, "soil")])),
                            "C1-BODY"),
    "unknown body": (dict(mutate=_particle("BODY_ID", lambda a: a + 5)), "C1-BODYREF"),
    "repeated id": (dict(mutate=_particle("ID", _repeat_id)), "C1-ID"),
    "position float32": (dict(mutate=_particle("POSITION", lambda a: a.astype("<f4"))),
                         "C1-DTYPE"),
    "stress 4 columns in 2D": (dict(mutate=_particle("STRESS", lambda a: np.c_[a, a[:, :1]])),
                               "C1-SHAPE"),
    "STRESS33 in 3D": (dict(mutate=_stress33_in_3d), "C1-STRESS33"),
    "stress NaN": (dict(mutate=_particle("STRESS", _first(np.nan, 0))), "C1-FINITE"),
    "negative volume": (dict(mutate=_particle("VOLUME", _first(-0.0625))), "C1-VOL"),
    "zero volume": (dict(mutate=_particle("VOLUME", _first(0.0))), "C1-VOL"),
    "half-length > h/2": (dict(mutate=_particle("GIMP_HALF_LENGTH", _first(0.3, 0))),
                          "C1-GIMP"),
    "outside the grid box": (dict(mutate=_particle("POSITION", _first(0.05, 0))), "C1-POS"),
    "material tag 0": (dict(mutate=_particle("MATERIAL_TAG", lambda a: a * 0)), "C1-MAT"),
    "flag bit set": (dict(mutate=_add_flags), "C1-FLAGS"),
    "unknown field": (dict(mutate=_set(["particles", "PORE_PRESSURE"], np.zeros(12))),
                      "C1-UNKNOWN"),
    "apeGmsh without session": (dict(mutate=_pop(["provenance", "session_id"])), "C1-PROV"),
    "scale-offset filter": (dict(h5_mutate=_scaleoffset), "C1-FILTER"),
}


@pytest.mark.parametrize("case", sorted(C1_CASES))
def test_contract1_refuses(tmp_path, case):
    kwargs, code = C1_CASES[case]
    codes = _c1(tmp_path, **kwargs)
    assert code in codes, (case, codes)


# --------------------------------------------------------------------------- Contract 2 refusals


def _c2(tmp_path, *, model=None, outputs=None, h5_mutate=None, **kw):
    m = model if model is not None else fx.deep_copy_model()
    outs = outputs if outputs is not None else fx.contract2_outputs(m)
    p = tmp_path / "bad.ladruno"
    fx.write_contract2(p, m, outs, **kw)
    if h5_mutate:
        with h5py.File(p, "r+") as f:
            h5_mutate(f)
    return mc.validate_particle_results(str(p))


def _outs(fn):
    outs = fx.contract2_outputs()
    fn(outs)
    return outs


def _reuse_removed(outs):
    """Particle 11, removed after output 1, comes back in a fourth output."""
    o = outs[1]
    outs.append({"step": 30, "time": 3e-3, "ids": o["ids"].copy(),
                 "fields": {k: v.copy() for k, v in o["fields"].items()}})


def _repeat_in_output(outs):
    outs[1]["ids"] = outs[1]["ids"].copy()
    outs[1]["ids"][1] = outs[1]["ids"][0]


def _foreign_id(outs):
    outs[0]["ids"] = outs[0]["ids"].copy()
    outs[0]["ids"][-1] = 99


def _drop_field(name):
    def f(outs):
        for o in outs:
            o["fields"].pop(name)
    return f


def _field(name, fn):
    def f(outs):
        for o in outs:
            o["fields"][name] = fn(o["fields"][name].copy())
    return f


def _h5(fn):
    return dict(h5_mutate=fn)


def _setattr(path, name, value):
    def f(h):
        h[path].attrs[name] = value
    return f


def _delattr(path, name):
    def f(h):
        del h[path].attrs[name]
    return f


def _replace(path, data, **kw):
    def f(h):
        parent, leaf = path.rsplit("/", 1)
        del h[path]
        h[parent].create_dataset(leaf, data=data, **kw)
    return f


def _bad_offsets(h):
    ds = h[f"{OP}/OFFSETS"]
    v = ds[...]
    v[1] = v[2] + 1
    ds[...] = v


def _neg_mass(h):
    ds = h[f"{MP}/MASS"]
    v = ds[...]
    v[0] = -1.0
    ds[...] = v


def _partitioned(h):
    h["INFO"].attrs["PARTITIONED"] = np.array([1], dtype="<i4")
    h["INFO"].attrs["NUM_PARTITIONS"] = np.array([2], dtype="<i4")


def _short_data(h):
    d = h[f"{OP}/VOLUME/DATA"]
    d.resize((d.shape[0] - 1, 1))


def _extra_dataset(h):
    h[OP].create_dataset("NOTES", data=np.zeros(3))


def _volume_scaleoffset(h):
    v = h[f"{OP}/VOLUME/DATA"][...]
    del h[f"{OP}/VOLUME/DATA"]
    h[f"{OP}/VOLUME"].create_dataset("DATA", data=v, chunks=v.shape, scaleoffset=6)


def _nodes_not_empty_shape(h):
    del h[f"{STAGE}/MODEL/NODES/ID"]
    h[f"{STAGE}/MODEL/NODES"].create_dataset("ID", data=np.zeros(0, dtype="<i4"))


def _epoch_back(h):
    _second_stage(h, 0, np.arange(12))
    h[STAGE].attrs["PARTICLE_EPOCH"] = np.array([1], dtype="<i4")


def _same_epoch_reuse(h):
    _second_stage(h, 0, np.arange(12))


C2_CASES = {
    "minor 1": (dict(minor=1), "C2-MINOR"),
    "partitioned": (_h5(_partitioned), "C2-PART"),
    "SPATIAL_DIM 1": (_h5(_setattr("INFO", "SPATIAL_DIM", np.array([1], dtype="<i4"))),
                      "C2-NDM"),
    "model without results": (_h5(lambda h: h.__delitem__(OP)), "C2-PAIR"),
    "eigen stage": (_h5(_setattr(STAGE, "KIND", np.array([b"eigen"]))), "C2-KIND"),
    "KIND absent": (_h5(_delattr(STAGE, "KIND")), "C2-KIND"),
    "KIND empty": (_h5(_setattr(STAGE, "KIND", np.array([b""]))), "C2-KIND"),
    "unknown mode": (dict(mode="dynamic"), "C2-MODE"),
    "explicit under a static stage": (_h5(_setattr(STAGE, "KIND", np.array([b"static"]))),
                                      "C2-MODEKIND"),
    "no PATTERN_TAG": (_h5(_delattr(MP, "PATTERN_TAG")), "C2-PATTERN"),
    "no PARTICLE_EPOCH": (_h5(_delattr(STAGE, "PARTICLE_EPOCH")), "C2-EPOCH"),
    "PARTICLE_EPOCH goes back": (_h5(_epoch_back), "C2-EPOCH"),
    "no EMPTY_PARTITION": (_h5(_delattr(STAGE, "EMPTY_PARTITION")), "C2-EMPTY"),
    "NODES/ID not [0 x 1]": (_h5(_nodes_not_empty_shape), "C2-EMPTY"),
    "bad SOURCE_SET_ID": (_h5(_setattr(MP, "SOURCE_SET_ID", np.array([b"nope"]))), "C2-SRC"),
    "no thickness in 2D": (_h5(_delattr(MP, "PLANE_STRAIN_THICKNESS")), "C2-THK"),
    "no UNIT_TIME": (_h5(_delattr(MP, "UNIT_TIME")), "C2-UNIT"),
    "no grid": (_h5(lambda h: h.__delitem__(f"{MP}/GRID")), "C2-GRID"),
    "duplicate body name": (_h5(_replace(f"{MP}/BODIES/NAME", np.array([b"soil", b"soil"]))),
                            "C2-BODY"),
    "negative registry mass": (_h5(_neg_mass), "C2-REG"),
    "PARTITION_REDUCTION SUM": (_h5(_setattr(OP, "PARTITION_REDUCTION", np.array([b"SUM"]))),
                                "C2-RED"),
    "STEP not increasing": (_h5(_replace(f"{OP}/STEP", np.array([0, 20, 10], dtype="<i4"))),
                            "C2-AXES"),
    "bad OFFSETS": (_h5(_bad_offsets), "C2-OFF"),
    "id repeated in an output": (dict(outputs=_outs(_repeat_in_output)), "C2-ID"),
    "id not in the registry": (dict(outputs=_outs(_foreign_id)), "C2-ID"),
    "removed id comes back": (dict(outputs=_outs(_reuse_removed)), "C2-REUSE"),
    "removed id comes back in the same epoch": (_h5(_same_epoch_reuse), "C2-REUSE"),
    "no POSITION": (dict(outputs=_outs(_drop_field("POSITION"))), "C2-POS"),
    "POSITION outside the grid box": (dict(outputs=_outs(_field("POSITION",
                                                                _first(5.0, 0)))), "C2-POS"),
    "POSITION float32": (_h5(_replace(f"{OP}/POSITION/DATA",
                                      np.zeros((35, 2), dtype="<f4"))), "C2-DTYPE"),
    "VOLUME DATA one row short": (_h5(_short_data), "C2-DATA"),
    "no DISPLAY_NAME": (_h5(_delattr(f"{OP}/VOLUME", "DISPLAY_NAME")), "C2-ATTR"),
    "repeated component label": (_h5(_setattr(f"{OP}/VELOCITY", "COMPONENTS",
                                              np.array([b"Vx,Vx"]))), "C2-COMP"),
    "STRESS labelled as the old 2D order": (_h5(_setattr(f"{OP}/STRESS", "COMPONENTS",
                                                         np.array([b"sigma11,sigma22,sigma33"]))),
                                            "C2-CAT"),
    "NaN stress": (dict(outputs=_outs(_field("STRESS", _first(np.nan, 0)))), "C2-FINITE"),
    "lowercase field name": (_h5(lambda h: h[OP].move("VOLUME", "volume")), "C2-NAME"),
    "stray dataset": (_h5(_extra_dataset), "C2-UNKNOWN"),
    "scale-offset filter": (_h5(_volume_scaleoffset), "C2-FILTER"),
}


@pytest.mark.parametrize("case", sorted(C2_CASES))
def test_contract2_refuses(tmp_path, case):
    kwargs, code = C2_CASES[case]
    problems = _c2(tmp_path, **kwargs)
    assert code in _codes(problems), (case, problems)


def test_contract2_major_version_bump_is_refused(tmp_path):
    """A particle file stamped FORMAT_VERSION 2 is refused: the change is a minor (232a §2.1)."""
    problems = _c2(tmp_path, h5_mutate=_setattr("INFO", "FORMAT_VERSION",
                                                np.array([2], dtype="<i4")))
    assert any("FORMAT_VERSION must be 1" in p for p in problems), problems


def test_every_rule_code_has_a_refusal_case():
    """Each C1-/C2- code the contract page names is fired by at least one case above."""
    named = set(re.findall(r"\bC[12]-[A-Z0-9]+\b", DOC.read_text(encoding="utf-8")))
    covered = {code for _, code in C1_CASES.values()} | {code for _, code in C2_CASES.values()}
    assert named, "no rule codes found in the contract page"
    assert named <= covered, sorted(named - covered)
    assert covered <= named, sorted(covered - named)


# --------------------------------------------------------------------------- the crash rule (§2.5)


def test_torn_fixture_reads_its_committed_outputs():
    """A kill after OFFSETS was extended but before its value was written: the last entry
    holds the fill value 0. The in-flight output is dropped; the 3 before it read."""
    with h5py.File(fx.PARTICLE_TORN, "r") as f:
        op = f[OP]
        assert op["OFFSETS"][...].tolist() == [0, 12, 24, 35, 0]
        assert op["STEP"].shape == (4,) and op["TIME"].shape == (4,)
        assert op["ID"].shape == (46,)
    assert mc.committed_outputs(fx.PARTICLE_TORN) == 3
    with lf.LadrunoReader(str(fx.PARTICLE_TORN)) as r:
        ids, _ = r.particle_output(STAGE, "POSITION", 2)
    assert ids.tolist() == list(range(11))


def test_validator_still_refuses_a_torn_file():
    """The validator checks closed files: a torn one is C2-OFF."""
    assert "C2-OFF" in _codes(mc.validate_particle_results(str(fx.PARTICLE_TORN)))


def test_torn_fixture_matches_the_generator(tmp_path):
    fx.write_torn(tmp_path / "torn.ladruno")
    assert _same(_content(fx.PARTICLE_TORN), _content(tmp_path / "torn.ladruno"))


@pytest.mark.parametrize("path,stage,n", [
    (fx.PARTICLE_RESULTS, None, 3),
    (fx.PARTICLE_MIXED, "MODEL_STAGE[1]", 3),
])
def test_closed_fixtures_are_unchanged_by_the_rule(path, stage, n):
    assert mc.committed_outputs(path, stage) == n
    with h5py.File(path, "r") as f:
        st = f[stage or STAGE]
        assert mc.committed_outputs(st) == st["RESULTS/ON_PARTICLES/STEP"].shape[0] == n
        assert mc.committed_outputs(st["RESULTS/ON_PARTICLES"]) == n


def test_a_decrease_before_the_last_entry_is_refused(tmp_path):
    p = tmp_path / "corrupt.ladruno"
    fx.write_contract2(p)
    with h5py.File(p, "r+") as f:
        f[f"{OP}/OFFSETS"][1] = 30                 # [0, 30, 24, 35]: 30 > 24 inside
    with pytest.raises(mc.TornStageError):
        mc.committed_outputs(p)


def test_all_empty_outputs_read_a_torn_entry_as_one_more_empty_output(tmp_path):
    """When every offset is 0, the torn 0 equals its predecessor: one more empty output,
    harmless because its STEP and TIME were written before OFFSETS (232a §2.5)."""
    p = tmp_path / "empty.ladruno"
    fx.write_contract2(p)
    with h5py.File(p, "r+") as f:
        f[f"{OP}/OFFSETS"][...] = 0
    fx.tear_last_output(p)
    with h5py.File(p, "r") as f:
        assert f[f"{OP}/OFFSETS"][...].tolist() == [0, 0, 0, 0, 0]
    assert mc.committed_outputs(p) == 4


def test_an_output_whose_rows_are_missing_is_not_committed(tmp_path):
    """The len(ID) / len(DATA) bound: an output whose rows did not all reach disk."""
    p = tmp_path / "short.ladruno"
    fx.write_contract2(p)
    with h5py.File(p, "r+") as f:
        d = f[f"{OP}/VELOCITY/DATA"]
        d.resize((d.shape[0] - 1,) + d.shape[1:])
    assert mc.committed_outputs(p) == 2


def test_a_stage_torn_while_opening_holds_no_committed_output(tmp_path):
    """Step 0: the registry is written once, before output 0. Registry columns of
    different lengths mean a kill while the stage opened: nothing is committed."""
    p = tmp_path / "opening.ladruno"
    fx.write_contract2(p)
    assert mc.committed_outputs(p) == 3
    with h5py.File(p, "r+") as f:
        mp = f[f"{STAGE}/MODEL/PARTICLES"]
        mass = mp["MASS"][...]
        del mp["MASS"]
        mp.create_dataset("MASS", data=mass[:-1])
    assert mc.committed_outputs(p) == 0


_SWMR_CHILD = r"""
import sys, h5py
f = h5py.File(sys.argv[1], "w", libver="latest")
d = f.create_dataset("D", shape=(0,), maxshape=(None,), chunks=(16,), dtype="i8")
f.flush()
f.swmr_mode = True
print("GO", flush=True)
while True:
    d.resize((d.shape[0] + 16,)); d[-16:] = 1; f.flush()
"""


def test_the_reader_rule_opens_a_killed_swmr_file(tmp_path):
    """232a §2.5 reader rule, against an h5py SWMR writer killed mid-run: a plain
    read-only open is refused, `open_particle_file` falls back to SWMR-read, with no
    HDF5_USE_FILE_LOCKING override; a cleanly closed file opens plainly."""
    import os
    import subprocess
    import sys
    import time
    p = tmp_path / "killed.h5"
    env = dict(os.environ)
    env.pop("HDF5_USE_FILE_LOCKING", None)
    w = subprocess.Popen([sys.executable, "-c", _SWMR_CHILD, str(p)], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, stdin=subprocess.DEVNULL, text=True, env=env)
    try:
        assert "GO" in w.stdout.readline()
        time.sleep(0.5)
    finally:
        w.kill()
        w.wait()
    with pytest.raises(OSError, match="already open for write"):
        h5py.File(p, "r")
    with mc.open_particle_file(p) as f:
        assert f.swmr_mode and f["D"].shape[0] > 0
        assert (f["D"][...] == 1).all()
    q = tmp_path / "closed.h5"
    with h5py.File(q, "w", libver="latest") as f:
        f.create_dataset("D", shape=(0,), maxshape=(None,), chunks=(16,), dtype="i8")
        f.swmr_mode = True
    with mc.open_particle_file(q) as f:
        assert not f.swmr_mode


def _layout_killed(path, extra_eoa=0, last_stage_partial=False, first_stage_broken=False):
    """A particle-like file in the latest format, two complete stages (copies of the
    fixture's), its superblock flags then set to WRITE_ACCESS only (0x01): the state a
    kill leaves between a reset's read-write reopen and H5Fstart_swmr_write."""
    with h5py.File(fx.PARTICLE_RESULTS, "r") as s, h5py.File(path, "w", libver="latest") as f:
        s.copy(s["INFO"], f, "INFO")
        for k in (0, 1):
            s.copy(s[STAGE], f, f"MODEL_STAGE[{k}]")
        if first_stage_broken:
            del f["MODEL_STAGE[0]/RESULTS"]
        if last_stage_partial:
            f.create_group("MODEL_STAGE[2]").create_group("MODEL")
    flags, _ = mc.set_superblock(path, flags=0x01)
    assert flags == 0x00
    if extra_eoa:
        mc.set_superblock(path, eoa=path.stat().st_size + extra_eoa)


def test_a_file_killed_mid_layout_is_refused_by_both_opens_and_recovers(tmp_path):
    """232a §2.5 reader rule step 3: flags 0x01 (WRITE_ACCESS only) are refused by the
    plain AND the SWMR-read open; open_particle_file names the recovery; the recovery
    works on a copy and returns every stage."""
    p, q = tmp_path / "killed.ladruno", tmp_path / "recovered.ladruno"
    _layout_killed(p)
    before = p.read_bytes()
    with pytest.raises(OSError, match="already open for write"):
        h5py.File(p, "r")
    with pytest.raises(OSError, match="not already open for SWMR writing"):
        h5py.File(p, "r", swmr=True)
    with pytest.raises(mc.StageLayoutKilledError, match="recover_particle_file"):
        mc.open_particle_file(p)
    assert mc.recover_particle_file(p, q) == [
        {"stage": "MODEL_STAGE[0]", "committed": 3},
        {"stage": "MODEL_STAGE[1]", "committed": 3}]
    assert p.read_bytes() == before                      # the source is never repaired
    with h5py.File(q, "r") as f:                          # the copy opens plainly
        assert not f.swmr_mode
    with pytest.raises(ValueError, match="copy"):
        mc.recover_particle_file(p, p)


def test_recovery_pads_a_file_shorter_than_its_stored_eoa(tmp_path):
    """eof < stored EOA (seen in 1 of the 21 measured refusals): the copy is padded with
    zeros up to the stored EOA, then opens plainly."""
    p, q = tmp_path / "short.ladruno", tmp_path / "recovered.ladruno"
    _layout_killed(p, extra_eoa=4096)
    size = p.stat().st_size
    with pytest.raises(OSError, match="truncated file"):  # checked before the flags
        h5py.File(p, "r")
    with pytest.raises(mc.StageLayoutKilledError):
        mc.open_particle_file(p)
    rep = mc.recover_particle_file(p, q)
    assert q.stat().st_size == size + 4096
    assert [r["committed"] for r in rep] == [3, 3]


def test_recovery_reports_the_stage_being_laid_out_and_refuses_an_earlier_broken_one(tmp_path):
    """Only the LAST stage may be torn while laid out; an earlier stage that will not
    read is corruption."""
    p, q = tmp_path / "partial.ladruno", tmp_path / "recovered.ladruno"
    _layout_killed(p, last_stage_partial=True)
    rep = mc.recover_particle_file(p, q)
    assert [r.get("committed") for r in rep[:2]] == [3, 3]
    assert rep[2]["stage"] == "MODEL_STAGE[2]" and rep[2]["torn_while_laid_out"]
    p2, q2 = tmp_path / "broken.ladruno", tmp_path / "recovered2.ladruno"
    _layout_killed(p2, first_stage_broken=True)
    with pytest.raises(mc.TornStageError, match="not the last stage"):
        mc.recover_particle_file(p2, q2)


def test_the_repair_reads_a_file_killed_in_swmr_mode_like_the_swmr_read_open(tmp_path):
    """232a §2.5: a reader without SWMR-read (e.g. h5wasm) may apply step 3's repair to
    ANY file step 1 refuses, flags 0x05 included. The image of a live SWMR writer (what a
    kill leaves) reads the same through the repair as through the SWMR-read open."""
    live, image, rec = tmp_path / "live.ladruno", tmp_path / "image.ladruno", tmp_path / "rec.ladruno"
    m = 5
    f = h5py.File(live, "w", libver="latest")
    try:
        with h5py.File(fx.PARTICLE_RESULTS, "r") as s:
            s.copy(s["INFO"], f, "INFO")
        st = f.create_group("MODEL_STAGE[0]")
        mp = st.create_group("MODEL/PARTICLES")
        for name, dt, val in (("ID", "<i8", np.arange(m)), ("BODY_ID", "<i4", np.ones(m)),
                              ("MATERIAL_TAG", "<i4", np.ones(m)), ("MASS", "<f8", np.full(m, 2.0))):
            mp.create_dataset(name, data=np.asarray(val, dtype=dt))
        op = st.create_group("RESULTS/ON_PARTICLES")
        ds = {k: op.create_dataset(k, shape=(0,), maxshape=(None,), chunks=(64,), dtype=dt)
              for k, dt in (("STEP", "<i4"), ("TIME", "<f8"), ("OFFSETS", "<i8"), ("ID", "<i8"))}
        pos = op.create_group("POSITION").create_dataset(
            "DATA", shape=(0, 2), maxshape=(None, 2), chunks=(m, 2), dtype="<f8")
        ds["OFFSETS"].resize((1,))
        ds["OFFSETS"][0] = 0
        f.flush()
        f.swmr_mode = True
        for k in range(4):                                # rows, barrier, OFFSETS, flush
            pos.resize((pos.shape[0] + m, 2))
            pos[-m:] = k + np.arange(2 * m).reshape(m, 2) / 10.0
            for name, val in (("ID", np.arange(m)), ("STEP", [k + 1]), ("TIME", [(k + 1) * 1e-3])):
                d = ds[name]
                d.resize((d.shape[0] + len(val),))
                d[-len(val):] = val
            f.flush()
            ds["OFFSETS"].resize((k + 2,))
            ds["OFFSETS"][k + 1] = (k + 1) * m
            f.flush()
        image.write_bytes(live.read_bytes())              # the image a kill leaves
    finally:
        f.close()
    assert mc._superblock_flags(image) == 0x05
    with pytest.raises(OSError, match="already open for write|truncated file"):
        h5py.File(image, "r")

    def rows(h):
        n = mc.committed_outputs(h["MODEL_STAGE[0]"])
        op = h["MODEL_STAGE[0]/RESULTS/ON_PARTICLES"]
        end = int(op["OFFSETS"][n])
        return n, op["ID"][:end], op["POSITION/DATA"][:end]

    with mc.open_particle_file(image) as h:
        assert h.swmr_mode
        n1, ids1, pos1 = rows(h)
    rep = mc.recover_particle_file(image, rec)
    with h5py.File(rec, "r") as h:
        n2, ids2, pos2 = rows(h)
    assert n1 == n2 == 4 and rep == [{"stage": "MODEL_STAGE[0]", "committed": 4}]
    assert np.array_equal(ids1, ids2) and pos1.tobytes() == pos2.tobytes()
