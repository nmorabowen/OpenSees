"""WP-232a: the MPM contracts -- fixtures, validators, malformed variants.

Pure h5py + numpy, no OpenSees build (run by the static-gates job of ladruno.yml). The
contracts are Ladruno_implementation/232a_mpm_contracts.md; the validators are
tests/mpm_contracts.py (Contract 1) and Ladruno_scripts/ladruno_recorder_tests/
ladruno_format.py (Contract 2, via mpm_contracts.validate_particle_results). Each malformed
variant must be refused by the rule code the contract names for it.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mpm_contract_fixtures as fx  # noqa: E402
import mpm_contracts as mc  # noqa: E402

lf = mc.lf


# --------------------------------------------------------------------------- helpers


def _content(path) -> dict:
    """Every object's dtype / shape / values and every attribute, for a content diff."""
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
    assert _same(_content(fx.PARTICLE_SET), _content(tmp_path / "c1.h5"))
    assert _same(_content(fx.PARTICLE_RESULTS), _content(tmp_path / "c2.ladruno"))


def test_contract2_is_an_additive_minor_not_a_major():
    """The orchestrator's ruling: FORMAT_VERSION stays 1, the minor is 2 (232a §2.1)."""
    with lf.LadrunoReader(str(fx.PARTICLE_RESULTS)) as r:
        info = r.info()
    assert int(info["FORMAT_VERSION"]) == lf.FORMAT_VERSION == 1
    assert int(info["FORMAT_MINOR_VERSION"]) == lf.FORMAT_MINOR_VERSION_PARTICLES == 2


def test_contract2_fixture_uses_the_fork_attribute_encodings():
    """1-element int arrays, fixed-length strings, empty strings absent (232a §0)."""
    with h5py.File(fx.PARTICLE_RESULTS, "r") as f:
        fv = f["INFO"].attrs["FORMAT_VERSION"]
        assert fv.shape == (1,) and fv.dtype == np.dtype("<i4")
        kind = f["MODEL_STAGE[1]"].attrs["KIND"]
        assert kind.shape == (1,) and kind.dtype.kind == "S"
        fg = f["MODEL_STAGE[1]/RESULTS/ON_PARTICLES/DEFORMATION_GRADIENT"]
        assert "DIMENSION" not in fg.attrs and "DESCRIPTION" not in fg.attrs
        assert f["MODEL_STAGE[1]/MODEL/PARTICLES/BODIES/NAME"].dtype.kind == "S"


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
        assert np.array_equal(s0, P["STRESS"][...])
        _, v0 = r.particle_output(stage, "VOLUME", 0)
        assert np.array_equal(v0[:, 0], P["VOLUME"][...])


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
    fx.write_contract1(tmp_path / "c1.h5", m)
    assert mc.validate_particle_set(str(tmp_path / "c1.h5")) == []


def test_contract1_three_dimensional_set(tmp_path):
    m = fx.deep_copy_model()
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
    fx.write_contract1(tmp_path / "c1.h5", m)
    assert mc.validate_particle_set(str(tmp_path / "c1.h5")) == []


def test_contract2_accepts_an_added_particle(tmp_path):
    """Addition: a new id above every used one, appended to the registry (232a §2.7)."""
    m = fx.deep_copy_model()
    outs = fx.contract2_outputs(m)
    p = m["particles"]
    for k in ("ID", "BODY_ID", "MATERIAL_TAG", "VOLUME", "DENSITY"):
        p[k] = np.r_[p[k], p[k][-1:] + (1 if k == "ID" else 0)]
    for k in ("POSITION", "GIMP_HALF_LENGTH", "STRESS"):
        p[k] = np.vstack([p[k], p[k][-1:]])
    last = outs[-1]
    last["ids"] = np.r_[last["ids"], 12]
    last["fields"] = {k: np.vstack([v, v[-1:]]) for k, v in last["fields"].items()}
    fx.write_contract2(tmp_path / "c2.ladruno", m, outs)
    assert mc.validate_particle_results(str(tmp_path / "c2.ladruno")) == []


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


def _repeat_id(a):
    a[1] = a[0]
    return a


def _neg_first(a):
    a[0] = -a[0]
    return a


def _zero_first(a):
    a[0] = 0.0
    return a


def _nan_first(a):
    a[0, 0] = np.nan
    return a


def _outside(a):
    a[0, 0] = 0.05                       # x - l = -0.075 < grid origin
    return a


def _too_big(a):
    a[0, 0] = 0.3                        # > CELL_SIZE / 2 = 0.25
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


C1_CASES = {
    "missing VOLUME": (dict(mutate=_pop(["particles", "VOLUME"])), "C1-REQ"),
    "missing BODY_ID": (dict(mutate=_pop(["particles", "BODY_ID"])), "C1-REQ"),
    "missing GRID": (dict(h5_mutate=lambda f: f.__delitem__("GRID")), "C1-REQ"),
    "negative volume": (dict(mutate=_particle("VOLUME", _neg_first)), "C1-VOL"),
    "zero volume": (dict(mutate=_particle("VOLUME", _zero_first)), "C1-VOL"),
    "repeated id": (dict(mutate=_particle("ID", _repeat_id)), "C1-ID"),
    "position float32": (dict(mutate=_particle("POSITION", lambda a: a.astype("<f4"))),
                         "C1-DTYPE"),
    "id int32": (dict(mutate=_particle("ID", lambda a: a.astype("<i4"))), "C1-DTYPE"),
    "stress NaN": (dict(mutate=_particle("STRESS", _nan_first)), "C1-FINITE"),
    "stress 6 columns in 2D": (dict(mutate=_particle("STRESS", lambda a: np.c_[a, a[:, :2]])),
                               "C1-SHAPE"),
    "outside the grid": (dict(mutate=_particle("POSITION", _outside)), "C1-POS"),
    "half-length > h/2": (dict(mutate=_particle("GIMP_HALF_LENGTH", _too_big)), "C1-GIMP"),
    "unknown body": (dict(mutate=_particle("BODY_ID", lambda a: a + 5)), "C1-BODYREF"),
    "material tag 0": (dict(mutate=_particle("MATERIAL_TAG", lambda a: a * 0)), "C1-MAT"),
    "flag bit set": (dict(mutate=_add_flags), "C1-FLAGS"),
    "mass and density": (dict(mutate=_both_mass), "C1-MASS"),
    "unknown field": (dict(mutate=_set(["particles", "PORE_PRESSURE"], np.zeros(12))),
                      "C1-UNKNOWN"),
    "2D without thickness": (dict(mutate=_pop(["attrs", "PLANE_STRAIN_THICKNESS"])), "C1-THK"),
    "wrong Voigt order": (dict(mutate=_set(["attrs", "VOIGT_ORDER"], "11,22,12,33")),
                          "C1-VOIGT"),
    "compression-negative sign": (dict(mutate=_set(["attrs", "STRESS_SIGN"],
                                                   "compression-positive")), "C1-SIGN"),
    "format version 2": (dict(mutate=_set(["attrs", "FORMAT_VERSION"], np.int32(2))), "C1-VER"),
    "bad set id": (dict(mutate=_set(["attrs", "SET_ID"], "not-a-uuid")), "C1-SETID"),
    "NUM_PARTICLES wrong": (dict(mutate=_set(["attrs", "NUM_PARTICLES"], np.int64(13))),
                            "C1-N"),
    "unknown basis": (dict(mutate=_set(["grid", "BASIS"], "bspline")), "C1-BASIS"),
    "zero cell size": (dict(mutate=_set(["grid", "CELL_SIZE"], np.array([0.5, 0.0]))),
                       "C1-GRID"),
    "duplicate body name": (dict(mutate=_set(["bodies"], [(0, "soil"), (1, "soil")])),
                            "C1-BODY"),
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


def _position_f32(f):
    ds = f["MODEL_STAGE[1]/RESULTS/ON_PARTICLES/POSITION/DATA"]
    v = ds[...]
    del f["MODEL_STAGE[1]/RESULTS/ON_PARTICLES/POSITION/DATA"]
    f["MODEL_STAGE[1]/RESULTS/ON_PARTICLES/POSITION"].create_dataset("DATA", data=v.astype("<f4"))


def _bad_offsets(f):
    ds = f["MODEL_STAGE[1]/RESULTS/ON_PARTICLES/OFFSETS"]
    v = ds[...]
    v[1] = v[2] + 1
    ds[...] = v


def _stress_comps(f):
    g = f["MODEL_STAGE[1]/RESULTS/ON_PARTICLES/STRESS"]
    g.attrs["COMPONENTS"] = np.array([b"s11,s22,s33,s12"])


def _neg_mass(f):
    ds = f["MODEL_STAGE[1]/MODEL/PARTICLES/MASS"]
    v = ds[...]
    v[0] = -1.0
    ds[...] = v


def _drop_results(f):
    del f["MODEL_STAGE[1]/RESULTS/ON_PARTICLES"]


def _partitioned(f):
    f["INFO"].attrs["PARTITIONED"] = np.array([1], dtype="<i4")
    f["INFO"].attrs["NUM_PARTITIONS"] = np.array([2], dtype="<i4")


def _major_2(f):
    f["INFO"].attrs["FORMAT_VERSION"] = np.array([2], dtype="<i4")


def _eigen(f):
    f["MODEL_STAGE[1]"].attrs["KIND"] = np.array([b"eigen"])


def _no_grid(f):
    del f["MODEL_STAGE[1]/MODEL/PARTICLES/GRID"]


def _moved_far(a):
    a[0, 0] = 5.0
    return a


def _nan_stress(a):
    a[0, 0] = np.nan
    return a


C2_CASES = {
    "removed id comes back": (dict(outputs=_outs(_reuse_removed)), "C2-REUSE"),
    "id repeated in an output": (dict(outputs=_outs(_repeat_in_output)), "C2-ID"),
    "id not in the registry": (dict(outputs=_outs(_foreign_id)), "C2-ID"),
    "no POSITION": (dict(outputs=_outs(_drop_field("POSITION"))), "C2-POS"),
    "POSITION outside the grid": (dict(outputs=_outs(_field("POSITION", _moved_far))),
                                  "C2-POS"),
    "NaN stress": (dict(outputs=_outs(_field("STRESS", _nan_stress))), "C2-FINITE"),
    "POSITION float32": (dict(h5_mutate=_position_f32), "C2-DTYPE"),
    "bad OFFSETS": (dict(h5_mutate=_bad_offsets), "C2-OFF"),
    "STRESS labels": (dict(h5_mutate=_stress_comps), "C2-CAT"),
    "negative registry mass": (dict(h5_mutate=_neg_mass), "C2-REG"),
    "model without results": (dict(h5_mutate=_drop_results), "C2-PAIR"),
    "minor 1": (dict(minor=1), "C2-MINOR"),
    "partitioned": (dict(h5_mutate=_partitioned), "C2-PART"),
    "eigen stage": (dict(h5_mutate=_eigen), "C2-KIND"),
    "unknown mode": (dict(mode="dynamic"), "C2-MODE"),
    "no grid": (dict(h5_mutate=_no_grid), "C2-GRID"),
}


@pytest.mark.parametrize("case", sorted(C2_CASES))
def test_contract2_refuses(tmp_path, case):
    kwargs, code = C2_CASES[case]
    problems = _c2(tmp_path, **kwargs)
    assert code in _codes(problems), (case, problems)


def test_contract2_major_version_bump_is_refused(tmp_path):
    """A particle file stamped FORMAT_VERSION 2 is refused: the change is a minor (232a §2.1)."""
    problems = _c2(tmp_path, h5_mutate=_major_2)
    assert any("FORMAT_VERSION must be 1" in p for p in problems), problems
