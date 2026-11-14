"""WP-232a: write the two MPM contract fixtures (python3.12, h5py + numpy).

    python3.12 tests/mpm_contract_fixtures.py        # rewrites tests/data/mpm/*

* ``tests/data/mpm/particle_set_v1.h5`` -- Contract 1, written the way apeGmsh (h5py) will
  write it: scalar attributes, variable-length UTF-8 strings.
* ``tests/data/mpm/particles_mixed_v1m2.ladruno`` -- Contract 2 after an FE stage:
  make_synthetic's static FE stage, then the same particle stage (232a §4).
* ``tests/data/mpm/particles_v1m2.ladruno`` -- Contract 2: a minimal .ladruno shell
  (FORMAT_VERSION 1, FORMAT_MINOR_VERSION 2) holding only the particle group, written the
  way the fork's C++ wrapper writes attributes (1-element arrays, fixed-length strings,
  empty strings absent) through ``make_synthetic.write_particle_group``.

One model seen from both sides (232a §3-§4): 2D plane strain, t = 1 m, kN_m_s; a 4 x 4 grid
of 0.5 m cells (GIMP); body ``soil`` (8 particles, tag 2, K0 stress, rho 2.0) and body
``block`` (4 particles, tag 1, rho 7.85), 2 x 2 particles per cell. Contract 2: one
transient stage (LadrunoRecorder's empty stage plus the particle groups), outputs at steps
0, 10, 20; particle 11 is removed after output 1.

The builders take the model as data so ``tests/test_mpm_contracts.py`` can write each
malformed variant from a mutated copy. Everything is deterministic (no clock, no uuid4()).
"""
from __future__ import annotations

import copy
import sys
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mpm_contracts as mc  # noqa: E402  (puts the fork's format tools on sys.path)
import make_synthetic  # noqa: E402

lf = mc.lf

DATA = HERE / "data" / "mpm"
PARTICLE_SET = DATA / "particle_set_v1.h5"
PARTICLE_RESULTS = DATA / "particles_v1m2.ladruno"
PARTICLE_MIXED = DATA / "particles_mixed_v1m2.ladruno"

SET_ID = "3f6c1d2e-8a4b-4c5d-9e7f-0a1b2c3d4e5f"
SESSION_ID = "9b2e7c41-5d3a-4f6e-8c1b-2a3d4e5f6a7b"
G = 9.81
K0 = 0.5
SOIL_TOP = 1.0


def _lattice(cells_x, cells_y, h=0.5, n=2):
    """Particle centres of an n x n lattice in each listed cell (row-major, y then x)."""
    off = (np.arange(n) + 0.5) * h / n
    xs = np.concatenate([i * h + off for i in cells_x])
    ys = np.concatenate([j * h + off for j in cells_y])
    return np.array([[x, y] for y in ys for x in xs])


def contract1_model() -> dict:
    """The particle set as plain data (Contract 1, §1)."""
    soil = _lattice([1, 2], [1])           # x 0.5..1.5, y 0.5..1.0 -> 8 particles
    block = _lattice([2], [2])             # x 1.0..1.5, y 1.0..1.5 -> 4 particles
    x = np.vstack([soil, block])
    n = x.shape[0]
    rho = np.r_[np.full(len(soil), 2.0), np.full(len(block), 7.85)]
    stress = np.zeros((n, 3))                           # plane strain: [11, 22, 12]
    stress33 = np.zeros(n)                              # sigma_33, its own dataset
    syy = -2.0 * G * (SOIL_TOP - soil[:, 1])            # tension positive: compression < 0
    stress[: len(soil)] = np.c_[K0 * syy, syy, np.zeros(len(soil))]
    stress33[: len(soil)] = K0 * syy                    # completes the K0 state
    return {
        "attrs": {
            "FORMAT": mc.FORMAT, "FORMAT_VERSION": np.int32(1), "GENERATOR": "apeGmsh",
            "SET_ID": SET_ID, "CREATED_UTC": "2026-11-13T00:00:00Z", "NDM": np.int32(2),
            "NUM_PARTICLES": np.int64(n), "PLANE_STRAIN_THICKNESS": np.float64(1.0),
            "UNIT_LENGTH": "m", "UNIT_FORCE": "kN", "UNIT_TIME": "s", "UNIT_SYSTEM": "kN_m_s",
            "MASS_SOURCE": "density", "STRESS_SIGN": "tension-positive",
            "VOIGT_ORDER": mc.VOIGT[2],
        },
        "grid": {"ORIGIN": np.array([0.0, 0.0]), "CELL_SIZE": np.array([0.5, 0.5]),
                 "NUM_CELLS": np.array([4, 4], dtype=np.int32), "BASIS": "gimp"},
        "bodies": [(0, "soil"), (1, "block")],
        "provenance": {"apegmsh_version": "fixture-wp232a", "apegmsh_commit": "unknown",
                       "session_id": SESSION_ID},
        "seeding": {"METHOD": "lattice", "PARTICLES_PER_CELL": np.array([2, 2], dtype=np.int32)},
        "particles": {
            "ID": np.arange(n, dtype="<i8"),
            "POSITION": x.astype("<f8"),
            "VOLUME": np.full(n, 0.25 * 0.25 * 1.0),
            "DENSITY": rho,
            "GIMP_HALF_LENGTH": np.full((n, 2), 0.125),
            "STRESS": stress,
            "STRESS33": stress33,
            "BODY_ID": np.r_[np.zeros(len(soil)), np.ones(len(block))].astype("<i4"),
            "MATERIAL_TAG": np.r_[np.full(len(soil), 2), np.full(len(block), 1)].astype("<i4"),
        },
    }


def write_contract1(path, model: dict | None = None, compress: int = 1) -> None:
    """Contract 1 as apeGmsh writes it (h5py: scalar attrs, vlen UTF-8 strings)."""
    m = contract1_model() if model is None else model
    with h5py.File(path, "w", track_order=True) as f:
        for k, v in m["attrs"].items():
            f.attrs[k] = v
        g = f.create_group("GRID")
        for k, v in m["grid"].items():
            g.attrs[k] = v
        b = f.create_group("BODIES")
        b.create_dataset("ID", data=np.array([i for i, _ in m["bodies"]], dtype="<i4"))
        b.create_dataset("NAME", data=[s for _, s in m["bodies"]], dtype=h5py.string_dtype())
        prov = f.create_group("PROVENANCE")
        for k, v in m["provenance"].items():
            prov.attrs[k] = v
        if m.get("seeding"):
            s = f.create_group("SEEDING")
            for k, v in m["seeding"].items():
                s.attrs[k] = v
        P = f.create_group("PARTICLES")
        for k, v in m["particles"].items():
            v = np.asarray(v)
            kw = {"compression": "gzip", "compression_opts": compress, "shuffle": True} \
                if compress else {}
            P.create_dataset(k, data=v, chunks=(max(1, v.shape[0]),) + v.shape[1:], **kw)


def contract2_outputs(model: dict | None = None) -> list[dict]:
    """Three outputs of the run the fork would make from ``model`` (232a §4)."""
    m = contract1_model() if model is None else model
    p = m["particles"]
    n = len(p["ID"])
    block = p["BODY_ID"] == 1
    soil = ~block
    outs = []
    for step, keep_all in ((0, True), (10, True), (20, False)):
        t = step * 1.0e-4
        x = p["POSITION"].copy()
        v = np.zeros((n, 2))
        x[block, 1] -= 0.5 * G * t * t                  # the block falls freely
        v[block, 1] = -G * t
        eps = np.full((n, 1), np.nan)
        eps[soil, 0] = 1.0e-5 * step
        fields = {
            "POSITION": x, "VELOCITY": v, "STRESS": p["STRESS"].copy(),
            "VOLUME": p["VOLUME"].reshape(-1, 1).copy(),
            "STRAIN": np.zeros((n, 3)),                  # eps_hat: unrotated, engineering shear
            "DEFORMATION_GRADIENT": np.tile([1.0, 0.0, 0.0, 1.0], (n, 1)),
            "STATE_eqPlasticStrain": eps,
        }
        rows = np.arange(n) if keep_all else np.arange(n - 1)   # particle 11 removed
        outs.append({"step": step, "time": t, "ids": p["ID"][rows],
                     "fields": {k: a[rows] for k, a in fields.items()}})
    return outs


def _mass(m: dict) -> np.ndarray:
    p = m["particles"]
    return p["MASS"] if "MASS" in p else p["DENSITY"] * p["VOLUME"]


def write_contract2(path, model: dict | None = None, outputs: list | None = None,
                    minor: int = 2, mode: str = "explicit", compress: int = 1,
                    epoch: int = 0, source: bool = True) -> None:
    """A minimal .ladruno shell with the particle group, as the fork writes it: the stage
    is LadrunoRecorder's empty stage (EMPTY_PARTITION = 1, NODES/ID int32 [0 x 1], empty
    ELEMENTS and ON_NODES / ON_ELEMENTS / ON_DOMAIN / ON_REGIONS) plus the particle groups."""
    m = contract1_model() if model is None else model
    outs = contract2_outputs(m) if outputs is None else outputs
    a = m["attrs"]
    ai, as_ = make_synthetic._ai, make_synthetic._as
    with h5py.File(path, "w", track_order=True) as f:
        info = f.create_group("INFO")
        as_(info, "GENERATOR", "Ladruno")
        ai(info, "FORMAT_VERSION", 1)
        ai(info, "FORMAT_MINOR_VERSION", minor)
        as_(info, "SOLVER_NAME", "OpenSees")
        ai(info, "SOLVER_VERSION", [3, 5, 1])
        ai(info, "SPATIAL_DIM", int(a["NDM"]))
        as_(info, "STORED_PRECISION", "f64")
        ai(info, "PARTITION_ID", 0)
        ai(info, "PARTITIONED", 0)
        ai(info, "NUM_PARTITIONS", 1)
        as_(info, "RUN_ID", "wp232a-fixture")
        as_(info, "RUN_ID_SCOPE", "user")
        info.create_group("PROVENANCE")

        write_particle_stage(f, "MODEL_STAGE[1]", m, outs, mode=mode, compress=compress,
                             epoch=epoch, source=source)


def write_particle_stage(f, name: str, m: dict, outs: list, mode: str = "explicit",
                         compress: int = 1, epoch: int = 0, source: bool = True,
                         step: int = 0) -> None:
    """One particle-only stage: LadrunoRecorder's empty stage plus the particle groups.
    KIND is pinned from MODE (232a §2.3)."""
    a, p = m["attrs"], m["particles"]
    ai, ad, as_ = make_synthetic._ai, make_synthetic._ad, make_synthetic._as
    st = f.create_group(name)
    ai(st, "STEP", step)
    ad(st, "TIME", 0.0)
    as_(st, "KIND", lf.MODE_KIND.get(mode, "transient"))
    ai(st, "EMPTY_PARTITION", 1)
    model = st.create_group("MODEL")
    nodes = model.create_group("NODES")
    ndm = int(a["NDM"])
    nodes.create_dataset("ID", data=np.zeros((0, 1), dtype="<i4"))
    nodes.create_dataset("COORDINATES", data=np.zeros((0, ndm), dtype="<f8"))
    model.create_group("ELEMENTS")
    res = st.create_group("RESULTS")
    for g in lf.EMPTY_RESULT_GROUPS:
        res.create_group(g)

    meta = {"PATTERN_TAG": 1, "MODE": mode,
            "UNIT_LENGTH": a["UNIT_LENGTH"], "UNIT_FORCE": a["UNIT_FORCE"],
            "UNIT_TIME": a["UNIT_TIME"]}
    if source:
        meta["SOURCE_SET_ID"] = a["SET_ID"]
    if "UNIT_SYSTEM" in a:
        meta["UNIT_SYSTEM"] = a["UNIT_SYSTEM"]
    if ndm == 2:
        meta["PLANE_STRAIN_THICKNESS"] = float(a["PLANE_STRAIN_THICKNESS"])
    make_synthetic.write_particle_group(
        st, ndm=ndm, grid=m["grid"], bodies=m["bodies"],
        registry={"ID": p["ID"], "BODY_ID": p["BODY_ID"],
                  "MATERIAL_TAG": p["MATERIAL_TAG"], "MASS": _mass(m)},
        outputs=outs, meta=meta,
        field_meta={"STATE_eqPlasticStrain": (
            "Equivalent plastic strain (illustrative token)", "epsPeq", "", 0)},
        compress=compress, epoch=epoch)


def deep_copy_model() -> dict:
    return copy.deepcopy(contract1_model())


def write_mixed(path, model: dict | None = None) -> None:
    """An FE stage, then a particle stage (232a §4): make_synthetic's static FE file as
    MODEL_STAGE[0], minor raised to 2, and the particle stage of write_contract2 as
    MODEL_STAGE[1]. Today's apeGmsh readers open it and ignore the particles."""
    m = contract1_model() if model is None else model
    make_synthetic.build(str(path))
    with h5py.File(path, "r+") as f:
        make_synthetic._ai(f["INFO"], "FORMAT_MINOR_VERSION", 2)
        write_particle_stage(f, "MODEL_STAGE[1]", m, contract2_outputs(m), step=1)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    write_contract1(PARTICLE_SET)
    write_contract2(PARTICLE_RESULTS)
    write_mixed(PARTICLE_MIXED)
    for path, check in ((PARTICLE_SET, mc.validate_particle_set),
                        (PARTICLE_RESULTS, mc.validate_particle_results),
                        (PARTICLE_MIXED, mc.validate_particle_results)):
        problems = check(str(path))
        print(f"wrote {path.relative_to(HERE.parent)}: "
              f"{'conformant' if not problems else problems}")


if __name__ == "__main__":
    main()
