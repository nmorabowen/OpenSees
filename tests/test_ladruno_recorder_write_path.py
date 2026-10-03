"""Ladruno recorder write path (WP-164): same results, written differently.

WP-164 keeps every group, dataset and shape a reader depends on and changes how they
are written (Ladruno_implementation/164_ladruno_recorder_write_path.md):

  * DATA/TIME/STEP stay open for the stage; the file is flushed on a wall-clock
    cadence (`-flush <s>`, default 10 s; `-flush 0` = every step, the old path).
    Gate: the two produce IDENTICAL arrays.
  * chunks target ~1 MiB; a slab above 1 MiB tiles the id axis (one entity's
    history no longer inflates the whole dataset).
  * `-compress <0..9>` sets the deflate level; 0 drops the shuffle+deflate filters.
  * `-envelope` creates its datasets once and overwrites them in place (it used to
    delete and recreate every group every step, so the file grew with the step
    count). Gate: the file size does not grow with the number of steps.

Speed is measured by Ladruno_scripts/ladruno_recorder_tests/bench/recorder_bench.py,
not here (wall time on a shared CI runner is not a gate).
"""
import os

import numpy as np
import pytest

from _testbed import ops

h5py = pytest.importorskip("h5py")
pytestmark = [pytest.mark.zone_a]

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")


def _brick_block(nx, ny, nz):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 3.0e7, 0.2)

    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                ops.node(nid(i, j, k), float(i), float(j), float(k))
                if k == 0:
                    ops.fix(nid(i, j, k), 1, 1, 1)
    e = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                e += 1
                ops.element("stdBrick", e, nid(i, j, k), nid(i + 1, j, k),
                            nid(i + 1, j + 1, k), nid(i, j + 1, k), nid(i, j, k + 1),
                            nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1),
                            nid(i, j + 1, k + 1), 1)
    return [nid(i, j, nz) for j in range(ny + 1) for i in range(nx + 1)]


def _run_static(top, steps):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in top:
        ops.load(n, 1.0, 0.0, -1.0)
    ops.constraints("Plain"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-8, 5); ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0 / steps); ops.analysis("Static")
    assert ops.analyze(steps) == 0
    ops.wipe()


def _stage(f):
    return f[[k for k in f if k.startswith("MODEL_STAGE")][0]]


def _arrays(path):
    out = {}
    with h5py.File(path, "r") as f:
        st = _stage(f)

        def grab(name, obj):
            if isinstance(obj, h5py.Dataset) and name.rsplit("/", 1)[-1] in ("DATA", "TIME", "STEP"):
                out[name] = obj[...]
        st["RESULTS"].visititems(grab)
    return out


def test_flush_cadence_writes_identical_data(tmp_path):
    a, b = str(tmp_path / "every.ladruno"), str(tmp_path / "cadence.ladruno")
    for path, flush in ((a, 0), (b, 10)):
        top = _brick_block(4, 4, 3)
        ops.recorder("ladruno", path, "-N", "displacement", "reactionForce",
                     "-E", "stresses", "-flush", flush)
        _run_static(top, 6)
    da, db = _arrays(a), _arrays(b)
    assert da.keys() == db.keys() and da, (da.keys(), db.keys())
    for k in da:
        assert da[k].shape[0] == 6, (k, da[k].shape)
        if k.endswith("/STEP"):
            # STEP is the domain commitTag, which keeps counting across `wipe`
            # in one interpreter (run 1: 0..5, run 2: 6..11) — compare increments.
            np.testing.assert_array_equal(np.diff(da[k]), np.diff(db[k]), err_msg=k)
        else:
            np.testing.assert_array_equal(da[k], db[k], err_msg=k)


def test_large_slab_tiles_the_id_axis(tmp_path):
    """14x14x14 bricks x 48 GP-stress columns x 8 B = 1.05 MiB per step > 1 MiB."""
    path = str(tmp_path / "tile.ladruno")
    top = _brick_block(14, 14, 14)
    ops.recorder("ladruno", path, "-E", "stresses")
    _run_static(top, 3)
    with h5py.File(path, "r") as f:
        data = next(iter(_stage(f)["RESULTS/ON_ELEMENTS/stresses"].values()))["DATA"]
        assert data.shape == (3, 14 ** 3, 48), data.shape
        assert data.chunks[1] < data.shape[1], data.chunks          # ids tiled
        assert data.chunks[0] * data.chunks[1] * data.chunks[2] * 8 <= (1 << 20) + 48 * 8
        k = data.shape[1] // 2
        hist = data[:, k, :]
        assert np.all(np.isfinite(hist)) and np.any(hist != 0.0)


@pytest.mark.parametrize("level,expect", [(0, None), (1, "gzip"), (9, "gzip")])
def test_compress_option_sets_the_filter(tmp_path, level, expect):
    path = str(tmp_path / f"c{level}.ladruno")
    top = _brick_block(2, 2, 2)
    ops.recorder("ladruno", path, "-N", "displacement", "-compress", level)
    _run_static(top, 2)
    with h5py.File(path, "r") as f:
        d = _stage(f)["RESULTS/ON_NODES/DISPLACEMENT/DATA"]
        assert d.compression == expect, d.compression
        if expect:
            assert d.compression_opts == level
            assert d.shuffle


def test_envelope_file_does_not_grow_with_steps(tmp_path):
    sizes = {}
    for steps in (3, 30):
        path = str(tmp_path / f"env{steps}.ladruno")
        top = _brick_block(6, 6, 3)
        ops.recorder("ladruno", path, "-N", "displacement", "-E", "stresses", "-envelope")
        _run_static(top, steps)
        sizes[steps] = os.path.getsize(path)
        with h5py.File(path, "r") as f:
            env = _stage(f)["RESULTS/ENVELOPES"]
            assert "ON_NODES/DISPLACEMENT/ABSMAX" in env
            assert any("COLUMN_MAP" in g for g in env["ON_ELEMENTS/stresses"].values())
    # in place: same datasets rewritten; the old delete/recreate grew ~linearly
    assert sizes[30] < 1.10 * sizes[3], sizes
