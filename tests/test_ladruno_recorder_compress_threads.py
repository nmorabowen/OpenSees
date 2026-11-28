"""Ladruno recorder, off-thread compression and threaded reactions (WP-251).

The compressed DATA layout stores one step per chunk (WP-233), so each recorded step is
a set of complete chunks. WP-251 compresses them on worker threads (byte shuffle + zlib,
the exact pipeline the dataset's filters declare) and the solver thread hands the bytes
to H5Dwrite_chunk; `-compressThreads 0` keeps the old in-line HDF5 filter path. Gates:

  * the stored chunk BYTES and the decoded arrays are identical to the old path, for
    f64 and `-precision f32`;
  * a deck that exits without `wipe` still writes every step (the at-exit drain);
  * a killed run reads by the WP-233 rule, every committed row equal to a full run's;
  * reactionForce under `ladrunoThreads 4` (the threaded reaction loop) is
    bit-identical to the serial loop.

Speed: scratch bench in the WP-251 PR (recorder 124 -> 44 ms/step serial on 8192
stdBricks), not gated here.
"""
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest

from _testbed import ops

h5py = pytest.importorskip("h5py")
pytestmark = [pytest.mark.zone_a]

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "Ladruno_scripts" / "ladruno_recorder_tests"))
import ladruno_format as lf  # noqa: E402

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

N = 6          # 6 x 6 x 6 stdBricks: the stress slab (216 x 48 f64) is compressed
RESULTS = ["-N", "displacement", "velocity", "reactionForce", "-E", "stresses", "strains"]


def _block(n=N):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 3.0e7, 0.2, 2.0)   # WP-192 thread-safe

    def nid(i, j, k):
        return 1 + i + (n + 1) * (j + (n + 1) * k)

    for k in range(n + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                ops.node(nid(i, j, k), float(i), float(j), float(k))
                if k == 0:
                    ops.fix(nid(i, j, k), 1, 1, 1)
    e = 0
    for k in range(n):
        for j in range(n):
            for i in range(n):
                e += 1
                ops.element("stdBrick", e, nid(i, j, k), nid(i + 1, j, k),
                            nid(i + 1, j + 1, k), nid(i, j + 1, k), nid(i, j, k + 1),
                            nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1),
                            nid(i, j + 1, k + 1), 1)
    return [nid(i, j, n) for j in range(n + 1) for i in range(n + 1)]


def _transient(top):
    ops.timeSeries("Sine", 1, 0.0, 1.0e9, 0.05)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, 10.0, 3.0, -2.0)
    ops.constraints("Plain"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-12, 10); ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25); ops.analysis("Transient")


def _deck(path, steps, *extra, wipe=True):
    _transient(_block())
    ops.recorder("ladruno", path, *RESULTS, *extra)
    assert ops.analyze(steps, 0.01) == 0
    if wipe:
        ops.wipe()


def _stored(path):
    """{dataset: (decoded array, [stored chunk bytes] or None)} for DATA/TIME/STEP."""
    out = {}
    with h5py.File(path, "r") as f:
        def grab(name, obj):
            if isinstance(obj, h5py.Dataset) and name.rsplit("/", 1)[-1] in ("DATA", "TIME", "STEP"):
                raw = None
                if obj.compression:
                    raw = [obj.id.read_direct_chunk(obj.id.get_chunk_info(i).chunk_offset)
                           for i in range(obj.id.get_num_chunks())]
                out[name] = (obj[...], raw, obj.dtype.str)
        f.visititems(grab)
    return out


def _assert_same(a, b, need_compressed=True):
    sa, sb = _stored(a), _stored(b)
    assert sa.keys() == sb.keys()
    n_raw = 0
    for k in sa:
        va, ra, ta = sa[k]
        vb, rb, tb = sb[k]
        assert ta == tb, k
        if k.endswith("/STEP"):
            # commit tags keep counting across in-process decks: compare the axis shape
            va, vb = va - va[:1], vb - vb[:1]
        assert np.array_equal(va, vb), k
        assert ra == rb, f"{k}: stored chunk bytes differ"
        n_raw += len(ra or [])
    if need_compressed:
        assert n_raw > 0, "no compressed chunk was compared (the slab fell below the 4 KiB cut)"


@pytest.mark.parametrize("precision", ["f64", "f32"])
def test_off_thread_chunks_equal_the_hdf5_filter_path(tmp_path, precision):
    sync, asy = str(tmp_path / "sync.ladruno"), str(tmp_path / "async.ladruno")
    _deck(sync, 40, "-precision", precision, "-compressThreads", 0)
    _deck(asy, 40, "-precision", precision, "-compressThreads", 3)
    _assert_same(sync, asy)


def test_flush_cadence_and_auto_threads(tmp_path):
    # -flush 0 drains every step; the default (auto threads, 10 s cadence) lags and
    # drains on close. Same file contents either way.
    a, b = str(tmp_path / "f0.ladruno"), str(tmp_path / "auto.ladruno")
    _deck(a, 30, "-flush", 0, "-compressThreads", 0)
    _deck(b, 30)
    _assert_same(a, b)


_CHILD = r"""
import os, sys
sys.path.insert(0, os.environ["LADRUNO_TESTS"])
from _testbed import ops
import test_ladruno_recorder_compress_threads as t
out, mode, steps = sys.argv[1], sys.argv[2], int(sys.argv[3])
t._transient(t._block())
ops.recorder("ladruno", out, *t.RESULTS, "-flush", 0 if mode == "kill" else 10)
if mode == "kill":
    k = 0
    while True:
        assert ops.analyze(1, 0.01) == 0
        k += 1
        if k == 50:
            print("STEP 50", flush=True)
else:
    assert ops.analyze(steps, 0.01) == 0
    if mode == "wipe":
        ops.wipe()
print("DONE", flush=True)
"""


def _child(out, mode, steps=0):
    env = dict(os.environ, LADRUNO_TESTS=str(ROOT / "tests"), LADRUNO_OPENSEES_QUIET="1")
    env.pop("HDF5_USE_FILE_LOCKING", None)
    return subprocess.Popen([sys.executable, "-c", _CHILD, out, mode, str(steps)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            stdin=subprocess.DEVNULL, text=True, env=env)


def test_exit_without_wipe_writes_every_step(tmp_path):
    """Rows still in the compression queue are not in HDF5's chunk cache, so HDF5's own
    atexit close cannot save them; the recorder's at-exit hook (registered after HDF5's,
    so it runs first) compresses them on the exiting thread and writes them."""
    ref, nowipe = str(tmp_path / "ref.ladruno"), str(tmp_path / "nowipe.ladruno")
    for path, mode in ((ref, "wipe"), (nowipe, "exit")):
        o, e = _child(path, mode, 60).communicate(timeout=300)
        assert "DONE" in o, e[-2000:]
    _assert_same(ref, nowipe)
    with h5py.File(nowipe, "r") as f:
        st = f[[k for k in f if k.startswith("MODEL_STAGE")][0]]
        assert st["RESULTS/ON_NODES/DISPLACEMENT/TIME"].shape[0] == 60


def test_a_killed_async_run_reads_by_the_rule(tmp_path):
    out = str(tmp_path / "killed.ladruno")
    p = _child(out, "kill")
    try:
        line = p.stdout.readline()
        assert "STEP 50" in line, line + p.stderr.read()[-2000:]
        time.sleep(0.3)
    finally:
        p.kill()
        p.wait()
    try:
        f = lf.open_ladruno(out)
    except lf.StageLayoutKilledError:
        rep = str(tmp_path / "repaired.ladruno")
        lf.repair_copy(out, rep)
        f = h5py.File(rep, "r")
    with f:
        report = lf.stage_report(f)
        n = report[0]["committed"]
        assert n >= 50 - 8, report           # at most the in-flight rows behind the kill
        ref = str(tmp_path / "ref.ladruno")
        o, e = _child(ref, "wipe", n + 10).communicate(timeout=300)
        assert "DONE" in o, e[-2000:]
        with h5py.File(ref, "r") as fr:
            st, rs = f[report[0]["stage"]], fr[report[0]["stage"]]
            for g in lf._streaming_results(st):
                rel = g.name.split("/RESULTS/", 1)[1]
                for ds in ("DATA", "TIME", "STEP"):
                    assert np.array_equal(g[ds][:n], rs[f"RESULTS/{rel}"][ds][:n]), (rel, ds)


def test_threaded_reactions_are_bit_identical(tmp_path):
    serial, threaded = str(tmp_path / "t1.ladruno"), str(tmp_path / "t4.ladruno")
    try:
        ops.ladrunoThreads(1)
        _deck(serial, 30, "-compressThreads", 0)
        ops.ladrunoThreads(4)
        _deck(threaded, 30, "-compressThreads", 0)
    finally:
        ops.ladrunoThreads(1)
    _assert_same(serial, threaded)


def test_compress_threads_needs_an_int(tmp_path, capfd):
    _block(2)
    with pytest.raises(Exception):
        ops.recorder("ladruno", str(tmp_path / "bad.ladruno"), "-N", "displacement",
                     "-compressThreads", "many")
    assert "-compressThreads requires an int" in capfd.readouterr().err
