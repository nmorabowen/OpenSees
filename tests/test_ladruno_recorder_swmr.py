"""LadrunoRecorder under HDF5 SWMR (WP-233, #1082): a killed run stays readable.

Schema §7b.2 (Ladruno_implementation/ladruno_schema_v1.md). The recorder writes the latest
HDF5 format, in SWMR-write mode except during a layout step, with a streaming layout whose
chunks never move (SWMR never reuses freed space). These tests pin, on the real binary:

  * the layout: `-compress 1` (default) gives one deflated step per chunk for a slab of
    >= 4 KiB, a smaller slab or `-compress 0` an unfiltered ~64 KiB chunk;
  * no SWMR bloat: a `-flush 0` run's file stays near its payload (the WP-164 layout under
    SWMR grew 2.6 -> 206 MB);
  * a closed file opens plainly (consistency flags 0x00, version-3 superblock);
  * a run killed mid-analysis (TerminateProcess) reads by the reader rule and the crash
    rule (`ladruno_format.open_ladruno` / `repair_copy` / `stage_report`), with its rows
    equal to a completed run's.

The exhaustive evidence is elsewhere: every byte-level kill point of an h5py mirror
(Ladruno_scripts/ladruno_recorder_tests/recorder_crash_replay.py) and a kill stress per
lifecycle window (recorder_binary_kill_stress.py).
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Ladruno_scripts" / "ladruno_recorder_tests"))
import ladruno_format as lf  # noqa: E402


def _chain(n, mode="displacement"):
    """A 2-D cantilever of n elasticBeamColumn (n+1 nodes), a sine tip load, Newmark."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for i in range(n + 1):
        ops.node(i + 1, 0.0, 0.1 * i)
        if i:
            ops.mass(i + 1, 1.0, 1.0, 0.01)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("Linear", 1)
    for i in range(n):
        ops.element("elasticBeamColumn", i + 1, i + 1, i + 2, 0.01, 2.0e8, 1.0e-5, 1)


def _transient(n, steps):
    ops.timeSeries("Sine", 1, 0.0, 1.0e9, 0.37)
    ops.pattern("Plain", 1, 1)
    ops.load(n + 1, 10.0, 0.0, 0.0)
    ops.constraints("Plain"); ops.numberer("RCM"); ops.system("BandGen")
    ops.test("NormDispIncr", 1.0e-12, 10); ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25); ops.analysis("Transient")
    assert ops.analyze(steps, 0.001) == 0
    ops.wipe()


def _data(f, name="DISPLACEMENT"):
    st = f[[k for k in f if k.startswith("MODEL_STAGE")][0]]
    return st[f"RESULTS/ON_NODES/{name}/DATA"]


@pytest.mark.parametrize("n,compress,deflated", [
    (300, 1, True),     # 301 nodes x 2 (2-D displacement) x 8 B = 4.8 KiB: one deflated step per chunk
    (20, 1, False),     # 21 nodes: 336 B, under 4 KiB: unfiltered ~64 KiB chunks
    (300, 0, False),    # -compress 0: unfiltered
])
def test_streaming_layout(tmp_path, n, compress, deflated):
    path = str(tmp_path / f"lay{n}_{compress}.ladruno")
    _chain(n)
    ops.recorder("ladruno", path, "-N", "displacement", "-compress", compress)
    _transient(n, 30)
    with h5py.File(path, "r") as f:
        d = _data(f)
        assert d.shape == (30, n + 1, 2)
        if deflated:
            assert d.compression == "gzip" and d.chunks[0] == 1, (d.compression, d.chunks)
        else:
            assert d.compression is None and d.chunks[0] > 1, (d.compression, d.chunks)
            assert d.chunks[0] * d.chunks[1] * d.chunks[2] * 8 <= (64 << 10)


def test_closed_file_opens_plainly_latest_format(tmp_path):
    path = str(tmp_path / "closed.ladruno")
    _chain(20)
    ops.recorder("ladruno", path, "-N", "displacement")
    _transient(20, 10)
    assert lf.superblock_flags(path) == 0
    with open(path, "rb") as fh:
        assert fh.read(9)[8] == 3           # version-3 superblock: the latest format
    with h5py.File(path, "r") as f:         # a plain open, no SWMR-read needed
        assert _data(f).shape[0] == 10
    assert lf.validate(path) == []


@pytest.mark.parametrize("compress", [1, 0])
def test_flush_every_step_does_not_bloat(tmp_path, compress):
    """The WP-164 layout under SWMR leaked the space of every recompressed partial chunk
    (a 300-step run at -flush 0: 2.6 -> 206 MB). Chunks that never move do not leak."""
    path = str(tmp_path / f"bloat{compress}.ladruno")
    n, steps = 300, 300
    _chain(n)
    ops.recorder("ladruno", path, "-N", "displacement", "-flush", 0, "-compress", compress)
    _transient(n, steps)
    payload = steps * (n + 1) * 2 * 8
    assert os.path.getsize(path) < 2.0 * payload + (1 << 20), (os.path.getsize(path), payload)


_KILL_CHILD = r"""
import os, sys
sys.path.insert(0, os.environ["LADRUNO_TESTS"])
from _testbed import ops
import test_ladruno_recorder_swmr as t
out, every, limit = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
t._chain(60)
ops.recorder("ladruno", out, "-N", "displacement", "velocity", "-E", "force", "-flush", 0)
ops.timeSeries("Sine", 1, 0.0, 1.0e9, 0.37); ops.pattern("Plain", 1, 1); ops.load(61, 10.0, 0.0, 0.0)
ops.constraints("Plain"); ops.numberer("RCM"); ops.system("BandGen")
ops.test("NormDispIncr", 1.0e-12, 10); ops.algorithm("Linear")
ops.integrator("Newmark", 0.5, 0.25); ops.analysis("Transient")
k = 0
while not limit or k < limit:
    assert ops.analyze(1, 0.001) == 0
    k += 1
    if k == 50:
        print("STEP 50", flush=True)
    if every and k % every == 0:
        ops.domainChange()
ops.wipe()
print("DONE", flush=True)
"""


def _child(out, every, limit):
    env = dict(os.environ, LADRUNO_TESTS=str(ROOT / "tests"), LADRUNO_OPENSEES_QUIET="1")
    env.pop("HDF5_USE_FILE_LOCKING", None)
    return subprocess.Popen([sys.executable, "-c", _KILL_CHILD, out, str(every), str(limit)],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            stdin=subprocess.DEVNULL, text=True, env=env)


def _reference(path, steps, every):
    """The same deck run to completion in a FRESH process (commit tags and time start
    where the killed child's did; the analysis is deterministic: a killed run is its
    prefix)."""
    o, e = _child(path, every, steps).communicate()
    assert "DONE" in o, e[-2000:]


@pytest.mark.parametrize("every", [0, 25])
def test_a_killed_run_reads_by_the_rule(tmp_path, every):
    """TerminateProcess mid-run (in one stage, or with a new stage every 25 steps): the
    file reads by schema §7b.2 and every committed row equals a completed run's."""
    out = str(tmp_path / "killed.ladruno")
    p = _child(out, every, 0)
    try:
        line = p.stdout.readline()
        assert "STEP 50" in line, line + p.stderr.read()[-2000:]
        time.sleep(0.4)
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
        assert report and all("committed" in r for r in report[:-1]), report
        ref = str(tmp_path / "ref.ladruno")
        n_total = sum(r["committed"] or 0 for r in report if "committed" in r)
        _reference(ref, n_total + 60, every)
        with h5py.File(ref, "r") as fr:
            for r in report:
                if "committed" not in r:
                    continue                      # the stage being laid out (torn)
                n = r["committed"]
                st, rs = f[r["stage"]], fr[r["stage"]]
                for g in lf._streaming_results(st):
                    rel = g.name.split("/RESULTS/", 1)[1]
                    rg = rs[f"RESULTS/{rel}"]
                    for ds in ("DATA", "TIME", "STEP"):
                        assert np.array_equal(g[ds][:n], rg[ds][:n]), (r["stage"], rel, ds)
        assert n_total >= 50
    if lf.superblock_flags(out) == 0x05:
        # the validator checks closed files: a killed one is refused, loudly
        assert any("FILE-KILLED" in e for e in lf.validate(out))
