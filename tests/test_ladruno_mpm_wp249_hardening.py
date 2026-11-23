"""WP-249 (MPM hardening before phase 7, program #1075).

* #1095 -- a mode-A step-cut retry must start from exactly the state a fresh step of the cut
  size starts from, also when the steps around it change the particles' connectivity. The case
  is #1095's: 3D, quasi-static LoadControl, DruckerPrager (non-associated), gravity, the default
  `renormalize` empty-node treatment at cutoff 5e-2, one thread. Three attempts of a step are
  forced to fail AFTER the trial (`mpm debug failAttempt <k> 3 afterTrial`: the materials carry
  a trial state, as after a Newton failure), on two steps whose own preparation re-points the
  particles (a connectivity change at the failed attempt; asserted on the reference run's first
  sub-step), so each call cuts three times and then takes eight sub-steps of dt/8 from the
  re-pointed state. The run must equal (rel 1e-12) the run that takes those sub-steps directly;
  it runs in a child process, so a crash is an exit code, and no size or host warning may
  appear. A guard: #1095 itself was not reproduced (WP-249).
* #1106 item 2 -- `mpm debug failDump` is a process-wide switch; it must not outlive its host
  (`wipe`, `mpm reset`) and is re-asserted from the live host's own setting.
* #1103 item 2 -- a version-1 restart pack (the ladruno-v20261116 layout: version 2 without the
  particle body force) restores and continues as the uninterrupted run.
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import pytest

TESTS = Path(__file__).resolve().parent
sys.path.insert(0, str(TESTS))
from _testbed.subprocess_run import run_python_script  # noqa: E402

pytestmark = pytest.mark.zone_a


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


# --------------------------------------------------------------------------- #
#  #1095 -- cut retries across connectivity changes (child process)           #
# --------------------------------------------------------------------------- #
_CHILD = r'''
import json, os, sys
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
sys.path.insert(0, sys.argv[1])
from _testbed import ops

B, RHO, G_ACC, E, NU, N_CELLS = 1.0, 1000.0, 10.0, 5.0e4, 0.3, 4
K = E / (3.0 * (1.0 - 2.0 * NU)); G = E / (2.0 * (1.0 + NU))
DT = 0.09375


def build():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    # #1095's block: DruckerPrager, non-associated (rhoBar 0), yields under K0 loading
    ops.nDMaterial("DruckerPrager", 1, K, G, 400.0, 0.21, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    h = B / N_CELLS
    ops.mpm("grid", 1, "-origin", 0.0, 0.0, 0.0, "-h", h, "-cells", N_CELLS, N_CELLS, N_CELLS + 2,
            "-basis", "gimp")
    ops.mpm("mode", "quasiStatic")
    ops.mpm("body", 1, "box", 0.0, 0.0, 0.0, B, B, B, "-ppc", 2, "-mat", 1, "-rho", RHO)
    ops.mpm("fix", "xmin", "x"); ops.mpm("fix", "xmax", "x")
    ops.mpm("fix", "ymin", "y"); ops.mpm("fix", "ymax", "y")
    ops.mpm("fix", "zmin", "xyz")
    ops.mpm("gravity", 0.0, 0.0, -G_ACC)
    ops.mpm("massCutoff", 5.0e-2)
    ops.mpm("emptyNodes", "renormalize")
    ops.mpm("cut", 2.0, 4)
    ops.timeSeries("Path", 7, "-time", 0.0, 1.0, 1.0e9, "-values", 0.0, 1.0, 1.0)
    ops.mpm("timeSeries", 7)
    ops.constraints("Plain"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 50); ops.algorithm("Newton")
    ops.integrator("LoadControl", DT); ops.analysis("Static")


def info(w):
    return int(ops.mpm("info", w))


def state():
    return {f: [float(v) for v in ops.mpm("particleData", f)] for f in ("x", "sig", "F")}


out = {}
# the retried run: two forced cut sequences, each followed by eight dt/8 sub-steps. The armed
# steps (2, then 1 + 8 + 1 + 1 = 11) are ones whose own prepareStep re-points the particles
build()
ops.mpm("analyze", 1, DT)
conn, cut_retries = [], []
for k in (2, 11):                         # epochSteps + 1 of the armed call
    ops.mpm("debug", "failAttempt", k, 3, "afterTrial")
    c0, r0 = info("connectivityChanges"), info("retries")
    ops.mpm("analyze", 1, DT)
    conn.append(info("connectivityChanges") - c0)
    cut_retries.append(info("retries") - r0)
    ops.mpm("analyze", 1, DT)
out["retry"] = state()
out["retries_in_cut_calls"] = cut_retries
out["retries"] = info("retries")
out["accepted"] = info("acceptedSteps")
out["conn_in_cut_calls"] = conn
# the reference: the same sub-steps, taken directly; the first sub-step's connectivity change
# is the armed step's own (its prepareStep runs from the same committed particles)
build()
ops.mpm("analyze", 1, DT)
first = []
for _ in range(2):
    c0 = info("connectivityChanges")
    ops.mpm("analyze", 1, DT / 8.0)
    first.append(info("connectivityChanges") - c0)
    ops.mpm("analyze", 7, DT / 8.0)
    ops.mpm("analyze", 1, DT)
out["reference_first_substep_conn"] = first
out["reference"] = state()
out["reference_retries"] = info("retries")
print("WP249-RESULT " + json.dumps(out))
'''

# text a #1095-class failure prints (Matrix / Vector size checks, the element's host-bug lines)
_BAD = ("incompatible sizes", "not of correct size", "MPM host bug", "out of range",
        "is not in the Domain", "access violation", "Fatal Python error")


def test_cut_retry_across_connectivity_changes_3d():
    _ops()                                     # skip cleanly where there is no build
    rc, text = run_python_script(_CHILD, cwd=str(TESTS), timeout=900, argv=(str(TESTS),))
    assert rc == 0, f"the child died (rc {rc}):\n{text[-3000:]}"
    for bad in _BAD:
        assert bad not in text, f"'{bad}' printed during the cut retries:\n{text[-3000:]}"
    line = [ln for ln in text.splitlines() if ln.startswith("WP249-RESULT ")]
    assert line, text[-3000:]
    out = json.loads(line[-1][len("WP249-RESULT "):])
    assert out["retries"] == 6, out["retries"]                 # 3 forced failures x 2 calls
    assert out["retries_in_cut_calls"] == [3, 3], out["retries_in_cut_calls"]   # armed where meant
    assert out["reference_retries"] == 0
    assert out["accepted"] == 1 + 2 * (8 + 1), out["accepted"]
    # "across a connectivity change": each armed step's own preparation re-points the particles,
    # so its failed attempts and its retries all start from a changed connectivity
    assert all(c >= 1 for c in out["reference_first_substep_conn"]), out["reference_first_substep_conn"]
    assert all(c >= 1 for c in out["conn_in_cut_calls"]), out["conn_in_cut_calls"]
    for f in ("x", "sig", "F"):
        a, b = out["retry"][f], out["reference"][f]
        assert len(a) == len(b)
        den = max(abs(v) for v in b) or 1.0
        assert max(abs(x - y) for x, y in zip(a, b)) / den <= 1e-12, f


# --------------------------------------------------------------------------- #
#  #1106 item 2 -- the failDump switch dies with its host                     #
# --------------------------------------------------------------------------- #
def _small_2d(ops):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e3, 0.3)
    ops.mpm("grid", 1, "-origin", 0.0, 0.0, "-h", 0.5, "-cells", 4, 4, "-basis", "gimp")
    ops.mpm("mode", "quasiStatic")
    ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, 1.0, "-ppc", 2, "-mat", 1, "-rho", 1.0)
    ops.mpm("fix", "ymin", "xy")
    ops.mpm("gravity", 0.0, -1.0)
    ops.timeSeries("Path", 7, "-time", 0.0, 1.0, 1.0e9, "-values", 0.0, 1.0, 1.0)
    ops.mpm("timeSeries", 7)
    ops.constraints("Plain"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 25); ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.25); ops.analysis("Static")


def test_faildump_switch_does_not_outlive_its_host():
    ops = _ops()
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
        _small_2d(ops)
        ops.mpm("debug", "failDump", str(Path(td) / "dump.txt"), 3)
        assert int(ops.mpm("info", "failDump")) == 1
        ops.mpm("analyze", 1, 0.25)
        assert int(ops.mpm("info", "failDump")) == 1
        ops.mpm("reset")                          # teardown clears it ...
        assert int(ops.mpm("info", "failDump")) == 0
        ops.mpm("analyze", 1, 0.25)               # ... and this host's analyze re-asserts it
        assert int(ops.mpm("info", "failDump")) == 1
        _small_2d(ops)                            # `wipe`: the next model pays nothing
        assert int(ops.mpm("info", "failDump")) == 0
        ops.mpm("analyze", 1, 0.25)
        assert int(ops.mpm("info", "failDump")) == 0
        ops.wipe()


# --------------------------------------------------------------------------- #
#  #1103 item 2 -- a version-1 restart pack restores                          #
# --------------------------------------------------------------------------- #
def _analysis_2d(ops):
    ops.constraints("Plain"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 25); ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.25); ops.analysis("Static")


def _save_restore(ops, db, tag):
    ops.database("File", db)
    ops.save(tag)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.database("File", db)
    ops.restore(tag)
    _analysis_2d(ops)


def _sig(ops):
    return [float(v) for v in ops.mpm("particleData", "sig")]


def _close(a, b, rel=1e-10):
    scale = max(abs(v) for v in b) or 1.0
    return len(a) == len(b) and max(abs(x - y) for x, y in zip(a, b)) <= rel * scale


def test_restart_pack_v1_restores_and_continues():
    ops = _ops()
    _small_2d(ops)
    for _ in range(4):
        ops.mpm("analyze", 1, 0.25)
    ref = _sig(ops)
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
        _small_2d(ops)
        for _ in range(2):
            ops.mpm("analyze", 1, 0.25)
        ops.mpm("debug", "packVersion", 1)
        _save_restore(ops, str(Path(td) / "v1db"), 2)
        for _ in range(2):
            ops.mpm("analyze", 1, 0.25)
        assert _close(_sig(ops), ref)
        ops.wipe()


def test_restart_pack_v1_falls_back_when_it_cannot_hold_the_state():
    """A particle body force has no place in version 1: the writer falls back to the current
    layout, so the restored run keeps the body force and continues as the uninterrupted one."""
    ops = _ops()

    def loaded():
        _small_2d(ops)
        n = int(ops.mpm("info", "numParticles"))
        ops.mpm("particleBodyForce", *([0.0, -0.5] * n))

    loaded()
    for _ in range(3):
        ops.mpm("analyze", 1, 0.25)
    ref = _sig(ops)
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
        loaded()
        ops.mpm("analyze", 1, 0.25)
        ops.mpm("debug", "packVersion", 1)
        _save_restore(ops, str(Path(td) / "v1bf"), 1)
        for _ in range(2):
            ops.mpm("analyze", 1, 0.25)
        assert _close(_sig(ops), ref)
        ops.wipe()
