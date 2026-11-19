"""WP-239 -- ASDSteel1D -buckling -implex on the WP-192 threaded loops.

ASDSteel1D opts into the threaded element loops whenever `-slip` is off. With
`-buckling -implex`, commitState re-runs the buckling RVE implicitly and keeps
the Newton history of that correction (U_iterative, U_iterative_el,
Rnorm_iterative) to pick the best iterate after 10 iterations. Before WP-239
those were plain function-scope statics: one array for every ASDSteel1D in the
process, so two points committing on different threads could copy back each
other's iterate. They are now thread_local (bit-identical serially).

Test: many zeroLength springs, each its own ASDSteel1D -buckling -implex
instance, driven through different cyclic histories with tolerances tight
enough that the implicit corrections run long and take the min-search branch.
The full response digest must be identical at 1 / 2 / 4 threads, and no
spring may be left serial in place.
"""
import hashlib
import os
import sys

import pytest

_DIST_BIN = os.environ.get("LADRUNO_DIST_BIN") or os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "bin")
if os.path.isdir(_DIST_BIN) and _DIST_BIN not in sys.path:
    sys.path.insert(0, _DIST_BIN)

from _testbed import ops                        # noqa: E402
from _testbed.subprocess_run import run_python_script   # noqa: E402


_OPENMP_PROBE = (
    "import sys; sys.path.insert(0, %r)\n"
    "import opensees as ops\n"
    "print('LADRUNO_PROBE_HW', ops.ladrunoThreads(1 << 20))\n"
)
_OPENMP_ABSENT_MARKER = "built WITHOUT LADRUNO_OPENMP"


def _probe():
    try:
        code, out = run_python_script(_OPENMP_PROBE % _DIST_BIN)
    except Exception:
        return True, ""
    if code != 0 or "LADRUNO_PROBE_HW" not in out:
        return True, ""
    if _OPENMP_ABSENT_MARKER in out:
        return False, "binary built with LADRUNO_OPENMP=OFF: nothing to thread"
    return True, ""


_OPENMP_IN, _SKIP_WHY = _probe()

pytestmark = [
    pytest.mark.zone_a,
    pytest.mark.skipif(not _OPENMP_IN, reason=_SKIP_WHY),
]

_NSPRING = 96
_HIST = [0.0, 0.004, -0.006, 0.008, -0.010, 0.006, -0.004, 0.0]


def _threads_available():
    hw = ops.ladrunoThreads(1024)
    ops.ladrunoThreads(1)
    return [n for n in (1, 2, 4) if n <= hw]


def _run(nthreads, steps_per_leg=6):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    # -tolU/-tolR at 1e-12: tight enough that the implicit corrections run
    # long and use the 10-slot history and its min-search, loose enough that
    # the explicit trial still converges. On the pre-WP-239 binary this deck
    # gave a different digest on every 4-thread run, none equal to serial.
    ops.uniaxialMaterial("ASDSteel1D", 1, 200000.0, 400.0, 500.0, 0.10,
                         "-implex", "-buckling", 120.0, 8.0,
                         "-max_iter", 12, "-tolU", 1.0e-12, "-tolR", 1.0e-12)
    for k in range(_NSPRING):
        a, b = 2 * k + 1, 2 * k + 2
        ops.node(a, 0.0)
        ops.node(b, 0.0)
        ops.fix(a, 1)
        ops.element("zeroLength", k + 1, a, b, "-mat", 1, "-dir", 1)
    ops.ladrunoThreads(nthreads)
    try:
        ops.timeSeries("Path", 1, "-time", *[float(i) for i in range(len(_HIST))],
                       "-values", *_HIST)
        ops.pattern("Plain", 1, 1)
        for k in range(_NSPRING):
            # a different history per spring (amplitude 0.6 .. 1.4)
            ops.sp(2 * k + 2, 1, 0.6 + 0.8 * k / (_NSPRING - 1))
        ops.constraints("Penalty", 1.0e12, 1.0e12)
        ops.numberer("Plain")
        ops.system("BandGeneral")
        ops.test("NormDispIncr", 1.0e-10, 20, 0)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0 / steps_per_leg)
        ops.analysis("Static")
        hh = hashlib.md5()
        rcs = []
        for _ in range(steps_per_leg * (len(_HIST) - 1)):
            rc = ops.analyze(1)
            rcs.append(rc)
            for e in range(1, _NSPRING + 1):
                hh.update(repr(list(ops.eleResponse(e, "force"))).encode())
            if rc != 0:
                break
        return hh.hexdigest(), rcs
    finally:
        ops.ladrunoThreads(1)


def test_asdsteel1d_buckling_implex_threaded_equal_serial(capfd):
    counts = _threads_available()
    if len(counts) < 2:
        pytest.skip("single hardware thread")
    capfd.readouterr()
    digests = {}
    for n in counts:
        d, rcs = _run(n)
        assert all(rc == 0 for rc in rcs), (n, rcs)
        digests[n] = d
    text = capfd.readouterr().err
    assert len(set(digests.values())) == 1, digests
    # the springs really ran on the threads: all 96 in the update loop, and all
    # 96 spring FEs in the tangent/residual loops (the serial-in-place FEs there
    # are the Penalty constraint FEs, not springs)
    assert "element update loop THREADED" in text and "(%d elements)" % _NSPRING in text, text
    assert "PARTLY threaded" not in text, text
    assert "%d of %d FE_Elements computed in parallel" % (_NSPRING, 3 * _NSPRING) in text, text
    assert "first: FE %d, a constraint FE" % _NSPRING in text, text
