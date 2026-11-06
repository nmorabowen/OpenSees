"""WP-195 (item 2) -- zeroLength springs and dashpots on the WP-192 opt-in.

The Tier 2 decks (2A-2D) carry their soil as zeroLength springs and dashpots:
`zeroLength -mat P P P -dir 1 2 3 -doRayleigh 0 -orient ...` with
P = Parallel(Elastic 1.0, Viscous 1.0 1.0) -factors k c. WP-192 left zeroLength
un-audited (serial in place). WP-195 audits it: its class-wide matrices and
vectors become per thread (read through ladrunoM()/ladrunoV()), and
ZeroLength, ParallelMaterial and ViscousMaterial opt in (Elastic already had).

Test: a slab of ASDShellQ4 on such springs/dashpots, static + transient, has
the SAME full-field digest (displacement, velocity, reaction of every node and a
spring's force, every step) at 1 / 2 / 4 threads, and the loops report every
element threaded -- no serial-in-place element.

See `Ladruno_implementation/195_ladruno_brick_tangent_cache_adr.md`.
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


def _threads_available():
    hw = ops.ladrunoThreads(1024)
    ops.ladrunoThreads(1)
    return [n for n in (1, 2, 4) if n <= hw]


def _shell_on_springs(nx=6, ny=4):
    """A 3D slab of elastic ASDShellQ4 on zeroLength springs and dashpots, as
    the 2C foundation: zeroLength -mat P P P -dir 1 2 3 -doRayleigh 0
    -orient ..., with P = Parallel(Elastic 1.0, Viscous 1.0 1.0) -factors k c."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.section("ElasticMembranePlateSection", 3, 25000.0, 0.2, 300.0, 2.4e-9)
    ops.uniaxialMaterial("Elastic", 9001, 1.0)
    ops.uniaxialMaterial("Viscous", 9002, 1.0, 1.0)
    h = 1000.0
    nid, gid = {}, {}
    t = 1
    for i in range(nx + 1):
        for j in range(ny + 1):
            nid[(i, j)] = t
            ops.node(t, i * h, j * h, 0.0)
            gid[(i, j)] = 100000 + t
            ops.node(100000 + t, i * h, j * h, 0.0)
            ops.fix(100000 + t, 1, 1, 1, 1, 1, 1)
            ops.mass(t, 0.02, 0.02, 0.02, 0.0, 0.0, 0.0)
            t += 1
    e = 1
    for i in range(nx):
        for j in range(ny):
            ops.element("ASDShellQ4", e, nid[(i, j)], nid[(i + 1, j)],
                        nid[(i + 1, j + 1)], nid[(i, j + 1)], 3)
            e += 1
    m = 10
    for (i, j), n in sorted(nid.items()):
        tags = []
        for d, (k, c) in enumerate(((2.2e4 + 10.0 * i, 412.0),
                                    (2.3e4 + 7.0 * j, 413.0),
                                    (1.2e5, 727.0 + i))):
            ops.uniaxialMaterial("Parallel", m, 9001, 9002, "-factors", k, c)
            tags.append(m)
            m += 1
        ops.element("zeroLength", 50000 + n, gid[(i, j)], n, "-mat", *tags,
                    "-dir", 1, 2, 3, "-doRayleigh", 0,
                    "-orient", 1.0, 0.0, 0.0, 0.0, 1.0, 0.0)
        # rotations restrained (the springs carry translations only)
        ops.fix(n, 0, 0, 0, 0, 0, 1)
    return nid


def _run_springs(nthreads):
    nid = _shell_on_springs()
    ops.ladrunoThreads(nthreads)
    try:
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for n in nid.values():
            ops.load(n, 0.0, 0.0, -50.0, 0.0, 0.0, 0.0)
        ops.constraints("Plain")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("NormDispIncr", 1.0e-8, 30, 0)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 0.5)
        ops.analysis("Static")
        hh = hashlib.md5()
        rcs = []

        def snap():
            ops.reactions()
            for k in sorted(nid.values()):
                hh.update(("%d|%s|%s|%s" % (k, repr(ops.nodeDisp(k)),
                                            repr(ops.nodeVel(k)),
                                            repr(ops.nodeReaction(k)))).encode())
            hh.update(repr(list(ops.eleResponse(50001, "force"))).encode())

        for _ in range(2):
            rcs.append(ops.analyze(1))
            snap()
        ops.loadConst("-time", 0.0)
        ops.wipeAnalysis()
        ops.rayleigh(0.1, 0.0, 0.002, 0.0)
        ops.timeSeries("Path", 2, "-dt", 0.01,
                       "-values", *[0.0, 500.0, -800.0, 1000.0, -500.0, 300.0, 0.0] * 2)
        ops.pattern("UniformExcitation", 2, 1, "-accel", 2)
        ops.pattern("UniformExcitation", 3, 2, "-accel", 2, "-fact", 0.6)
        ops.constraints("Plain")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("NormDispIncr", 1.0e-8, 30, 0)
        ops.algorithm("KrylovNewton")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        for _ in range(12):
            rc = ops.analyze(1, 0.01)
            rcs.append(rc)
            snap()
            if rc != 0:
                break
        return hh.hexdigest(), rcs
    finally:
        ops.ladrunoThreads(1)


def test_zerolength_springs_dashpots_threaded_equal_serial(capfd):
    counts = _threads_available()
    if len(counts) < 2:
        pytest.skip("single hardware thread")
    capfd.readouterr()
    digests, rcs_all = {}, {}
    for n in counts:
        d, rcs = _run_springs(n)
        assert all(rc == 0 for rc in rcs), (n, rcs)
        digests[n] = d
        rcs_all[n] = rcs
    text = capfd.readouterr().err
    assert len(set(digests.values())) == 1, digests
    # every element threaded: no zeroLength left serial in place
    assert "PARTLY threaded" not in text, text
    assert "run serially in place" not in text, text
    assert "element update loop THREADED" in text, text
    assert "element tangent loop threaded" in text, text
    assert "element residual loop threaded" in text, text
