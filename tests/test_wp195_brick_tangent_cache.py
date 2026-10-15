"""WP-195 (item 1) -- stdBrick element tangent cache.

1. stdBrick recomputed its 24x24 stiffness (8 Gauss points x 8 x 8 B^T D B
   products) at every tangent call, although in the Paper 1 decks the soil is
   ElasticIsotropic. The element now caches the stiffness together with the
   RAW BITS of the eight material tangents and of the nodal coordinates that
   built it, and reuses it only when, at the next tangent call, those bits are
   unchanged. Elasticity is observed through the material (its getTangent()),
   never assumed, so the cache is bit-identical to the uncached path BY
   CONSTRUCTION -- the tests below hold it to equality of full-field digests,
   never to a tolerance:

     * elastic block, static + transient            cache on == cache off
     * plastic block (J2, material state changes)    cache on == cache off,
                                                     and the cache rebuilt
     * parameter update both ways (updateParameter on the material, and
       setParameter -ele through the element)        on == off, rebuilt
     * geometry change (setNodeCoord)                on == off, rebuilt
     * revertToStart (reset) then the same history   on == off == first run,
                                                     cache state 0 after reset
     * LADRUNO_BRICK_TANGENT_CACHE=0                 the cache is off (the
                                                     element response says so)

   The switch is read at element construction from the environment.

See `Ladruno_implementation/195_ladruno_brick_tangent_cache_adr.md`.
"""
import ast
import hashlib
import os
import subprocess
import sys

import pytest

_DIST_BIN = os.environ.get("LADRUNO_DIST_BIN") or os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "bin")
if os.path.isdir(_DIST_BIN) and _DIST_BIN not in sys.path:
    sys.path.insert(0, _DIST_BIN)

from _testbed import ops                        # noqa: E402
from _testbed.subprocess_run import pinned_child   # noqa: E402

pytestmark = [pytest.mark.zone_a]

_ENV = "LADRUNO_BRICK_TANGENT_CACHE"
_HERE = os.path.dirname(os.path.abspath(__file__))


def _set_cache(on):
    """No-op in-process: the switch is read with getenv() at element
    construction, and on Windows Python's os.environ does not reach the C
    runtime's environment block. Each setting runs in a fresh child (_child)."""


def _child(fn, *args, on=True):
    """Run module function `fn(*args)` in a fresh interpreter on the PARENT's
    engine (WP-176 pinning) with the cache switched on or off, and return its
    (repr-able) result."""
    env = dict(os.environ)
    if on:
        env.pop(_ENV, None)
    else:
        env[_ENV] = "0"
    code = ("import sys; sys.path.insert(0, %r); sys.path.insert(0, %r)\n"
            "import test_wp195_brick_tangent_cache as m\n"
            "print('WP195RESULT', repr(m.%s(*%r)))\n" % (_HERE, _DIST_BIN, fn, args))
    argv, env = pinned_child(["-u", "-c", code], env)
    p = subprocess.run(argv, env=env, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, text=True, encoding="utf-8",
                       errors="replace", timeout=600, cwd=_HERE)
    for line in p.stdout.splitlines():
        if line.startswith("WP195RESULT "):
            return ast.literal_eval(line[len("WP195RESULT "):])
    raise AssertionError("child %s failed (rc %d):\n%s\n%s"
                         %(fn, p.returncode, p.stdout[-3000:], p.stderr[-3000:]))


def _cache_info(eid):
    """[enabled, state, hits, builds, consecutive misses]"""
    return list(ops.eleResponse(eid, "ladrunoTangentCache"))


# ---------------------------------------------------------------------------
#  a stdBrick block (n^3 elements), base fixed, top loaded
# ---------------------------------------------------------------------------
def _block(mat_cmd, n=3, h=1000.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial(*mat_cmd)
    nid = {}
    t = 1
    for i in range(n + 1):
        for j in range(n + 1):
            for k in range(n + 1):
                nid[(i, j, k)] = t
                # a slightly distorted mesh: Jacobians differ per element
                ops.node(t, i * h + 37.0 * j * k, j * h + 11.0 * i, k * h + 5.0 * i * j)
                t += 1
    eid = 1
    for i in range(n):
        for j in range(n):
            for k in range(n):
                c = [nid[(i, j, k)], nid[(i + 1, j, k)], nid[(i + 1, j + 1, k)],
                     nid[(i, j + 1, k)], nid[(i, j, k + 1)], nid[(i + 1, j, k + 1)],
                     nid[(i + 1, j + 1, k + 1)], nid[(i, j + 1, k + 1)]]
                ops.element("stdBrick", eid, *c, 1, 0.0, 0.0, -2.0e-5)
                eid += 1
    for i in range(n + 1):
        for j in range(n + 1):
            ops.fix(nid[(i, j, 0)], 1, 1, 1)
    top = [nid[(i, j, n)] for i in range(n + 1) for j in range(n + 1)]
    return nid, top, eid - 1


def _digest(h, nid):
    ops.reactions()
    for k in sorted(nid.values()):
        h.update(("%d|%s|%s" % (k, repr(ops.nodeDisp(k)),
                                repr(ops.nodeReaction(k)))).encode())
    h.update(repr(list(ops.eleResponse(1, "stresses"))).encode())


def _static(nid, top, px, steps, h, hook=None):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in top:
        ops.load(n, px, 0.3 * px, -0.5 * px)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / steps)
    ops.analysis("Static")
    rcs = []
    for s in range(steps):
        if hook is not None:
            hook(s)
        rc = ops.analyze(1)
        rcs.append(rc)
        _digest(h, nid)
        if rc != 0:
            break
    return rcs


def _transient(nid, steps, h):
    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()
    ops.rayleigh(0.2, 0.0, 0.001, 0.0)
    ops.timeSeries("Path", 2, "-dt", 0.01,
                   "-values", *[0.0, 500.0, -800.0, 1000.0, -500.0, 300.0, 0.0] * 2)
    ops.pattern("UniformExcitation", 2, 1, "-accel", 2)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    rcs = []
    for _ in range(steps):
        rc = ops.analyze(1, 0.01)
        rcs.append(rc)
        _digest(h, nid)
        if rc != 0:
            break
    return rcs


_ELASTIC = ("ElasticIsotropic", 1, 340.0, 0.262, 2.4e-9)
# J2: yields under the push below, so the Gauss-point tangents change
_PLASTIC = ("J2Plasticity", 1, 300.0, 120.0, 0.05, 0.08, 10.0, 1.0)


def _run_elastic(on):
    _set_cache(on)
    nid, top, ne = _block(_ELASTIC)
    h = hashlib.md5()
    rcs = _static(nid, top, 20.0, 4, h)
    rcs += _transient(nid, 12, h)
    info = [_cache_info(e) for e in (1, ne)]
    return h.hexdigest(), rcs, info


def test_elastic_block_cache_equals_uncached_static_and_transient():
    d_off, r_off, i_off = _child("_run_elastic", True, on=False)
    d_on, r_on, i_on = _child("_run_elastic", True, on=True)
    assert all(rc == 0 for rc in r_on), r_on
    assert r_on == r_off
    assert d_on == d_off
    # off really is off; on really hits (one build per element, the rest hits)
    for inf in i_off:
        assert inf[0] == 0.0 and inf[2] == 0.0 and inf[3] == 0.0, inf
    for inf in i_on:
        assert inf[0] == 1.0, inf
        assert inf[1] == 1.0, inf        # valid, uniform D
        assert inf[3] == 1.0, inf        # built once
        assert inf[2] >= 10.0, inf       # served from the cache afterwards


def _run_plastic(on):
    _set_cache(on)
    nid, top, ne = _block(_PLASTIC)
    h = hashlib.md5()
    rcs = _static(nid, top, 3000.0, 8, h)
    info = [_cache_info(e) for e in range(1, ne + 1)]
    return h.hexdigest(), rcs, info


def test_material_state_change_rebuilds_and_matches_uncached():
    d_off, r_off, _ = _child("_run_plastic", True, on=False)
    d_on, r_on, info = _child("_run_plastic", True, on=True)
    assert r_on == r_off
    assert sum(1 for rc in r_on if rc == 0) >= 4, r_on
    assert d_on == d_off
    # the yielding elements had to rebuild: some element built more than once
    assert max(inf[3] for inf in info) > 1.0, info


def _run_param(on, mode):
    _set_cache(on)
    nid, top, ne = _block(_ELASTIC)
    if mode == "updateParameter":
        # registers the MATERIALS of element 5 directly: Brick::updateParameter
        # is never called on this route -- only the raw-bits compare sees it
        ops.parameter(1, "element", 5, "E")

    def hook(s):
        if s == 2:
            if mode == "updateParameter":
                ops.updateParameter(1, 900.0)
            else:
                ops.setParameter("-val", 900.0, "-ele", 5, "E")

    h = hashlib.md5()
    rcs = _static(nid, top, 20.0, 5, h, hook)
    return h.hexdigest(), rcs, _cache_info(5)


@pytest.mark.parametrize("mode", ["updateParameter", "setParameter"])
def test_parameter_update_invalidates_and_matches_uncached(mode):
    d_off, r_off, _ = _child("_run_param", True, mode, on=False)
    d_on, r_on, info = _child("_run_param", True, mode, on=True)
    d_ref, _, _ = _child("_run_elastic_no_param", on=True)
    assert all(rc == 0 for rc in r_on), r_on
    assert r_on == r_off
    assert d_on == d_off
    assert d_on != d_ref            # the update did change the answer
    assert info[3] == 2.0, info     # built, then rebuilt after the update


def _run_elastic_no_param():
    _set_cache(True)
    nid, top, ne = _block(_ELASTIC)
    h = hashlib.md5()
    rcs = _static(nid, top, 20.0, 5, h)
    return h.hexdigest(), rcs, None


def _run_geometry(on):
    _set_cache(on)
    nid, top, ne = _block(_ELASTIC)
    corner = nid[(1, 1, 1)]

    def hook(s):
        if s == 2:
            ops.setNodeCoord(corner, 3, 1000.0 + 120.0)

    h = hashlib.md5()
    rcs = _static(nid, top, 20.0, 5, h, hook)
    return h.hexdigest(), rcs, _cache_info(1)


def test_geometry_change_invalidates_and_matches_uncached():
    d_off, r_off, _ = _child("_run_geometry", True, on=False)
    d_on, r_on, info = _child("_run_geometry", True, on=True)
    assert all(rc == 0 for rc in r_on), r_on
    assert r_on == r_off
    assert d_on == d_off
    assert info[3] == 2.0, info     # element 1 touches the moved node


def _run_reset(on):
    _set_cache(on)
    nid, top, ne = _block(_PLASTIC)
    h1 = hashlib.md5()
    r1 = _static(nid, top, 3000.0, 6, h1)
    ops.reset()                     # Domain::revertToStart
    after = _cache_info(1)
    ops.remove("loadPattern", 1)
    ops.remove("timeSeries", 1)
    ops.setTime(0.0)
    h2 = hashlib.md5()
    r2 = _static(nid, top, 3000.0, 6, h2)
    return h1.hexdigest(), h2.hexdigest(), r1, r2, after


def test_revert_to_start_invalidates_and_replays_identically():
    a1, a2, ra1, ra2, _ = _child("_run_reset", True, on=False)
    b1, b2, rb1, rb2, after = _child("_run_reset", True, on=True)
    assert ra1 == rb1 and ra2 == rb2
    assert a1 == b1 and a2 == b2
    assert b1 == b2                 # the replay after reset is the first run
    assert after[1] == 0.0, after   # the cache was emptied by revertToStart
