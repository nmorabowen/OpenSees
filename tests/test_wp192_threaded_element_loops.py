"""WP-192 -- threaded element loops for the Paper 1 element set.

WP-107 threaded ONE loop (`Domain::update`) for ONE element (LadrunoQuad) and
ran the whole loop serial as soon as any un-audited element was present.
WP-192 adds:

  * the integrators' element TANGENT and RESIDUAL loops (ADR-75b loops B/C),
    as a chunked compute + ORDERED serial scatter into the SOE, so the
    assembled system is bit-identical to serial at every thread count;
  * the element COMMIT loop of `Domain::commit`;
  * a PER-TYPE opt-in: audited elements run on threads, every other element
    runs alone in its serial place (`Element::ladrunoThreadSafeUpdate`,
    `::ladrunoThreadSafeForm`, `::ladrunoThreadSafeCommit`);
  * thread-safe ASDShellQ4 (+ LayeredShell / PlateFromPlaneStress / PlateFiber /
    PlaneStress wrappers, PlateRebar, ASDConcrete3D, ElasticMembranePlate),
    ForceBeamColumn3d (+ FiberSection3d, LinearCrdTransf3d, ASDConcrete1D,
    Hysteretic, Elastic), stdBrick (+ ElasticIsotropic 3D), ASDEmbeddedNode,
    ASDAbsorbingBoundary3D (update/forms; commit stays serial);
  * per-thread copies of the class-wide pools on those paths (FE_Element,
    TransformationFE, Element's Rayleigh pool, Matrix's solve work area).

THE ORACLE IS THE FULL FIELD AT EVERY STEP (WP-107 red-team S6): every node's
displacement and reaction plus element responses, at repr() precision, md5'd
per step, and the digests must be EQUAL across 1 / 2 / 4 threads -- equality,
not a tolerance, because both the update loop (no reduction) and the scatter
(serial, in serial order) are bit-identical by construction.

WHICH TEST CATCHES WHICH MUTATION:
  * revert ASDShellQ4Globals::instance() to a plain `static` -> the shell
    digests split across thread counts (test_shell_frame_*);
  * revert FE_Element::ladrunoTan()/ladrunoRes() to the shared pool -> every
    multi-thread digest differs (every FE builds into one buffer);
  * revert Brick::stiff/resid to plain static -> test_brick_block_* splits;
  * drop the serial re-insertion of un-audited elements -> test_mixed_* splits
    or crashes.

See `Ladruno_implementation/192_ladruno_threaded_element_loops_adr.md`.
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


# ---------------------------------------------------------------------------
#  materials of the Paper 1 decks (rev3 material 5 / 9 / Hysteretic 8)
# ---------------------------------------------------------------------------
_TE = [0.0, 0.00011553, 0.00019255, 0.00439755, 0.0216028, 0.216028]
_TS = [0.0, 3.276, 3.64, 0.728, 0.00364, 0.00364]
_TD = [0.0, 0.0, 0.259259, 0.992231, 0.999992, 0.999999]
_CE = [0.0, 0.000641834, 0.000903853, 0.00116587, 0.00142789, 0.00168991,
       0.00195193, 0.00221394, 0.00247596, 0.00273798, 0.003, 0.28206, 0.28262]
_CS = [0.0, 18.2, 24.0932, 28.0248, 30.8124, 32.8264, 34.2723, 35.278, 35.9297,
       36.2886, 36.4, 7.28, 7.28]
_CD = [0.0, 0.0, 0.0, 0.0, 0.0283541, 0.109709, 0.184949, 0.252924, 0.314092,
       0.369252, 0.419201, 0.998847, 0.998849]
_HYST = [491.0, 0.002455, 492.18599, 0.02, 652.29469, 0.1,
         -491.0, -0.002455, -492.18599, -0.02, -652.29469, -0.1, 0.5, 0.8, 0.0, 0.0]


def _materials(implex):
    flags = ["-implex"] if implex else []
    ops.nDMaterial("ASDConcrete3D", 1, 28356.2, 0.2, "-Te", *_TE, "-Ts", *_TS,
                   "-Td", *_TD, "-Ce", *_CE, "-Cs", *_CS, "-Cd", *_CD,
                   "-rho", 2.4e-9, "-Kc", 2.0 / 3.0, *flags,
                   "-crackPlanes", 4, 4, 45.0, "-autoRegularization", 6.37755)
    ops.uniaxialMaterial("Hysteretic", 2, *_HYST)
    ops.nDMaterial("PlateFromPlaneStress", 3, 1, 11910.0)
    ops.nDMaterial("PlateRebar", 4, 2, 0.0)
    ops.nDMaterial("PlateRebar", 5, 2, 90.0)
    ops.uniaxialMaterial("ASDConcrete1D", 6, 28356.2, "-Te", *_TE, "-Ts", *_TS,
                         "-Td", *_TD, "-Ce", *_CE, "-Cs", *_CS, "-Cd", *_CD,
                         *flags, "-autoRegularization", 6.37755)
    # LayeredShell: cover, rebar x/y, core (PlateFromPlaneStress) and one layer
    # of the raw 3D material (-> PlateFiber wrapper), as section 19 of rev3 1D.
    ops.section("LayeredShell", 10, 7, 3, 25.0, 4, 1.5, 5, 1.5, 3, 50.0, 1, 50.0,
                5, 1.5, 4, 1.5)


# ---------------------------------------------------------------------------
#  model 1: RC wall of layered ASDShellQ4 + two fiber forceBeamColumns +
#           a rigid diaphragm on top (TransformationFE path)
# ---------------------------------------------------------------------------
def _shell_frame(nx=6, nz=6, implex=False, unaudited=False):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    _materials(implex)
    L, H = 3000.0, 3000.0
    nid = {}
    t = 1
    for i in range(nx + 1):
        for k in range(nz + 1):
            nid[(i, k)] = t
            ops.node(t, i * L / nx, 0.0, k * H / nz)
            ops.mass(t, 0.5, 0.5, 0.5, 0.0, 0.0, 0.0)
            t += 1
    eid = 1
    for i in range(nx):
        for k in range(nz):
            ops.element("ASDShellQ4", eid, nid[(i, k)], nid[(i + 1, k)],
                        nid[(i + 1, k + 1)], nid[(i, k + 1)], 10,
                        "-noeas", "-drillingNL", "-local", 1.0, 0.0, 0.0)
            eid += 1
    # two columns beside the wall, fiber section of ASDConcrete1D + Hysteretic
    ops.section("Fiber", 20, "-GJ", 1.0e12)
    for y in (-150.0, -50.0, 50.0, 150.0):
        for z in (-150.0, -50.0, 50.0, 150.0):
            ops.fiber(y, z, 10000.0, 6)
    for y, z in ((-120.0, -120.0), (120.0, -120.0), (-120.0, 120.0), (120.0, 120.0)):
        ops.fiber(y, z, 380.0, 2)
    ops.geomTransf("Linear", 1, 1.0, 0.0, 0.0)
    ops.beamIntegration("Lobatto", 1, 20, 5)
    cols = []
    for x in (-1000.0, L + 1000.0):
        b, tp = t, t + 1
        ops.node(b, x, 0.0, 0.0)
        ops.node(tp, x, 0.0, H)
        ops.mass(tp, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0)
        ops.fix(b, 1, 1, 1, 1, 1, 1)
        ops.element("forceBeamColumn", eid, b, tp, 1, 1)
        cols.append(tp)
        eid += 1
        t += 2
    for i in range(nx + 1):
        ops.fix(nid[(i, 0)], 1, 1, 1, 1, 1, 1)
    master = t
    ops.node(master, L / 2.0, 0.0, H)
    ops.fix(master, 0, 0, 1, 1, 1, 0)
    slaves = [nid[(i, nz)] for i in range(nx + 1)] + cols
    ops.rigidDiaphragm(3, master, *slaves)
    if unaudited:
        ops.uniaxialMaterial("Elastic", 99, 1.0e3)
        ops.node(990001, L / 2.0, 0.0, -10.0)
        ops.fix(990001, 1, 1, 1, 1, 1, 1)
        # un-audited element: a Truss (zeroLength until WP-195 opted it in)
        ops.element("Truss", 990000, 990001, nid[(nx // 2, 0)], 1.0, 99)
    return nid, master


def _gravity_push_and_shake(nid, master, steps_push, steps_dyn, tag):
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(master, 100000.0, 15000.0, 0.0, 0.0, 0.0, 0.0)
    for (i, k), n in nid.items():
        ops.load(n, 0.0, 0.0, -2000.0, 0.0, 0.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-6, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / steps_push)
    ops.analysis("Static")
    h = hashlib.md5()
    rcs = []

    def snap():
        ops.reactions()
        for n in sorted(nid.values()) + [master]:
            h.update(("%d|%s|%s" % (n, repr(ops.nodeDisp(n)),
                                    repr(ops.nodeReaction(n)))).encode())
        h.update(repr(list(ops.eleResponse(1, "stresses"))).encode())

    for _ in range(steps_push):
        rc = ops.analyze(1)
        rcs.append(rc)
        snap()
        if rc != 0:
            return h.hexdigest(), rcs
    ops.loadConst("-time", 0.0)
    ops.wipeAnalysis()
    ops.rayleigh(0.1, 0.0, 0.002, 0.0)
    ops.timeSeries("Path", 2, "-dt", 0.01,
                   "-values", *[0.0, 500.0, -800.0, 1000.0, -500.0, 300.0, 0.0] * 3)
    ops.pattern("UniformExcitation", 2, 1, "-accel", 2)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-6, 100, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    for _ in range(steps_dyn):
        rc = ops.analyze(1, 0.01)
        rcs.append(rc)
        snap()
        if rc != 0:
            break
    return h.hexdigest(), rcs


def _run_shell(n, implex=False, unaudited=False):
    nid, master = _shell_frame(implex=implex, unaudited=unaudited)
    ops.ladrunoThreads(n)
    try:
        return _gravity_push_and_shake(nid, master, 6, 12, "shell")
    finally:
        ops.ladrunoThreads(1)


@pytest.mark.parametrize("implex", [False, True])
def test_shell_frame_bit_identity_across_thread_counts(implex, capfd):
    counts = _threads_available()
    if len(counts) < 2:
        pytest.skip("single hardware thread")
    digests, rcss = {}, {}
    for n in counts:
        d, rcs = _run_shell(n, implex=implex)
        # Equality is the claim, convergence is not: a step that fails must fail
        # identically at every thread count. The static push (6 steps) must
        # converge, or the deck tests nothing.
        assert sum(1 for rc in rcs if rc == 0) >= 6, (n, rcs)
        digests[n], rcss[n] = d, rcs
    assert len(set(map(tuple, rcss.values()))) == 1, rcss
    text = capfd.readouterr().err + capfd.readouterr().out
    assert len(set(digests.values())) == 1, digests
    # the threaded loops must actually have run (fail-loud in both directions)
    assert "element update loop THREADED" in text, text
    assert "element tangent loop threaded" in text, text
    assert "element residual loop threaded" in text, text


def test_mixed_unaudited_element_runs_serial_in_place(capfd):
    """A Truss (not on the allowlist) in the middle of the domain: the
    others thread, it runs alone in its serial place, and the field is still
    bit-identical to 1 thread."""
    counts = _threads_available()
    if len(counts) < 2:
        pytest.skip("single hardware thread")
    d1, r1 = _run_shell(1, unaudited=True)
    dn, rn = _run_shell(counts[-1], unaudited=True)
    text = capfd.readouterr().err
    assert sum(1 for rc in r1 if rc == 0) >= 6, r1
    assert r1 == rn
    assert d1 == dn
    assert "PARTLY threaded" in text, text
    assert "element 990000 " in text, text
    assert "run serially in place" in text, text


def test_repeated_runs_hash_identical_race_probe():
    """The race probe without ThreadSanitizer: the same deck, the most threads
    this box allows, five times. A schedule-dependent race shows up as a digest
    that differs between repetitions even when the 1/2/4 comparison happened to
    pass. (TSan is not available on the Windows/MSVC build; see the ADR.)"""
    counts = _threads_available()
    if len(counts) < 2:
        pytest.skip("single hardware thread")
    n = counts[-1]
    ref, rref = _run_shell(1)
    for _ in range(5):
        d, rcs = _run_shell(n)
        assert rcs == rref
        assert d == ref


# ---------------------------------------------------------------------------
#  model 2: stdBrick soil block with ASDEmbeddedNode ties to a beam
# ---------------------------------------------------------------------------
def _brick_block(n=5):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 340.0, 0.262, 2.4e-9)
    h = 1000.0
    nid = {}
    t = 1
    for i in range(n + 1):
        for j in range(n + 1):
            for k in range(n + 1):
                nid[(i, j, k)] = t
                ops.node(t, i * h, j * h, k * h)
                t += 1
    eid = 1
    for i in range(n):
        for j in range(n):
            for k in range(n):
                c = [nid[(i, j, k)], nid[(i + 1, j, k)], nid[(i + 1, j + 1, k)],
                     nid[(i, j + 1, k)], nid[(i, j, k + 1)], nid[(i + 1, j, k + 1)],
                     nid[(i + 1, j + 1, k + 1)], nid[(i, j + 1, k + 1)]]
                ops.element("stdBrick", eid, *c, 1)
                eid += 1
    for i in range(n + 1):
        for j in range(n + 1):
            ops.fix(nid[(i, j, 0)], 1, 1, 1)
    # embedded nodes (free 3-dof nodes tied inside bricks) carrying point loads
    emb = []
    for m in range(4):
        en = 900000 + m
        x = (0.35 + 0.1 * m) * h * n
        ops.node(en, x, 0.5 * h * n, 0.75 * h * n)
        i = int(x // h)
        ops.element("ASDEmbeddedNodeElement", 800000 + m, en,
                    nid[(i, 2, 3)], nid[(i + 1, 2, 3)], nid[(i, 3, 3)], nid[(i, 2, 4)],
                    "-K", 1.0e6)
        emb.append(en)
    return nid, emb


def _run_brick(nthreads):
    nid, emb = _brick_block()
    ops.ladrunoThreads(nthreads)
    try:
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for en in emb:
            ops.load(en, 1.0e4, -2.0e4, 3.0e3)
        ops.constraints("Transformation")
        ops.numberer("RCM")
        ops.system("UmfPack")
        ops.test("NormDispIncr", 1.0e-8, 20, 0)
        ops.algorithm("Newton")
        ops.rayleigh(0.2, 0.0, 0.001, 0.0)
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        hh = hashlib.md5()
        rcs = []
        for _ in range(10):
            rc = ops.analyze(1, 0.01)
            rcs.append(rc)
            ops.reactions()
            for k in sorted(nid.values()) + emb:
                hh.update(("%d|%s|%s" % (k, repr(ops.nodeDisp(k)),
                                         repr(ops.nodeReaction(k)))).encode())
            if rc != 0:
                break
        return hh.hexdigest(), rcs
    finally:
        ops.ladrunoThreads(1)


def test_brick_block_bit_identity_across_thread_counts(capfd):
    counts = _threads_available()
    if len(counts) < 2:
        pytest.skip("single hardware thread")
    digests = {}
    for n in counts:
        d, rcs = _run_brick(n)
        assert all(rc == 0 for rc in rcs), (n, rcs)
        digests[n] = d
    text = capfd.readouterr().err
    assert len(set(digests.values())) == 1, digests
    assert "element tangent loop threaded" in text, text


def test_hybrid_flag_is_accepted_and_one_thread_stays_default():
    assert ops.ladrunoThreads() == 1
    assert ops.ladrunoThreads(1, "-hybrid") == 1
    assert ops.ladrunoThreads() == 1
