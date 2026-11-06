"""WP-196 -- the commit phase's EnergyBalance sweep: same numbers, less work.

What changed in the EnergyBalance recorder (``SRC/recorder/EnergyBalanceKernel.h``,
``EnergyBalanceRecorder.cpp``):

  * the per-element terms are formed by ``ebkernel::elementEnergyTerms`` into one
    slot per element and summed afterwards in the live element order -- so the
    element sweep can run on the WP-192 threads and the sums still see the same
    operands in the same order;
  * an element that declares ``ladrunoEnergyRayleighOnly()`` (ASDShellQ4 without a
    Damping object) skips the WP-184 "damping outside F" check, which for such an
    element always answers "no" (C v == F_inc - F - M a to round-off) and costs a
    second residual and a second tangent per element per step;
  * the element mass is copied into a per-thread scratch instead of a fresh heap
    Matrix.

The shells are `-noeas -drillingNL`, as in the Paper 1 decks.

THE ORACLE IS THE RECORDER FILE, BYTE FOR BYTE: ``-v2`` with ``-region`` blocks and
``-precision 17``. Three comparisons:

  1. ``LADRUNO_EB_LEGACY=1`` (the WP-184 check on every element, the pre-WP-196
     work) vs the default -- equal bytes, on a Rayleigh shell + fiber-beam model
     AND on one with a zeroLength Viscous dashpot (the element the check exists
     for, which must still take the check path and book v.g);
  2. 1 thread vs 2 / 4 threads (the threaded sweep) -- equal bytes, with the
     announcement asserted so a sweep that silently stayed serial fails;
  3. the dashpot model threaded: the unaudited zeroLength runs alone in its
     serial place and its DW is unchanged.

Mutations this catches: making ``elementEnergyTerms`` skip the check for every
element (the dashpot's DW changes: comparison 1, dashpot arm); reducing the slots
in another order or summing per thread (comparison 2 at 4 threads); dropping the
``thread_local`` on the kernel's scratch (comparison 2, race).

See ``Ladruno_implementation/196_ladruno_commit_phase_cost_adr.md``.
"""
import hashlib
import os
import sys
import tempfile

import pytest

_DIST_BIN = os.environ.get("LADRUNO_DIST_BIN") or os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "bin")
if os.path.isdir(_DIST_BIN) and _DIST_BIN not in sys.path:
    sys.path.insert(0, _DIST_BIN)

from _testbed.subprocess_run import run_python_script   # noqa: E402

pytestmark = [pytest.mark.zone_a]

# The child builds the model, runs it, and prints the md5 of the EnergyBalance
# file plus the recorder's announcement lines. argv: threads dashpot outfile
_CHILD = r'''
import hashlib, sys
sys.path.insert(0, %(dist)r)
try:
    import opensees as ops
except ImportError:
    import openseespy.opensees as ops
nthreads, dashpot, out = int(sys.argv[1]), sys.argv[2] == "1", sys.argv[3]
eas = sys.argv[4] == "1"
shellflags = [] if eas else ["-noeas", "-drillingNL"]
hw = ops.ladrunoThreads(1024); ops.ladrunoThreads(1)
if nthreads > hw:
    print("SKIP_HW", hw); sys.exit(0)

ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 6)
ops.nDMaterial("ElasticIsotropic", 1, 30000.0, 0.2, 2.4e-9)
ops.nDMaterial("PlateFiber", 2, 1)
ops.section("LayeredShell", 10, 3, 2, 60.0, 2, 80.0, 2, 60.0)
nx, nz, L, H = 6, 6, 3000.0, 3000.0
nid = {}; t = 1
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
                    nid[(i + 1, k + 1)], nid[(i, k + 1)], 10, *shellflags, "-local", 1.0, 0.0, 0.0)
        eid += 1
ops.uniaxialMaterial("Elastic", 5, 30000.0)
ops.section("Fiber", 20, "-GJ", 1.0e12)
for y in (-150.0, -50.0, 50.0, 150.0):
    for z in (-150.0, -50.0, 50.0, 150.0):
        ops.fiber(y, z, 10000.0, 5)
ops.geomTransf("Linear", 1, 1.0, 0.0, 0.0)
ops.beamIntegration("Lobatto", 1, 20, 5)
for x in (-1000.0, L + 1000.0):
    ops.node(t, x, 0.0, 0.0); ops.node(t + 1, x, 0.0, H)
    ops.mass(t + 1, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0)
    ops.fix(t, 1, 1, 1, 1, 1, 1)
    ops.element("forceBeamColumn", eid, t, t + 1, 1, 1)
    eid += 1; t += 2
for i in range(nx + 1):
    ops.fix(nid[(i, 0)], 1, 1, 1, 1, 1, 1)
if dashpot:
    # the element the WP-184 check exists for: its viscous force is inside F
    ops.uniaxialMaterial("Viscous", 98, 50.0, 1.0)
    ops.node(990001, L / 2.0, 0.0, H + 10.0)
    ops.fix(990001, 1, 1, 1, 1, 1, 1)
    ops.element("zeroLength", 990000, nid[(nx // 2, nz)], 990001, "-mat", 98, "-dir", 1)
ops.region(1, "-eleRange", 1, nx * nz)
ops.region(2, "-eleRange", nx * nz + 1, eid - 1)
ops.recorder("EnergyBalance", "-v2", "-region", 1, "-region", 2, "-time",
             "-precision", 17, "-file", out)
ops.rayleigh(0.3, 0.0, 0.002, 0.0)
ops.timeSeries("Path", 2, "-dt", 0.01,
               "-values", *[0.0, 500.0, -800.0, 1000.0, -500.0, 300.0, 0.0] * 3)
ops.pattern("UniformExcitation", 2, 1, "-accel", 2)
ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("UmfPack")
ops.test("NormDispIncr", 1.0e-8, 50, 0)
ops.algorithm("Newton")
ops.integrator("Newmark", 0.5, 0.25)
ops.analysis("Transient")
ops.ladrunoThreads(nthreads)
rc = sum(ops.analyze(1, 0.01) != 0 for _ in range(15))
ops.ladrunoThreads(1)
ops.wipe()
print("RC", rc)
print("MD5", hashlib.md5(open(out, "rb").read()).hexdigest())
'''


def _run(tmp, threads, dashpot, legacy, eas=False):
    out = os.path.join(tmp, "eb_%d_%d_%d_%d.txt" % (threads, dashpot, legacy, eas))
    old = os.environ.get("LADRUNO_EB_LEGACY")
    if legacy:
        os.environ["LADRUNO_EB_LEGACY"] = "1"
    else:
        os.environ.pop("LADRUNO_EB_LEGACY", None)
    try:
        code, text = run_python_script(_CHILD % {"dist": _DIST_BIN},
                                       argv=(str(threads), "1" if dashpot else "0", out, "1" if eas else "0"),
                                       merge_stderr=True)
    finally:
        if old is None:
            os.environ.pop("LADRUNO_EB_LEGACY", None)
        else:
            os.environ["LADRUNO_EB_LEGACY"] = old
    assert code == 0, text
    if "SKIP_HW" in text:
        return None, text, out
    assert "RC 0" in text, text
    md5 = [ln.split()[1] for ln in text.splitlines() if ln.startswith("MD5 ")]
    assert md5, text
    return md5[0], text, out


@pytest.mark.parametrize("dashpot,eas", [(False, False), (True, False), (False, True)])
def test_legacy_check_and_default_write_identical_bytes(dashpot, eas):
    """eas=True: the AGQI shell does NOT opt in (its condensation state keeps the
    kind of the last evaluation, ASDShellQ4::ladrunoEnergyRayleighOnly), so it
    takes the check like before and the bytes must still match."""
    with tempfile.TemporaryDirectory() as tmp:
        a, _, fa = _run(tmp, 1, dashpot, legacy=True, eas=eas)
        b, _, fb = _run(tmp, 1, dashpot, legacy=False, eas=eas)
        assert a == b, (open(fa).read()[-400:], open(fb).read()[-400:])
        # the file is not trivially empty: the model block's DW_ele is nonzero
        last = open(fb).read().strip().splitlines()[-1].split()
        assert any(float(x) != 0.0 for x in last[1:]), last


@pytest.mark.parametrize("dashpot", [False, True])
def test_threaded_sweep_writes_identical_bytes(dashpot):
    with tempfile.TemporaryDirectory() as tmp:
        ref, _, _ = _run(tmp, 1, dashpot, legacy=False)
        ran = 0
        for n in (2, 4):
            d, text, _ = _run(tmp, n, dashpot, legacy=False)
            if d is None:
                continue
            ran += 1
            assert d == ref, (n, text[-2000:])
            assert "EnergyBalanceRecorder: element sweep threaded on %d threads" % n in text, text
        if ran == 0:
            pytest.skip("fewer than 2 hardware threads")


def test_dashpot_check_still_runs_under_threads():
    """The dashpot zeroLength keeps the WP-184 check under threads; the legacy
    single-thread file (which takes the check on every element) is the oracle.
    Since WP-195 (#983) zeroLength is thread-safe and opted in, so all 39
    elements (36 shells, 2 beams, the zeroLength) are formed in parallel."""
    with tempfile.TemporaryDirectory() as tmp:
        ref, _, _ = _run(tmp, 1, True, legacy=True)
        d, text, _ = _run(tmp, 4, True, legacy=False)
        if d is None:
            pytest.skip("fewer than 4 hardware threads")
        assert d == ref
        # 36 shells + 2 beams + the zeroLength (WP-195) of 39 elements are formed in parallel
        assert "39 of 39 elements formed in parallel" in text, text
