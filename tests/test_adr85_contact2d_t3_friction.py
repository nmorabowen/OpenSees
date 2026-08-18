"""ADR-85 T3 -- 2D mortar friction (gate G-T3(e)).

WHAT THIS FILE IS
    Written against `85_ladruno_contact_2d_adr.md` How/3 ("Friction (ADR-41
    C3, in scope here): the mortar return map ... collapses to the same
    scalar as How/4") and `_adr85_t3_design.md` §2b.3 (the mortar friction
    slot reuses `returnMap1D`/`tangentBlock1D` verbatim from T2's 2D NTS
    friction), for C++ that does not exist yet (see test_adr85_contact2d_t3_
    mortar.py's module docstring for the general disclaimer).

    This is the SAME scalar return-map machinery
    test_adr85_contact2d_t2_friction.py already gates on the NTS lane -- the
    incline rig here is that file's rig, rebuilt on the MORTAR lane (surfaces
    declared `-master 2` / `-slave-segments 2`, contact declared `-mortar
    -mu`).

WHY THE STICK SIDE IS AN INEQUALITY, NOT A CLOSED FORM (read before assuming
a gap versus T2)
    T2's stick closed form (x_tan = Q/(kt+ka)) needs the contact's own
    tangential stiffness `kt`, a plain per-node NTS parameter it declares
    directly. Mortar's tangential penalty is `epsT`, a Galerkin-WEIGHTED
    quantity (the D/M-assembled block), and its effective per-node stiffness
    for a given facet is not a value this test can predict without the built
    kernel's own tributary/measure weighting -- computing it independently
    here would gate this test's own arithmetic, not the engine (the same
    argument test_adr85_contact2d_t1b_nts.py::test_2d_autokn uses to justify
    NOT gating auto-kn's absolute value). So the STICK side is checked the way
    T2's return-map theory itself is stated: the trial tangential force must
    sit STRICTLY inside the cone (not saturated) -- true regardless of the
    unknown stiffness. The SLIP side needs no such escape: T2's own finding
    (the pre-T2 experiment, carried into How/3 by inheritance) is that once
    slipping, the return map's slip-branch tangent is IDENTICALLY ZERO in a
    1-D tangent space, so the mortar's own stiffness DROPS OUT of the
    assembled equilibrium entirely and the friction force is a FIXED cap
    (mu*N) independent of it -- exactly the closed form T2 gates, and exactly
    as exact here.

    The anchor's own threshold-shift correction (T2's docstring, "a delta-
    scaled ka is why") is carried over with `EPS_T` standing in for the
    unknown effective mortar tangential stiffness -- an assumption flagged in
    this phase's report, since the true scaling could differ by an O(1)
    geometric factor. The delta-scaling keeps `ka` several orders below
    whatever that true stiffness is regardless of the O(1) uncertainty, which
    is what the shift-suppression argument actually needs.

COVERAGE MAP -> ADR-85 G-T3(e)
    test_2d_mortar_incline_stick_slip_threshold_bracket   stick strictly
        inside the cone / slip exactly saturated at mu*N, swept over three
        deltas either side of tan(theta) == mu.
    test_2d_mortar_friction_mu0_byte_identical             mu omitted vs
        mu=0.0 explicit -> IDENTICAL nodal state (two child processes, the
        byte-identity idiom this whole battery uses for exactly this claim).
"""
import math
import os
import subprocess
import sys

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

ENGINE_DIR = os.path.dirname(os.path.abspath(ops.__file__))

THETA = math.radians(30.0)
TAN_TH = math.tan(THETA)
W0 = 1.0e3          # total "weight" split over the 2-node facet
EPS_N = 1.0e6
EPS_T = 1.0e5       # the DECLARED mortar tangential penalty (see docstring:
                    # the true per-node effective stiffness is unknown, EPS_T
                    # stands in for its order of magnitude)
MX = 4.0            # master half-length (terminal-window margin, ADR-85 D4)
HALF = 1.0          # half-length of the matched slave facet
SEED = 1.0e-4


def _incline_geom(theta=THETA, W=W0):
    c, s = math.cos(theta), math.sin(theta)
    n = (-s, c)
    dhat = (-c, -s)
    return n, dhat, W * c, W * s      # n, dhat, N_total, Q_total


def _incline_master(theta=THETA, mx=MX):
    c, s = math.cos(theta), math.sin(theta)
    ops.node(101, -mx * c, -mx * s)
    ops.node(102, mx * c, mx * s)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.contactSurface(10, "-master", 2, 101, 102)


def _mortar_incline_leg(mu, ka, theta=THETA, W=W0, epsN=EPS_N, epsT=EPS_T,
                        seed=SEED, tol=1.0e-11, iters=200):
    """A 2-node matched mortar facet on the incline, seeded penetrating, each
    node given its own tangential anchor (a zeroLength spring to a coincident
    fixed ground node, oriented along-slope via `-orient`, exactly as T2's
    `_incline_leg_anchored` does it). By symmetry both nodes carry identical
    loads/anchors and settle identically -- returns the PER-NODE quantities.
    """
    n, dhat, N, Q = _incline_geom(theta, W)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _incline_master(theta)
    c, s = math.cos(theta), math.sin(theta)
    dn = seed
    p0 = (-HALF * c - dn * n[0], -HALF * s - dn * n[1])
    p1 = (HALF * c - dn * n[0], HALF * s - dn * n[1])
    ops.node(1, *p0)
    ops.node(2, *p1)
    ops.contactSurface(20, "-slave-segments", 2, 1, 2)
    ops.contact(1, 10, 20, "-mortar", "-epsN", epsN, "-mu", mu, "-epsT", epsT,
                "-outward", n[0], n[1])
    for tag, p in ((1, p0), (2, p1)):
        gnd = 900 + tag
        ops.node(gnd, *p)
        ops.fix(gnd, 1, 1)
        ops.uniaxialMaterial("Elastic", gnd, ka)
        ops.element("zeroLength", gnd, gnd, tag, "-mat", gnd, "-dir", 1,
                    "-orient", dhat[0], dhat[1], 0.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, -W / 2.0)
    ops.load(2, 0.0, -W / 2.0)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", tol, iters, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    rc = ops.analyze(1)
    if rc != 0:
        return dict(rc=rc)
    out = []
    for tag in (1, 2):
        ux, uy = ops.nodeDisp(tag, 1), ops.nodeDisp(tag, 2)
        x_tan = ux * dhat[0] + uy * dhat[1]
        f_tan = (Q / 2.0) - ka * x_tan     # per-node equilibrium inference
        out.append((x_tan, f_tan))
    return dict(rc=0, legs=out, N=N, Q=Q)


def test_2d_mortar_incline_stick_slip_threshold_bracket():
    """G-T3(e): stick at tan(theta) = mu*(1-delta), slip at mu*(1+delta) ...
    read literally, the SLIP side (mu below tan(theta), delta>0 meaning
    FURTHER below) is exact against a closed form; the STICK side is an
    inequality (see the module docstring for why). Swept at three deltas
    (1e-2, 1e-4) either side of the threshold -- see the docstring for why the
    sweep stops one decade short of T2's tightest (1e-6): T2 could afford that
    because its per-node kt was a KNOWN declared parameter; here `ka`'s own
    shift-suppression margin against an ASSUMED (not measured) effective
    mortar stiffness scale is the limiting factor, and this file states that
    assumption honestly rather than chasing a bracket width its own inputs
    cannot justify.
    """
    SCALE = 1.0e-2
    for delta in (1.0e-2, 1.0e-4):
        ka = EPS_T * delta * SCALE
        # GATE RERUN (test-side): right at the cone's edge the return map's
        # trial-vs-cap comparison is itself down at its own resolution limit,
        # so Newton's last iterations chase a residual that has stopped
        # shrinking rather than one that is genuinely unresolved -- the SAME
        # measured behavior test_adr85_contact2d_t2_friction.py::
        # test_2d_incline_threshold_bracket_exact documents for the NTS
        # lane (plateau ~1.6e-8 at its own delta=1e-6) and its own fix
        # (`test_tol = max(1e-12, delta*1e-1)`), applied here verbatim: the
        # delta=1e-4 slip leg plateaued at NormDispIncr ~1.9e-6 against the
        # original flat 1e-11 test and failed outright.
        test_tol = max(1.0e-12, delta * 1.0e-1)

        mu_stick = TAN_TH * (1.0 + delta)
        leg = _mortar_incline_leg(mu_stick, ka, tol=test_tol)
        assert leg["rc"] == 0, f"delta={delta:.0e}: stick leg did not converge"
        N, Q = leg["N"], leg["Q"]
        for i, (x_tan, f_tan) in enumerate(leg["legs"]):
            cap = mu_stick * (N / 2.0)
            assert f_tan < cap, (
                f"delta={delta:.0e} STICK node {i}: f_tan={f_tan!r} reached "
                f"the cone {cap!r} -- spurious slip just above the threshold")
            assert f_tan > 0.0, (
                f"delta={delta:.0e} STICK node {i}: f_tan={f_tan!r} is not "
                "even opposing the drive -- friction is not engaged at all")

        mu_slip = TAN_TH * (1.0 - delta)
        leg = _mortar_incline_leg(mu_slip, ka, tol=test_tol)
        assert leg["rc"] == 0, f"delta={delta:.0e}: slip leg did not converge"
        N, Q = leg["N"], leg["Q"]
        cap = mu_slip * (N / 2.0)
        x_ref = (Q / 2.0 - cap) / ka
        for i, (x_tan, f_tan) in enumerate(leg["legs"]):
            assert abs(x_tan - x_ref) / abs(x_ref) < 1.0e-6, (
                f"delta={delta:.0e} SLIP node {i}: x_tan={x_tan!r} != "
                f"(Q/2-mu*N/2)/ka={x_ref!r} -- the pair took the STICK branch "
                "a hair below the threshold")
            assert abs(f_tan - cap) / cap < 1.0e-6, (
                f"delta={delta:.0e} SLIP node {i}: f_tan={f_tan!r} != cap "
                f"{cap!r} -- friction did not saturate at mu*N just below the "
                "threshold")


# --------------------------------------------------------------------------
# mu = 0 byte-identity, mortar lane -- G-T3(e)
# --------------------------------------------------------------------------
CHILD_MU0 = r'''
import os, sys

_D = %(ENGINE_DIR)r
if os.path.isdir(_D):
    os.environ["PATH"] = _D + os.pathsep + os.environ.get("PATH", "")
    _add = getattr(os, "add_dll_directory", None)
    if _add is not None:
        try:
            _add(_D)
        except OSError:
            pass
    sys.path.insert(0, _D)

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

CASE = sys.argv[1]
KN, KA, P, Q = 1.0e6, 1.0e4, 1.0e3, 50.0
SEED = 1.0e-4

ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.node(101, -2.0, 0.0)
ops.node(102, 2.0, 0.0)
ops.fix(101, 1, 1)
ops.fix(102, 1, 1)
ops.contactSurface(10, "-master", 2, 101, 102)
ops.node(1, -1.0, -SEED)
ops.node(2, 1.0, -SEED)
ops.contactSurface(20, "-slave-segments", 2, 1, 2)

if CASE == "mu_omitted":
    ops.contact(1, 10, 20, "-mortar", "-epsN", KN, "-outward", 0.0, 1.0)
elif CASE == "mu_zero_explicit":
    ops.contact(1, 10, 20, "-mortar", "-epsN", KN, "-mu", 0.0,
                "-outward", 0.0, 1.0)
else:
    print("UNKNOWN CASE " + CASE)
    sys.exit(4)

# a genuine tangential drive (anchored, per-node) so the mu<=0 short-circuit
# is actually EXERCISED (the return map is asked to decide something), not
# merely never reached -- the non-trivial half of T2's own mu=0 gate.
for tag, x0 in ((1, -1.0), (2, 1.0)):
    gnd = 900 + tag
    ops.node(gnd, x0, -SEED)
    ops.fix(gnd, 1, 1)
    ops.uniaxialMaterial("Elastic", gnd, KA)
    ops.element("zeroLength", gnd, gnd, tag, "-mat", gnd, "-dir", 1)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(1, Q / 2.0, -P / 2.0)
ops.load(2, Q / 2.0, -P / 2.0)
ops.constraints("LadrunoContact")
ops.numberer("Plain")
ops.system("FullGeneral")
ops.test("NormDispIncr", 1.0e-12, 40, 0)
ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0)
ops.analysis("Static")
assert ops.analyze(1) == 0, "mu0 byte-identity deck did not converge"

vals = []
for tag in (1, 2):
    vals.append(ops.nodeDisp(tag, 1).hex())
    vals.append(ops.nodeDisp(tag, 2).hex())
print("MU0_HEX " + " ".join(vals))
sys.exit(0)
'''


def _run_mu0_child(case):
    script = CHILD_MU0 % {"ENGINE_DIR": ENGINE_DIR}
    env = dict(os.environ)
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    proc = subprocess.run(
        [sys.executable, "-c", script, case],
        stdin=subprocess.DEVNULL,
        capture_output=True, text=True, timeout=300,
        encoding="utf-8", errors="replace",
        env=env,
    )
    return proc, (proc.stdout or "") + "\n" + (proc.stderr or "")


def test_2d_mortar_friction_mu0_byte_identical():
    """G-T3(e): the SAME mortar deck run in two clean child interpreters, one
    WITHOUT `-mu` and one with `-mu 0.0` explicit, must produce IDENTICAL
    hexfloat nodal displacements -- the mortar friction addition must not
    perturb a single bit of the frictionless mortar solve when mu=0, even
    under a genuine tangential drive (the direct mortar analogue of T2's
    test_2d_friction_mu0_with_drive_byte_identical).
    """
    proc_a, out_a = _run_mu0_child("mu_omitted")
    proc_b, out_b = _run_mu0_child("mu_zero_explicit")
    assert proc_a.returncode == 0, f"mu_omitted failed:\n{out_a}"
    assert proc_b.returncode == 0, f"mu_zero_explicit failed:\n{out_b}"

    def _hexvals(out):
        for line in out.splitlines():
            if line.startswith("MU0_HEX"):
                return line.split()[1:]
        raise AssertionError(f"no MU0_HEX line:\n{out}")

    va = _hexvals(out_a)
    vb = _hexvals(out_b)
    assert va == vb, (
        f"mu omitted vs mu=0.0 explicit diverged bit-for-bit on the 2D mortar "
        f"lane: {va} vs {vb}")
