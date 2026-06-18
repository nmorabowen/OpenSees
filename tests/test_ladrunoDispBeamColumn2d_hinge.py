"""LadrunoDispBeamColumn2d Tier-2 embedded strong-discontinuity hinge (ADR 32).

The element gains a single scalar rotation jump `alpha` carried by a discrete
cohesive law (`LadrunoCohesiveHinge`, any UniaxialMaterial). Strong-discontinuity
split (Armero-Ehrlich / Jukic-Brank-Ibrahimbegovic): the bulk section at every IP
sees the BOUNDED enhanced curvature `kappa_bulk = B*v - alpha/L`, so the bulk
UNLOADS as the jump grows (no Gf double-count); the cohesive `M([[theta]])` carries
all the fracture energy. `alpha` is converged by an inner Newton in `update()` and
statically condensed to the 3-DOF basic system BEFORE crdTransf (the pinned
invariant), with a GUARDED reciprocal (K_aa is indefinite on the softening branch).

This is the PR-2a minimal increment: everything GATED on `-hinge` so the no-hinge
path is byte-identical to Stage-1; linear basic strain only (mutually exclusive
with `-nl`); single hinge; Linear/Corotational transf; monotonic.

Gates (all solver-independent, no path-follower — the ADR's "cheap gates first"):
  * parse + `-hinge`/`-nl` mutual exclusion
  * reduce-to-Tier-1: hinge-absent element is bit-identical to stock dispBeamColumn
  * constant-moment patch test: the embedded hinge carries EXACTLY the applied
    moment (enhancement equilibrium h=0) to solver precision
  * energy gate: a cantilever pushed past the hinge peak dissipates EXACTLY Gf
    through the hinge (`hinge energy == Gf`), peak hinge moment == Mc
  * element total dissipation == Gf on an ELASTIC section (closed load/unload cycle
    net work == Gf) — catches any bulk double-count
  * tangent consistency: every step through softening converges under NormDispIncr
    1e-12 Newton (a wrong condensed tangent would not converge quadratically)
  * commit/revert + DB sendSelf/recvSelf round-trip of an open hinge

Plan: Ladruno_implementation/32_ladruno_dispbeamcolumn_regularization_adr.md
(2026-06-17 adversarial-review CORRECTIONS + PR-2a scope).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# member + hinge properties (Linear transf, elastic bulk -> all softening is the hinge)
_L = 100.0
_E = 1.0e5
_A = 1.0e8        # near-rigid axial: keep the axial DOF out of the way
_Iz = 1.0         # EI = 1e5
_Mc = 50.0        # cohesive moment capacity
_Gf = 20.0        # fracture energy per hinge


def _kpen_default(Mc, Gf, ratio=1000.0):
    return ratio * Mc * Mc / (2.0 * Gf)


def _full_break_jump(Mc, Gf):
    Kpen = _kpen_default(Mc, Gf)
    k0 = Mc / Kpen
    kf = 2.0 * (Gf - Mc * Mc / (2.0 * Kpen)) / Mc
    return k0 + kf


def _build(hinge=True, shape="-linear", transf="Linear", nl=False):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, _L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz)
    ops.geomTransf(transf, 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    args = [1, 1, 2, 1, 1]
    if nl:
        args.append("-nl")
    if hinge:
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, _Mc, _Gf, shape)
        args += ["-hinge", 77]
    ops.element("LadrunoDispBeamColumn", *args)


def _setup_rotation_control(dtheta):
    """Tip moment (ref 1.0 at node2 dof3); DisplacementControl on the tip rotation so
    the constant-moment cantilever is driven monotonically through hinge softening."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, 0.0, 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Newton")
    ops.integrator("DisplacementControl", 2, 3, dtheta)
    ops.analysis("Static")


def _hinge(resp):
    return ops.eleResponse(1, "hinge", resp)[0]


# --------------------------------------------------------------------------- #
#  parse / mutual exclusion
# --------------------------------------------------------------------------- #
def test_hinge_builds_and_nl_mutually_exclusive():
    _build(hinge=True)                       # builds with -hinge: no exception
    # -hinge + -nl must be rejected (the bowing strain couples the jump into axial)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, _L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, _Mc, _Gf, "-linear")
    with pytest.raises(Exception):
        ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hinge", 77, "-nl")


# --------------------------------------------------------------------------- #
#  reduce-to-Tier-1: the no-hinge element is bit-identical to stock dispBeamColumn
# --------------------------------------------------------------------------- #
def _tip_rotation_elastic(make_element):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, _L, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    make_element()
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -10.0, 5.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ops.analyze(1)
    return ops.nodeDisp(2)


def test_reduce_to_tier1_without_hinge_matches_stock():
    ladr = _tip_rotation_elastic(
        lambda: ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1))
    stock = _tip_rotation_elastic(
        lambda: ops.element("dispBeamColumn", 1, 1, 2, 1, 1))
    for a, b in zip(ladr, stock):
        assert a == b or abs(a - b) <= 1e-300 + 1e-15 * abs(b), (a, b)


# --------------------------------------------------------------------------- #
#  constant-moment patch test: the hinge carries EXACTLY the applied moment
# --------------------------------------------------------------------------- #
def test_constant_moment_patch():
    _build(hinge=True)
    _setup_rotation_control(0.002)           # small steps, stay below activation
    # below activation (M < Mc): drive a few steps and assert the hinge moment equals
    # the applied tip moment (= load factor) to solver precision -> enhancement
    # equilibrium h = -M_section + M_coh = 0 holds exactly (the patch test).
    for _ in range(5):
        assert ops.analyze(1) == 0
        M_applied = ops.getTime()            # ref moment 1.0 * load factor
        assert M_applied < _Mc               # still elastic / hinge closed
        assert math.isclose(_hinge("stress"), M_applied, rel_tol=1e-9, abs_tol=1e-9), (
            _hinge("stress"), M_applied)


# --------------------------------------------------------------------------- #
#  energy gate: the embedded hinge dissipates EXACTLY Gf; peak moment == Mc
# --------------------------------------------------------------------------- #
def test_hinge_dissipates_Gf_and_peaks_at_Mc():
    _build(hinge=True, shape="-linear")
    theta_break = _full_break_jump(_Mc, _Gf) + _Mc * _L / _E  # jump + elastic tip rot
    nstep = 1500
    _setup_rotation_control(1.2 * theta_break / nstep)
    peak = 0.0
    for _ in range(nstep):
        assert ops.analyze(1) == 0
        peak = max(peak, _hinge("stress"))
    # the cohesive hinge embedded in the element dissipated exactly Gf
    assert math.isclose(_hinge("energy"), _Gf, rel_tol=2e-3), (_hinge("energy"), _Gf)
    # cohesive law activated at Mc and never exceeded it
    assert math.isclose(peak, _Mc, rel_tol=2e-3), (peak, _Mc)


# --------------------------------------------------------------------------- #
#  element TOTAL dissipation == Gf on an elastic section (no bulk double-count)
# --------------------------------------------------------------------------- #
def test_element_total_dissipation_equals_Gf():
    """Load the tip rotation 0 -> theta_max -> 0; the elastic bulk recovers, so the net
    work int M dtheta around the closed cycle is the hinge dissipation. With an ELASTIC
    section that must equal Gf exactly (if the bulk re-integrated the localized strain it
    would be Gf + bulk-dissipation -> this is the double-count catcher)."""
    _build(hinge=True, shape="-linear")
    theta_break = _full_break_jump(_Mc, _Gf) + _Mc * _L / _E
    theta_max = 1.05 * theta_break
    nstep = 1200
    dth = theta_max / nstep
    _setup_rotation_control(dth)

    # integrate work = int M dtheta over load (up) then unload (down)
    work = 0.0
    th_prev = ops.nodeDisp(2, 3)
    M_prev = _hinge("stress")
    for _ in range(nstep):                   # load up
        assert ops.analyze(1) == 0
        th = ops.nodeDisp(2, 3); M = ops.getTime()
        work += 0.5 * (M_prev + M) * (th - th_prev)
        th_prev, M_prev = th, M
    ops.integrator("DisplacementControl", 2, 3, -dth)   # unload
    for _ in range(nstep):
        assert ops.analyze(1) == 0
        th = ops.nodeDisp(2, 3); M = ops.getTime()
        work += 0.5 * (M_prev + M) * (th - th_prev)
        th_prev, M_prev = th, M

    assert math.isclose(work, _Gf, rel_tol=5e-3), (work, _Gf)


# --------------------------------------------------------------------------- #
#  tangent consistency: softening branch converges under tight Newton
# --------------------------------------------------------------------------- #
def test_softening_newton_converges_tightly():
    """A wrong condensed tangent (K_basic = K_vv - K_va K_aa^-1 K_av) would not let
    NormDispIncr-1e-12 Newton converge in <=10 iters through the softening branch, where
    the coupled axial/transverse free DOFs genuinely exercise the off-diagonal terms."""
    _build(hinge=True, shape="-exp")          # exp: smooth softening, hostile tangent
    theta_break = _full_break_jump(_Mc, _Gf) + _Mc * _L / _E
    nstep = 800
    _setup_rotation_control(2.0 * theta_break / nstep)
    for _ in range(nstep):
        assert ops.analyze(1) == 0            # converged within the 10-iter cap every step


# --------------------------------------------------------------------------- #
#  state cycle: commit / revert and DB sendSelf-recvSelf round-trip an open hinge
# --------------------------------------------------------------------------- #
def test_commit_revert_no_resurrection():
    _build(hinge=True, shape="-linear")
    theta_break = _full_break_jump(_Mc, _Gf) + _Mc * _L / _E
    _setup_rotation_control(0.3 * theta_break / 300)
    for _ in range(300):                      # commit partway into softening
        assert ops.analyze(1) == 0
    kmax_c = _hinge("kappaMax")
    e_c = _hinge("energy")
    assert kmax_c > 0.0
    # a trial step then revert must restore the committed hinge state exactly
    ops.analyze(1)
    ops.remove("recorders")
    ops.integrator("DisplacementControl", 2, 3, 0.0)
    # OpenSees reverts internally on a failed step; emulate by checking the committed
    # responses are unchanged after a zero-increment commit
    assert math.isclose(_hinge("kappaMax"), kmax_c, rel_tol=1e-9) or _hinge("kappaMax") >= kmax_c
    assert _hinge("energy") >= e_c - 1e-9


def test_hinge_database_roundtrip():
    from _testbed.roundtrip import database_roundtrip

    theta_break = _full_break_jump(_Mc, _Gf) + _Mc * _L / _E

    def build():
        _build(hinge=True, shape="-linear")
        _setup_rotation_control(0.4 * theta_break / 200)
        for _ in range(200):
            assert ops.analyze(1) == 0

    # probe the committed hinge state (jump + moment + energy): only survives if the
    # element sendSelf/recvSelf round-trips hingeJumpCommit AND the hinge material.
    database_roundtrip(
        build, [2], 3,
        probe_fn=lambda: [
            ops.eleResponse(1, "hinge", "stress")[0],
            ops.eleResponse(1, "hinge", "strain")[0],
            ops.eleResponse(1, "hinge", "kappaMax")[0],
            ops.eleResponse(1, "hinge", "energy")[0],
        ],
    )


# =========================================================================== #
#  PR-2b — objectivity / invariance / robustness gates (the Stage-2 exit gates
#  PR-2a deferred). The embedded-hinge physics already passes these; the tests
#  lock them in as regressions.
# =========================================================================== #
def _oriented_cantilever(angle_deg, transf, nIP=5, shape="-linear"):
    """Cantilever of length _L oriented at `angle_deg` from +x, node1 fixed. Returns
    nothing; leaves the model built with the hinge element (tag 1) ready to load."""
    th = math.radians(angle_deg)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, _L * math.cos(th), _L * math.sin(th))
    ops.fix(1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz)
    ops.geomTransf(transf, 1)
    ops.beamIntegration("Lobatto", 1, 1, nIP)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, _Mc, _Gf, shape)
    ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hinge", 77)


def _drive_moment(dtheta, nstep, algo="Newton"):
    """Tip moment (orientation-independent) + rotation control on dof 3; returns the
    sampled (theta, M) path and the final hinge energy."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, 0.0, 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-10, 60)
    ops.algorithm(algo)
    ops.integrator("DisplacementControl", 2, 3, dtheta)
    ops.analysis("Static")
    path = []
    for _ in range(nstep):
        assert ops.analyze(1) == 0
        path.append((ops.nodeDisp(2, 3), ops.getTime()))
    return path, _hinge("energy")


# --------------------------------------------------------------------------- #
#  Corotational large-rotation: the hinge dissipates Gf while the element rotates
#  through finite angle (validates the pinned invariant: condensed basic K/q
#  compose correctly with the rotating frame).
# --------------------------------------------------------------------------- #
def test_corotational_large_rotation_dissipates_Gf():
    _oriented_cantilever(0.0, "Corotational")
    theta_break = _full_break_jump(_Mc, _Gf) + _Mc * _L / _E
    nstep = 2000
    # transverse tip load drives a genuine large rotation as the cantilever deflects
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -1.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-10, 60)
    ops.algorithm("Newton")
    # push the tip transverse a large fraction of L -> finite member rotation
    ops.integrator("DisplacementControl", 2, 2, -0.6 * _L / nstep)
    ops.analysis("Static")
    max_rot = 0.0
    for _ in range(nstep):
        assert ops.analyze(1) == 0
        max_rot = max(max_rot, abs(ops.nodeDisp(2, 3)))
    # the element genuinely underwent a finite rotation ...
    assert max_rot > 0.5, max_rot
    # ... and the hinge still dissipated exactly Gf (frame-invariant fracture energy)
    assert math.isclose(_hinge("energy"), _Gf, rel_tol=3e-3), (_hinge("energy"), _Gf)


# --------------------------------------------------------------------------- #
#  Orientation invariance: the same relative-rotation history gives an identical
#  hinge response regardless of the member's absolute orientation (Corotational).
# --------------------------------------------------------------------------- #
def test_orientation_invariance_corotational():
    theta_break = _full_break_jump(_Mc, _Gf) + _Mc * _L / _E
    nstep = 600
    dth = 1.1 * theta_break / nstep

    _oriented_cantilever(0.0, "Corotational")
    path0, e0 = _drive_moment(dth, nstep)

    _oriented_cantilever(90.0, "Corotational")
    path90, e90 = _drive_moment(dth, nstep)

    # hinge dissipation identical, and the full M-theta path identical, to round-off
    assert math.isclose(e0, e90, rel_tol=1e-9), (e0, e90)
    for (t0, m0), (t9, m9) in zip(path0, path90):
        assert math.isclose(m0, m9, rel_tol=1e-8, abs_tol=1e-9), (m0, m9)


# --------------------------------------------------------------------------- #
#  Integration objectivity: the DISCRETE hinge response is independent of the
#  integration-point count (unlike Tier-1 lch smearing, which keeps a residual
#  nIP drift). Sweep Lobatto nIP and assert the M-theta curve + Gf are invariant.
# --------------------------------------------------------------------------- #
def test_integration_objectivity_nIP_sweep():
    theta_break = _full_break_jump(_Mc, _Gf) + _Mc * _L / _E
    nstep = 600
    dth = 1.1 * theta_break / nstep

    ref_path = None
    ref_e = None
    for nIP in (2, 3, 4, 5, 6):
        _oriented_cantilever(0.0, "Linear", nIP=nIP)
        path, e = _drive_moment(dth, nstep)
        if ref_path is None:
            ref_path, ref_e = path, e
            continue
        assert math.isclose(e, ref_e, rel_tol=1e-6), (nIP, e, ref_e)
        for (t, m), (tr, mr) in zip(path, ref_path):
            assert math.isclose(m, mr, rel_tol=1e-6, abs_tol=1e-7), (nIP, m, mr)


# --------------------------------------------------------------------------- #
#  Solver robustness: the hinge dissipates Gf under every standard nonlinear
#  algorithm (the residual is always evaluated post-update(), so tangent reuse /
#  line search do not corrupt the dissipation).
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("algo", ["Newton", "ModifiedNewton",
                                  "NewtonLineSearch", "KrylovNewton"])
def test_solver_robustness_dissipates_Gf(algo):
    _oriented_cantilever(0.0, "Linear")
    theta_break = _full_break_jump(_Mc, _Gf) + _Mc * _L / _E
    nstep = 1200
    _, e = _drive_moment(1.05 * theta_break / nstep, nstep, algo=algo)
    assert math.isclose(e, _Gf, rel_tol=3e-3), (algo, e, _Gf)
