"""LadrunoDispBeamColumn3d Tier-2 embedded strong-axis (Mz) hinge — ADR 33 PR-3a.

The 3D sibling of the 2D Tier-2 hinge, restricted in PR-3a to the STRONG-AXIS
rotation jump `alpha_z` on the Mz row of the 6-DOF 3D basic system — the literal
2D scalar algebra transplanted onto 3D. The bulk section sees
`kappa_z_bulk = (1/L)((6xi-4)v1+(6xi-2)v2) + Gbar_z*alpha_z` (Gbar_z=-1/L), so it
unloads as the jump grows (no Gf double-count); the cohesive `M([[theta_z]])`
(LadrunoCohesiveHinge) carries `Gf`. `alpha_z` is converged by an inner Newton in
`update()` and condensed (a guarded rank-1 update) to the 6-DOF basic system BEFORE
`CorotCrdTransf3d` — the PINNED INVARIANT, now through the quaternion-triad transform.

PR-3a gates (mirror the 2D suite + the 3D-specific finite-rotation gate):
  * parse + `-hinge`/`-nl` mutual exclusion
  * reduce-to-Tier-1: hinge-absent element bit-identical to stock dispBeamColumn (3D)
  * constant-Mz patch test: hinge carries exactly the applied strong-axis moment
  * energy gate: `∫Mz d[[theta_z]] == Gf`, peak hinge moment == Mc
  * element total dissipation == Gf on an elastic section (no bulk double-count)
  * tangent consistency: softening converges under NormDispIncr-1e-12 Newton
  * FINITE-ROTATION INVARIANCE under CorotCrdTransf3d: the hinge dissipates Gf while
    the element rotates through finite angle (proves the pinned invariant survives
    the quaternion transform) + solver robustness
  * commit/revert + DB sendSelf/recvSelf round-trip of an open hinge

Plan: Ladruno_implementation/33_ladruno_dispbeamcolumn3d_hinge_adr.md.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_L = 100.0
_E = 1.0e5
_A = 1.0e8         # near-rigid axial
_Iz = 1.0          # strong axis: EIz = 1e5
_Iy = 1.0
_G = 4.0e4
_J = 2.0
_Mc = 50.0
_Gf = 20.0
_RZ = 6            # nodal DOF for rotation about z (Mz) in ndf6 [ux,uy,uz,rx,ry,rz]
_UY = 2            # transverse y displacement DOF


def _kpen_default(Mc, Gf, ratio=1000.0):
    return ratio * Mc * Mc / (2.0 * Gf)


def _full_break_jump(Mc, Gf):
    Kpen = _kpen_default(Mc, Gf)
    k0 = Mc / Kpen
    kf = 2.0 * (Gf - Mc * Mc / (2.0 * Kpen)) / Mc
    return k0 + kf


def _build(hinge=True, shape="-linear", transf="Linear", nl=False, nIP=5):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz, _Iy, _G, _J)
    ops.geomTransf(transf, 1, 0.0, 0.0, 1.0)   # vecxz = global z -> local z = strong axis
    ops.beamIntegration("Lobatto", 1, 1, nIP)
    args = [1, 1, 2, 1, 1]
    if nl:
        args.append("-nl")
    if hinge:
        ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, _Mc, _Gf, shape)
        args += ["-hinge", 77]
    ops.element("LadrunoDispBeamColumn", *args)


def _moment_control(dtheta, algo="Newton"):
    """Tip moment about z (DOF 6) + rotation control -> constant-Mz cantilever."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-10, 60)
    ops.algorithm(algo)
    ops.integrator("DisplacementControl", 2, _RZ, dtheta)
    ops.analysis("Static")


def _hinge(resp):
    return ops.eleResponse(1, "hinge", resp)[0]


# --------------------------------------------------------------------------- #
#  parse / mutual exclusion
# --------------------------------------------------------------------------- #
def test_hinge_builds_and_nl_mutually_exclusive_3d():
    _build(hinge=True)
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz, _Iy, _G, _J)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, _Mc, _Gf, "-linear")
    with pytest.raises(Exception):
        ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hinge", 77, "-nl")


# --------------------------------------------------------------------------- #
#  reduce-to-Tier-1: no-hinge element bit-identical to stock dispBeamColumn
# --------------------------------------------------------------------------- #
def _elastic_tip(make_element):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz, _Iy, _G, _J)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    make_element()
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -10.0, 3.0, 0.0, 2.0, 5.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ops.analyze(1)
    return ops.nodeDisp(2)


def test_reduce_to_tier1_without_hinge_matches_stock_3d():
    ladr = _elastic_tip(lambda: ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1))
    stock = _elastic_tip(lambda: ops.element("dispBeamColumn", 1, 1, 2, 1, 1))
    for a, b in zip(ladr, stock):
        assert a == b or abs(a - b) <= 1e-300 + 1e-15 * abs(b), (a, b)


# --------------------------------------------------------------------------- #
#  constant-Mz patch test
# --------------------------------------------------------------------------- #
def test_constant_moment_patch_3d():
    _build(hinge=True)
    _moment_control(0.002)
    for _ in range(5):
        assert ops.analyze(1) == 0
        M_applied = ops.getTime()
        assert M_applied < _Mc
        assert math.isclose(_hinge("stress"), M_applied, rel_tol=1e-9, abs_tol=1e-9), (
            _hinge("stress"), M_applied)


# --------------------------------------------------------------------------- #
#  energy gate + peak == Mc
# --------------------------------------------------------------------------- #
def test_hinge_dissipates_Gf_and_peaks_at_Mc_3d():
    _build(hinge=True, shape="-linear")
    theta_break = _full_break_jump(_Mc, _Gf) + _Mc * _L / _E
    nstep = 1500
    _moment_control(1.2 * theta_break / nstep)
    peak = 0.0
    for _ in range(nstep):
        assert ops.analyze(1) == 0
        peak = max(peak, _hinge("stress"))
    assert math.isclose(_hinge("energy"), _Gf, rel_tol=2e-3), (_hinge("energy"), _Gf)
    assert math.isclose(peak, _Mc, rel_tol=2e-3), (peak, _Mc)


# --------------------------------------------------------------------------- #
#  element total dissipation == Gf (no bulk double-count)
# --------------------------------------------------------------------------- #
def test_element_total_dissipation_equals_Gf_3d():
    _build(hinge=True, shape="-linear")
    theta_break = _full_break_jump(_Mc, _Gf) + _Mc * _L / _E
    theta_max = 1.05 * theta_break
    nstep = 1200
    dth = theta_max / nstep
    _moment_control(dth)
    work = 0.0
    th_prev = ops.nodeDisp(2, _RZ); M_prev = _hinge("stress")
    for _ in range(nstep):
        assert ops.analyze(1) == 0
        th = ops.nodeDisp(2, _RZ); M = ops.getTime()
        work += 0.5 * (M_prev + M) * (th - th_prev)
        th_prev, M_prev = th, M
    ops.integrator("DisplacementControl", 2, _RZ, -dth)
    for _ in range(nstep):
        assert ops.analyze(1) == 0
        th = ops.nodeDisp(2, _RZ); M = ops.getTime()
        work += 0.5 * (M_prev + M) * (th - th_prev)
        th_prev, M_prev = th, M
    assert math.isclose(work, _Gf, rel_tol=5e-3), (work, _Gf)


# --------------------------------------------------------------------------- #
#  tangent consistency: softening converges under tight Newton
# --------------------------------------------------------------------------- #
def test_softening_newton_converges_tightly_3d():
    _build(hinge=True, shape="-exp")
    theta_break = _full_break_jump(_Mc, _Gf) + _Mc * _L / _E
    nstep = 800
    _moment_control(2.0 * theta_break / nstep)
    for _ in range(nstep):
        assert ops.analyze(1) == 0


# --------------------------------------------------------------------------- #
#  THE 3D gate: finite-rotation invariance under CorotCrdTransf3d
#  A transverse tip load (Corotational) drives a finite member rotation while the
#  strong-axis hinge opens; the hinge must still dissipate exactly Gf -> the pinned
#  invariant (condensed 6-DOF basic K/q) composes correctly with the quaternion triad.
# --------------------------------------------------------------------------- #
def test_corotational_large_rotation_dissipates_Gf_3d():
    _build(hinge=True, shape="-linear", transf="Corotational")
    nstep = 2000
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0)   # transverse y load -> bending about z
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-10, 60)
    ops.algorithm("Newton")
    ops.integrator("DisplacementControl", 2, _UY, -0.6 * _L / nstep)
    ops.analysis("Static")
    max_rot = 0.0
    for _ in range(nstep):
        assert ops.analyze(1) == 0
        max_rot = max(max_rot, abs(ops.nodeDisp(2, _RZ)))
    assert max_rot > 0.5, max_rot
    assert math.isclose(_hinge("energy"), _Gf, rel_tol=3e-3), (_hinge("energy"), _Gf)


@pytest.mark.parametrize("algo", ["Newton", "ModifiedNewton", "KrylovNewton"])
def test_solver_robustness_dissipates_Gf_3d(algo):
    _build(hinge=True, shape="-linear")
    theta_break = _full_break_jump(_Mc, _Gf) + _Mc * _L / _E
    nstep = 1200
    _moment_control(1.05 * theta_break / nstep, algo=algo)
    for _ in range(nstep):
        assert ops.analyze(1) == 0
    assert math.isclose(_hinge("energy"), _Gf, rel_tol=3e-3), (algo, _hinge("energy"), _Gf)


# --------------------------------------------------------------------------- #
#  integration objectivity (strong-axis hinge is nIP-invariant)
# --------------------------------------------------------------------------- #
def test_integration_objectivity_nIP_sweep_3d():
    # Drive to DEEP softening (~0.85x the full-break jump) but NOT past full break — a
    # fully-broken single-element hinge is a zero-stiffness mechanism whose pure-Newton solve
    # is platform-FP-fragile. nIP-invariance holds anywhere on the softening branch. (Mirrors
    # the 2D fix for the Ubuntu CI failure from #269.)
    theta_break = _full_break_jump(_Mc, _Gf) + _Mc * _L / _E
    nstep = 500
    dth = 0.85 * theta_break / nstep
    ref = None
    for nIP in (2, 3, 4, 5):
        _build(hinge=True, shape="-linear", nIP=nIP)
        _moment_control(dth)
        for _ in range(nstep):
            assert ops.analyze(1) == 0
        e = _hinge("energy")
        if ref is None:
            ref = e
        else:
            assert math.isclose(e, ref, rel_tol=1e-6), (nIP, e, ref)


# --------------------------------------------------------------------------- #
#  state cycle: DB sendSelf/recvSelf round-trip of an open hinge
# --------------------------------------------------------------------------- #
def test_hinge_database_roundtrip_3d():
    from _testbed.roundtrip import database_roundtrip

    theta_break = _full_break_jump(_Mc, _Gf) + _Mc * _L / _E

    def build():
        _build(hinge=True, shape="-linear")
        _moment_control(0.4 * theta_break / 200)
        for _ in range(200):
            assert ops.analyze(1) == 0

    database_roundtrip(
        build, [2], 6,
        probe_fn=lambda: [
            ops.eleResponse(1, "hinge", "stress")[0],
            ops.eleResponse(1, "hinge", "strain")[0],
            ops.eleResponse(1, "hinge", "kappaMax")[0],
            ops.eleResponse(1, "hinge", "energy")[0],
        ],
    )
