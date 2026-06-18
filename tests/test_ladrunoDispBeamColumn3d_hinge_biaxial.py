"""LadrunoDispBeamColumn3d Tier-2 BIAXIAL embedded hinge — ADR 33 PR-3b.

Adds the weak-axis (My) rotation jump `alpha_y` to the 3D strong-axis hinge (PR-3a),
making the embedded hinge biaxial: `alpha = [alpha_z, alpha_y]`, one jump per bending
axis, carried by two block-diagonal cohesive laws (`-hinge $matZ -hingeY $matY`). The
two jumps are solved by ONE unified inner Newton (a single setTrialSectionDeformation
per IP sets BOTH bulk curvatures `kappa_*_bulk = B_*v - alpha_*/L`) and condensed by a
TRUE coupled 2x2 `K_aa` whose off-diagonal `(1/L)sum wt*ks(MZ,MY)` is exact when the
fiber section couples bending axes (product of inertia Iyz != 0, generic for any
axially-loaded biaxially-bent section) — using an EIGENVALUE-FLOORED 2x2 inverse, not
two independent rank-1 updates, not a det-floored adjugate.

PR-3b gates:
  * parse: -hingeY requires -hinge; -hingeY/-nl mutually exclusive
  * reduce-to-PR-3a: biaxial element with only Mz loaded -> alpha_y ~ 0, Mz dissipates Gf_z
  * independent weak axis: pure My load -> alpha_z ~ 0, My dissipates Gf_y
  * constant biaxial patch test: each hinge carries exactly its applied moment
  * energy partition: drive z then y -> hingeZ == Gf_z, hingeY == Gf_y
  * THE coupled-2x2 gate: a product-of-inertia fiber section driven through STAGGERED
    biaxial softening converges under tight Newton (the consistent coupled tangent is
    the oracle; two independent rank-1 updates lose convergence when ks(MZ,MY) != 0)
  * finite-rotation invariance under CorotCrdTransf3d (biaxial, both hinges open)
  * solver robustness; integration-objectivity nIP sweep; commit/revert + DB roundtrip

Plan: Ladruno_implementation/33_ladruno_dispbeamcolumn3d_hinge_adr.md (PR-3b).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_L = 100.0
_E = 1.0e5
_A = 1.0e8         # near-rigid axial
_Iz = 1.0          # EIz = 1e5
_Iy = 1.0          # EIy = 1e5
_G = 4.0e4
_J = 2.0
_McZ, _GfZ = 50.0, 20.0
_McY, _GfY = 35.0, 12.0
_RY = 5            # nodal DOF for rotation about y (My) in ndf6 [ux,uy,uz,rx,ry,rz]
_RZ = 6            # nodal DOF for rotation about z (Mz)
_UY = 2
_UZ = 3

# Centred fiber pattern with a genuine product of inertia Iyz != 0 (rho ~ 0.43) so the
# coupled 2x2 K_aa has a non-trivial off-diagonal even with an elastic fiber material.
#   Iz = sum A*y^2 = 8.32, Iy = sum A*z^2 = 12.5, Iyz = sum A*y*z = 4.4, centroid at origin.
_FIBERS = [(2.0, 1.5, 1.0), (-2.0, -1.5, 1.0), (0.4, -2.0, 1.0), (-0.4, 2.0, 1.0)]


def _kpen_default(Mc, Gf, ratio=1000.0):
    return ratio * Mc * Mc / (2.0 * Gf)


def _full_break_jump(Mc, Gf):
    Kpen = _kpen_default(Mc, Gf)
    k0 = Mc / Kpen
    kf = 2.0 * (Gf - Mc * Mc / (2.0 * Kpen)) / Mc
    return k0 + kf


def _theta_break(Mc, Gf):
    return _full_break_jump(Mc, Gf) + Mc * _L / _E


def _build(transf="Linear", nIP=5):
    """Biaxial hinge on an elastic (uncoupled) section — clean per-axis behaviour."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz, _Iy, _G, _J)
    ops.geomTransf(transf, 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, nIP)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, _McZ, _GfZ, "-linear")
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 78, _McY, _GfY, "-linear")
    ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hinge", 77, "-hingeY", 78)


def _build_coupled(transf="Linear", nIP=4, shape="-linear"):
    """Biaxial hinge on a product-of-inertia fiber section -> ks(MZ,MY) != 0."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 9, _E)
    ops.section("Fiber", 1, "-GJ", _G * _J)
    for (y, z, a) in _FIBERS:
        ops.fiber(y, z, a, 9)
    ops.geomTransf(transf, 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, nIP)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, _McZ, _GfZ, shape)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 78, _McY, _GfY, shape)
    ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hinge", 77, "-hingeY", 78)


def _moment_control(dofs_loads, control_dof, dtheta, algo="Newton", tol=1.0e-10):
    """Apply tip moment load(s) (dict dof->magnitude) + rotation control on control_dof."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    load = [0.0] * 6
    for dof, mag in dofs_loads.items():
        load[dof - 1] = mag
    ops.load(2, *load)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", tol, 60)
    ops.algorithm(algo)
    ops.integrator("DisplacementControl", 2, control_dof, dtheta)
    ops.analysis("Static")


def _prescribe(sp_map, nstep, tol=1.0e-11, maxiter=10, algo="Newton"):
    """Ramp PRESCRIBED nodal displacements (dof->target) on node 2 via sp constraints.

    Displacement-prescribing the softening DOFs removes the single-element snap a free
    DOF under a dead load suffers at hinge activation, so a single element traces deep
    (even simultaneous biaxial) softening robustly — and tight Newton through it is the
    coupled-tangent oracle.
    """
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for dof, target in sp_map.items():
        ops.sp(2, dof, target)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", tol, maxiter)
    ops.algorithm(algo)
    ops.integrator("LoadControl", 1.0 / nstep)
    ops.analysis("Static")


def _hz(resp):
    return ops.eleResponse(1, "hinge", resp)[0]


def _hy(resp):
    return ops.eleResponse(1, "hingeY", resp)[0]


# --------------------------------------------------------------------------- #
#  parse / mutual exclusion
# --------------------------------------------------------------------------- #
def _model_skeleton():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz, _Iy, _G, _J)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, _McZ, _GfZ, "-linear")
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 78, _McY, _GfY, "-linear")


def test_hingeY_requires_hinge():
    _model_skeleton()
    with pytest.raises(Exception):
        ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hingeY", 78)


def test_hingeY_and_nl_mutually_exclusive():
    _model_skeleton()
    with pytest.raises(Exception):
        ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hinge", 77, "-hingeY", 78, "-nl")


def test_biaxial_builds():
    _build()
    assert ops.eleResponse(1, "hinge", "stress")[0] == 0.0
    assert ops.eleResponse(1, "hingeY", "stress")[0] == 0.0


# --------------------------------------------------------------------------- #
#  reduce-to-PR-3a: pure Mz load -> weak-axis jump stays closed
# --------------------------------------------------------------------------- #
def test_pure_Mz_leaves_weak_axis_closed():
    _build()
    tb = _theta_break(_McZ, _GfZ)
    nstep = 1500
    _moment_control({_RZ: 1.0}, _RZ, 1.2 * tb / nstep)
    peak = 0.0
    for _ in range(nstep):
        assert ops.analyze(1) == 0
        peak = max(peak, _hz("stress"))
    assert math.isclose(_hz("energy"), _GfZ, rel_tol=2e-3), (_hz("energy"), _GfZ)
    assert math.isclose(peak, _McZ, rel_tol=2e-3), (peak, _McZ)
    # the weak-axis hinge never saw a moment -> closed, zero jump, zero dissipation
    assert abs(_hy("strain")) <= 1e-10, _hy("strain")
    assert abs(_hy("energy")) <= 1e-6 * _GfY, _hy("energy")


# --------------------------------------------------------------------------- #
#  independent weak axis: pure My load -> strong-axis jump stays closed
# --------------------------------------------------------------------------- #
def test_pure_My_dissipates_GfY_and_leaves_strong_axis_closed():
    _build()
    tb = _theta_break(_McY, _GfY)
    nstep = 1500
    _moment_control({_RY: 1.0}, _RY, 1.2 * tb / nstep)
    peak = 0.0
    for _ in range(nstep):
        assert ops.analyze(1) == 0
        peak = max(peak, _hy("stress"))
    assert math.isclose(_hy("energy"), _GfY, rel_tol=2e-3), (_hy("energy"), _GfY)
    assert math.isclose(peak, _McY, rel_tol=2e-3), (peak, _McY)
    assert abs(_hz("strain")) <= 1e-10, _hz("strain")
    assert abs(_hz("energy")) <= 1e-6 * _GfZ, _hz("energy")


# --------------------------------------------------------------------------- #
#  constant biaxial patch test: each hinge carries exactly its applied moment
# --------------------------------------------------------------------------- #
def test_constant_biaxial_patch():
    _build()
    # both moments held below their peaks via LoadControl (structure stays ~elastic)
    Mz_load, My_load = 0.6 * _McZ, 0.6 * _McY
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, 0.0, 0.0, 0.0, My_load, Mz_load)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-10, 60)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / 5.0)
    ops.analysis("Static")
    for _ in range(5):
        assert ops.analyze(1) == 0
        lam = ops.getTime()
        assert math.isclose(_hz("stress"), lam * Mz_load, rel_tol=1e-8, abs_tol=1e-8), (
            _hz("stress"), lam * Mz_load)
        assert math.isclose(_hy("stress"), lam * My_load, rel_tol=1e-8, abs_tol=1e-8), (
            _hy("stress"), lam * My_load)


# --------------------------------------------------------------------------- #
#  energy partition: drive z to break, then y to break -> each hinge == its Gf
# --------------------------------------------------------------------------- #
def test_biaxial_energy_partition():
    _build()
    tbz, tby = _theta_break(_McZ, _GfZ), _theta_break(_McY, _GfY)
    nz, ny = 1400, 1400
    _moment_control({_RZ: 1.0}, _RZ, 1.05 * tbz / nz)
    for _ in range(nz):
        assert ops.analyze(1) == 0
    # remove the Mz pattern and drive the weak axis with a FRESH My pattern from zero load
    # factor (the z-hinge stays broken/committed; a single load factor can't trace both axes).
    ops.remove("loadPattern", 1)
    ops.setTime(0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(2, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
    ops.integrator("DisplacementControl", 2, _RY, 1.05 * tby / ny)
    for _ in range(ny):
        assert ops.analyze(1) == 0
    assert math.isclose(_hz("energy"), _GfZ, rel_tol=3e-3), (_hz("energy"), _GfZ)
    assert math.isclose(_hy("energy"), _GfY, rel_tol=3e-3), (_hy("energy"), _GfY)


# --------------------------------------------------------------------------- #
#  THE coupled-2x2 gate: staggered biaxial softening on a product-of-inertia
#  section converges under tight Newton. With ks(MZ,MY) != 0 the consistent
#  condensed tangent REQUIRES the true coupled 2x2; two independent rank-1
#  updates are the exact tangent of a different residual and lose convergence
#  through coupled softening. Tight Newton through it is the tangent oracle.
# --------------------------------------------------------------------------- #
def test_coupled_staggered_softening_converges_tightly():
    _build_coupled(shape="-exp")
    # Product-of-inertia section (ks(MZ,MY) != 0): the condensed tangent REQUIRES the true
    # coupled 2x2 (two independent rank-1 updates are the exact tangent of a different residual
    # and lose convergence through coupled softening). Prescribe BOTH rotations to deep softening
    # (z and y activate at different lambda -> STAGGERED biaxial); tight Newton (maxIter 10) must
    # converge every step. Convergence through the coupled descending branch is the tangent oracle.
    tbz, tby = _theta_break(_McZ, _GfZ), _theta_break(_McY, _GfY)
    nstep = 1500
    _prescribe({_RZ: 1.3 * tbz, _RY: 1.3 * tby}, nstep, tol=1.0e-11, maxiter=10)
    for _ in range(nstep):
        assert ops.analyze(1) == 0
    # both jumps fully opened and dissipated their fracture energy (exp tail truncation ~6%)
    assert math.isclose(_hz("energy"), _GfZ, rel_tol=8e-2), (_hz("energy"), _GfZ)
    assert math.isclose(_hy("energy"), _GfY, rel_tol=8e-2), (_hy("energy"), _GfY)


# --------------------------------------------------------------------------- #
#  finite-rotation invariance under CorotCrdTransf3d (biaxial, both hinges open)
# --------------------------------------------------------------------------- #
def test_corotational_large_rotation_biaxial():
    # Prescribe SKEW tip transverse displacements (uy, uz) under CorotCrdTransf3d: the member
    # rotates through finite angle (> 1 rad) while BOTH bending hinges open. The condensed biaxial
    # kb passes through the quaternion-mean triad -> the pinned invariant must hold biaxially.
    _build(transf="Corotational")
    nstep = 2500
    _prescribe({_UY: 0.45 * _L, _UZ: 0.30 * _L}, nstep, tol=1.0e-10, maxiter=40)
    max_rot = 0.0
    for _ in range(nstep):
        assert ops.analyze(1) == 0
        max_rot = max(max_rot, math.hypot(ops.nodeDisp(2, _RY), ops.nodeDisp(2, _RZ)))
    assert max_rot > 0.5, max_rot
    # both bending hinges fully dissipated their Gf while the element rotated through finite angle
    assert math.isclose(_hz("energy"), _GfZ, rel_tol=3e-3), (_hz("energy"), _GfZ)
    assert math.isclose(_hy("energy"), _GfY, rel_tol=5e-2), (_hy("energy"), _GfY)


@pytest.mark.parametrize("algo", ["Newton", "ModifiedNewton", "KrylovNewton"])
def test_solver_robustness_biaxial(algo):
    # Both hinges driven to full break together (prescribed rotations) under three algorithms.
    _build()
    tbz, tby = _theta_break(_McZ, _GfZ), _theta_break(_McY, _GfY)
    nstep = 1500
    _prescribe({_RZ: 1.3 * tbz, _RY: 1.3 * tby}, nstep, tol=1.0e-10, maxiter=20, algo=algo)
    for _ in range(nstep):
        assert ops.analyze(1) == 0
    assert math.isclose(_hz("energy"), _GfZ, rel_tol=8e-2), (algo, _hz("energy"), _GfZ)
    assert math.isclose(_hy("energy"), _GfY, rel_tol=8e-2), (algo, _hy("energy"), _GfY)


# --------------------------------------------------------------------------- #
#  integration objectivity: biaxial energy partition is nIP-invariant
# --------------------------------------------------------------------------- #
def test_integration_objectivity_biaxial_nIP():
    tbz = _theta_break(_McZ, _GfZ)
    nstep = 500
    dth = 0.85 * tbz / nstep
    ref = None
    for nIP in (2, 3, 4, 5):
        _build(nIP=nIP)
        _moment_control({_RZ: 1.0}, _RZ, dth)
        for _ in range(nstep):
            assert ops.analyze(1) == 0
        e = _hz("energy")
        if ref is None:
            ref = e
        else:
            assert math.isclose(e, ref, rel_tol=1e-6), (nIP, e, ref)


# --------------------------------------------------------------------------- #
#  state cycle: DB sendSelf/recvSelf round-trip with BOTH hinges open
# --------------------------------------------------------------------------- #
def test_biaxial_database_roundtrip():
    from _testbed.roundtrip import database_roundtrip

    tbz, tby = _theta_break(_McZ, _GfZ), _theta_break(_McY, _GfY)

    def build():
        # both hinges partially open (into softening but not fully broken) for a non-trivial state
        _build()
        _prescribe({_RZ: 0.45 * tbz, _RY: 0.45 * tby}, 300, tol=1.0e-11, maxiter=10)
        for _ in range(300):
            assert ops.analyze(1) == 0

    database_roundtrip(
        build, [2], 6,
        probe_fn=lambda: [
            ops.eleResponse(1, "hinge", "stress")[0],
            ops.eleResponse(1, "hinge", "strain")[0],
            ops.eleResponse(1, "hinge", "energy")[0],
            ops.eleResponse(1, "hingeY", "stress")[0],
            ops.eleResponse(1, "hingeY", "strain")[0],
            ops.eleResponse(1, "hingeY", "energy")[0],
        ],
    )
