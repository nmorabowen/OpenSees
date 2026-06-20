"""LadrunoDispBeamColumn3d + LadrunoCohesiveHingeBiaxial — coupled Mz-My hinge (ADR 34).

The coupled biaxial cohesive law (`-hingeBiaxial $ndMatTag`, an order-2 NDMaterial
`LadrunoCohesiveHingeBiaxial`) replaces the two block-diagonal scalar hinges (PR-3b)
with a single law carrying an **Mz-My interaction surface** (ellipse). It reduces
EXACTLY to the scalar `LadrunoCohesiveHinge` on each pure bending axis, and for
combined loading activates on the ellipse `√((Mz/Mcz)²+(My/Mcy)²)=1` (earlier than
the independent per-axis capacities) with a mixed-mode fracture energy.

Gates:
  * parse + exclusivity (-hingeBiaxial vs -hinge/-hingeY/-nl)
  * on-axis reduction: pure-Mz coupled element == scalar -hinge element (Mz path, Gf_z),
    weak axis closed; symmetric pure-My
  * elliptical interaction onset: combined proportional load activates on the ellipse
  * mixed-mode dissipation: combined radial break dissipates Gf between Gf_y and Gf_z
  * coupled tangent: product-of-inertia section, deep biaxial softening, tight Newton converges
  * finite-rotation invariance under CorotCrdTransf3d; DB roundtrip

Plan: Ladruno_implementation/34_ladruno_cohesive_hinge_biaxial_adr.md.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_L = 100.0
_E = 1.0e5
_A = 1.0e8
_Iz = 1.0
_Iy = 1.0
_G = 4.0e4
_J = 2.0
_McZ, _GfZ = 50.0, 20.0
_McY, _GfY = 35.0, 12.0
_RY, _RZ, _UY, _UZ = 5, 6, 2, 3

_FIBERS = [(2.0, 1.5, 1.0), (-2.0, -1.5, 1.0), (0.4, -2.0, 1.0), (-0.4, 2.0, 1.0)]


def _kpen(Mc, Gf, ratio=1000.0):
    return ratio * Mc * Mc / (2.0 * Gf)


def _theta_break(Mc, Gf):
    Kpen = _kpen(Mc, Gf)
    k0 = Mc / Kpen
    kf = 2.0 * (Gf - Mc * Mc / (2.0 * Kpen)) / Mc
    return k0 + kf + Mc * _L / _E


def _model(section="elastic"):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    if section == "fiber":
        ops.uniaxialMaterial("Elastic", 9, _E)
        ops.section("Fiber", 1, "-GJ", _G * _J)
        for (y, z, a) in _FIBERS:
            ops.fiber(y, z, a, 9)
    else:
        ops.section("Elastic", 1, _E, _A, _Iz, _Iy, _G, _J)
    ops.geomTransf("Linear" if section != "corot" else "Corotational", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)


def _build_coupled(section="elastic", shape="-linear", transf="Linear"):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    if section == "fiber":
        ops.uniaxialMaterial("Elastic", 9, _E)
        ops.section("Fiber", 1, "-GJ", _G * _J)
        for (y, z, a) in _FIBERS:
            ops.fiber(y, z, a, 9)
    else:
        ops.section("Elastic", 1, _E, _A, _Iz, _Iy, _G, _J)
    ops.geomTransf(transf, 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 79, _McZ, _GfZ, _McY, _GfY, shape)
    ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hingeBiaxial", 79)


def _build_blockdiag(shape="-linear"):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz, _Iy, _G, _J)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, _McZ, _GfZ, shape)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 78, _McY, _GfY, shape)
    ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hinge", 77, "-hingeY", 78)


def _prescribe(sp_map, nstep, tol=1.0e-11, maxiter=10, algo="Newton"):
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


def _hc(resp):
    return ops.eleResponse(1, "hingeBiaxial", resp)[0]


# --------------------------------------------------------------------------- #
#  parse / exclusivity
# --------------------------------------------------------------------------- #
def test_coupled_builds():
    _build_coupled()
    assert _hc("momentZ") == 0.0 and _hc("momentY") == 0.0


def test_coupled_exclusive_with_blockdiag_and_nl():
    _model()
    ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 79, _McZ, _GfZ, _McY, _GfY, "-linear")
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 77, _McZ, _GfZ, "-linear")
    with pytest.raises(Exception):
        ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hinge", 77, "-hingeBiaxial", 79)
    _model()
    ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 79, _McZ, _GfZ, _McY, _GfY, "-linear")
    with pytest.raises(Exception):
        ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hingeBiaxial", 79, "-nl")


# --------------------------------------------------------------------------- #
#  on-axis reduction: coupled element == scalar -hinge element on a pure axis
# --------------------------------------------------------------------------- #
def _drive_pure(dof, target, nstep, builder):
    builder()
    _prescribe({dof: target}, nstep)
    Ms, ths = [], []
    for _ in range(nstep):
        assert ops.analyze(1) == 0
        ths.append(ops.nodeDisp(2, dof))
        # tip reaction moment about the driven axis = element basic moment carried by the hinge
        Ms.append(ops.eleForce(1, 5 if dof == _RY else 6))
    return ths, Ms


def test_on_axis_reduction_pure_Mz_matches_scalar():
    tbz = _theta_break(_McZ, _GfZ)
    n = 1000
    th_c, M_c = _drive_pure(_RZ, 1.1 * tbz, n, lambda: _build_coupled(shape="-linear"))
    coupled_energy = _hc("energy")
    coupled_jumpY = _hc("jumpY")
    th_s, M_s = _drive_pure(_RZ, 1.1 * tbz, n, lambda: _build_blockdiag(shape="-linear"))
    # identical Mz-theta path (the coupled law reduces exactly to the scalar on the pure axis)
    for a, b in zip(M_c, M_s):
        assert abs(a - b) <= 1e-6 * (abs(b) + 1.0), (a, b)
    # coupled hinge dissipated Gf_z and the weak axis never opened
    assert math.isclose(coupled_energy, _GfZ, rel_tol=3e-3), (coupled_energy, _GfZ)
    assert abs(coupled_jumpY) <= 1e-10, coupled_jumpY


def test_on_axis_reduction_pure_My_dissipates_GfY():
    tby = _theta_break(_McY, _GfY)
    n = 1000
    _build_coupled(shape="-linear")
    _prescribe({_RY: 1.1 * tby}, n)
    for _ in range(n):
        assert ops.analyze(1) == 0
    assert math.isclose(_hc("energy"), _GfY, rel_tol=3e-3), (_hc("energy"), _GfY)
    assert abs(_hc("jumpZ")) <= 1e-10, _hc("jumpZ")


# --------------------------------------------------------------------------- #
#  elliptical interaction onset: combined load activates ON the ellipse, earlier
#  than the independent per-axis capacities (the whole point of the coupled law)
# --------------------------------------------------------------------------- #
def test_elliptical_interaction_onset():
    _build_coupled(shape="-linear")
    # Prescribe BOTH rotations proportionally with rz:ry = Mcz:Mcy (a 45deg-normalized path, no
    # dead-load snap): the hinge activates on the ellipse. At onset (r=1, D=0) the cohesive moments
    # satisfy exactly sqrt((Mz/Mcz)^2+(My/Mcy)^2) = 1, EARLIER than either independent capacity
    # (Mz=Mcz/sqrt2, My=Mcy/sqrt2) — the whole point of the coupled interaction surface.
    nstep = 4000
    _prescribe({_RZ: 0.12 * _McZ, _RY: 0.12 * _McY}, nstep, tol=1.0e-11, maxiter=20)
    onset = None
    for _ in range(nstep):
        assert ops.analyze(1) == 0
        if _hc("damage") > 1.0e-6:
            onset = (_hc("momentZ"), _hc("momentY"))
            break
    assert onset is not None, "hinge never activated"
    Mz, My = onset
    ellipse = math.hypot(Mz / _McZ, My / _McY)
    assert math.isclose(ellipse, 1.0, rel_tol=2e-2), (ellipse, Mz, My)
    # strictly inside the independent box (activated before either axis hit its own capacity)
    assert Mz < 0.9 * _McZ and My < 0.9 * _McY, (Mz, My)


# --------------------------------------------------------------------------- #
#  mixed-mode dissipation: combined radial break dissipates Gf between Gf_y, Gf_z
# --------------------------------------------------------------------------- #
def test_mixed_mode_dissipation_between_GfY_and_GfZ():
    _build_coupled(shape="-linear")
    tbz, tby = _theta_break(_McZ, _GfZ), _theta_break(_McY, _GfY)
    n = 1500
    _prescribe({_RZ: 1.3 * tbz, _RY: 1.3 * tby}, n)
    for _ in range(n):
        assert ops.analyze(1) == 0
    e = _hc("energy")
    assert min(_GfY, _GfZ) - 1e-6 <= e <= max(_GfY, _GfZ) + 1e-6, (e, _GfY, _GfZ)
    # both jumps opened
    assert abs(_hc("jumpZ")) > 0.0 and abs(_hc("jumpY")) > 0.0


# --------------------------------------------------------------------------- #
#  coupled tangent: product-of-inertia section, deep biaxial softening converges
# --------------------------------------------------------------------------- #
def test_coupled_section_plus_coupled_hinge_converges_tightly():
    _build_coupled(section="fiber", shape="-exp")
    tbz, tby = _theta_break(_McZ, _GfZ), _theta_break(_McY, _GfY)
    n = 1500
    _prescribe({_RZ: 1.3 * tbz, _RY: 1.3 * tby}, n, tol=1.0e-11, maxiter=12)
    for _ in range(n):
        assert ops.analyze(1) == 0
    assert _hc("energy") > 0.2 * min(_GfY, _GfZ), _hc("energy")
    assert abs(_hc("jumpZ")) > 0.0 and abs(_hc("jumpY")) > 0.0


# --------------------------------------------------------------------------- #
#  finite-rotation invariance under CorotCrdTransf3d
# --------------------------------------------------------------------------- #
def test_corotational_large_rotation_coupled():
    # Strong-axis-dominant skew transverse displacement under CorotCrdTransf3d: the member rotates
    # through finite angle (> 0.5 rad) while the coupled hinge opens on BOTH axes (jz, jy != 0). The
    # condensed biaxial kb passes through the quaternion-mean triad -> the pinned invariant holds for
    # the coupled law too. (Deeply weak-axis-dominant indirect-control paths sit on the unstable
    # elliptical onset where the frozen-mix tangent is at the edge — prescribe rotations there.)
    _build_coupled(shape="-linear", transf="Corotational")
    n = 2500
    _prescribe({_UY: 0.45 * _L, _UZ: 0.10 * _L}, n, tol=1.0e-10, maxiter=40)
    max_rot = 0.0
    for _ in range(n):
        assert ops.analyze(1) == 0
        max_rot = max(max_rot, math.hypot(ops.nodeDisp(2, _RY), ops.nodeDisp(2, _RZ)))
    assert max_rot > 0.5, max_rot
    assert _hc("energy") > 0.0 and abs(_hc("jumpZ")) > 0.0 and abs(_hc("jumpY")) > 0.0


# --------------------------------------------------------------------------- #
#  Benzeggagh-Kenane mode-mix (ADR 34 #2): a nonlinear Gf interpolation makes the
#  damage mix-dependent so the consistent off-radial tangent term is live; the
#  element's coupled condensation must still converge through deep biaxial softening.
# --------------------------------------------------------------------------- #
def test_bk_modemix_coupled_converges_and_reduces_on_axis():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.section("Elastic", 1, _E, _A, _Iz, _Iy, _G, _J)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    # eta = 2 -> nonlinear B-K fracture-energy mix (live off-radial mode-mix tangent)
    ops.nDMaterial("LadrunoCohesiveHingeBiaxial", 79, _McZ, _GfZ, _McY, _GfY, "-exp", "-bk", 2.0)
    ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hingeBiaxial", 79)
    tbz, tby = _theta_break(_McZ, _GfZ), _theta_break(_McY, _GfY)
    n = 1500
    _prescribe({_RZ: 1.3 * tbz, _RY: 1.3 * tby}, n, tol=1.0e-11, maxiter=12)
    for _ in range(n):
        assert ops.analyze(1) == 0
    # combined break dissipates a mixed-mode energy strictly between the two axes
    e = _hc("energy")
    assert min(_GfY, _GfZ) - 1e-6 <= e <= max(_GfY, _GfZ) + 1e-6, (e, _GfY, _GfZ)
    assert abs(_hc("jumpZ")) > 0.0 and abs(_hc("jumpY")) > 0.0


def test_bk_pure_axis_unchanged_by_eta():
    # on a pure axis Gf_mix = Gf_z regardless of eta -> the B-K element matches the
    # default (eta=1) element's Mz path bit-for-bit.
    tbz = _theta_break(_McZ, _GfZ)
    n = 800

    def build(eta):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, _L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.section("Elastic", 1, _E, _A, _Iz, _Iy, _G, _J)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.beamIntegration("Lobatto", 1, 1, 5)
        a = [79, _McZ, _GfZ, _McY, _GfY, "-linear"]
        if eta != 1.0:
            a += ["-bk", eta]
        ops.nDMaterial("LadrunoCohesiveHingeBiaxial", *a)
        ops.element("LadrunoDispBeamColumn", 1, 1, 2, 1, 1, "-hingeBiaxial", 79)

    def run(eta):
        build(eta)
        _prescribe({_RZ: 1.1 * tbz}, n)
        out = []
        for _ in range(n):
            assert ops.analyze(1) == 0
            out.append(ops.eleForce(1, 6))
        return out

    m1, m2 = run(1.0), run(3.0)
    for a, b in zip(m1, m2):
        assert abs(a - b) <= 1e-9 * (abs(b) + 1.0), (a, b)


# --------------------------------------------------------------------------- #
#  DB roundtrip with the coupled hinge open
# --------------------------------------------------------------------------- #
def test_coupled_database_roundtrip():
    from _testbed.roundtrip import database_roundtrip

    tbz, tby = _theta_break(_McZ, _GfZ), _theta_break(_McY, _GfY)

    def build():
        _build_coupled(shape="-linear")
        _prescribe({_RZ: 0.5 * tbz, _RY: 0.5 * tby}, 300)
        for _ in range(300):
            assert ops.analyze(1) == 0

    database_roundtrip(
        build, [2], 6,
        probe_fn=lambda: [
            ops.eleResponse(1, "hingeBiaxial", "momentZ")[0],
            ops.eleResponse(1, "hingeBiaxial", "momentY")[0],
            ops.eleResponse(1, "hingeBiaxial", "jumpZ")[0],
            ops.eleResponse(1, "hingeBiaxial", "jumpY")[0],
            ops.eleResponse(1, "hingeBiaxial", "energy")[0],
            ops.eleResponse(1, "hingeBiaxial", "rMax")[0],
        ],
    )
