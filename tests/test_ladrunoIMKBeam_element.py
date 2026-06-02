"""LadrunoIMKBeam (Ladruno concentrated-plasticity IMK beam, classTag 33003) —
Zone-A element battery.

The element is a 3D beam-column with UNCOUPLED moment-rotation rotational hinges
at the ends (elastic interior in series with an IMK uniaxial law per bending
axis). v1 correctness anchors:

  * **elastic-degenerate**: with NO hinge materials the element is a pure elastic
    beam and must reproduce ``elasticBeamColumn`` to ~1e-9 (same A,E,G,Jx,Iy,Iz,
    same CrdTransf) on a 6-dof mixed load. This proves the basic-system assembly
    and the basic->global transform.

  * **single-hinge series**: a strong-axis hinge at end j, driven under rotation
    control, must reproduce an INDEPENDENT python model of "elastic rotational
    spring a = 4EIz/L in series with the IMK backbone M(theta)". This proves the
    exact series F = F_elastic + F_hinge and the 2x2 internal Newton (no
    n-factor). The check is non-tautological: python solves a(theta - thj) =
    M(thj) for the hinge rotation thj and compares BOTH thj and the moment to the
    element's reported hingeRotation / hingeMoment.

Plan: Ladruno_implementation/14_ladruno_imk_beam.md.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# member elastic properties (consistent SI-ish units)
_A = 1.0e-2
_E = 2.0e11
_G = 8.0e10
_Jx = 1.0e-5
_Iy = 1.0e-5
_Iz = 2.0e-5
_L = 3.0


def _new_model():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, _L, 0.0, 0.0)
    # local x = global x; vecxz = global z -> local y = global y
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)


# ---------------------------------------------------------------------------
# Test 1: elastic-degenerate == elasticBeamColumn
# ---------------------------------------------------------------------------
def _solve_cantilever_tip(make_element):
    """Fix node 1, build one element via make_element(), apply a mixed 6-dof
    load at node 2, linear-solve, return the node-2 displacement (6)."""
    _new_model()
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    make_element()

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    # mixed load so the whole 6x6 basic stiffness is exercised
    ops.load(2, 1.0e4, 2.0e4, -1.5e4, 3.0e3, -4.0e3, 5.0e3)

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return [ops.nodeDisp(2, dof) for dof in range(1, 7)]


def test_elastic_degenerate_matches_elasticBeamColumn():
    def make_imk():
        ops.element("LadrunoIMKBeam", 1, 1, 2,
                    _A, _E, _G, _Jx, _Iy, _Iz, 1)  # no -matZ/-matY -> elastic

    def make_elastic():
        ops.element("elasticBeamColumn", 1, 1, 2,
                    _A, _E, _G, _Jx, _Iy, _Iz, 1)

    d_imk = _solve_cantilever_tip(make_imk)
    d_ref = _solve_cantilever_tip(make_elastic)

    for a, b in zip(d_imk, d_ref):
        assert abs(a - b) <= 1.0e-9 * (abs(b) + 1.0e-12), (a, b)

    # sanity: the tip actually moved (load really applied)
    assert any(abs(x) > 1.0e-9 for x in d_ref)


# ---------------------------------------------------------------------------
# Test 2: single strong-axis hinge reproduces the elastic-series + IMK backbone
# ---------------------------------------------------------------------------
# Steel01 monotonic backbone (pure bilinear kinematic): M(theta).
_My = 1000.0     # yield moment
_k0 = 1.0e6      # initial rotational stiffness of the hinge
_bh = 0.02       # strain-hardening ratio
_theta_y = _My / _k0


def _imk_backbone(theta):
    if theta <= _theta_y:
        return _k0 * theta
    return _My + _bh * _k0 * (theta - _theta_y)


def _series_predict(theta_imposed, a):
    """Solve a*(theta - thj) = M(thj) for the hinge rotation thj (monotonic,
    thj >= 0), return (thj, moment)."""
    thj = 0.0
    for _ in range(100):
        if thj <= _theta_y:
            M = _k0 * thj
            dM = _k0
        else:
            M = _My + _bh * _k0 * (thj - _theta_y)
            dM = _bh * _k0
        g = a * (theta_imposed - thj) - M
        if abs(g) <= 1.0e-12 * (a * theta_imposed + 1.0):
            break
        # dg/dthj = -a - dM
        thj -= g / (-a - dM)
    return thj, _imk_backbone(thj)


def test_single_hinge_series_matches_independent_model():
    # Iz chosen so a = 4 E Iz / L is a round number is unnecessary; use _Iz.
    a = 4.0 * _E * _Iz / _L

    _new_model()
    ops.uniaxialMaterial("Steel01", 10, _My, _k0, _bh)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    # node 2: only the strong-axis rotation (dof 6) is free
    ops.fix(2, 1, 1, 1, 1, 1, 0)

    ops.element("LadrunoIMKBeam", 1, 1, 2,
                _A, _E, _G, _Jx, _Iy, _Iz, 1,
                "-hinge", "j", "-matZ", 10)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0)  # reference moment on dof 6

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    dtheta = 1.0e-4
    ops.integrator("DisplacementControl", 2, 6, dtheta)
    ops.test("NormDispIncr", 1.0e-12, 50)
    ops.algorithm("Newton")
    ops.analysis("Static")

    n_steps = 45  # spans elastic (theta_y=1e-3) well into hardening (~4.5e-3)
    for step in range(1, n_steps + 1):
        assert ops.analyze(1) == 0, f"step {step} failed"

        theta_imposed = ops.nodeDisp(2, 6)
        hingeRot = ops.eleResponse(1, "hingeRotation")    # [z_i,z_j,y_i,y_j]
        hingeMom = ops.eleResponse(1, "hingeMoment")      # [Mz_i,Mz_j,My_i,My_j]
        thj_ele = hingeRot[1]
        M_ele = hingeMom[1]

        thj_pred, M_pred = _series_predict(theta_imposed, a)

        assert abs(thj_ele - thj_pred) <= 1.0e-7 * (abs(thj_pred) + 1.0e-9), (
            f"step {step}: hinge rot {thj_ele} vs {thj_pred}")
        assert abs(M_ele - M_pred) <= 1.0e-6 * (abs(M_pred) + 1.0e-6), (
            f"step {step}: hinge moment {M_ele} vs {M_pred}")

    # we must have crossed yield for the test to be meaningful
    assert ops.nodeDisp(2, 6) > _theta_y * 2.0


def test_no_hinge_when_material_omitted_is_pure_elastic():
    """A LadrunoIMKBeam with -hinge but no -matZ/-matY must still be elastic
    (the hinge exists only where a material is supplied)."""
    def make_imk_hinge_no_mat():
        ops.element("LadrunoIMKBeam", 1, 1, 2,
                    _A, _E, _G, _Jx, _Iy, _Iz, 1, "-hinge", "both")

    def make_elastic():
        ops.element("elasticBeamColumn", 1, 1, 2,
                    _A, _E, _G, _Jx, _Iy, _Iz, 1)

    d_imk = _solve_cantilever_tip(make_imk_hinge_no_mat)
    d_ref = _solve_cantilever_tip(make_elastic)
    for a, b in zip(d_imk, d_ref):
        assert abs(a - b) <= 1.0e-9 * (abs(b) + 1.0e-12), (a, b)
