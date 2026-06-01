"""LadrunoBrick ``-geom corot`` (v2) — element-level corotational battery.

``corot`` strips the element rigid rotation R = polar(H) (Procrustes best fit of
the nodal cloud), feeds the SMALL deformational strain in the co-rotating frame
to the EXISTING small-strain material library (``ElasticIsotropic`` here, via
``setTrialStrain``), and rotates the core force / stiffness back to global,
adding the geometric (load) stiffness. No finite-strain material is involved —
that is the whole point of corot (see Ladruno_implementation/10_solid_corotational_adr.md).

Gates, in priority order:
  * RIGID ROTATION ⇒ ZERO STRESS / ZERO FORCE — the defining corotational
    property; an exact algebraic identity (R = Q exactly ⇒ u_d = 0), so it must
    hold to polar-convergence precision regardless of rotation magnitude.
  * SMALL DEFORMATION ⇒ corot reduces to ``linear`` (same std/bbar kernel).
  * OBJECTIVITY: a small strain carried through a large rotation gives the same
    force as the unrotated strain, rotated by R.
  * CONSISTENT TANGENT vs finite difference — v2.0 tangent is the rotated
    material stiffness + the symmetric dominant geometric term (G1 + G1ᵀ); the
    polar-Hessian term is deferred, so the FD match is held to a (documented)
    looser tolerance than the finite path while the deformation stays in the
    corot (small-strain) regime.
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# A distorted (positive-Jacobian) reference hex, so R = polar(H) and the spin
# map Γ are genuinely exercised (not a trivial cube).
_NODES = {
    1: (0.00, 0.00, 0.00),
    2: (1.00, 0.10, 0.05),
    3: (1.10, 1.00, 0.00),
    4: (0.05, 0.95, 0.10),
    5: (0.00, 0.05, 1.00),
    6: (1.00, 0.00, 1.05),
    7: (1.05, 1.00, 1.10),
    8: (0.00, 1.00, 0.95),
}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]
_E = 200.0
_NU = 0.3

_WIGGLE = [
    (0.000, 0.000, 0.000), (0.013, -0.009, 0.006), (-0.007, 0.011, -0.010),
    (0.009, 0.005, 0.012), (-0.011, 0.008, -0.006), (0.006, -0.012, 0.009),
    (-0.010, 0.007, 0.011), (0.012, -0.006, -0.008),
]


def _build(geom, formulation="std"):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    ops.element("LadrunoBrick", 1, *_CONN, 1,
                "-formulation", formulation, "-geom", geom)
    return 1


def _impose_and_solve(u, geom, formulation="std"):
    _build(geom, formulation)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag in _NODES:
        base = (tag - 1) * 3
        for d in range(3):
            ops.sp(tag, d + 1, float(u[base + d]))
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-11, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    return ops.analyze(1)


def _resisting_force(u, geom="corot", formulation="std"):
    assert _impose_and_solve(u, geom, formulation) == 0, "imposed-displacement solve failed"
    return np.array(ops.eleForce(1), dtype=float)


def _element_tangent(u, geom="corot", formulation="std"):
    assert _impose_and_solve(u, geom, formulation) == 0, "imposed-displacement solve failed"
    K = ops.eleResponse(1, "stiff")
    assert K and len(K) == 24 * 24, "LadrunoBrick must answer eleResponse('stiff') with 24x24"
    return np.array(K, dtype=float).reshape(24, 24)


def _affine_disp(Fbar, wiggle=False):
    u = np.zeros(24)
    I = np.eye(3)
    for tag, (x, y, z) in _NODES.items():
        X = np.array([x, y, z])
        d = (np.asarray(Fbar) - I) @ X
        if wiggle:
            d = d + np.array(_WIGGLE[tag - 1])
        u[(tag - 1) * 3:(tag - 1) * 3 + 3] = d
    return u


def _rot(thz, thy):
    cz, sz = math.cos(thz), math.sin(thz)
    Rz = np.array([[cz, -sz, 0.0], [sz, cz, 0.0], [0.0, 0.0, 1.0]])
    cy, sy = math.cos(thy), math.sin(thy)
    Ry = np.array([[cy, 0.0, sy], [0.0, 1.0, 0.0], [-sy, 0.0, cy]])
    return Ry @ Rz


# --------------------------------------------------------------------------- #
#  THE defining gate: pure rigid rotation -> zero stress and zero force.      #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("thz,thy", [(0.35, 0.25), (1.20, -0.80), (2.50, 0.40)])
def test_corot_rigid_rotation_is_stress_free(thz, thy):
    R = _rot(thz, thy)
    u = _affine_disp(R)  # x_a = R X_a  => polar(H) = R exactly => u_d = 0
    assert _impose_and_solve(u, "corot") == 0

    s = ops.eleResponse(1, "stresses")
    assert s and len(s) == 48, "LadrunoBrick 'stresses' response missing"
    smax = max(abs(v) for v in s)
    assert smax <= 1.0e-7 * _E, f"rigid rotation induced stress {smax:.3e}"

    f = np.array(ops.eleForce(1), dtype=float)
    assert np.abs(f).max() <= 1.0e-7 * _E, f"rigid rotation induced force {np.abs(f).max():.3e}"


def test_corot_rigid_rotation_bbar_is_stress_free():
    R = _rot(0.9, 0.5)
    u = _affine_disp(R)
    assert _impose_and_solve(u, "corot", "bbar") == 0
    f = np.array(ops.eleForce(1), dtype=float)
    assert np.abs(f).max() <= 1.0e-7 * _E, f"corot+bbar rigid rotation force {np.abs(f).max():.3e}"


# --------------------------------------------------------------------------- #
#  Small deformation: corot reduces to the linear (small-strain) kernel.      #
# --------------------------------------------------------------------------- #
def test_corot_reduces_to_linear_at_small_deformation():
    eps = 1.0e-5
    Fbar = [[1.0 + eps, 0.5 * eps, 0.0],
            [0.5 * eps, 1.0 - 0.3 * eps, 0.0],
            [0.0, 0.0, 1.0 + 0.2 * eps]]
    u = _affine_disp(Fbar, wiggle=False)
    f_cor = _resisting_force(u, "corot")
    f_lin = _resisting_force(u, "linear")
    scale = max(np.abs(f_lin).max(), 1.0e-30)
    assert np.abs(f_cor - f_lin).max() <= 1.0e-3 * scale, (
        "corot did not reduce to linear at small deformation: "
        f"max diff {np.abs(f_cor - f_lin).max():.3e} vs scale {scale:.3e}"
    )


# --------------------------------------------------------------------------- #
#  Objectivity: a small strain carried through a large rotation gives the     #
#  same internal force as the unrotated small strain, rotated by R.           #
# --------------------------------------------------------------------------- #
def test_corot_objectivity_rotated_small_strain():
    # small uniaxial-ish strain (in the reference frame)
    E = np.array([[1.02, 0.0, 0.0],
                  [0.0, 0.99, 0.0],
                  [0.0, 0.0, 1.005]])
    u_ref = _affine_disp(E.tolist(), wiggle=False)
    f_ref = _resisting_force(u_ref, "corot")           # unrotated

    R = _rot(0.7, -0.4)
    u_rot = _affine_disp((R @ E).tolist(), wiggle=False)  # same strain, then rotate
    f_rot = _resisting_force(u_rot, "corot")

    # the rotated force should equal R applied node-by-node to the unrotated one
    f_ref_rot = np.zeros(24)
    for a in range(8):
        f_ref_rot[3*a:3*a+3] = R @ f_ref[3*a:3*a+3]

    scale = max(np.abs(f_ref).max(), 1.0e-30)
    err = np.abs(f_rot - f_ref_rot).max()
    assert err <= 1.0e-4 * scale, (
        f"corot not objective: rotated force != R·(unrotated force), "
        f"max err {err:.3e} vs scale {scale:.3e}"
    )


# --------------------------------------------------------------------------- #
#  Consistent tangent vs finite difference (v2.0: rotated material + dominant  #
#  symmetric geometric term; polar-Hessian deferred -> looser-but-bounded).    #
# --------------------------------------------------------------------------- #
def test_corot_tangent_matches_finite_difference():
    # rotation + a genuine but small strain (corot regime), no wiggle so the
    # state stays small-strain and the deferred polar-Hessian term stays small.
    R = _rot(0.30, -0.20)
    strain = np.array([[1.03, 0.015, -0.010],
                       [0.015, 0.975, 0.012],
                       [-0.010, 0.012, 1.020]])
    u = _affine_disp((R @ strain).tolist(), wiggle=False)

    K = _element_tangent(u, "corot")
    assert np.isfinite(K).all()

    h = 1.0e-6
    Kfd = np.zeros((24, 24))
    for d in range(24):
        up = u.copy(); up[d] += h
        um = u.copy(); um[d] -= h
        Kfd[:, d] = (_resisting_force(up, "corot") - _resisting_force(um, "corot")) / (2.0 * h)

    scale = np.abs(K).max()
    err = np.abs(K - Kfd).max()
    # v2.0 tolerance: the deferred polar-Hessian term ~ O(strain·force) shows
    # up here. In the corot (small-strain) regime it is bounded well below 1e-2
    # relative; tighten in v2.1 when the Hessian term lands.
    assert err <= 1.0e-2 * scale, (
        f"corot tangent vs FD: max abs err {err:.3e} vs scale {scale:.3e}"
    )
    assert np.abs(K - K.T).max() <= 1.0e-6 * scale, "corot tangent not symmetric"
