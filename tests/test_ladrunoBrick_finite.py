"""LadrunoBrick ``-geom finite`` (v3) — updated-Lagrangian finite-strain battery.

v3 adds the geometry-method axis ``-geom <linear|corot|finite>`` to LadrunoBrick.
``finite`` drives a FiniteStrainNDMaterial (e.g. ``LogStrain``) via ``setTrialF(F)``
and assembles, on the current configuration,

    f = ∫ Bᵀ σ dv ,   K = ∫ Bᵀ c B dv + ∫ Gᵀ Σ G dv

where the material returns the Cauchy stress σ and the spatial CONSTITUTIVE
modulus c, and the ELEMENT owns the geometric/initial-stress stiffness
∫ Gᵀ Σ G dv (the LOCKED seam-3 contract).

The headline gate is the **consistent-tangent finite-difference test**: the
element tangent K must equal d(resisting force)/d(nodal disp) — this is the
arbiter that the push-forward AND the geometric stiffness are correct. On top:

  * pure rigid rotation → zero stress / zero resisting force (objectivity);
  * small strain → ``finite`` reduces to ``linear`` (same std kernel);
  * ``-geom finite`` is refused for a non-finite material and for bbar/uri.

Tangent contract: Ladruno_implementation/09_finite_strain_material_wrapper.md (§4).
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# A deliberately distorted (positive-Jacobian) reference hex, so the F
# push-forward and the geometric stiffness are genuinely exercised.
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

# Small fixed per-node "wiggle" (deterministic, no RNG) to excite bending /
# hourglass modes on top of an affine deformation.
_WIGGLE = [
    (0.000, 0.000, 0.000), (0.013, -0.009, 0.006), (-0.007, 0.011, -0.010),
    (0.009, 0.005, 0.012), (-0.011, 0.008, -0.006), (0.006, -0.012, 0.009),
    (-0.010, 0.007, 0.011), (0.012, -0.006, -0.008),
]


def _build(geom, inner_small_strain=False):
    """Build the single-element model. Returns the element/material tags.

    finite -> LogStrain(ElasticIsotropic);  linear -> ElasticIsotropic directly.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    if geom == "finite":
        ops.nDMaterial("LogStrain", 2, 1)
        mtag = 2
    else:
        mtag = 1
    ops.element("LadrunoBrick", 1, *_CONN, mtag, "-formulation", "std", "-geom", geom)
    return mtag


def _impose_and_solve(u, geom):
    """Impose the full 24-dof nodal displacement vector u via single-point
    constraints and solve one step, leaving the element committed at that state.
    The Lagrange handler keeps the DOFs and adds multiplier equations (so the
    system is non-empty — a fully sp-constrained element gives 0 free DOFs, which
    the Transformation handler cannot solve) while imposing u exactly.
    Returns the analyze() return code."""
    _build(geom)
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


def _resisting_force(u, geom="finite"):
    assert _impose_and_solve(u, geom) == 0, "imposed-displacement solve failed"
    return np.array(ops.eleForce(1), dtype=float)


def _element_tangent(u, geom="finite"):
    assert _impose_and_solve(u, geom) == 0, "imposed-displacement solve failed"
    K = ops.eleResponse(1, "stiff")
    assert K and len(K) == 24 * 24, "LadrunoBrick must answer eleResponse('stiff') with 24x24"
    return np.array(K, dtype=float).reshape(24, 24)


def _affine_disp(Fbar, wiggle=False):
    """Nodal displacements for a uniform deformation gradient Fbar applied to the
    reference coords: u_a = (Fbar - I) X_a (+ optional fixed wiggle)."""
    u = np.zeros(24)
    I = np.eye(3)
    for tag, (x, y, z) in _NODES.items():
        X = np.array([x, y, z])
        d = (np.asarray(Fbar) - I) @ X
        if wiggle:
            d = d + np.array(_WIGGLE[tag - 1])
        u[(tag - 1) * 3:(tag - 1) * 3 + 3] = d
    return u


# --------------------------------------------------------------------------- #
#  The headline gate: consistent tangent == finite difference of the force.   #
# --------------------------------------------------------------------------- #
def test_finite_consistent_tangent_matches_finite_difference():
    # A moderate, genuinely 3D deformation: stretch + shear + per-node wiggle.
    Fbar = [[1.15, 0.08, -0.05],
            [0.04, 0.92, 0.06],
            [-0.03, 0.05, 1.10]]
    u = _affine_disp(Fbar, wiggle=True)

    K = _element_tangent(u, "finite")
    f0 = _resisting_force(u, "finite")
    assert np.isfinite(K).all() and np.isfinite(f0).all()

    h = 1.0e-6
    Kfd = np.zeros((24, 24))
    for d in range(24):
        up = u.copy(); up[d] += h
        um = u.copy(); um[d] -= h
        Kfd[:, d] = (_resisting_force(up, "finite") - _resisting_force(um, "finite")) / (2.0 * h)

    scale = np.abs(K).max()
    err = np.abs(K - Kfd).max()
    # central FD truncation ~ O(h^2) plus the element's own conditioning; the
    # consistent tangent (material + geometric) must match well within 1e-4.
    assert err <= 1.0e-4 * scale, (
        f"finite tangent != FD of resisting force: max abs err {err:.3e} "
        f"vs scale {scale:.3e} (geometric stiffness / push-forward bug?)"
    )
    # symmetry (elastic inner + conservative geometric stiffness)
    assert np.abs(K - K.T).max() <= 1.0e-6 * scale, "finite tangent not symmetric"


# --------------------------------------------------------------------------- #
#  Objectivity: a pure rigid rotation produces no stress and no force.        #
# --------------------------------------------------------------------------- #
def test_finite_rigid_rotation_is_stress_free():
    th = 0.35  # rad, about the z axis then tilted — a finite rotation
    cz, sz = math.cos(th), math.sin(th)
    Rz = np.array([[cz, -sz, 0.0], [sz, cz, 0.0], [0.0, 0.0, 1.0]])
    cy, sy = math.cos(0.25), math.sin(0.25)
    Ry = np.array([[cy, 0.0, sy], [0.0, 1.0, 0.0], [-sy, 0.0, cy]])
    R = Ry @ Rz
    u = _affine_disp(R)  # pure rotation => F = R at every GP, no stretch

    assert _impose_and_solve(u, "finite") == 0
    # element-level 'stresses' response (8 GP x 6 Cauchy components = 48)
    s = ops.eleResponse(1, "stresses")
    assert s and len(s) == 48, "LadrunoBrick 'stresses' response missing"
    smax = max(abs(v) for v in s)
    assert smax <= 1.0e-7 * _E, f"rigid rotation induced stress {smax:.3e} (not objective)"

    f = np.array(ops.eleForce(1), dtype=float)
    assert np.abs(f).max() <= 1.0e-7 * _E, f"rigid rotation induced force {np.abs(f).max():.3e}"


# --------------------------------------------------------------------------- #
#  Small strain: finite reduces to the linear (small-strain) std kernel.      #
# --------------------------------------------------------------------------- #
def test_finite_reduces_to_linear_at_small_strain():
    eps = 1.0e-5
    Fbar = [[1.0 + eps, 0.5 * eps, 0.0],
            [0.5 * eps, 1.0 - 0.3 * eps, 0.0],
            [0.0, 0.0, 1.0 + 0.2 * eps]]
    # NO wiggle here: the wiggle (~1e-2) would dominate eps (~1e-5) and the state
    # would not be small-strain. A pure tiny affine deformation is.
    u = _affine_disp(Fbar, wiggle=False)

    f_fin = _resisting_force(u, "finite")
    f_lin = _resisting_force(u, "linear")

    scale = max(np.abs(f_lin).max(), 1.0e-30)
    # difference is O(strain) relative; at eps=1e-5 a few x1e-4 relative is ample.
    assert np.abs(f_fin - f_lin).max() <= 1.0e-3 * scale, (
        "finite did not reduce to linear at small strain: "
        f"max diff {np.abs(f_fin - f_lin).max():.3e} vs scale {scale:.3e}"
    )


# --------------------------------------------------------------------------- #
#  Homogeneous-deformation patch test: constant F at every GP => constant      #
#  Cauchy stress, matching the closed-form Hencky-elastic value.               #
# --------------------------------------------------------------------------- #
def test_finite_homogeneous_deformation_patch():
    import logstrain_reference as o  # the numpy oracle (tests/)
    Fm = np.array([[1.10, 0.05, -0.02],
                   [0.03, 0.96, 0.04],
                   [-0.01, 0.02, 1.07]])
    u = _affine_disp(Fm.tolist(), wiggle=False)   # affine => F = Fm at every GP
    assert _impose_and_solve(u, "finite") == 0

    s = np.array(ops.eleResponse(1, "stresses"), dtype=float).reshape(8, 6)
    s0 = s[0]
    # (1) constant stress across all 8 Gauss points (no spurious gradients)
    assert np.abs(s - s0).max() <= 1.0e-7 * max(np.abs(s0).max(), 1.0), (
        f"finite patch: Cauchy stress not constant across GPs (max dev "
        f"{np.abs(s - s0).max():.3e})"
    )
    # (2) matches the closed-form Hencky-elastic Cauchy stress for Fm
    Kbulk = _E / (3.0 * (1.0 - 2.0 * _NU))
    G = _E / (2.0 * (1.0 + _NU))
    sig, tau, J, eps = o.cauchy_from_F_elastic(Fm, Kbulk, G)
    sig_v = np.array([sig[0, 0], sig[1, 1], sig[2, 2],
                      sig[0, 1], sig[1, 2], sig[2, 0]])
    tol = 1.0e-6 * max(np.abs(sig_v).max(), 1.0)
    assert np.allclose(s0, sig_v, rtol=1.0e-6, atol=tol), (
        f"finite patch: GP Cauchy {s0} != analytic {sig_v}"
    )


# --------------------------------------------------------------------------- #
#  API guards.                                                                #
# --------------------------------------------------------------------------- #
def _ele_absent_after(make):
    """Run a (rejected) element() call and assert no element was created,
    whether openseespy raises or just returns on the factory's null."""
    try:
        make()
    except Exception:
        pass
    tags = ops.getEleTags()
    if tags is None:
        return
    if isinstance(tags, int):
        tags = [tags]
    assert 1 not in tags


def test_finite_requires_a_finite_strain_material():
    # -geom finite with a plain small-strain material must be refused.
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    _ele_absent_after(
        lambda: ops.element("LadrunoBrick", 1, *_CONN, 1, "-formulation", "std", "-geom", "finite")
    )


def test_finite_only_with_std_formulation():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    ops.nDMaterial("LogStrain", 2, 1)
    _ele_absent_after(
        lambda: ops.element("LadrunoBrick", 1, *_CONN, 2, "-formulation", "bbar", "-geom", "finite")
    )
