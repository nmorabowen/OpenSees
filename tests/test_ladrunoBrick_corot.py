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


def _make_material(material, tag=1):
    """Create the requested small-strain nDMaterial; corot reuses it unchanged."""
    if material == "elastic":
        ops.nDMaterial("ElasticIsotropic", tag, _E, _NU)
    elif material == "j2":
        # von Mises (J2) with a little linear hardening — a path-dependent
        # material, to prove corot reuses the small-strain library as-is.
        K = _E / (3.0 * (1.0 - 2.0 * _NU))
        G = _E / (2.0 * (1.0 + _NU))
        sigY = 0.30
        ops.nDMaterial("J2Plasticity", tag, K, G, sigY, sigY, 0.0, 2.0)
    else:
        raise ValueError(material)
    return tag


def _build(geom, formulation="std", material="elastic"):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    _make_material(material, 1)
    ops.element("LadrunoBrick", 1, *_CONN, 1,
                "-formulation", formulation, "-geom", geom)
    return 1


def _impose_and_solve(u, geom, formulation="std", material="elastic"):
    _build(geom, formulation, material)
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


def _resisting_force(u, geom="corot", formulation="std", material="elastic"):
    assert _impose_and_solve(u, geom, formulation, material) == 0, "imposed-displacement solve failed"
    return np.array(ops.eleForce(1), dtype=float)


def _element_tangent(u, geom="corot", formulation="std", material="elastic"):
    assert _impose_and_solve(u, geom, formulation, material) == 0, "imposed-displacement solve failed"
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
# includes cases that cross π (179°, 180°, 270°) — the regime where R has an
# eigenvalue near −1 and the polar iteration / spin-map inverse are most stressed
# (ADR §5 Gate 1 asks for 179°; §7 deferred >π tracking — exercise it here).
@pytest.mark.parametrize("thz,thy", [
    (0.35, 0.25), (1.20, -0.80), (2.50, 0.40),
    (math.radians(179.0), 0.0), (math.pi, 0.0), (1.5 * math.pi, 0.0),
])
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
    R = _rot(0.7, -0.4)

    assert _impose_and_solve(_affine_disp(E.tolist()), "corot") == 0
    f_ref = np.array(ops.eleForce(1), dtype=float)                       # unrotated
    s_ref = np.array(ops.eleResponse(1, "stresses"), dtype=float)        # 48 = 8 GP x 6

    assert _impose_and_solve(_affine_disp((R @ E).tolist()), "corot") == 0
    f_rot = np.array(ops.eleForce(1), dtype=float)                       # same strain, rotated
    s_rot = np.array(ops.eleResponse(1, "stresses"), dtype=float)

    # (1) force: rotated == R applied node-by-node to the unrotated force
    f_ref_rot = np.zeros(24)
    for a in range(8):
        f_ref_rot[3*a:3*a+3] = R @ f_ref[3*a:3*a+3]
    scale = max(np.abs(f_ref).max(), 1.0e-30)
    assert np.abs(f_rot - f_ref_rot).max() <= 1.0e-4 * scale, (
        f"corot not objective: rotated force != R·(unrotated force), "
        f"max err {np.abs(f_rot - f_ref_rot).max():.3e} vs scale {scale:.3e}"
    )
    # (2) the corotated GP stress is frame-INVARIANT: the material sees the same
    #     u_d in both cases, so eleResponse('stresses') must be bit-identical to
    #     polar precision (an exact identity — tight tolerance).
    sscale = max(np.abs(s_ref).max(), 1.0e-30)
    assert np.abs(s_rot - s_ref).max() <= 1.0e-6 * sscale, (
        f"corot stress not frame-invariant: rotated GP stress != unrotated "
        f"(max diff {np.abs(s_rot - s_ref).max():.3e})"
    )


# --------------------------------------------------------------------------- #
#  Consistent tangent vs finite difference.                                    #
#  v2.0 omits two O(strain) geometric terms (polar-Hessian ∂Γ/∂u·m AND the     #
#  material-frame ∂R coupling), so the FD-tangent gap scales ~linearly with    #
#  strain. We therefore (a) prove the IMPLEMENTED terms are correct by showing  #
#  the gap VANISHES as strain→0 (a wrong R̄k_dR̄ᵀ or G1+G1ᵀ would not vanish),   #
#  and (b) keep a documented loose bound at finite strain. Symmetry is checked  #
#  on the FD tangent (the in-code symmetriser makes asserting it on K vacuous). #
# --------------------------------------------------------------------------- #
def _fd_tangent(u, geom="corot", formulation="std", h=1.0e-6):
    Kfd = np.zeros((24, 24))
    for d in range(24):
        up = u.copy(); up[d] += h
        um = u.copy(); um[d] -= h
        Kfd[:, d] = (_resisting_force(up, geom, formulation)
                     - _resisting_force(um, geom, formulation)) / (2.0 * h)
    return Kfd


@pytest.mark.parametrize("formulation", ["std", "bbar"])
def test_corot_tangent_gap_vanishes_with_strain(formulation):
    # gap = |K_impl - K_fd| / |K|; must shrink ~linearly as strain -> 0.
    R = _rot(0.30, -0.20)
    base = np.array([[1.0, 0.4, -0.3], [0.4, -0.6, 0.2], [-0.3, 0.2, 0.8]])
    rel = {}
    for eps in (2.0e-3, 1.0e-3, 5.0e-4):
        u = _affine_disp((R @ (np.eye(3) + eps * base)).tolist(), wiggle=False)
        K = _element_tangent(u, "corot", formulation)
        assert np.isfinite(K).all()
        Kfd = _fd_tangent(u, "corot", formulation)
        rel[eps] = np.abs(K - Kfd).max() / max(np.abs(K).max(), 1e-30)
        # symmetry on the FD tangent (bypasses the in-code symmetriser)
        fds = max(np.abs(Kfd).max(), 1e-30)
        assert np.abs(Kfd - Kfd.T).max() <= 1.0e-5 * fds, (
            f"corot FD tangent asymmetric @eps={eps}: geometric term mis-assembled"
        )
    # the deferred terms are O(strain) -> the gap must roughly halve as eps halves
    assert rel[1.0e-3] < 0.75 * rel[2.0e-3], f"gap not vanishing: {rel}"
    assert rel[5.0e-4] < 0.75 * rel[1.0e-3], f"gap not vanishing: {rel}"
    # and be small at the smallest strain (implemented terms correct)
    assert rel[5.0e-4] < 5.0e-4, f"corot small-strain tangent gap too large: {rel}"


def test_corot_tangent_finite_strain_bound():
    # documented deferred-term bound (NOT a tight correctness gate): at ~3% strain
    # the omitted O(strain) geometric terms appear; assert bounded and log.
    R = _rot(0.30, -0.20)
    strain = np.array([[1.03, 0.015, -0.010],
                       [0.015, 0.975, 0.012],
                       [-0.010, 0.012, 1.020]])
    u = _affine_disp((R @ strain).tolist(), wiggle=False)
    K = _element_tangent(u, "corot")
    Kfd = _fd_tangent(u, "corot")
    scale = np.abs(K).max()
    err = np.abs(K - Kfd).max()
    print(f"corot finite-strain FD-tangent rel err = {err/scale:.3e} (deferred O(strain) terms)")
    assert err <= 2.0e-2 * scale, f"corot tangent gap {err/scale:.3e} exceeds deferred-term bound"


# --------------------------------------------------------------------------- #
#  ADR Gate 6 — material reuse: J2 plasticity works in the corotated frame.    #
#  A path-dependent (von Mises) material, driven by setTrialStrain, must yield  #
#  correctly under a large rigid rotation: the material sees the same          #
#  deformational strain u_d whether or not the body is rotated, so its         #
#  corotated-frame response is frame-invariant.                                #
# --------------------------------------------------------------------------- #
def test_corot_j2_plasticity_objectivity():
    # a strain large enough to yield (sigY chosen small in _make_material)
    Ey = np.array([[1.04, 0.0, 0.0],
                   [0.0, 0.98, 0.0],
                   [0.0, 0.0, 1.01]])
    R = _rot(0.8, -0.5)

    assert _impose_and_solve(_affine_disp(Ey.tolist()), "corot", "std", "j2") == 0
    f_ref = np.array(ops.eleForce(1), dtype=float)
    s_ref = np.array(ops.eleResponse(1, "stresses"), dtype=float)
    assert np.isfinite(f_ref).all() and np.isfinite(s_ref).all()

    assert _impose_and_solve(_affine_disp((R @ Ey).tolist()), "corot", "std", "j2") == 0
    f_rot = np.array(ops.eleForce(1), dtype=float)
    s_rot = np.array(ops.eleResponse(1, "stresses"), dtype=float)

    # force rotates by R node-by-node; corotated GP stress is frame-invariant
    f_ref_rot = np.zeros(24)
    for a in range(8):
        f_ref_rot[3*a:3*a+3] = R @ f_ref[3*a:3*a+3]
    fscale = max(np.abs(f_ref).max(), 1.0e-30)
    assert np.abs(f_rot - f_ref_rot).max() <= 1.0e-4 * fscale, (
        "corot+J2 not objective: rotated force != R·(unrotated force)"
    )
    sscale = max(np.abs(s_ref).max(), 1.0e-30)
    assert np.abs(s_rot - s_ref).max() <= 1.0e-6 * sscale, (
        "corot+J2 GP stress not frame-invariant under rigid rotation"
    )
    # sanity: the state actually yielded (deviatoric stress capped near sigY),
    # i.e. J2 is genuinely active in the corotated frame, not silently elastic.
    s0 = s_ref[:6]
    mean = (s0[0] + s0[1] + s0[2]) / 3.0
    dev = np.array([s0[0]-mean, s0[1]-mean, s0[2]-mean, s0[3], s0[4], s0[5]])
    vm = math.sqrt(1.5 * (dev[0]**2 + dev[1]**2 + dev[2]**2
                          + 2.0*(dev[3]**2 + dev[4]**2 + dev[5]**2)))
    assert vm > 0.30, f"corot+J2 did not yield (von Mises {vm:.3f} <= sigY); test is vacuous"


# --------------------------------------------------------------------------- #
#  ADR Gate 5 — load-driven large-displacement cantilever (drives Newton).     #
#  A slender stack of bricks, base fixed, transverse tip load. This is the     #
#  only test that solves for its OWN equilibrium (residual + tangent through a  #
#  real Newton loop) rather than probing prescribed kinematics: convergence    #
#  across load steps catches a wrong-sign / mis-scaled geometric stiffness,     #
#  and agreement with a -geom finite reference catches a wrong force globalize. #
# --------------------------------------------------------------------------- #
def _build_cantilever(geom, nz=6, P=0.0):
    """1x1xnz column of LadrunoBricks, base (z=0) fully fixed, transverse (x)
    tip load P split over the 4 top nodes. corot reuses the small-strain
    ElasticIsotropic directly; the finite reference needs a LogStrain wrapper
    (driven by setTrialF). Returns the top-node tags."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    # nodes: level k (z=k) has 4 nodes  id = 4*k + {1,2,3,4}
    corners = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    for k in range(nz + 1):
        for c, (x, y) in enumerate(corners):
            ops.node(4*k + c + 1, x, y, float(k))
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    if geom == "finite":
        ops.nDMaterial("LogStrain", 2, 1)
        mtag = 2
    else:
        mtag = 1
    for k in range(nz):
        b = 4*k
        conn = [b+1, b+2, b+3, b+4, b+5, b+6, b+7, b+8]
        ops.element("LadrunoBrick", k + 1, *conn, mtag, "-formulation", "std", "-geom", geom)
    for c in range(4):                       # fix the base
        ops.fix(c + 1, 1, 1, 1)
    top = [4*nz + c + 1 for c in range(4)]
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, P / 4.0, 0.0, 0.0)       # transverse (x) tip load
    return top


def _solve_cantilever(geom, nz=6, P=0.6, nsteps=20):
    top = _build_cantilever(geom, nz, P)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-8, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            return None
    return float(np.mean([ops.nodeDisp(t, 1) for t in top]))   # mean tip x-disp


def test_corot_cantilever_converges_and_matches_finite():
    # P=0.6 drives the tip to ~26% of the length (large rotation, small strain) —
    # the corot regime. Newton converges through the load steps despite the
    # deferred tangent terms (the residual is the exact gradient).
    ux_cor = _solve_cantilever("corot")
    assert ux_cor is not None, "corot cantilever Newton failed to converge"
    assert ux_cor > 0.0 and math.isfinite(ux_cor), f"bad corot tip disp {ux_cor}"

    ux_fin = _solve_cantilever("finite")
    assert ux_fin is not None, "finite cantilever reference failed to converge"

    # corot (large-rotation/small-strain) must agree closely with the full
    # finite-strain solution for this slender elastic column.
    rel = abs(ux_cor - ux_fin) / abs(ux_fin)
    print(f"cantilever tip ux: corot={ux_cor:.5e} finite={ux_fin:.5e} rel={rel:.3e}")
    assert rel <= 2.0e-2, (
        f"corot tip disp {ux_cor:.4e} disagrees with finite {ux_fin:.4e} (rel {rel:.3e})"
    )
