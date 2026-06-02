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


def _build(geom, form="std", nu=_NU):
    """Build the single-element model. Returns the element/material tags.

    finite -> LogStrain(ElasticIsotropic);  linear -> ElasticIsotropic directly.
    form selects the formulation axis ("std" or, for finite, "bbar" = F-bar).
    nu lets a test exercise the near-incompressible regime (where F-bar matters).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, nu)
    if geom == "finite":
        ops.nDMaterial("LogStrain", 2, 1)
        mtag = 2
    else:
        mtag = 1
    ops.element("LadrunoBrick", 1, *_CONN, mtag, "-formulation", form, "-geom", geom)
    return mtag


def _impose_and_solve(u, geom, form="std", nu=_NU):
    """Impose the full 24-dof nodal displacement vector u via single-point
    constraints and solve one step, leaving the element committed at that state.
    The Lagrange handler keeps the DOFs and adds multiplier equations (so the
    system is non-empty — a fully sp-constrained element gives 0 free DOFs, which
    the Transformation handler cannot solve) while imposing u exactly.
    Returns the analyze() return code."""
    _build(geom, form, nu)
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


def _resisting_force(u, geom="finite", form="std", nu=_NU):
    assert _impose_and_solve(u, geom, form, nu) == 0, "imposed-displacement solve failed"
    return np.array(ops.eleForce(1), dtype=float)


def _element_tangent(u, geom="finite", form="std", nu=_NU):
    assert _impose_and_solve(u, geom, form, nu) == 0, "imposed-displacement solve failed"
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


def test_finite_accepts_std_and_bbar_rejects_uri_eas():
    # std (plain F) and bbar (F-bar) are both valid with -geom finite; uri/eas
    # finite are reserved and must be refused.
    def _mk(form):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        for tag, (x, y, z) in _NODES.items():
            ops.node(tag, x, y, z)
        ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
        ops.nDMaterial("LogStrain", 2, 1)
        ops.element("LadrunoBrick", 1, *_CONN, 2, "-formulation", form, "-geom", "finite")

    for form in ("std", "bbar"):
        _mk(form)
        tags = ops.getEleTags()
        tags = [tags] if isinstance(tags, int) else (tags or [])
        assert 1 in tags, f"-geom finite -formulation {form} should be accepted"

    for form in ("uri", "eas"):
        _ele_absent_after(lambda f=form: _mk(f))


# --------------------------------------------------------------------------- #
#  F-bar (bbar + finite): the headline gate — consistent tangent == FD of the  #
#  resisting force. F-bar's tangent is GENERALLY UNSYMMETRIC (dSNPO eq 15.10),  #
#  so this test does NOT assert symmetry; the FD match is the sole arbiter of   #
#  the eq 15.11 coupling coefficient q = (1/3) a:(I⊗I) − (2/3) σ⊗I.             #
# --------------------------------------------------------------------------- #
def test_fbar_consistent_tangent_matches_finite_difference():
    # Non-affine state (wiggle) so F varies across GPs => G0 != G and the F-bar
    # coupling term is genuinely exercised (an affine state would zero it out).
    Fbar = [[1.15, 0.08, -0.05],
            [0.04, 0.92, 0.06],
            [-0.03, 0.05, 1.10]]
    u = _affine_disp(Fbar, wiggle=True)

    K = _element_tangent(u, "finite", "bbar")
    f0 = _resisting_force(u, "finite", "bbar")
    assert np.isfinite(K).all() and np.isfinite(f0).all()

    h = 1.0e-6
    Kfd = np.zeros((24, 24))
    for d in range(24):
        up = u.copy(); up[d] += h
        um = u.copy(); um[d] -= h
        Kfd[:, d] = (_resisting_force(up, "finite", "bbar")
                     - _resisting_force(um, "finite", "bbar")) / (2.0 * h)

    scale = np.abs(K).max()
    err = np.abs(K - Kfd).max()
    assert err <= 1.0e-4 * scale, (
        f"F-bar tangent != FD of resisting force: max abs err {err:.3e} vs scale "
        f"{scale:.3e} (eq 15.11 coupling q_ij = (1/3)a_ijpp - (2/3)sigma_ij?)"
    )
    # F-bar carries a genuinely unsymmetric coupling on this non-affine state:
    # confirm the test is actually exercising it (else it proves nothing).
    asym = np.abs(K - K.T).max()
    assert asym > 1.0e-6 * scale, (
        "F-bar tangent came out symmetric on a non-affine state — the eq 15.10 "
        "coupling term is not being exercised (G0 == G?)"
    )


def test_fbar_consistent_tangent_matches_fd_near_incompressible():
    # Defense-in-depth on the eq 15.11 coupling q = (1/3)a:(I⊗I) − (2/3)σ⊗I. The
    # easily-dropped −(2/3)σ⊗I term feeds the volumetric block (1/3)a:(I⊗I), whose
    # magnitude blows up as ν→1/2 (bulk modulus ~ E/(3(1−2ν))). Re-running the FD
    # gate at ν=0.499 makes it directly guard that term where it matters most: a
    # regression to the "(1/3)c" spatial-shortcut coefficient would be glaring here.
    nu = 0.499
    Fbar = [[1.15, 0.08, -0.05],
            [0.04, 0.92, 0.06],
            [-0.03, 0.05, 1.10]]
    u = _affine_disp(Fbar, wiggle=True)

    K = _element_tangent(u, "finite", "bbar", nu)
    assert np.isfinite(K).all()

    h = 1.0e-6
    Kfd = np.zeros((24, 24))
    for d in range(24):
        up = u.copy(); up[d] += h
        um = u.copy(); um[d] -= h
        Kfd[:, d] = (_resisting_force(up, "finite", "bbar", nu)
                     - _resisting_force(um, "finite", "bbar", nu)) / (2.0 * h)

    scale = np.abs(K).max()
    err = np.abs(K - Kfd).max()
    assert err <= 1.0e-4 * scale, (
        f"near-incompressible F-bar tangent != FD: max abs err {err:.3e} vs scale "
        f"{scale:.3e} (the −(2/3)σ⊗I term in eq 15.11 is most consequential at ν→1/2)"
    )


# --------------------------------------------------------------------------- #
#  F-bar RESIDUAL matches std finite under HOMOGENEOUS deformation: F constant  #
#  => J0 == J => Fbar = (J0/J)^(1/3) F == F at every GP => sigma_bar == sigma,  #
#  and the residual carries no coupling term, so f_bar == f_std exactly.        #
#                                                                               #
#  NOTE the TANGENT does NOT reduce to std here. F-bar uses the element-CENTROID #
#  dilatation gradient G0 (dSNPO eq 15.4), not the volume-average B-bar gradient #
#  — so (G0 - g) != 0 per Gauss point even for an affine state on a distorted    #
#  hex, and the eq 15.10 coupling term is genuinely nonzero. F-bar is a distinct #
#  element (different volumetric stiffness), not std-with-F-replaced; its        #
#  tangent is validated by the FD-tangent gate above, not by equality with std. #
# --------------------------------------------------------------------------- #
def test_fbar_residual_matches_std_on_homogeneous_deformation():
    Fm = [[1.10, 0.05, -0.02],
          [0.03, 0.96, 0.04],
          [-0.01, 0.02, 1.07]]
    u = _affine_disp(Fm, wiggle=False)            # affine => homogeneous F

    f_std = _resisting_force(u, "finite", "std")
    f_bar = _resisting_force(u, "finite", "bbar")
    scale = max(np.abs(f_std).max(), 1.0e-30)
    assert np.abs(f_bar - f_std).max() <= 1.0e-9 * scale, (
        "F-bar residual did not reduce to std finite on a homogeneous deformation "
        f"(Fbar should equal F): max force diff {np.abs(f_bar - f_std).max():.3e}"
    )

    # The F-bar tangent legitimately DIFFERS from std even here (centroid vs
    # local dilatation): confirm the coupling is actually present (not a no-op).
    K_std = _element_tangent(u, "finite", "std")
    K_bar = _element_tangent(u, "finite", "bbar")
    ks = np.abs(K_std).max()
    assert np.abs(K_bar - K_std).max() > 1.0e-6 * ks, (
        "F-bar tangent is identical to std on a distorted hex — the eq 15.10 "
        "centroid coupling term appears to be missing"
    )


# --------------------------------------------------------------------------- #
#  T4: F-bar relieves volumetric locking. A slender cantilever (one element     #
#  through the thickness) bent under a tip shear locks for the standard hex as  #
#  nu -> 1/2; F-bar stays compliant. Reference-free signature: compare the      #
#  tip deflection of bbar vs std at nu=0.499 (F-bar must be much softer) and at  #
#  nu=0.30 (both close, F-bar barely changes the compressible answer).          #
# --------------------------------------------------------------------------- #
def _cantilever_tip_disp(form, nu, nz=6, tip_force=0.02):
    """Build a 1x1xnz column of LadrunoBricks (LogStrain/ElasticIsotropic at the
    given nu), fix the base, apply a transverse (x) tip force split over the 4 top
    nodes, solve, and return the mean tip x-displacement."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    # node grid: (i,j) in {0,1}^2 across the section, k=0..nz up the column.
    def nid(i, j, k):
        return 1 + i + 2 * j + 4 * k

    for k in range(nz + 1):
        for j in range(2):
            for i in range(2):
                ops.node(nid(i, j, k), float(i), float(j), float(k))
    # fix the base layer (k=0) fully
    for j in range(2):
        for i in range(2):
            ops.fix(nid(i, j, 0), 1, 1, 1)

    ops.nDMaterial("ElasticIsotropic", 1, _E, nu)
    ops.nDMaterial("LogStrain", 2, 1)
    for k in range(nz):
        # bottom face (k) nodes 1-4, top face (k+1) nodes 5-8, hex node order
        conn = [nid(0, 0, k), nid(1, 0, k), nid(1, 1, k), nid(0, 1, k),
                nid(0, 0, k + 1), nid(1, 0, k + 1), nid(1, 1, k + 1), nid(0, 1, k + 1)]
        ops.element("LadrunoBrick", k + 1, *conn, 2, "-formulation", form, "-geom", "finite")

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(2):
        for i in range(2):
            ops.load(nid(i, j, nz), tip_force / 4.0, 0.0, 0.0)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")          # F-bar tangent is unsymmetric
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.25)
    ops.analysis("Static")
    assert ops.analyze(4) == 0, f"cantilever solve failed (form={form}, nu={nu})"

    return float(np.mean([ops.nodeDisp(nid(i, j, nz), 1)
                          for j in range(2) for i in range(2)]))


def test_fbar_relieves_volumetric_locking():
    # Tip deflection at a moderate (nu=0.30) and a near-incompressible (nu=0.499)
    # Poisson ratio, for std (plain F) and bbar (F-bar). The volumetric-locking
    # signature is reference-free and lives in the CROSS-nu behaviour:
    #   - std stiffens sharply (deflection collapses) as nu -> 1/2  => locks;
    #   - F-bar deflection stays essentially constant                => no lock.
    # (A 1-element-through-thickness bending mesh already locks partially at
    #  nu=0.30, so F-bar is somewhat softer there too — that is expected and not
    #  asserted; the discriminating evidence is the nu-trend below.)
    d_std_c = _cantilever_tip_disp("std", 0.30)
    d_bar_c = _cantilever_tip_disp("bbar", 0.30)
    d_std_i = _cantilever_tip_disp("std", 0.499)
    d_bar_i = _cantilever_tip_disp("bbar", 0.499)
    assert min(d_std_c, d_bar_c, d_std_i, d_bar_i) > 0.0

    # (1) standard hex LOCKS: tip deflection collapses going nu 0.30 -> 0.499.
    std_lock = d_std_i / d_std_c
    assert std_lock < 0.5, (
        f"standard hex did not show volumetric locking: nu-ratio {std_lock:.3f} "
        "(expected << 1 as the deflection collapses toward incompressibility)"
    )
    # (2) F-bar stays compliant: deflection barely changes with nu.
    bar_lock = d_bar_i / d_bar_c
    assert 0.7 < bar_lock < 1.5, (
        f"F-bar deflection should be ~nu-insensitive; nu-ratio {bar_lock:.3f}"
    )
    # (3) and so at nu=0.499 F-bar is dramatically softer than the locked std hex.
    unlock = d_bar_i / d_std_i
    assert unlock > 3.0, (
        f"F-bar did not relieve volumetric locking: bbar/std tip-deflection at "
        f"nu=0.499 is only {unlock:.3f} (expected >> 1)"
    )
