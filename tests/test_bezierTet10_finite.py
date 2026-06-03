"""BezierTet10 ``-geom finite`` — updated-Lagrangian finite-strain battery (3D tet).

Tet sibling of test_ladrunoBrick_finite.py, std only. ``finite`` drives a
FiniteStrainNDMaterial (``LogStrain``) via ``setTrialF(F)`` and assembles, on the
current configuration,

    f = ∫ σ_ij ∂Nₐ/∂x_j dv ,
    K = ∫ (∂Nₐ/∂x_j) a_ijkl (∂N_b/∂x_l) dv ,   a_ijkl = c_ijkl − σ_il δ_jk

where the material returns the Cauchy stress σ and the FULL 4th-order spatial
modulus c (getSpatialTangentTensor — the 6×6 getTangent is lossy in (k,l)), and
the −σδ geometric/initial-stress term is the element's. F is built per GP from the
REFERENCE Bernstein gradients (F = I + Σ uₐ⊗∂Nₐ/∂X); the spatial gradients use F⁻¹.

The headline gate is the **consistent-tangent finite-difference test**: K must
equal d(resisting force)/d(nodal disp) — the arbiter that the push-forward AND the
geometric stiffness are correct. On top:

  * pure rigid rotation → zero stress / zero force (Hencky/LogStrain objectivity);
  * small strain → ``finite`` reduces to ``linear`` (same kernel at vanishing F−I);
  * homogeneous F → constant Cauchy at every GP, matching the closed-form oracle;
  * API guards: finite requires a finite-strain material; a finite-strain material
    requires ``-geom finite``; ``-bbar``/``-pressure`` + finite are refused (F-bar
    is a step-2 follow-up).

This is STD only — near-incompressible un-locking is NOT asserted (std finite
locks volumetrically; F-bar is the step-2 cure, tested separately when it lands).

Tangent contract: Ladruno_implementation/09_finite_strain_material_wrapper.md (§4).
"""
import math

import numpy as np
import pytest

from _testbed import ops
from _testbed.roundtrip import database_roundtrip

pytestmark = [pytest.mark.zone_a]

_E = 200.0
_NU = 0.3

# mid-edge (vertexA, vertexB) order for nodes 5..10 (TenNodeTetrahedron / element)
_EDGE_V = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]

# A distorted (positive-volume) reference tet, so the F push-forward and the
# geometric stiffness are genuinely exercised (not a regular simplex).
_VERTS = {
    1: (0.00, 0.00, 0.00),
    2: (1.00, 0.12, 0.05),
    3: (0.18, 1.00, 0.08),
    4: (0.07, 0.15, 1.00),
}

# small per-node wiggle (10 nodes) to break any accidental symmetry in FD tests
_WIGGLE = [
    (0.000, 0.000, 0.000), (0.011, -0.008, 0.006), (-0.007, 0.010, -0.009),
    (0.008, 0.005, 0.011), (-0.010, 0.007, -0.006), (0.006, -0.011, 0.008),
    (-0.009, 0.006, 0.010), (0.012, -0.005, -0.007), (-0.006, 0.009, 0.011),
    (0.010, -0.007, -0.005),
]


def _tet_nodes():
    """Full 10-node coordinate map: 4 vertices + 6 straight-edge midpoints."""
    coord = dict(_VERTS)
    for e, (a, b) in enumerate(_EDGE_V):
        va, vb = _VERTS[a + 1], _VERTS[b + 1]
        coord[5 + e] = tuple(0.5 * (va[i] + vb[i]) for i in range(3))
    return coord


_COORD = _tet_nodes()


def _build(geom, nu=_NU):
    """Build the single-element model.

    finite -> LogStrain(ElasticIsotropic);  linear -> ElasticIsotropic directly.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _COORD.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, nu)
    if geom == "finite":
        ops.nDMaterial("LogStrain", 2, 1)
        mtag = 2
    else:
        mtag = 1
    args = ["BezierTet10", 1, *range(1, 11), mtag]
    if geom is not None:
        args += ["-geom", geom]
    ops.element(*args)
    return mtag


def _impose_and_solve(u, geom, nu=_NU):
    """Impose the full 30-dof nodal displacement vector u via single-point
    constraints and solve one step, leaving the element committed at that state.
    Lagrange keeps the DOFs (a fully sp-constrained element has 0 free DOFs, which
    the Transformation handler cannot solve) while imposing u exactly."""
    _build(geom, nu)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag in _COORD:
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


def _resisting_force(u, geom="finite", nu=_NU):
    assert _impose_and_solve(u, geom, nu) == 0, "imposed-displacement solve failed"
    return np.array(ops.eleForce(1), dtype=float)


def _element_tangent(u, geom="finite", nu=_NU):
    assert _impose_and_solve(u, geom, nu) == 0, "imposed-displacement solve failed"
    K = ops.eleResponse(1, "stiff")
    assert K and len(K) == 30 * 30, "BezierTet10 must answer eleResponse('stiff') with 30x30"
    return np.array(K, dtype=float).reshape(30, 30)


def _affine_disp(Fbar, wiggle=False):
    """u_a = (F - I) X_a at every node. The field is AFFINE, so the Bernstein
    control-point coefficient at each node equals the field value there (the
    degree-2 blossom of a linear field) — true for vertices AND mid-edges, so an
    affine F is reproduced EXACTLY at every Gauss point."""
    u = np.zeros(30)
    I = np.eye(3)
    F = np.asarray(Fbar)
    for tag, (x, y, z) in _COORD.items():
        X = np.array([x, y, z])
        d = (F - I) @ X
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
    Kfd = np.zeros((30, 30))
    for d in range(30):
        up = u.copy(); up[d] += h
        um = u.copy(); um[d] -= h
        Kfd[:, d] = (_resisting_force(up, "finite") - _resisting_force(um, "finite")) / (2.0 * h)

    scale = np.abs(K).max()
    err = np.abs(K - Kfd).max()
    # central FD truncation ~ O(h^2) plus conditioning; the consistent tangent
    # (full a_ijkl = c − σδ via getSpatialTangentTensor) must match within 1e-4.
    assert err <= 1.0e-4 * scale, (
        f"finite tangent != FD of resisting force: max abs err {err:.3e} "
        f"vs scale {scale:.3e} (geometric stiffness / F⁻¹ push-forward bug?)"
    )
    # std finite tangent is symmetric (a_ijkl = c − σδ is symmetric under
    # (a,i)<->(b,k) for symmetric σ); the lossy 6×6 getTangent path would not be.
    assert np.abs(K - K.T).max() <= 1.0e-6 * scale, "finite tangent not symmetric"


# --------------------------------------------------------------------------- #
#  Objectivity: a pure rigid rotation produces no stress and no force.        #
#  Includes cases crossing π/2 and π, where the log-strain map is most         #
#  stressed.                                                                   #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("thz,thy", [
    (0.35, 0.25), (1.20, -0.80), (2.50, 0.40),
    (math.radians(179.0), 0.0), (math.pi, 0.0), (1.5 * math.pi, 0.0),
])
def test_finite_rigid_rotation_is_stress_free(thz, thy):
    R = _rot(thz, thy)
    u = _affine_disp(R)  # x_a = R X_a  =>  F = R exactly => zero log-strain
    assert _impose_and_solve(u, "finite") == 0

    s = ops.eleResponse(1, "stresses")
    assert s and len(s) == 24, "BezierTet10 'stresses' response missing (4 GP x 6)"
    smax = max(abs(v) for v in s)
    assert smax <= 1.0e-7 * _E, f"rigid rotation induced stress {smax:.3e} (not objective)"

    f = np.array(ops.eleForce(1), dtype=float)
    assert np.abs(f).max() <= 1.0e-7 * _E, f"rigid rotation induced force {np.abs(f).max():.3e}"


# --------------------------------------------------------------------------- #
#  Small strain: finite reduces to the linear (small-strain) kernel.          #
# --------------------------------------------------------------------------- #
def test_finite_reduces_to_linear_at_small_strain():
    eps = 1.0e-5
    Fbar = [[1.0 + eps, 0.5 * eps, 0.0],
            [0.5 * eps, 1.0 - 0.3 * eps, 0.0],
            [0.0, 0.0, 1.0 + 0.2 * eps]]
    # NO wiggle: the wiggle (~1e-2) would dominate eps (~1e-5) and the state would
    # not be small-strain. A pure tiny affine deformation is.
    u = _affine_disp(Fbar, wiggle=False)

    f_fin = _resisting_force(u, "finite")
    f_lin = _resisting_force(u, "linear")

    scale = max(np.abs(f_lin).max(), 1.0e-30)
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

    s = np.array(ops.eleResponse(1, "stresses"), dtype=float).reshape(4, 6)
    s0 = s[0]
    # (1) constant stress across all 4 Gauss points (no spurious gradients)
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


def _nodes_only(nu=_NU):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _COORD.items():
        ops.node(tag, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, _E, nu)
    ops.nDMaterial("LogStrain", 2, 1)


def test_finite_requires_a_finite_strain_material():
    # -geom finite with a plain small-strain material must be refused.
    _nodes_only()
    _ele_absent_after(
        lambda: ops.element("BezierTet10", 1, *range(1, 11), 1, "-geom", "finite")
    )


def test_finite_material_requires_geom_finite():
    # Converse: a finite-strain material (LogStrain) under -geom linear|corot
    # must be refused (it disables setTrialStrain -> zero-stress phantom).
    for geom in ("linear", "corot"):
        _nodes_only()
        _ele_absent_after(
            lambda g=geom: ops.element("BezierTet10", 1, *range(1, 11), 2, "-geom", g)
        )


def test_finite_accepts_std_and_bbar():
    # Both std (plain F) and bbar (F-bar, shipped step 2) are valid under
    # -geom finite. (Detailed F-bar coverage lives in test_bezierTet10_fbar.py.)
    for extra in ([], ["-bbar"]):
        _nodes_only()
        ops.element("BezierTet10", 1, *range(1, 11), 2, *extra, "-geom", "finite")
        tags = ops.getEleTags()
        tags = [tags] if isinstance(tags, int) else (tags or [])
        assert 1 in tags, f"-geom finite {extra} should be accepted"


def test_finite_rejects_pressure():
    # The +z pressure hack is unvalidated under finite and must be refused.
    _nodes_only()
    _ele_absent_after(
        lambda: ops.element("BezierTet10", 1, *range(1, 11), 2,
                            "-pressure", 1.0, "-geom", "finite")
    )


# --------------------------------------------------------------------------- #
#  Robustness: det F <= 0 is guarded — the step is cut, not silently wrong.    #
# --------------------------------------------------------------------------- #
def test_finite_detF_guard_propagates():
    # A reflection-like affine field gives det F < 0 at every GP; updateFinite
    # must return failure so the analysis step-cuts (analyze != 0), rather than
    # assembling a negative-volume contribution.
    Fbad = [[-0.5, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0]]          # det = -0.5 < 0
    u = _affine_disp(Fbad, wiggle=False)
    assert _impose_and_solve(u, "finite") != 0, (
        "det F <= 0 was not caught — the element assembled a negative-volume state"
    )

    # positive control: a healthy J > 0 affine field converges.
    Fok = [[1.05, 0.02, 0.0],
           [0.0, 1.03, 0.0],
           [0.0, 0.0, 1.04]]
    assert _impose_and_solve(_affine_disp(Fok, wiggle=False), "finite") == 0


# --------------------------------------------------------------------------- #
#  Serialization: a -geom finite element must round-trip through the           #
#  FE_Datastore (broker reconstruction via classTag 33001 + iData slot 15      #
#  carrying geomMethodID). nodeDisp alone is VACUOUS here — it is pinned by the #
#  sp-constraints and survives even a recvSelf that silently reverted to        #
#  linear — so we ALSO probe ops.eleResponse(1,'stiff'): the finite tangent at  #
#  the committed (stretched+sheared) state differs from the small-strain K, so  #
#  matching it before save vs after restore proves the geometry-method +        #
#  committed F genuinely survived recvSelf.                                     #
# --------------------------------------------------------------------------- #
def test_finite_serialization_roundtrip():
    Fbar = [[1.10, 0.04, -0.03],
            [0.02, 0.95, 0.05],
            [-0.02, 0.03, 1.08]]
    u = _affine_disp(Fbar, wiggle=True)

    def build():
        assert _impose_and_solve(u, "finite") == 0, "roundtrip build failed to converge"

    build()
    probe = list(_COORD)
    dmax = max(abs(ops.nodeDisp(t, 1)) for t in probe)
    assert dmax > 1.0e-6, f"roundtrip build produced ~zero state ({dmax:.2e}); test would be vacuous"
    # the element-owned 'stiff' is the part that actually depends on recvSelf
    # reconstructing the FINITE element (a reverted linear element gives a
    # different tangent at this deformed state).
    database_roundtrip(build, probe, ndf=3,
                       probe_fn=lambda: ops.eleResponse(1, "stiff"))


# --------------------------------------------------------------------------- #
#  Load-driven, multi-step Newton: the geometric stiffness must work as a       #
#  genuine TANGENT (not just an FD-matched matrix at one frozen state). A       #
#  slender cantilever of stacked tet-filled cubes, fixed at the base, is pushed #
#  transversely over several LoadControl increments; convergence at every step  #
#  exercises K (incl. the −σδ geometric term) through the solver. Reference-    #
#  free: assert convergence + a finite, positive, monotonically growing tip     #
#  deflection (no finite oracle exists for Tet10 — see corot test).             #
# --------------------------------------------------------------------------- #
def _kuhn_tets(n0, n1, n2, n3, n4, n5, n6, n7):
    """Split a hex (corner nodes in std hex order) into 6 tets (Kuhn/diagonal
    triangulation through node n0..n6). Returns 6 vertex-quadruples (1-based
    vertex ids in TenNodeTetrahedron v1..v4 convention)."""
    return [
        (n0, n1, n2, n6), (n0, n2, n3, n6),
        (n0, n3, n7, n6), (n0, n7, n4, n6),
        (n0, n4, n5, n6), (n0, n5, n1, n6),
    ]


def test_finite_load_driven_cantilever_converges():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    nz = 4
    side = 0.5

    def vid(i, j, k):
        return 1 + i + 2 * j + 4 * k

    for k in range(nz + 1):
        for j in range(2):
            for i in range(2):
                ops.node(vid(i, j, k), side * i, side * j, float(k))
    for j in range(2):
        for i in range(2):
            ops.fix(vid(i, j, 0), 1, 1, 1)   # fixed base

    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    ops.nDMaterial("LogStrain", 2, 1)

    # straight-sided tets: BezierTet10 needs the 10-node connectivity, but for a
    # straight tet the 6 mid-edge nodes sit at edge midpoints — create them on the
    # fly per element so the geometry is exact.
    etag = 1
    nodecount = vid(1, 1, nz) + 1
    for k in range(nz):
        hexc = [vid(0, 0, k), vid(1, 0, k), vid(1, 1, k), vid(0, 1, k),
                vid(0, 0, k + 1), vid(1, 0, k + 1), vid(1, 1, k + 1), vid(0, 1, k + 1)]
        for tet in _kuhn_tets(*hexc):
            v = list(tet)
            mids = []
            for (a, b) in _EDGE_V:
                xa = np.array(ops.nodeCoord(v[a]))
                xb = np.array(ops.nodeCoord(v[b]))
                xm = 0.5 * (xa + xb)
                ops.node(nodecount, *[float(c) for c in xm])
                mids.append(nodecount)
                nodecount += 1
            ops.element("BezierTet10", etag, *v, *mids, 2, "-geom", "finite")
            etag += 1

    top = [vid(i, j, nz) for j in range(2) for i in range(2)]
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, 0.01, 0.0, 0.0)          # transverse tip shear, split over 4 nodes

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-9, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.25)
    ops.analysis("Static")

    tip_prev = 0.0
    for step in range(4):
        assert ops.analyze(1) == 0, f"finite cantilever step {step} failed to converge"
        tip = float(np.mean([ops.nodeDisp(t, 1) for t in top]))
        assert np.isfinite(tip), "tip displacement is not finite"
        assert tip > tip_prev, (
            f"tip deflection should grow monotonically under increasing load "
            f"(step {step}: {tip:.4e} <= prev {tip_prev:.4e})"
        )
        tip_prev = tip

    assert tip_prev > 0.0, "final tip deflection should be positive (transverse +x load)"
