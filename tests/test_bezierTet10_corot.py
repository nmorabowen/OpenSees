"""BezierTet10 ``-geom corot`` — element-level corotational battery (3D tet).

Tet sibling of test_ladrunoBrick_corot.py. ``corot`` strips the element rigid
rotation R = polar(H) (Procrustes best fit of the 10-node cloud), feeds the SMALL
deformational strain in the co-rotating frame to the EXISTING small-strain
material library (``ElasticIsotropic``, via ``setTrialStrain``), and rotates the
core force / stiffness back to global, adding the geometric (load) stiffness via
the shared SolidTransformation layer. No finite-strain material is involved.

LOAD-FRAME CONTRACT (Ladruno, Tet10 corot — see Ladruno_implementation/
14_bezier_tet10_corot.md §3): the INTERNAL force ∫Bᵀσ is rotated by R; fixed-
direction external loads (body force, the +z pressure hack, Q) are applied in the
GLOBAL frame AFTER globalize, so the fCore feeding K_geo is pure ∫Bᵀσ and the f/K
paths are trivially consistent. (-pressure under corot is rejected in v1.)

Gates, in priority order:
  * RIGID ROTATION ⇒ ZERO STRESS / ZERO FORCE (std + bbar) — the defining
    corotational property; exact to polar precision at any rotation magnitude.
  * SMALL DEFORMATION ⇒ corot reduces to ``linear`` (std + bbar kernel).
  * ``-geom linear`` is BIT-IDENTICAL to the default (no -geom) — the seam proof.
  * OBJECTIVITY: a small strain through a large rotation gives R·(unrotated force)
    and a frame-invariant corotated GP stress.
  * CONSISTENT TANGENT vs finite difference (THE arbiter of the pure-Bᵀσ Option-A
    assembly): the gap VANISHES ~linearly as strain→0 (implemented R k_d Rᵀ +
    G1+G1ᵀ correct), with a documented loose bound at finite strain.
  * bbar near-incompressible under corot stays more flexible than std.
  * Load-driven cantilever (a tet cube-stack) drives Newton through the geometric
    stiffness and converges.
  * sendSelf/recvSelf round-trip preserves -geom corot (no silent revert).
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

# A distorted (positive-volume) reference tet, so R = polar(H) and the spin map Γ
# are genuinely exercised (not a regular simplex).
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


def _make_material(material, tag=1, nu=_NU):
    if material == "elastic":
        ops.nDMaterial("ElasticIsotropic", tag, _E, nu)
    else:
        raise ValueError(material)
    return tag


def _build(geom, bbar=False, nu=_NU):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _COORD.items():
        ops.node(tag, x, y, z)
    _make_material("elastic", 1, nu)
    args = ["BezierTet10", 1, *range(1, 11), 1]
    if bbar:
        args.append("-bbar")
    if geom is not None:
        args += ["-geom", geom]
    ops.element(*args)
    return 1


def _impose_and_solve(u, geom, bbar=False, nu=_NU):
    _build(geom, bbar, nu)
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


def _resisting_force(u, geom="corot", bbar=False, nu=_NU):
    assert _impose_and_solve(u, geom, bbar, nu) == 0, "imposed-displacement solve failed"
    return np.array(ops.eleForce(1), dtype=float)


def _element_tangent(u, geom="corot", bbar=False):
    assert _impose_and_solve(u, geom, bbar) == 0, "imposed-displacement solve failed"
    K = ops.eleResponse(1, "stiff")
    assert K and len(K) == 30 * 30, "BezierTet10 must answer eleResponse('stiff') with 30x30"
    return np.array(K, dtype=float).reshape(30, 30)


def _affine_disp(Fbar, wiggle=False):
    """u_a = (F - I) X_a at every node. The field is AFFINE, so the Bernstein
    control-point coefficient at each node equals the field value there (the
    degree-2 blossom of a linear field) — true for vertices AND mid-edges."""
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
#  THE defining gate: pure rigid rotation -> zero stress and zero force.      #
#  Includes cases crossing π (179°, 180°, 270°), where R has an eigenvalue    #
#  near −1 and the polar iteration / spin-map inverse are most stressed.      #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("thz,thy", [
    (0.35, 0.25), (1.20, -0.80), (2.50, 0.40),
    (math.radians(179.0), 0.0), (math.pi, 0.0), (1.5 * math.pi, 0.0),
])
def test_corot_rigid_rotation_is_stress_free(thz, thy):
    R = _rot(thz, thy)
    u = _affine_disp(R)  # x_a = R X_a  =>  polar(H) = R exactly  =>  u_d = 0
    assert _impose_and_solve(u, "corot") == 0

    s = ops.eleResponse(1, "stresses")
    assert s and len(s) == 24, "BezierTet10 'stresses' response missing (4 GP x 6)"
    smax = max(abs(v) for v in s)
    assert smax <= 1.0e-7 * _E, f"rigid rotation induced stress {smax:.3e}"

    f = np.array(ops.eleForce(1), dtype=float)
    assert np.abs(f).max() <= 1.0e-7 * _E, f"rigid rotation induced force {np.abs(f).max():.3e}"


def test_corot_rigid_rotation_bbar_is_stress_free():
    R = _rot(0.9, 0.5)
    u = _affine_disp(R)
    assert _impose_and_solve(u, "corot", bbar=True) == 0
    f = np.array(ops.eleForce(1), dtype=float)
    assert np.abs(f).max() <= 1.0e-7 * _E, f"corot+bbar rigid rotation force {np.abs(f).max():.3e}"


# --------------------------------------------------------------------------- #
#  Seam proof: -geom linear is BIT-IDENTICAL to the default (no -geom).        #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("bbar", [False, True])
def test_linear_geom_is_bit_identical_to_default(bbar):
    Fbar = [[1.004, 0.002, -0.001],
            [0.002, 0.997, 0.003],
            [-0.001, 0.003, 1.002]]
    u = _affine_disp(Fbar, wiggle=True)
    f_default = _resisting_force(u, None, bbar)     # no -geom flag at all
    f_linear = _resisting_force(u, "linear", bbar)  # explicit -geom linear
    assert np.abs(f_default - f_linear).max() <= 1.0e-12 * max(np.abs(f_default).max(), 1.0), (
        "-geom linear is not bit-identical to the default (identity seam broken): "
        f"max diff {np.abs(f_default - f_linear).max():.3e}"
    )


# --------------------------------------------------------------------------- #
#  Small deformation: corot reduces to the linear (small-strain) kernel.      #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("bbar", [False, True])
def test_corot_reduces_to_linear_at_small_deformation(bbar):
    eps = 1.0e-5
    Fbar = [[1.0 + eps, 0.5 * eps, 0.0],
            [0.5 * eps, 1.0 - 0.3 * eps, 0.0],
            [0.0, 0.0, 1.0 + 0.2 * eps]]
    u = _affine_disp(Fbar)
    f_cor = _resisting_force(u, "corot", bbar)
    f_lin = _resisting_force(u, "linear", bbar)
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
    E = np.array([[1.02, 0.0, 0.0],
                  [0.0, 0.99, 0.0],
                  [0.0, 0.0, 1.005]])
    R = _rot(0.7, -0.4)

    assert _impose_and_solve(_affine_disp(E.tolist()), "corot") == 0
    f_ref = np.array(ops.eleForce(1), dtype=float)                # unrotated
    s_ref = np.array(ops.eleResponse(1, "stresses"), dtype=float)  # 24 = 4 GP x 6

    assert _impose_and_solve(_affine_disp((R @ E).tolist()), "corot") == 0
    f_rot = np.array(ops.eleForce(1), dtype=float)                # same strain, rotated
    s_rot = np.array(ops.eleResponse(1, "stresses"), dtype=float)

    # (1) force: rotated == R applied node-by-node to the unrotated force
    f_ref_rot = np.zeros(30)
    for a in range(10):
        f_ref_rot[3*a:3*a+3] = R @ f_ref[3*a:3*a+3]
    scale = max(np.abs(f_ref).max(), 1.0e-30)
    assert np.abs(f_rot - f_ref_rot).max() <= 1.0e-4 * scale, (
        f"corot not objective: rotated force != R·(unrotated force), "
        f"max err {np.abs(f_rot - f_ref_rot).max():.3e} vs scale {scale:.3e}"
    )
    # (2) the corotated GP stress is frame-INVARIANT (same u_d both cases)
    sscale = max(np.abs(s_ref).max(), 1.0e-30)
    assert np.abs(s_rot - s_ref).max() <= 1.0e-6 * sscale, (
        f"corot stress not frame-invariant: max diff {np.abs(s_rot - s_ref).max():.3e}"
    )


# --------------------------------------------------------------------------- #
#  Consistent tangent vs finite difference (the Option-A arbiter).             #
#  Like the brick, v2.0 omits two O(strain) geometric terms, so the FD gap     #
#  vanishes ~linearly with strain. Also a direct cross-check that the fCore     #
#  reconstructed inside getTangentStiff equals getResistingForce's internal     #
#  force (caught here as a non-vanishing or asymmetric tangent).                #
# --------------------------------------------------------------------------- #
def _fd_tangent(u, geom="corot", bbar=False, h=1.0e-6):
    Kfd = np.zeros((30, 30))
    for d in range(30):
        up = u.copy(); up[d] += h
        um = u.copy(); um[d] -= h
        Kfd[:, d] = (_resisting_force(up, geom, bbar)
                     - _resisting_force(um, geom, bbar)) / (2.0 * h)
    return Kfd


@pytest.mark.parametrize("bbar", [False, True])
def test_corot_tangent_gap_vanishes_with_strain(bbar):
    R = _rot(0.30, -0.20)
    base = np.array([[1.0, 0.4, -0.3], [0.4, -0.6, 0.2], [-0.3, 0.2, 0.8]])
    rel = {}
    for eps in (2.0e-3, 1.0e-3, 5.0e-4):
        u = _affine_disp((R @ (np.eye(3) + eps * base)).tolist())
        K = _element_tangent(u, "corot", bbar)
        assert np.isfinite(K).all()
        Kfd = _fd_tangent(u, "corot", bbar)
        rel[eps] = np.abs(K - Kfd).max() / max(np.abs(K).max(), 1e-30)
        fds = max(np.abs(Kfd).max(), 1e-30)
        assert np.abs(Kfd - Kfd.T).max() <= 1.0e-5 * fds, (
            f"corot FD tangent asymmetric @eps={eps}: geometric term mis-assembled"
        )
    assert rel[1.0e-3] < 0.75 * rel[2.0e-3], f"gap not vanishing: {rel}"
    assert rel[5.0e-4] < 0.75 * rel[1.0e-3], f"gap not vanishing: {rel}"
    assert rel[5.0e-4] < 5.0e-4, f"corot small-strain tangent gap too large: {rel}"


def test_corot_tangent_finite_strain_bound():
    R = _rot(0.30, -0.20)
    strain = np.array([[1.03, 0.015, -0.010],
                       [0.015, 0.975, 0.012],
                       [-0.010, 0.012, 1.020]])
    u = _affine_disp((R @ strain).tolist())
    K = _element_tangent(u, "corot")
    Kfd = _fd_tangent(u, "corot")
    scale = np.abs(K).max()
    err = np.abs(K - Kfd).max()
    print(f"corot finite-strain FD-tangent rel err = {err/scale:.3e} (deferred O(strain) terms)")
    assert err <= 2.0e-2 * scale, f"corot tangent gap {err/scale:.3e} exceeds deferred-term bound"


# --------------------------------------------------------------------------- #
#  bbar near-incompressible under corot: the volumetric split survives the     #
#  rotation strip. An AFFINE field gives constant strain where B̄ ≡ B (the      #
#  volume average equals the pointwise value), so it cannot distinguish bbar    #
#  from std — we must impose a genuinely NON-AFFINE (quadratic bending) field,  #
#  which produces a through-depth-varying volumetric strain that std over-      #
#  stiffens (locks) at ν→0.5 and bbar relieves. Imposed via the degree-2        #
#  Bernstein blossom (vertex DOF = field; mid-edge DOF = 2·f(mid) − ½f(a) − ½f(b))#
#  and superposed on a large rigid rotation R to prove it works in the          #
#  corotated frame.                                                             #
# --------------------------------------------------------------------------- #
def _bending_field(X, kappa):
    # pure bending about y: ε_xx = −κ z (linear through depth), ε_zz=γ_xz=0.
    x, y, z = X
    return np.array([-kappa * x * z, 0.0, 0.5 * kappa * x * x])


def _quadratic_corot_disp(field_fn, R):
    """Nodal DOFs for the deformed control net R·(X + u_def), with u_def imposed
    via the degree-2 blossom so the prescribed field is exactly quadratic."""
    u = np.zeros(30)
    for tag, xyz in _COORD.items():
        Xn = np.array(xyz)
        if tag <= 4:                                   # vertex: blossom = f(vertex)
            udef = field_fn(Xn)
        else:                                          # mid-edge: degree-2 blossom
            a, b = _EDGE_V[tag - 5]
            Xa, Xb = np.array(_VERTS[a + 1]), np.array(_VERTS[b + 1])
            udef = 2.0 * field_fn(0.5 * (Xa + Xb)) - 0.5 * field_fn(Xa) - 0.5 * field_fn(Xb)
        u[(tag - 1) * 3:(tag - 1) * 3 + 3] = R @ (Xn + udef) - Xn
    return u


def test_corot_bbar_relieves_volumetric_locking_under_rotation():
    nu = 0.4999
    R = _rot(0.6, -0.35)
    kappa = 0.02                                       # small curvature: small strain
    u = _quadratic_corot_disp(lambda X: _bending_field(X, kappa), R)
    f_std = _resisting_force(u, "corot", bbar=False, nu=nu)
    f_bbar = _resisting_force(u, "corot", bbar=True, nu=nu)
    n_std = np.abs(f_std).max()
    n_bbar = np.abs(f_bbar).max()
    print(f"corot near-incompressible bending: |f_std|={n_std:.4e} |f_bbar|={n_bbar:.4e}")
    # (1) bbar is genuinely active under corot (differs from std on a non-affine field)
    assert abs(n_bbar - n_std) > 1.0e-3 * n_std, (
        "corot+bbar == corot+std on a quadratic field — volumetric split inactive"
    )
    # (2) the locking signature: std over-stiffens the locked bending mode at ν→0.5,
    #     so its reaction is strictly larger than the bbar-relieved one.
    assert n_bbar < n_std, (
        f"corot+bbar not more flexible than corot+std at nu={nu}: "
        f"|f_bbar|={n_bbar:.4e} >= |f_std|={n_std:.4e} (volumetric split lost under rotation?)"
    )


# --------------------------------------------------------------------------- #
#  Load-driven cantilever: a tet cube-stack drives Newton through K_geo.       #
#  Each unit cube is split into 6 tets (Freudenthal, shared 0-6 diagonal);     #
#  midedges are deduped across the whole mesh. Base fixed, transverse tip load.#
# --------------------------------------------------------------------------- #
class _Mesh:
    """Vertex + deduped-midedge registry -> 10-node tet connectivities."""
    def __init__(self):
        self.coord = {}
        self._vkey = {}
        self._ekey = {}
        self._next = 1

    def _new(self, xyz):
        t = self._next
        self._next += 1
        self.coord[t] = tuple(float(c) for c in xyz)
        return t

    def vertex(self, xyz):
        k = tuple(round(c, 9) for c in xyz)
        if k not in self._vkey:
            self._vkey[k] = self._new(xyz)
        return self._vkey[k]

    def _mid(self, ta, tb):
        k = frozenset((ta, tb))
        if k not in self._ekey:
            a, b = self.coord[ta], self.coord[tb]
            self._ekey[k] = self._new(tuple(0.5 * (a[i] + b[i]) for i in range(3)))
        return self._ekey[k]

    def tet(self, v):
        conn = list(v)
        for (a, b) in _EDGE_V:
            conn.append(self._mid(v[a], v[b]))
        return conn

    def define(self):
        for t, (x, y, z) in self.coord.items():
            ops.node(t, x, y, z)


# Freudenthal 6-tet split of a cube (corner index -> bit coords), all sharing
# the 0-6 main diagonal. Corner order: 0=(0,0,0) 1=(1,0,0) 2=(1,1,0) 3=(0,1,0)
# 4=(0,0,1) 5=(1,0,1) 6=(1,1,1) 7=(0,1,1).
_CUBE_CORNERS = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
                 (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
_CUBE_TETS = [(0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6),
              (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6)]


def _build_cantilever(geom, nz=4, P=0.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mesh = _Mesh()
    ele_conn = []
    for k in range(nz):                       # stack nz unit cubes along z
        cv = [mesh.vertex((cx, cy, cz + k)) for (cx, cy, cz) in _CUBE_CORNERS]
        for tet in _CUBE_TETS:
            ele_conn.append(mesh.tet([cv[i] for i in tet]))
    mesh.define()
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    for e, conn in enumerate(ele_conn):
        ops.element("BezierTet10", e + 1, *conn, 1, "-geom", geom)

    base = [t for t, (x, y, z) in mesh.coord.items() if abs(z) < 1e-9]
    top = [t for t, (x, y, z) in mesh.coord.items() if abs(z - nz) < 1e-9]
    for t in base:
        ops.fix(t, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in top:
        ops.load(t, P / len(top), 0.0, 0.0)   # transverse (x) tip load
    return top


def _solve_cantilever(geom, nz=8, P=0.5, nsteps=30):
    # slender 1x1xnz column (aspect ~8) so a modest transverse load produces a
    # LARGE ROTATION at SMALL STRAIN — the corot regime. tol 1e-7 (not 1e-8):
    # the v2.0 corot tangent omits two O(strain) geometric terms, so late-step
    # Newton converges a touch more slowly than the exact-tangent finite path.
    top = _build_cantilever(geom, nz, P)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-7, 200, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            return None
    return float(np.mean([ops.nodeDisp(t, 1) for t in top]))   # mean tip x-disp


def test_corot_cantilever_converges_and_shows_geometric_effect():
    ux_cor = _solve_cantilever("corot")
    assert ux_cor is not None, "corot cantilever Newton failed to converge"
    assert math.isfinite(ux_cor) and ux_cor > 0.0, f"bad corot tip disp {ux_cor}"

    ux_lin = _solve_cantilever("linear")
    assert ux_lin is not None, "linear cantilever failed to converge"

    # Physical anchor (no -geom finite reference exists for Tet10): under a
    # transverse tip DEAD load the large-rotation response is geometrically
    # STIFFENING, so corot must give a SMALLER tip deflection than linear-geom,
    # and the reduction must sit in a sane band — a wrong-sign or mis-scaled
    # K_geo lands outside it (ratio>1 = softening, ratio→0 = runaway stiffening).
    ratio = ux_cor / ux_lin
    print(f"cantilever tip ux: corot={ux_cor:.5e} linear={ux_lin:.5e} ratio={ratio:.4f}")
    assert ux_cor < ux_lin, (
        f"corot tip {ux_cor:.4e} not stiffer than linear {ux_lin:.4e}: K_geo sign/scale wrong"
    )
    assert 0.4 < ratio < 0.98, (
        f"corot/linear tip ratio {ratio:.3f} outside the physical large-rotation band "
        "(0.4, 0.98) — geometric stiffness mis-scaled"
    )


# --------------------------------------------------------------------------- #
#  Serialization: sendSelf/recvSelf must preserve -geom corot (else a parallel  #
#  worker silently reverts to linear). Round-trips committed state via the      #
#  FE_Datastore (exercises the widened ID buffer + the broker).                 #
# --------------------------------------------------------------------------- #
def test_corot_serialization_roundtrip():
    def build():
        top = _build_cantilever("corot", nz=2, P=0.2)
        ops.constraints("Transformation")
        ops.numberer("RCM")
        ops.system("BandGeneral")
        ops.test("NormDispIncr", 1.0e-7, 100, 0)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        # MUST converge to a NON-ZERO state, else the round-trip compares 0→0 and
        # is vacuous (it would "pass" even if recvSelf silently reverted corot).
        assert ops.analyze(1) == 0, "serialization-roundtrip build failed to converge"
        build.top = top

    build.top = []
    # prime to discover the probe nodes, then hand the same builder to the helper
    build()
    probe = list(build.top)
    # guard against a vacuous 0→0 round-trip: the committed state must be non-trivial
    dmax = max(abs(ops.nodeDisp(t, 1)) for t in probe)
    assert dmax > 1.0e-6, f"roundtrip build produced ~zero state ({dmax:.2e}); test would be vacuous"
    database_roundtrip(build, probe, ndf=3)
