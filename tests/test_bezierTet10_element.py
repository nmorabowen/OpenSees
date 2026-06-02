"""BezierTet10 (Ladruno quadratic Bézier tetrahedron, classTag 33001) — T0 + T1.

3D sibling of test_bezierTri6_element.py. Element class merge gate = T0 + T1
(+ T2a convergence, deferred — see ADR 06_bezier_tet10, geometry-limited for
straight-sided curved benchmarks).

T0:
  * rank / zero-energy: a free element has EXACTLY 6 zero modes (3D rigid body).
  * constant-strain patch: a 4-tet "central node" patch (one interior vertex P +
    4 interior mid-edges, all free), linear field on the boundary -> every
    Gauss-point stress is the closed-form constant stress.
  * quadratic completeness: the distinguishing test for a *quadratic* element.
    A quadratic displacement field is in the element's span, so with the
    consistent body force b = -div(sigma) the patch reproduces it exactly and
    every Gauss-point strain equals the analytic (linear) strain at that GP's
    physical location.
T1:
  * uniaxial bar (displacement-controlled): total support reaction = E*eps*A and
    the lateral contraction = -nu*eps*L, both exact for a constant-strain-
    reproducing element.

CRITICAL Bezier subtlety (same as BezierTri6): the DOFs are Bernstein
CONTROL-POINT COEFFICIENTS, not nodal field values — the basis is
non-interpolatory at mid-edge nodes. To impose a field p(x) you prescribe the
control net, not point samples:
  * vertex DOF  = p(vertex)
  * mid-edge DOF on edge (va,vb) = 2*p(mid) - 0.5*p(va) - 0.5*p(vb)  (degree-2
    blossom / polar form). For a LINEAR p this reduces to p(mid).

Command: element('BezierTet10', tag, n1..n10, matTag, '-bodyForce', bx,by,bz).
Node order = TenNodeTetrahedron (D2'): 4 vertices, then mid-edges
(1-2)(2-3)(1-3)(1-4)(3-4)(2-4). Straight-sided only. NGAUSS = 4, ndf = 3.
"""
import math

import pytest

from _testbed import ops
from _testbed.fem_checks import assert_zero_energy, check_constant_stress
from _testbed.equilibrium import assert_equilibrium

pytestmark = [pytest.mark.zone_a]

E = 1000.0
NU = 0.25
NGAUSS = 4

# 3D isotropic elasticity, Voigt {xx,yy,zz,xy,yz,zx}
LAM = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU))
MU = E / (2.0 * (1.0 + NU))

# mid-edge (vertexA, vertexB) order for nodes 5..10 (D2', matches the element)
EDGE_V = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]


def _D_times(eps):
    """sigma = D : eps for {exx,eyy,ezz,gxy,gyz,gzx} (engineering shears)."""
    exx, eyy, ezz, gxy, gyz, gzx = eps
    tr = exx + eyy + ezz
    return [LAM * tr + 2 * MU * exx,
            LAM * tr + 2 * MU * eyy,
            LAM * tr + 2 * MU * ezz,
            MU * gxy, MU * gyz, MU * gzx]


# ─────────────────────────────────────────────────────────────────────────
#  mesh helper: global node registry (vertices + deduped mid-edges)
# ─────────────────────────────────────────────────────────────────────────
class _Mesh:
    def __init__(self):
        self.coord = {}        # tag -> (x,y,z)
        self._vkey = {}        # rounded vertex coord -> tag
        self._ekey = {}        # frozenset(vtagA,vtagB) -> mid tag
        self._next = 1

    def _new(self, xyz):
        t = self._next
        self._next += 1
        self.coord[t] = xyz
        return t

    def vertex(self, xyz):
        k = tuple(round(c, 9) for c in xyz)
        if k not in self._vkey:
            self._vkey[k] = self._new(tuple(float(c) for c in xyz))
        return self._vkey[k]

    def _mid(self, ta, tb):
        k = frozenset((ta, tb))
        if k not in self._ekey:
            a, b = self.coord[ta], self.coord[tb]
            self._ekey[k] = self._new(tuple(0.5 * (a[i] + b[i]) for i in range(3)))
        return self._ekey[k]

    def tet(self, v):
        """v = 4 vertex tags -> 10-node connectivity in D2' order."""
        conn = list(v)
        for (a, b) in EDGE_V:
            conn.append(self._mid(v[a], v[b]))
        return conn

    def define(self):
        for t, (x, y, z) in self.coord.items():
            ops.node(t, x, y, z)


# distorted "central node" patch: outer tet ABCD + interior P -> 4 sub-tets
# (P + each face).  P and the four P-edge midsides are interior (free); the
# outer corners and outer-edge midsides are the boundary.
_A = (0.0, 0.0, 0.0)
_B = (2.0, 0.0, 0.0)
_C = (0.3, 1.8, 0.0)
_D = (0.4, 0.5, 1.6)
_P = (0.65, 0.55, 0.45)


def _build_patch(body=(0.0, 0.0, 0.0)):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    m = _Mesh()
    a, b, c, d = (m.vertex(_A), m.vertex(_B), m.vertex(_C), m.vertex(_D))
    p = m.vertex(_P)
    tets = [(p, a, b, c), (p, a, b, d), (p, a, c, d), (p, b, c, d)]
    conns = [m.tet(t) for t in tets]
    m.define()
    for e, conn in enumerate(conns, start=1):
        ops.element("BezierTet10", e, *conn, 1, "-bodyForce", *body)
    # boundary = nodes ON a face of ABCD; interior = P and P-edge midsides.
    corners = {a, b, c, d}
    interior = {p}
    for (ta, tb) in [(p, a), (p, b), (p, c), (p, d)]:
        interior.add(m._ekey[frozenset((ta, tb))])
    boundary = [t for t in m.coord if t not in interior]
    return m, boundary, list(interior), len(conns)


# ───────────────────────────────────────────────────────── T0: rank
@pytest.mark.t0
def test_rank_zero_energy():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    m = _Mesh()
    v = [m.vertex(_A), m.vertex(_B), m.vertex(_C), m.vertex(_D)]
    conn = m.tet(v)
    m.define()
    ops.element("BezierTet10", 1, *conn, 1)
    assert_zero_energy(list(m.coord), ndf=3, ndm=3)   # exactly 6 rigid-body modes


# ───────────────────────────────────────────────────────── T0: patch
@pytest.mark.t0
def test_constant_strain_patch():
    eps = [1.0e-3, 4.0e-4, -3.0e-4, 2.0e-4, -1.5e-4, 1.0e-4]  # xx,yy,zz,xy,yz,zx
    exx, eyy, ezz, gxy, gyz, gzx = eps

    def ux(x, y, z): return exx * x + 0.5 * gxy * y + 0.5 * gzx * z
    def uy(x, y, z): return 0.5 * gxy * x + eyy * y + 0.5 * gyz * z
    def uz(x, y, z): return 0.5 * gzx * x + 0.5 * gyz * y + ezz * z

    m, boundary, interior, nele = _build_patch()
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in boundary:                       # linear field: control coeff = value
        x, y, z = m.coord[n]
        ops.sp(n, 1, ux(x, y, z))
        ops.sp(n, 2, uy(x, y, z))
        ops.sp(n, 3, uz(x, y, z))
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static", "-noWarnings")
    ops.analyze(1)

    target = _D_times(eps)
    sref = max(abs(v) for v in eps)
    for e in range(1, nele + 1):
        check_constant_stress(e, NGAUSS, target, E=E, strain_ref=sref)


# ─────────────────────────────────────────── T0: quadratic completeness
def _blossom(fn, mid, va, vb):
    """degree-2 Bernstein control coeff for a mid-edge DOF of a polynomial fn."""
    return 2.0 * fn(*mid) - 0.5 * fn(*va) - 0.5 * fn(*vb)


@pytest.mark.t0
def test_quadratic_completeness_patch():
    a, c, g = 1.0e-4, -0.6e-4, 0.4e-4       # u_x=a x^2, u_y=c y^2, u_z=g z^2
    def ux(x, y, z): return a * x * x
    def uy(x, y, z): return c * y * y
    def uz(x, y, z): return g * z * z
    # sigma* linear -> div sigma* = ((lam+2mu)*2a, (lam+2mu)*2c, (lam+2mu)*2g)
    k = LAM + 2 * MU
    body = (-k * 2 * a, -k * 2 * c, -k * 2 * g)

    m, boundary, interior, nele = _build_patch(body=body)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    # classify each boundary node: vertex (4 corners) -> field value;
    # mid-edge -> blossom from its two endpoints.
    corners = set(m._vkey.values())
    mid_of = {v: tuple(ab) for ab, v in m._ekey.items()}   # mid tag -> (vtagA, vtagB)
    for n in boundary:
        if n in corners:
            x, y, z = m.coord[n]
            ops.sp(n, 1, ux(x, y, z)); ops.sp(n, 2, uy(x, y, z)); ops.sp(n, 3, uz(x, y, z))
        else:
            ta, tb = mid_of[n]
            va, vb, mid = m.coord[ta], m.coord[tb], m.coord[n]
            ops.sp(n, 1, _blossom(ux, mid, va, vb))
            ops.sp(n, 2, _blossom(uy, mid, va, vb))
            ops.sp(n, 3, _blossom(uz, mid, va, vb))
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static", "-noWarnings")
    ops.analyze(1)

    # every GP strain == analytic linear strain (2a x, 2c y, 2g z, 0,0,0)
    for e in range(1, nele + 1):
        gp = ops.eleResponse(e, "gpCoord")      # [x,y,z]*4
        st = ops.eleResponse(e, "strain")       # [exx..gzx]*4
        assert gp and st, f"element {e} missing gpCoord/strain response"
        for kgp in range(NGAUSS):
            xg, yg, zg = gp[3*kgp], gp[3*kgp+1], gp[3*kgp+2]
            s = st[6*kgp:6*kgp+6]
            assert abs(s[0] - 2*a*xg) <= 1e-8 * abs(2*a*xg) + 1e-11
            assert abs(s[1] - 2*c*yg) <= 1e-8 * abs(2*c*yg) + 1e-11
            assert abs(s[2] - 2*g*zg) <= 1e-8 * abs(2*g*zg) + 1e-11
            assert max(abs(v) for v in s[3:]) <= 1e-9


# ───────────────────────────────────────────────────────── T1: uniaxial
def test_uniaxial_reaction_and_poisson():
    # a box [0,Lx]x[0,Ly]x[0,Lz] meshed as ONE hex split into 6 tets (Kuhn),
    # uniaxial strain eps along x via prescribed end displacement.
    Lx, Ly, Lz = 3.0, 1.0, 1.0
    delta = 0.009
    eps = delta / Lx

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    m = _Mesh()
    # 8 hex corners, id = i + 2j + 4k
    corner = {}
    for kz in (0, 1):
        for jy in (0, 1):
            for ix in (0, 1):
                corner[(ix, jy, kz)] = m.vertex((ix*Lx, jy*Ly, kz*Lz))
    cid = lambda i, j, kk: corner[(i, j, kk)]
    # Kuhn 6-tet subdivision around the main diagonal (0,0,0)-(1,1,1)
    kuhn = [
        ((0,0,0),(1,0,0),(1,1,0),(1,1,1)),
        ((0,0,0),(1,1,0),(0,1,0),(1,1,1)),
        ((0,0,0),(0,1,0),(0,1,1),(1,1,1)),
        ((0,0,0),(0,1,1),(0,0,1),(1,1,1)),
        ((0,0,0),(0,0,1),(1,0,1),(1,1,1)),
        ((0,0,0),(1,0,1),(1,0,0),(1,1,1)),
    ]
    conns = [m.tet([cid(*c) for c in tet]) for tet in kuhn]
    m.define()
    for e, conn in enumerate(conns, start=1):
        ops.element("BezierTet10", e, *conn, 1)

    tol = 1e-9
    x0 = [t for t, (x, y, z) in m.coord.items() if abs(x) < tol]          # x=0 face
    xL = [t for t, (x, y, z) in m.coord.items() if abs(x - Lx) < tol]     # x=Lx face
    for n in x0:
        ops.fix(n, 1, 0, 0)                       # roller: u_x=0
    # pin rigid body in y,z at one corner and y at another
    o = corner[(0, 0, 0)]
    ops.fix(o, 0, 1, 1)
    ops.fix(corner[(0, 1, 0)], 0, 0, 1)
    ops.fix(corner[(0, 0, 1)], 0, 1, 0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in xL:
        ops.sp(n, 1, delta)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static", "-noWarnings")
    ops.analyze(1)
    ops.reactions()

    # total axial reaction on the fixed face = E*eps*A  (uniaxial stress state)
    Rx = sum(ops.nodeReaction(n, 1) for n in x0)
    A = Ly * Lz
    assert abs(abs(Rx) - E * eps * A) <= 1e-6 * E * eps * A, \
        f"Rx={Rx} vs E*eps*A={E*eps*A}"

    # Poisson lateral contraction at the far corner: u_y(y=Ly) = -nu*eps*Ly
    far = corner[(1, 1, 1)]
    assert abs(ops.nodeDisp(far, 2) - (-NU * eps * Ly)) <= 1e-7 + 1e-7 * NU * eps * Ly

    # global equilibrium closes (no external x-load: end reactions balance)
    assert_equilibrium(list(m.coord), n_force_dofs=3,
                       applied_total=[0.0, 0.0, 0.0], atol=1e-6)


# ─────────────────────────────────────── T0: characteristic length
#
# Crack-band materials (ASDConcrete3D) regularize by
# ops_TheActiveElement->getCharacteristicLength(), read once on the first
# setTrialStrain. BezierTet10 OVERRIDES Element's min-inter-node-distance
# default (which on a quadratic element collapses to ~½ the edge length:
# corner-to-mid-edge) with a volume-based element size lch = cbrt(6·V).
#
# For a right-corner tetrahedron with three mutually perpendicular legs of
# length L, V = L³/6, so the override returns lch = cbrt(6·L³/6) = L exactly.
# The inherited default would return the min node spacing L/2 — so this test
# pins the formula AND fails loudly if the override is reverted (½× regression).

@pytest.mark.t0
@pytest.mark.parametrize("L", [1.0, 2.5, 0.3])
def test_characteristic_length(L):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    m = _Mesh()
    # right-corner tet: legs along +x,+y,+z from the origin
    v = [m.vertex((0.0, 0.0, 0.0)), m.vertex((L, 0.0, 0.0)),
         m.vertex((0.0, L, 0.0)), m.vertex((0.0, 0.0, L))]
    conn = m.tet(v)
    m.define()
    ops.element("BezierTet10", 1, *conn, 1)

    lch = ops.eleResponse(1, "charLength")
    assert lch, "BezierTet10 missing 'charLength' response"
    lch = lch[0]
    # volume-based override: lch = cbrt(6·V) = L for the right-corner tet
    assert abs(lch - L) <= 1e-12 * L + 1e-14, f"lch={lch}, expected {L}"
    # regression guard: must NOT be the inherited min-distance default (L/2)
    assert lch > 0.75 * L, f"lch={lch} looks like the ½-edge min-distance default"
