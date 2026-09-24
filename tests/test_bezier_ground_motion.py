"""Ground-motion (UniformExcitation) sign + magnitude for BezierTri6 / BezierTet10 -- WP-117.

Both elements built their ground-motion load as Q += +M*R*a_g while
getResistingForce() subtracts Q, so the element mass was driven by -a_g: a
UniformExcitation run shook a Bezier mesh the wrong way, and a mixed model
(Bezier soil + nodal-mass structure) in opposite directions. Every other fork
element and vanilla OpenSees use Q += -M*R*a_g. No test ever ran a Bezier
element under UniformExcitation.

Gates:
  t1  rigid-body probe -- all nodes free in x, fixed otherwise; a rigid
      x-translation has zero stiffness, so under a constant ground
      acceleration a_g every node's RELATIVE acceleration must be exactly
      -a_g. Legs: {lumped, consistent (-cMass)} x {element -rho, material rho}.
  t1  deformable differential -- a fixed-base column under UniformExcitation
      must reproduce, step for step, the same column driven by the equivalent
      nodal loads -m_i*a_g (lumped mass: rho*A*t/6 per Tri6 node, rho*V/10 per
      Tet10 node). Pins the magnitude as well as the sign.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, NU, RHO, THK = 1000.0, 0.3, 1.5, 1.0
AG = 2.0
_TET_EDGES = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]   # e12 e23 e13 e14 e34 e24


class _Mesh:
    def __init__(self, dim):
        self.dim, self.coord, self.key, self.mids = dim, {}, {}, {}

    def node(self, xyz):
        k = tuple(round(c, 9) for c in xyz)
        if k not in self.key:
            t = len(self.coord) + 1
            self.key[k], self.coord[t] = t, tuple(float(c) for c in xyz)
            ops.node(t, *self.coord[t])
        return self.key[k]

    def mid(self, a, b):
        k = frozenset((a, b))
        if k not in self.mids:
            ca, cb = self.coord[a], self.coord[b]
            self.mids[k] = self.node(tuple(0.5 * (ca[i] + cb[i]) for i in range(self.dim)))
        return self.mids[k]


def _material(mat_rho):
    if mat_rho:
        ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
    else:
        ops.nDMaterial("ElasticIsotropic", 1, E, NU)


def _tri6_mesh(ncell, cmass, mat_rho):
    """1 x ncell column of unit cells, 2 BezierTri6 per cell. Returns (mesh, lumped node masses)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    _material(mat_rho)
    m, lumped, e = _Mesh(2), {}, 0
    for j in range(ncell):
        n1, n2 = m.node((0.0, float(j))), m.node((1.0, float(j)))
        n3, n4 = m.node((1.0, j + 1.0)), m.node((0.0, j + 1.0))
        for c in ((n1, n2, n3), (n1, n3, n4)):
            conn = [*c, m.mid(c[0], c[1]), m.mid(c[1], c[2]), m.mid(c[2], c[0])]
            e += 1
            args = ["BezierTri6", e, *conn, THK, "PlaneStrain", 1]
            if not mat_rho:
                args += ["-rho", RHO]
            if cmass:
                args += ["-cMass"]
            ops.element(*args)
            for t in conn:                      # lumped: rho * A_e * t / 6 per node (A_e = 1/2)
                lumped[t] = lumped.get(t, 0.0) + RHO * 0.5 * THK / 6.0
    return m, lumped


def _tet10_mesh(ncell, cmass, mat_rho):
    """1 x 1 x ncell column of unit cubes, 6 Kuhn BezierTet10 per cube."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _material(mat_rho)
    m, lumped, e = _Mesh(3), {}, 0

    def vol6(v):
        p = [m.coord[t] for t in v]
        a, b, c = ([p[i][k] - p[0][k] for k in range(3)] for i in (1, 2, 3))
        return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))

    for k in range(ncell):
        for perm in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
            p, path = [0, 0, 0], [m.node((0, 0, k))]
            for ax in perm:
                p[ax] = 1
                path.append(m.node((p[0], p[1], p[2] + k)))
            if vol6(path) < 0.0:
                path[1], path[2] = path[2], path[1]
            conn = path + [m.mid(path[a], path[b]) for a, b in _TET_EDGES]
            e += 1
            args = ["BezierTet10", e, *conn, 1]
            if not mat_rho:
                args += ["-rho", RHO]
            if cmass:
                args += ["-cMass"]
            ops.element(*args)
            for t in conn:                      # lumped: rho * V_e / 10 per node (V_e = 1/6)
                lumped[t] = lumped.get(t, 0.0) + RHO * (1.0 / 6.0) / 10.0
    return m, lumped


_BUILD = {"BezierTri6": _tri6_mesh, "BezierTet10": _tet10_mesh}


def _transient():
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-12, 30)
    ops.analysis("Transient")


# --------------------------------------------------------------------------
# gate 1: rigid-body probe -- relative acceleration must be exactly -a_g
# --------------------------------------------------------------------------
@pytest.mark.parametrize("mat_rho", [False, True], ids=["element_rho", "material_rho"])
@pytest.mark.parametrize("cmass", [False, True], ids=["lumped", "consistent"])
@pytest.mark.parametrize("kind", ["BezierTri6", "BezierTet10"])
def test_rigid_body_relative_accel_is_minus_ag(kind, cmass, mat_rho):
    m, _ = _BUILD[kind](1, cmass, mat_rho)
    for t in m.coord:
        ops.fix(t, 0, *([1] * (m.dim - 1)))      # free in x only
    ops.timeSeries("Constant", 1, "-factor", AG)
    ops.pattern("UniformExcitation", 1, 1, "-accel", 1)
    _transient()
    assert ops.analyze(1, 0.01) == 0
    acc = [ops.nodeAccel(t, 1) for t in m.coord]
    worst = max(abs(a + AG) for a in acc)
    assert worst <= 1e-9 * AG, (
        f"[{kind}] relative x-acceleration {min(acc):+.6f}..{max(acc):+.6f}, expected {-AG} "
        "(+a_g means the ground-motion inertia load has the wrong sign)")


# --------------------------------------------------------------------------
# gate 2: UniformExcitation == equivalent nodal loads -m_i * a_g
# --------------------------------------------------------------------------
def _column_history(kind, via_ground):
    m, lumped = _BUILD[kind](2, False, False)
    base = [t for t, c in m.coord.items() if abs(c[m.dim - 1]) < 1e-12]
    for t in base:
        ops.fix(t, *([1] * m.dim))
    ops.timeSeries("Constant", 1, "-factor", AG)
    if via_ground:
        ops.pattern("UniformExcitation", 1, 1, "-accel", 1)
    else:
        ops.pattern("Plain", 1, 1)
        for t, mi in lumped.items():
            if t not in base:
                ops.load(t, -mi, *([0.0] * (m.dim - 1)))   # x 'Constant' factor AG => -m_i*a_g
    _transient()
    lam = ops.eigen("-fullGenLapack", 1)[0]
    nstep = 100
    dt = 1.5 * (2.0 * math.pi / math.sqrt(lam)) / nstep
    top = max(m.coord, key=lambda t: (m.coord[t][m.dim - 1], -sum(m.coord[t][:m.dim - 1])))
    out = []
    for _ in range(nstep):
        assert ops.analyze(1, dt) == 0
        out.append(ops.nodeDisp(top, 1))
    return out


@pytest.mark.parametrize("kind", ["BezierTri6", "BezierTet10"])
def test_uniform_excitation_equals_equivalent_nodal_loads(kind):
    ground = _column_history(kind, via_ground=True)
    loads = _column_history(kind, via_ground=False)
    scale = max(abs(u) for u in loads)
    assert scale > 0.0
    assert ground[-1] * loads[-1] > 0.0 or abs(ground[-1]) < 1e-12 * scale, "opposite sign"
    worst = max(abs(a - b) for a, b in zip(ground, loads))
    assert worst <= 1e-9 * scale, (
        f"[{kind}] UniformExcitation departs from the equivalent -m*a_g nodal-load run by "
        f"{worst / scale:.3e} of peak")
