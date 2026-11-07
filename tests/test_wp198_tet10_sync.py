"""WP-198: TenNodeTetrahedron synced with jaabell/OpenSees `fix/tet10-integration`.

Jose Abell fixed these TenNodeTetrahedron defects on jaabell/OpenSees
`ladruño` (upstream PR OpenSees/OpenSees#1833) after the fork's own Tet10
fixes had gone the other way. Each case below is his check, rebuilt here,
and fails on the fork before WP-198:

  * inertia load sized for 4 nodes / 1 Gauss point (ac9e5f2dc): a 48-element
    unit cube, all nodes fixed, rho = 2, constant ground acceleration 3 gave a
    base reaction of 7.042 (static) / 4.392 (transient) instead of 6;
  * shape-function derivatives from the corner coordinates instead of the
    full Jacobian (1d46bff1c): a curved element failed the patch test by up
    to 0.5 in stress;
  * setParameter reached Gauss point 1 only; setDispInit was unreachable;
  * `-doInitDisp` was rejected unless all three body forces were given;
  * a deactivated element still returned its initial stiffness;
  * an inverted element was accepted with no warning;
  * `print -JSON` wrote 3 of the 10 nodes;
  * `eleResponse $e material` without a point number read past argv.

TenNodeTetrahedronStd (node 9 on edge 2-4, node 10 on edge 3-4) is the same
element with the other mid-edge order; it must give the same answers when its
connectivity lists the same physical nodes in its own order.
"""

import itertools
import json
import os
import subprocess

import pytest

from _testbed import ops
from _testbed.subprocess_run import pinned_child


pytestmark = [pytest.mark.zone_a]

_PENALTY = 1.0e18

# Gmsh order (TenNodeTetrahedron): v1 v2 v3 v4, e12 e23 e13 e14 e34 e24.
# Standard order (TenNodeTetrahedronStd): the last two swapped, e24 e34.
_EDGES_GMSH = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]
_EDGES_STD = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
_EDGES = {'TenNodeTetrahedron': _EDGES_GMSH, 'TenNodeTetrahedronStd': _EDGES_STD}
_TYPES = list(_EDGES)

_REF = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]


def _signed6v(p):
    a = [p[1][i] - p[0][i] for i in range(3)]
    b = [p[2][i] - p[0][i] for i in range(3)]
    c = [p[3][i] - p[0][i] for i in range(3)]
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def _mid(pa, pb):
    return tuple(0.5 * (pa[j] + pb[j]) for j in range(3))


# --------------------------------------------------------------------------
# 1. Inertia load (the confirmed fork bug)
# --------------------------------------------------------------------------

def _build_cube(ele_type, rho):
    """Unit cube, 2x2x2 sub-cubes, 6 Kuhn tets each = 48 quadratic tets."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.nDMaterial('ElasticIsotropic', 1, 1000.0, 0.25, rho)

    tags = {}

    def node(xyz):
        key = tuple(round(v, 12) for v in xyz)
        if key not in tags:
            tags[key] = len(tags) + 1
            ops.node(tags[key], *key)
        return tags[key]

    n, h = 2, 0.5
    ele = 0
    for i, j, k in itertools.product(range(n), repeat=3):
        o = (i * h, j * h, k * h)
        for perm in itertools.permutations(range(3)):
            p = [list(o)]
            for axis in perm:
                q = list(p[-1])
                q[axis] += h
                p.append(q)
            p = [tuple(v) for v in p]
            if _signed6v(p) < 0.0:
                p[1], p[2] = p[2], p[1]
            corners = [node(v) for v in p]
            mids = [node(_mid(p[a], p[b])) for a, b in _EDGES[ele_type]]
            ele += 1
            ops.element(ele_type, ele, *corners, *mids, 1)
    assert ele == 48
    for t in tags.values():
        ops.fix(t, 1, 1, 1)
    return list(tags.values())


def _base_reaction(ele_type, transient):
    rho, ag = 2.0, 3.0
    nodes = _build_cube(ele_type, rho)
    ops.timeSeries('Constant', 1, '-factor', ag)
    ops.pattern('UniformExcitation', 1, 1, '-accel', 1)
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-12, 10, 0)
    ops.algorithm('Linear')
    if transient:
        ops.integrator('Newmark', 0.5, 0.25)
        ops.analysis('Transient')
        assert ops.analyze(1, 0.01) == 0
        ops.reactions('-dynamic')
    else:
        ops.integrator('LoadControl', 1.0)
        ops.analysis('Static')
        assert ops.analyze(1) == 0
        ops.reactions()
    return sum(ops.nodeReaction(t, 1) for t in nodes), rho * 1.0 * ag


@pytest.mark.parametrize('ele_type', _TYPES)
@pytest.mark.parametrize('transient', [False, True], ids=['static', 'transient'])
def test_ground_motion_base_reaction_is_rho_v_ag(ele_type, transient):
    """Total base reaction of a fully fixed cube = rho * V * ag = 6.

    Before WP-198 addInertiaLoadToUnbalance filled 12 entries (4 nodes,
    1 Gauss point) and getResistingForceIncInertia returned a 12-entry
    vector: 7.042 static, 4.392 transient, with no warning.
    """
    total, expected = _base_reaction(ele_type, transient)
    assert abs(total) == pytest.approx(expected, rel=1e-9), total


# --------------------------------------------------------------------------
# 2. Single-element patch tests (straight and curved edges)
# --------------------------------------------------------------------------

_E, _NU, _EPS = 1000.0, 0.25, 1.0e-3
_SIG = (1.2, 0.4, 0.4, 0.0, 0.0, 0.0)   # exact stress for u = eps*x, nu = 0.25

# Deterministic off-edge offsets for the six mid-edge nodes (gmsh edge order).
_OFFSETS = [(0.3, -0.5, 0.8), (-0.6, 0.2, 0.4), (0.5, 0.7, -0.3),
            (-0.4, -0.2, 0.6), (0.2, -0.8, -0.5), (0.7, 0.4, 0.1)]


def _geometry(curve):
    """Corner and mid-edge coordinates keyed by physical edge (a, b)."""
    mids = {}
    for (a, b), off in zip(_EDGES_GMSH, _OFFSETS):
        m = _mid(_REF[a], _REF[b])
        mids[(a, b)] = tuple(m[j] + curve * off[j] for j in range(3))
    return mids


def _single(ele_type, curve=0.0, extra=()):
    """One element, u = (eps*x, 0, 0) prescribed (Penalty) at all 10 nodes."""
    mids = _geometry(curve)
    pts = list(_REF) + [mids[e] for e in _EDGES[ele_type]]
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.nDMaterial('ElasticIsotropic', 1, _E, _NU)
    for t, xyz in enumerate(pts, start=1):
        ops.node(t, *xyz)
    ops.element(ele_type, 1, *range(1, 11), 1, *extra)
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for t, xyz in enumerate(pts, start=1):
        ops.sp(t, 1, _EPS * xyz[0])
        ops.sp(t, 2, 0.0)
        ops.sp(t, 3, 0.0)
    ops.constraints('Penalty', _PENALTY, _PENALTY)
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.test('NormDispIncr', 1.0e-12, 20, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')


def _stresses():
    s = ops.eleResponse(1, 'stresses')
    assert len(s) == 24
    return [s[6 * g:6 * g + 6] for g in range(4)]


@pytest.mark.parametrize('ele_type', _TYPES)
@pytest.mark.parametrize('curve', [0.0, 0.1, 0.3])
def test_curved_element_passes_the_patch_test(ele_type, curve):
    """Uniform strain is exact on a curved quadratic tet.

    Before WP-198 the derivative cofactors came from the corner coordinates
    while the determinant came from the full Jacobian: max stress error
    0.14-0.22 (offsets ~0.1) and 0.34-0.50 (~0.3); straight edges were exact.
    """
    _single(ele_type, curve)
    assert ops.analyze(1) == 0
    for g, sig in enumerate(_stresses()):
        for c in range(6):
            assert sig[c] == pytest.approx(_SIG[c], abs=1e-8), (g, c, sig)


def test_both_orders_give_the_same_element():
    """Gmsh and Std order, same physical nodes: identical stresses."""
    res = {}
    for ele_type in _TYPES:
        mids = _geometry(0.2)
        pts = list(_REF) + [mids[e] for e in _EDGES[ele_type]]
        ops.wipe()
        ops.model('basic', '-ndm', 3, '-ndf', 3)
        ops.nDMaterial('ElasticIsotropic', 1, _E, _NU)
        for t, xyz in enumerate(pts, start=1):
            ops.node(t, *xyz)
        ops.element(ele_type, 1, *range(1, 11), 1)
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 0, 1, 1)
        ops.fix(3, 1, 0, 1)
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        ops.load(4, 0.3, -0.2, 1.0)           # corner node 4: same in both orders
        ops.constraints('Plain')
        ops.numberer('Plain')
        ops.system('FullGeneral')
        ops.test('NormDispIncr', 1.0e-14, 10, 0)
        ops.algorithm('Linear')
        ops.integrator('LoadControl', 1.0)
        ops.analysis('Static')
        assert ops.analyze(1) == 0
        res[ele_type] = (ops.nodeDisp(4), _stresses())
    (u_g, s_g), (u_s, s_s) = res['TenNodeTetrahedron'], res['TenNodeTetrahedronStd']
    assert u_s == pytest.approx(u_g, rel=1e-10)
    for g in range(4):
        assert s_s[g] == pytest.approx(s_g[g], rel=1e-9, abs=1e-12)


# --------------------------------------------------------------------------
# 3. Parameters
# --------------------------------------------------------------------------

def test_set_parameter_reaches_gauss_point_4():
    """`material 4 E` used to be refused (pointNum <= 1)."""
    _single('TenNodeTetrahedron')
    assert ops.analyze(1) == 0
    ops.parameter(1, 'element', 1, 'material', 4, 'E')
    ops.updateParameter(1, 2.0 * _E)
    assert ops.analyze(1) == 0
    sig = _stresses()
    for g in range(3):
        assert sig[g][0] == pytest.approx(1.2, rel=1e-9)
    assert sig[3][0] == pytest.approx(2.4, rel=1e-9)


def test_set_parameter_all_points_form_reaches_every_gauss_point():
    """`E` with no point number used to change Gauss point 1 only."""
    _single('TenNodeTetrahedron')
    assert ops.analyze(1) == 0
    ops.parameter(1, 'element', 1, 'E')
    ops.updateParameter(1, 3.0 * _E)
    assert ops.analyze(1) == 0
    for g, sig in enumerate(_stresses()):
        assert sig[0] == pytest.approx(3.6, rel=1e-9), g


def test_set_disp_init_is_reachable():
    """`setDispInit` matched no string; now it zeroes the strain reference."""
    _single('TenNodeTetrahedron')
    assert ops.analyze(1) == 0
    assert _stresses()[0][0] == pytest.approx(1.2, rel=1e-9)
    ops.parameter(1, 'element', 1, 'setDispInit')
    ops.updateParameter(1, 1)
    assert ops.analyze(1) == 0
    for sig in _stresses():
        assert abs(sig[0]) < 1e-9


@pytest.mark.parametrize('extra', [('-doInitDisp', 1), (0.0, '-doInitDisp', 1),
                                   (0.0, 0.0, 0.0, '-doInitDisp', 1)])
def test_do_init_disp_parses_with_any_body_force_prefix(extra):
    """`... $matTag -doInitDisp 1` was rejected as invalid double data."""
    _single('TenNodeTetrahedron', extra=extra)
    assert ops.getEleTags() == [1]
    assert ops.analyze(1) == 0
    assert _stresses()[0][0] == pytest.approx(1.2, rel=1e-9)


def _deactivation_model():
    """Tet with only node 4 free in z, plus a z spring k = 100 on node 4."""
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.nDMaterial('ElasticIsotropic', 1, _E, _NU)
    ops.uniaxialMaterial('Elastic', 2, 100.0)
    pts = list(_REF) + [_mid(_REF[a], _REF[b]) for a, b in _EDGES_GMSH]
    for t, xyz in enumerate(pts, start=1):
        ops.node(t, *xyz)
        ops.fix(t, 1, 1, 0 if t == 4 else 1)
    ops.node(11, *_REF[3])
    ops.fix(11, 1, 1, 1)
    ops.element('TenNodeTetrahedron', 1, *range(1, 11), 1)
    ops.element('zeroLength', 2, 11, 4, '-mat', 2, '-dir', 3)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(4, 0.0, 0.0, 1.0)
    ops.constraints('Plain')
    ops.numberer('Plain')
    ops.system('FullGeneral')
    ops.integrator('LoadControl', 1.0)


@pytest.mark.parametrize('initial', [False, True], ids=['Linear', 'Linear-initial'])
def test_deactivated_element_contributes_no_initial_stiffness(initial):
    """`-initial` used to solve with the full Ki of a switched-off element."""
    _deactivation_model()
    ops.parameter(1, 'element', 1, 'update')
    ops.updateParameter(1, 0)
    if initial:
        ops.algorithm('Linear', '-initial')
    else:
        ops.algorithm('Linear')
    ops.analysis('Static')
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(4, 3) == pytest.approx(1.0 / 100.0, rel=1e-12)


# --------------------------------------------------------------------------
# 4. Diagnostics, printing, responses
# --------------------------------------------------------------------------

def _child(body):
    """Run `body` in a fresh interpreter on this engine; (returncode, output).

    Used where the pre-WP-198 failure is a crash (the interpreter carrying an
    in-process assert would die with it) or an `opserr` line.
    """
    script = ('from _testbed import ops\n'
              'ops.wipe()\n'
              "ops.model('basic', '-ndm', 3, '-ndf', 3)\n"
              "ops.nDMaterial('ElasticIsotropic', 1, 1000.0, 0.25)\n" + body)
    argv, env = pinned_child(['-u', '-c', script])
    p = subprocess.run(argv, env=env, cwd=os.path.dirname(os.path.abspath(__file__)),
                       stdin=subprocess.DEVNULL, capture_output=True, text=True,
                       encoding='utf-8', errors='replace', timeout=120)
    return p.returncode, p.stdout + p.stderr


def _nodes_script(pts):
    return ''.join("ops.node(%d, %r, %r, %r)\n" % ((t,) + tuple(xyz))
                   for t, xyz in enumerate(pts, start=1))


def test_inverted_element_warns():
    """Left-handed corners used to be accepted with negative volume."""
    pts = [_REF[0], _REF[2], _REF[1], _REF[3]]
    pts += [_mid(pts[a], pts[b]) for a, b in _EDGES_GMSH]
    rc, out = _child(_nodes_script(pts)
                     + "ops.element('TenNodeTetrahedron', 1, *range(1, 11), 1)\n"
                     + "print('CHILD-DONE')\n")
    assert rc == 0 and 'CHILD-DONE' in out, out
    assert 'non-positive Jacobian' in out, out


def test_right_handed_element_does_not_warn():
    pts = list(_REF) + [_mid(_REF[a], _REF[b]) for a, b in _EDGES_GMSH]
    rc, out = _child(_nodes_script(pts)
                     + "ops.element('TenNodeTetrahedron', 1, *range(1, 11), 1)\n"
                     + "print('CHILD-DONE')\n")
    assert rc == 0 and 'CHILD-DONE' in out, out
    assert 'non-positive Jacobian' not in out, out


def test_print_json_lists_all_ten_nodes(tmp_path):
    """`print -JSON` wrote [n1, n2, n4]."""
    _single('TenNodeTetrahedron')
    f = tmp_path / 'model.json'
    ops.printModel('-JSON', '-file', str(f))
    ops.wipe()
    data = json.loads(f.read_text())
    eles = data['StructuralAnalysisModel']['geometry']['elements']
    tet = [e for e in eles if e['type'] == 'TenNodeTetrahedron']
    assert len(tet) == 1
    assert tet[0]['nodes'] == list(range(1, 11))


def test_material_response_without_point_number_is_refused():
    """`eleResponse $e material` read argv[1] past the end (access violation)."""
    pts = list(_REF) + [_mid(_REF[a], _REF[b]) for a, b in _EDGES_GMSH]
    rc, out = _child(_nodes_script(pts)
                     + "ops.element('TenNodeTetrahedron', 1, *range(1, 11), 1)\n"
                     + "print('NOPOINT', list(ops.eleResponse(1, 'material')))\n"
                     + "print('POINT4', len(ops.eleResponse(1, 'material', 4, 'stress')))\n")
    assert rc == 0, out
    assert 'NOPOINT []' in out, out
    assert 'POINT4 6' in out, out
