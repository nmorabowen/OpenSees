"""WP-204: LadrunoMembrane as a real diaphragm.

ADR-98 condensed the in-plane CST stiffness onto the OUTER polygon only and
slaved every other boundary vertex to straight edges, so hole-edge nodes and
interior support nodes (columns, walls inside the panel) got ZERO in-plane
stiffness (YOO probe t1: interior column ux = 0; membrane_raw2: no in-plane
reaction at hole/support nodes). WP-204 condenses an LST mesh exactly onto
every external node (outer loop, holes, -supportNodes, -diaphragmNodes), adds
ETABS-style -f11/-f22/-f12 modifiers and membrane-force outputs.

Oracles:
  * exactness: the same panel rebuilt with explicit 6-node triangles (tri6n)
    on the element's own diaphragm mesh — displacements, reactions and the
    recovered membrane forces to round-off;
  * physics: a refined ASDShellQ4 membrane mesh (tests/_testbed/
    membrane_diaphragm_ref.py), compared on WELL-POSED quantities. A point
    load or a point support is logarithmically singular in 2-D elasticity, so
    displacement AT a point connection never converges with refinement (the
    h = 0.1 -> 0.05 reference moves 8 % on the deep-beam case below) — in
    ETABS as here. The reaction SPLIT between supports converges (0.4-2.6 %
    over -meshSize 1.0 .. 0.25); displacements are compared at matched
    resolution (-meshSize 0.25 vs h = 0.1: 2.1 %).
"""
import numpy as np
import pytest

from _testbed import ops
from _testbed.membrane_diaphragm_ref import Panel, rel_err

pytestmark = [pytest.mark.zone_a]

H_REF = 0.1        # reference ASDShellQ4 grid spacing
MESH = ('-meshSize', 0.25)   # matched resolution for the displacement checks


def _deep_beam(extra_hole_pts=()):
    """8 x 6 panel, opening (2,1.5)-(4,3), a wall line on x = 0 (nodes every
    0.5 m, fixed), an interior point column at (6, 4.5) (fixed), and a uniform
    shear q = 10 on the x = 8 edge (tributary nodal forces every 0.5 m)."""
    ys = [0.5 * i for i in range(1, 12)]
    wall = [(0.0, y) for y in ys]
    edge = [(8.0, y) for y in ys]
    p = Panel(8.0, 6.0, hole=(2.0, 1.5, 4.0, 3.0), extra=[(6.0, 4.5)] + wall + edge)
    col = 8
    iwall = [0, 3] + list(range(9, 9 + len(wall)))
    iedge = [1, 2] + list(range(9 + len(wall), 9 + 2 * len(wall)))
    fix = {k: (1, 1) for k in iwall + [col]}
    loads = {k: (0.0, 5.0) for k in iedge}
    loads[1] = (0.0, 2.5)
    loads[2] = (0.0, 2.5)
    dia = list(range(9, len(p.points)))       # wall/edge nodes: in-plane only
    return p, fix, loads, dia, col, iwall, iedge


def test_reaction_split_converges_to_refined_shell_mesh():
    """wall + interior column share the diaphragm shear like the fine mesh,
    at every -meshSize; the stock element gave the column ZERO."""
    p, fix, loads, dia, col, iwall, _ = _deep_beam()
    _, rr = p.reference(loads, fix, h=H_REF)
    sup = iwall + [col]
    errs = []
    for ms in (1.0, 0.5, 0.25):
        _, rm = p.membrane(loads, fix, extra_args=('-meshSize', ms), diaphragm=dia)
        assert rm[sup].sum(axis=0) == pytest.approx(rr[sup].sum(axis=0), rel=1e-9, abs=1e-9)
        errs.append(rel_err(rm[sup], rr[sup]))
        assert abs(rm[col, 1]) > 0.5 * abs(rr[col, 1])    # the column takes its share
    assert errs[0] < 0.03 and errs[-1] < 0.01, errs
    assert errs[-1] < errs[0]


def test_deep_beam_displacement_matched_resolution():
    p, fix, loads, dia, col, iwall, iedge = _deep_beam()
    ur, rr = p.reference(loads, fix, h=H_REF)
    um, rm = p.membrane(loads, fix, extra_args=MESH, diaphragm=dia)
    assert rel_err(um[iedge], ur[iedge]) < 0.03
    assert rel_err(rm[iwall + [col]], rr[iwall + [col]]) < 0.01
    # the hole corners move with the panel (stock: unconnected, zero stiffness)
    assert np.abs(um[4:8]).max() > 0.1 * np.abs(um[iedge]).max()


def test_core_around_shaft_takes_the_diaphragm_shear():
    """12 x 12 slab around a 4 x 4 shaft whose walls (hole edges, nodes every
    1 m) are fixed, plus 4 corner columns; floor force on the perimeter
    (nodes every 1 m): the core takes its share through the diaphragm, split
    like the fine mesh. The stock element passed NOTHING to the core."""
    hole = (4.0, 4.0, 8.0, 8.0)
    core_pts = ([(x, 4.0) for x in (5.0, 6.0, 7.0)] + [(8.0, y) for y in (5.0, 6.0, 7.0)] +
                [(x, 8.0) for x in (5.0, 6.0, 7.0)] + [(4.0, y) for y in (5.0, 6.0, 7.0)])
    per_pts = ([(float(x), 0.0) for x in range(1, 12)] + [(12.0, float(y)) for y in range(1, 12)] +
               [(float(x), 12.0) for x in range(1, 12)] + [(0.0, float(y)) for y in range(1, 12)])
    p = Panel(12.0, 12.0, hole=hole, extra=core_pts + per_pts)
    core = [4, 5, 6, 7] + list(range(8, 8 + len(core_pts)))
    cols = [0, 1, 2, 3]
    perimeter = list(range(8 + len(core_pts), len(p.points)))
    fix = {k: (1, 1) for k in core + cols}
    loads = {k: (10.0, 4.0) for k in perimeter}
    dia = list(range(8, len(p.points)))
    um, rm = p.membrane(loads, fix, extra_args=MESH, diaphragm=dia)
    ur, rr = p.reference(loads, fix, h=H_REF)
    total = np.array([10.0, 4.0]) * len(perimeter)
    assert (rm[core].sum(axis=0) + rm[cols].sum(axis=0)) == pytest.approx(-total, rel=1e-9)
    share_m = rm[core, 0].sum() / -total[0]
    share_r = rr[core, 0].sum() / -total[0]
    assert share_m > 0.5                                   # the core carries most of it
    assert share_m == pytest.approx(share_r, rel=0.01)
    # per-wall resultants (corner nodes split half to each wall)
    sides = [[4, 5] + list(range(8, 11)), [5, 6] + list(range(11, 14)),
             [6, 7] + list(range(14, 17)), [7, 4] + list(range(17, 20))]
    def walls(r):
        return np.array([r[sd[2:]].sum(axis=0) + 0.5 * (r[sd[0]] + r[sd[1]]) for sd in sides])
    assert rel_err(walls(rm), walls(rr)) < 0.02      # 0.9 % measured
    assert rel_err(rm[cols], rr[cols]) < 0.02        # 0.2 % measured
    assert rel_err(um[perimeter], ur[perimeter]) < 0.03   # mean drift 0.8 %


@pytest.mark.parametrize("mods", [(0.5, 0.5, 0.5), (0.5, 0.3, 0.25), (1.0, 0.35, 0.7)])
def test_inplane_modifiers_match_shellmodifier(mods):
    f11, f22, f12 = mods
    p, fix, loads, dia, col, iwall, iedge = _deep_beam()
    # add a horizontal load component so f11 and f22 both work
    loads = {k: (2.0, v[1]) for k, v in loads.items()}
    mod = ('-f11', f11, '-f22', f22, '-f12', f12)
    um, rm = p.membrane(loads, fix, extra_args=MESH + mod, diaphragm=dia)

    def section():
        ops.section('ElasticMembranePlateSection', 2, p.E, p.nu, p.t, 0.0)
        ops.section('LadrunoShellModifier', 1, 2, *mod)
    ur, rr = p.reference(loads, fix, h=H_REF, section=section)
    assert rel_err(um[iedge], ur[iedge]) < 0.03
    assert rel_err(rm[iwall + [col]], rr[iwall + [col]]) < 0.01
    # and the modifiers matter
    u1, _ = p.membrane(loads, fix, extra_args=MESH, diaphragm=dia)
    assert rel_err(um[iedge], u1[iedge]) > 0.2


def test_unit_modifiers_change_nothing():
    p, fix, loads, dia, *_ = _deep_beam()
    u0, _ = p.membrane(loads, fix, diaphragm=dia)
    u1, _ = p.membrane(loads, fix, extra_args=('-f11', 1.0, '-f22', 1.0, '-f12', 1.0), diaphragm=dia)
    assert np.array_equal(u0, u1)


# ------------------------------------------------------------------ exactness
def _explicit_lst(verts, tris, points, loads, fix, E, nu, t):
    """the panel rebuilt with explicit 6-node triangles on the element's own
    diaphragm mesh (2-D plane stress); retained points by coordinates."""
    ops.wipe(); ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.nDMaterial('ElasticIsotropic', 1, E, nu)
    used = sorted({int(v) for tri in tris for v in tri})
    for v in used:
        ops.node(v + 1, float(verts[v, 0]), float(verts[v, 1]))
    mid = {}
    nxt = len(verts) + 1
    e = 1
    for a, b, c in tris.astype(int).tolist():
        ms = []
        for i, j in ((a, b), (b, c), (c, a)):
            key = (min(i, j), max(i, j))
            if key not in mid:
                mid[key] = nxt
                ops.node(nxt, float(0.5 * (verts[i, 0] + verts[j, 0])), float(0.5 * (verts[i, 1] + verts[j, 1])))
                nxt += 1
            ms.append(mid[key])
        ops.element('tri6n', e, a + 1, b + 1, c + 1, *ms, t, 'PlaneStress', 1)
        e += 1
    tags = []
    for x, y in points:
        d = np.hypot(verts[:, 0] - x, verts[:, 1] - y)
        tags.append(int(np.argmin(d)) + 1)
        assert d.min() < 1e-9
    for k, tg in enumerate(tags):
        fx, fy = fix.get(k, (0, 0))
        if fx or fy:
            ops.fix(tg, fx, fy)
    ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1)
    for k, (fx, fy) in loads.items():
        ops.load(tags[k], fx, fy)
    ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0); ops.algorithm('Linear'); ops.analysis('Static')
    assert ops.analyze(1) == 0
    ops.reactions()
    u = np.array([[ops.nodeDisp(tg, 1), ops.nodeDisp(tg, 2)] for tg in tags])
    r = np.array([[ops.nodeReaction(tg, 1), ops.nodeReaction(tg, 2)] for tg in tags])
    # LST stress is linear: the mean of the three symmetric Gauss points is the
    # centroid value
    n = np.array([t * np.mean(np.reshape(ops.eleResponse(k, 'stresses'), (3, 3)), axis=0)
                  for k in range(1, e)])
    return u, r, n


@pytest.mark.parametrize("mesh", [(), ('-meshSize', 0.6)])
def test_condensation_is_exact_on_its_own_mesh(mesh):
    """Guyan condensation is exact for loads at retained nodes: the element must
    equal the full LST mesh it condensed, to round-off (outer, hole, interior
    support and an in-plane-only diaphragm node all retained)."""
    p = Panel(8.0, 6.0, hole=(2.0, 1.5, 4.0, 3.0), extra=[(6.0, 4.5), (7.0, 1.0)])
    fix = {0: (1, 1), 1: (0, 1), 6: (1, 0)}
    loads = {2: (100.0, 50.0), 4: (-20.0, 35.0), 8: (-60.0, 30.0), 9: (15.0, -40.0)}
    um, rm = p.membrane(loads, fix, extra_args=mesh, diaphragm=[9])
    verts = np.array(ops.eleResponse(100, 'diaphragmVertices')).reshape(-1, 3)[:, :2]
    tris = np.array(ops.eleResponse(100, 'diaphragmTriangles')).reshape(-1, 3)
    nm = np.array(ops.eleResponse(100, 'membraneForces')).reshape(-1, 3)
    mean = np.array(ops.eleResponse(100, 'meanMembraneForce'))
    ue, re, ne = _explicit_lst(verts, tris, p.points, loads, fix, p.E, p.nu, p.t)
    assert rel_err(um, ue) < 1e-9
    assert rel_err(rm, re) < 1e-9
    # force recovery: per-triangle centroid forces and the area mean
    assert nm.shape == ne.shape
    assert rel_err(nm, ne) < 1e-8
    v = verts[tris.astype(int)]
    area = 0.5 * np.abs((v[:, 1, 0] - v[:, 0, 0]) * (v[:, 2, 1] - v[:, 0, 1]) -
                        (v[:, 2, 0] - v[:, 0, 0]) * (v[:, 1, 1] - v[:, 0, 1]))
    assert rel_err(mean, area @ ne / area.sum()) < 1e-8


# --------------------------------------------------------------------- outputs
def test_diaphragm_forces_and_virial_identity():
    """nodal diaphragm forces equilibrate the loads, and the panel integral of
    N equals the virial of the nodal forces (exact for an FE field containing
    linear displacements): int N12 dA = 1/2 sum(f_x y + f_y x)."""
    ops.wipe(); ops.model('basic', '-ndm', 3, '-ndf', 6)
    pts = [(0.0, 0.0), (6.0, 0.0), (6.0, 4.0), (0.0, 4.0), (1.0, 1.0), (3.0, 1.0), (3.0, 2.5),
           (1.0, 2.5), (4.5, 3.0)]
    for k, (x, y) in enumerate(pts, start=1):
        ops.node(k, x, y, 0.0)
        ops.fix(k, int(k == 1), int(k in (1, 2)), 1, 1, 1, 1)
    ops.element('LadrunoMembrane', 100, '-nodes', 1, 2, 3, 4, '-hole', 5, 6, 7, 8,
                '-diaphragmNodes', 9, '-stiff', 25.0e6, 0.2, 0.2)
    ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1)
    ops.load(3, 10.0, 40.0, 0, 0, 0, 0)
    ops.load(9, -25.0, 15.0, 0, 0, 0, 0)
    ops.load(6, 5.0, -8.0, 0, 0, 0, 0)
    ops.system('UmfPack'); ops.numberer('Plain'); ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0); ops.algorithm('Linear'); ops.analysis('Static')
    assert ops.analyze(1) == 0
    f = np.array(ops.eleResponse(100, 'diaphragmForces')).reshape(-1, 3)
    assert np.allclose(f.sum(axis=0), 0.0, atol=1e-8)        # self-equilibrated
    assert f[2, :2] == pytest.approx([10.0, 40.0], rel=1e-8)  # free loaded node
    area = 6.0 * 4.0 - 2.0 * 1.5
    mean = np.array(ops.eleResponse(100, 'meanMembraneForce'))
    xy = np.array(pts)
    virial = np.array([np.sum(f[:, 0] * xy[:, 0]), np.sum(f[:, 1] * xy[:, 1]),
                       0.5 * np.sum(f[:, 0] * xy[:, 1] + f[:, 1] * xy[:, 0])])
    assert mean * area == pytest.approx(virial, rel=1e-8, abs=1e-8)


def test_save_restore_keeps_diaphragm_nodes_and_modifiers(tmp_path):
    p = Panel(8.0, 6.0, hole=(2.0, 1.5, 4.0, 3.0), extra=[(6.0, 4.5)])
    fix = {0: (1, 1), 1: (0, 1)}
    loads = {2: (100.0, 50.0), 8: (-60.0, 30.0)}
    args = ('-f11', 0.5, '-f22', 0.7, '-f12', 0.4)
    ref, _ = p.membrane(loads, fix, extra_args=args, diaphragm=[8])
    p.membrane(loads, fix, extra_args=args, diaphragm=[8])
    ops.database('File', str(tmp_path / 'mem'))
    ops.save(3)
    ops.restore(3)
    ops.analyze(1)
    got = np.array([[ops.nodeDisp(t, 1), ops.nodeDisp(t, 2)] for t in range(1, 10)])
    assert rel_err(got, 2.0 * ref) < 1e-9   # restored state + one more load step
