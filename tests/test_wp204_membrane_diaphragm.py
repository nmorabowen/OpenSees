"""WP-204 (+ WP-208): LadrunoMembrane as a real diaphragm.

ADR-98 condensed the in-plane CST stiffness onto the OUTER polygon only, so
hole-edge nodes and interior support nodes (columns, walls inside the panel)
got ZERO in-plane stiffness (YOO probe t1: interior column ux = 0;
membrane_raw2: no in-plane reaction at hole/support nodes). WP-204 condenses an
LST mesh exactly onto every external node (outer loop, holes, -supportNodes,
-diaphragmNodes), adds -f11/-f22/-f12 and membrane-force outputs. WP-208 makes
the slab edge straight between consecutive element nodes on a loop or frame
line (the linear-edge elements it connects to), which WP-204 had left free:
that failed the constant-strain patch test (tests/test_wp208_*).

Oracle: a structured ASDShellQ4 membrane mesh of spacing h, with the
membrane's loops carrying every reference boundary node (``boundary_step = h``)
and ``-meshSize h``: the same connection semantics, as ETABS slab and wall
meshes share their nodes. The remaining difference is Q4 (EAS) vs LST.
Interior point connections (a column) are mesh dependent in both (2-D point
singularity), hence matched resolution.
"""
import numpy as np
import pytest

from _testbed import ops
from _testbed.membrane_diaphragm_ref import Panel, rel_err

pytestmark = [pytest.mark.zone_a]

H = 0.25           # reference spacing = boundary node spacing = -meshSize
MESH = ('-meshSize', H)


def _deep_beam(h=H, column=True):
    """8 x 6 panel, opening (2,1.5)-(4,3), wall line x = 0 (fixed), optional
    interior point column (6, 4.5) (fixed), uniform shear q = 10 (+ axial 2)
    on the x = 8 edge as tributary nodal forces."""
    p = Panel(8.0, 6.0, hole=(2.0, 1.5, 4.0, 3.0), extra=[(6.0, 4.5)], boundary_step=h)
    pts = np.array(p.points)
    col = 8
    wall = [k for k in range(len(pts)) if abs(pts[k, 0]) < 1e-9]
    edge = [k for k in range(len(pts)) if abs(pts[k, 0] - 8.0) < 1e-9]
    fix = {k: (1, 1) for k in wall + ([col] if column else [])}
    loads = {}
    for k in edge:
        trib = 0.5 * h if pts[k, 1] in (0.0, 6.0) else h
        loads[k] = (2.0 * trib, 10.0 * trib)
    return p, fix, loads, wall, edge, col


def test_deep_beam_matches_shell_mesh():
    """wall-supported deep beam with an opening (no point connection)."""
    p, fix, loads, wall, edge, col = _deep_beam(column=False)
    um, rm = p.membrane(loads, fix, extra_args=MESH)
    ur, rr = p.reference(loads, fix, h=H)
    assert rel_err(um[edge], ur[edge]) < 0.02                # 0.5 % measured (1.3 % at h = 0.5)
    for c in (0, 1):
        assert um[edge, c].mean() == pytest.approx(ur[edge, c].mean(), rel=0.01)
    assert rel_err(rm[wall], rr[wall]) < 0.02                # 1.2 % measured
    # the hole corners and the free interior node move with the panel
    assert np.abs(um[4:9]).max() > 0.1 * np.abs(um[edge]).max()


def test_deep_beam_with_point_column_splits_the_shear():
    """add the interior point column: the shear split between wall and column
    matches; the drift does not converge in EITHER model (2-D point support:
    the reference uy moves 15 % from h 0.5 to 0.125), LST/Q4 differ by ~7 %."""
    p, fix, loads, wall, edge, col = _deep_beam()
    um, rm = p.membrane(loads, fix, extra_args=MESH)
    ur, rr = p.reference(loads, fix, h=H)
    sup = wall + [col]
    assert rm[sup].sum(axis=0) == pytest.approx(rr[sup].sum(axis=0), rel=1e-9, abs=1e-9)
    assert rel_err(rm[sup], rr[sup]) < 0.03                 # 2.0 % measured
    assert rm[col, 1] == pytest.approx(rr[col, 1], rel=0.03)  # the column takes its share
    assert um[edge, 1].mean() == pytest.approx(ur[edge, 1].mean(), rel=0.10)  # 7.0 %


def test_core_around_shaft_takes_the_diaphragm_shear():
    """12 x 12 slab around a 4 x 4 shaft whose walls (the hole loop) are fixed,
    plus 4 corner columns; floor force on the perimeter nodes. The stock
    element passed NOTHING to the core."""
    p = Panel(12.0, 12.0, hole=(4.0, 4.0, 8.0, 8.0), boundary_step=H)
    hole = list(p.hole_loop)
    cols = [0, 1, 2, 3]
    per = [k for k in p.outer_loop if k not in cols]
    fix = {k: (1, 1) for k in hole + cols}
    loads = {k: (10.0 * H, 4.0 * H) for k in per}
    um, rm = p.membrane(loads, fix, extra_args=MESH)
    ur, rr = p.reference(loads, fix, h=H)
    total = np.array([10.0, 4.0]) * H * len(per)
    assert (rm[hole].sum(axis=0) + rm[cols].sum(axis=0)) == pytest.approx(-total, rel=1e-9)
    share_m = rm[hole, 0].sum() / -total[0]
    share_r = rr[hole, 0].sum() / -total[0]
    assert share_m > 0.5
    assert share_m == pytest.approx(share_r, rel=0.01)       # 0.3 % measured
    assert rel_err(rm[cols], rr[cols]) < 0.04                # 2.5 % measured
    assert um[per, 0].mean() == pytest.approx(ur[per, 0].mean(), rel=0.01)  # 0.2 %


@pytest.mark.parametrize("mods", [(0.5, 0.5, 0.5), (0.5, 0.3, 0.25), (1.0, 0.35, 0.7)])
def test_inplane_modifiers_match_shellmodifier(mods):
    f11, f22, f12 = mods
    p, fix, loads, wall, edge, col = _deep_beam(column=False)
    mod = ('-f11', f11, '-f22', f22, '-f12', f12)
    um, rm = p.membrane(loads, fix, extra_args=MESH + mod)

    def section():
        ops.section('ElasticMembranePlateSection', 2, p.E, p.nu, p.t, 0.0)
        ops.section('LadrunoShellModifier', 1, 2, *mod)
    ur, rr = p.reference(loads, fix, h=H, section=section)
    assert rel_err(um[edge], ur[edge]) < 0.02
    assert rel_err(rm[wall], rr[wall]) < 0.02
    u1, _ = p.membrane(loads, fix, extra_args=MESH)
    assert rel_err(um[edge], u1[edge]) > 0.2                 # and the modifiers matter


def test_unit_modifiers_change_nothing():
    p, fix, loads, *_ = _deep_beam(h=0.5)
    u0, _ = p.membrane(loads, fix)
    u1, _ = p.membrane(loads, fix, extra_args=('-f11', 1.0, '-f22', 1.0, '-f12', 1.0))
    assert np.array_equal(u0, u1)


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
