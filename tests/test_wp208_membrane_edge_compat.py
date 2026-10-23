"""WP-208: LadrunoMembrane's slab edge is straight between consecutive element
nodes (the constant-strain patch test).

WP-204 left every mesh node on a loop segment between two element nodes free
(LST midsides, and mesh vertices when -meshSize is finer than the node
spacing), so the slab edge bowed between connection nodes. Reported by the
apeGmsh SLE port (PR #2495 xfail): a 2 x 2 panel with 16 boundary nodes every
0.5 in -nodes, ux fixed at x = 0 and ux = eps S imposed at x = S gave an edge
reaction 0.740 of E t eps S (0.624 at -meshSize 0.25), ux(0.5, 0) 1.30 x the
linear field and uy(2, 2) = -3.7e-5 instead of 0 at nu = 0.

Oracle: the exact constant-strain field. Every mesh node strictly inside a
connection line (loop segment split at element nodes, or a supporting frame
whose ends are element nodes) is slaved linearly to the line's two ends; the
interior is then exactly linear, so reactions, displacements AND the recovered
membrane forces are exact to round-off.
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

S, E, T = 2.0, 30.0e6, 0.25


def _loop(step):
    n = int(round(S / step))
    pts = [(i * step, 0.0) for i in range(n)] + [(S, j * step) for j in range(n)]
    pts += [(S - i * step, S) for i in range(n)] + [(0.0, S - j * step) for j in range(n)]
    return pts


def _apgmsh_case(nu, mesh, step=0.5):
    """the apeGmsh reproducer: uniaxial pull, edges every `step` in -nodes."""
    ops.wipe(); ops.model('basic', '-ndm', 3, '-ndf', 6)
    pts = _loop(step)
    for k, (x, y) in enumerate(pts, start=1):
        ops.node(k, x, y, 0.0)
        ops.fix(k, int(abs(x) < 1e-12), int(k == 1), 1, 1, 1, 1)
    ops.element('LadrunoMembrane', 100, '-nodes', *range(1, len(pts) + 1),
                '-stiff', E, nu, T, '-meshSize', mesh)
    eps = 1.0e-4
    ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1)
    right = [k for k, (x, _) in enumerate(pts, start=1) if abs(x - S) < 1e-12]
    for k in right:
        ops.sp(k, 1, eps * S)
    ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints('Transformation')
    ops.integrator('LoadControl', 1.0); ops.algorithm('Linear'); ops.analysis('Static')
    assert ops.analyze(1) == 0
    ops.reactions()
    return pts, right, eps


@pytest.mark.parametrize("nu", [0.0, 0.2])
@pytest.mark.parametrize("mesh", [0.5, 0.25])
def test_apegmsh_uniaxial_patch(nu, mesh):
    pts, right, eps = _apgmsh_case(nu, mesh)
    R = sum(ops.nodeReaction(k, 1) for k in right)
    assert R == pytest.approx(E * T * eps * S, rel=1e-9)     # was 0.740 / 0.624
    for k, (x, y) in enumerate(pts, start=1):
        assert ops.nodeDisp(k, 1) == pytest.approx(eps * x, abs=1e-12)
        assert ops.nodeDisp(k, 2) == pytest.approx(-nu * eps * y, abs=1e-12)
    n = np.array(ops.eleResponse(100, 'membraneForces')).reshape(-1, 3)
    assert np.allclose(n, [E * T * eps, 0.0, 0.0], rtol=1e-8, atol=1e-6)


def _general_patch(eps, mesh, mods=()):
    """hole + interior support + interior diaphragm node + an edge diaphragm
    node + an interior supporting frame, every element node driven by the
    linear field u = [exx x + gxy/2 y, gxy/2 x + eyy y]."""
    ops.wipe(); ops.model('basic', '-ndm', 3, '-ndf', 6)
    outer = _loop(0.5)
    hole = [(0.5, 0.5), (1.0, 0.5), (1.0, 1.0), (0.5, 1.0)]
    extra = [(1.5, 1.5), (1.5, 0.75), (0.25, 0.0)]   # support, diaphragm, on-edge diaphragm
    pts = outer + hole + extra
    for k, (x, y) in enumerate(pts, start=1):
        ops.node(k, x, y, 0.0)
        ops.fix(k, 0, 0, 1, 1, 1, 1)
    no, nh = len(outer), len(hole)
    # an interior beam from (1.5, 1.5) to (1.5, 0.75) carrying the slab
    ops.geomTransf('Linear', 1, 1.0, 0.0, 0.0)
    ops.element('elasticBeamColumn', 50, no + nh + 1, no + nh + 2, 1e-6, 1.0, 1.0, 1e-6, 1e-6, 1e-6, 1)
    ops.element('LadrunoMembrane', 100, '-nodes', *range(1, no + 1),
                '-hole', *range(no + 1, no + nh + 1), '-supportNodes', no + nh + 1,
                '-diaphragmNodes', no + nh + 2, no + nh + 3, '-frames', 50,
                '-stiff', E, 0.2, T, '-meshSize', mesh, *mods)
    exx, eyy, gxy = eps
    ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1)
    for k, (x, y) in enumerate(pts, start=1):
        ops.sp(k, 1, exx * x + 0.5 * gxy * y)
        ops.sp(k, 2, 0.5 * gxy * x + eyy * y)
    ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints('Transformation')
    ops.integrator('LoadControl', 1.0); ops.algorithm('Linear'); ops.analysis('Static')
    assert ops.analyze(1) == 0


@pytest.mark.parametrize("mesh", [0.5, 0.2])
@pytest.mark.parametrize("mods", [(), ('-f11', 0.5, '-f22', 0.3, '-f12', 0.4)])
def test_general_patch_forces_exact(mesh, mods):
    eps = (1.0e-4, -3.0e-5, 6.0e-5)
    _general_patch(eps, mesh, mods)
    f = dict(zip(mods[0::2], mods[1::2]))
    s = np.sqrt([f.get('-f11', 1.0), f.get('-f22', 1.0), f.get('-f12', 1.0)])
    c = E / (1.0 - 0.2 ** 2)
    D = np.array([[c, 0.2 * c, 0.0], [0.2 * c, c, 0.0], [0.0, 0.0, 0.4 * c]])
    expect = T * (np.outer(s, s) * D) @ np.array(eps)
    n = np.array(ops.eleResponse(100, 'membraneForces')).reshape(-1, 3)
    assert n.shape[0] > 8
    assert np.allclose(n, expect, rtol=1e-7, atol=1e-7 * np.abs(expect).max())
    assert ops.eleResponse(100, 'meanMembraneForce') == pytest.approx(list(expect), rel=1e-7)

