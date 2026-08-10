"""Regression: BezierTet10 symmetrized an unsymmetric plastic tangent.

BezierTet10::formCore (and getInitialStiff, and both BezierTri6 stiffness
loops) assembled only the upper triangle of B'DB and mirrored it into the
lower triangle.  The mirror writes the B'(I)DB(J) contraction into K(J,I),
whose true value is the D-transpose contraction — exact while D is symmetric,
silently wrong the moment any Gauss point yields on a material with
non-associated flow (DruckerPrager with rho_bar != rho, ManzariDafalias
always).  On a mesh where some Gauss points yield and others stay elastic,
Newton iterates against the corrupted tangent and stalls with no error from
any element or material (TIMs campaign defect report, item 1: elastic OK,
sigma_y large OK, first yield FAILS; TenNodeTetrahedron — which assembles the
full product — OK throughout).

Gates, mirroring the reporters' discriminator table:

  * mixed-yield convergence: a shear-plus-compression cube of BezierTet10
    (std and -bbar) on non-associated DruckerPrager must converge past first
    yield.  Before the fix this fails; TenNodeTetrahedron on the identical
    deck converges (kept here as the deck-feasibility control).
  * the bug's fingerprint: eleResponse 'stiffness' on a yielded element must
    be measurably UNSYMMETRIC.  The old code returned an exactly symmetric
    matrix by construction, so max|K - K'| == 0 whatever the material did.
  * elastic no-regression: with ElasticIsotropic the rewritten full-product
    loop must stay symmetric to round-off, and the corner-node displacements
    must match TenNodeTetrahedron (straight-sided quadratic Bezier and
    Lagrange tets span the same space, and quadratic Bernstein functions are
    interpolatory at vertices, so the corner values of the Galerkin solution
    coincide).
"""

import pytest

from _testbed import ops


pytestmark = [pytest.mark.zone_a]

# ---- geometry: unit cube, Kuhn 6-tet split, quadratic nodes --------------

_VERTS = {
    1: (0.0, 0.0, 0.0), 2: (1.0, 0.0, 0.0), 3: (1.0, 1.0, 0.0), 4: (0.0, 1.0, 0.0),
    5: (0.0, 0.0, 1.0), 6: (1.0, 0.0, 1.0), 7: (1.0, 1.0, 1.0), 8: (0.0, 1.0, 1.0),
}
# All six tets share the body diagonal 1-7 and are positively oriented for
# TenNodeTetrahedron's signed-Jacobian convention (checked analytically).
_TETS = [(1, 2, 3, 7), (1, 3, 4, 7), (1, 4, 8, 7), (1, 8, 5, 7), (1, 5, 6, 7), (1, 6, 2, 7)]
# Mid-edge local order shared by TenNodeTetrahedron and BezierTet10:
# e12 e23 e13 e14 e34 e24 (0-based local vertex pairs).
_EDGE_ORDER = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]

_NELD = 30      # 10 nodes x 3 DOF
_NELEM = len(_TETS)

# ---- material: kPa-scale DruckerPrager (PEER example constants) ----------

_K, _G = 27777.78, 9259.26
_RHO = 0.398
_RHO_BAR = 0.1          # rho_bar != rho -> NON-associated -> unsymmetric Cep
_SIG_Y = 5.0            # yields under the applied load
_SIG_Y_ELASTIC = 500.0  # never yields (reporters' "OK" row)
# Linear isotropic hardening (theta = 1): perfectly plastic corner flow makes
# even the TenNodeTetrahedron control limp near collapse (a DruckerPrager
# apex/corner-return trait, nothing to do with the element defect); modest
# hardening keeps the TRUE-tangent Newton quadratic while Cep stays
# unsymmetric — the asymmetry comes from rho_bar != rho, not from H.
_HARD_H = 1000.0

# Total top-corner loads, calibrated on the TenNodeTetrahedron control:
# yield onset is at H ~ 0.9 and a genuine flip-flop stall (material-level,
# element-independent) appears at H ~ 1.5, so the ramp ends at H = 1.35 —
# comfortably inside the mixed elastic/plastic regime the reporters
# identified, short of the deck's own limit.  Compression raises the DP
# margin; the corner concentration keeps the yield partial.
_H_TOTAL = 1.35
_V_TOTAL = 2.0
_NSTEPS = 16


def _mesh():
    """(node tag -> coords, per-element 10-tag connectivity)."""
    nodes = dict(_VERTS)
    mid = {}
    next_tag = 9
    conn = []
    for tet in _TETS:
        elem = list(tet)
        for a, b in _EDGE_ORDER:
            key = tuple(sorted((tet[a], tet[b])))
            if key not in mid:
                pa, pb = nodes[key[0]], nodes[key[1]]
                nodes[next_tag] = tuple(0.5 * (pa[j] + pb[j]) for j in range(3))
                mid[key] = next_tag
                next_tag += 1
            elem.append(mid[key])
        conn.append(elem)
    return nodes, conn


def _build(ele_type, sig_y=None, bbar=False, elastic=False):
    """Cube model: base fixed, shear + compression at the four top corners."""
    nodes, conn = _mesh()

    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)

    if elastic:
        _E = 9.0 * _K * _G / (3.0 * _K + _G)
        _NU = (3.0 * _K - 2.0 * _G) / (2.0 * (3.0 * _K + _G))
        ops.nDMaterial('ElasticIsotropic', 1, _E, _NU)
    else:
        ops.nDMaterial('DruckerPrager', 1, _K, _G, sig_y, _RHO, _RHO_BAR,
                       0.0, 0.0, 0.0, 0.0, _HARD_H, 1.0, 0.0)

    for tag, (x, y, z) in nodes.items():
        ops.node(tag, x, y, z)
    for etag, elem in enumerate(conn, start=1):
        args = ['-bbar'] if bbar else []
        ops.element(ele_type, etag, *elem, 1, *args)

    for tag, (_x, _y, z) in nodes.items():
        if z == 0.0:
            ops.fix(tag, 1, 1, 1)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for tag in (5, 6, 7, 8):          # top corners
        ops.load(tag, _H_TOTAL / 4.0, 0.0, -_V_TOTAL / 4.0)

    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGen')             # unsymmetric-capable LU
    ops.test('NormDispIncr', 1.0e-8, 50, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / _NSTEPS)
    ops.analysis('Static')


def _max_rel_asym():
    """max over elements of max|K - K'| / max|K| from eleResponse 'stiffness'."""
    worst = 0.0
    for etag in range(1, _NELEM + 1):
        flat = ops.eleResponse(etag, 'stiffness')
        assert len(flat) == _NELD * _NELD, (etag, len(flat))
        scale = max(abs(v) for v in flat)
        dev = max(
            abs(flat[i * _NELD + j] - flat[j * _NELD + i])
            for i in range(_NELD) for j in range(i)
        )
        if scale > 0.0:
            worst = max(worst, dev / scale)
    return worst


def _run(ele_type, sig_y, bbar=False, probe=True):
    """Full ramp; returns (loaded-corner displacement, max relative asymmetry).

    probe=False skips the stiffness probe (TenNodeTetrahedron registers no
    'stiffness' response — the control only needs to converge).
    """
    _build(ele_type, sig_y=sig_y, bbar=bbar)
    for step in range(_NSTEPS):
        assert ops.analyze(1) == 0, (
            f'{ele_type} (bbar={bbar}, sigY={sig_y}) failed to converge at '
            f'step {step + 1}/{_NSTEPS} — the mixed elastic/plastic tangent '
            f'regime is exactly where the symmetrized-B\'DB defect bites'
        )
    ux = ops.nodeDisp(7, 1)
    return ux, (_max_rel_asym() if probe else None)


# ---- the regression ------------------------------------------------------

@pytest.mark.parametrize('bbar', [False, True], ids=['std', 'bbar'])
def test_mixed_yield_converges_and_tangent_goes_unsymmetric(bbar):
    """Past first yield: converge, AND report a genuinely unsymmetric K.

    The second assertion is the bug's fingerprint: the old mirror produced
    max|K - K'| == 0 identically, so it doubles as proof that the deck really
    entered the non-associated plastic regime rather than staying elastic.
    """
    ux_plastic, asym_plastic = _run('BezierTet10', _SIG_Y, bbar=bbar)
    ux_elastic, asym_elastic = _run('BezierTet10', _SIG_Y_ELASTIC, bbar=bbar)

    assert asym_elastic < 1.0e-9, (
        'never-yielding run must keep a symmetric tangent', asym_elastic)
    assert asym_plastic > 1.0e-8, (
        'yielded run must expose the unsymmetric consistent tangent '
        '(exactly 0 under the pre-fix mirror)', asym_plastic)
    assert abs(ux_plastic - ux_elastic) > 5.0e-3 * abs(ux_elastic), (
        'plastic and never-yield displacements coincide — the load never '
        'caused yield and the test lost its discriminating power',
        ux_plastic, ux_elastic)


def test_control_lagrange_tet10_converges():
    """Deck-feasibility control — the reporters' passing column."""
    ux, _ = _run('TenNodeTetrahedron', _SIG_Y, probe=False)
    assert ux != 0.0


def test_elastic_stiffness_stays_symmetric():
    """The rewritten full-product loop must not break the symmetric case."""
    _build('BezierTet10', elastic=True)
    assert ops.analyze(1) == 0
    assert _max_rel_asym() < 1.0e-10


def test_elastic_corner_displacements_match_lagrange_tet10():
    """Straight-sided Bezier and Lagrange tet10 give one Galerkin solution.

    Corner (vertex) nodes only: quadratic Bernstein functions are
    interpolatory there, so corner DOFs are field values in both elements;
    mid-edge Bezier DOFs are control values and legitimately differ.
    """
    corner_disp = {}
    for ele_type in ('BezierTet10', 'TenNodeTetrahedron'):
        _build(ele_type, elastic=True)
        for _ in range(_NSTEPS):
            assert ops.analyze(1) == 0
        corner_disp[ele_type] = {
            tag: [ops.nodeDisp(tag, d) for d in (1, 2, 3)]
            for tag in (5, 6, 7, 8)
        }

    for tag in (5, 6, 7, 8):
        for d in range(3):
            a = corner_disp['BezierTet10'][tag][d]
            b = corner_disp['TenNodeTetrahedron'][tag][d]
            assert a == pytest.approx(b, rel=1.0e-7, abs=1.0e-12), (tag, d, a, b)
