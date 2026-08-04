"""Regression: TenNodeTetrahedron::getResponse wrote 24 doubles into a 6-double buffer.

Upstream bug (present in OpenSees/OpenSees master, fixed in this fork):
`getResponse` declared `static Vector stresses(6)` — copy-pasted from
FourNodeTetrahedron, which has a SINGLE Gauss point so 6 is correct there —
but both the "stresses" (responseID 3) and "strains" (responseID 4) branches
loop over all NumGaussPoints=4 points writing 6 components each, i.e. 24
doubles.  `Vector::operator()` bounds-checks only under _G3DEBUG, so release
builds silently wrote 18 doubles (144 bytes) past the end of that static heap
block on every recorder step -> heap corruption / segfault when recording
tet10 element fields.

Two independent observables, both gated here:

  * eleResponse returns 6*4 = 24 values, not 6.  `Information::setVector` does
    `*theVector = newVector` and `Vector::operator=` REALLOCATES on a size
    mismatch, so the undersized scratch also shrank the response Vector that
    `setResponse` had advertised as Vector(6*4) — visible without any crash.
  * an Element recorder emits 24 columns per element.  ElementRecorder sizes
    its columns from the advertised size at setup but copies eleData.Size()
    per element at record time, so the shrink desynchronised the column layout
    across elements.

The model is a single tet10 with a uniform-strain displacement field imposed
on all 10 nodes.  For a quadratic tet a linear displacement field is exactly
representable, so every Gauss point must report the SAME state — which makes
the 4 GP blocks individually checkable against the analytic answer rather than
merely non-empty.
"""

import pytest

from _testbed import ops


pytestmark = [pytest.mark.zone_a]

_E = 1.0e7
_NU = 0.0
_EPS = 1.0e-3          # imposed axial strain
_NGP = 4               # TenNodeTetrahedron::NumGaussPoints
_NCOMP = 6             # Voigt {xx, yy, zz, xy, yz, zx}
_PENALTY = 1.0e18

_ELE = 1

# Vertices of the reference tet, then mid-edge nodes in the OpenSees
# TenNodeTetrahedron order: v1 v2 v3 v4, e12 e23 e13 e14 e34 e24.
_VERTS = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
]
_EDGES = [(1, 2), (2, 3), (1, 3), (1, 4), (3, 4), (2, 4)]


def _coords():
    """Node tag -> (x, y, z) for the 10 nodes, in element connectivity order."""
    pts = {i + 1: _VERTS[i] for i in range(4)}
    for k, (a, b) in enumerate(_EDGES, start=5):
        pa, pb = pts[a], pts[b]
        pts[k] = tuple(0.5 * (pa[j] + pb[j]) for j in range(3))
    return pts


def _build_and_solve():
    """One tet10 driven by u = (eps*x, 0, 0) prescribed at every node.

    With nu = 0 the exact answer is a uniform uniaxial state:
    strain = (eps, 0, 0, 0, 0, 0), stress = (E*eps, 0, 0, 0, 0, 0) at all 4 GP.

    The field is imposed with sp-constraints under a Penalty handler rather
    than fix(), so the DOFs stay in the system (a fully eliminated model would
    leave neq = 0) and the element is genuinely assembled and solved.
    """
    pts = _coords()

    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.nDMaterial('ElasticIsotropic', 1, _E, _NU)
    for tag, (x, y, z) in pts.items():
        ops.node(tag, x, y, z)
    ops.element('TenNodeTetrahedron', _ELE, *range(1, 11), 1)

    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for tag, (x, _y, _z) in pts.items():
        ops.sp(tag, 1, _EPS * x)
        ops.sp(tag, 2, 0.0)
        ops.sp(tag, 3, 0.0)

    ops.constraints('Penalty', _PENALTY, _PENALTY)
    ops.numberer('RCM')
    ops.system('BandGen')
    ops.test('NormDispIncr', 1.0e-12, 20, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    assert ops.analyze(1) == 0, 'tet10 uniform-strain patch step failed to converge'


def _gauss_blocks(values):
    """Split a flat response into per-Gauss-point 6-component blocks."""
    return [values[i * _NCOMP:(i + 1) * _NCOMP] for i in range(_NGP)]


@pytest.mark.parametrize('kind', ['stresses', 'strains'])
def test_eleresponse_returns_all_four_gauss_points(kind):
    """The core regression: 24 values, not 6.

    Before the fix this returned 6 (the response Vector was reallocated down
    to the undersized scratch), while 18 doubles were written out of bounds.
    """
    _build_and_solve()
    out = ops.eleResponse(_ELE, kind)
    assert len(out) == _NCOMP * _NGP, (
        f"eleResponse({_ELE}, {kind!r}) returned {len(out)} values; expected "
        f"{_NCOMP * _NGP} (= 6 components x {_NGP} Gauss points, which is what "
        f"setResponse advertises as Vector(6*4))."
    )


def test_all_gauss_points_carry_the_analytic_uniaxial_state():
    """Every GP block must hold the exact answer — not just GP1.

    With the undersized buffer only the first block ever landed inside the
    Vector, so this also fails if the tail blocks are garbage rather than
    simply absent.
    """
    _build_and_solve()

    strains = _gauss_blocks(ops.eleResponse(_ELE, 'strains'))
    stresses = _gauss_blocks(ops.eleResponse(_ELE, 'stresses'))

    for gp in range(_NGP):
        eps = strains[gp]
        sig = stresses[gp]
        assert eps[0] == pytest.approx(_EPS, rel=1e-8), (gp, eps)
        assert sig[0] == pytest.approx(_E * _EPS, rel=1e-8), (gp, sig)
        # Off-axis components are zero only up to the penalty solve's round-off
        # (measured ~1e-12 relative to the axial component), so bound them
        # against the axial magnitude rather than against absolute zero.
        for comp in range(1, _NCOMP):
            assert abs(eps[comp]) <= 1e-9 * _EPS, (gp, comp, eps)
            assert abs(sig[comp]) <= 1e-9 * _E * _EPS, (gp, comp, sig)


def test_element_recorder_emits_all_24_columns(tmp_path):
    """Gate the ElementRecorder path — the one that actually crashed.

    The recorder both drives getResponse repeatedly (each call re-ran the
    144-byte overrun against the same static block) and lays out its columns
    from the advertised response size, so a shrunk Vector desynchronised the
    file even when the heap damage stayed latent.
    """
    out_file = tmp_path / 'tet10_stresses.out'

    pts = _coords()
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.nDMaterial('ElasticIsotropic', 1, _E, _NU)
    for tag, (x, y, z) in pts.items():
        ops.node(tag, x, y, z)
    ops.element('TenNodeTetrahedron', _ELE, *range(1, 11), 1)

    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for tag, (x, _y, _z) in pts.items():
        ops.sp(tag, 1, _EPS * x)
        ops.sp(tag, 2, 0.0)
        ops.sp(tag, 3, 0.0)

    ops.recorder('Element', '-file', str(out_file), '-time',
                 '-ele', _ELE, 'stresses')

    ops.constraints('Penalty', _PENALTY, _PENALTY)
    ops.numberer('RCM')
    ops.system('BandGen')
    ops.test('NormDispIncr', 1.0e-12, 20, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')

    # Several steps: the overrun was per-record, so repeated recording is what
    # turned latent heap damage into a crash.
    for _ in range(5):
        assert ops.analyze(1) == 0
    ops.remove('recorders')

    rows = [ln.split() for ln in out_file.read_text().splitlines() if ln.strip()]
    assert rows, 'recorder wrote no rows'
    for row in rows:
        # -time prepends one column.
        assert len(row) == 1 + _NCOMP * _NGP, (
            f"recorder row has {len(row)} columns; expected "
            f"{1 + _NCOMP * _NGP} (time + 6 x {_NGP} Gauss points)."
        )
