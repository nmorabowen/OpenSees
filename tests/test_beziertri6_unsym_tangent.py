"""Regression: BezierTri6 symmetrized an unsymmetric plastic tangent.

Companion to test_beziertet10_unsym_tangent.py — commit acc9cd53d fixed the
identical upper-triangle-and-mirror BᵀDB assembly in BOTH BezierTri6
stiffness loops (getTangentStiff and getInitialStiff), but the tet test
cannot gate the 2D element; the adversarial review flagged that a regression
re-introducing the mirror in BezierTri6 alone would pass the whole suite.

This file carries the watertight gate: `eleResponse 'stiffness'` on a
yielded element must be measurably UNSYMMETRIC (non-associated
plane-strain DruckerPrager), while a never-yielding run stays symmetric to
round-off.  Under the pre-fix mirror the reported matrix satisfies
max|K - K'| == 0 bit-exactly, whatever the material does, so the plastic
assertion fails with certainty on a regressed binary.

Deck: unit square, two plane-strain BezierTri6 (diagonal 1-3), base fixed,
shear + compression at the top corners — a scaled-down 2D version of the
tet test's mixed elastic/plastic regime.
"""

import pytest

from _testbed import ops


pytestmark = [pytest.mark.zone_a]

# corners 1-4, mid-edge nodes 5-9; element order: c1 c2 c3 mid12 mid23 mid31
_NODES = {
    1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0),
    5: (0.5, 0.0),   # edge 1-2
    6: (1.0, 0.5),   # edge 2-3
    7: (0.5, 0.5),   # edge 1-3 (shared diagonal)
    8: (0.5, 1.0),   # edge 3-4
    9: (0.0, 0.5),   # edge 1-4
}
_CONN = [
    (1, 2, 3, 5, 6, 7),   # CCW, straight-sided
    (1, 3, 4, 7, 8, 9),
]

_NELD = 12       # 6 nodes x 2 DOF
_THICK = 1.0

# Same kPa-scale non-associated DruckerPrager as the tet test; the material
# is created 3D and the element requests the 'PlaneStrain' copy.
_K, _G = 27777.78, 9259.26
_RHO, _RHO_BAR = 0.398, 0.1
_SIG_Y = 5.0
_SIG_Y_ELASTIC = 500.0
_HARD_H = 1000.0

# Per-variant lateral load, calibrated by sweeping H under THIS proportional
# ramp (H and V reach full value together): std first yields in H = (1.5,
# 1.6] and still converges at 1.7; -bbar yields in (1.0, 1.1] and hits the
# material-level flip-flop stall (as in the tet test) in (1.5, 1.6].  Each
# ramp ends inside its own mixed elastic/plastic window.
_H_TOTAL = {False: 1.7, True: 1.3}
_V_TOTAL = 1.6
_NSTEPS = 12


def _build(sig_y, bbar=False):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.nDMaterial('DruckerPrager', 1, _K, _G, sig_y, _RHO, _RHO_BAR,
                   0.0, 0.0, 0.0, 0.0, _HARD_H, 1.0, 0.0)

    for tag, (x, y) in _NODES.items():
        ops.node(tag, x, y)
    for etag, elem in enumerate(_CONN, start=1):
        args = ['-bbar'] if bbar else []
        ops.element('BezierTri6', etag, *elem, _THICK, 'PlaneStrain', 1, *args)

    for tag, (_x, y) in _NODES.items():
        if y == 0.0:
            ops.fix(tag, 1, 1)

    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for tag in (3, 4):
        ops.load(tag, _H_TOTAL[bbar] / 2.0, -_V_TOTAL / 2.0)

    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandGen')
    ops.test('NormDispIncr', 1.0e-8, 50, 0)
    ops.algorithm('Newton')
    ops.integrator('LoadControl', 1.0 / _NSTEPS)
    ops.analysis('Static')


def _max_rel_asym():
    worst = 0.0
    for etag in (1, 2):
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


def _run(sig_y, bbar=False):
    _build(sig_y, bbar=bbar)
    for step in range(_NSTEPS):
        assert ops.analyze(1) == 0, (
            f'BezierTri6 (bbar={bbar}, sigY={sig_y}) failed to converge at '
            f'step {step + 1}/{_NSTEPS}'
        )
    return _max_rel_asym()


@pytest.mark.parametrize('bbar', [False, True], ids=['std', 'bbar'])
def test_yielded_tangent_goes_unsymmetric(bbar):
    """The fingerprint pair: plastic K unsymmetric, elastic K symmetric.

    Pre-fix the mirror returned max|K - K'| == 0 identically, so the first
    assertion fails with certainty on a regressed binary; the second proves
    the rewritten full-product loop did not break the symmetric case and
    that the asymmetry really is the plastic non-associated signature.
    """
    asym_plastic = _run(_SIG_Y, bbar=bbar)
    asym_elastic = _run(_SIG_Y_ELASTIC, bbar=bbar)

    assert asym_elastic < 1.0e-9, (
        'never-yielding BezierTri6 run must keep a symmetric tangent',
        asym_elastic)
    assert asym_plastic > 1.0e-8, (
        'yielded BezierTri6 run must expose the unsymmetric consistent '
        'tangent (exactly 0 under the pre-fix mirror)', asym_plastic)
