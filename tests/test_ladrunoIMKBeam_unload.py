"""WP-261: LadrunoIMKBeam (2D and 3D) hinges UNLOADING from a plastic state.

Regression for the local hinge Newton of ``ladrunoIMKSolveAxis``. With bilinear
Steel01 hinges (b = 1e-4, k0 = 100 * 6EI/L, so k0 ~ 150 * 4EI/L) the hinge
residual a(v - th) - M(th) is monotone but S-shaped. A hinge that must unload
from a plastic committed state sent the plain Newton into the arctan 2-cycle
across the narrow elastic band. The step then failed or cycled, and in frame
analyses a line-search recovery committed a reverse-yield iterate silently.

One beam, both ends hinged, node 1 fixed, node 2 rotation-controlled (so
v = (0, phi)): load into the plastic range, then reverse in coarse steps and
cycle. The controlled rotation is the only free DOF, so the global test cannot
hide anything: the committed element state is exactly what the local solve
returns. At every step:

  * the step converges with plain Newton;
  * the hinge rotations and moments equal an independent oracle: the exact
    series solution, found by enumerating the Steel01 regions of each hinge
    from its committed state (unique by monotonicity);
  * every reported moment lies on the Steel01 clip of its own rotation (the
    hinge law holds in the state about to be committed).

Before WP-261 the first unloading step (step 20) "converged" in one global
iteration onto th_j = 0.0177 with M_j = 1.4e3, while the Steel01 clip of that
rotation is ~1.0e5 (the exact answer is th_j = 1.19e-3, M_j = 8.95e4): the
silent reverse-yield commit seen in the 20-story frame, on one element.
"""
import itertools

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_A, _E, _Iz, _L = 1.0e-2, 2.0e11, 2.0e-5, 3.0
_a, _b = 4 * _E * _Iz / _L, 2 * _E * _Iz / _L
_k0 = 100.0 * 6.0 * _E * _Iz / _L
_My, _bh = 1.0e5, 1.0e-4
_kp, _cb = _bh * _k0, _My * (1.0 - _bh)
_EPS = np.finfo(float).eps


class _Clip:
    """Steel01 (default a1..a4) trial moment from a committed (theta_n, M_n)."""

    def __init__(self):
        self.thn, self.Mn = 0.0, 0.0

    def __call__(self, th):
        c = self.Mn + _k0 * (th - self.thn)
        return max(_kp * th - _cb, min(_kp * th + _cb, c))

    def region_line(self, r):
        # M = k * th + m on region r (-1, 0, +1)
        return (_k0, self.Mn - _k0 * self.thn) if r == 0 else (_kp, r * _cb)


def _exact(clips, v):
    K = np.array([[_a, _b], [_b, _a]])
    for regs in itertools.product((-1, 0, 1), repeat=2):
        kk, m0 = zip(*(c.region_line(r) for c, r in zip(clips, regs)))
        th = np.linalg.solve(K + np.diag(kk), K @ v - np.array(m0))
        q = K @ (v - th)
        M = np.array([c(t) for c, t in zip(clips, th)])
        if np.abs(q - M).max() <= 1e-9 * (np.abs(q).max() + np.abs(M).max() + 1.0):
            return th, M
    raise AssertionError("no consistent Steel01 region pair")


def _protocol():
    up = [1.0e-3] * 20                  # to 0.02 rad (~20 theta_y of the hinge)
    back = [-2.0e-3] * 20               # unload, reverse-yield, to -0.02
    again = [4.0e-3] * 5 + [-4.0e-3] * 3  # coarse re-reversals
    return up + back + again


def _build(ndm=2):
    """Return the controlled rotation DOF of node 2 (strong axis)."""
    ops.wipe()
    if ndm == 2:
        ops.model("basic", "-ndm", 2, "-ndf", 3)
        ops.node(1, 0.0, 0.0)
        ops.node(2, _L, 0.0)
        ops.fix(1, 1, 1, 1)
        ops.fix(2, 1, 1, 0)
        ops.geomTransf("Linear", 1)
        ops.uniaxialMaterial("Steel01", 10, _My, _k0, _bh)
        ops.element("LadrunoIMKBeam", 1, 1, 2, _A, _E, _Iz, 1, "-hinge", "both", "-matZ", 10)
        dof, ref = 3, (0.0, 0.0, 1.0)
    else:   # the 3D element: same strong-axis kernel, its own update()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, _L, 0.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.fix(2, 1, 1, 1, 1, 1, 0)
        ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
        ops.uniaxialMaterial("Steel01", 10, _My, _k0, _bh)
        ops.element("LadrunoIMKBeam", 1, 1, 2, _A, _E, 8.0e10, 1.0e-5, 1.0e-5, _Iz, 1,
                    "-hinge", "both", "-matZ", 10)
        dof, ref = 6, (0.0, 0.0, 0.0, 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, *ref)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    return dof


@pytest.mark.parametrize("ndm", [2, 3])
def test_unloading_from_plastic_converges_and_matches_oracle(ndm):
    dof = _build(ndm)
    clips = [_Clip(), _Clip()]
    max_iters, reversals = 0, 0
    prev = 0.0
    for k, dphi in enumerate(_protocol()):
        if dphi * prev < 0:
            reversals += 1
        prev = dphi
        ops.integrator("DisplacementControl", 2, dof, dphi)
        ops.analysis("Static")
        assert ops.analyze(1) == 0, f"step {k} (dphi {dphi:+.0e}) failed"
        max_iters = max(max_iters, ops.testIter())

        phi = ops.nodeDisp(2, dof)
        th = np.array(ops.eleResponse(1, "hingeRotation")[:2])
        M = np.array(ops.eleResponse(1, "hingeMoment")[:2])
        th_x, M_x = _exact(clips, np.array([0.0, phi]))
        thy = _My / _k0
        assert np.abs(th - th_x).max() <= 1e-9 * thy + 1e-9 * np.abs(th_x).max(), (
            f"step {k}: hinge rotations {th} vs exact {th_x}")
        assert np.abs(M - M_x).max() <= 1e-9 * _My, (
            f"step {k}: hinge moments {M} vs exact {M_x}")
        # the hinge law in the state about to be committed
        for c, t, m in zip(clips, th, M):
            assert abs(m - c(t)) <= 1e-9 * _My, f"step {k}: M {m} off the clip {c(t)}"
        for c, t, m in zip(clips, th, M):
            c.thn, c.Mn = t, m

    assert reversals >= 3
    assert max_iters <= 8, f"global Newton needed {max_iters} iterations"
    # the protocol really did reverse-yield both hinges
    assert min(c.Mn for c in clips) < 0.0


def test_unloading_step_is_the_reported_two_cycle_case():
    """The first unloading step from 0.02 rad by -0.002 is the case the
    g++ harness (test_ladrunoIMK_safe_solve_cpp.py) shows the plain local
    Newton cannot solve. Here it must land on the exact (elastic-unloading)
    state."""
    _build()
    ops.integrator("DisplacementControl", 2, 3, 1.0e-3)
    ops.analysis("Static")
    assert ops.analyze(20) == 0
    clips = [_Clip(), _Clip()]
    for c, t, m in zip(clips, ops.eleResponse(1, "hingeRotation")[:2],
                       ops.eleResponse(1, "hingeMoment")[:2]):
        c.thn, c.Mn = t, m
    ops.integrator("DisplacementControl", 2, 3, -2.0e-3)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    assert abs(ops.nodeDisp(2, 3) - 0.018) < 1e-12
    _, M_x = _exact(clips, np.array([0.0, 0.018]))
    M = np.array(ops.eleResponse(1, "hingeMoment")[:2])
    assert np.abs(M - M_x).max() <= 1e-9 * _My, f"hinge moments {M} vs exact {M_x}"


@pytest.mark.parametrize("ndm", [2, 3])
def test_tiny_first_step_in_steel01_dead_zone_does_not_fail(ndm):
    """Steel01 ignores a strain change below DBL_EPSILON (absolute), so near
    zero M(th) jumps by k0 * eps and, for eps < v < (1 + k0/a) eps, the hinge
    equation has no root. The first steps of a ground motion sit there. The
    residual is the law's resolution (~k0 * 2e-16), so the step must pass,
    as it did before WP-261, and not be rejected as an unconverged hinge."""
    dof = _build(ndm)
    for dphi in (1.0e-14, 1.0e-14, -3.0e-14, 1.0e-3):
        ops.integrator("DisplacementControl", 2, dof, dphi)
        ops.analysis("Static")
        assert ops.analyze(1) == 0, f"step dphi {dphi} rejected"
