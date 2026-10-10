"""WP-139: LadrunoBrick `-lumped` uses ONE mass model for the inertia residual and the tangent.

Before WP-139 the residual inertia was the CONSISTENT mass (inherited from upstream Brick)
while getMass() -- the Newton tangent, alphaM Rayleigh and the ground load -- was the row-sum
LUMPED mass: implicit dynamics integrated a hybrid and Newton converged only linearly (an
elastic step did not reach 1e-10 in 40 iterations; WP-124 gap C8). Oracles, on a regular unit
cube where row-sum lumping gives every node exactly rho V / 8:

  * the inertia residual is M_L a, node by node, for an ARBITRARY acceleration field (the
    consistent mass couples the nodes, so it cannot pass);
  * Newton converges quadratically again (<= 3 iterations per Newmark step);
  * without -lumped nothing changed: the consistent residual still couples the nodes and its
    x-sum is rho V times the mean acceleration.
"""
import pytest

import opensees as ops

pytestmark = [pytest.mark.zone_a]

E, NU, RHO = 1000.0, 0.3, 2.0
CUBE = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0)]
ACCEL = [(0.7, -0.2, 0.1), (0.0, 0.3, 0.0), (-0.4, 0.0, 0.9), (0.0, 0.0, 0.0),
         (1.1, 0.5, -0.3), (0.0, -0.6, 0.0), (0.2, 0.0, 0.4), (-0.8, 0.1, 0.0)]


def _cube(lumped, fixed=()):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, xyz in enumerate(CUBE):
        ops.node(i + 1, *xyz)
    for n in fixed:
        ops.fix(n, 1, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
    ops.element("LadrunoBrick", 1, *range(1, 9), 1, *(("-lumped",) if lumped else ()))


def _inertial(lumped):
    _cube(lumped)
    for n, acc in enumerate(ACCEL, start=1):
        for d, a in enumerate(acc, start=1):
            ops.setNodeAccel(n, d, a, "-commit")
    return ops.eleResponse(1, "inertialForce")


def test_lumped_inertia_residual_is_M_L_a():
    f = _inertial(lumped=True)
    m = RHO * 1.0 / 8.0
    want = [m * a for acc in ACCEL for a in acc]
    assert len(f) == 24
    for c, (got, w) in enumerate(zip(f, want)):
        assert abs(got - w) <= 1.0e-12, f"dof {c}: inertia {got!r}, lumped M_L a = {w!r}"


def test_consistent_mass_is_unchanged():
    """Without -lumped: sum_x (M_c a) = sum_j (column sum j) a_jx = (rho V / 8) sum_j a_jx on the
    regular cube -- the same TOTAL as lumped -- but node by node the consistent mass couples."""
    f = _inertial(lumped=False)
    m = RHO * 1.0 / 8.0
    fx = sum(f[i] for i in range(0, 24, 3))
    want = m * sum(acc[0] for acc in ACCEL)
    assert abs(fx - want) <= 1.0e-12 * max(1.0, abs(want))
    coupled = any(abs(f[3 * n + d] - m * ACCEL[n][d]) > 1.0e-6 for n in range(8) for d in range(3))
    assert coupled, "the consistent mass must still couple the nodes"


@pytest.mark.parametrize("alphaM", [0.0, 0.5])
def test_lumped_newton_converges_quadratically(alphaM):
    _cube(lumped=True, fixed=(1, 2, 3, 4))
    ops.rayleigh(alphaM, 0.0, 0.0, 0.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (5, 6, 7, 8):
        ops.load(n, 1.0, 0.3, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 4, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient", "-noWarnings")
    for step in range(10):
        assert ops.analyze(1, 0.01) == 0, f"step {step}: Newton did not converge in 4 iterations"
        assert ops.testIter() <= 3, f"step {step}: {ops.testIter()} iterations (linear convergence)"
