"""Non-finite acceleration breaker of the explicit central difference.

Vector::pNorm(0) implements the maximum with `(data > value) ? data : value`.
Every comparison against NaN is false, so NaN entries are skipped and a guard on
pNorm(0) only fires on +/-Inf. An all-NaN acceleration must stop the step.

Fixture for NaN: a NaN committed displacement makes F_int = k NaN, so the residual
and the acceleration are NaN on the first step.
"""
import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

NAN = float("nan")


def _sdof(d0):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0)
    ops.fix(2, 0, 1)
    ops.mass(2, 1.0, 0.0)
    ops.uniaxialMaterial("Elastic", 1, 100.0)  # omega = 10, dt_cr = 0.2
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.setNodeDisp(2, 1, d0, "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceExplicit")
    ops.analysis("Transient")


def test_nan_acceleration_trips_breaker(capfd):
    _sdof(NAN)
    capfd.readouterr()
    rc = ops.analyze(1, 1.0e-3)
    err = capfd.readouterr().err
    assert rc != 0, "analyze must fail on an all-NaN acceleration"
    assert "non-finite" in err, "breaker message missing; stderr:\n%s" % err
    assert ops.getTime() == pytest.approx(0.0, abs=1e-15), (
        "time advanced past the aborted step (%r)" % ops.getTime()
    )


def test_inf_acceleration_trips_breaker(capfd):
    """Divergence above dt_cr overflows into a genuine +/-Inf acceleration.

    At omega dt = 2.1 the amplification factor is the larger root of
    z^2 - (2 - 2.1^2) z + 1 = 0, |z| = 1.9 per step (rounded), so the 3000
    steps exceed the range of a double well before the end of the loop.
    """
    _sdof(0.1)
    capfd.readouterr()
    rc = 0
    for _ in range(3000):
        rc = ops.analyze(1, 0.21)
        if rc != 0:
            break
    err = capfd.readouterr().err
    assert rc != 0, "divergent run never aborted"
    assert "non-finite" in err, "breaker message missing; stderr:\n%s" % err
