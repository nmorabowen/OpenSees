"""Transient integrators reject a zero time step.

Model: one-DOF oscillator, truss of stiffness k = E*A/L = 100 between a
fixed node and a free node of mass m = 1, under a constant load P = 1.

TRBDF2, TRBDF3, Houbolt and BackwardEuler form their integration constants
as 2/dt and 4/dt^2. With dt = 0.0 the constants are inf, the predicted
acceleration -4/dt * 0 is NaN, and with algorithm Linear analyze(1, 0.0)
returns 0 and commits a NaN state. Each integrator must instead fail the
step (negative return code), leave the committed state finite (zero), and
accept a following step with dt > 0.
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


K = 100.0
M = 1.0
P = 1.0


def _build_sdof(integrator, algorithm):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.mass(2, M)
    ops.uniaxialMaterial("Elastic", 1, K)
    ops.element("truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, P)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-12, 20)
    ops.algorithm(*algorithm)
    ops.integrator(*integrator)
    ops.analysis("Transient")


def _state_is_finite(tag):
    values = (ops.nodeDisp(tag, 1), ops.nodeVel(tag, 1), ops.nodeAccel(tag, 1))
    return all(math.isfinite(v) for v in values)


@pytest.mark.parametrize("integrator", [
    ("TRBDF2",),
    ("TRBDF3",),
    ("Houbolt",),
    ("BackwardEuler",),
])
def test_zero_time_step_is_rejected(integrator):
    try:
        _build_sdof(integrator, ("Linear",))
        assert ops.analyze(1, 0.0) < 0
        assert _state_is_finite(2)
        assert ops.nodeDisp(2, 1) == 0.0
        assert ops.analyze(1, 0.01) == 0
        assert _state_is_finite(2)
        assert ops.nodeDisp(2, 1) > 0.0
    finally:
        ops.wipe()
