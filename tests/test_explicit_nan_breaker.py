"""NaN-capable circuit breaker on the Ladruno explicit integrators.

The pre-fix breaker computed `A_max = U.pNorm(0)` and tested `A_max != A_max`.
Vector::pNorm(0) implements max via `(data > value) ? data : value`, and every
comparison against NaN is FALSE — NaN entries are silently SKIPPED, so A_max is
never NaN and the breaker only ever fired on +Inf. An all-NaN acceleration
(inf - inf residual, a NaN material state, a poisoned initial condition) sailed
straight through and was committed into the node state.

Fixture: a NaN committed displacement makes Fint = k*NaN = NaN -> the residual and
the solved acceleration are NaN on the very first step. Pre-fix, analyze() returns
0 and the NaN silently propagates into the committed node state; post-fix the
breaker aborts with the existing "non-finite acceleration" message.

Sites covered: CentralDifferenceLadruno, unified ExplicitBathe, and
LadrunoDynamicRelaxation (same pNorm(0) pattern in all three).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

NAN = float("nan")


def _sdof(integrator_args):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0)
    ops.fix(2, 0, 1)
    ops.mass(2, 1.0, 0.0)
    ops.uniaxialMaterial("Elastic", 1, 100.0)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.setNodeDisp(2, 1, NAN, "-commit")           # NaN state -> NaN residual/accel
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")


@pytest.mark.parametrize(
    "integrator_args",
    [
        ("CentralDifferenceLadruno",),
        ("ExplicitBathe", 0.54),
        ("LadrunoDynamicRelaxation",),
    ],
    ids=["CDL", "ExplicitBathe", "DynamicRelaxation"],
)
def test_nan_acceleration_trips_breaker(integrator_args, capfd):
    _sdof(integrator_args)
    capfd.readouterr()                              # drop setup chatter
    rc = ops.analyze(1, 1.0e-3)
    err = capfd.readouterr().err

    assert rc != 0, (
        "%s: analyze() must FAIL on an all-NaN acceleration (got rc=0 — the pNorm(0) "
        "breaker is NaN-blind and the NaN was silently committed)" % (integrator_args[0],)
    )
    assert "non-finite" in err, (
        "%s: the non-finite-acceleration ABORT message must fire; stderr was:\n%s"
        % (integrator_args[0], err)
    )
    # the abort must have prevented the NaN from being committed: after the
    # framework's revertToLastCommit the committed displacement is the pre-step
    # value (NaN was the SEED here, so just require the run STOPPED — the real
    # guard is rc != 0 above; this leg pins that no extra step ran).
    assert ops.getTime() == pytest.approx(0.0, abs=1e-15), (
        "time advanced past the aborted step (%r)" % ops.getTime()
    )


def test_inf_acceleration_still_trips(capfd):
    """Regression lock: the old Inf path must keep working after the isfinite swap.

    A SEEDED Inf degenerates to NaN before reaching the breaker (inf + dt*(-inf)),
    so the honest Inf producer is natural divergence: run just above dt_cr and let
    the amplitude overflow — the first overflow step hands the breaker a genuine
    +/-Inf acceleration (this is exactly how the pre-fix breaker ever fired).
    Passes pre-fix AND post-fix (regression lock, not a change-detector)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0)
    ops.fix(2, 0, 1)
    ops.mass(2, 1.0, 0.0)
    ops.uniaxialMaterial("Elastic", 1, 100.0)     # omega=10 -> dt_cr = 0.2
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.setNodeDisp(2, 1, 0.1, "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno")
    ops.analysis("Transient")
    capfd.readouterr()
    rc = 0
    for _ in range(3000):                          # |rho| ~ 1.9/step at dt=0.21
        rc = ops.analyze(1, 0.21)
        if rc != 0:
            break
    err = capfd.readouterr().err
    assert rc != 0, "divergent run never aborted — the Inf breaker path is dead"
    assert "non-finite" in err, (
        "the overflow abort must be the non-finite-acceleration breaker; stderr:\n%s" % err
    )
