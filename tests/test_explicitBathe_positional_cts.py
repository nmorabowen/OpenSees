"""ExplicitBathe: the documented <computeCriticalTimestep> positional (r6 parity).

`integrator ExplicitBathe <p> <computeCriticalTimestep>` is the documented grammar
(and the upstream one), but the parser read only p and passed the flag to the option
loop as "unknown option 1 (ignored)", so criticalTimeStep() stayed -1. Checked on a
one-element 2-D Truss, E = 100, A = 1, L = 1, rho = 1 (le/c = 0.1):

  * `0.54 1` gives criticalTimeStep = 0.1 (1e-9), equal to `0.54 -cfl`;
  * `0.54 0` and `0.54` give <= 0;
  * the positional composes with options (`0.54 1 -lump diagonal` on consistent mass:
    2/sqrt(600));
  * the deprecated ExplicitBatheLNVD alias honours its positional after <p> <alpha>.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def _truss(integrator_args, cmass=0):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 1.0, 0.0)
    ops.fix(2, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 100.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1, "-rho", 1.0, "-cMass", cmass)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1.0e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    assert ops.analyze(1, 1.0e-3) == 0
    return ops.criticalTimeStep()


def test_positional_flag_turns_the_estimate_on():
    dt = _truss(("ExplicitBathe", 0.54, 1))
    assert dt == pytest.approx(0.1, rel=1e-9)
    assert dt == _truss(("ExplicitBathe", 0.54, "-cfl"))


@pytest.mark.parametrize("args", [("ExplicitBathe", 0.54, 0), ("ExplicitBathe", 0.54)])
def test_no_estimate_without_the_flag(args):
    assert _truss(args) <= 0.0


def test_positional_flag_composes_with_options():
    dt = _truss(("ExplicitBathe", 0.54, 1, "-lump", "diagonal"), cmass=1)
    assert dt == pytest.approx(2.0 / math.sqrt(600.0), rel=1e-9)


def test_lnvd_alias_positional_flag():
    assert _truss(("ExplicitBatheLNVD", 0.54, 0.8, 1)) == pytest.approx(0.1, rel=1e-9)
