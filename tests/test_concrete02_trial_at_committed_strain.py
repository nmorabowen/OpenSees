"""Concrete02: a trial strain back at the committed one must return the committed response.

The upstream early return on |deps| < DBL_EPSILON left sig and e at the previous TRIAL, so
the stress jumped by a whole iterate at that one strain. Fibers of a section whose strain
stays put between steps (constant axial load, near-zero moment: the column end away from
the hinge) land there, and the force-based element's internal iteration then cycles on the
jump forever, just above whatever tolerance it has. This was the "class E" stop of the blind
RC column recipe, on LadrunoForceBeamColumn and vanilla forceBeamColumn alike.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def _truss(e_committed):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Concrete02", 1, -35.0, -0.0025, -12.0, -0.018, 0.1, 2.5, 1500.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-12, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("DisplacementControl", 2, 1, e_committed)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return ops.nodeDisp(2, 1)


def _trial(u):
    ops.ladrunoSetNodeTrial(2, u, 0.0, 0.0)
    ops.ladrunoTrialResidualNorm()
    return ops.eleResponse(1, "material", "stress")[0], ops.eleResponse(1, "material", "tangent")[0]


@pytest.mark.parametrize("offset", [+1.0e-6, -1.0e-6])
def test_returning_to_the_committed_strain_returns_the_committed_stress(offset):
    u0 = _truss(-1.0e-4)
    s0, k0 = _trial(u0)
    s1, _ = _trial(u0 + offset)
    assert abs(s1 - s0) > 1.0e-3          # the excursion really moved the stress
    s2, k2 = _trial(u0)
    assert s2 == s0
    assert k2 == k0


def _column(element):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 1600.0)
    ops.fix(1, 1, 1, 1)
    ops.geomTransf("PDelta", 1)
    ops.uniaxialMaterial("Concrete02", 1, -40.0, -0.004, -20.0, -0.03, 0.1, 3.0, 1500.0)
    ops.uniaxialMaterial("Concrete02", 2, -35.0, -0.0025, -7.0, -0.008, 0.1, 3.0, 1500.0)
    ops.uniaxialMaterial("Steel02", 3, 438.0, 200000.0, 0.01, 20.0, 0.925, 0.15)
    ops.section("Fiber", 1)
    ops.patch("rect", 1, 12, 1, -130.0, -130.0, 130.0, 130.0)
    for y in (-130.0, 130.0):
        ops.patch("rect", 2, 2, 1, min(y, y + (20.0 if y > 0 else -20.0)), -150.0,
                  max(y, y + (20.0 if y > 0 else -20.0)), 150.0)
    ops.layer("straight", 3, 4, 500.0, -127.0, -127.0, -127.0, 127.0)
    ops.layer("straight", 3, 4, 500.0, 127.0, -127.0, 127.0, 127.0)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.element(element, 1, 1, 2, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -600.0e3, 0.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-7, 40, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    assert ops.analyze(10) == 0
    ops.loadConst("-time", 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(2, 1.0, 0.0, 0.0)


def _cycle(amps, du):
    u, done = 0.0, 0
    for a in amps:
        while abs(a - u) > 1.0e-9:
            d = max(-du, min(du, a - u))
            ops.integrator("DisplacementControl", 2, 1, d)
            rc = -1
            for alg in (("Newton",), ("KrylovNewton",)):
                ops.algorithm(*alg)
                rc = ops.analyze(1)
                if rc == 0:
                    break
            if rc != 0:
                return done, False
            u = ops.nodeDisp(2, 1)
            done += 1
    return done, True


@pytest.mark.parametrize("element", ["LadrunoForceBeamColumn", "forceBeamColumn"])
def test_cyclic_column_does_not_stall_on_the_kink(element):
    _column(element)
    amps = []
    for a in (4.0, 8.0, 12.0, 16.0, 24.0, 32.0):
        amps += [a, -a]
    n, ok = _cycle(amps, 0.25)
    assert ok, f"stopped after {n} increments"
