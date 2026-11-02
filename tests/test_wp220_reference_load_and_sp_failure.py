"""WP-220: the path-control reference load, and element failures under imposed displacements.

phat was R(lambda+1), which assumes no unbalance now. A pattern added at t != 0, or an
integrator re-issued with a residual left over, put that residual in phat: at a drive DOF
r_new = r*R_prev/(q + R_prev), and with q = 1 N it fed on itself until an interior node of two
trusses in series carried 76 kN. phat is now R(lambda+1) - R(lambda), except on a Constant
series, which has no such difference and keeps R(lambda+1).

TransformationConstraintHandler::enforceSPs dropped the update code of every element. A
transient step now fails at newStep. A static step already failed in Domain::update.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_PATH_CONTROLS = {
    "DisplacementControl": ("DisplacementControl", 2, 1, 0.5),
    "MinUnbalDispNorm": ("MinUnbalDispNorm", 0.5),
    "ArcLength": ("ArcLength", 0.5, 1.0),
    "ArcLength1": ("ArcLength1", 0.5, 1.0),
    "HSConstraint": ("HSConstraint", 0.5),
    "LadrunoIndirectControl": ("LadrunoIndirectControl", 0.5, "-dof", 2, 1, 1.0),
    "LadrunoArcLength": ("LadrunoArcLength", 0.5, 1.0),
}


@pytest.mark.parametrize("name", sorted(_PATH_CONTROLS))
@pytest.mark.parametrize("t0", [1.0, 3.0])
def test_pattern_added_at_nonzero_time_converges_at_once(name, t0):
    # Elastic, one free DOF: with a true reference load the first iteration balances it. The
    # old phat left r*t/(1+t) per iteration: 36 and 87 iterations for DisplacementControl.
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1000.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Elastic", 1, 200.0e3)
    ops.element("Truss", 1, 1, 2, 1000.0, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormUnbalance", 1.0e-6, 200, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", t0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(2, 1000.0)
    ops.integrator(*_PATH_CONTROLS[name])
    assert ops.analyze(1) == 0
    assert ops.testIter() <= 3


def _series_trusses(q, tol):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1000.0)
    ops.node(3, 2000.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Steel02", 1, 400.0, 200.0e3, 0.01, 15.0, 0.925, 0.15)
    ops.element("Truss", 1, 1, 2, 1000.0, 1)
    ops.element("Truss", 2, 2, 3, 1100.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, q)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", tol, 40, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    u, worst = 0.0, 0.0
    for a in (2.0, -3.0, 4.0, -5.0, 6.0):
        while abs(a - u) > 1.0e-9:
            du = max(-0.1, min(0.1, a - u))
            ops.integrator("DisplacementControl", 3, 1, du)
            assert ops.analyze(1) == 0
            u = ops.nodeDisp(3, 1)
            worst = max(worst, abs(ops.eleForce(1)[1] + ops.eleForce(2)[0]))
    return worst


def test_reissued_integrator_leaves_no_interior_force():
    assert _series_trusses(1.0, 1.0e-4) < 1.0


def _sp_driven(element_args, step):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 1000.0)
    ops.fix(1, 1, 1, 1)
    ops.uniaxialMaterial("Steel02", 1, 400.0, 200.0e3, 0.01, 15.0, 0.925, 0.15)
    ops.section("Fiber", 1)
    ops.patch("rect", 1, 8, 1, -100.0, -100.0, 100.0, 100.0)
    ops.geomTransf("Linear", 1)
    ops.beamIntegration("Lobatto", 1, 1, 5)
    ops.element(*element_args)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(2, 1, 1.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-8, 10, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", step)
    ops.analysis("Static")


def test_constant_series_reference_load_still_drives():
    # The load does not follow lambda, so R(lambda+1) - R(lambda) is zero. Such decks keep working.
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1000.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Elastic", 1, 200.0e3)
    ops.element("Truss", 1, 1, 2, 1000.0, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1000.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-8, 10, 0)
    ops.algorithm("Newton")
    ops.integrator("DisplacementControl", 2, 1, 0.5)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


def test_static_element_failure_under_imposed_displacement_commits_nothing():
    _sp_driven(("forceBeamColumn", 1, 1, 2, 1, 1, "-iter", 1, 1.0e-30), 200.0)
    assert ops.analyze(1) == -3
    assert ops.nodeDisp(2, 1) == 0.0


def test_transient_element_failure_under_imposed_displacement_fails_at_newstep():
    _sp_driven(("forceBeamColumn", 1, 1, 2, 1, 1, "-iter", 1, 1.0e-30), 200.0)
    ops.wipeAnalysis()
    ops.remove("sp", 2, 1)
    ops.timeSeries("Linear", 5, "-factor", 40000.0)
    ops.pattern("Plain", 5, 5)
    ops.sp(2, 1, 1.0)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-8, 10, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    assert ops.analyze(1, 0.005) == -2
    assert ops.nodeDisp(2, 1) == 0.0


def test_imposed_displacement_that_updates_still_runs():
    _sp_driven(("forceBeamColumn", 1, 1, 2, 1, 1), 0.5)
    assert ops.analyze(1) == 0
    assert ops.nodeDisp(2, 1) == pytest.approx(0.5)
