"""WP-240 -- `ladrunoDR relax` and `ladrunoDR maxDisplacement`.

Ported from jaabell/OpenSees 9a429237e (José Abell).

    ladrunoDR relax -tolResidual r <-tolKE k> -maxSteps n <-runaway f>
                    <-uRef u> <-checkEvery m> <-dt dt>

runs one-step transient analyses until rest, runaway or n steps, and returns
[status, steps, residualNorm, kineticEnergy, maxDisplacement] with status
0 rest, 1 maxSteps, 2 runaway, 3 failed step.

  X-1  rest: status 0, and the march is the one a script loop of `analyze(1)`
       with the same test would do, step for step and bit for bit.
  X-2  -tolKE delays the stop until the kinetic energy is below it too.
  X-3  -maxSteps: status 1 after exactly n steps of -dt.
  X-4  -runaway with an explicit -uRef: status 2, decided on a -checkEvery
       multiple; with the default uRef (max|u| after the first step) a
       relaxation to a far equilibrium is flagged, and a large factor is not.
  X-5  a failed step: status 3, and the step is counted.
  X-6  maxDisplacement is max|u| over the free DOFs.
  X-7  bad input is an error, not a silent default.

The classic-Tcl bridge is covered by tests/tcl/ladruno_solver_queries.tcl.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


def _arch(load):
    """Von Mises truss (snap-through arch), apex load `load` downwards. The
    limit load is about 3.80; 3.0 relaxes on the near branch, 4.4 snaps to the
    far one (uy about -0.22)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, -1.0, 0.0)
    ops.node(2, 0.0, 0.1)
    ops.node(3, 1.0, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(3, 1, 1)
    ops.fix(2, 1, 0)
    ops.uniaxialMaterial("Elastic", 1, 1.0e4)
    ops.element("corotTruss", 1, 1, 2, 1.0, 1)
    ops.element("corotTruss", 2, 2, 3, 1.0, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.0, -load)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormUnbalance", 1.0e30, 1, 0)
    ops.algorithm("Linear")
    ops.integrator("LadrunoDynamicRelaxation")
    ops.analysis("Transient")


def _relax(*args):
    out = ops.ladrunoDR("relax", *args)
    assert len(out) == 5, out
    return int(out[0]), int(out[1]), out[2], out[3], out[4]


def test_X1_rest_equals_the_script_loop():
    tol, cap = 1.0e-8, 50000
    _arch(3.0)
    status, steps, res, ke, umax = _relax("-tolResidual", tol, "-maxSteps", cap)
    assert status == 0
    assert 1 < steps < cap
    assert 0.0 <= res < tol
    uy = ops.nodeDisp(2, 2)
    assert umax == abs(uy)
    assert res == ops.ladrunoDR("residualNorm")
    assert ke == ops.ladrunoDR("kineticEnergy")
    assert ops.getTime() == pytest.approx(steps * 1.0)

    _arch(3.0)
    n = 0
    while n < cap:
        assert ops.analyze(1, 1.0) == 0
        n += 1
        if ops.ladrunoDR("residualNorm") < tol:
            break
    assert n == steps
    assert ops.nodeDisp(2, 2) == uy


def test_X2_tolKE_holds_the_stop_until_the_kinetic_energy_is_small():
    tol = 1.0e-3
    _arch(3.0)
    s0, n0, _, ke0, _ = _relax("-tolResidual", tol, "-maxSteps", 50000)
    assert s0 == 0 and ke0 > 0.0
    tol_ke = 1.0e-6 * ke0
    _arch(3.0)
    s1, n1, res1, ke1, _ = _relax("-tolResidual", tol, "-tolKE", tol_ke,
                                  "-maxSteps", 50000)
    assert s1 == 0
    assert n1 > n0
    assert res1 < tol and ke1 < tol_ke


def test_X3_maxSteps_is_status_1_after_exactly_n_steps():
    _arch(3.0)
    status, steps, res, _, _ = _relax("-tolResidual", 1.0e-30, "-maxSteps", 7,
                                      "-dt", 1.0)
    assert (status, steps) == (1, 7)
    assert res > 0.0
    assert ops.getTime() == pytest.approx(7.0)
    # a second call goes on from where the first stopped
    status, steps, _, _, _ = _relax("-tolResidual", 1.0e-30, "-maxSteps", 5)
    assert (status, steps) == (1, 5)
    assert ops.getTime() == pytest.approx(12.0)


def test_X4_runaway_is_status_2_on_a_checkEvery_multiple():
    # super-critical load: the apex snaps to |uy| about 0.22
    _arch(4.4)
    status, steps, _, _, umax = _relax("-tolResidual", 1.0e-10, "-maxSteps", 50000,
                                       "-runaway", 5.0, "-uRef", 0.02,
                                       "-checkEvery", 10)
    assert status == 2
    assert steps % 10 == 0 and steps > 0
    assert umax > 5.0 * 0.02
    assert umax == abs(ops.nodeDisp(2, 2))

    # the same run with the threshold above the far branch comes to rest
    _arch(4.4)
    status, _, res, _, umax = _relax("-tolResidual", 1.0e-10, "-maxSteps", 200000,
                                     "-runaway", 5.0, "-uRef", 0.1,
                                     "-checkEvery", 10)
    assert status == 0 and res < 1.0e-10
    assert 0.15 < umax < 0.5


def test_X4_default_uRef_is_the_displacement_after_the_first_step():
    _arch(3.0)
    assert ops.analyze(1, 1.0) == 0
    u1 = abs(ops.nodeDisp(2, 2))
    assert u1 > 0.0

    # the equilibrium is far beyond 3 x the first step: flagged
    _arch(3.0)
    status, steps, _, _, umax = _relax("-tolResidual", 1.0e-10, "-maxSteps", 50000,
                                       "-runaway", 3.0, "-checkEvery", 1)
    assert status == 2
    assert umax > 3.0 * u1
    # ... and not one step later than it had to be
    _arch(3.0)
    assert ops.analyze(steps - 1, 1.0) == 0
    assert abs(ops.nodeDisp(2, 2)) <= 3.0 * u1

    # a large factor does not fire on a relaxation that converges
    _arch(3.0)
    status, _, _, _, _ = _relax("-tolResidual", 1.0e-8, "-maxSteps", 50000,
                                "-runaway", 1.0e6)
    assert status == 0


def test_X5_a_failed_step_is_status_3():
    _arch(3.0)
    status, steps, _, _, _ = _relax("-tolResidual", 1.0e-30, "-maxSteps", 20)
    assert (status, steps) == (1, 20)
    t = ops.getTime()
    # DR refuses a second update() in one step (WP-153's injection)
    ops.test("NormUnbalance", 1.0e-300, 2, 0)
    ops.algorithm("Newton")
    status, steps, _, _, _ = _relax("-tolResidual", 1.0e-30, "-maxSteps", 20)
    assert (status, steps) == (3, 1)
    assert ops.getTime() == pytest.approx(t)


def test_X6_maxDisplacement_is_the_largest_free_dof_displacement():
    _arch(3.0)
    assert ops.ladrunoDR("maxDisplacement") == 0.0
    assert ops.analyze(25, 1.0) == 0
    assert ops.ladrunoDR("maxDisplacement") == abs(ops.nodeDisp(2, 2)) > 0.0


@pytest.mark.parametrize("args", [
    (),                                             # no -tolResidual, no -maxSteps
    ("-maxSteps", 10),                              # no -tolResidual
    ("-tolResidual", 1.0e-6),                       # no -maxSteps
    ("-tolResidual", -1.0, "-maxSteps", 10),        # not positive
    ("-tolResidual", 1.0e-6, "-maxSteps", 0),
    ("-tolResidual", 1.0e-6, "-maxSteps", 10, "-notAnOption", 1.0),
    ("-tolResidual", 1.0e-6, "-maxSteps"),          # value missing
])
def test_X7_bad_input_is_an_error(args):
    _arch(3.0)
    with pytest.raises(Exception):
        ops.ladrunoDR("relax", *args)
    assert ops.getTime() == 0.0                     # nothing marched


def test_X7_relax_needs_the_transient_analysis():
    _arch(3.0)
    ops.wipeAnalysis()
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.algorithm("Linear")
    ops.integrator("LadrunoDynamicRelaxation")
    with pytest.raises(Exception):
        ops.ladrunoDR("relax", "-tolResidual", 1.0e-6, "-maxSteps", 10)
