"""`constraints Auto` + a non-homogeneous `sp` was a SILENT WRONG ANSWER. CI gate.

`AutoConstraintHandler::applyLoad()` ran the two `enforceSPs` loops but omitted
the `updateElement()` loop that `TransformationConstraintHandler::enforceSPs`
has (`TransformationConstraintHandler.cpp:518-521`). Both handlers route
SP-constrained nodes through a `TransformationDOF_Group`, so the constrained DOF
is genuinely eliminated (`setID(dof,-1)`) and there is no `K*du_prescribed`
column: the ONLY way the prescribed motion reaches the RHS is via the state of
the elements attached to the driven node. Skip the element update and the first
residual of the step is stale.

Failure mode it guards (measured on the pre-fix binary, ADR-80 gate G4):

    two equal trusses in series, node 3 driven to 3.0, node 2 free
    Transformation : u2 = 1.5   (correct)     2 iterations
    Auto           : u2 = 0.0   (WRONG)       1 iteration, analyze() == 0

`analyze()` returns success. Nothing warns. The step commits with only the driven
node moved, because a stale residual makes dU ~ 0 and a dU-based convergence test
declares victory at iteration 1 -- `CTestNormDispIncr` has no minimum-iteration
guard (`CTestNormDispIncr.cpp:160-175`).

TEST-DEPENDENCE IS THE MECHANISM, and it is why this stayed hidden: the bug is
only fatal for convergence tests that look at dU. `NormDispIncr` and `EnergyIncr`
both commit the wrong answer; `NormUnbalance` looks at the (large, stale)
residual instead, refuses to converge at iteration 1, and the next
`LoadControl::update` -> `updateDomain()` silently repairs the state -- so the
same model gives the right answer purely by choice of test. Both branches are
asserted below; a fix that repaired only the dU tests would be incomplete.

Contrast worth keeping: `constraints Plain` is ALSO wrong here (it cannot do
non-homogeneous SPs) but it is LOUD -- it prints "non-homogeneos constraint for
node 3 homogeneous constraint assumed". That is documented by-design behaviour,
not this bug. Auto was wrong and silent.

No MKL, no fork elements -- upstream Truss/stdBrick + Elastic, so Zone-A gates it.
Adjacent precedent, same code path: `tests/test_sp_subtract_init.py`.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, A, L = 200000.0, 100.0, 1000.0
DELTA = 3.0
EXPECTED = 0.5 * DELTA          # two equal springs in series => midpoint halves

# Tests whose convergence measure is dU-based, i.e. the ones the stale residual
# fools. NormUnbalance is deliberately NOT in this list -- see module docstring.
DU_BASED_TESTS = ["NormDispIncr", "EnergyIncr"]
ALL_TESTS = DU_BASED_TESTS + ["NormUnbalance"]


def _truss_chain(handler, testtype):
    """Node 1 fixed, node 3 driven to DELTA, node 2 free. Returns (ok, u2, iters)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.node(3, 2.0 * L, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.fix(3, 0, 1)                 # dof 1 left free for the sp
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.element("Truss", 1, 1, 2, A, 1)
    ops.element("Truss", 2, 2, 3, A, 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(3, 1, DELTA)

    ops.constraints(handler)
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test(testtype, 1e-10, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ok = ops.analyze(1)
    return ok, ops.nodeDisp(2, 1), ops.testIter()


def _brick_stack(handler, testtype):
    """Two stacked stdBricks, ux/uy fixed => a 1D chain. Returns (ok, uz_mid, iters)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    n = 0
    for z in (0.0, 1.0, 2.0):
        for x, y in ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)):
            n += 1
            ops.node(n, x, y, z)
    ops.nDMaterial("ElasticIsotropic", 1, 200000.0, 0.25)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.element("stdBrick", 2, 5, 6, 7, 8, 9, 10, 11, 12, 1)
    for i in range(1, 13):
        ops.fix(i, 1, 1, 1) if i <= 4 else ops.fix(i, 1, 1, 0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for i in range(9, 13):
        ops.sp(i, 3, DELTA)

    ops.constraints(handler)
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test(testtype, 1e-10, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ok = ops.analyze(1)
    return ok, ops.nodeDisp(5, 3), ops.testIter()


MODELS = {"truss": _truss_chain, "brick": _brick_stack}


@pytest.mark.parametrize("model", sorted(MODELS))
def test_transformation_is_the_reference(model):
    """Guards the premise: the reference handler must give the analytic answer.

    If this ever fails the rest of the file is meaningless -- it would mean the
    baseline moved, not that Auto regressed.
    """
    ok, u, _ = MODELS[model]("Transformation", "NormDispIncr")
    assert ok == 0
    assert u == pytest.approx(EXPECTED, rel=1e-9)


@pytest.mark.parametrize("model", sorted(MODELS))
@pytest.mark.parametrize("testtype", ALL_TESTS)
def test_auto_matches_transformation(model, testtype):
    """THE REGRESSION. Auto must agree with Transformation for every test type.

    Pre-fix this failed for NormDispIncr and EnergyIncr with u = 0.0 exactly.
    """
    ok_a, u_auto, _ = MODELS[model]("Auto", testtype)
    ok_t, u_ref, _ = MODELS[model]("Transformation", testtype)
    assert ok_a == 0 and ok_t == 0
    assert u_auto == pytest.approx(u_ref, rel=1e-9), (
        f"Auto disagrees with Transformation under {testtype}: "
        f"{u_auto} vs {u_ref} -- the updateElement() loop in "
        f"AutoConstraintHandler::applyLoad() is missing again")


@pytest.mark.parametrize("model", sorted(MODELS))
@pytest.mark.parametrize("testtype", DU_BASED_TESTS)
def test_auto_does_not_leave_the_interior_unmoved(model, testtype):
    """The bug's exact signature, stated positively so it cannot pass vacuously.

    The pre-fix answer was not merely inaccurate -- the interior DOF never moved
    at all. Asserting `> 0` catches the failure even if EXPECTED were ever
    retuned for a different mesh.
    """
    _, u, _ = MODELS[model]("Auto", testtype)
    assert abs(u) > 1e-6, (
        "interior node did not move: the step committed with only the driven "
        "layer displaced (stale first residual, false convergence)")
    assert u == pytest.approx(EXPECTED, rel=1e-9)


@pytest.mark.parametrize("model", sorted(MODELS))
def test_auto_no_longer_converges_at_iteration_one(model):
    """The tell that made this diagnosable, kept as a gate.

    A single Newton iteration under a dU-based test means the first residual
    carried no information about the prescribed motion. The corrected path needs
    at least two.
    """
    _, _, iters = MODELS[model]("Auto", "NormDispIncr")
    assert iters >= 2, (
        "converged in one iteration -- the first residual was stale again")


def test_homogeneous_only_model_is_unaffected():
    """Blast-radius guard: the fix must not perturb ordinary fixed-only models.

    The added loop runs for every element touching ANY SP-constrained node,
    including plain `fix`. For homogeneous SPs enforceSPs writes an unchanged
    value, so updateElement() is numerically idempotent -- this pins that.
    """
    def _loaded_chain(handler):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.node(1, 0.0, 0.0)
        ops.node(2, L, 0.0)
        ops.node(3, 2.0 * L, 0.0)
        ops.fix(1, 1, 1)
        ops.fix(2, 0, 1)
        ops.fix(3, 0, 1)
        ops.uniaxialMaterial("Elastic", 1, E)
        ops.element("Truss", 1, 1, 2, A, 1)
        ops.element("Truss", 2, 2, 3, A, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(3, 1.0e4, 0.0)          # Neumann only -- no non-homogeneous sp
        ops.constraints(handler)
        ops.numberer("Plain")
        ops.system("BandGeneral")
        ops.test("NormDispIncr", 1e-10, 30, 0)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        return ops.nodeDisp(3, 1)

    assert _loaded_chain("Auto") == pytest.approx(_loaded_chain("Transformation"),
                                                  rel=1e-12)
