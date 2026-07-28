"""`sp ... -subtractInit` — the openseespy path was a NO-OP. CI gate for the fix.

`OPS_SP()` in `SRC/interpreter/OpenSeesPatternCommands.cpp` handled
`-subtractInit` by setting `retZeroInitValue = true` — which is *already* the
initialiser a few lines above, so the branch did nothing and incremental SP was
unreachable from Python. Both Tcl parsers set it to `false` correctly
(`SRC/modelbuilder/tcl/TclModelBuilder.cpp:3808-3810`,
`SRC/runtime/commands/modeling/constraint.cpp:392-394`), so this was a
Python-vs-Tcl behaviour split, not a missing feature.

Failure mode it guards: a staged deck (gravity, then an imposed displacement)
silently YANKS the node to the absolute value instead of moving it incrementally
from its staged state — a wrong answer with no warning, in the direction of
"looks plausible". Measured on the fix commit: 3.0 (broken) vs 5.5 (correct).

HANDLER CAVEAT (why this test uses Transformation): only handlers that read
`SP_Constraint::getInitialValue()` honour the flag — `TransformationDOF_Group`
(:1068), `PenaltySP_FE` (:125), `LagrangeSP_FE` (:143). The base `DOF_Group`
never calls it, so under `constraints Plain` the flag is inert BY DESIGN, fixed
or not. Testing this under Plain would show no difference and read as a
regression; it is not one.

No MKL, no fork elements — pure upstream Truss + Elastic, so Zone-A gates it.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, A, L = 200000.0, 100.0, 1000.0
PX = 5.0e4          # stage-1 load  => u0 = PX*L/(E*A) = 2.5
DELTA = 3.0         # stage-2 imposed displacement
U0 = PX * L / (E * A)


def _run(flag, handler="Transformation"):
    """Stage 1 loads the node to U0; stage 2 imposes DELTA on the moved node."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, L, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)                    # leave dof 1 free for the SP
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.element("Truss", 1, 1, 2, A, 1)

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, PX, 0.0)
    ops.constraints(handler)
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-10, 20)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    u0 = ops.nodeDisp(2, 1)

    ops.loadConst("-time", 0.0)
    if flag:
        ops.sp(2, 1, DELTA, "-subtractInit")
    else:
        ops.sp(2, 1, DELTA)
    assert ops.analyze(1) == 0
    return u0, ops.nodeDisp(2, 1)


def test_stage1_reaches_expected_state():
    """Guards the premise: if u0 were 0 the whole test would pass vacuously."""
    u0, _ = _run(False)
    assert u0 == pytest.approx(U0, rel=1e-9)
    assert abs(u0) > 1e-6


def test_without_flag_is_absolute():
    """No flag => the node is driven to DELTA itself, discarding the staged state."""
    _, u = _run(False)
    assert u == pytest.approx(DELTA, rel=1e-9)


def test_subtract_init_is_incremental():
    """The fix: DELTA is applied ON TOP of the staged displacement.

    `TransformationDOF_Group::enforceSPs` does `setTrialDisp(value + initial)`,
    so the expected answer is DELTA + u0 = 5.5, not DELTA = 3.0.
    """
    u0, u = _run(True)
    assert u == pytest.approx(DELTA + u0, rel=1e-9)


def test_flag_actually_changes_the_answer():
    """The regression that would catch a re-introduction of the no-op.

    Stated as a DIFFERENCE rather than an absolute so it stays meaningful even
    if the sign convention of the initial-value subtraction is ever revisited:
    whatever `-subtractInit` means, it must not mean 'nothing'.
    """
    _, u_plain = _run(False)
    _, u_flag = _run(True)
    assert abs(u_flag - u_plain) > 1e-9, (
        "-subtractInit changed nothing — the OPS_SP() branch is a no-op again")


def test_inert_under_plain_handler_is_expected():
    """Documents the caveat so it is never mistaken for this bug returning.

    `constraints Plain` routes SPs through the base DOF_Group, which never calls
    getInitialValue(), so the flag has no effect there — before OR after the fix.
    """
    _, u_plain = _run(False, handler="Plain")
    _, u_flag = _run(True, handler="Plain")
    assert u_flag == pytest.approx(u_plain, rel=1e-9)
