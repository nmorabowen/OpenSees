"""ADR-80 P3 — `LadrunoLoadControl -tangentPredictor` (the Kratos `b -= K*du_D` route).

WHY IT EXISTS. Under `constraints Transformation` a non-homogeneous `sp` has its
DOF genuinely eliminated (`setID(dof,-1)`), so the prescribed motion can only
reach the RHS by pre-updating the elements attached to the driven node and
reading their internal force back
(`TransformationConstraintHandler.cpp:518-521`). That forces the FIRST
constitutive evaluation of every increment to happen with the driven face
advanced by the whole increment and the interior still at the last commit -- an
overstrain of order L/h in one element layer. For a path-dependent law that
layer yields SPURIOUSLY and its consistent tangent collapses.

S1's `-extrapolate` reduced the overstrain and failed its acceptance gate
(cutbacks 23 -> 23; see `Ladruno_implementation/80c_s1_extrapolate_verdict_*`).
`-tangentPredictor` ELIMINATES it instead: `newStep` applies the domain loads
(which sets every sp value) but does NOT let the handler enforce them, so
iteration 1 forms K at the committed state and takes its prescribed-motion
forcing from `-K*du_D` assembled directly -- with no constitutive evaluation in
the lagging state anywhere. The sps are enforced in `update()`, on an interior
the predictor solve has already distributed.

WHAT IS GATED HERE
  1. the answer does not move (a predictor cannot change a converged solution);
  2. it is off by default -- an unflagged LadrunoLoadControl is stock;
  3. it is NOT a silent-wrong-answer machine. This is the sharp one: the route
     deliberately withholds the prescribed motion from iteration 1's residual.
     If the `K*du_D` term ever fails to assemble, dU comes back ~0 and a dU-based
     convergence test declares victory on a step that never moved -- which is
     EXACTLY the AutoConstraintHandler defect ADR-80 S2 fixed
     (`tests/test_auto_handler_sp_update.py`). The fallback in
     `LadrunoLoadControl::formUnbalance()` exists for that, and is gated below on
     the model that triggers it (no `sp` at all).
  4. `-extrapolate` and `-tangentPredictor` refuse to compose.

Upstream elements + Elastic/J2 only, no MKL -- Zone-A gates it.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, A, L = 200000.0, 100.0, 1000.0
DELTA = 3.0
EXPECTED = 0.5 * DELTA          # two equal springs in series => midpoint halves


def _truss_chain(integrator_args, testtype="NormDispIncr", drive="disp"):
    """Node 1 fixed, node 3 driven (or loaded), node 2 free.

    Returns (ok, u2, iters).
    """
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
    if drive == "disp":
        ops.sp(3, 1, DELTA)
    else:
        # the equivalent traction: same converged u2 = DELTA/2
        ops.load(3, E * A * DELTA / (2.0 * L), 0.0)

    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("BandGeneral")
    ops.test(testtype, 1e-10, 30, 0)
    ops.algorithm("Newton")
    ops.integrator(*integrator_args)
    ops.analysis("Static")
    ok = ops.analyze(1)
    return ok, ops.nodeDisp(2, 1), ops.testIter()


STOCK = ("LoadControl", 1.0)
LADRUNO = ("LadrunoLoadControl", 1.0)
TANGENT = ("LadrunoLoadControl", 1.0, "-tangentPredictor")


def test_premise_stock_gives_the_analytic_answer():
    """Guards the baseline; without it every comparison below is vacuous."""
    ok, u, _ = _truss_chain(STOCK)
    assert ok == 0
    assert u == pytest.approx(EXPECTED, rel=1e-9)


@pytest.mark.parametrize("testtype", ["NormDispIncr", "EnergyIncr", "NormUnbalance"])
def test_tangent_predictor_does_not_move_the_answer(testtype):
    """A predictor changes the PATH, never the converged solution."""
    ok_t, u_tan, _ = _truss_chain(TANGENT, testtype)
    ok_s, u_ref, _ = _truss_chain(STOCK, testtype)
    assert ok_t == 0 and ok_s == 0
    assert u_tan == pytest.approx(u_ref, rel=1e-9)
    assert u_tan == pytest.approx(EXPECTED, rel=1e-9)


def test_default_is_stock():
    """`-tangentPredictor` is opt-in: unflagged LadrunoLoadControl must match."""
    ok_l, u_lad, it_lad = _truss_chain(LADRUNO)
    ok_s, u_ref, it_ref = _truss_chain(STOCK)
    assert ok_l == 0 and ok_s == 0
    assert u_lad == pytest.approx(u_ref, rel=1e-12)
    assert it_lad == it_ref


@pytest.mark.parametrize("testtype", ["NormDispIncr", "EnergyIncr"])
def test_tangent_predictor_never_leaves_the_interior_unmoved(testtype):
    """THE SHARP ONE -- the S2 failure mode this route could re-create.

    Iteration 1 deliberately runs before the sps are enforced. If `-K*du_D` did
    not assemble, dU would be ~0, a dU-based test would converge immediately,
    and the step would commit with only the driven node displaced. Assert the
    interior actually moved, not merely that a number matched.
    """
    ok, u, iters = _truss_chain(TANGENT, testtype)
    assert ok == 0
    assert abs(u) > 1e-6, (
        "interior node did not move -- the K*du_prescribed term failed to "
        "assemble and the step committed on a false convergence")
    assert u == pytest.approx(EXPECTED, rel=1e-9)
    assert iters >= 1


def test_fallback_on_a_model_with_no_sp():
    """No non-homogeneous `sp` => nothing to assemble => must degrade to stock.

    This is the path that exercises `formUnbalance()`'s fallback. It must give
    the stock answer, not a step that silently never moved.
    """
    ok_t, u_tan, _ = _truss_chain(TANGENT, drive="load")
    ok_s, u_ref, _ = _truss_chain(STOCK, drive="load")
    assert ok_t == 0 and ok_s == 0
    assert u_tan == pytest.approx(u_ref, rel=1e-9)
    assert u_tan == pytest.approx(EXPECTED, rel=1e-6)


def test_extrapolate_and_tangent_predictor_refuse_to_compose():
    """Stacking them would double-count the interior motion; -extrapolate loses."""
    _truss_chain(("LadrunoLoadControl", 1.0, "-extrapolate", 1.0,
                  "-tangentPredictor"))
    assert ops.ladrunoLoadControl("extrapolate") == pytest.approx(0.0)
    assert ops.ladrunoLoadControl("tangentPredictor") == pytest.approx(1.0)


def test_tangent_predictor_reports_how_often_it_fired():
    """`armed` is meaningless for a stateless route; `tangentPredicts` is not.

    ADR-80 80c §2 records that a diagnostic counter is the only thing that
    distinguishes "the predictor helped little" from "the predictor never ran".
    """
    _truss_chain(TANGENT)
    assert ops.ladrunoLoadControl("tangentPredicts") >= 1


def _j2_chain(integrator_args, nsteps=5, maxiter=30):
    """1D stdBrick chain with LadrunoJ2, driven face. Returns (ok, u_mid, iters).

    Sized like the ADR-80 synthetic gates: converged sigma = 300 MPa against
    fy = 379.5, so NOTHING yields physically and the whole cost is spurious.
    """
    Lc, ax, emod, fy, hiso, delta, N = 100.0, 10.0, 200000.0, 379.5, 2000.0, 0.15, 10
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    h = Lc / N
    for k in range(N + 1):
        z, b = k * h, 4 * k
        ops.node(b + 1, 0.0, 0.0, z)
        ops.node(b + 2, ax, 0.0, z)
        ops.node(b + 3, ax, ax, z)
        ops.node(b + 4, 0.0, ax, z)
    ops.nDMaterial("LadrunoJ2", 1, emod / 3.0, emod / 2.0,
                   "-iso", "voce", fy, 0.0, 0.0, hiso)
    for k in range(N):
        b = 4 * k
        ops.element("stdBrick", k + 1, *[b + j for j in range(1, 9)], 1)
    for k in range(N + 1):
        b = 4 * k
        for j in range(1, 5):
            ops.fix(b + j, 1, 1, 1) if k == 0 else ops.fix(b + j, 1, 1, 0)

    top = 4 * N
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(1, 5):
        ops.sp(top + j, 3, delta)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1e-10, maxiter, 0)
    ops.algorithm("KrylovNewton")
    ops.analysis("Static")

    iters = 0
    for _ in range(nsteps):
        ops.integrator(*[a if a != 1.0 else 1.0 / nsteps for a in integrator_args])
        if ops.analyze(1) != 0:
            return -1, 0.0, iters
        iters += ops.testIter()
    # THE MIDPOINT, not the driven face. Reporting u at the sp'd node would
    # report the prescribed value back to itself and would pass even for a step
    # that moved nothing but the boundary layer -- the S2 failure signature.
    return 0, ops.nodeDisp(4 * (N // 2) + 1, 3), iters


def test_j2_answer_is_unchanged_and_the_path_is_cheaper():
    """The ADR-80 mechanism itself, at CI scale.

    Two claims, and the second one is the point of the whole phase: the answer
    must be identical, and the tangent route must cost STRICTLY FEWER Newton
    iterations than stock on a path-dependent law under prescribed displacement.
    Stated as an inequality, not a tuned count, so a mesh or tolerance change
    cannot turn it into a false alarm.
    """
    ok_s, u_ref, it_ref = _j2_chain(list(STOCK))
    ok_t, u_tan, it_tan = _j2_chain(list(TANGENT))
    assert ok_s == 0 and ok_t == 0
    # exact for this chain: the midpoint of a uniform bar at delta = 0.15
    assert u_ref == pytest.approx(0.075, rel=1e-6)
    assert u_tan == pytest.approx(u_ref, rel=1e-8)
    assert it_tan < it_ref, (
        f"-tangentPredictor cost {it_tan} iterations against stock's {it_ref}: "
        "the b -= K*du_D term is not doing its job")
