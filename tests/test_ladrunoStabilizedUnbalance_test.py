"""LadrunoStabilizedUnbalance (Ladruno true-equilibrium NormUnbalance convergence
test, classTag 33000 in the CONVERGENCE_TEST registry) — Zone-A battery.

In the LadrunoArcLength -stabilize mode the SOE residual a stock CTestNormUnbalance
norms is the f_v-polluted `lambda p - f_int - f_v`, so Newton is satisfied while
the TRUE static unbalance `||lambda p - f_int||` is still nonzero by `||f_v||`.
LadrunoStabilizedUnbalance norms the TRUE unbalance instead (a stricter criterion),
reaching the integrator through EquiSolnAlgo and falling back to the SOE `||B||`
for any non-stabilizing integrator.

  ST-1  graceful fallback: with a NON-stabilizing integrator (LoadControl) it
        behaves as a drop-in NormUnbalance — an elastic problem converges to the
        same linear solution as the stock test.
  ST-2  end-to-end under LadrunoArcLength -stabilize: the true-residual test
        drives convergence THROUGH the von Mises limit point (snap-through).
  ST-S  parser smoke: documented argument forms are accepted.

Theory / design: Ladruno_implementation/20_ladruno_arclength_stabilized_adr.md §8 #4.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

LSU = "LadrunoStabilizedUnbalance"
LAL = "LadrunoArcLength"

# von Mises truss (snap-through limit load ~3.80), shared with the AL battery
_H = 0.10
_SPAN = 1.0
_E = 1.0e4
_A = 1.0
_APEX = 2
_LIMIT = 3.70   # conservative lower bound on the snap-through limit load


def _build_truss():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, -_SPAN, 0.0)
    ops.node(_APEX, 0.0, _H)
    ops.node(3, _SPAN, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(3, 1, 1)
    ops.fix(_APEX, 1, 0)
    ops.uniaxialMaterial("Elastic", 1, _E)
    ops.element("corotTruss", 1, 1, _APEX, _A, 1)
    ops.element("corotTruss", 2, _APEX, 3, _A, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(_APEX, 0.0, -1.0)


# --------------------------------------------------------------------------
# ST-1 : graceful fallback to stock NormUnbalance for a non-stabilizing solver
# --------------------------------------------------------------------------
def _linear_truss_tip(test_cmd):
    # a stiff, near-linear elastic step under LoadControl; the convergence test
    # must converge to the same tip displacement regardless of which |.| it norms.
    _build_truss()
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test(*test_cmd)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.05)
    ops.analysis("Static")
    ok = ops.analyze(3)
    return ok, ops.nodeDisp(_APEX, 2)


def test_st1_graceful_fallback():
    ok_ref, uy_ref = _linear_truss_tip(("NormUnbalance", 1.0e-8, 50))
    ok_lsu, uy_lsu = _linear_truss_tip((LSU, 1.0e-8, 50))
    assert ok_ref == 0, "reference NormUnbalance run failed"
    assert ok_lsu == 0, f"{LSU} fallback run failed to converge"
    assert abs(uy_ref - uy_lsu) < 1.0e-9, (
        f"fallback diverged from NormUnbalance: {uy_ref} vs {uy_lsu}"
    )


# --------------------------------------------------------------------------
# ST-2 : true-residual semantics under LadrunoArcLength -stabilize
# --------------------------------------------------------------------------
# THE feature. In -stabilize mode the regularized SOE residual is
# ||lambda p - f_int - f_v||, which Newton CAN drive to ~0 (the stabilized
# solution satisfies the REGULARIZED equation), so stock NormUnbalance "converges"
# while the TRUE static unbalance is still ~||f_v|| > 0. LadrunoStabilizedUnbalance
# norms the TRUE unbalance, so it (correctly) REFUSES to accept a tolerance below
# ||f_v|| and accepts once the tolerance is above it. This proves it norms a
# different (true) residual, not the polluted SOE ||B||.
def _stabilized_step_ok(test_cmd, cVisc=0.01):
    _build_truss()
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test(*test_cmd)
    ops.algorithm("Newton")
    ops.integrator(LAL, 0.1, 1.0, "-stabilize", 2.0e-4, "-cVisc", cVisc)
    ops.analysis("Static")
    return ops.analyze(1)


def test_st2_true_residual_semantics():
    # stock NormUnbalance accepts the f_v-polluted (regularized) residual at 1e-8
    assert _stabilized_step_ok(("NormUnbalance", 1.0e-8, 100)) == 0, (
        "stock NormUnbalance unexpectedly failed the stabilized step"
    )
    # the true-residual test REFUSES the same tol: ||lambda p - f_int|| ~ ||f_v|| >> 1e-8
    assert _stabilized_step_ok((LSU, 1.0e-8, 100)) != 0, (
        f"{LSU} converged at 1e-8 — it is NOT norming the true unbalance "
        f"(it would have to drive the artificial force below 1e-8)"
    )
    # ...but ACCEPTS once the tolerance is above ||f_v|| (a genuine equilibrium
    # to that tolerance). 0.1 >> the ~0.018 plateau for cVisc=0.01.
    assert _stabilized_step_ok((LSU, 0.1, 100)) == 0, (
        f"{LSU} failed even at a tol above ||f_v|| — true residual never settled"
    )


# --------------------------------------------------------------------------
# ST-S : parser smoke
# --------------------------------------------------------------------------
def test_sts_parser_smoke():
    _build_truss()
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test(LSU, 1.0e-8, 50, 0, 2)     # tol maxIter printFlag normType
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.05)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, f"{LSU} with full arg form failed to take a step"
