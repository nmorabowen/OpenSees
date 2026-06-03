"""LadrunoArcLength (Ladruno adaptive arc-length static integrator, classTag
33004) — Zone-A battery.

PR-1 scope = stock-ArcLength identity + Ramm desired-iteration radius
adaptation ("-adapt Jd ellMin ellMax"). The viscous "-stabilize" mode and the
Layer-B reduceStep/revert mutators are deferred to follow-up PRs.

Canonical model: a shallow von Mises truss (two corotational bars, apex loaded
downward). E*A is tuned so the snap-through limit load is ~3.8. This is the
smallest problem that exhibits a limit point; it is the standing prototype
fixture from tests/_proto_arch_snapthrough.py, promoted here.

  AL-1  LadrunoArcLength (no flags) reproduces stock ArcLength step-for-step.
  AL-3  -adapt grows the radius (desired-iteration rule) => reaches a target
        load in fewer steps than the same fixed radius.
  AL-S  smoke: parser accepts the documented argument forms.

Theory / design: Ladruno_implementation/20_ladruno_arclength_stabilized_adr.md.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

LAL = "LadrunoArcLength"

# von Mises truss parameters (limit load ~3.80)
_H = 0.10
_SPAN = 1.0
_E = 1.0e4
_A = 1.0
_APEX = 2


def _build_truss():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, -_SPAN, 0.0)
    ops.node(_APEX, 0.0, _H)
    ops.node(3, _SPAN, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(3, 1, 1)
    ops.fix(_APEX, 1, 0)            # apex: X fixed by symmetry, Y free
    ops.uniaxialMaterial("Elastic", 1, _E)
    ops.element("corotTruss", 1, 1, _APEX, _A, 1)
    ops.element("corotTruss", 2, _APEX, 3, _A, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(_APEX, 0.0, -1.0)


def _setup_solver():
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")


def _trace(integrator_cmd, nsteps):
    """Run a static analysis, return the list of (lambda, uy) for each
    converged step. integrator_cmd is the full integrator(...) arg tuple."""
    _build_truss()
    _setup_solver()
    ops.integrator(*integrator_cmd)
    ops.analysis("Static")
    path = []
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        path.append((ops.getTime(), ops.nodeDisp(_APEX, 2)))
    return path


# --------------------------------------------------------------------------
# AL-1 : identity — LadrunoArcLength (no flags) == stock ArcLength
# --------------------------------------------------------------------------
def test_al1_identity_vs_arclength():
    arc, alpha, nsteps = 0.01, 0.1, 25
    ref = _trace(("ArcLength", arc, alpha), nsteps)
    lad = _trace((LAL, arc, alpha), nsteps)

    assert len(ref) > 10, "reference ArcLength produced too few steps"
    assert len(lad) == len(ref), (
        "LadrunoArcLength traced a different number of steps "
        f"({len(lad)}) than ArcLength ({len(ref)})"
    )
    for i, ((lr, ur), (ll, ul)) in enumerate(zip(ref, lad)):
        assert abs(ll - lr) <= 1e-9 + 1e-7 * abs(lr), f"lambda mismatch at step {i}: {ll} vs {lr}"
        assert abs(ul - ur) <= 1e-9 + 1e-7 * abs(ur), f"uy mismatch at step {i}: {ul} vs {ur}"


# --------------------------------------------------------------------------
# AL-3 : -adapt grows the radius => fewer steps to a target load
# --------------------------------------------------------------------------
def _steps_to_reach(path, target_lambda):
    for i, (lam, _uy) in enumerate(path):
        if lam >= target_lambda:
            return i + 1
    return None


def test_al3_adapt_reaches_target_in_fewer_steps():
    # both start from the SAME small radius; the adaptive one is allowed to
    # grow it up to ellMax via the desired-iteration rule.
    arc, alpha, target = 0.004, 0.1, 3.0

    fixed = _trace((LAL, arc, alpha), 400)
    adapt = _trace((LAL, arc, alpha, "-adapt", 8, 1.0e-3, 5.0e-2), 400)

    n_fixed = _steps_to_reach(fixed, target)
    n_adapt = _steps_to_reach(adapt, target)

    assert n_fixed is not None, "fixed-radius run never reached the target load"
    assert n_adapt is not None, "adaptive run never reached the target load"
    # desired-iteration adaptation (Jd=8 > the 2-4 iters Newton actually needs)
    # grows the radius, so it should reach the target in strictly fewer steps.
    assert n_adapt < n_fixed, (
        f"adaptive run did not accelerate: {n_adapt} steps vs fixed {n_fixed}"
    )


# --------------------------------------------------------------------------
# AL-S : parser smoke — all documented forms construct without error
# --------------------------------------------------------------------------
@pytest.mark.parametrize("cmd", [
    (LAL, 0.01, 1.0),                                   # identity form
    (LAL, 0.01, 0.1, "-adapt", 5, 1.0e-3, 1.0e-1),     # adaptive form
    (LAL, 0.01, 0.1, "-adapt", 5, 1.0e-3, 1.0e-1, "-p", 0.5),  # + Crisfield exp
])
def test_als_parser_smoke(cmd):
    _build_truss()
    _setup_solver()
    ops.integrator(*cmd)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, f"first step failed for {cmd}"
    assert math.isfinite(ops.nodeDisp(_APEX, 2))
