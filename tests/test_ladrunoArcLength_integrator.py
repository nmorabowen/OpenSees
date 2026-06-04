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
  AL-5  -stabilize clears the limit point where plain LoadControl fails.
  AL-9  Layer-B reduceStep + revertToLastStep restores the COMMITTED per-step
        state ({deltaLambdaStep, deltaUstep, sign, currentLambda, arcLength})
        through a FAILED step, then the script-owned cut-and-retry recovers.
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


_LIMIT = 3.80   # von Mises snap-through limit load


# --------------------------------------------------------------------------
# AL-5 : -stabilize clears the limit point where plain LoadControl fails
# --------------------------------------------------------------------------
def test_al5_stabilize_clears_limit_point():
    # Stock LoadControl fails at lambda ~ 3.80 (no equilibrium on the near branch
    # above the limit). -stabilize regularizes K_T (K + cOverDt*M*) so the
    # load-incrementing Newton converges THROUGH the limit and the apex snaps to
    # the far branch. dLambda0 = sqrt(arc^2) = arc (load increment per step).
    arc, cV = 0.1, 0.01
    _build_truss()
    # load-incrementing snap-through needs a few more corrector iterations than the
    # smooth-path default (the regularized jump near the limit is nonlinear).
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-8, 100, 0)
    ops.algorithm("Newton")
    ops.integrator(LAL, arc, 1.0, "-stabilize", 2.0e-4, "-cVisc", cV)
    ops.analysis("Static")
    lam_reached, uy_final, n = 0.0, 0.0, 0
    for _ in range(200):
        if ops.analyze(1) != 0:
            break
        lam_reached = ops.getTime()
        uy_final = ops.nodeDisp(_APEX, 2)
        n += 1
    assert n > 0, "-stabilize produced no converged steps"
    # cleared the limit: converged at load factors stock LoadControl cannot reach.
    assert lam_reached > _LIMIT, (
        f"-stabilize stalled at lambda={lam_reached:.3f} (<= limit {_LIMIT})"
    )
    # and the regularized path reached the far stable branch (snapped through).
    assert uy_final < -2.0 * _H, f"apex did not reach the far branch: uy={uy_final:.4f}"


# --------------------------------------------------------------------------
# AL-9 : Layer-B reduceStep + revertToLastStep restores the COMMITTED state
#        through a failed step (ADR-20 §4.3 / review §2.10).
# --------------------------------------------------------------------------
def test_al9_layerb_reduce_revert_restores_committed_state():
    # The crux: stock ArcLength overwrites deltaLambdaStep/deltaUstep in newStep()
    # BEFORE its own b24ac<0 bail, so a failed step leaves them POLLUTED — "state
    # preserved for free" is false. LadrunoArcLength snapshots the per-step state
    # at commit() and an overridden revertToLastStep() (which StaticAnalysis::analyze
    # already calls on a failed step) restores it. This test drives a few converged
    # steps, forces a genuine failure with an oversized radius, and asserts the
    # COMMITTED state is restored — then that reduceStep()+retry recovers.
    arc, alpha = 0.05, 0.1
    _build_truss()
    _setup_solver()
    ops.integrator(LAL, arc, alpha)
    ops.analysis("Static")

    # Climb the rising branch toward the limit point (~3.80). A fixed radius this
    # large eventually fails to converge near the limit (Newton can't clear the
    # over-long step) — that natural failure IS the failed step we test. Record
    # the committed per-step state after every converged step.
    committed = None
    failed = False
    for _ in range(2000):
        ok = ops.analyze(1)
        if ok != 0:
            failed = True
            break
        committed = dict(
            lam=ops.getTime(),
            uy=ops.nodeDisp(_APEX, 2),
            dlds=ops.ladrunoArcLength("deltaLambdaStep"),
            sign=ops.ladrunoArcLength("sign"),
            nrm=ops.ladrunoArcLength("deltaUstepNorm"),
            al=ops.ladrunoArcLength("arcLength"),
        )
    assert failed, "fixed-radius arc-length never hit a failed step"
    assert committed is not None, "no converged step before the failure"
    assert committed["lam"] > 3.0, (
        f"failed too early to be the limit point (lambda={committed['lam']:.3f})"
    )
    assert abs(committed["al"] - arc) <= 1e-12, "radius drifted with no -adapt flag"
    assert abs(committed["nrm"]) > 0.0, "committed deltaUstep should be non-zero"

    # After the FAILED step StaticAnalysis::analyze already called
    # revertToLastCommit()+revertToLastStep(). The committed per-step state must be
    # restored — NOT the in-flight values newStep()/update() left polluted before
    # the bail (the ADR-20 §2.10 crux).
    assert ops.getTime() == pytest.approx(committed["lam"]), "domain load factor not reverted"
    assert ops.nodeDisp(_APEX, 2) == pytest.approx(committed["uy"]), "apex disp not reverted"
    assert ops.ladrunoArcLength("deltaLambdaStep") == pytest.approx(committed["dlds"]), (
        "deltaLambdaStep NOT restored — left polluted by the failed step"
    )
    assert ops.ladrunoArcLength("sign") == committed["sign"], (
        "signLastDeltaLambdaStep NOT restored"
    )
    assert ops.ladrunoArcLength("deltaUstepNorm") == pytest.approx(committed["nrm"]), (
        "deltaUstep NOT restored — left polluted by the failed step"
    )

    # Script-owned cut-and-retry (ADR-20 §4.3): shrink the radius and re-attempt
    # until the step clears. Each failed retry reverts the per-step state again but
    # leaves the (script-owned) radius alone, so reduceStep() accumulates.
    recovered = False
    for _ in range(14):
        ops.ladrunoArcLength("reduceStep", 0.5)
        if ops.analyze(1) == 0:
            recovered = True
            break
    assert recovered, "reduceStep cut-and-retry never cleared the failed step"
    # the radius genuinely shrank below the committed value (reduceStep accumulated)
    assert ops.ladrunoArcLength("arcLength") < committed["al"], (
        "reduceStep did not shrink the arc radius across the retry loop"
    )
    # and the recovered step advanced the path past the last committed state
    advanced = (ops.getTime() != pytest.approx(committed["lam"])) or (
        ops.nodeDisp(_APEX, 2) != pytest.approx(committed["uy"])
    )
    assert advanced, "retry converged but did not advance the path"


# --------------------------------------------------------------------------
# AL-S : parser smoke — all documented forms construct without error
# --------------------------------------------------------------------------
@pytest.mark.parametrize("cmd", [
    (LAL, 0.01, 1.0),                                   # identity form
    (LAL, 0.01, 0.1, "-adapt", 5, 1.0e-3, 1.0e-1),     # adaptive form
    (LAL, 0.01, 0.1, "-adapt", 5, 1.0e-3, 1.0e-1, "-p", 0.5),  # + Crisfield exp
    (LAL, 0.1, 1.0, "-stabilize", 2.0e-4),             # viscous (auto c)
    (LAL, 0.1, 1.0, "-stabilize", 2.0e-4, "-cVisc", 0.5),   # explicit c
    (LAL, 0.1, 1.0, "-stabilize", 2.0e-4, "-mass", "gershgorin"),
])
def test_als_parser_smoke(cmd):
    _build_truss()
    _setup_solver()
    ops.integrator(*cmd)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, f"first step failed for {cmd}"
    assert math.isfinite(ops.nodeDisp(_APEX, 2))
