"""revertToLastStep() on the Ladruno explicit integrators (CentralDifferenceLadruno
33003 + unified ExplicitBathe 33000).

The analysis framework calls Domain::revertToLastCommit() then
Integrator::revertToLastStep() when a step fails. Both explicit classes carry
private leap-frog / sub-step state that advances BEFORE their failure exits (their
own -divergence / NaN circuit breakers included), so without an override the
inherited no-op leaves the integrator desynced from the reverted Domain and an
adaptive driver that halves dt and retries continues from garbage.

The fix re-seeds the integrator from the COMMITTED node state and (for CD) re-arms
the first-step starter, making retry-after-abort exactly equal to a fresh restart
from the committed state:

  * CD leg is the change-detector: pre-fix, Ut/Vhalf hold the dt_bad-advanced
    garbage and firstStep stays consumed, so the continuation diverges from the
    clean reference immediately.
  * ExplicitBathe defers its cross-step advance to commit(), so its newStep/update
    aborts were already revert-safe; its leg is a REGRESSION LOCK (plus the
    getTime() assertion that the mid-step (1-p)dt advance is undone by
    revertToLastCommit). The genuinely-broken EB windows (commit-failure ordering,
    stale prevKE on a late updateDomain failure) are not scriptable from Zone A
    and are covered by the commit() reorder + code inspection.

Tolerances are TIGHT on purpose: for the undamped SDOF, restart == continuation up
to one floating-point round-trip on the half-step velocity (v_{n+1/2} reconstructed
as (v_{n+1/2} + dt/2 a) - dt/2 a), so rel 1e-9 over a few hundred steps is the
honest bound — a physically-loose band would pass on a half-broken revert.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"
EB = "ExplicitBathe"

K, M = 100.0, 1.0                 # omega = 10 -> CD dt_cr = 2/omega = 0.2
OMEGA = math.sqrt(K / M)
DT_OK = 0.005
DT_BAD = 0.5                      # >> dt_cr: one step trips the -divergence breaker
N_BEFORE, N_AFTER = 200, 400

# CIRCULAR-MOTION fixture (load-bearing): a plain free-vibration SDOF passes
# through v ~ 0 at every velocity trough, where the -divergence KE proxy's
# per-step ratio ke/prevKE spikes unboundedly (sin^2 growth off a near-zero
# floor) and FALSE-trips the breaker mid-run. Equal springs in x and y with the
# quadrature seed (v_x0 = A*omega, u_y0 = A) give vx = A*omega*cos, vy =
# -A*omega*sin -> the TOTAL KE is CONSTANT, the stable-phase ratio is ~1, and
# only the genuine dt_bad fault step (KE jump ~O(100x)) trips the breaker.
A0 = 0.1                          # motion amplitude; v amplitude = A0*omega = 1


def _build(integrator_args, seed=None):
    """seed: (ux, uy, vx, vy, ax, ay) committed initial state; default = circular
    motion (accel = -omega^2 * u, consistent). The accel seed matters for the
    ExplicitBathe reference runs: Noh-Bathe is SELF-STARTING from (u, v, a), so a
    fresh run meant to reproduce a restart must carry the committed a too (CD's
    starter re-derives a0 from the residual, so it ignores the seed — harmless)."""
    if seed is None:
        seed = (0.0, A0, A0 * OMEGA, 0.0, 0.0, -A0 * OMEGA * OMEGA)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0)
    ops.mass(2, M, M)
    ops.uniaxialMaterial("Elastic", 1, K)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, 1, "-dir", 1, 2)
    ux, uy, vx, vy, ax, ay = seed
    for dof, u, v, a in ((1, ux, vx, ax), (2, uy, vy, ay)):
        if u != 0.0:
            ops.setNodeDisp(2, dof, u, "-commit")
        if v != 0.0:
            ops.setNodeVel(2, dof, v, "-commit")
        if a != 0.0:
            ops.setNodeAccel(2, dof, a, "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")


def _state():
    return (ops.nodeDisp(2, 1), ops.nodeDisp(2, 2),
            ops.nodeVel(2, 1), ops.nodeVel(2, 2),
            ops.nodeAccel(2, 1), ops.nodeAccel(2, 2))


def _run(nsteps, dt):
    """Returns list of (t, (ux,uy,vx,vy)) AFTER each successful step; stops on failure."""
    out = []
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        out.append((ops.getTime(), _state()))
    return out


def _assert_matches(cont, ref, tol, what):
    assert len(cont) == len(ref), (
        "%s: continuation ran %d steps vs reference %d (a post-revert integrator "
        "poisoned enough to diverge aborts early)" % (what, len(cont), len(ref))
    )
    names = ("ux", "uy", "vx", "vy", "ax", "ay")
    for i, ((tf, sf), (tr, sr)) in enumerate(zip(cont, ref)):
        for got, exp, name in zip(sf, sr, names):
            # abs tolerance scaled to the motion amplitude (u ~ A0, v ~ A0*omega = 1):
            # components pass through zero, so a pure relative check is meaningless there.
            assert abs(got - exp) < tol, (
                "%s: %s diverges from the reference at continuation step %d: got %.12g "
                "expected %.12g (|diff| %.2e >= %g) — revertToLastStep left stale "
                "integrator state" % (what, name, i, got, exp, abs(got - exp), tol)
            )


def _integrator(name):
    if name == CDL:
        return (CDL, "-divergence", 10.0)
    return (EB, 0.54, "-divergence", 10.0)


@pytest.mark.parametrize("name", [CDL, EB])
def test_abort_then_retry_matches_clean_run(name):
    """One supra-stable step trips the KE breaker; the continuation at the original
    dt must match a never-failed clean run step-for-step."""
    args = _integrator(name)

    _build(args)
    before = _run(N_BEFORE, DT_OK)
    assert len(before) == N_BEFORE, "fixture: the pre-fault phase must run clean"
    t_committed = before[-1][0]
    rc_bad = ops.analyze(1, DT_BAD)
    assert rc_bad != 0, (
        "fixture: the dt=%g step must trip the -divergence breaker (analyze<0); "
        "got %r — the fault leg never happened" % (DT_BAD, rc_bad)
    )
    # Domain::revertToLastCommit restored the committed time (for ExplicitBathe this
    # also locks the mid-step p*dt / (1-p)*dt advances being undone).
    assert ops.getTime() == pytest.approx(t_committed, abs=1e-12), (
        "domain time %.12g != committed %.12g after the failed step — "
        "revertToLastCommit did not restore time" % (ops.getTime(), t_committed)
    )
    cont = _run(N_AFTER, DT_OK)

    _build(args)
    clean = _run(N_BEFORE + N_AFTER, DT_OK)
    assert len(clean) == N_BEFORE + N_AFTER, "fixture: the clean reference must run"

    _assert_matches(cont, clean[N_BEFORE:], 1e-9, "%s retry-at-same-dt" % name)


@pytest.mark.parametrize("name", [CDL, EB])
def test_abort_then_retry_at_smaller_dt(name):
    """Adaptive-driver shape: after the abort, continue at dt/2. The reference is a
    FRESH analysis seeded with the committed state — exactly what a correct revert
    reduces the retry to (the CD starter re-derives a_n and re-centers the half-step
    for the NEW dt, so restart == the seeded fresh run, not the never-failed
    variable-dt continuation)."""
    args = _integrator(name)
    dt_retry = DT_OK / 2.0

    _build(args)
    before = _run(N_BEFORE, DT_OK)
    assert len(before) == N_BEFORE, "fixture: the pre-fault phase must run clean"
    committed = _state()
    assert ops.analyze(1, DT_BAD) != 0, "fixture: fault leg must abort"
    cont = _run(N_AFTER, dt_retry)

    _build(args, seed=committed)              # fresh run from the committed state
    ref = _run(N_AFTER, dt_retry)
    assert len(ref) == N_AFTER, "fixture: the seeded reference must run"

    _assert_matches(cont, ref, 1e-9, "%s retry-at-dt/2" % name)


