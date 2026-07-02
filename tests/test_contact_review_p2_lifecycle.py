"""Contact review fix P2 (2026-07) — friction-state lifecycle: revert + reset.

Two lifecycle holes in the NTS friction path-state, both review-verified:

1. REVERT (the C3.2 MAJOR-2 pattern, never back-ported from mortar/edge to the NTS
   lane it was copied from): `FrictionState.gT0/engaged` are captured in getResidual
   (mid-iteration), but `revertToLastCommit()` only restored `gpTtrial=gpT` — a
   REJECTED implicit trial step latched an engagement origin captured at the rejected
   configuration, permanently offsetting every later stick traction. The fix adds the
   `gT0committed/engagedCommitted` double-buffer (commit promotes, revert restores).
   Gate: a step that FAILS while first-engaging friction, then retried at dt/2, must
   produce the BIT-SAME trajectory as a reference run that never saw the failed
   attempt (post-revert state is bitwise the committed state).

2. RESET: `Domain::revertToStart` (ops.reset) rewound nodes/elements but the contact
   engine had NO revertToStart at all — committed friction slip gpT, engagement
   origins, ALM multipliers, edge signs and re-emit state all survived, so a re-run
   after reset started from the previous run's friction history (a large spurious
   backward stick force at first contact). The fix drops all path-state maps (slots
   lazily recreate zeroed — by construction the pristine state of a fresh run).
   Gate: run a frictional slide (accumulates gpT), reset, re-run the identical
   loading — the two trajectories must match bit-tight.
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6
KT = 1.0e5
MU = 0.3


def _flat_master(L=20.0):
    """One long flat quad master strip in the z=0 plane, all nodes fixed."""
    ops.node(100, -L, -L, 0.0); ops.fix(100, 1, 1, 1)
    ops.node(101,  L, -L, 0.0); ops.fix(101, 1, 1, 1)
    ops.node(102,  L,  L, 0.0); ops.fix(102, 1, 1, 1)
    ops.node(103, -L,  L, 0.0); ops.fix(103, 1, 1, 1)
    return [100, 101, 102, 103]


# --------------------------------------------------------------- 1. failed-step revert
def _implicit_drag(fail_first):
    """Slave starts SEPARATED (gap +1e-4) so friction first ENGAGES mid-step at a
    laterally-moved trial config; a combined lateral+downward load drives it into
    contact. If fail_first: attempt one dt=1e-3 step under an impossible tolerance
    (friction engages during the rejected trials -> revert fires), then retry at
    dt/2 with the workable tolerance. Else: run dt/2 with the workable tolerance from
    the start. Returns the slave x-trajectory over 20 dt/2 steps."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _flat_master()
    ops.node(1, 0.0, 0.0, 1.0e-4)          # separated: engagement happens MID-STEP
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, MU, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 600.0, 0.0, -1000.0)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Newton")
    ops.analysis("Transient")

    if fail_first:
        ops.test("NormDispIncr", 1e-30, 1, 0)      # impossible => the dt step aborts
        rc = ops.analyze(1, 1.0e-3)
        assert rc != 0, "the poisoned step should have failed to converge"
        info = ops.ladrunoContactInfo()            # [numContacts, numCommits, numReverts]
        assert info[2] >= 1, f"revert hook must fire on the aborted step, contactInfo={info}"

    ops.test("NormDispIncr", 1e-10, 50, 0)
    traj = []
    for _ in range(20):
        rc = ops.analyze(1, 0.5e-3)
        assert rc == 0, "retry/reference step failed to converge"
        traj.append(ops.nodeDisp(1, 1))
    # the run must actually engage friction (anti-vacuous: slave in contact + dragged)
    z_dist = 1.0e-4 + ops.nodeDisp(1, 3)
    assert z_dist < 0.0, f"slave never engaged contact (z-dist={z_dist})"
    return traj


def test_failed_step_does_not_latch_stale_engagement_origin():
    """A rejected implicit step that FIRST-ENGAGES friction must leave no trace: the
    dt/2 retry must match, bit-tight, a reference run that never saw the failed
    attempt. Pre-fix, `engaged/gT0` survive the revert latched at the rejected trial
    config, so the retry's stick traction carries a permanent spurious offset."""
    traj_ref = _implicit_drag(fail_first=False)
    traj_retry = _implicit_drag(fail_first=True)
    err = max(abs(a - b) for a, b in zip(traj_retry, traj_ref))
    assert err < 1.0e-12, (
        f"stale engagement origin after a rejected step: max trajectory diff {err}")


# --------------------------------------------------------------------------- 2. reset
def _build_slide_model():
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mt = _flat_master()
    ops.node(1, 0.0, 0.0, -1.0e-3)          # penetrating: engaged from step 1
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *mt)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, MU, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 600.0, 0.0, -1000.0)         # slides: Fx > mu*Fz


def _build_analysis():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")


def _run_slide(nsteps=500):
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)
    traj = []
    for _ in range(nsteps):
        rc = ops.analyze(1, dt)
        assert rc == 0, "explicit slide step failed"
        traj.append(ops.nodeDisp(1, 1))
    return traj


def test_reset_drops_friction_history():
    """ops.reset() must return contact to the pristine t=0 state: a re-run of the
    identical frictional slide must reproduce the first run bit-tight. Pre-fix the
    committed plastic slip gpT (accumulated over the slide) survived the reset, so
    the re-run started with a large spurious backward stick force."""
    ops.wipe()
    _build_slide_model()
    _build_analysis()
    traj1 = _run_slide()
    # the slide must actually have slid (anti-vacuous: gpT accumulated)
    assert traj1[-1] > 1.0e-3, f"run 1 never slid (x_final={traj1[-1]})"

    ops.reset()                              # Domain::revertToStart
    ops.wipeAnalysis()                       # fresh integrator/handler; the contact
    _build_analysis()                        # ENGINE lives on the Domain and survives
    traj2 = _run_slide()

    err = max(abs(a - b) for a, b in zip(traj2, traj1))
    assert err < 1.0e-12, (
        f"contact state survived ops.reset(): max trajectory diff {err} "
        f"(run1 x_final={traj1[-1]}, run2 x_final={traj2[-1]})")
