"""ADR-31 v1 acceptance gate — the robust-solve torture battery as pytest.

Asserts BOTH directions of the robustness claim on the live build:
  * plain / adaptive LOAD control cannot pass a strength peak or a limit point;
  * the rung-3 constraint switch (and the Layer-0 robust_drive) does.

Runs against the built opensees.pyd (via the fixtures' LADRUNO_DIST). Skips
cleanly if the build is absent, so the suite is safe in CI without a dist.
"""
import os

import pytest

DIST = os.environ.get("LADRUNO_DIST",
                      r"C:\Users\nmb\Documents\Github\OpenSees\dist\bin")
pytestmark = pytest.mark.skipif(
    not os.path.isdir(DIST),
    reason=f"no built opensees.pyd at {DIST} (set LADRUNO_DIST)")

import robust_drive                       # noqa: E402
import torture_softening as soft          # noqa: E402
import torture_snapthrough as snap        # noqa: E402


# --- softening: the clean load-control killer -----------------------------
def _soft_targets():
    lam_T = 1.3 * soft.PEAK
    return lam_T, lam_T / 200.0


def test_softening_loadcontrol_cannot_pass_peak():
    lam_T, dlam = _soft_targets()
    peak, eps_min, reached_soft = soft.run_loadcontrol(dlam, lam_T, 1000)
    assert not reached_soft                       # never got onto the softening branch
    assert peak <= soft.PEAK * 1.01               # stuck at the strength peak
    assert eps_min > 1.5 * soft.EPSC0             # strain never passed the peak


def test_softening_adaptive_loadcontrol_also_dies():
    lam_T, dlam = _soft_targets()
    _, eps, _ = soft.run_adaptive_loadcontrol(dlam, lam_T)
    assert eps > 1.5 * soft.EPSC0                 # halving+ladder cannot pass either


def test_softening_dispcontrol_traces_branch():
    eps_min, reached_soft, _ = soft.run_dispcontrol(-5.0e-5, 800)
    assert reached_soft and eps_min <= -0.009


def test_softening_robust_traces_branch():
    lam_T, dlam = _soft_targets()
    eps_min, reached_soft, switches = soft.run_robust(dlam, -5.0e-5, 800)
    assert reached_soft and switches == 1 and eps_min <= -0.009


# --- snap-through: load control stuck at the limit; rung-3 crosses ---------
def _snap_calibrate():
    _, _, lam_lim, _ = snap.run_dispcontrol(-1.0e-3, 600)
    lam_T = 1.5 * lam_lim
    return lam_T, lam_T / 300.0


def test_snapthrough_loadcontrol_stuck_at_limit():
    lam_T, dlam = _snap_calibrate()
    _, a_miny, _ = snap.run_loadcontrol(dlam, lam_T, 1200)
    assert not snap.crossed(a_miny)               # cannot drive the apex past flat


def test_snapthrough_robust_crosses_with_switch():
    _, dlam = _snap_calibrate()
    _, d_miny, d_rev, d_sw = snap.run_robust(dlam, -1.0e-3, 600)
    assert snap.crossed(d_miny) and d_rev and d_sw == 1


# --- the Layer-0 driver's own self-test -----------------------------------
def test_robust_drive_selftest_passes():
    assert robust_drive._selftest() == 0


# --- the driver must REFUSE rung-4 until the C++ getters ship (R-OBS) ------
def test_robust_drive_refuses_stabilize():
    with pytest.raises(NotImplementedError):
        robust_drive.robust_drive(
            None, done=lambda: True, load_increment=1.0, stabilize=True)
