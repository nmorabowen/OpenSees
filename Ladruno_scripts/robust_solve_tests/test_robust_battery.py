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
import torture_stabilize as stab          # noqa: E402
import torture_dynamics as dyn            # noqa: E402


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


# --- rung-4 (auto-stabilization) is wired; verdict discipline is honest -----
# There is no clean-WIN fixture for rung-4 in this battery: pure softening defeats
# stabilize too (it IS load control), and a snap-through is dynamically jumped by
# adaptive load control before rung-4 is reached (R-LOG-MASK). So we assert the
# ESCALATION + honest verdict, and that rung-3 is preferred when a control DOF
# exists. The seam getters/gate/ramp-down are exercised by test_stabilize_* above.
def test_robust_drive_rung4_engages_on_softening():
    import opensees as ops
    soft.build()
    ops.integrator("LoadControl", soft.PEAK / 200.0)
    res = robust_drive.robust_drive(
        ops, done=lambda: soft.strain() <= -0.010,
        load_increment=soft.PEAK / 200.0, control=None, stabilize=True,
        stab_f=1.0e-3, max_substeps=4000)
    assert res.switches >= 1                 # escalated past plain load control
    assert res.mode == "Stabilized"          # reached rung-4
    assert not res.completed                 # stabilize cannot pass pure softening
    assert res.verdict in ("incomplete", "unverified")
    assert not bool(res)                     # NEVER a clean equilibrium


def test_robust_drive_rung3_preferred_over_rung4():
    # With a control DOF, rung-3 wins cleanly even though stabilize is enabled.
    import opensees as ops
    soft.build()
    ops.integrator("LoadControl", soft.PEAK / 200.0)
    res = robust_drive.robust_drive(
        ops, done=lambda: soft.strain() <= -0.010,
        load_increment=soft.PEAK / 200.0, control=(2, 1, -5.0e-5), stabilize=True)
    assert bool(res) and res.verdict == "equilibrium"
    assert res.mode == "DisplacementControl"   # rung-3, not rung-4


# --- rung-4 c-reduction diffusion bound (R-DIFFUSION) ------------------------
# The only trustworthy fidelity certificate for a stabilized run: re-run at half
# the viscosity and bound the peak-load drift. `robust_drive(verify_rebuild=…)`
# computes this internally on a rung-4 success and upgrades `unverified` →
# `regularized` iff drift ≤ verify_tol; `diffusion_drift` is the standalone form.
def test_diffusion_drift_discriminates():
    # pure math: an f-insensitive metric → ~0 drift; an f-sensitive one → large.
    flat = robust_drive.diffusion_drift(lambda f: 30.0, 1e-3)
    assert flat == 0.0
    sensitive = robust_drive.diffusion_drift(lambda f: 1.0 + 10.0 * f, 0.1)
    assert sensitive > 0.1                     # halving f moved the metric a lot


def test_stabilize_diffusion_bound_trustworthy_on_softening():
    # the stabilized strength peak is f-independent → drift below the 2% gate, so a
    # verify_rebuild=peak_load_stabilized rung-4 run would stamp `regularized`.
    drift = robust_drive.diffusion_drift(stab.peak_load_stabilized, 1e-3)
    assert drift <= 0.02


# --- rung-5 (dynamics fall-through) + the ladrunoDR command + R-HANDOFF ------
def test_ladrunoDR_command_surface():
    import opensees as ops
    rn, ke = dyn.dr_command_surface()
    assert rn >= 0.0 and ke >= 0.0           # the getters resolve on active DR
    # the command must reject a bogus subcommand (openseespy raises on stderr)
    with pytest.raises(ops.OpenSeesError):
        ops.ladrunoDR("bogus")


def test_robust_drive_rung5_dynamics_handoff():
    res, t_after = dyn.run_rung5()
    assert res.switches >= 1                  # escalated all the way to rung-5
    assert res.mode == "Dynamics"
    assert res.dr_settled                     # DR reached the mass-free settling gate
    # R-HANDOFF contract: lambda is frozen and restored EXACTLY across the excursion
    assert res.dr_lambda is not None
    assert t_after == res.dr_lambda
    # a DR rest state is regularized, never a clean equilibrium
    assert res.verdict == "regularized" and not bool(res)


# --- Layer-1.5: the stabilization-energy observability seam (rung-4 gate) ---
# These verify the C++ getters (LadrunoArcLength 33004) the driver reads to gate
# rung-4 and ramp the viscous coefficient (ADR-31 "gate vs audit", R-RAMPDOWN).
def test_stabilize_seam_reports_sane_energy():
    import opensees as ops
    r = stab.run_stabilized(adaptStab=True)
    assert r["committed"] > 0
    assert r["ref"] > 0.0                          # Estrain0 calibrated at commit
    assert r["dissip"] >= 0.0
    # identity: dissipationRatio == dissipatedEnergy / referenceEnergy
    assert r["ratio"] == pytest.approx(r["dissip"] / r["ref"], rel=1e-9)


def test_stabilize_adaptstab_bounds_dissipation():
    # -adaptStab holds the cumulative ratio near fTarget; without it the ratio
    # accumulates well above (the watchdog signal the driver gates on is real).
    a = stab.run_stabilized(adaptStab=True)
    b = stab.run_stabilized(adaptStab=False)
    assert a["ratio"] < 100 * stab.FTARGET        # O(fTarget), tight margin
    assert b["ratio"] > a["ratio"]                # cumulative drifts up


def test_stabilize_scaleCVisc_ramps_down():
    # R-RAMPDOWN: scaling cVisc by 0.1 cuts the per-window artificial dissipation.
    inc_full, inc_ramped = stab.measure_rampdown(factor=0.1)
    assert inc_full > 0.0
    assert inc_ramped < inc_full


def test_stabilize_scaleCVisc_rejects_nonpositive():
    import opensees as ops
    stab.build()
    stab._arm(adaptStab=False)
    ops.analyze(1)
    # factor must be > 0; the C++ guard warns to stderr, which openseespy raises.
    with pytest.raises(ops.OpenSeesError):
        ops.ladrunoArcLength("scaleCVisc", 0.0)
    with pytest.raises(ops.OpenSeesError):
        ops.ladrunoArcLength("scaleCVisc", -2.0)
