"""ADR-52 W1-I1b — half-increment-residual accuracy gate (runtime validation).

The W1-I1a transient lane sizes dt from the iteration count (a COST criterion); W1-I1b
adds the Abaqus-style half-increment-residual ACCURACY gate. OpenSeesPy alone cannot form
the residual at an injected half-step state (setNodeDisp triggers no element->update(), so
reactions()/printB() report a residual that ignores the displacement-dependent internal
force). Two fork commands close that gap:

  * ladrunoSetNodeTrial nodeTag <ndof disp><ndof vel><ndof accel>  -- full-vector trial setter
    (the per-dof setNodeDisp can't build a multi-dof trial state; it resets to committed each call)
  * ladrunoTrialResidualNorm <loadTime>  -- Domain::update() then formUnbalance(), returns the
    inf-norm of the free-DOF dynamic unbalance at the current trial state, loads optionally
    re-applied at loadTime (the step midpoint). Committed state is never touched.

Decisive checks:
  1. COMMANDS EXIST; the residual at a converged step is ~0.
  2. ORACLE: the assembled half-increment residual matches an independent numpy computation
     (internal k*u_h + inertia m*a_h, load at the MIDPOINT) -- validates update(),
     formUnbalance(), and the load-at-midpoint applyLoad path.
  3. FULL-VECTOR SETTER: on a 2-dof node ladrunoSetNodeTrial sets BOTH components; a per-dof
     setNodeDisp loses the dominant one (proving why the full-vector setter is needed).
  4. COMMITTED STATE UNTOUCHED: nodeDisp is identical before/after a gate evaluation.
  5. GATE EFFICACY: where convergence-only sizing lets dt grow past accuracy, the error gate
     holds the solution close to a fine-dt reference (and far closer than gate-off).
"""
import os
import sys

import pytest

from _testbed import ops

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "Ladruno_scripts"))
from robust_drive import robust_transient, _half_increment_residual  # noqa: E402

pytestmark = [pytest.mark.zone_a]

E, A, L, M = 1000.0, 1.0, 1.0, 1.0      # k = EA/L = 1000, omega = sqrt(k/m) ~ 31.62, T ~ 0.1987
K = E * A / L


def _build_sdof(load_kind="constant", slope=1.0, p0=1.0):
    """A 1-dof mass-on-truss SDOF, Newmark average-acceleration, undamped.

    load_kind 'constant' -> Constant series factor p0; 'ramp' -> Linear series (P(t)=slope*t)."""
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, L)
    ops.fix(1, 1)
    ops.mass(2, M)
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.element("Truss", 1, 1, 2, A, 1)
    if load_kind == "constant":
        ops.timeSeries("Constant", 1, "-factor", p0)
    else:
        ops.timeSeries("Linear", 1, "-factor", slope)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGen")
    ops.test("NormDispIncr", 1e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")


def test_commands_exist():
    assert hasattr(ops, "ladrunoTrialResidualNorm")
    assert hasattr(ops, "ladrunoSetNodeTrial")
    _build_sdof()
    ops.analyze(1, 0.005)
    r = float(ops.ladrunoTrialResidualNorm())    # at the converged state -> ~0
    assert abs(r) < 1e-6


def test_halfres_matches_numpy_oracle():
    """Per-step assembled half-increment residual vs an independent numpy computation.

    A ramp (Linear) load makes the load value differ between the midpoint and t_{n+1}, so this
    also proves the helper applies the load at the MIDPOINT (the loadTime argument), not t_{n+1}."""
    slope = 5.0
    _build_sdof(load_kind="ramp", slope=slope)
    dt = 0.01
    max_abs = 0.0
    for _ in range(40):
        u0, v0, a0 = ops.nodeDisp(2, 1), ops.nodeVel(2, 1), ops.nodeAccel(2, 1)
        assert ops.analyze(1, dt) == 0
        u1, v1, a1 = ops.nodeDisp(2, 1), ops.nodeVel(2, 1), ops.nodeAccel(2, 1)
        # Independent check that the integrator IS average-acceleration (so the midpoint formula
        # below is the right interpolation family): the avg-accel full-step identities must hold
        # for the states OpenSees actually produced -- this validates the kinematics separately
        # from the formula the driver shares with this test.
        assert abs(u1 - (u0 + dt * v0 + 0.25 * dt * dt * (a0 + a1))) < 1e-9 * (abs(u1) + 1e-6)
        assert abs(v1 - (v0 + 0.5 * dt * (a0 + a1))) < 1e-9 * (abs(v1) + 1e-6)
        # constant-average-acceleration half-step state (same formula the driver uses)
        a_h = 0.5 * (a0 + a1)
        v_h = v0 + 0.25 * dt * (a0 + a1)
        u_h = u0 + 0.5 * dt * v0 + (dt * dt / 16.0) * (a0 + a1)
        t_mid = ops.getTime() - 0.5 * dt
        p_mid = slope * t_mid                      # Linear series value at the midpoint
        r_oracle = abs(p_mid - (K * u_h + M * a_h))
        # the fork path: inject the half state, read the assembled residual at the midpoint
        pre = {1: ([0.0], [0.0], [0.0]), 2: ([u0], [v0], [a0])}
        post = {1: ([0.0], [0.0], [0.0]), 2: ([u1], [v1], [a1])}
        r_fork = _half_increment_residual(ops, pre, post, dt)
        max_abs = max(max_abs, abs(r_fork - r_oracle))
    assert max_abs < 1e-6, f"half-increment residual oracle mismatch: max_abs={max_abs:.2e}"


def test_full_vector_trial_setter_multidof():
    """ladrunoSetNodeTrial sets a FULL multi-dof trial vector; per-dof setNodeDisp cannot.

    kx*ux is the DOMINANT term and dof-1 is the one the per-dof sequence loses -> the per-dof
    residual collapses to ky*|uy|, decisively smaller than the full-vector residual."""
    kx, ky = 4000.0, 1000.0
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)                          # free node, 2 dof
    ops.node(2, 1.0, 0.0)
    ops.node(3, 0.0, 1.0)
    ops.fix(2, 1, 1)
    ops.fix(3, 1, 1)
    ops.uniaxialMaterial("Elastic", 11, kx)
    ops.uniaxialMaterial("Elastic", 12, ky)
    ops.element("Truss", 1, 1, 2, 1.0, 11)         # horizontal -> stiffness kx in x
    ops.element("Truss", 2, 1, 3, 1.0, 12)         # vertical   -> stiffness ky in y
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGen")
    ops.test("NormDispIncr", 1e-12, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.0)
    ops.analysis("Static")
    ops.analyze(1)                                 # establish DOF numbering (no load -> u=0)

    ux, uy = 0.003, -0.002
    # FULL-VECTOR setter: both components land in the trial state.
    ops.ladrunoSetNodeTrial(1, ux, uy, 0.0, 0.0, 0.0, 0.0)
    r_full = float(ops.ladrunoTrialResidualNorm())  # static: ||P - F_int||inf = max(kx|ux|, ky|uy|)
    expected = max(kx * abs(ux), ky * abs(uy))      # = kx*ux = 12.0 (dominant)
    assert abs(r_full - expected) / expected < 1e-8, (r_full, expected)

    # per-dof setNodeDisp can express only ONE component: the second call resets dof-1 to committed.
    ops.ladrunoSetNodeTrial(1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)   # reset trial to committed (0)
    ops.setNodeDisp(1, 1, ux)
    ops.setNodeDisp(1, 2, uy)                       # resets dof-1 (the dominant term) back to 0
    r_perdof = float(ops.ladrunoTrialResidualNorm())
    assert abs(r_perdof - ky * abs(uy)) / (ky * abs(uy)) < 1e-8     # only the y term survives
    assert r_perdof < 0.5 * r_full                  # the full-vector setter captured strictly more


def test_committed_state_untouched():
    """A gate evaluation must not perturb the committed trajectory."""
    _build_sdof(load_kind="constant", p0=1.0)
    dt = 0.01
    for _ in range(5):
        ops.analyze(1, dt)
    u_before, v_before, a_before = ops.nodeDisp(2, 1), ops.nodeVel(2, 1), ops.nodeAccel(2, 1)
    pre = {1: ([0.0], [0.0], [0.0]), 2: ([u_before], [v_before], [a_before])}
    _half_increment_residual(ops, pre, pre, dt)    # bogus-but-harmless eval (post == pre)
    assert ops.nodeDisp(2, 1) == u_before          # committed disp untouched
    assert ops.nodeVel(2, 1) == v_before
    assert ops.analyze(1, dt) == 0                 # analysis still advances afterwards


def _fixed_dt_disp(dt_fine, t_end, p0):
    """Fine fixed-dt Newmark reference: return the displacement at t_end."""
    _build_sdof(load_kind="constant", p0=p0)
    n = int(round(t_end / dt_fine))
    for _ in range(n):
        ops.analyze(1, dt_fine)
    return ops.nodeDisp(2, 1)


def test_error_gate_improves_accuracy():
    """The gate keeps the solution near a fine reference where convergence-only sizing drifts."""
    p0, t_end = 1.0, 0.15
    amp = 2.0 * p0 / K                              # oscillation amplitude (1-cos response)
    u_ref = _fixed_dt_disp(1.0e-4, t_end, p0)

    # gate OFF (W1-I1a): linear problem converges in ~1 iter -> dt grows to the cap -> inaccurate
    _build_sdof(load_kind="constant", p0=p0)
    off = robust_transient(ops, t_end, 2.0e-3, dt_max=2.0e-3 * 64, target_iters=5)
    u_off = ops.nodeDisp(2, 1)

    # gate ON (W1-I1b): the half-increment residual holds dt where accuracy demands
    _build_sdof(load_kind="constant", p0=p0)
    on = robust_transient(ops, t_end, 2.0e-3, dt_max=2.0e-3 * 64, target_iters=5,
                          error_gate=True, haftol=0.02)
    u_on = ops.nodeDisp(2, 1)

    assert bool(off) and bool(on)
    err_off, err_on = abs(u_off - u_ref), abs(u_on - u_ref)
    assert err_on < err_off, (err_on, err_off, u_ref)        # gate is materially more accurate
    assert err_on < 0.25 * amp, (err_on, amp)                # and bounded well within amplitude
    assert on.max_dt_used < off.max_dt_used                  # gate held dt below the cap
    assert on.halfres_max > 0.0                              # the gate actually measured residuals
    assert on.verdict in ("accuracy_gated", "integrated")


def test_rate_dependent_material_not_corrupted():
    """Guards the ops_Dt fix: a rate-dependent material (Maxwell reads the global ops_Dt inside
    update()) must produce a FINITE, sane half-increment residual and integrate under the gate
    without blow-up. The bug this guards: the midpoint applyLoad sets ops_Dt = t_mid - t_{n+1} < 0,
    so update() would integrate the relaxation with a negative dt (exp(-ops_Dt/tR) -> growing)."""
    K, C, alpha, Lvar = 1000.0, 50.0, 1.0, 0.0     # relaxation time tR = C/K = 0.05 ~ dt scale
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.mass(2, 1.0)
    ops.uniaxialMaterial("Maxwell", 1, K, C, alpha, Lvar)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Constant", 1, "-factor", 1.0)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGen")
    ops.test("NormDispIncr", 1e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")

    # direct helper probe: one step, then a half-residual must be finite & positive (not garbage)
    for _ in range(3):
        ops.analyze(1, 0.02)
    u0, v0, a0 = ops.nodeDisp(2, 1), ops.nodeVel(2, 1), ops.nodeAccel(2, 1)
    ops.analyze(1, 0.02)
    pre = {1: ([0.0], [0.0], [0.0]), 2: ([u0], [v0], [a0])}
    post = {1: ([0.0], [0.0], [0.0]),
            2: ([ops.nodeDisp(2, 1)], [ops.nodeVel(2, 1)], [ops.nodeAccel(2, 1)])}
    r = _half_increment_residual(ops, pre, post, 0.02)
    assert r == r and r != float("inf") and 0.0 <= r < 1e6, r   # finite, not NaN/inf/garbage

    # full gate run on the rate-dependent model: completes with a bounded response
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.node(2, 1.0); ops.fix(1, 1); ops.mass(2, 1.0)
    ops.uniaxialMaterial("Maxwell", 1, K, C, alpha, Lvar)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Constant", 1, "-factor", 1.0)
    ops.pattern("Plain", 1, 1); ops.load(2, 1.0)
    ops.constraints("Plain"); ops.numberer("Plain"); ops.system("BandGen")
    ops.test("NormDispIncr", 1e-10, 50, 0); ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25); ops.analysis("Transient")
    res = robust_transient(ops, 0.4, 5.0e-3, error_gate=True, haftol=0.05)
    u = ops.nodeDisp(2, 1)
    assert bool(res)
    assert u == u and abs(u) < 1.0, u                # finite, bounded (static ~ P/K = 1e-3)
    assert res.halfres_max == res.halfres_max and res.halfres_max < 1e6


def test_gate_off_matches_w1i1a():
    """error_gate=False is the unchanged W1-I1a lane (completes, verdict integrated)."""
    _build_sdof(load_kind="constant", p0=1.0)
    res = robust_transient(ops, 0.2, 5.0e-3)
    assert bool(res)
    assert res.verdict == "integrated"
    assert res.halfres_max == 0.0 and res.n_overtol == 0
