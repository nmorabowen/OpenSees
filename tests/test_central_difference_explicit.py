"""CentralDifferenceExplicit: closed forms and dynamics checks.

Scheme on an undamped oscillator, with Omega = omega dt:

    u_{n+1} - 2 u_n + u_{n-1} = -Omega^2 u_n,   u_1 = u_0 + dt v_0 + dt^2/2 a_0.

Substituting u_n = exp(i n theta) gives cos(theta) = 1 - Omega^2/2, real for
Omega <= 2. With sin(theta) != 0 the exact discrete solution is

    u_n = u_0 cos(n theta) + (dt v_0 / sin(theta)) sin(n theta),

and the discrete period is T_d = 2 pi dt / theta = 2 pi dt / (2 asin(Omega/2)).

With viscous damping lagged at v_{n+1/2}, c = 2 xi omega m and gamma = 2 xi Omega,
the characteristic polynomial is z^2 - (2 - gamma - Omega^2) z + (1 - gamma). The
Jury conditions reduce to Omega^2 + 4 xi Omega - 4 < 0, that is

    dt_cr = (2/omega) (sqrt(1 + xi^2) - xi).
"""
import math

import pytest

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

CDE = "CentralDifferenceExplicit"
TWO_PI = 2.0 * math.pi


def _rms(xs):
    return math.sqrt(sum(v * v for v in xs) / len(xs))


def _linfit_slope(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    return sxy / sxx


def _all_finite(xs):
    return all(math.isfinite(v) for v in xs)


def _max_abs(xs):
    return max(abs(v) for v in xs)


def _sdof(k, m, integrator_args, dt, nsteps, d0=0.0, v0=0.0, alphaM=None, system=("Diagonal",)):
    """Free vibration on a zeroLength spring with nodal mass. Returns t, u, v."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0)
    ops.fix(2, 0, 1)
    ops.mass(2, m, 0.0)
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    if d0 != 0.0:
        ops.setNodeDisp(2, 1, d0, "-commit")
    if v0 != 0.0:
        ops.setNodeVel(2, 1, v0, "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system(*system)
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    if alphaM is not None:
        ops.rayleigh(alphaM, 0.0, 0.0, 0.0)
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    t, u, v = [ops.getTime()], [ops.nodeDisp(2, 1)], [ops.nodeVel(2, 1)]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        t.append(ops.getTime())
        u.append(ops.nodeDisp(2, 1))
        v.append(ops.nodeVel(2, 1))
    return t, u, v


def _build_bar(N, L, E, A, rho, fixed=True):
    """Chain of N lumped-mass trusses along x; returns the element length."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    le = L / N
    for i in range(N + 1):
        ops.node(i + 1, i * le, 0.0)
    ops.fix(1, 1 if fixed else 0, 1)
    for i in range(1, N + 1):
        ops.fix(i + 1, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    for i in range(N):
        ops.element("Truss", i + 1, i + 1, i + 2, A, 1, "-rho", rho)
    return le


def _analysis(*integrator_args):
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")


def test_exact_discrete_solution():
    """The whole trajectory matches the closed-form discrete solution to round-off."""
    m, w = 2.0, 3.0
    k = m * w * w
    u0, v0 = 0.8, -1.1
    Omega = 0.7
    dt = Omega / w
    _, u, _ = _sdof(k, m, (CDE,), dt, 500, d0=u0, v0=v0)
    assert len(u) == 501
    theta = math.acos(1.0 - 0.5 * Omega * Omega)
    B = dt * v0 / math.sin(theta)
    err = max(
        abs(ui - (u0 * math.cos(n * theta) + B * math.sin(n * theta)))
        for n, ui in enumerate(u)
    )
    assert err < 1e-11, "max error against the closed form %.3e" % err


def test_discrete_period():
    """Zero crossings give T_d = 2 pi dt / (2 asin(Omega/2)), shorter than 2 pi/omega."""
    m, w = 1.0, TWO_PI
    k = m * w * w
    dt = 0.05
    t, u, _ = _sdof(k, m, (CDE,), dt, 2000, d0=1.0)
    crossings = []
    for i in range(1, len(u)):
        if u[i - 1] < 0.0 <= u[i]:
            crossings.append(t[i - 1] + dt * (-u[i - 1]) / (u[i] - u[i - 1]))
    measured = (crossings[-1] - crossings[0]) / (len(crossings) - 1)
    exact = TWO_PI * dt / (2.0 * math.asin(0.5 * w * dt))
    assert abs(measured - exact) / exact < 1e-4, "T_d %.6f vs %.6f" % (measured, exact)
    assert measured < TWO_PI / w


def test_second_order_convergence():
    m, w = 1.0, TWO_PI
    k = m * w * w
    v0 = 1.0
    log_dts, log_errs = [], []
    for N in (40, 80, 160, 320, 640):
        dt = 1.0 / N
        t, u, _ = _sdof(k, m, (CDE,), dt, N, v0=v0)
        err = [ui - (v0 / w) * math.sin(w * ti) for ti, ui in zip(t, u)]
        log_dts.append(math.log(dt))
        log_errs.append(math.log(_rms(err)))
    slope = _linfit_slope(log_dts, log_errs)
    assert 1.7 <= slope <= 2.3, "convergence slope %.3f" % slope


def _stable_at(Omega, m=1.0, w=TWO_PI):
    k = m * w * w
    _, u, _ = _sdof(k, m, (CDE,), Omega / w, 200, d0=1.0)
    if not _all_finite(u) or len(u) < 201:
        return False
    return _max_abs(u[-40:]) <= 5.0 * max(_max_abs(u[:40]), 1e-12)


def test_undamped_stability_limit():
    assert _stable_at(1.99)
    assert not _stable_at(2.01)


@pytest.mark.parametrize("xi", [0.02, 0.1, 0.5])
def test_damped_stability_limit(xi):
    m, w = 1.0, TWO_PI
    k = m * w * w
    alphaM = 2.0 * xi * w
    dt_cr = (2.0 / w) * (math.sqrt(1.0 + xi * xi) - xi)
    _, u_lo, _ = _sdof(k, m, (CDE,), 0.98 * dt_cr, 400, d0=1.0, alphaM=alphaM)
    _, u_hi, _ = _sdof(k, m, (CDE,), 1.02 * dt_cr, 400, d0=1.0, alphaM=alphaM)
    assert len(u_lo) == 401 and _all_finite(u_lo)
    assert _max_abs(u_lo[-20:]) <= _max_abs(u_lo[:20])
    assert (
        len(u_hi) < 401
        or not _all_finite(u_hi)
        or _max_abs(u_hi[-20:]) > 5.0 * _max_abs(u_hi[:20])
    )


def test_damped_free_vibration_decay():
    """Amplitude ratio per step is |z| = sqrt(1 - gamma), gamma = 2 xi Omega."""
    m, w, xi = 1.0, TWO_PI, 0.05
    k = m * w * w
    Omega = 0.2
    dt = Omega / w
    n = 2000
    gamma = 2.0 * xi * Omega
    _, u, _ = _sdof(k, m, (CDE,), dt, n, d0=1.0, alphaM=2.0 * xi * w)
    # reference: the two-term recurrence of the characteristic polynomial
    b = 2.0 - gamma - Omega * Omega
    c = 1.0 - gamma
    # u_1 from the starter: a_0 = -(c0 v0 + k u0)/m with v0 = 0
    u_ref = [1.0, 1.0 - 0.5 * Omega * Omega]
    for _ in range(n - 1):
        u_ref.append(b * u_ref[-1] - c * u_ref[-2])
    err = max(abs(a - r) for a, r in zip(u, u_ref))
    assert err < 1e-10, "max error against the recurrence %.3e" % err
    # envelope decays with rho = sqrt(1 - gamma)
    env_ratio = _max_abs(u[1500:1600]) / _max_abs(u[500:600])
    assert env_ratio == pytest.approx(math.sqrt(1.0 - gamma) ** 1000, rel=0.1)


def test_crosscheck_newmark_and_central_difference():
    m, w = 1.0, TWO_PI
    k = m * w * w
    dt = 1.0 / 200
    _, ucde, _ = _sdof(k, m, (CDE,), dt, 400, d0=1.0)
    _, unm, _ = _sdof(k, m, ("Newmark", 0.5, 0.25), dt, 400, d0=1.0)
    _, ucd, _ = _sdof(k, m, ("CentralDifference",), dt, 400, d0=1.0)
    assert max(abs(a - b) for a, b in zip(ucde, unm)) < 0.03
    assert max(abs(a - b) for a, b in zip(ucde, ucd)) < 0.03


def test_first_step_with_initial_velocity_and_acceleration():
    """u_1 = u_0 + dt v_0 + dt^2/2 a_0 exactly; first steps beat the older class."""
    m, w = 1.0, TWO_PI
    k = m * w * w
    u0, v0 = 1.0, 0.7
    dt = 1.0 / 200
    t, u, v = _sdof(k, m, (CDE,), dt, 5, d0=u0, v0=v0)
    a0 = -w * w * u0
    assert u[1] == pytest.approx(u0 + dt * v0 + 0.5 * dt * dt * a0, abs=1e-14)
    ue = [u0 * math.cos(w * ti) + (v0 / w) * math.sin(w * ti) for ti in t]
    err_cde = max(abs(a - b) for a, b in zip(u, ue))
    _, ucd, _ = _sdof(k, m, ("CentralDifference",), dt, 5, d0=u0, v0=v0)
    err_cd = max(abs(a - b) for a, b in zip(ucd, ue))
    assert err_cde < 1e-3
    assert err_cde <= err_cd


def test_full_step_velocity_output():
    """Nodal velocity is the centered difference (u_{n+1} - u_{n-1}) / (2 dt)."""
    m, w = 1.0, TWO_PI
    k = m * w * w
    dt = 0.01
    _, u, v = _sdof(k, m, (CDE,), dt, 100, d0=1.0, v0=0.3)
    for n in range(1, 100):
        assert v[n] == pytest.approx((u[n + 1] - u[n - 1]) / (2.0 * dt), abs=1e-10)


def test_wave_speed():
    N, L, E, A, rho, P = 100, 20.0, 100.0, 1.0, 1.0, 1.0
    c = math.sqrt(E / rho)
    le = _build_bar(N, L, E, A, rho)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(N + 1, P, 0.0)
    _analysis(CDE)
    dt = 0.5 * le / c
    vth = 0.5 * P / (A * rho * c)
    dists = [5.0, 8.0, 11.0, 14.0]
    nodes = [int(round((L - d) / le)) + 1 for d in dists]
    arrival = {nd: None for nd in nodes}
    tcur = 0.0
    for _ in range(int(2.2 * L / c / dt)):
        ops.analyze(1, dt)
        tcur += dt
        for nd in nodes:
            if arrival[nd] is None and abs(ops.nodeVel(nd, 1)) > vth:
                arrival[nd] = tcur
    pts = [(d, arrival[nd]) for d, nd in zip(dists, nodes) if arrival[nd]]
    assert len(pts) >= 2, "wave front not detected at two stations"
    c_meas = 1.0 / _linfit_slope([p[0] for p in pts], [p[1] for p in pts])
    assert abs(c_meas - c) / c < 0.05, "wave speed %.3f vs %.3f" % (c_meas, c)


def test_rigid_body_momentum():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 0, 1)
    m, v0, dt, N = 2.0, 3.0, 0.01, 100
    ops.mass(1, m, 0.0)
    ops.setNodeVel(1, 1, v0, "-commit")
    _analysis(CDE)
    for _ in range(N):
        ops.analyze(1, dt)
    assert ops.nodeVel(1, 1) == pytest.approx(v0, abs=1e-9)
    assert ops.nodeDisp(1, 1) == pytest.approx(v0 * N * dt, abs=1e-6)


def test_mechanical_energy_closure():
    """Undamped distributed-mass bar: KE and IE exchange, KE + IE drifts below 2%."""
    N, L, E, A, rho, v0 = 20, 10.0, 100.0, 1.0, 1.0, 2.0
    le = _build_bar(N, L, E, A, rho)
    for i in range(1, N + 1):
        ops.setNodeVel(i + 1, 1, v0, "-commit")
    _analysis(CDE)
    dt = 0.5 * le / math.sqrt(E / rho)
    kl = E * A / le
    node_mass = [0.5 * rho * A * le] + [rho * A * le] * (N - 1) + [0.5 * rho * A * le]
    ke, ie = [], []
    for _ in range(400):
        assert ops.analyze(1, dt) == 0
        vel = [ops.nodeVel(i + 1, 1) for i in range(N + 1)]
        dsp = [ops.nodeDisp(i + 1, 1) for i in range(N + 1)]
        ke.append(0.5 * sum(mi * vi * vi for mi, vi in zip(node_mass, vel)))
        ie.append(0.5 * kl * sum((dsp[i + 1] - dsp[i]) ** 2 for i in range(N)))
    mech = [a + b for a, b in zip(ke, ie)]
    assert max(ie) > 0.05 * ke[0]
    assert (max(mech) - min(mech)) / (sum(mech) / len(mech)) < 0.02


@pytest.mark.parametrize("flag", [(), ("-noMassCache",)])
def test_mass_cache_equivalence(flag):
    """Reused and reassembled mass tangents give the same trajectory, bit for bit."""
    m, w = 1.0, TWO_PI
    k = m * w * w
    _, a, _ = _sdof(k, m, (CDE,), 0.01, 300, d0=1.0, v0=0.2)
    _, b, _ = _sdof(k, m, (CDE,) + flag, 0.01, 300, d0=1.0, v0=0.2)
    assert a == b


def test_mass_cache_with_a_factoring_solver():
    """A reused tangent stays valid under a solver that factors in place."""
    m, w = 1.0, TWO_PI
    k = m * w * w
    _, a, _ = _sdof(k, m, (CDE,), 0.01, 300, d0=1.0, v0=0.2)
    _, b, _ = _sdof(k, m, (CDE,), 0.01, 300, d0=1.0, v0=0.2, system=("ProfileSPD",))
    assert max(abs(x - y) for x, y in zip(a, b)) < 1e-12


def test_retry_after_divergence_abort():
    """A failed step restores the committed state; the retry is a clean restart."""
    m, w = 1.0, TWO_PI
    k = m * w * w
    Omega = 0.2
    dt = Omega / w

    def run(with_failure):
        _sdof(k, m, (CDE, "-cflAbort"), dt, 0, d0=1.0)
        for _ in range(50):
            assert ops.analyze(1, dt) == 0
        if with_failure:
            assert ops.analyze(1, 3.0 / w) != 0    # Omega = 3 > 2: refused by -cflAbort
        u50, v50 = ops.nodeDisp(2, 1), ops.nodeVel(2, 1)
        return u50, v50, ops.getTime()

    u_ok, v_ok, t_ok = run(False)
    u_f, v_f, t_f = run(True)
    assert (u_f, v_f, t_f) == (u_ok, v_ok, t_ok)
    # the restart centers v_{-1/2} from the committed (u, v, a), so the next steps
    # keep the accuracy of the scheme
    for _ in range(200):
        assert ops.analyze(1, dt) == 0
    t = ops.getTime()
    assert abs(ops.nodeDisp(2, 1) - math.cos(w * t)) < 0.05
