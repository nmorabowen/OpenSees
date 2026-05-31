"""CentralDifferenceLadruno (Ladruno explicit leap-frog integrator, classTag
33003) — Zone-A dynamics battery (CDL-1..10).

Zone-A port of the CDL portion of ``Ladruno_scripts/_verify_explicit.py``. Each
``_verify_explicit`` CDL check becomes one pytest. The script depended on numpy;
Zone-A CI installs only ``conan pytest pyyaml`` (no numpy), so the linear-algebra
bits (RMS, log-log slope, peak finding) are reimplemented in pure ``math`` here.

The integrator is the explicit leap-frog central difference "done right":
  * correct first-step starter: a0 = M^-1 (P0 - C v0 - Fint(u0)),
    v_{-1/2} = v0 - dt/2 a0 on the first newStep;
  * built-in critical time step (CriticalTimeStep, factor 1.0);
  * clean full-step velocity output;
  * data-driven betaK guard (Rayleigh stiffness damping shrinks dt_cr).

CDL-8 (first-step correctness with nonzero v0/a0) and CDL-9 (betaK collapses
dt_cr ~quadratically while alphaM stays benign) are the differentiators versus
the legacy CentralDifference class. CDL-10 (energy closure via the
EnergyBalanceRecorder) was link-gated until the CriticalTimeStep dsygvx_ link
fix (PR #28) restored the full Linux build — it now runs.

Theory: Ladruno_implementation/05_robust_central_difference.md.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"
TWO_PI = 2.0 * math.pi


# --------------------------------------------------------------------------
# pure-python numeric helpers (numpy-free)
# --------------------------------------------------------------------------
def _rms(xs):
    return math.sqrt(sum(v * v for v in xs) / len(xs))


def _linfit_slope(xs, ys):
    """Least-squares slope of ys vs xs (both 1-D, equal length)."""
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


# --------------------------------------------------------------------------
# shared models (ported from _verify_explicit.py, numpy removed)
# --------------------------------------------------------------------------
def _sdof(k, m, integrator_args, dt, nsteps, d0=0.0, v0=0.0):
    """Free-vibration SDOF; returns (t[], u[]). integrator_args = tuple."""
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
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    t = [ops.getTime()]
    u = [ops.nodeDisp(2, 1)]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        t.append(ops.getTime())
        u.append(ops.nodeDisp(2, 1))
    return t, u


def _stable_at(integrator_args, omega_dt, m=1.0, w=TWO_PI):
    """Run an SDOF at Omega = w*dt; True if it does not blow up."""
    k = m * w * w
    dt = omega_dt / w
    _, u = _sdof(k, m, integrator_args, dt, 200, d0=1.0)
    if not _all_finite(u):
        return False
    early = _max_abs(u[:40])
    late = _max_abs(u[-40:])
    return late <= 5.0 * max(early, 1e-12)


def _build_bar(N, L, E, A, rho):
    """Fixed-free chain of N truss elements along x; returns element length."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    le = L / N
    for i in range(N + 1):
        ops.node(i + 1, i * le, 0.0)
    ops.fix(1, 1, 1)
    for i in range(1, N + 1):
        ops.fix(i + 1, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    for i in range(N):
        ops.element("Truss", i + 1, i + 1, i + 2, A, 1, "-rho", rho)
    return le


# ==========================================================================
# CDL-1  order of accuracy ~ 2.0
# ==========================================================================
def test_cdl_order():
    m = 1.0
    w = TWO_PI
    k = m * w * w
    v0 = 1.0
    Tend = 1.0
    log_dts, log_errs = [], []
    for N in (40, 80, 160, 320, 640):
        dt = Tend / N
        t, u = _sdof(k, m, (CDL,), dt, N, v0=v0)
        err = [ui - (v0 / w) * math.sin(w * ti) for ti, ui in zip(t, u)]
        log_dts.append(math.log(dt))
        log_errs.append(math.log(_rms(err)))
    slope = _linfit_slope(log_dts, log_errs)
    assert 1.7 <= slope <= 2.3, "convergence slope = %.3f (expect ~2.0)" % slope


# ==========================================================================
# CDL-2  undamped stability boundary: stable for w*dt < 2, unstable just above
# ==========================================================================
def test_cdl_stability():
    assert _stable_at((CDL,), 1.9), "should be stable at Omega=1.9 (< 2.0)"
    assert not _stable_at((CDL,), 2.1), "should be unstable at Omega=2.1 (> 2.0)"


# ==========================================================================
# CDL-3  damped stability limit (2/w)(sqrt(1+xi^2)-xi), xi = 0.02/0.1/0.5
# ==========================================================================
def _sdof_damped(k, m, alphaM, dt, nsteps, d0=1.0):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0)
    ops.fix(2, 0, 1)
    ops.mass(2, m, 0.0)
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.setNodeDisp(2, 1, d0, "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.rayleigh(alphaM, 0.0, 0.0, 0.0)
    ops.integrator(CDL)
    ops.analysis("Transient")
    u = [ops.nodeDisp(2, 1)]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return None
        u.append(ops.nodeDisp(2, 1))
    return u


@pytest.mark.parametrize("xi", [0.02, 0.1, 0.5])
def test_cdl_damped(xi):
    m = 1.0
    w = TWO_PI
    k = m * w * w
    alphaM = 2.0 * xi * w
    dt_cr = (2.0 / w) * (math.sqrt(1.0 + xi * xi) - xi)  # damped CD limit (ADR T3)
    u_lo = _sdof_damped(k, m, alphaM, 0.9 * dt_cr, 200)  # under -> stable
    u_hi = _sdof_damped(k, m, alphaM, 1.1 * dt_cr, 200)  # over  -> blows up
    stable_lo = (
        u_lo is not None
        and _all_finite(u_lo)
        and _max_abs(u_lo[-20:]) <= _max_abs(u_lo[:20])
    )
    unstable_hi = (
        u_hi is None
        or not _all_finite(u_hi)
        or _max_abs(u_hi[-20:]) > 5.0 * _max_abs(u_hi[:20])
    )
    assert stable_lo, "xi=%.2f: should be stable below dt_cr=%.4f" % (xi, dt_cr)
    assert unstable_hi, "xi=%.2f: should be unstable above dt_cr=%.4f" % (xi, dt_cr)


# ==========================================================================
# CDL-4  cross-check vs Newmark (avg accel) and CentralDifference at small dt
# ==========================================================================
def test_cdl_crosscheck():
    m = 1.0
    w = TWO_PI
    k = m * w * w
    d0 = 1.0
    dt = (TWO_PI / w) / 200
    N = 400
    _, ucdl = _sdof(k, m, (CDL,), dt, N, d0=d0)
    _, unm = _sdof(k, m, ("Newmark", 0.5, 0.25), dt, N, d0=d0)
    _, ucd = _sdof(k, m, ("CentralDifference",), dt, N, d0=d0)
    e_nm = max(abs(a - b) for a, b in zip(ucdl, unm)) / d0
    e_cd = max(abs(a - b) for a, b in zip(ucdl, ucd)) / d0
    assert e_nm < 0.03, "max|CDL-Newmark| = %.3f%% (expect < 3%%)" % (100 * e_nm)
    assert e_cd < 0.03, "max|CDL-CD| = %.3f%% (expect < 3%%)" % (100 * e_cd)


# ==========================================================================
# CDL-5  1-D wave speed exact = sqrt(E/rho)
# ==========================================================================
def test_cdl_wave():
    N = 100
    L = 20.0
    E = 100.0
    A = 1.0
    rho = 1.0
    P = 1.0
    c = math.sqrt(E / rho)
    le = _build_bar(N, L, E, A, rho)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(N + 1, P, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")
    dt = 0.5 * le / c
    v_particle = P / (A * rho * c)
    vth = 0.5 * v_particle
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
    assert len(pts) >= 2, "wave front not detected at >= 2 stations"
    slope = _linfit_slope([p[0] for p in pts], [p[1] for p in pts])  # = 1/speed
    c_meas = 1.0 / slope
    err = abs(c_meas - c) / c
    assert err < 0.05, "measured wave speed %.3f vs sqrt(E/rho)=%.3f (err %.1f%%)" % (
        c_meas,
        c,
        100 * err,
    )


# ==========================================================================
# CDL-6  rigid-body momentum conserved exactly (zero stiffness, nonzero v0)
# ==========================================================================
def test_cdl_rigid_body():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 0, 1)
    m = 2.0
    v0 = 3.0
    ops.mass(1, m, 0.0)
    ops.setNodeVel(1, 1, v0, "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")
    dt = 0.01
    N = 100
    for _ in range(N):
        ops.analyze(1, dt)
    u = ops.nodeDisp(1, 1)
    v = ops.nodeVel(1, 1)
    u_exact = v0 * (N * dt)  # free particle: u = v0 t, v = v0
    assert abs(v - v0) < 1e-9, "velocity drift: v=%.9f vs v0=%.1f" % (v, v0)
    assert abs(u - u_exact) < 1e-6, "position: u=%.6f vs v0*t=%.6f" % (u, u_exact)


# ==========================================================================
# CDL-7  criticalTimeStep() = le/c exact on a 1-element bar (2/omega_max)
# ==========================================================================
def test_cdl_dtcr_bar():
    N = 1
    L = 1.0
    E = 100.0
    A = 1.0
    rho = 1.0
    le = _build_bar(N, L, E, A, rho)
    c = math.sqrt(E / rho)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl")
    ops.analysis("Transient")
    dtcr = ops.criticalTimeStep()
    target = le / c  # single free DOF: omega_max = 2c/le -> dt_cr = le/c
    err = abs(dtcr - target) / target
    assert err < 0.02, "criticalTimeStep()=%.5f vs le/c=%.5f (err %.1f%%)" % (
        dtcr,
        target,
        100 * err,
    )


# ==========================================================================
# CDL-8  first-step correctness (THE differentiator): nonzero v0 AND nonzero
#        a0 (= -w^2 u0). The starter reproduces u(t)=u0 cos wt + (v0/w) sin wt
#        from step 1, where the legacy CentralDifference carries a 1st-step error.
# ==========================================================================
def test_cdl_first_step():
    m = 1.0
    w = TWO_PI
    k = m * w * w
    u0 = 1.0
    v0 = 0.7  # both nonzero -> a0 = -w^2 u0 != 0
    dt = (TWO_PI / w) / 200
    N = 5  # look at the FIRST few steps only
    t, u = _sdof(k, m, (CDL,), dt, N, d0=u0, v0=v0)
    ue = [u0 * math.cos(w * ti) + (v0 / w) * math.sin(w * ti) for ti in t]
    err_cdl = max(abs(a - b) for a, b in zip(u, ue))
    # legacy class on the same problem: expected larger first-step error
    _, ucd = _sdof(k, m, ("CentralDifference",), dt, N, d0=u0, v0=v0)
    err_cd = max(abs(a - b) for a, b in zip(ucd, ue))
    assert err_cdl < 1e-3, "max|CDL-exact| = %.2e (expect < 1e-3)" % err_cdl
    assert err_cdl <= err_cd, (
        "CDL first-step error %.2e should be <= legacy CD's %.2e" % (err_cdl, err_cd)
    )


# ==========================================================================
# CDL-9  betaK trap + guard (T4): on a stiff bar dt_cr collapses ~quadratically
#        as Rayleigh stiffness damping betaK rises, while the same nominal
#        damping applied as mass-proportional alphaM stays benign.
# ==========================================================================
def _bar_dtcr_cdl(extra_rayleigh):
    N = 20
    L = 1.0
    E = 1.0e6
    A = 1.0
    rho = 1.0
    _build_bar(N, L, E, A, rho)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    if extra_rayleigh is not None:
        ops.rayleigh(*extra_rayleigh)
    ops.integrator(CDL, "-cfl")
    ops.analysis("Transient")
    return ops.criticalTimeStep()


def test_cdl_betak_trap():
    dt0 = _bar_dtcr_cdl(None)                      # undamped reference
    dt_bk1 = _bar_dtcr_cdl((0.0, 1e-5, 0.0, 0.0))  # small betaK
    dt_bk2 = _bar_dtcr_cdl((0.0, 4e-5, 0.0, 0.0))  # 4x betaK -> ~1/4 the step
    dt_am = _bar_dtcr_cdl((50.0, 0.0, 0.0, 0.0))   # alphaM -> benign
    # betaK collapses dt_cr ~quadratically (4x betaK -> ~4x smaller dt_cr);
    # alphaM leaves it essentially unchanged.
    assert dt_bk1 < 0.95 * dt0, "betaK should shrink dt_cr (%.3e vs %.3e)" % (dt_bk1, dt0)
    assert dt_bk2 < 0.6 * dt_bk1, (
        "4x betaK should ~quarter dt_cr (%.3e vs %.3e)" % (dt_bk2, dt_bk1)
    )
    assert dt_am > 0.9 * dt0, "alphaM should leave dt_cr ~unchanged (%.3e vs %.3e)" % (
        dt_am,
        dt0,
    )


# ==========================================================================
# CDL-10  energy closure via EnergyBalanceRecorder (multi-DOF). This was the
#         link-gated skip in _verify_explicit.py; the CriticalTimeStep dsygvx_
#         link fix (PR #28) restored the full build, so it now runs. Undamped
#         free vibration of a distributed-mass bar: KE <-> IE exchange with
#         total mechanical energy (KE+IE) conserved to ~1%.
# ==========================================================================
def test_cdl_energy_closure(tmp_path):
    N = 20
    L = 10.0
    E = 100.0
    A = 1.0
    rho = 1.0
    v0 = 2.0
    efile = str(tmp_path / "cdl_energy.txt")
    le = _build_bar(N, L, E, A, rho)
    for i in range(1, N + 1):  # initial x-velocity on every free node
        ops.setNodeVel(i + 1, 1, v0, "-commit")
    # columns (with -time): time, KE, IE, DW, ULW, RES, ERR%
    ops.recorder("EnergyBalance", "-file", efile, "-time")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")
    c = math.sqrt(E / rho)
    dt = 0.5 * le / c
    for _ in range(400):
        ops.analyze(1, dt)
    ops.wipe()  # flush/close the recorder

    rows = []
    with open(efile) as fh:
        for line in fh:
            vals = line.split()
            if vals:
                rows.append([float(x) for x in vals])
    assert len(rows) > 10, "recorder produced too few rows (%d)" % len(rows)

    ke = [r[1] for r in rows]
    ie = [r[2] for r in rows]
    mech = [a + b for a, b in zip(ke, ie)]  # total mechanical energy KE + IE
    assert ke[0] > 0.0, "initial KE not recorded (>0 expected with v0 on every node)"
    # genuine exchange must occur (otherwise IE column is dead): some IE built up
    assert max(ie) > 0.05 * ke[0], "no KE->IE exchange recorded (max IE = %.4g)" % max(ie)
    e_mean = sum(mech) / len(mech)
    drift = (max(mech) - min(mech)) / e_mean
    assert drift < 0.02, "mechanical energy drift = %.3f%% (expect < 2%%)" % (100 * drift)
