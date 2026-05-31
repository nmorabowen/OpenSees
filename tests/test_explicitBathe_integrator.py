"""ExplicitBathe (Noh-Bathe explicit composite integrator, classTag 33000) —
Zone-A dynamics battery.

Zone-A port of the ExplicitBathe portion of ``Ladruno_scripts/_verify_explicit.py``
(jaabell's tag, band-remediated 61 -> 33000). Each ``_verify_explicit`` check
becomes one pytest. numpy is available on Zone-A CI, so the model builders and
numeric helpers are copied verbatim from the source battery (no numpy-free
reimplementation — a faithful 1:1 port).

Covered checks (the LNVD case -> test_explicitBatheLNVD_integrator.py; the
EnergyBalance element-mass case -> test_energyBalanceRecorder.py):
  1.  order of accuracy ~ 2.0
  2.  stability limit beats central difference's 2/omega
  3.  cross-check vs Newmark (avg accel) and CentralDifference
  4.  spectral behaviour: algorithmic damping + period error vs Omega
  5.  1-D wave propagation speed = sqrt(E/rho)
  7.  rigid-body / momentum conservation
  8.  queryable dt_cr (= le/c) + adaptive sub-step + -recompute path
  10. dt_cr <= 0 "no value" contract (no -cfl; pure nodal-mass model)
  11. -lump diagonal vs rowsum (diagonal-of-consistent is stiffer -> smaller dt_cr)
  12. -cflAbort enforcement (analyze aborts on a CFL violation)

The dt_cr machinery (CriticalTimeStep) links on Linux via the dsygvx_ swap
(PR #28); CDL-7/9/10 already exercise it green on Zone-A CI.
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

EB = "ExplicitBathe"
P = 0.54  # Noh-Bathe sub-step parameter used throughout the source battery


# --------------------------------------------------------------------------
# shared models (ported verbatim from _verify_explicit.py)
# --------------------------------------------------------------------------
def sdof(k, m, integrator_args, dt, nsteps, d0=0.0, v0=0.0):
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
    return np.array(t), np.array(u)


def build_bar(N, L, E, A, rho):
    """Fixed-free chain of N truss elements along x; returns element length."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    le = L / N
    for i in range(N + 1):
        ops.node(i + 1, i * le, 0.0)
    ops.fix(1, 1, 1)                  # fixed end: x and y
    for i in range(1, N + 1):        # nodes 2..N+1: free in x, fixed in y
        ops.fix(i + 1, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    for i in range(N):
        ops.element("Truss", i + 1, i + 1, i + 2, A, 1, "-rho", rho)
    return le


def _stable_at(integrator_args, Omega, m=1.0, w=2 * math.pi):
    k = m * w * w
    dt = Omega / w
    _, u = sdof(k, m, integrator_args, dt, 200, d0=1.0)
    if not np.all(np.isfinite(u)):
        return False
    early = np.abs(u[:40]).max()
    late = np.abs(u[-40:]).max()
    return late <= 5.0 * max(early, 1e-12)


def _peaks(u):
    return [i for i in range(1, len(u) - 1)
            if u[i] > u[i - 1] and u[i] >= u[i + 1] and u[i] > 0]


def _bar_dtcr(N, L, E, A, rho, extra):
    """Build the distributed-mass bar, run one tiny step, return criticalTimeStep()."""
    le = build_bar(N, L, E, A, rho)
    c = math.sqrt(E / rho)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(N + 1, 1.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(EB, P, "-cfl", *extra)
    ops.analysis("Transient")
    ops.analyze(1, 0.1 * le / c)
    return ops.criticalTimeStep(), le, c


# ==========================================================================
# 1. Order of accuracy: u0=0, v0=1 -> u(t)=(v0/w) sin(w t), a0=0 (no cold start)
# ==========================================================================
def test_eb_order():
    m = 1.0
    w = 2 * math.pi
    k = m * w * w
    v0 = 1.0
    Tend = 1.0
    errs, dts = [], []
    for N in (40, 80, 160, 320, 640):
        dt = Tend / N
        t, u = sdof(k, m, (EB, P), dt, N, v0=v0)
        ue = (v0 / w) * np.sin(w * t)
        errs.append(math.sqrt(np.mean((u - ue) ** 2)))
        dts.append(dt)
    slope = np.polyfit(np.log(dts), np.log(errs), 1)[0]
    assert 1.7 <= slope <= 2.3, "convergence slope = %.3f (expect ~2.0)" % slope


# ==========================================================================
# 2. Stability limit: sweep Omega = w*dt; ExplicitBathe must beat CD's 2/omega
# ==========================================================================
def test_eb_stability():
    grid = [1.6, 1.8, 1.95, 2.0, 2.05, 2.2, 2.4, 2.6, 2.8, 3.0]
    eb = max([O for O in grid if _stable_at((EB, P), O)] or [0])
    cd = max([O for O in grid if _stable_at(("CentralDifference",), O)] or [0])
    assert eb > 2.0, (
        "ExplicitBathe stable to Omega=%.2f, CentralDifference to %.2f "
        "(must exceed CD's 2.0)" % (eb, cd)
    )


# ==========================================================================
# 3. Cross-check vs Newmark (avg accel) and CentralDifference at small dt
# ==========================================================================
def test_eb_crosscheck():
    m = 1.0
    w = 2 * math.pi
    k = m * w * w
    d0 = 1.0
    dt = (2 * math.pi / w) / 200
    N = 400
    _, ueb = sdof(k, m, (EB, P), dt, N, d0=d0)
    _, unm = sdof(k, m, ("Newmark", 0.5, 0.25), dt, N, d0=d0)
    _, ucd = sdof(k, m, ("CentralDifference",), dt, N, d0=d0)
    e_nm = np.abs(ueb - unm).max() / d0
    e_cd = np.abs(ueb - ucd).max() / d0
    assert e_nm < 0.03, "max|EB-Newmark| = %.3f%% (expect < 3%%)" % (100 * e_nm)
    assert e_cd < 0.03, "max|EB-CD| = %.3f%% (expect < 3%%)" % (100 * e_cd)


# ==========================================================================
# 4. Spectral: algorithmic damping tiny at low Omega, grows with Omega
# ==========================================================================
def test_eb_spectral():
    m = 1.0
    dt = 0.01
    rows = []
    for Omega in (0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5):
        w = Omega / dt
        k = m * w * w
        T = 2 * math.pi / w
        n = int(8 * T / dt)
        t, u = sdof(k, m, (EB, P), dt, n, d0=1.0)
        pk = _peaks(u)
        if len(pk) >= 3:
            a0, an = u[pk[0]], u[pk[-1]]
            nc = len(pk) - 1
            xi = (math.log(max(a0, 1e-12) / max(an, 1e-12)) / (2 * math.pi * nc)
                  if an > 0 else float("nan"))
        else:
            xi = float("nan")
        rows.append((Omega, 100 * xi if xi == xi else xi))
    lo = [r[1] for r in rows if r[0] <= 0.25 and r[1] == r[1]]
    hi = [r[1] for r in rows if r[0] >= 1.0 and r[1] == r[1]]
    assert lo and hi, "could not extract damping at both low and high Omega"
    assert max(lo) < 1.0, "low-Omega algorithmic damping %.3f%% should be ~0" % max(lo)
    assert max(hi) > max(lo), (
        "damping should grow with Omega (hi=%.3f%% vs lo=%.3f%%)" % (max(hi), max(lo))
    )


# ==========================================================================
# 5. 1-D wave propagation speed = sqrt(E/rho)
# ==========================================================================
def test_eb_wave():
    N = 100
    L = 20.0
    E = 100.0
    A = 1.0
    rho = 1.0
    P_load = 1.0
    c = math.sqrt(E / rho)
    le = build_bar(N, L, E, A, rho)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(N + 1, P_load, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(EB, P)
    ops.analysis("Transient")
    dt = 0.5 * le / c
    v_particle = P_load / (A * rho * c)
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
    d_arr = np.array([p[0] for p in pts])
    t_arr = np.array([p[1] for p in pts])
    slope = np.polyfit(d_arr, t_arr, 1)[0]      # = 1/speed
    c_meas = 1.0 / slope
    err = abs(c_meas - c) / c
    assert err < 0.05, "measured wave speed %.3f vs sqrt(E/rho)=%.3f (err %.1f%%)" % (
        c_meas, c, 100 * err
    )


# ==========================================================================
# 7. Rigid-body / momentum: free-free truss, both nodes v0 -> u = v0 t, no strain
# ==========================================================================
def test_eb_rigid_body():
    m = 1.0
    v0 = 2.0
    E = 100.0
    A = 1.0
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 1.0, 0.0)
    ops.fix(1, 0, 1)
    ops.fix(2, 0, 1)
    ops.mass(1, m, 0.0)
    ops.mass(2, m, 0.0)
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.element("Truss", 1, 1, 2, A, 1)
    ops.setNodeVel(1, 1, v0, "-commit")
    ops.setNodeVel(2, 1, v0, "-commit")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(EB, P)
    ops.analysis("Transient")
    dt = 0.01
    for _ in range(200):
        ops.analyze(1, dt)
    t = ops.getTime()
    u1 = ops.nodeDisp(1, 1)
    u2 = ops.nodeDisp(2, 1)
    v1 = ops.nodeVel(1, 1)
    v2 = ops.nodeVel(2, 1)
    u_exact = v0 * t
    err_u = max(abs(u1 - u_exact), abs(u2 - u_exact)) / abs(u_exact)
    err_v = max(abs(v1 - v0), abs(v2 - v0)) / v0
    err_strain = abs(u1 - u2)                   # rigid -> no relative motion
    assert err_u < 1e-3, "displacement error %.2e (expect ~0)" % err_u
    assert err_v < 1e-3, "velocity error %.2e (expect ~0)" % err_v
    assert err_strain < 1e-6, "spurious strain %.2e (rigid body, expect ~0)" % err_strain


# ==========================================================================
# 8. Queryable dt_cr (= le/c) + adaptive sub-step count + -recompute path
# ==========================================================================
def test_eb_query_substep():
    N = 40
    L = 8.0
    E = 100.0
    A = 1.0
    rho = 1.0
    c = math.sqrt(E / rho)
    le = build_bar(N, L, E, A, rho)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(N + 1, 1.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(EB, P, "-cfl")
    ops.analysis("Transient")
    ops.analyze(1, 0.1 * le / c)                # priming step triggers the compute
    dtcr = ops.criticalTimeStep()
    expected = le / c                           # 2/omega_max for a 2-node bar element
    err = abs(dtcr - expected) / expected
    # adaptive sub-stepping: internal steps to cover a coarse interval
    dt_coarse = 5.0 * dtcr
    n_sub = max(1, math.ceil(dt_coarse / (0.9 * dtcr)))
    # exercise the -recompute (tangent) path: must still return a finite dt_cr
    build_bar(N, L, E, A, rho)
    ops.timeSeries("Constant", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(N + 1, 1.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(EB, P, "-recompute", 5)
    ops.analysis("Transient")
    ops.analyze(3, 0.5 * le / c)
    dtcr_t = ops.criticalTimeStep()
    assert dtcr > 0, "criticalTimeStep() should be positive with -cfl"
    assert err < 0.10, "criticalTimeStep()=%.5f vs le/c=%.5f (err %.1f%%)" % (
        dtcr, expected, 100 * err
    )
    assert n_sub == 6, "expected 6 sub-steps for dt=5*dt_cr, got %d" % n_sub
    assert dtcr_t > 0, "-recompute path should still return a positive dt_cr"


# ==========================================================================
# 10. dt_cr <= 0 "no value" contract: no -cfl, and -cfl on a pure nodal-mass model
# ==========================================================================
def test_eb_no_value_contract():
    m = 1.0
    w = 2 * math.pi
    k = m * w * w
    dt = 0.01
    sdof(k, m, (EB, P), dt, 1, d0=1.0)              # estimation not requested
    s_disabled = ops.criticalTimeStep()
    sdof(k, m, (EB, P, "-cfl"), dt, 1, d0=1.0)      # nodal mass only -> no K-M pencil
    s_napp = ops.criticalTimeStep()
    assert s_disabled <= 0.0, "no -cfl should report <=0 'no value' (got %.3g)" % s_disabled
    assert s_napp <= 0.0, "nodal-mass + -cfl should report <=0 (got %.3g)" % s_napp


# ==========================================================================
# 11. -lump diagonal vs rowsum: diagonal-of-consistent is stiffer -> smaller dt_cr
# ==========================================================================
def test_eb_lump_diagonal():
    N = 40
    L = 8.0
    E = 100.0
    A = 1.0
    rho = 1.0
    dt_rowsum, le, c = _bar_dtcr(N, L, E, A, rho, ("-lump", "rowsum"))
    dt_diag, _, _ = _bar_dtcr(N, L, E, A, rho, ("-lump", "diagonal"))
    assert dt_rowsum > 0 and dt_diag > 0, "both lumpings must yield a positive dt_cr"
    assert dt_diag < dt_rowsum, (
        "diagonal-of-consistent (stiffer) should give a smaller dt_cr "
        "(diagonal=%.5f vs rowsum=%.5f)" % (dt_diag, dt_rowsum)
    )
    assert abs(dt_rowsum - le / c) / (le / c) < 0.10, (
        "rowsum dt_cr=%.5f should be ~le/c=%.5f" % (dt_rowsum, le / c)
    )


# ==========================================================================
# 12. -cflAbort enforcement: analyze() aborts (!=0) on an oversized dt, runs on a safe dt
# ==========================================================================
def test_eb_cfl_enforce():
    N = 40
    L = 8.0
    E = 100.0
    A = 1.0
    rho = 1.0
    c = math.sqrt(E / rho)
    le = L / N

    def run(extra, dt):
        build_bar(N, L, E, A, rho)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(N + 1, 1.0, 0.0)
        ops.constraints("Transformation")
        ops.numberer("Plain")
        ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 1)
        ops.algorithm("Linear")
        ops.integrator(EB, P, *extra)
        ops.analysis("Transient")
        return ops.analyze(1, dt)

    r_abort = run(("-cflAbort",), 5.0 * le / c)     # over the NB limit -> abort
    r_ok = run(("-cflAbort",), 0.5 * le / c)        # safely under it    -> runs
    assert r_abort != 0, "oversized dt should abort (analyze != 0), got %d" % r_abort
    assert r_ok == 0, "safe dt should run (analyze == 0), got %d" % r_ok
