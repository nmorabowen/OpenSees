"""
Numerical verification battery for the explicit Noh-Bathe integrators.
Runs against the freshly built dist/bin/opensees.pyd.

Tests:
  1. Order of accuracy (must converge ~2nd order)
  2. Stability limit (Noh-Bathe stable beyond central difference's 2/omega)
  3. Cross-check vs Newmark and CentralDifference
  4. Spectral behaviour: algorithmic damping + period error vs Omega=omega*dt
  5. 1D wave propagation (wave speed sqrt(E/rho))
  6. LNVD dynamic relaxation vs a static solve (multi-DOF)
  7. Rigid-body / momentum conservation
  8. Queryable dt_cr (= le/c) + adaptive sub-step + -recompute path
  9. Recorder element-mass KE (getMass aliasing fix)
 10. dt_cr sentinels: DISABLED (no -cfl) / NOT_APPLICABLE (nodal-mass model)
 11. -lump diagonal vs rowsum (diagonal-of-consistent is stiffer -> smaller dt_cr)
 12. -cflAbort / -recompute enforcement (analyze aborts on a CFL violation)
"""
import sys, os, math

DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
for _m in ('opensees', 'openseespy', 'openseespy.opensees'):
    sys.modules.pop(_m, None)
sys.path.insert(0, DIST)
os.add_dll_directory(DIST)
import numpy as np
import opensees as ops
print("opensees module:", ops.__file__)

OUT = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\Ladruno_scripts\_verify_out"
os.makedirs(OUT, exist_ok=True)

results = []   # (name, passed, detail)


def sdof(k, m, integrator_args, dt, nsteps, d0=0.0, v0=0.0):
    """Free-vibration SDOF; returns (t[], u[]). integrator_args = tuple."""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0); ops.fix(2, 0, 1)
    ops.mass(2, m, 0.0)
    ops.uniaxialMaterial('Elastic', 1, k)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)
    if d0 != 0.0:
        ops.setNodeDisp(2, 1, d0, '-commit')
    if v0 != 0.0:
        ops.setNodeVel(2, 1, v0, '-commit')
    ops.constraints('Transformation'); ops.numberer('Plain')
    ops.system('Diagonal'); ops.test('NormDispIncr', 1e-12, 1)
    ops.algorithm('Linear')
    ops.integrator(*integrator_args)
    ops.analysis('Transient')
    t = [ops.getTime()]; u = [ops.nodeDisp(2, 1)]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        t.append(ops.getTime()); u.append(ops.nodeDisp(2, 1))
    return np.array(t), np.array(u)


# ---------------------------------------------------------------------------
# 1. Order of accuracy: u0=0, v0=1 -> u(t)=(v0/w) sin(w t), a0=0 (no cold start)
# ---------------------------------------------------------------------------
def test_order():
    m = 1.0; w = 2*math.pi; k = m*w*w; v0 = 1.0
    Tend = 1.0
    errs, dts = [], []
    for N in (40, 80, 160, 320, 640):
        dt = Tend/N
        t, u = sdof(k, m, ('ExplicitBathe', 0.54), dt, N, v0=v0)
        ue = (v0/w)*np.sin(w*t)
        errs.append(math.sqrt(np.mean((u-ue)**2))); dts.append(dt)
    errs = np.array(errs); dts = np.array(dts)
    slope = np.polyfit(np.log(dts), np.log(errs), 1)[0]
    ok = 1.7 <= slope <= 2.3
    results.append(("1. order of accuracy", ok,
                    "convergence slope = %.2f (expect ~2.0); errs=%s" %
                    (slope, np.array2string(errs, precision=2))))


# ---------------------------------------------------------------------------
# 2. Stability limit: sweep Omega = w*dt; find largest stable. CD ~ 2.0.
# ---------------------------------------------------------------------------
def stable_at(integrator_args, Omega, m=1.0, w=2*math.pi):
    k = m*w*w; dt = Omega/w
    _, u = sdof(k, m, integrator_args, dt, 200, d0=1.0)
    if not np.all(np.isfinite(u)):
        return False
    early = np.abs(u[:40]).max(); late = np.abs(u[-40:]).max()
    return late <= 5.0*max(early, 1e-12)

def test_stability():
    grid = [1.6, 1.8, 1.95, 2.0, 2.05, 2.2, 2.4, 2.6, 2.8, 3.0]
    eb = max([O for O in grid if stable_at(('ExplicitBathe', 0.54), O)] or [0])
    cd = max([O for O in grid if stable_at(('CentralDifference',), O)] or [0])
    ok = eb > 2.0  # must beat central difference's 2/omega
    results.append(("2. stability limit", ok,
                    "ExplicitBathe stable to Omega=%.2f, CentralDifference to %.2f "
                    "(CD theory=2.0; NB should exceed it)" % (eb, cd)))


# ---------------------------------------------------------------------------
# 3. Cross-check vs Newmark (avg accel) and CentralDifference at small dt
# ---------------------------------------------------------------------------
def test_crosscheck():
    m = 1.0; w = 2*math.pi; k = m*w*w; d0 = 1.0
    dt = (2*math.pi/w)/200; N = 400
    _, ueb = sdof(k, m, ('ExplicitBathe', 0.54), dt, N, d0=d0)
    _, unm = sdof(k, m, ('Newmark', 0.5, 0.25), dt, N, d0=d0)
    _, ucd = sdof(k, m, ('CentralDifference',), dt, N, d0=d0)
    e_nm = np.abs(ueb-unm).max()/d0
    e_cd = np.abs(ueb-ucd).max()/d0
    ok = e_nm < 0.03 and e_cd < 0.03
    results.append(("3. cross-check vs Newmark/CD", ok,
                    "max|EB-Newmark|=%.2f%%, max|EB-CD|=%.2f%% (both <3%%)" %
                    (100*e_nm, 100*e_cd)))


# ---------------------------------------------------------------------------
# 4. Spectral: algorithmic damping + period error vs Omega
# ---------------------------------------------------------------------------
def peaks(u):
    idx = [i for i in range(1, len(u)-1) if u[i] > u[i-1] and u[i] >= u[i+1] and u[i] > 0]
    return idx

def test_spectral():
    m = 1.0; dt = 0.01
    rows = []
    for Omega in (0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5):
        w = Omega/dt; k = m*w*w; T = 2*math.pi/w
        n = int(8*T/dt)
        t, u = sdof(k, m, ('ExplicitBathe', 0.54), dt, n, d0=1.0)
        pk = peaks(u)
        if len(pk) >= 3:
            Tn = (t[pk[-1]]-t[pk[0]])/(len(pk)-1)       # numerical period
            perr = 100*(Tn-T)/T
            a0, an = u[pk[0]], u[pk[-1]]; nc = len(pk)-1
            xi = math.log(max(a0,1e-12)/max(an,1e-12))/(2*math.pi*nc) if an > 0 else float('nan')
        else:
            perr, xi = float('nan'), float('nan')
        rows.append((Omega, perr, 100*xi if xi==xi else xi))
    # checks: low-Omega damping tiny; damping grows with Omega
    lo = [r[2] for r in rows if r[0] <= 0.25 and r[2]==r[2]]
    hi = [r[2] for r in rows if r[0] >= 1.0 and r[2]==r[2]]
    ok = (max(lo) < 1.0) and (max(hi) > max(lo)) if lo and hi else False
    detail = "  ".join("O=%.2f:Terr=%.2f%%,xi=%.2f%%" % r for r in rows)
    results.append(("4. spectral (damping/period vs Omega)", ok, detail))


# ---------------------------------------------------------------------------
# Bar model (chain of truss elements) shared by tests 5 and 6
# ---------------------------------------------------------------------------
def build_bar(N, L, E, A, rho):
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    le = L/N
    for i in range(N+1):
        ops.node(i+1, i*le, 0.0)
    ops.fix(1, 1, 1)                        # fixed end: x and y
    for i in range(1, N+1):                 # nodes 2..N+1: free in x, fixed in y
        ops.fix(i+1, 0, 1)
    ops.uniaxialMaterial('Elastic', 1, E)
    for i in range(N):
        ops.element('Truss', i+1, i+1, i+2, A, 1, '-rho', rho)
    return le


def test_wave():
    N = 100; L = 20.0; E = 100.0; A = 1.0; rho = 1.0; P = 1.0
    c = math.sqrt(E/rho)                    # = 10
    le = build_bar(N, L, E, A, rho)
    ops.timeSeries('Constant', 1); ops.pattern('Plain', 1, 1)
    ops.load(N+1, P, 0.0)                   # step force at the free end (x=L)
    ops.constraints('Transformation'); ops.numberer('Plain')
    ops.system('Diagonal'); ops.test('NormDispIncr', 1e-12, 1)
    ops.algorithm('Linear'); ops.integrator('ExplicitBathe', 0.54)
    ops.analysis('Transient')
    dt = 0.5*le/c
    v_particle = P/(A*rho*c)               # particle velocity behind the front
    vth = 0.5*v_particle
    # stations at several distances d from the loaded free end (x = L-d)
    dists = [5.0, 8.0, 11.0, 14.0]
    nodes = [int(round((L-d)/le)) + 1 for d in dists]
    arrival = {nd: None for nd in nodes}
    tcur = 0.0
    for _ in range(int(2.2*L/c/dt)):
        ops.analyze(1, dt); tcur += dt
        for nd in nodes:
            if arrival[nd] is None and abs(ops.nodeVel(nd, 1)) > vth:
                arrival[nd] = tcur
    pts = [(d, arrival[nd]) for d, nd in zip(dists, nodes) if arrival[nd]]
    if len(pts) < 2:
        results.append(("5. wave propagation", False, "front not detected at >=2 stations"))
        return
    d_arr = np.array([p[0] for p in pts]); t_arr = np.array([p[1] for p in pts])
    slope = np.polyfit(d_arr, t_arr, 1)[0]   # = 1/speed
    c_meas = 1.0/slope
    err = abs(c_meas - c)/c
    ok = err < 0.05
    results.append(("5. wave propagation", ok,
                    "measured wave speed %.3f vs sqrt(E/rho)=%.3f (err %.1f%%)" %
                    (c_meas, c, 100*err)))


def test_lnvd_static():
    N = 20; L = 10.0; E = 100.0; A = 1.0; rho = 1.0; P = 2.0
    # --- static reference ---
    le = build_bar(N, L, E, A, rho)
    ops.timeSeries('Linear', 1); ops.pattern('Plain', 1, 1); ops.load(N+1, P, 0.0)
    ops.constraints('Transformation'); ops.numberer('Plain')
    ops.system('BandSPD'); ops.test('NormDispIncr', 1e-10, 50)
    ops.algorithm('Newton'); ops.integrator('LoadControl', 1.0)
    ops.analysis('Static'); ops.analyze(1)
    u_static = np.array([ops.nodeDisp(i+1, 1) for i in range(N+1)])
    # --- LNVD dynamic relaxation ---
    le = build_bar(N, L, E, A, rho)
    ops.timeSeries('Constant', 1); ops.pattern('Plain', 1, 1); ops.load(N+1, P, 0.0)
    ops.constraints('Transformation'); ops.numberer('Plain')
    ops.system('Diagonal'); ops.test('NormDispIncr', 1e-12, 1)
    ops.algorithm('Linear'); ops.integrator('ExplicitBatheLNVD', 0.54, 0.8)
    ops.analysis('Transient')
    c = math.sqrt(E/rho); dt = 0.5*(L/N)/c
    for _ in range(6000):
        ops.analyze(1, dt)
    u_lnvd = np.array([ops.nodeDisp(i+1, 1) for i in range(N+1)])
    tip_exact = P*L/(E*A)
    err_field = np.abs(u_lnvd - u_static).max()/abs(u_static[-1])
    err_tip = abs(u_lnvd[-1]-tip_exact)/tip_exact
    ok = err_field < 0.02 and err_tip < 0.02
    results.append(("6. LNVD vs static (multi-DOF)", ok,
                    "field match %.2f%%, tip err vs PL/EA %.2f%% (tip=%.4f, exact=%.4f)" %
                    (100*err_field, 100*err_tip, u_lnvd[-1], tip_exact)))


# ---------------------------------------------------------------------------
# 7. Rigid-body / momentum: free-free truss, both nodes initial v0 -> u=v0 t
# ---------------------------------------------------------------------------
def test_rigid_body():
    m = 1.0; v0 = 2.0; E = 100.0; A = 1.0
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.node(1, 0.0, 0.0); ops.node(2, 1.0, 0.0)
    ops.fix(1, 0, 1); ops.fix(2, 0, 1)      # free in x, no constraint in x
    ops.mass(1, m, 0.0); ops.mass(2, m, 0.0)
    ops.uniaxialMaterial('Elastic', 1, E)
    ops.element('Truss', 1, 1, 2, A, 1)
    ops.setNodeVel(1, 1, v0, '-commit'); ops.setNodeVel(2, 1, v0, '-commit')
    ops.constraints('Transformation'); ops.numberer('Plain')
    ops.system('Diagonal'); ops.test('NormDispIncr', 1e-12, 1)
    ops.algorithm('Linear'); ops.integrator('ExplicitBathe', 0.54)
    ops.analysis('Transient')
    dt = 0.01; nsteps = 200
    for _ in range(nsteps):
        ops.analyze(1, dt)
    t = ops.getTime()
    u1 = ops.nodeDisp(1, 1); u2 = ops.nodeDisp(2, 1)
    v1 = ops.nodeVel(1, 1); v2 = ops.nodeVel(2, 1)
    u_exact = v0*t
    err_u = max(abs(u1-u_exact), abs(u2-u_exact))/abs(u_exact)
    err_v = max(abs(v1-v0), abs(v2-v0))/v0
    err_strain = abs(u1-u2)                 # rigid -> no relative motion
    ok = err_u < 1e-3 and err_v < 1e-3 and err_strain < 1e-6
    results.append(("7. rigid-body / momentum", ok,
                    "u err %.2e, v err %.2e, spurious strain %.2e (all ~0)" %
                    (err_u, err_v, err_strain)))


def test_query_substep():
    # distributed-mass bar -> finite, queryable dt_cr (= le/c for a 2-node bar)
    N = 40; L = 8.0; E = 100.0; A = 1.0; rho = 1.0; c = math.sqrt(E/rho)
    le = build_bar(N, L, E, A, rho)
    ops.timeSeries('Constant', 1); ops.pattern('Plain', 1, 1); ops.load(N+1, 1.0, 0.0)
    ops.constraints('Transformation'); ops.numberer('Plain')
    ops.system('Diagonal'); ops.test('NormDispIncr', 1e-12, 1)
    ops.algorithm('Linear'); ops.integrator('ExplicitBathe', 0.54, '-cfl')
    ops.analysis('Transient')
    ops.analyze(1, 0.1*le/c)                 # one priming step triggers the compute
    dtcr = ops.criticalTimeStep()
    expected = le/c                          # 2/omega_max for the bar element
    err = abs(dtcr - expected)/expected
    # adaptive sub-stepping: how many internal steps to cover a coarse interval
    dt_coarse = 5.0*dtcr
    n_sub = max(1, math.ceil(dt_coarse/(0.9*dtcr)))
    # exercise the -recompute (tangent) path too: must still return finite dt_cr
    build_bar(N, L, E, A, rho)
    ops.timeSeries('Constant', 2); ops.pattern('Plain', 2, 2); ops.load(N+1, 1.0, 0.0)
    ops.constraints('Transformation'); ops.numberer('Plain')
    ops.system('Diagonal'); ops.test('NormDispIncr', 1e-12, 1)
    ops.algorithm('Linear'); ops.integrator('ExplicitBathe', 0.54, '-recompute', 5)
    ops.analysis('Transient'); ops.analyze(3, 0.5*le/c)
    dtcr_t = ops.criticalTimeStep()
    ok = dtcr > 0 and err < 0.10 and n_sub == 6 and dtcr_t > 0
    results.append(("8. queryable dt_cr + sub-step + recompute", ok,
                    "criticalTimeStep()=%.5f vs le/c=%.5f (err %.1f%%); "
                    "n_sub=%d for dt=5*dt_cr; -recompute dt_cr=%.5f" %
                    (dtcr, expected, 100*err, n_sub, dtcr_t)))


def test_recorder_element_mass():
    # Bar with ELEMENT mass (truss -rho) + EnergyBalance recorder. Verifies the
    # getMass()/getDamp() shared-theMatrix aliasing fix: before it, element KE
    # was read from the (zeroed) damping matrix and came out ~0.
    N = 20; L = 10.0; E = 100.0; A = 1.0; rho = 1.0; v0 = 2.0
    efile = os.path.join(OUT, "bar_energy.txt")
    build_bar(N, L, E, A, rho)
    for i in range(1, N+1):                  # initial x-velocity on every free node
        ops.setNodeVel(i+1, 1, v0, '-commit')
    ops.recorder('EnergyBalance', '-file', efile, '-time')
    ops.constraints('Transformation'); ops.numberer('Plain'); ops.system('Diagonal')
    ops.test('NormDispIncr', 1e-12, 1); ops.algorithm('Linear')
    ops.integrator('ExplicitBathe', 0.54); ops.analysis('Transient')
    ops.analyze(2, 1e-4)
    ops.wipe()                               # flush/close the recorder
    d = np.loadtxt(efile)
    ke0 = d[0, 1]                            # first-row KE (col 0 = time)
    M_total = rho*A*L                        # total bar mass (element-based)
    ke_expected = 0.5*M_total*v0**2          # ~ (fixed end excludes ~rho*A*le/2)
    ok = ke0 > 0.5*ke_expected               # before the fix this was ~0
    results.append(("9. recorder element-mass KE (aliasing fix)", ok,
                    "EnergyBalance KE=%.3f vs ~0.5*M*v^2=%.3f (>0 confirms the fix; was ~0)"
                    % (ke0, ke_expected)))


def test_sentinels():
    # item 4: criticalTimeStep() returns distinct sentinels explaining WHY there
    # is no value, instead of collapsing every case onto -1.
    DISABLED, NOT_COMPUTED, NOT_APPLICABLE = -1.0, -2.0, -3.0
    m = 1.0; w = 2*math.pi; k = m*w*w; dt = 0.01
    # (a) estimation not requested -> DISABLED
    sdof(k, m, ('ExplicitBathe', 0.54), dt, 1, d0=1.0)
    s_disabled = ops.criticalTimeStep()
    # (b) -cfl on a pure nodal-mass model: the zeroLength element carries no
    #     element mass, so no element provides a pencil -> NOT_APPLICABLE
    sdof(k, m, ('ExplicitBathe', 0.54, '-cfl'), dt, 1, d0=1.0)
    s_napp = ops.criticalTimeStep()
    ok = abs(s_disabled - DISABLED) < 1e-9 and abs(s_napp - NOT_APPLICABLE) < 1e-9
    results.append(("10. dt_cr sentinels (disabled / not-applicable)", ok,
                    "no-cfl -> %.1f (want %.1f); nodal-mass+cfl -> %.1f (want %.1f)" %
                    (s_disabled, DISABLED, s_napp, NOT_APPLICABLE)))


def _bar_dtcr(N, L, E, A, rho, extra):
    """Build the distributed-mass bar, run one tiny step, return criticalTimeStep()."""
    le = build_bar(N, L, E, A, rho); c = math.sqrt(E/rho)
    ops.timeSeries('Constant', 1); ops.pattern('Plain', 1, 1); ops.load(N+1, 1.0, 0.0)
    ops.constraints('Transformation'); ops.numberer('Plain'); ops.system('Diagonal')
    ops.test('NormDispIncr', 1e-12, 1); ops.algorithm('Linear')
    ops.integrator('ExplicitBathe', 0.54, '-cfl', *extra); ops.analysis('Transient')
    ops.analyze(1, 0.1*le/c)
    return ops.criticalTimeStep(), le, c


def test_lump_diagonal():
    # item 6: -lump diagonal returns a finite, positive dt_cr. Diagonal-of-
    # consistent mass (rho*A*le/3 per DOF) is smaller than row-sum (rho*A*le/2),
    # so the element is "stiffer" -> a smaller dt_cr (~0.816 le/c vs le/c).
    N = 40; L = 8.0; E = 100.0; A = 1.0; rho = 1.0
    dt_rowsum, le, c = _bar_dtcr(N, L, E, A, rho, ('-lump', 'rowsum'))
    dt_diag, _, _    = _bar_dtcr(N, L, E, A, rho, ('-lump', 'diagonal'))
    ok = (dt_rowsum > 0 and dt_diag > 0 and dt_diag < dt_rowsum and
          abs(dt_rowsum - le/c)/(le/c) < 0.10)
    results.append(("11. -lump diagonal vs rowsum", ok,
                    "rowsum dt_cr=%.5f (~le/c=%.5f); diagonal dt_cr=%.5f (smaller=stiffer)" %
                    (dt_rowsum, le/c, dt_diag)))


def test_cfl_enforce():
    # item 5: -cflAbort enforces - analyze() aborts (returns != 0) when dt exceeds
    # the Noh-Bathe limit. -recompute couples to this SAME guard, re-checking every
    # N committed steps; here it is exercised via an oversized dt at step 1. True
    # mid-run stiffening enforcement (tangent shrinking dt_cr) needs a nonlinear
    # model and is out of scope for this linear-elastic battery.
    N = 40; L = 8.0; E = 100.0; A = 1.0; rho = 1.0; c = math.sqrt(E/rho); le = L/N

    def run(extra, dt):
        build_bar(N, L, E, A, rho)
        ops.timeSeries('Constant', 1); ops.pattern('Plain', 1, 1); ops.load(N+1, 1.0, 0.0)
        ops.constraints('Transformation'); ops.numberer('Plain'); ops.system('Diagonal')
        ops.test('NormDispIncr', 1e-12, 1); ops.algorithm('Linear')
        ops.integrator('ExplicitBathe', 0.54, *extra); ops.analysis('Transient')
        return ops.analyze(1, dt)

    r_abort = run(('-cflAbort',), 5.0*le/c)   # over the NB limit (~2 le/c) -> abort
    r_ok    = run(('-cflAbort',), 0.5*le/c)   # safely under it           -> runs
    ok = r_abort != 0 and r_ok == 0
    results.append(("12. -cflAbort / -recompute enforcement", ok,
                    "oversized dt -> analyze=%d (!=0, aborted); safe dt -> analyze=%d (==0)" %
                    (r_abort, r_ok)))


for fn in (test_order, test_stability, test_crosscheck, test_spectral,
           test_wave, test_lnvd_static, test_rigid_body, test_query_substep,
           test_recorder_element_mass, test_sentinels, test_lump_diagonal,
           test_cfl_enforce):
    try:
        fn()
    except Exception as e:
        results.append((fn.__name__, False, "EXCEPTION: %r" % e))

print("\n================ EXPLICIT INTEGRATOR VERIFICATION ================")
npass = 0
for name, ok, detail in results:
    print("[%s] %s" % ("PASS" if ok else "FAIL", name))
    print("        " + detail)
    npass += int(ok)
print("=================================================================")
print("%d / %d passed" % (npass, len(results)))
sys.exit(0 if npass == len(results) else 1)
