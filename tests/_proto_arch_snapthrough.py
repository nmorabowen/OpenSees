"""ADR-20 first-move prototype — snap-through oracle harness.

Minimal limit-point model: a shallow von Mises truss (two corotational bars,
apex loaded downward). As the apex is pushed past the flat configuration the
structure snaps through, so the load-displacement curve has a limit point where
the tangent stiffness vanishes and then goes negative.

Establishes the three things ADR-20 §6 needs as test fixtures BEFORE any C++:
  A. LoadControl  -> FAILS at the limit point          (the problem is real)
  B. ArcLength    -> traces the full equilibrium path   (AL-5 reference path,
                     including the unstable descending branch)
  C. Dynamic      -> snaps PAST the instability to the   (AL-6 oracle for the
     relaxation       far stable branch under damping        future -stabilize mode)

Run:
  $env:PYTHONPATH=<main>/dist/bin; py -3.12 tests/_proto_arch_snapthrough.py
"""
try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops


# ---- model ---------------------------------------------------------------
# supports at (-1,0),(1,0); apex at (0,h) shallow. apex free in Y only.
# E*A chosen so the snap-through limit load is O(1):
#   P_cr ~ E*A*(H/L)^3 ~ 3.8 for these numbers.
H = 0.10          # apex rise (shallow)
SPAN = 1.0        # half-span
E = 1.0e4
A = 1.0
APEX = 2          # node tag of the loaded apex


def build_geom():
    """Geometry + material + elements only — each runner adds its own load."""
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.node(1, -SPAN, 0.0)
    ops.node(APEX, 0.0, H)
    ops.node(3,  SPAN, 0.0)
    ops.fix(1, 1, 1)
    ops.fix(3, 1, 1)
    ops.fix(APEX, 1, 0)          # apex: X fixed by symmetry, Y free
    ops.uniaxialMaterial('Elastic', 1, E)
    ops.element('corotTruss', 1, 1, APEX, A, 1)
    ops.element('corotTruss', 2, APEX, 3, A, 1)


def add_ref_load():
    """Unit downward load at apex, scaled by the load factor lambda (static)."""
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(APEX, 0.0, -1.0)


def add_constant_load(P):
    """Time-invariant downward load of magnitude P (for the transient run)."""
    ops.timeSeries('Constant', 2)
    ops.pattern('Plain', 2, 2)
    ops.load(APEX, 0.0, -P)


def _common_static():
    ops.system('BandGeneral')
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.test('NormDispIncr', 1.0e-8, 50, 0)
    ops.algorithm('Newton')


# ---- A. load control: expect failure at the limit point ------------------
def run_load_control():
    build_geom()
    add_ref_load()
    _common_static()
    dLam = 0.05
    ops.integrator('LoadControl', dLam)
    ops.analysis('Static')
    last_lam, last_uy = 0.0, 0.0
    failed_at = None
    for i in range(200):
        ok = ops.analyze(1)
        if ok != 0:
            failed_at = (ops.getTime(), ops.nodeDisp(APEX, 2))
            break
        last_lam, last_uy = ops.getTime(), ops.nodeDisp(APEX, 2)
    return dict(last_converged_lambda=last_lam, last_uy=last_uy,
                failed=failed_at is not None, fail_state=failed_at)


# ---- B. reference path: DisplacementControl traces the FULL curve --------
# apex-Y is monotonic through the snap-through, so controlling it parametrizes
# the whole equilibrium path incl. the unstable (negative-stiffness) branch.
def run_reference_path(du=-0.005, nsteps=120):
    build_geom()
    add_ref_load()
    _common_static()
    ops.integrator('DisplacementControl', APEX, 2, du)
    ops.analysis('Static')
    # the snap-through LIMIT load is the FIRST local maximum of lambda(u) — past
    # it the path turns down (unstable branch) before the far branch re-stiffens
    # (where lambda climbs again, so the global max is NOT the limit load).
    path = []
    limit_lam, limit_uy = -1e30, 0.0
    frozen = False
    saw_negative_stiffness = False
    prev = (0.0, 0.0)
    for i in range(nsteps):
        if ops.analyze(1) != 0:
            break
        lam, uy = ops.getTime(), ops.nodeDisp(APEX, 2)
        path.append((lam, uy))
        if not frozen:
            if lam < prev[0] - 1e-9:        # first decrease -> limit point passed
                frozen = True
            elif lam > limit_lam:
                limit_lam, limit_uy = lam, uy
        if lam < prev[0] - 1e-9 and uy < prev[1]:   # load drops as disp grows
            saw_negative_stiffness = True
        prev = (lam, uy)
    final = path[-1] if path else (0.0, 0.0)
    return dict(nsteps=len(path), peak_lambda=limit_lam, peak_uy=limit_uy,
                final_lambda=final[0], final_uy=final[1],
                traced_unstable_branch=saw_negative_stiffness)


# ---- B'. stock ArcLength, tuned (load-factor de-weighted via alpha) -------
def run_arc_length(arc=0.01, alpha=0.1, nsteps=600):
    build_geom()
    add_ref_load()
    _common_static()
    ops.integrator('ArcLength', arc, alpha)
    ops.analysis('Static')
    peak_lam, nsuccess = -1e30, 0
    saw_negative_stiffness = False
    prev = (0.0, 0.0)
    for i in range(nsteps):
        if ops.analyze(1) != 0:
            break
        nsuccess += 1
        lam, uy = ops.getTime(), ops.nodeDisp(APEX, 2)
        peak_lam = max(peak_lam, lam)
        if lam < prev[0] - 1e-9 and uy < prev[1]:
            saw_negative_stiffness = True
        prev = (lam, uy)
    return dict(nsteps=nsuccess, peak_lambda=peak_lam,
                traced_unstable_branch=saw_negative_stiffness)


# ---- C. dynamic relaxation: snap PAST the instability --------------------
def run_dynamic_relaxation(P_over_crit, m=1.0, damp=8.0, dt=2e-3, nsteps=4000):
    build_geom()
    # constant super-critical load: no equilibrium on the near branch, so the
    # apex must dynamically snap through and damp out on the far stable branch.
    ops.mass(APEX, m, m)
    add_constant_load(P_over_crit)
    _common_static()
    ops.rayleigh(damp, 0.0, 0.0, 0.0)          # mass-proportional damping
    ops.integrator('Newmark', 0.5, 0.25)
    ops.analysis('Transient')
    uy_hist = []
    for i in range(nsteps):
        ok = ops.analyze(1, dt)
        if ok != 0:
            return dict(converged=False, step=i, uy=ops.nodeDisp(APEX, 2))
        if i % 200 == 0:
            uy_hist.append(ops.nodeDisp(APEX, 2))
    uy_final = ops.nodeDisp(APEX, 2)
    return dict(converged=True, uy_final=uy_final,
                snapped_through=uy_final < -2.0 * H,   # below the far flat
                uy_samples=[round(u, 4) for u in uy_hist[:8]])


if __name__ == '__main__':
    print("=" * 64)
    print("ADR-20 snap-through oracle harness  (H=%.3f span=%.1f E=%.0e)" % (H, SPAN, E))
    print("=" * 64)

    a = run_load_control()
    print("\n[A] LoadControl (the problem):")
    print("    last converged lambda = %.5f at uy=%.5f" % (a['last_converged_lambda'], a['last_uy']))
    print("    FAILED at limit point = %s   fail_state=%s" % (a['failed'], a['fail_state']))

    b = run_reference_path()
    print("\n[B] DisplacementControl reference path (AL-5):")
    print("    steps traced           = %d" % b['nsteps'])
    print("    PEAK (limit) lambda    = %.5f at uy=%.5f" % (b['peak_lambda'], b['peak_uy']))
    print("    final state            = lambda=%.5f uy=%.5f" % (b['final_lambda'], b['final_uy']))
    print("    traced unstable branch = %s" % b['traced_unstable_branch'])

    bp = run_arc_length()
    print("\n[B'] stock ArcLength, tuned (alpha=0.1, arc=0.01):")
    print("    steps traced           = %d" % bp['nsteps'])
    print("    PEAK (limit) lambda    = %.5f" % bp['peak_lambda'])
    print("    traced unstable branch = %s" % bp['traced_unstable_branch'])

    crit = b['peak_lambda']
    c = run_dynamic_relaxation(P_over_crit=1.15 * crit)
    print("\n[C] Dynamic relaxation @1.15*crit (AL-6 stabilization oracle):")
    print("    %s" % c)

    print("\n" + "-" * 64)
    print("SUMMARY")
    print("  A load control fails at limit  : %s" % a['failed'])
    print("  B reference limit load (DispC) : %.5f" % b['peak_lambda'])
    print("  B traced unstable branch       : %s" % b['traced_unstable_branch'])
    print("  B' ArcLength limit load        : %.5f  (unstable=%s)" % (bp['peak_lambda'], bp['traced_unstable_branch']))
    print("  C snapped through to far side  : %s" % c.get('snapped_through'))
    print("-" * 64)
