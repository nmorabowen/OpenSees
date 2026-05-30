"""
Explicit Noh-Bathe integrator + energy-balance diagnostics -- worked example.

Demonstrates the "explicit recipe" and pairs ExplicitBathe with the
EnergyBalanceRecorder so spurious energy growth is visible. Run with the
project venv:

    & "C:\\Users\\nmora\\venv\\opensees_venv\\Scripts\\python.exe" explicit_bathe_example.py

Explicit recipe (what makes these integrators behave):
  * lumped / nodal mass            (the RHS is just M; consistent mass wastes it)
  * system Diagonal                (trivial M^-1; never a sparse solver here)
  * algorithm Linear               (the scheme does exactly 2 solves per step)
  * numberer Plain, test on 1 iter
  * pick dt below the reported critical step (use -cfl to print it)
"""

import openseespy.opensees as ops


# ---------------------------------------------------------------------------
# 1) Free vibration of an undamped SDOF -- energy should be conserved
#    (KE <-> IE trade off each cycle, RES ~ 0).
# ---------------------------------------------------------------------------
def sdof_free_vibration():
    k, m, d0 = 100.0, 1.0, 0.10        # omega = sqrt(k/m) = 10 rad/s
    # central-difference limit dt_cd = 2/omega = 0.2 ; Noh-Bathe ~ 0.4

    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)

    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0); ops.fix(2, 0, 1)        # free in x only
    ops.mass(2, m, 0.0)

    ops.uniaxialMaterial('Elastic', 1, k)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)

    # initial displacement, released from rest
    ops.setNodeDisp(2, 1, d0, '-commit')

    # whole-model + (optional) per-region energy balance, written as a sidecar
    ops.recorder('EnergyBalance', '-file', 'sdof_energy.txt', '-time')

    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('Diagonal')                 # trivial diagonal solve
    ops.test('NormDispIncr', 1e-12, 1)     # single pass (Linear)
    ops.algorithm('Linear')
    # p = 0.54 (recommended); -cfl prints the conservative + Noh-Bathe dt limits
    ops.integrator('ExplicitBathe', 0.54, '-cfl')
    ops.analysis('Transient')

    dt = 0.05                               # < dt_cd = 0.2, comfortably stable
    for _ in range(126):                    # ~2 periods (T = 0.628 s)
        if ops.analyze(1, dt) != 0:
            print("analysis failed/aborted (check energy residual)")
            break

    print("wrote sdof_energy.txt -- columns: time KE IE DW ULW RES ERR%")
    print("for a stable run RES stays ~0 and ERR% stays small.")


# ---------------------------------------------------------------------------
# 2) Quasi-static equilibrium via LNVD dynamic relaxation.
#    FLAC local damping bleeds kinetic energy; watch |unbalance| -> 0.
#    (Use -verbose to print the per-step unbalance norm.)
# ---------------------------------------------------------------------------
def lnvd_dynamic_relaxation():
    k, m, P = 100.0, 1.0, 5.0              # static answer: u = P/k = 0.05

    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)
    ops.node(1, 0.0, 0.0); ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0); ops.fix(2, 0, 1)
    ops.mass(2, m, 0.0)
    ops.uniaxialMaterial('Elastic', 1, k)
    ops.element('zeroLength', 1, 1, 2, '-mat', 1, '-dir', 1)

    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, P, 0.0)

    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('Diagonal')
    ops.test('NormDispIncr', 1e-12, 1)
    ops.algorithm('Linear')
    # p = 0.54, alpha = 0.8 (classic FLAC); damps toward static equilibrium
    ops.integrator('ExplicitBatheLNVD', 0.54, 0.8)
    ops.analysis('Transient')

    for _ in range(400):
        ops.analyze(1, 0.05)
    print("LNVD relaxed u(2) = %.5f  (target P/k = %.5f)" %
          (ops.nodeDisp(2, 1), P / k))


# ---------------------------------------------------------------------------
# 3) Adaptive sub-stepping for a time history.
#    You keep YOUR dt (e.g. the ground-motion sampling step) for loads and
#    output, but query the integrator's critical step and sub-divide each
#    interval so every internal step is stable. Each sub-step is a real,
#    fully-committed step (correct for nonlinear materials) -- this is the
#    right way to "lower dt if it's below critical".
# ---------------------------------------------------------------------------
def adaptive_substep_skeleton(dt_record, n_intervals, safety=0.9):
    # ... build model with distributed (element) mass, system Diagonal,
    #     algorithm Linear, integrator ExplicitBathe with -cfl or -recompute ...
    ops.integrator('ExplicitBathe', 0.54, '-recompute', 50)   # refresh dt_cr
    ops.analysis('Transient')
    ops.analyze(1, 1e-9)                       # prime: triggers the dt_cr compute
    import math
    for _ in range(n_intervals):
        dt_cr = ops.criticalTimeStep()         # queryable (feature a)
        if dt_cr <= 0:                         # nodal-mass-only model: no element limit
            n_sub = 1
        else:
            n_sub = max(1, math.ceil(dt_record / (safety * dt_cr)))
        ops.analyze(n_sub, dt_record / n_sub)  # stable, fully-committed sub-steps
        # Notes:
        # - n_sub is re-queried EACH interval, so it adapts as the model softens
        #   or stiffens -- but it is FIXED within an interval, lagged by the
        #   -recompute period. A sudden mid-interval stiffening (e.g. contact
        #   impact) is not caught until the next interval; pair with -cflAbort if
        #   you need the run to stop rather than go unstable on such an event.
        # - criticalTimeStep() returns the CONSERVATIVE central-difference limit
        #   (2/omega). Noh-Bathe is stable to ~2x that, so safety=0.9 here is
        #   ~2x extra-safe; if you trust the Noh-Bathe factor you can raise it.


if __name__ == '__main__':
    sdof_free_vibration()
    lnvd_dynamic_relaxation()
