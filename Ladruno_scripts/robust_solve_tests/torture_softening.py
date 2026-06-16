"""Torture fixture #2 - material softening (the clean load-control killer).

A single truss of Concrete02 pushed in compression past its peak strength. The
load-displacement envelope rises to a peak (A*|fpc|) then SOFTENS monotonically
to a residual - there is NO equilibrium at a load above the peak and NO far
branch to jump to. So LOAD control dies hard at the peak (and step-halving + the
whole algorithm ladder cannot help: the problem is not stiffness, it is that the
target state does not exist under load control). Only a constraint that follows
displacement (rung 3) traces the descending branch.

This is the RC-relevant motivator (concrete/UC softening) and the cleanest
binary in the battery. Needs NO energy recorder. Run with the dist build:
    python Ladruno_scripts/robust_solve_tests/torture_softening.py
"""
import os
import sys

DIST = os.environ.get("LADRUNO_DIST",
                      r"C:\Users\nmb\Documents\Github\OpenSees\dist\bin")
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)
    sys.path.insert(0, DIST)
import opensees as ops  # noqa: E402

L, A = 1.0, 1.0
FPC, EPSC0, FPCU, EPSU = -30.0, -0.002, -6.0, -0.012
PEAK = A * abs(FPC)                 # peak compressive load = 30
EPS_TARGET = -0.010                 # push well onto the softening branch


def build():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, L)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Concrete02", 1, FPC, EPSC0, FPCU, EPSU, 0.1, 3.0, 300.0)
    ops.element("Truss", 1, 1, 2, A, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, -1.0)              # compressive reference load
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGen")
    ops.test("NormDispIncr", 1e-8, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def strain():
    return ops.nodeDisp(2, 1) / L


def on_softening(eps, lam):
    return eps < 1.5 * EPSC0 and lam < 0.95 * PEAK   # past peak AND load dropped


def run_loadcontrol(dlam, lam_target, max_steps):
    build()
    ops.integrator("LoadControl", dlam)
    peak = 0.0
    soft = False
    eps_min = 0.0
    for _ in range(max_steps):
        if ops.analyze(1) != 0:
            break
        lam = ops.getLoadFactor(1)
        peak = max(peak, lam)
        eps_min = min(eps_min, strain())
        if on_softening(strain(), lam):
            soft = True
        if lam >= lam_target - 1e-9:
            break
    return peak, eps_min, soft


def run_adaptive_loadcontrol(dlam, lam_target):
    from ladruno_solve import adaptive_static
    build()
    n = max(1, round(lam_target / dlam))
    res = adaptive_static(
        ops, lambda s: ops.integrator("LoadControl", dlam * s), n_steps=n)
    return ops.getLoadFactor(1), strain(), bool(res)


def run_dispcontrol(du, max_steps):
    build()
    ops.integrator("DisplacementControl", 2, 1, du)
    soft = False
    eps_min = 0.0
    lam_at_soft = None
    for _ in range(max_steps):
        if ops.analyze(1) != 0:
            break
        lam = ops.getLoadFactor(1)
        eps_min = min(eps_min, strain())
        if on_softening(strain(), lam):
            soft = True
            lam_at_soft = lam
        if strain() <= EPS_TARGET:
            break
    return eps_min, soft, lam_at_soft


def run_robust(dlam, du, max_steps):
    build()
    ops.integrator("LoadControl", dlam)
    mode = "LoadControl"
    soft = False
    eps_min = 0.0
    switches = 0
    for _ in range(max_steps):
        if ops.analyze(1) == 0:
            lam = ops.getLoadFactor(1)
            eps_min = min(eps_min, strain())
            if on_softening(strain(), lam):
                soft = True
            if strain() <= EPS_TARGET:
                break
            continue
        if mode == "LoadControl":                 # rung-3 constraint switch
            ops.integrator("DisplacementControl", 2, 1, du)
            mode = "DisplacementControl"
            switches += 1
            continue
        break
    return eps_min, soft, switches


if __name__ == "__main__":
    lam_T = 1.3 * PEAK
    dlam = lam_T / 200.0
    DU = -5.0e-5
    MAXS = 800

    a_peak, a_eps, a_soft = run_loadcontrol(dlam, lam_T, 1000)
    b_lam, b_eps, b_ok = run_adaptive_loadcontrol(dlam, lam_T)
    c_eps, c_soft, c_lam = run_dispcontrol(DU, MAXS)
    d_eps, d_soft, d_sw = run_robust(dlam, DU, MAXS)

    print(f"peak strength A*|fpc| = {PEAK:.1f}; load-control target = {lam_T:.1f} "
          f"(above peak => no equilibrium exists)\n")
    print(f"{'scenario':<36}{'min strain':>12}{'peak load':>11}  note")
    print("-" * 82)
    print(f"{'(a) LoadControl + Newton':<36}{a_eps:>12.5f}{a_peak:>11.3f}  "
          f"{'reached softening' if a_soft else 'DIED at peak (no state beyond)'}")
    print(f"{'(b) adaptive halving / LoadControl':<36}{b_eps:>12.5f}{b_lam:>11.3f}  "
          f"{'reached softening' if (b_eps < 1.5*EPSC0) else 'STILL DIES (no far branch to find)'}")
    print(f"{'(c) DisplacementControl [reference]':<36}{c_eps:>12.5f}"
          f"{(c_lam or 0):>11.3f}  {'traced softening branch' if c_soft else 'no'}")
    print(f"{'(d) robust driver (load->disp)':<36}{d_eps:>12.5f}{'-':>11}  "
          f"{'traced softening branch' if d_soft else 'no'} ({d_sw} switch)")
    print("-" * 82)
    gap = not a_soft and (b_eps > 1.5 * EPSC0)            # load control cannot soften
    fix = d_soft and c_soft                              # disp + robust trace softening
    ok = gap and fix
    print("\nVERDICT:", "PASS - load control (with the full algorithm ladder + "
          "adaptive halving) cannot pass the strength peak; the rung-3 constraint "
          "switch traces the softening branch."
          if ok else "FAIL - inspect above")
    raise SystemExit(0 if ok else 1)
