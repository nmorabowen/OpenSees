"""Torture fixture #1 - von Mises shallow-truss snap-through (rung-3 motivator).

A 2-bar shallow truss snaps through as the apex is pushed past the flat line.
At the limit load the tangent is singular, so LOAD control (and step-halving on
load control) cannot pass it - only a constraint SWITCH (load -> displacement /
arc-length, the driver's rung 3) carries the full equilibrium path.

Success metric (uniform across scenarios): did the method drive the apex PAST
the flat position (min apex_y < 0, i.e. through the snap) with the load
reversing sign? Load control cannot; the rung-3 switch can.

Needs NO energy recorder, so it runs on today's build. Run with the dist build:
    python Ladruno_scripts/robust_solve_tests/torture_snapthrough.py
"""
import os
import sys

DIST = os.environ.get("LADRUNO_DIST",
                      r"C:\Users\nmb\Documents\Github\OpenSees\dist\bin")
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)
    sys.path.insert(0, DIST)
import opensees as ops  # noqa: E402

E, A = 1.0e4, 1.0
a, b = 1.0, 0.10                    # half-span, rise (shallow => snap-through)
CTRL_NODE, CTRL_DOF = 3, 2          # apex vertical DOF (the single free DOF)
APEX_Y0 = b
Y_TARGET = -0.20                    # push apex from +0.10 to -0.20 (through snap)


def build():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 2 * a, 0.0)
    ops.node(3, a, b)
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    ops.fix(3, 1, 0)               # symmetry: lock apex X -> single DOF (Y)
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.element("corotTruss", 1, 1, 3, A, 1)
    ops.element("corotTruss", 2, 2, 3, A, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 0.0, -1.0)         # downward unit reference load
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-8, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def apex_y():
    return APEX_Y0 + ops.nodeDisp(CTRL_NODE, CTRL_DOF)


# --- displacement control: the reference path-follower (passes limit pts) -----
def run_dispcontrol(du, max_steps):
    build()
    ops.integrator("DisplacementControl", CTRL_NODE, CTRL_DOF, du)
    done = 0
    lam_lim = 0.0          # snap-through limit load = max lam while apex still above flat
    reversed_ = False
    min_y = APEX_Y0
    for _ in range(max_steps):
        if ops.analyze(1) != 0:
            break
        done += 1
        lam = ops.getLoadFactor(1)
        y = apex_y()
        min_y = min(min_y, y)
        if y >= 0.0:
            lam_lim = max(lam_lim, lam)
        if lam < -1e-9:
            reversed_ = True
        if y <= Y_TARGET:
            break
    return done, min_y, lam_lim, reversed_


# --- load control to a target ABOVE the limit load: must die at the limit -----
def run_loadcontrol(dlam, lam_target, max_steps):
    build()
    ops.integrator("LoadControl", dlam)
    done = 0
    peak = 0.0
    min_y = APEX_Y0
    for _ in range(max_steps):
        if ops.analyze(1) != 0:
            break
        done += 1
        peak = max(peak, ops.getLoadFactor(1))
        min_y = min(min_y, apex_y())
        if ops.getLoadFactor(1) >= lam_target - 1e-9:
            break
    return done, min_y, peak


# --- adaptive step-halving ON load control (still cannot pass a limit point) --
def run_adaptive_loadcontrol(dlam, lam_target):
    from ladruno_solve import adaptive_static
    build()
    n = max(1, round(lam_target / dlam))
    res = adaptive_static(
        ops, lambda s: ops.integrator("LoadControl", dlam * s), n_steps=n)
    return res.completed / n, apex_y(), ops.getLoadFactor(1), bool(res)


# --- the robust driver: load control, rung-3 switch to disp control on stall --
def run_robust(dlam, du, max_steps):
    build()
    ops.integrator("LoadControl", dlam)
    mode = "LoadControl"
    done = 0
    switches = 0
    reversed_ = False
    min_y = APEX_Y0
    for _ in range(max_steps):
        if ops.analyze(1) == 0:
            done += 1
            min_y = min(min_y, apex_y())
            if ops.getLoadFactor(1) < -1e-9:
                reversed_ = True
            if apex_y() <= Y_TARGET:
                break
            continue
        if mode == "LoadControl":                 # rung-3 constraint switch
            ops.integrator("DisplacementControl", CTRL_NODE, CTRL_DOF, du)
            mode = "DisplacementControl"
            switches += 1
            continue
        break
    return done, min_y, reversed_, switches


def crossed(min_y):
    return min_y < -0.05            # apex driven well past the flat line


if __name__ == "__main__":
    DU = -1.0e-3
    MAXS = 600

    # 1) calibrate the snap-through limit load from a displacement-control trace
    c_done, c_miny, lam_lim, c_rev = run_dispcontrol(DU, MAXS)

    # 2) set a load target 50% ABOVE the limit so load control MUST fail at it
    lam_T = 1.5 * lam_lim
    dlam = lam_T / 300.0

    a_done, a_miny, a_peak = run_loadcontrol(dlam, lam_T, 1200)
    b_frac, b_y, b_peak, b_ok = run_adaptive_loadcontrol(dlam, lam_T)
    d_done, d_miny, d_rev, d_sw = run_robust(dlam, DU, MAXS)

    print(f"snap-through limit load (calibrated): lam_lim = {lam_lim:.4f}; "
          f"load-control target lam_T = {lam_T:.4f}\n")
    print(f"{'scenario':<36}{'min apex_y':>11}{'peak lam':>11}  note")
    print("-" * 82)
    print(f"{'(a) LoadControl + Newton':<36}{a_miny:>11.4f}{a_peak:>11.3f}  "
          f"{'crossed' if crossed(a_miny) else 'STUCK above flat (died at limit)'}")
    print(f"{'(b) adaptive halving / LoadControl':<36}{b_y:>11.4f}{b_peak:>11.3f}  "
          f"{'crossed' if crossed(b_y) else 'STUCK (halving cannot pass a singular tangent)'}")
    print(f"{'(c) DisplacementControl [reference]':<36}{c_miny:>11.4f}{lam_lim:>11.3f}  "
          f"{'crossed + load reversed' if crossed(c_miny) and c_rev else 'no'}")
    print(f"{'(d) robust driver (load->disp)':<36}{d_miny:>11.4f}{'-':>11}  "
          f"{'crossed + load reversed' if crossed(d_miny) and d_rev else 'no'} "
          f"({d_sw} switch)")
    print("-" * 82)
    gap = not crossed(a_miny) and not crossed(b_y)        # load control cannot cross
    fix = crossed(d_miny) and d_rev                       # robust crosses + snaps
    ref = crossed(c_miny) and c_rev
    ok = gap and fix and ref
    print("\nVERDICT:", "PASS - load control (even with adaptive halving + the "
          "algorithm ladder) cannot pass the limit point; the rung-3 constraint "
          "switch carries the full snap-through path."
          if ok else "FAIL - inspect above")
    raise SystemExit(0 if ok else 1)
