"""ADR-65 Route B P0 headroom oracle (methodology dry-run).

Measures how much an adaptive-dt driver could save by logging the tangent-tracked
criticalTimeStep() per step (integrator CentralDifferenceLadruno -recompute 1) on
two 1-D truss-bar cases where dt_cr evolves:

  A) Concrete01 damage-softening bar, prescribed tip displacement: ramp into the
     softening branch, unload, then hold/wiggle at degraded stiffness. Probes
     (1) whether -recompute 1 actually tracks the tangent, (2) what the pencil
     returns on the NEGATIVE-tangent branch, (3) post-damage headroom.
  B) Steel01 (b=0.02) yield pulse: tangent dt_cr should rise ~sqrt(1/b) during
     plastic flow and snap back on elastic unload. Probes the driver's
     shrink-back requirement + the elastic-unload safety caveat (running at the
     plastic-tangent dt is UNSAFE against elastic unloading waves).

Adaptive-step estimate: N_adapt = sum_i dt_fixed / min(s*dtcr_i, DT_GROW_CAP...)
vs N_fixed actual steps, same physical duration. safety s = 0.9.
"""
import os
import sys
import math

DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
sys.path.insert(0, DIST)
os.add_dll_directory(DIST)
import opensees as ops  # noqa: E402

assert DIST.lower() in ops.__file__.lower(), "wrong opensees loaded: %s" % ops.__file__
print("[env] opensees =", ops.__file__)

CDL = "CentralDifferenceLadruno"
SAFETY = 0.9


def _bar(mat_cmd, N=10, L_e=1.0, A=1.0, rho=1.0):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    for i in range(N + 1):
        ops.node(i + 1, i * L_e)
    ops.fix(1, 1)
    mat_cmd()
    for e in range(1, N + 1):
        ops.element("Truss", e, e, e + 1, A, 1, "-rho", rho)
    return N + 1  # tip node


def _analysis():
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-recompute", 1)
    ops.analysis("Transient")


def _run(tip, dt, nsteps, label):
    trace = []
    bad = 0
    for i in range(nsteps):
        ok = ops.analyze(1, dt)
        if ok != 0:
            print("  [%s] analyze FAILED at step %d" % (label, i))
            break
        d = ops.criticalTimeStep()
        if d is None or d <= 0.0 or d != d:
            bad += 1
            d = float("nan")
        trace.append((ops.getTime(), d, ops.nodeDisp(tip, 1)))
    return trace, bad


def _report(label, trace, dt_fixed, dtcr_elastic, bad):
    ds = [d for (_, d, _) in trace if d == d]  # drop NaN
    if not ds:
        print("  [%s] no valid dt_cr samples!" % label)
        return
    d0, dmin, dmax, dend = ds[0], min(ds), max(ds), ds[-1]
    # TWO adaptive-step estimates over the same physical duration:
    #  (i) UNSAFE  — grow to s*dtcr_i using the instantaneous (tangent) dt_cr.
    #      This is what a naive driver that trusts criticalTimeStep() would do.
    #  (ii) SAFE   — clamp growth to the UNDAMAGED/elastic dt_cr, because a
    #      reloadable material (plasticity, damage-with-reload) can stiffen back
    #      to E0 in one step and a wave hitting it reloads elastic -> the only
    #      dt safe against reload is the elastic one.
    # Honest baseline = a careful user's CONSTANT safe step = SAFETY*dtcr_elastic
    # (they know the material reloads elastic, so they pick the undamaged limit).
    # Compare adaptive step-counts against THAT, at the SAME safety factor, so the
    # speedup isolates nonlinear headroom (not the safety-factor choice).
    duration = sum(dt_fixed for _ in trace)
    n_base = duration / (SAFETY * dtcr_elastic)
    n_unsafe = n_safe = 0.0
    last = d0
    for (_, d, _) in trace:
        if d == d:
            last = d
        n_unsafe += dt_fixed / (SAFETY * last)                     # trust tangent
        n_safe += dt_fixed / (SAFETY * min(last, dtcr_elastic))    # clamp to elastic
    n_fixed = len(trace)
    print("  [%s] steps=%d  dt_run=%.3e  (equal-safety baseline = %.0f steps)"
          % (label, n_fixed, dt_fixed, n_base))
    print("        dt_cr: t0=%.3e  min=%.3e  max=%.3e  end=%.3e   (elastic ref %.3e)"
          % (d0, dmin, dmax, dend, dtcr_elastic))
    print("        invalid/nonpositive dt_cr samples: %d / %d" % (bad, n_fixed))
    print("        max/elastic ratio = %.2f   (tangent OVER-reports safe headroom by "
          "this factor)" % (dmax / dtcr_elastic))
    print("        UNSAFE grow-to-tangent: %.0f steps -> %.2fx vs baseline  <-- MIRAGE "
          "(blows up on elastic reload)" % (n_unsafe, n_base / max(n_unsafe, 1.0)))
    print("        SAFE  grow-clamped-elastic: %.0f steps -> %.2fx vs baseline  <-- real "
          "growth headroom" % (n_safe, n_base / max(n_safe, 1.0)))


# ---------------------------------------------------------------- Case A
def case_a():
    print("\n=== Case A: Concrete01 damage-softening bar, prescribed tip disp ===")
    fpc, epsc0, fpcu, epsU = -30.0, -0.002, -6.0, -0.01
    E0 = 2.0 * fpc / epsc0                    # 30000
    rho = 1.0
    c = math.sqrt(E0 / rho)
    L_e = 1.0
    N = 10
    tip = _bar(lambda: ops.uniaxialMaterial("Concrete01", 1, fpc, epsc0, fpcu, epsU),
               N=N, L_e=L_e, rho=rho)
    dtcr_el = L_e / c
    dt = 0.5 * dtcr_el

    # prescribed tip displacement: ramp to -0.04 (avg strain -0.004, past peak),
    # unload to -0.025, hold. Slow vs wave time so it's quasi-static-ish.
    L = N * L_e
    T1, T2, T3 = 400 * dt, 200 * dt, 400 * dt
    u1, u2 = -0.004 * L, -0.0025 * L
    ops.timeSeries("Path", 1,
                   "-time", 0.0, T1, T1 + T2, T1 + T2 + T3,
                   "-values", 0.0, u1, u2, u2)
    ops.pattern("Plain", 1, 1)
    ops.sp(tip, 1, 1.0)
    _analysis()
    nsteps = int((T1 + T2 + T3) / dt)
    trace, bad = _run(tip, dt, nsteps, "A")
    _report("A", trace, dt, dtcr_el, bad)
    # phase medians to see the shape
    third = len(trace) // 3
    for name, seg in (("ramp/soften", trace[:third]),
                      ("unload", trace[third:2 * third]),
                      ("hold(damaged)", trace[2 * third:])):
        ds = sorted(d for (_, d, _) in seg if d == d)
        if ds:
            print("        phase %-14s median dt_cr = %.3e  (%.2fx elastic)"
                  % (name, ds[len(ds) // 2], ds[len(ds) // 2] / dtcr_el))


# ---------------------------------------------------------------- Case B
def case_b():
    print("\n=== Case B: Steel01 (b=0.02) yield pulse ===")
    Fy, E0, b = 300.0, 30000.0, 0.02
    rho = 1.0
    c = math.sqrt(E0 / rho)
    L_e = 1.0
    N = 10
    tip = _bar(lambda: ops.uniaxialMaterial("Steel01", 1, Fy, E0, b),
               N=N, L_e=L_e, rho=rho)
    dtcr_el = L_e / c
    dt = 0.5 * dtcr_el

    # tip force pulse: ramp to 1.2*Fy*A (yields whole bar), back to 0, hold.
    P = 1.2 * Fy
    T1, T2, T3 = 400 * dt, 400 * dt, 400 * dt
    ops.timeSeries("Path", 1,
                   "-time", 0.0, T1, T1 + T2, T1 + T2 + T3,
                   "-values", 0.0, P, 0.0, 0.0)
    ops.pattern("Plain", 1, 1)
    ops.load(tip, 1.0)
    _analysis()
    nsteps = int((T1 + T2 + T3) / dt)
    trace, bad = _run(tip, dt, nsteps, "B")
    _report("B", trace, dt, dtcr_el, bad)
    third = len(trace) // 3
    for name, seg in (("load->yield", trace[:third]),
                      ("unload", trace[third:2 * third]),
                      ("hold(elastic)", trace[2 * third:])):
        ds = sorted(d for (_, d, _) in seg if d == d)
        if ds:
            print("        phase %-14s median dt_cr = %.3e  (%.2fx elastic)"
                  % (name, ds[len(ds) // 2], ds[len(ds) // 2] / dtcr_el))
    print("        sqrt(1/b) = %.2f  (expected plastic-phase dt_cr inflation)"
          % math.sqrt(1.0 / b))


if __name__ == "__main__":
    case_a()
    case_b()
    print("\n[done]")
