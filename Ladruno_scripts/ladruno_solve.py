"""Ladruno blessed adaptive-stepping driver (v1b of the RC-3D plan).

Both RC-3D validation gates (Ladruno_implementation/21_rc3d_validation_gates.md)
needed step-halving + a non-default algorithm to converge through softening; the
adversarial reviews of ADR 20 agreed that adaptive stepping belongs in a *tested
Python wrapper*, NOT a new C++ integrator (halve-on-fail is an outer loop over
`analyze()`; only auto-arc-length-radius control would justify C++).

This module is that wrapper. It is integrator-agnostic: you pass a `set_step(scale)`
callback that (re)issues `ops.integrator(...)` at `scale * nominal_increment`, and the
driver advances a target number of nominal steps, halving the increment and walking an
algorithm ladder on non-convergence, growing back on success.

    from ladruno_solve import adaptive_static
    ops.integrator("DisplacementControl", ctrl, dof, du)          # nominal
    ok = adaptive_static(ops, lambda s: ops.integrator(
                             "DisplacementControl", ctrl, dof, du * s),
                         n_steps=400)

Returns an AdaptiveResult; `bool(result)` is True iff the full target was reached.
Pure stdlib + an `ops` handle (local build or openseespy). No OpenSees source change.
"""
from __future__ import annotations


class AdaptiveResult:
    __slots__ = ("completed", "target", "substeps", "min_scale_used",
                 "algo_used", "failed_at")

    def __init__(self):
        self.completed = 0.0          # nominal steps actually advanced
        self.target = 0.0
        self.substeps = 0             # total ops.analyze(1) calls
        self.min_scale_used = 1.0     # smallest scale that was needed
        self.algo_used = {}           # {algorithm_name: success_count}
        self.failed_at = None         # nominal-progress where it gave up (or None)

    def __bool__(self):
        return self.failed_at is None and self.completed >= self.target - 1e-9

    def __repr__(self):
        s = "OK" if self else f"INCOMPLETE@{self.failed_at:.3f}"
        return (f"<AdaptiveResult {s} {self.completed:.3f}/{self.target:.0f} "
                f"substeps={self.substeps} min_scale={self.min_scale_used:g} "
                f"algos={self.algo_used}>")


# (name, *args) passed straight to ops.algorithm(...). Ordered easy->robust.
_DEFAULT_LADDER = (
    ("KrylovNewton",),
    ("Newton",),
    ("NewtonLineSearch",),
    ("ModifiedNewton",),
)


def adaptive_static(ops, set_step, n_steps, *,
                    algorithms=_DEFAULT_LADDER,
                    min_scale=1.0 / 1024,
                    grow=2.0,
                    on_step=None,
                    verbose=False):
    """Advance `n_steps` nominal increments with adaptive sub-stepping.

    ops        : the OpenSees module handle (opensees / openseespy.opensees).
    set_step   : callable(scale) -> re-issues ops.integrator(...) at scale*nominal.
                 Called before each attempt; scale in (min_scale, 1].
    n_steps    : number of NOMINAL increments to advance (the analysis covers the
                 same total as `ops.analyze(n_steps)` at the nominal increment).
    algorithms : ladder of ops.algorithm(*spec) tuples, tried in order on failure.
    min_scale  : give up once the increment must shrink below this (returns INCOMPLETE).
    grow       : after a clean success, multiply scale by this (capped at 1.0).
    on_step    : optional callable(progress_fraction) called after each committed
                 sub-step (use to record history).
    verbose    : print a line per scale change / algorithm switch.

    Returns an AdaptiveResult (truthy iff the full target was reached).
    """
    res = AdaptiveResult()
    res.target = float(n_steps)
    done = 0.0
    scale = 1.0
    algo_i = 0

    def _set_algo(i):
        spec = algorithms[i % len(algorithms)]
        ops.algorithm(*spec)
        return spec[0]

    algo_name = _set_algo(algo_i)

    while done < res.target - 1e-9:
        # don't overshoot the target on the last sub-step
        scale = min(scale, res.target - done)
        set_step(scale)
        res.substeps += 1
        ok = (ops.analyze(1) == 0)
        if ok:
            done += scale
            res.completed = done
            res.min_scale_used = min(res.min_scale_used, scale)
            res.algo_used[algo_name] = res.algo_used.get(algo_name, 0) + 1
            if on_step is not None:
                on_step(min(1.0, done / res.target))
            # recover toward nominal, and fall back to the easy algorithm
            if scale < 1.0:
                scale = min(1.0, scale * grow)
            if algo_i != 0:
                algo_i = 0
                algo_name = _set_algo(algo_i)
            continue
        # failure: first try the next algorithm at the SAME scale, then halve
        if algo_i < len(algorithms) - 1:
            algo_i += 1
            algo_name = _set_algo(algo_i)
            if verbose:
                print(f"  [adaptive] fail @scale={scale:g} -> algorithm {algo_name}")
            continue
        # exhausted the ladder at this scale: halve and reset to the easy algorithm
        scale *= 0.5
        algo_i = 0
        algo_name = _set_algo(algo_i)
        if verbose:
            print(f"  [adaptive] fail -> halve to scale={scale:g}")
        if scale < min_scale:
            res.failed_at = done
            return res
    return res


# --------------------------------------------------------------------------
# self-test: the Gate-2 confined column at p/fc=0.20 needs KrylovNewton +
# step-halving. Show the wrapper drives it to completion starting from the
# WORST default (plain Newton, fixed step), where a naive loop diverges.
# --------------------------------------------------------------------------
def _selftest():
    import os
    import sys
    try:
        import opensees as ops
    except ModuleNotFoundError:
        dist = os.environ.get(
            "LADRUNO_DIST",
            r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin")
        if os.path.isdir(dist):
            os.add_dll_directory(dist)
            sys.path.insert(0, dist)
        import opensees as ops

    E, NU, FC, FT = 30000.0, 0.2, 30.0, 3.0
    CE = [-0.4 * FC / E, -0.002, -0.005, -0.015]
    CS = [-0.4 * FC, -FC, -0.5 * FC, -0.2 * FC]
    CD = [1 - abs(CS[i]) / (E * abs(CE[i])) for i in range(4)]
    TE, TS, TD = [FT / E, 1e-3], [FT, 0.3], [0.0, 1 - 0.3 / (E * 1e-3)]
    C = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
    p = 6.0

    def build():
        ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
        for t, c in C.items():
            ops.node(t, *c)
        ops.nDMaterial("ASDConcrete3D", 1, E, NU, "-Te", *TE, "-Ts", *TS,
                       "-Td", *TD, "-Ce", *CE, "-Cs", *CS, "-Cd", *CD)
        ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
        for n in (1, 4, 5, 8): ops.fix(n, 1, 0, 0)
        for n in (1, 2, 5, 6): ops.fix(n, 0, 1, 0)
        for n in (1, 2, 3, 4): ops.fix(n, 0, 0, 1)
        ops.equalDOF(2, 3, 1); ops.equalDOF(2, 6, 1); ops.equalDOF(2, 7, 1)
        ops.equalDOF(4, 3, 2); ops.equalDOF(4, 7, 2); ops.equalDOF(4, 8, 2)
        ops.equalDOF(5, 6, 3); ops.equalDOF(5, 7, 3); ops.equalDOF(5, 8, 3)
        ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
        for n in (2, 3, 6, 7): ops.load(n, -p / 4.0, 0, 0)
        for n in (3, 4, 7, 8): ops.load(n, 0, -p / 4.0, 0)
        ops.constraints("Transformation"); ops.numberer("RCM")
        ops.system("FullGeneral"); ops.test("NormDispIncr", 1e-6, 50, 0)
        ops.integrator("LoadControl", 1.0 / 100); ops.algorithm("KrylovNewton")
        ops.analysis("Static"); ops.analyze(100); ops.loadConst("-time", 0.0)
        ops.timeSeries("Linear", 2); ops.pattern("Plain", 2, 2)
        for n in (5, 6, 7, 8): ops.load(n, 0, 0, -0.25)

    nstep, umax = 600, -12e-3
    du = umax / nstep

    # (a) naive: plain Newton, fixed step, no adaptation
    build(); ops.algorithm("Newton"); ops.integrator("DisplacementControl", 5, 3, du)
    naive = 0
    for _ in range(nstep):
        if ops.analyze(1) != 0:
            break
        naive += 1

    # (b) the wrapper from the same worst-case start
    build(); ops.algorithm("Newton")
    res = adaptive_static(
        ops, lambda s: ops.integrator("DisplacementControl", 5, 3, du * s),
        n_steps=nstep, verbose=False)

    print(f"naive Newton fixed-step : reached {naive}/{nstep} steps "
          f"({'DIVERGED' if naive < nstep else 'ok'})")
    print(f"adaptive_static         : {res}")
    ok = (naive < nstep) and bool(res)
    print("\nSELF-TEST:", "PASS — wrapper completes where naive diverges" if ok
          else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(_selftest())
