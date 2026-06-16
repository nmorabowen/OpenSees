"""Ladruno robust nonlinear-solution driver - Layer 0 (ADR 31).

The increment-sizing SPINE (Abaqus-faithful automatic incrementation) plus the
rung-3 constraint switch, in pure Python over ops.analyze(1). No C++ change.

Rungs implemented at Layer 0:
  0  increment cutback / grow          (target-iteration sizing, via adaptive_static)
  1  line search                       (in the algorithm ladder)
  2  modified / Krylov / BFGS Newton   (algorithm ladder)
  3  CONSTRAINT SWITCH  load -> displacement control on a control DOF
  5  dynamics fall-through             (HOOK only - refused until handoff primitive ships)

Rung 4 (auto-stabilization) is deliberately REFUSED here: per ADR-31 R-OBS the
dissipation hard-gate is a no-op until the Layer-1.5 integrator getters land, so
the driver must fail loud rather than diffuse silently. `stabilize=True` raises.

Decision log: every cutback / algorithm switch / constraint switch / verdict is
appended as one JSON object per line (JSONL) so a run is auditable (R-LOG-MASK).

    from robust_drive import robust_drive
    res = robust_drive(
        ops,
        done=lambda: ops.nodeDisp(2, 1) <= -0.010,     # termination predicate
        load_increment=0.2,                            # nominal LoadControl step
        control=(2, 1, -5e-5),                         # (node, dof, du) for rung-3
        log_path="run.jsonl")
    assert res                                         # truthy iff target reached
"""
from __future__ import annotations

import json

from ladruno_solve import _DEFAULT_LADDER


class RobustResult:
    __slots__ = ("completed", "substeps", "switches", "mode", "verdict",
                 "min_scale_used", "algos")

    def __init__(self):
        self.completed = False
        self.substeps = 0
        self.switches = 0
        self.mode = "LoadControl"
        self.verdict = "incomplete"      # equilibrium | incomplete | aborted | refused
        self.min_scale_used = 1.0
        self.algos = {}

    def __bool__(self):
        return self.completed and self.verdict == "equilibrium"

    def __repr__(self):
        return (f"<RobustResult {self.verdict} completed={self.completed} "
                f"substeps={self.substeps} switches={self.switches} "
                f"mode={self.mode} min_scale={self.min_scale_used:g}>")


def robust_drive(ops, done, *,
                 load_increment,
                 control=None,
                 algorithms=_DEFAULT_LADDER,
                 min_scale=1.0 / 1024,
                 grow=2.0,
                 max_substeps=20000,
                 stabilize=False,
                 dynamics=False,
                 log_path=None,
                 verbose=False):
    """Drive a static analysis to `done()` with the rung-0..3 ladder.

    done           : callable() -> bool, True when the analysis target is reached.
    load_increment : nominal LoadControl load-factor step (scaled by the spine).
    control        : (node, dof, du) enabling the rung-3 switch to DisplacementControl;
                     None disables the switch (load control only).
    max_substeps   : global budget (R-THRASH/R-TERMINATION); INCOMPLETE if exceeded.
    stabilize      : rung-4 auto-stabilization - REFUSED at Layer 0 (raises).
    dynamics       : rung-5 fall-through - HOOK only; logs and aborts (not built).
    """
    if stabilize:
        raise NotImplementedError(
            "rung-4 auto-stabilization is refused at Layer 0 (ADR-31 R-OBS): the "
            "dissipation hard-gate is a no-op until the Layer-1.5 LadrunoArcLength "
            "getters ship. Fail loud, not diffuse silent.")

    res = RobustResult()
    _log_fh = open(log_path, "w", encoding="utf-8") if log_path else None

    def log(event, **kw):
        if _log_fh is not None:
            rec = {"event": event, "substep": res.substeps, "mode": res.mode}
            rec.update(kw)
            _log_fh.write(json.dumps(rec) + "\n")
            _log_fh.flush()
        if verbose:
            print(f"  [robust] {event} {kw}")

    scale = 1.0
    algo_i = 0

    def set_algo(i):
        spec = algorithms[i % len(algorithms)]
        ops.algorithm(*spec)
        return spec[0]

    def set_step():
        if res.mode == "LoadControl":
            ops.integrator("LoadControl", load_increment * scale)
        else:                                       # DisplacementControl
            node, dof, du = control
            ops.integrator("DisplacementControl", node, dof, du * scale)

    algo_name = set_algo(algo_i)
    log("start", load_increment=load_increment, has_control=control is not None)

    try:
        while not done():
            if res.substeps >= max_substeps:
                res.verdict = "incomplete"
                log("budget_exhausted", max_substeps=max_substeps)
                return res
            set_step()
            res.substeps += 1
            if ops.analyze(1) == 0:                 # converged
                res.min_scale_used = min(res.min_scale_used, scale)
                res.algos[algo_name] = res.algos.get(algo_name, 0) + 1
                if scale < 1.0:                     # grow back toward nominal
                    scale = min(1.0, scale * grow)
                if algo_i != 0:                     # fall back to the easy algorithm
                    algo_i = 0
                    algo_name = set_algo(algo_i)
                continue
            # --- failure: rung 1-2 (walk the algorithm ladder at this scale) ---
            if algo_i < len(algorithms) - 1:
                algo_i += 1
                algo_name = set_algo(algo_i)
                log("algorithm_switch", to=algo_name, scale=scale)
                continue
            # --- rung 0 (halve), reset to the easy algorithm -------------------
            scale *= 0.5
            algo_i = 0
            algo_name = set_algo(algo_i)
            if scale >= min_scale:
                log("cutback", scale=scale)
                continue
            # --- rung 3 (constraint switch load -> displacement) ---------------
            if res.mode == "LoadControl" and control is not None:
                res.mode = "DisplacementControl"
                res.switches += 1
                scale = 1.0
                algo_i = 0
                algo_name = set_algo(algo_i)
                log("constraint_switch", to="DisplacementControl",
                    reason="load-control exhausted at limit/peak (rung 3)")
                continue
            # --- rung 5 (dynamics) - hook only --------------------------------
            if dynamics:
                res.verdict = "aborted"
                log("dynamics_requested",
                    note="rung-5 fall-through not built at Layer 0 (handoff "
                         "primitive pending, ADR-31 R-HANDOFF)")
                return res
            res.verdict = "incomplete"
            log("give_up", scale=scale)
            return res

        res.completed = True
        res.verdict = "equilibrium"
        log("done", verdict=res.verdict)
        return res
    finally:
        if _log_fh is not None:
            _log_fh.close()


# --------------------------------------------------------------------------
# self-test: a single softening Concrete02 truss. Load control cannot pass the
# strength peak (no equilibrium beyond it); the rung-3 switch traces the
# softening branch. Mirrors robust_solve_tests/torture_softening.py.
# --------------------------------------------------------------------------
def _selftest():
    import os
    import sys
    DIST = os.environ.get("LADRUNO_DIST",
                          r"C:\Users\nmb\Documents\Github\OpenSees\dist\bin")
    if os.path.isdir(DIST):
        os.add_dll_directory(DIST)
        sys.path.insert(0, DIST)
    import opensees as ops

    FPC, EPSC0, FPCU, EPSU = -30.0, -0.002, -6.0, -0.012
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("Concrete02", 1, FPC, EPSC0, FPCU, EPSU, 0.1, 3.0, 300.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, -1.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGen")
    ops.test("NormDispIncr", 1e-8, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")

    res = robust_drive(
        ops,
        done=lambda: ops.nodeDisp(2, 1) <= -0.010,
        load_increment=39.0 / 200.0,
        control=(2, 1, -5.0e-5),
        log_path=None,
        verbose=True)

    eps = ops.nodeDisp(2, 1)
    softened = eps <= -0.003 and ops.getLoadFactor(1) < 0.95 * 30.0
    ok = bool(res) and res.switches == 1 and softened
    print(f"\n{res}")
    print(f"final strain={eps:.5f}  load={ops.getLoadFactor(1):.3f}")
    print("SELF-TEST:", "PASS - rung-3 switch traced the softening branch"
          if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(_selftest())
