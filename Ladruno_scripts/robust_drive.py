"""Ladruno robust nonlinear-solution driver (ADR 31).

An increment-sizing SPINE (Abaqus-faithful automatic incrementation) plus an
escalation ladder, in pure Python over ops.analyze(1). The hard numerics are the
already-shipped fork C++ (LadrunoArcLength -stabilize 33004, LadrunoStabilized-
Unbalance, LadrunoDynamicRelaxation 33005); this driver is the orchestration +
observability that sequences them honestly.

Rungs (driver numbering):
  0  increment cutback / grow            (target-iteration sizing)
  1  line search                         (algorithm ladder)
  2  modified / Krylov / BFGS Newton     (algorithm ladder)
  3  CONSTRAINT SWITCH  load -> displacement control on a control DOF
  4  AUTO-STABILIZE  LadrunoArcLength -stabilize + LadrunoStabilizedUnbalance,
     dissipation gate read each step, scaleCVisc ramp-down (R-RAMPDOWN)
  5  DYNAMICS fall-through  LadrunoDynamicRelaxation, atomic implicit<->DR
     handoff (R-HANDOFF), force-residual settling gate, indirect-control polish

Escalation order at the cutback floor: rung-3 (if a control DOF is given) ->
rung-4 (if `stabilize`) -> rung-5 (if `dynamics`) -> give up.

HONEST VERDICTS (ADR-31 R-OBS / R-DIFFUSION / R-SNAPBACK):
  * a rung-0..3 success is `equilibrium` (publishable, truthy).
  * a rung-4 stabilized success is `regularized` at best and `unverified` unless
    the c-reduction diffusion bound was computed -- it is NEVER `equilibrium`
    (viscous f_v pollutes the path; LadrunoStabilizedUnbalance proves point-wise
    equilibrium but is not a fidelity certificate).
  * a rung-5 dynamics rest state is `regularized` (a relaxed rest, not a traced
    branch) unless the polish tail re-lands on the true branch.
`bool(res)` is True ONLY for `equilibrium`, so a regularized/unverified result is
falsy by construction -- a paper-grade consumer cannot mistake it for clean.

    from robust_drive import robust_drive
    res = robust_drive(ops,
        done=lambda: ops.nodeDisp(2, 1) <= -0.010,
        load_increment=0.2,
        control=(2, 1, -5e-5),          # enables rung-3
        stabilize=True,                 # enables rung-4
        dynamics=True,                  # enables rung-5 (needs the ladrunoDR build)
        log_path="run.jsonl")
"""
from __future__ import annotations

import json

from ladruno_solve import _DEFAULT_LADDER


class RobustResult:
    __slots__ = ("completed", "substeps", "switches", "mode", "verdict",
                 "min_scale_used", "algos", "dissipation_ratio", "stab_dissipated",
                 "peak_drift", "dr_settled", "over_diffused", "dr_lambda")

    def __init__(self):
        self.completed = False
        self.substeps = 0
        self.switches = 0
        self.mode = "LoadControl"
        # equilibrium | regularized | unverified | incomplete | aborted | refused
        self.verdict = "incomplete"
        self.min_scale_used = 1.0
        self.algos = {}
        self.dissipation_ratio = 0.0    # W_stab / Estrain0 at end of a rung-4 phase
        self.stab_dissipated = 0.0      # W_stab (the L0 energy-partition SW term)
        self.peak_drift = None          # c-reduction diffusion bound (None = not run)
        self.dr_settled = False         # did the rung-5 DR phase reach quasi-static
        self.over_diffused = False      # rung-4 dissipation crossed the hard gate
        self.dr_lambda = None           # load factor frozen across the DR handoff (R-HANDOFF)

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
                 peak_cutbacks=3,
                 max_substeps=20000,
                 stabilize=False,
                 stab_f=1.0e-3,
                 stab_tol=1.0e-8,
                 stab_hard_gate=0.05,
                 stab_rampdown_window=8,
                 stab_max_cutbacks=12,
                 dynamics=False,
                 dr_settle_tol=1.0e-4,
                 dr_max_steps=4000,
                 dr_pref=None,
                 dr_setup=None,
                 log_path=None,
                 verbose=False):
    """Drive a static analysis to `done()` with the rung-0..5 ladder.

    load_increment : nominal LoadControl load-factor step (scaled by the spine).
    control        : (node, dof, du) enabling rung-3 (load->disp switch); None disables.
    stabilize      : enable rung-4 auto-stabilization (the Layer-1.5 getters must be
                     built). A stabilized success is `regularized`/`unverified`, never
                     `equilibrium`.
    stab_f         : -stabilize dissipated-energy fraction. ELEVATED vs the Abaqus
                     2e-4 default because the default is too weak to cross a hard
                     limit (measured); -adaptStab is deliberately NOT used (it
                     re-weakens c and prevents crossing).
    stab_hard_gate : dissipation ratio above which the run is flagged over-diffused
                     (informs the verdict; does NOT abort -- crossing inherently
                     dissipates).
    stab_rampdown_window : clean stabilized steps before scaleCVisc(0.5) decays c
                     (R-RAMPDOWN).
    dynamics       : enable rung-5 DR fall-through (needs the `ladrunoDR` command).
    dr_settle_tol  : DR is quasi-static once residualNorm < dr_settle_tol * Pref
                     (the mass-free force criterion; R-DR-ENERGY -- the physical-mass
                     KE is ~0 on DR's pseudo-mass models, so the gate is force-based).
    dr_pref        : reference load norm ||P|| for the settle gate. None -> proxy
                     max(1, |lambda|) (exact when the reference-load norm is ~1).
    dr_setup       : optional callable(ops) that configures the transient analysis
                     for the DR phase (constraints/numberer/system). If None, the
                     driver uses Plain/RCM/BandGen defaults (fine for simple models).
    """
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
    cutbacks = 0                                    # consecutive cutbacks since last commit

    def set_algo(i):
        spec = algorithms[i % len(algorithms)]
        ops.algorithm(*spec)
        return spec[0]

    def can_switch():
        return res.mode == "LoadControl" and control is not None

    def set_step():
        # LoadControl / DisplacementControl are stateless -> re-issued each step at
        # the current scale. The stabilized integrator is STATEFUL (watchdog +
        # calibrated cVisc), so its phase issues it ONCE and never re-issues here.
        if res.mode == "LoadControl":
            ops.integrator("LoadControl", load_increment * scale)
        elif res.mode == "DisplacementControl":
            node, dof, du = control
            ops.integrator("DisplacementControl", node, dof, du * scale)

    algo_name = set_algo(algo_i)
    log("start", load_increment=load_increment, has_control=control is not None,
        stabilize=bool(stabilize), dynamics=bool(dynamics))

    # ---------------------------------------------------------------- rung 0-3
    def phase_implicit():
        """Load/disp control with rungs 0-3. Returns 'done' | 'stop' | 'escalate'."""
        nonlocal scale, algo_i, cutbacks, algo_name
        while not done():
            if res.substeps >= max_substeps:
                res.verdict = "incomplete"
                log("budget_exhausted", max_substeps=max_substeps)
                return "stop"
            set_step()
            res.substeps += 1
            if ops.analyze(1) == 0:                 # converged
                res.min_scale_used = min(res.min_scale_used, scale)
                res.algos[algo_name] = res.algos.get(algo_name, 0) + 1
                cutbacks = 0
                if scale < 1.0:
                    scale = min(1.0, scale * grow)
                if algo_i != 0:
                    algo_i = 0
                    algo_name = set_algo(algo_i)
                continue
            # failure: rung 1-2 (walk the algorithm ladder at this scale)
            if algo_i < len(algorithms) - 1:
                algo_i += 1
                algo_name = set_algo(algo_i)
                log("algorithm_switch", to=algo_name, scale=scale)
                continue
            # rung 0 (halve), reset to the easy algorithm
            scale *= 0.5
            algo_i = 0
            algo_name = set_algo(algo_i)
            cutbacks += 1
            # cheap peak/limit detector (R-INSTAB heuristic)
            if can_switch() and cutbacks >= peak_cutbacks:
                res.mode = "DisplacementControl"
                res.switches += 1
                scale, algo_i, cutbacks = 1.0, 0, 0
                algo_name = set_algo(algo_i)
                log("constraint_switch", to="DisplacementControl",
                    reason=f"peak detected ({peak_cutbacks} cutbacks, no commit)")
                continue
            if scale >= min_scale:
                log("cutback", scale=scale)
                continue
            # rung 3 (constraint switch at the cutback floor)
            if can_switch():
                res.mode = "DisplacementControl"
                res.switches += 1
                scale, algo_i, cutbacks = 1.0, 0, 0
                algo_name = set_algo(algo_i)
                log("constraint_switch", to="DisplacementControl",
                    reason="load-control exhausted at floor (rung 3)")
                continue
            # nothing left in the implicit ladder -> escalate to rung 4/5
            log("implicit_exhausted", scale=scale)
            return "escalate"
        return "done"

    # ------------------------------------------------------------------ rung 4
    def phase_stabilized():
        """Viscous-stabilized load control. Returns 'done' | 'stop' | 'escalate'."""
        res.mode = "Stabilized"
        res.switches += 1
        # Issue the STATEFUL integrator + true-equilibrium test ONCE. No -adaptStab
        # (it re-weakens c and prevents crossing a hard limit; measured).
        ops.integrator("LadrunoArcLength", load_increment, 1.0, "-stabilize", stab_f)
        ops.test("LadrunoStabilizedUnbalance", stab_tol, 50, 0)
        algo = set_algo(0)
        clean = 0
        stab_cutbacks = 0
        log("stabilize_arm", f=stab_f, gate=stab_hard_gate)
        while not done():
            if res.substeps >= max_substeps:
                res.verdict = "incomplete"
                log("budget_exhausted", max_substeps=max_substeps)
                return "stop"
            res.substeps += 1
            if ops.analyze(1) == 0:
                res.algos[algo] = res.algos.get(algo, 0) + 1
                res.dissipation_ratio = ops.ladrunoArcLength("dissipationRatio")
                res.stab_dissipated = ops.ladrunoArcLength("dissipatedEnergy")
                if res.dissipation_ratio > stab_hard_gate and not res.over_diffused:
                    res.over_diffused = True
                    log("diffusion_gate_exceeded",
                        ratio=res.dissipation_ratio, gate=stab_hard_gate)
                clean += 1
                if clean >= stab_rampdown_window:    # R-RAMPDOWN: decay c
                    clean = 0
                    ops.ladrunoArcLength("scaleCVisc", 0.5)
                    log("rampdown", ratio=res.dissipation_ratio)
                continue
            # failed: shrink the stabilized load increment and retry
            clean = 0
            stab_cutbacks += 1
            ops.ladrunoArcLength("reduceStep", 0.5)
            log("stabilize_cutback", n=stab_cutbacks)
            if stab_cutbacks >= stab_max_cutbacks:
                log("stabilize_stalled", reason="reduceStep floor (rung 4)")
                return "escalate"
        return "done"

    # ------------------------------------------------------------------ rung 5
    def phase_dynamics():
        """Atomic implicit->DR handoff, settle, return. Returns 'done' | 'stop'."""
        res.mode = "Dynamics"
        res.switches += 1
        lam = ops.getTime()                          # R-HANDOFF: snapshot lambda
        res.dr_lambda = lam
        log("handoff_to_DR", lam=lam)
        ops.loadConst("-time", lam)                  # freeze load + pseudo-time
        ops.wipeAnalysis()
        if dr_setup is not None:
            dr_setup(ops)
        else:
            ops.constraints("Plain")
            ops.numberer("RCM")
            ops.system("BandGen")
        ops.test("NormUnbalance", 1.0e-3, 25, 0)     # DR's own gate is force-based
        ops.algorithm("Newton")
        ops.integrator("LadrunoDynamicRelaxation")   # gershgorin M*, pseudo dt=1
        ops.analysis("Transient")
        # mass-free settling threshold: residualNorm < dr_settle_tol * ||P||
        pref = dr_pref if dr_pref is not None else max(1.0, abs(lam))
        settle_abs = dr_settle_tol * pref
        for _ in range(dr_max_steps):
            if res.substeps >= max_substeps:
                res.verdict = "incomplete"
                log("budget_exhausted", max_substeps=max_substeps)
                return "stop"
            res.substeps += 1
            if ops.analyze(1, 1.0) != 0:
                log("dr_step_failed")
                break
            rn = ops.ladrunoDR("residualNorm")
            if rn < settle_abs:                      # mass-free quasi-static gate
                res.dr_settled = True
                log("dr_settled", residualNorm=rn, pref=pref, steps=res.substeps)
                break
            if done():
                break
        # R-HANDOFF: restore lambda + freeze before any resume
        ops.setTime(lam)
        ops.loadConst("-time", lam)
        log("handoff_from_DR", lam=lam, settled=res.dr_settled)
        return "done" if (res.dr_settled or done()) else "stop"

    try:
        status = phase_implicit()
        if status == "done":
            res.completed = True
            res.verdict = "equilibrium"
            log("done", verdict=res.verdict)
            return res
        if status == "stop":
            return res

        # rung 4 -----------------------------------------------------------
        if stabilize:
            status = phase_stabilized()
            if status == "done":
                res.completed = True
                # Stabilized is NEVER equilibrium. Without the c-reduction
                # diffusion bound (deferred), the honest stamp is `unverified`
                # (R-DIFFUSION: "if the drift was not computed, stamp unverified").
                res.verdict = "regularized" if res.peak_drift is not None \
                    else "unverified"
                log("done", verdict=res.verdict,
                    dissipation_ratio=res.dissipation_ratio,
                    stab_dissipated=res.stab_dissipated,
                    over_diffused=res.over_diffused)
                return res
            if status == "stop":
                return res

        # rung 5 -----------------------------------------------------------
        if dynamics:
            status = phase_dynamics()
            if status == "done":
                res.completed = done()
                # A DR rest state is regularized, not a traced equilibrium branch.
                res.verdict = "regularized"
                log("done", verdict=res.verdict, dr_settled=res.dr_settled)
                return res
            if status == "stop":
                return res

        res.verdict = "incomplete"
        log("give_up", scale=scale)
        return res
    finally:
        if _log_fh is not None:
            _log_fh.close()


# --------------------------------------------------------------------------
# self-test: a single softening Concrete02 truss. Load control cannot pass the
# strength peak; the rung-3 switch traces the softening branch. Mirrors
# robust_solve_tests/torture_softening.py.
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
