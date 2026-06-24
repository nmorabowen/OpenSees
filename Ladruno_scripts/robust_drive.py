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
                 "peak_drift", "dr_settled", "over_diffused", "dr_lambda", "peak_load")

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
        self.peak_load = 0.0            # peak |load factor| during a rung-4 stabilized phase

    def __bool__(self):
        return self.completed and self.verdict == "equilibrium"

    def __repr__(self):
        return (f"<RobustResult {self.verdict} completed={self.completed} "
                f"substeps={self.substeps} switches={self.switches} "
                f"mode={self.mode} min_scale={self.min_scale_used:g}>")


def diffusion_drift(run_at_f, f, factor=0.5):
    """The c-reduction diffusion bound (ADR-31 R-DIFFUSION), standalone.

    A viscous-stabilized solution depends on the artificial dissipation fraction
    `f`; `LadrunoStabilizedUnbalance` proves a point is in equilibrium but NOT that
    the path is `f`-independent. The only trustworthy fidelity certificate is to
    re-run with HALF the viscosity and check the answer barely moved.

    run_at_f : callable(f) -> float. Rebuilds the model, drives it stabilized at the
               dissipated-energy fraction `f`, and returns a comparable scalar metric
               (typically the peak load factor, or a control-DOF response).
    f        : the fraction the headline run used.
    factor   : viscosity reduction for the check (default 0.5 -> half).

    Returns the RELATIVE drift |m(f) - m(factor*f)| / max(|m(f)|, |m(factor*f)|).
    Small (e.g. < 0.02) => the stabilized answer is essentially `f`-independent and
    can be trusted (`regularized`); large => the viscosity is shaping the answer
    (`unverified`). This is what `robust_drive(verify_rebuild=...)` computes
    internally; use it directly to audit any `-stabilize` run.
    """
    m1 = float(run_at_f(f))
    m2 = float(run_at_f(f * factor))
    return abs(m1 - m2) / max(abs(m1), abs(m2), 1e-30)


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
                 verify_rebuild=None,
                 verify_tol=0.02,
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
    verify_rebuild : optional callable(f) -> peak_load. The c-reduction diffusion
                     bound (R-DIFFUSION): if a rung-4 phase reaches done(), the
                     driver calls verify_rebuild(stab_f/2) to re-run the analysis at
                     HALF the viscosity and compares the peak load factor. Drift <=
                     verify_tol -> verdict `regularized` (f-insensitive, trustworthy);
                     otherwise `unverified`. None -> the bound is NOT computed and the
                     stamp stays `unverified` (never silently trusted). See the
                     module-level `diffusion_drift` helper for the standalone form.
    verify_tol     : peak-load drift threshold for the c-reduction bound (default 2%).
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
                res.peak_load = max(res.peak_load, abs(ops.getTime()))  # |λ| for the c-bound
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
                # Stabilized is NEVER equilibrium. The c-reduction diffusion bound
                # (R-DIFFUSION) is the only trustworthy fidelity certificate: re-run
                # the stabilized analysis at HALF stab_f and bound the peak-load
                # drift. `verify_rebuild(f)` rebuilds the model, drives it stabilized
                # at fraction f, and returns the peak |load factor| it reached.
                #   - drift <= verify_tol  -> `regularized` (f-insensitive: trustworthy)
                #   - drift  > verify_tol  -> `unverified`  (the viscosity moved the answer)
                #   - no verify_rebuild    -> `unverified`  ("drift not computed")
                if verify_rebuild is not None:
                    try:
                        peak_half = float(verify_rebuild(stab_f * 0.5))
                        denom = max(abs(res.peak_load), abs(peak_half), 1e-30)
                        res.peak_drift = abs(res.peak_load - peak_half) / denom
                        res.verdict = "regularized" if res.peak_drift <= verify_tol \
                            else "unverified"
                    except Exception as exc:                 # verify must never crash the run
                        res.verdict = "unverified"
                        log("verify_failed", error=str(exc))
                else:
                    res.verdict = "unverified"
                log("done", verdict=res.verdict,
                    dissipation_ratio=res.dissipation_ratio,
                    stab_dissipated=res.stab_dissipated,
                    peak_drift=res.peak_drift, over_diffused=res.over_diffused)
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


# ==========================================================================
# W1-I1a (ADR-52): TRANSIENT adaptive-dt lane.
#
# A Python re-host of OpenSees' built-in VariableTimeStepDirectIntegrationAnalysis
# (the 5-arg `ops.analyze(n, dt, dtMin, dtMax, Jd)`): it reproduces that
# iteration-count step sizing -- grow dt when the last step converged in fewer than
# `target_iters`, shrink it when it took more -- and adds this driver's
# algorithm-ladder escalation + honest observability on top, stepping via
# `ops.analyze(1, dt)` so every decision is visible/loggable.
# ==========================================================================
class TransientResult:
    __slots__ = ("completed", "substeps", "cutbacks", "min_dt_used", "max_dt_used",
                 "verdict", "algos", "final_time", "halfres_max", "n_overtol")

    def __init__(self):
        self.completed = False
        self.substeps = 0
        self.cutbacks = 0
        self.min_dt_used = float("inf")
        self.max_dt_used = 0.0
        # accuracy_gated | integrated | incomplete | aborted
        self.verdict = "incomplete"
        self.algos = {}
        self.final_time = 0.0
        self.halfres_max = 0.0          # W1-I1b: max half-increment residual over accepted steps
        self.n_overtol = 0              # W1-I1b: # accepted steps whose half-res exceeded haftol

    def __bool__(self):
        # Truthy = the integration ran to t_end. `accuracy_gated` (W1-I1b) is the strong
        # verdict: every committed step met the half-increment-residual tolerance. `integrated`
        # is the weak one (W1-I1a convergence-only, OR gate-on but >=1 step exceeded haftol --
        # see n_overtol). Neither is an accuracy certificate unless it is `accuracy_gated`.
        return self.completed and self.verdict in ("integrated", "accuracy_gated")

    def __repr__(self):
        extra = (f" halfres_max={self.halfres_max:g} n_overtol={self.n_overtol}"
                 if self.halfres_max or self.n_overtol else "")
        return (f"<TransientResult {self.verdict} completed={self.completed} "
                f"substeps={self.substeps} cutbacks={self.cutbacks} "
                f"dt=[{self.min_dt_used:g},{self.max_dt_used:g}] t={self.final_time:g}{extra}>")


def _transient_node_state(ops):
    """Snapshot the COMMITTED nodal disp/vel/accel of every node {tag: (d, v, a)}.

    nodeDisp/Vel/Accel return the committed vectors, so this is the converged-step state
    (untouched by the trial-state injection the half-increment gate performs)."""
    state = {}
    for tag in ops.getNodeTags():
        state[tag] = (list(ops.nodeDisp(tag)),
                      list(ops.nodeVel(tag)),
                      list(ops.nodeAccel(tag)))
    return state


def _half_increment_residual(ops, pre, post, dt):
    """The Abaqus-style half-increment residual after a converged step t_n -> t_n+dt.

    Builds the constant-average-acceleration half-step state (exact for the trapezoidal rule,
    a faithful approximation for HHT/generalized-alpha):
        a_half = (a_n + a_{n+1}) / 2
        v_half = v_n + dt/4 (a_n + a_{n+1})
        u_half = u_n + dt/2 v_n + dt^2/16 (a_n + a_{n+1})
    injects it into the node TRIAL vectors (`ladrunoSetNodeTrial`, full-vector so multi-dof
    nodes work), then reads the free-DOF dynamic unbalance norm at that state
    (`ladrunoTrialResidualNorm`, which runs the element update() OpenSeesPy otherwise omits).
    The residual is nonzero because the integrator only enforced equilibrium at the step ENDS,
    not the midpoint -- so its size is a local-truncation-error proxy. Trial state is restored
    to the converged t_{n+1} values; committed state is never touched.

    FIDELITY (be honest): exact for rate- and path-INDEPENDENT (elastic) materials. For inelastic
    or rate-dependent (plastic/viscous/creep) materials it is APPROXIMATE -- the residual is
    formed AFTER the step commits, so the element reference (committed) state is t_{n+1} not t_n,
    and the helper imposes a representative positive half-step dt for rate terms. Step 1 also
    leans on the committed initial acceleration a0, which OpenSees leaves at 0 unless the caller
    established it (no Abaqus-style self-start solve). NOTE: between accepted steps the element
    trial state is parked at the half-step until the next analyze() refreshes it -- harmless for
    the committed trajectory and node-level done()/results (nodeDisp etc. read COMMITTED state),
    but an inter-step eleForce/recorder read would see the half-step forces."""
    for tag, (d0, v0, a0) in pre.items():
        d1, v1, a1 = post[tag]
        n = len(d0)
        uh = [d0[i] + 0.5 * dt * v0[i] + (dt * dt / 16.0) * (a0[i] + a1[i]) for i in range(n)]
        vh = [v0[i] + 0.25 * dt * (a0[i] + a1[i]) for i in range(n)]
        ah = [0.5 * (a0[i] + a1[i]) for i in range(n)]
        ops.ladrunoSetNodeTrial(tag, *uh, *vh, *ah)
    t_mid = ops.getTime() - 0.5 * dt                 # step midpoint (post-step time is t_{n+1})
    r = float(ops.ladrunoTrialResidualNorm(t_mid))   # load evaluated at the midpoint
    for tag, (d1, v1, a1) in post.items():           # restore trial -> converged t_{n+1}
        ops.ladrunoSetNodeTrial(tag, *d1, *v1, *a1)
    return r


def robust_transient(ops, t_end, dt, *,
                     dt_min=None,
                     dt_max=None,
                     target_iters=5,
                     done=None,
                     algorithms=_DEFAULT_LADDER,
                     grow=2.0,
                     shrink=0.5,
                     max_substeps=1_000_000,
                     time_tol=None,
                     error_gate=False,
                     haftol=None,
                     error_order=2.0,
                     log_path=None,
                     verbose=False):
    """Integrate a TRANSIENT analysis to `t_end` with convergence-driven adaptive dt.

    Re-hosts the built-in variable-step transient analysis in Python so the robust
    driver can adapt dt across the transient integrators (Newmark/HHT/GeneralizedAlpha
    /...). On each step it calls `ops.analyze(1, dt)`; on success it sizes the next dt
    from the iteration count (`ops.testIter()`): factor = target_iters / lastIters,
    clamped to [shrink, grow], then clamped to [dt_min, dt_max]. On failure it walks
    the algorithm ladder, then halves dt down to dt_min before aborting.

    IMPORTANT (ADR-52 W1-I1a caveat): iteration-count control is a COST / convergence
    criterion, NOT an accuracy one. A stiff-but-stable response can converge cheaply
    and let dt grow PAST where the solution is accurate. The verdict is therefore
    `integrated` (the integration completed), never an accuracy certificate -- the
    half-increment-residual ERROR gate (W1-I1b) is what makes a transient run
    accuracy-grade. `bool(res)` is True for `integrated`; treat it as "ran to t_end",
    not "is accurate".

    t_end        : integrate until ops.getTime() >= t_end (or done() returns True).
    dt           : initial / nominal time step.
    dt_min,dt_max: adaptive bounds (defaults dt/1024 and dt*64).
    target_iters : Jd-style target iteration count used for sizing (default 5).
    done         : optional callable() -> bool early-stop predicate (e.g. a DOF target).
    The C++ equivalent is `ops.analyze(n, dt, dt_min, dt_max, target_iters)`; this lane
    adds the algorithm ladder + per-step observability + the honest verdict on top.

    W1-I1b -- THE ACCURACY GATE (error_gate=True):
    Upgrades the sizing from convergence-only to accuracy-grade with the Abaqus-style
    half-increment residual. After each converged step the driver evaluates the dynamic
    unbalance at the constant-average-acceleration MIDPOINT of the step (see
    `_half_increment_residual`); its norm is a local-truncation-error proxy that, unlike the
    iteration count, actually senses when dt has outrun accuracy. dt sizing then takes the
    MOST CONSERVATIVE of the iteration-count factor and the accuracy factor
    `(haftol / r_half) ** (1/error_order)` (residual ~ O(dt^error_order); error_order=2 default).

    Because OpenSees commits a converged step before the residual can be seen (no Python
    pre-commit hook), this is a FEED-FORWARD gate: a too-large step is committed, and its
    residual shrinks the NEXT step -- it does NOT reject/redo the step the way Abaqus does.
    The verdict tells the truth about that:
      * `accuracy_gated` (truthy): completed AND every committed step met haftol (n_overtol==0).
      * `integrated`     (truthy): completed but >=1 committed step exceeded haftol
                                   (n_overtol>0) -- the gate flagged inaccuracy it could not
                                   pre-empt. Start with a conservative dt to keep n_overtol==0.
    `res.halfres_max` / `res.n_overtol` carry the evidence. Requires the fork build
    (`ladrunoTrialResidualNorm` / `ladrunoSetNodeTrial`).

    error_gate   : enable the half-increment-residual accuracy gate (default off -> W1-I1a).
    haftol       : half-increment-residual tolerance (force units, problem-scaled). REQUIRED
                   when error_gate is on. A step whose r_half exceeds it is accuracy-suspect.
    error_order  : assumed convergence order of r_half in dt for the accuracy factor (default 2).
    """
    res = TransientResult()
    if dt_min is None:
        dt_min = dt / 1024.0
    if dt_max is None:
        dt_max = dt * 64.0
    if time_tol is None:
        time_tol = dt_min * 0.5

    if error_gate:
        if haftol is None or haftol <= 0.0:
            raise ValueError("robust_transient(error_gate=True) requires a positive haftol "
                             "(half-increment-residual tolerance, force units).")
        if not hasattr(ops, "ladrunoTrialResidualNorm") or not hasattr(ops, "ladrunoSetNodeTrial"):
            raise RuntimeError("error_gate needs the fork build: ladrunoTrialResidualNorm / "
                               "ladrunoSetNodeTrial commands are missing (rebuild the Ladruno fork).")

    _log_fh = open(log_path, "w", encoding="utf-8") if log_path else None

    def log(event, **kw):
        if _log_fh is not None:
            rec = {"event": event, "substep": res.substeps, "t": ops.getTime()}
            rec.update(kw)
            _log_fh.write(json.dumps(rec) + "\n")
            _log_fh.flush()
        if verbose:
            print(f"  [transient] {event} {kw}")

    def set_algo(i):
        spec = algorithms[i % len(algorithms)]
        ops.algorithm(*spec)
        return spec[0]

    def reached():
        if done is not None and done():
            return True
        return ops.getTime() >= t_end - time_tol

    cur = dt
    algo_i = 0
    algo_name = set_algo(0)
    log("start", t_end=t_end, dt=dt, dt_min=dt_min, dt_max=dt_max,
        target_iters=target_iters, error_gate=bool(error_gate), haftol=haftol)

    try:
        while not reached():
            if res.substeps >= max_substeps:
                res.verdict = "incomplete"
                log("budget_exhausted", max_substeps=max_substeps)
                break
            step = min(cur, t_end - ops.getTime())   # never overshoot t_end
            if step <= time_tol:
                break
            pre = _transient_node_state(ops) if error_gate else None  # t_n committed state
            res.substeps += 1
            if ops.analyze(1, step) == 0:            # converged
                res.algos[algo_name] = res.algos.get(algo_name, 0) + 1
                res.min_dt_used = min(res.min_dt_used, step)
                res.max_dt_used = max(res.max_dt_used, step)
                it = 0
                try:
                    it = int(ops.testIter())
                except Exception:
                    it = 0
                factor = (target_iters / it) if it > 0 else grow
                if error_gate:
                    # W1-I1b: half-increment-residual accuracy factor, combined conservatively.
                    r_half = _half_increment_residual(ops, pre, _transient_node_state(ops), step)
                    res.halfres_max = max(res.halfres_max, r_half)
                    over = r_half > haftol
                    if over:
                        res.n_overtol += 1
                    log("halfres", r=r_half, haftol=haftol, dt=step, over=over)
                    acc_factor = (haftol / r_half) ** (1.0 / error_order) if r_half > 0.0 else grow
                    factor = min(factor, acc_factor)
                factor = max(shrink, min(grow, factor))
                cur = min(dt_max, max(dt_min, cur * factor))
                if algo_i != 0:                      # recovered -> back to easy algo
                    algo_i = 0
                    algo_name = set_algo(0)
                continue
            # failure: walk the algorithm ladder at this dt (rung 1-2)
            if algo_i < len(algorithms) - 1:
                algo_i += 1
                algo_name = set_algo(algo_i)
                log("algorithm_switch", to=algo_name, dt=step)
                continue
            # rung 0: halve dt, reset to the easy algorithm
            algo_i = 0
            algo_name = set_algo(0)
            cur *= shrink
            res.cutbacks += 1
            if cur >= dt_min:
                log("cutback", dt=cur)
                continue
            res.verdict = "aborted"
            log("dt_floor", dt=cur, dt_min=dt_min)
            break

        res.final_time = ops.getTime()
        if reached() and res.verdict != "aborted":
            res.completed = True
            # accuracy_gated only if the gate ran AND every committed step met haftol.
            res.verdict = ("accuracy_gated" if (error_gate and res.n_overtol == 0)
                           else "integrated")
        log("done", verdict=res.verdict, substeps=res.substeps,
            cutbacks=res.cutbacks, final_time=res.final_time,
            halfres_max=res.halfres_max, n_overtol=res.n_overtol)
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


# --------------------------------------------------------------------------
# self-test: a linear SDOF (mass on an elastic truss spring) under a constant
# load, integrated with Newmark via the W1-I1a transient lane. The linear
# problem converges in one iteration every step, so the iteration-count sizing
# grows dt to the cap -- exercising the adaptive-growth path and illustrating
# the W1-I1a caveat (growth is convergence-driven, not accuracy-driven). The
# integrator is unconditionally stable, so the run completes regardless.
# --------------------------------------------------------------------------
def _selftest_transient():
    import os
    import sys
    DIST = os.environ.get("LADRUNO_DIST",
                          r"C:\Users\nmb\Documents\Github\OpenSees\dist\bin")
    if os.path.isdir(DIST):
        os.add_dll_directory(DIST)
        sys.path.insert(0, DIST)
    import opensees as ops

    E, A, L, m, P = 1000.0, 1.0, 1.0, 1.0, 1.0   # k=1000, omega~31.6, T~0.199
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, L)
    ops.fix(1, 1)
    ops.mass(2, m)
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.element("Truss", 1, 1, 2, A, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, P)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("BandGen")
    ops.test("NormDispIncr", 1e-10, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)         # unconditionally stable
    ops.analysis("Transient")

    dt0, t_end = 0.005, 1.0
    res = robust_transient(ops, t_end, dt0, verbose=True)

    t = ops.getTime()
    u = ops.nodeDisp(2, 1)
    u_static = P / (E * A / L)
    grew = res.max_dt_used > dt0                  # adaptive growth happened
    bounded = -0.1 * u_static <= u <= 2.1 * u_static
    ok = bool(res) and abs(t - t_end) < 1e-6 and grew and bounded
    print(f"\n{res}")
    print(f"final t={t:.5f}  u={u:.6e}  u_static={u_static:.6e}  "
          f"max_dt={res.max_dt_used:g} (dt0={dt0})")
    print("TRANSIENT SELF-TEST:",
          "PASS - integrated to t_end with adaptive dt growth"
          if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    rc = 0
    if which in ("all", "static"):
        rc |= _selftest()
    if which in ("all", "transient"):
        rc |= _selftest_transient()
    raise SystemExit(rc)
