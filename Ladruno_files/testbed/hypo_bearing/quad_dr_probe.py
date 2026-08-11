"""quad_dr_probe -- the TANGENT-FREE path probe for the quadratic-element wall.

THE QUESTION.  Notes 81 and 82 established that on the Prandtl-Reissner collapse
problem (exact oracle q_u = q0*N_q) the LINEAR relieved hex reaches the exact
answer on a genuine plateau, while EVERY quadratic leg -- hex and tet, Lagrange
and Bernstein, relieved and locked -- stops early on a path-controller limit at
s/B ~ 0.001-0.011 against the linear element's 0.031-0.074.  Four explanations
are dead (mesh alignment; the constraint ratio r; "the tail slope splits on
order"; and span/mechanism deficiency, refuted directly at MATCHED settlement).
Note 82 sec. 7.3 measured the sharpest surviving clue: the quadratic tangent is
FLAT right up to the wall and then loses ~1000x in sigma_min across 0.06 mm of
settlement -- abrupt, not asymptotic, i.e. a numerical event and not a limit
point.  Note 82 sec. 8 item 3 names the experiment that settles the remaining
fork in the road, and this file is it:

    DOES A QUADRATIC ELEMENT FINISH THIS PROBLEM WHEN THE PATH DOES NOT DEPEND
    ON A TANGENT?

`LadrunoDynamicRelaxation` (33005, ADR 21) is quasi-static DR on an explicit
central-difference core with a fictitious Gershgorin mass and Cundall kinetic
damping.  No tangent is formed, factorized or trusted anywhere in the stepping
loop: `formTangent` pokes the diagonal M* onto the SOE and the solve degenerates
to an M*^-1 apply.  So a leg that completes under DR removes the whole
Newton/conditioning family of explanations; a leg that walls under DR removes it
the other way and points at the discretization or the material return map.

THE CONTROL COMES FIRST AND IS NOT OPTIONAL.  Leg `h8bbar` at h0 = 0.5 is the
KNOWN-GOOD leg: the merged R3 gate reads 0.9938 there and note 82's driver reads
0.9977, both plateaued.  If DR cannot reproduce a result we already know, then
nothing measured afterwards means anything and the DR-validation failure IS the
finding.  Run it before any quadratic leg.

WHAT "FINISHES" MEANS HERE -- AND A CORRECTION TO THE OBVIOUS FRAMING.  The
tempting DR acceptance rule is "at collapse the structure never settles, it just
keeps accelerating".  That rule is correct under LOAD control and WRONG here,
and the difference is structural, not a detail:

    this problem is DISPLACEMENT-controlled (a rigid footing, sp() on dof 3 of
    every footing node).  The collapse mechanism carries footing settlement, so
    under displacement control the mechanism is DRIVEN, not free -- it is not in
    the span of the free DOFs at all.  The constrained operator therefore stays
    non-singular through collapse, which is exactly what note 82 sec. 7.3
    MEASURED for the linear element (sigma_min flat, 1.2e-3, all the way to a
    fully formed mechanism).  So every increment settles, at and past collapse,
    and "does it settle" is NOT the plateau test here.

Keeping the loading identical to the established legs is worth more than
inheriting a slogan, so the acceptance is built on what displacement control
actually offers:

  * per increment, the imposed settlement is HELD and DR is relaxed to static
    rest -- so every recorded (s, q) point is a genuine STATIC equilibrium
    point, directly comparable with the static legs' converged points;
  * "settled" = the TRUE STATIC UNBALANCE ||f_ext - f_int||_inf has fallen
    below 1e-6 of the total surcharge (3.0e-4 kN), i.e. ten times tighter in
    absolute force than the static ladder's own NormUnbalance acceptance.  The
    residual is the whole gate on purpose: a per-chunk "the reaction stopped
    moving" test is chunk-length dependent -- shorten the chunk and the same
    state passes -- so it is reported in the CSV and decides nothing;
  * the COLLAPSE signature remains the q(s) PLATEAU -- unchanged from the
    static campaign, so the two are on the same footing;
  * an increment that runs out of DR iterations is the DR analogue of the static
    controller's step floor.  It is an ALLOWANCE, not a capacity, and the leg is
    classified accordingly.

CAPACITY RULE (the merged R3 gate's three clauses, transposed to DR):
  1. plateau -- tail dq/ds < 2 % of the initial tangent;
  2. termination mode TARGET (reached the s/B target).  ITERCAP / DIVERGED /
     WALL / TRUNCATED are seizure modes and are never a capacity;
  3. free advance -- over the final tenth of the run no increment spent more
     than FREE_ADVANCE_FRACTION of its iteration allowance.
A leg failing any clause is a CONTROLLER ALLOWANCE and the allowance is named.

THE dt NULL-CONTROL (a prediction, then a measurement).  With `-mass gershgorin`
the fictitious mass is m_i = (dt^2/4)*sum_j|K_ij|, so the leap-frog update is
  du = dt^2 * a = dt^2 * R/m = 4*R / sum_j|K_ij|
and dt CANCELS EXACTLY.  KE = 1/2 v^T M* v likewise scales as du^2, so even the
kinetic-damping decisions are dt-free.  The march must therefore be BIT-
IDENTICAL for any dt, provided the integrator's `-dt` matches the analyze dt.
`--knob dt` runs that pair; anything other than an exact tie is a bug in the
harness or the integrator, not a physical result.  It is a self-consistency
control with a predicted exact null, which is the cheapest kind worth having.

EVERY PROBLEM PARAMETER IS IMPORTED, NOT RE-IMPLEMENTED.  The mesh, the
consistent Q8 surcharge and its two verifications, the 1-D elastic stress patch
check, the spurious-mode census, the material, the cone, the oracle and the
surcharge step all come from `h20_prandtl.py` by import.  phi_txc = 20, q0 = 10,
gamma = 0, nu = 0.45, non-associated rho_bar = 0, Chen & Han plane-strain match.

Run:
    py -3.12 quad_dr_probe.py --elem h8bbar --h0 0.5            # THE CONTROL
    py -3.12 quad_dr_probe.py --elem h20uri --h0 0.5
    py -3.12 quad_dr_probe.py --elem h8bbar --h0 0.5 --dt 0.1 --suffix _dt01
"""
import argparse
import csv
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h20_prandtl as H                                          # noqa: E402

ops = H.ops
log = H.log
HERE = os.path.dirname(os.path.abspath(__file__))

ELEMS = {"h8bbar": (1, "bbar"), "h8std": (1, "std"),
         "h20uri": (2, "uri"), "h20std": (2, "std")}

DROP_TRUNCATE = 0.02          # h20_prandtl's truncation rule, unchanged
FREE_ADVANCE_FRACTION = 0.50  # clause 3: iteration spend over the final tenth
STABLE_CHUNKS = 2             # consecutive quiet chunks before "settled"


# --------------------------------------------------------------------------
def dr_options(args):
    """The integrator argument list.  `-dt` MUST equal the analyze dt (see the
    dt null-control in the docstring)."""
    o = ["-mass", args.mass]
    if args.mass == "lumped":
        o.append(args.mass_scale)
    o += ["-dt", args.dt]
    if args.damping == "viscous":
        o += ["-damping", "viscous", args.zeta]
    else:
        o += ["-damping", "kinetic"]
    if args.no_auto_refresh:
        o.append("-noAutoRefresh")
    if args.recompute > 0:
        o += ["-recompute", args.recompute]
    return o


def set_dr_analysis(args):
    ops.wipeAnalysis()
    ops.constraints("Transformation")   # SP'd dofs leave the SOE; the prescribed
    ops.numberer("Plain")               # value is enforced exactly (enforceSPs)
    ops.system("Diagonal")              # makes the M*-solve a diagonal apply
    ops.test("NormUnbalance", 1.0e-12, 1, 0)   # unused by `Linear`; 1 iteration
    ops.algorithm("Linear")             # exactly one solve per DR step
    ops.integrator("LadrunoDynamicRelaxation", *dr_options(args))
    ops.analysis("Transient")


def set_push(foot, uz_target, first):
    """Impose the TOTAL footing displacement uz_target through a Constant-series
    pattern.  `enforceSPs` sets the node's total trial disp to the SP value
    (TransformationDOF_Group.cpp:1058), so the value is an ABSOLUTE
    displacement, exactly as in the static harness where pseudo-time IS the
    imposed settlement.  Re-imposing means dropping and rebuilding the pattern:
    a second sp() on the same dof would ADD a duplicate SP_Constraint."""
    if not first:
        ops.remove("loadPattern", 2)
    ops.pattern("Plain", 2, 2)
    for t in foot:
        ops.sp(t, 3, float(uz_target))


def ramp_push(args, foot, uz_from, uz_to):
    """Take the footing from uz_from to uz_to as `--ramp` equal sub-jumps
    separated by `--rampevery` DR steps, instead of one discontinuity.

    WHY THIS IS NOT OPTIONAL, AND WHAT IT COST TO LEARN.  Stepping the whole
    increment at once and then relaxing gives an exact static equilibrium -- the
    residual really does go to 1e-13 kN -- but it is the equilibrium of a
    DIFFERENT LOAD HISTORY.  A displacement discontinuity at the footing
    launches a transient through a mesh whose mass is FICTITIOUS, so the
    transient has no physical wave speed and the near-field elements absorb the
    whole jump in a thin layer.  On an elastic material that rings and relaxes
    away harmlessly (measured: DR and Newton agree to 1.000000 on the elastic
    model, `dr_elastic_check.py`).  On an ELASTIC-PLASTIC material the overshoot
    leaves PLASTIC STRAIN behind, and the settled state is permanently softer:
    measured, a single 1 mm step increment landed 29 % below the static Newton
    answer at the same settlement while reporting a perfectly converged
    residual.  A settled residual proves the state is in equilibrium; it proves
    nothing about which path got there."""
    n = max(1, args.ramp)
    for i in range(1, n + 1):
        set_push(foot, uz_from + (uz_to - uz_from) * i / n, first=False)
        if i < n and args.rampevery > 0:
            if ops.analyze(args.rampevery, args.dt_march) != 0:
                return False
    return True


def continuous_push(args, foot, uz0):
    """Drive the footing down at a CONSTANT rate with no pattern rebuilds at all.

    This is the cheap, shock-free limit of `ramp_push`, and it exists because of
    the dt null-control (module docstring sec. 'THE dt NULL-CONTROL'): with
    `-mass gershgorin` the march depends only on the RATIO of the analyze dt to
    the integrator's `-dt`, never on their absolute size.  So the pseudo-time
    step is free to be repurposed as the loading clock.  Set `-dt rate/dtfac`
    and march at `dt = rate`; the domain clock then advances by exactly `rate`
    per DR step, and a `Linear` time series carrying `sp(node, 3, -1.0)` turns
    that clock directly into the imposed footing displacement:

        uz = -1.0 * t,   t0 = -uz0  =>  uz = uz0 - n*rate  after n steps.

    The dynamics are BIT-IDENTICAL to the same run at dt = 1.0 (measured), so
    nothing about the relaxation is traded away; what is bought is that the
    displacement is now continuous in the DR step rather than a staircase of
    sub-jumps, and that no `domainChanged` is triggered on the way -- which is
    where the ramped version spends most of its wall clock."""
    ops.remove("loadPattern", 2)          # the Constant-series pattern DR-0 used
    ops.remove("timeSeries", 2)
    _arm_ramp(foot, uz0)


def _arm_ramp(foot, uz_now):
    """(Re)install the constant-rate clock so that uz = -1.0 * t = uz_now now."""
    ops.setTime(-uz_now)
    ops.timeSeries("Linear", 4)
    ops.pattern("Plain", 2, 4)
    for t in foot:
        ops.sp(t, 3, -1.0)


def hold_check(args, foot, r_corr, res_abs, uz_now):
    """Freeze the ramp at the current settlement and relax to static rest.

    A continuously ramped point is quasi-static, not static: it carries whatever
    inertial lag the rate leaves behind.  This measures that lag directly rather
    than assuming it away -- swap the Linear-series pattern for a Constant one
    at the CURRENT displacement, relax to the same residual gate every held leg
    uses, and report how far q moved.  The ramp is then re-armed at the clock
    value the held displacement corresponds to, so the leg continues on the same
    path it was on."""
    q_before = reaction(foot, r_corr)
    ops.remove("loadPattern", 2)
    ops.remove("timeSeries", 4)
    ops.timeSeries("Constant", 3)
    ops.pattern("Plain", 2, 3)
    for t in foot:
        ops.sp(t, 3, float(uz_now))
    ok, nit, R, ke, res, why, dR = relax(args, foot, r_corr, res_abs)
    ops.remove("loadPattern", 2)
    ops.remove("timeSeries", 3)
    _arm_ramp(foot, uz_now)
    return q_before, R, nit, ok, why


def reaction(foot, r_corr):
    ops.reactions()          # flag 0 = STATIC reaction (no inertia, no damping)
    return -sum(ops.nodeReaction(t, 3) for t in foot) + r_corr


def relax(args, foot, r_corr, res_abs):
    """Hold the imposed displacement and relax to static rest.

    Returns (settled, iters, R, ke, res, reason, dR).  The gate must hold for
    STABLE_CHUNKS consecutive chunks: kinetic damping zeroes velocities at every
    KE peak, so a momentarily quiet chunk means nothing on its own."""
    n, stable, R_prev = 0, 0, None
    ke = res = float("nan")
    R = float("nan")
    while n < args.maxit:
        if ops.analyze(args.chunk, args.dt_march) != 0:
            return False, n, R, ke, res, "DR analyze failed (divergence)",                 float("nan")
        n += args.chunk
        ke = float(ops.ladrunoDR("kineticEnergy"))
        res = float(ops.ladrunoDR("residualNorm"))
        R = reaction(foot, r_corr)
        if not (math.isfinite(R) and math.isfinite(res)):
            return False, n, R, ke, res, "non-finite state", float("nan")
        if args.trace:
            log(f"        [trace] {n:>7d} steps  R = {R:12.6f} kN  "
                f"KE = {ke:.4e}  res = {res:.4e} kN")
        # The settle gate is the TRUE STATIC UNBALANCE and nothing else.  An
        # earlier version also required |dR| over a chunk to be small, which is
        # wrong for a subtle reason worth recording: a per-chunk change gate is
        # chunk-length dependent -- shorten the chunk and less happens inside
        # it, so the same state passes.  A gate whose verdict moves with a
        # reporting parameter is measuring the reporting parameter.  The
        # residual is chunk-free.  |dR| is still WRITTEN to the CSV so the
        # reader can see stationarity; it just does not decide anything.
        dR = abs(R - R_prev) if R_prev is not None else float("nan")
        stable = stable + 1 if res <= res_abs else 0
        R_prev = R
        if stable >= STABLE_CHUNKS:
            return True, n, R, ke, res, "settled", dR
    return False, n, R, ke, res, f"iteration cap {args.maxit} spent", \
        (abs(R - R_prev) if R_prev is not None else float("nan"))


def classify(args, rows, q_exact, mode, verdict, wall, ndof, nel, tag, path):
    a = np.array([r[:4] for r in rows])
    s, q = a[:, 0], a[:, 3]
    its = np.array([r[4] for r in rows], dtype=float)
    peak = np.maximum.accumulate(q)
    bad = np.where(q < (1.0 - DROP_TRUNCATE) * peak)[0]
    if len(bad):
        ntr = len(q) - bad[0]
        s, q, its = s[:bad[0]], q[:bad[0]], its[:bad[0]]
        verdict += f"; curve TRUNCATED at s/B = {s[-1] / H.B_FOOT:.4f} ({ntr} cut)"
        mode = "TRUNCATED"
    msk = s >= 0.9 * s[-1]
    t_last = np.polyfit(s[msk], q[msk], 1)[0] if msk.sum() > 2 else float("nan")
    n0f = max(4, len(s) // 50)
    t_init = np.polyfit(s[:n0f], q[:n0f], 1)[0]
    plateau = bool(abs(t_last) < 0.02 * abs(t_init))
    tail_it = float(its[msk].max()) if msk.sum() else float(its.max())
    free = bool(args.mode == "continuous"
                or tail_it <= FREE_ADVANCE_FRACTION * args.maxit)
    qmax = float(q.max())
    capacity = bool(plateau and mode == "TARGET" and free)
    log(f"    {len(rows)} points, {int(its.sum())} DR steps, {wall:.0f} s   "
        f"[{verdict}]")
    log(f"    q_max = {qmax:.2f} kPa at s/B = "
        f"{s[int(q.argmax())] / H.B_FOOT:.4f}; end s/B = {s[-1] / H.B_FOOT:.4f}")
    log(f"    tail dq/ds = {t_last:.2f} kPa/m = {100 * t_last / t_init:.3f} % "
        f"of initial  ->  {'PLATEAU' if plateau else 'STILL HARDENING'}")
    log(f"    MODE = {mode}")
    log(f"    q_num / q_exact = {qmax / q_exact:.4f}   -->  "
        f"{'CAPACITY' if capacity else 'CONTROLLER ALLOWANCE (' + mode + ')'}")
    return dict(leg=tag, elem=args.elem, h0=args.h0, ndof=ndof, nel=nel,
                qmax=qmax, q_exact=q_exact, ratio=qmax / q_exact,
                plateau=plateau, t_frac=100 * t_last / t_init,
                s_end=float(s[-1]) / H.B_FOOT, ninc=len(rows),
                iters=int(its.sum()), tail_it=tail_it, mode=mode, free=free,
                capacity=capacity, wall=wall, verdict=verdict, csv=path)


def push_continuous(args, foot, r_corr, uz0, res_abs, smax, area, q_exact,
                    nodes, cells):
    """The constant-rate push: one uninterrupted DR march, sampled every
    `--chunk` steps, with periodic HOLD checks that measure the inertial lag."""
    rate = args.rate
    nsteps = int(round(smax / rate))
    log(f"=== leg {args.elem} h0 = {args.h0}: CONTINUOUS push at "
        f"{rate:.3e} m per DR step to s/B = {args.sfrac} "
        f"({nsteps} DR steps), sampled every {args.chunk}; "
        f"hold check every {args.holdevery} samples")
    log(f"    ORACLE q_u = {q_exact:.3f} kPa")
    args.dt = rate / args.dtfac
    args.dt_march = rate
    set_dr_analysis(args)                     # re-sized -dt, same dynamics
    continuous_push(args, foot, uz0)

    tag = f"dr_{args.elem}_h{str(args.h0).replace('0.', '')}{args.suffix}"
    path = os.path.join(HERE, f"{tag}.csv")
    fh = open(path, "w", newline="")
    wr = csv.writer(fh)
    wr.writerow(["s_m", "s_over_B", "R_kN", "q_kPa", "dr_iters", "settled",
                 "ke", "res_kN", "dR_kN", "wall_s"])
    rows, verdict, mode, holds = [], "reached the target", "TARGET", []
    t0, n, k = time.time(), 0, 0
    while n < nsteps:
        if time.time() - t0 > args.tmax:
            verdict, mode = f"WALL-CLOCK CAP of {args.tmax:.0f} s hit", "WALL"
            break
        m = min(args.chunk, nsteps - n)
        if ops.analyze(m, args.dt_march) != 0:
            verdict, mode = f"DR analyze failed at step {n}", "DIVERGED"
            break
        n += m
        k += 1
        s = uz0 - ops.nodeDisp(foot[0], 3)
        R = reaction(foot, r_corr)
        ke = float(ops.ladrunoDR("kineticEnergy"))
        res = float(ops.ladrunoDR("residualNorm"))
        if not (math.isfinite(R) and math.isfinite(res)):
            verdict, mode = f"non-finite state at step {n}", "DIVERGED"
            break
        rows.append((s, s / H.B_FOOT, R, R / area, m, int(res <= res_abs), ke,
                     res, float("nan"), time.time() - t0))
        wr.writerow([f"{v:.9g}" for v in rows[-1]])
        fh.flush()
        if args.holdevery > 0 and k % args.holdevery == 0:
            qb, Rh, nit, okh, whyh = hold_check(args, foot, r_corr, res_abs,
                                                ops.nodeDisp(foot[0], 3))
            holds.append((s / H.B_FOOT, qb / area, Rh / area, nit, okh))
            log(f"    [HOLD] s/B = {s / H.B_FOOT:.5f}: ramped q = "
                f"{qb / area:8.3f} -> relaxed {Rh / area:8.3f} kPa "
                f"({100 * (qb / Rh - 1):+.3f} % inertial lag), {nit} steps, "
                f"{'settled' if okh else 'NOT SETTLED: ' + whyh}")
        if k % args.report == 0:
            log(f"    n = {n:>9d}/{nsteps}  s/B = {s / H.B_FOOT:.5f}  "
                f"q = {R / area:8.3f} kPa ({R / area / q_exact:.4f} of exact)  "
                f"KE = {ke:.2e}  res = {res:.2e}  [{time.time() - t0:.0f}s]")
    fh.close()
    if not rows:
        raise SystemExit("no points recorded")
    # a FINAL hold check: the headline must be a static point, not a ramped one
    qb, Rh, nit, okh, whyh = hold_check(args, foot, r_corr, res_abs,
                                        ops.nodeDisp(foot[0], 3))
    log(f"    [HOLD, final] ramped q = {qb / area:.3f} -> relaxed "
        f"{Rh / area:.3f} kPa ({100 * (qb / Rh - 1):+.3f} % inertial lag), "
        f"{nit} steps, {'settled' if okh else 'NOT SETTLED: ' + whyh}")
    for h in holds:
        log(f"    [HOLD summary] s/B {h[0]:.5f}: {h[1]:.3f} -> {h[2]:.3f} kPa "
            f"({100 * (h[1] / h[2] - 1):+.3f} %)")
    return classify(args, rows, q_exact, mode, verdict, time.time() - t0,
                    3 * len(nodes), len(cells), tag, path)


# --------------------------------------------------------------------------
def run_leg(args):
    order, form = ELEMS[args.elem]
    q_exact = H.Q0 * H.n_q(H.mc_from_cone(H.alpha_from_phi_txc(H.PHI_TXC))["ps"])
    area = H.B_FOOT * H.THICK

    nodes, cells, sets, x, y, z = H.strip_mesh(args.h0, order)
    w, nfaces = H.consistent_surcharge(nodes, cells, order)
    H.verify_surcharge(nodes, cells, w, nfaces, order)          # controls 2a/2b
    log(f"[mesh] {len(nodes)} nodes / {3 * len(nodes)} DOF, {len(cells)} "
        f"{cells.shape[1]}-node hexes, footing set {len(sets['footing'])} nodes")

    if args.patch:
        H.leg_patch(nodes, cells, sets, w, form)                # control 1
    if args.modes:
        H.leg_modes(nodes, cells, sets, w, form)                # control 4

    H.build_model(nodes, cells, sets, w, form, assoc=False)
    ftot = H.Q0 * float(w.sum())
    tol = 1.0e-5 * max(ftot, 1.0)
    H.surcharge_step(nodes, sets, w, tol)                       # control 3

    foot = [int(n) + 1 for n in sets["footing"]]
    r_corr = H.Q0 * float(w[sets["footing"]].sum())
    uz0 = ops.nodeDisp(foot[0], 3)
    R_static0 = reaction(foot, r_corr)
    res_abs = args.tol_res * ftot
    log(f"    [DR setup] uz0 = {uz0:.9g} m, static R at s = 0 is "
        f"{R_static0:.6f} kN; settle gate = "
        f"||f_ext-f_int||_inf < {res_abs:.4e} kN "
        f"(= {args.tol_res:.0e} x the {ftot:.1f} kN surcharge, i.e. 10x "
        f"TIGHTER in absolute force than the static ladder's own "
        f"NormUnbalance acceptance of {1e-5 * ftot:.4e} kN)")

    ops.loadConst("-time", 0.0)
    ops.timeSeries("Constant", 2)
    set_dr_analysis(args)
    log(f"    [DR setup] integrator LadrunoDynamicRelaxation "
        f"{' '.join(str(v) for v in dr_options(args))}; analyze dt = "
        f"{args.dt_march} (dtfac = {args.dtfac}, so omega_max*dt_march <= "
        f"{2 * args.dtfac:.3f} against the central-difference bound of 2)")

    # --- CONTROL DR-0: the zero-push relaxation ---------------------------
    # Re-impose the CURRENT footing displacement and relax.  The state is
    # already a static equilibrium, so DR must not move it.  This is the
    # cheapest test that the SP is really being enforced, that M* is sane, and
    # that the DR path reproduces the state the static path handed it.  If DR
    # drifts here, every number below is void.
    # It is deliberately run for a FIXED, long hold rather than stopped the
    # moment it looks quiet: `-mass gershgorin` sizes M* so that
    # omega_max*dt <= 2 by the Gershgorin bound, and for an FE stiffness that
    # bound is very nearly ATTAINED by the highest mode -- i.e. the march sits
    # on the central-difference stability boundary, where round-off is
    # amplified rather than damped.  A short hold cannot see that; a long one
    # can.  This is the control that made the failure visible at all.
    set_push(foot, uz0, first=True)
    t0 = time.time()
    hold = max(args.dr0_steps, args.chunk)
    n_h, ke_h, res_h = 0, 0.0, 0.0
    while n_h < hold:
        if ops.analyze(args.chunk, args.dt_march) != 0:
            log("    [control DR-0] *** FAIL: the zero-push hold DIVERGED ***")
            raise SystemExit("DR-0 control failed (analyze)")
        n_h += args.chunk
        ke_h = float(ops.ladrunoDR("kineticEnergy"))
        res_h = float(ops.ladrunoDR("residualNorm"))
        if args.trace:
            log(f"        [DR-0 hold] {n_h:>7d} steps  KE = {ke_h:.4e}  "
                f"res = {res_h:.4e} kN  R = {reaction(foot, r_corr):.6f} kN")
    log(f"    [control DR-0] held {n_h} steps at zero push: KE = {ke_h:.3e}, "
        f"res = {res_h:.3e} kN")
    ok0, n0, R0, ke0, res0, why0, _ = relax(args, foot, r_corr, res_abs)
    n0 += n_h
    drift_u = abs(ops.nodeDisp(foot[0], 3) - uz0)
    drift_R = abs(R0 - R_static0) / max(abs(R_static0), 1.0)
    log(f"    [control DR-0] zero-push relaxation: {n0} DR steps, {why0}; "
        f"footing |du| = {drift_u:.3e} m, R = {R0:.6f} vs static "
        f"{R_static0:.6f} kN ({100 * drift_R:+.6f} %), KE = {ke0:.3e}, "
        f"res = {res0:.3e} kN   [{time.time() - t0:.1f}s]")
    if not ok0 or drift_u > 1.0e-9 or drift_R > 1.0e-4 or res_h > res_abs:
        log("    [control DR-0] *** FAIL: DR does not reproduce the static "
            "state it was handed; every number below would be void ***")
        raise SystemExit("DR-0 control failed")
    log("    [control DR-0] PASS -- DR holds the handed-over static equilibrium")

    # --- the push ---------------------------------------------------------
    smax = args.sfrac * H.B_FOOT
    if args.mode == "continuous":
        return push_continuous(args, foot, r_corr, uz0, res_abs, smax, area,
                               q_exact, nodes, cells)
    ninc = int(round(smax / args.ds))
    log(f"=== leg {args.elem} h0 = {args.h0}: {ninc} held increments of "
        f"ds = {args.ds * 1e3:.4f} mm to s/B = {args.sfrac}, each RAMPED in "
        f"{args.ramp} sub-jumps of {args.ds / args.ramp * 1e6:.2f} um every "
        f"{args.rampevery} DR steps, then relaxed (<= {args.maxit} steps)")
    log(f"    ORACLE q_u = {q_exact:.3f} kPa")

    tag = f"dr_{args.elem}_h{str(args.h0).replace('0.', '')}{args.suffix}"
    path = os.path.join(HERE, f"{tag}.csv")
    fh = open(path, "w", newline="")
    wr = csv.writer(fh)
    wr.writerow(["s_m", "s_over_B", "R_kN", "q_kPa", "dr_iters", "settled",
                 "ke", "res_kN", "dR_kN", "wall_s"])
    rows, verdict, mode = [], "reached the target", "TARGET"
    t0 = time.time()
    for k in range(1, ninc + 1):
        s = k * args.ds
        if time.time() - t0 > args.tmax:
            verdict = f"WALL-CLOCK CAP of {args.tmax:.0f} s hit"
            mode = "WALL"
            break
        if not ramp_push(args, foot, uz0 - (k - 1) * args.ds, uz0 - s):
            verdict, mode = f"increment {k} DIVERGED during the ramp", "DIVERGED"
            break
        ok, nit, R, ke, res, why, dR = relax(args, foot, r_corr, res_abs)
        nit += max(0, args.ramp - 1) * args.rampevery
        rows.append((s, s / H.B_FOOT, R, R / area, nit, int(ok), ke, res,
                     dR, time.time() - t0))
        wr.writerow([f"{v:.9g}" for v in rows[-1]])
        fh.flush()
        if k % args.report == 0 or not ok:
            log(f"    inc {k:>4}/{ninc}  s/B = {s / H.B_FOOT:.5f}  "
                f"q = {R / area:8.3f} kPa ({R / area / q_exact:.4f} of exact)  "
                f"{nit:>6} DR steps  {'settled' if ok else 'NOT SETTLED: ' + why}"
                f"  KE = {ke:.2e}  res = {res:.2e}  [{time.time() - t0:.0f}s]")
        if not ok:
            verdict = f"increment {k} did not settle: {why}"
            mode = "DIVERGED" if "failed" in why or "finite" in why else "ITERCAP"
            break
    fh.close()
    wall = time.time() - t0
    if not rows:
        raise SystemExit("no increments completed")

    a = np.array([r[:4] for r in rows])
    s, q = a[:, 0], a[:, 3]
    its = np.array([r[4] for r in rows], dtype=float)
    peak = np.maximum.accumulate(q)
    bad = np.where(q < (1.0 - DROP_TRUNCATE) * peak)[0]
    if len(bad):
        ntr = len(q) - bad[0]
        s, q, its = s[:bad[0]], q[:bad[0]], its[:bad[0]]
        verdict += f"; curve TRUNCATED at s/B = {s[-1] / H.B_FOOT:.4f} ({ntr} cut)"
        mode = "TRUNCATED"

    msk = s >= 0.9 * s[-1]
    t_last = np.polyfit(s[msk], q[msk], 1)[0] if msk.sum() > 2 else float("nan")
    n0f = max(4, len(s) // 50)
    t_init = np.polyfit(s[:n0f], q[:n0f], 1)[0]
    plateau = bool(abs(t_last) < 0.02 * abs(t_init))
    tail_it = float(its[msk].max()) if msk.sum() else float(its.max())
    free = bool(tail_it <= FREE_ADVANCE_FRACTION * args.maxit)
    qmax = float(q.max())
    capacity = bool(plateau and mode == "TARGET" and free)

    log(f"    {len(rows)} increments, {int(its.sum())} DR steps total, "
        f"{wall:.0f} s   [{verdict}]")
    log(f"    q_max = {qmax:.2f} kPa at s/B = "
        f"{s[int(q.argmax())] / H.B_FOOT:.4f}; end s/B = {s[-1] / H.B_FOOT:.4f}")
    log(f"    tail dq/ds = {t_last:.2f} kPa/m = {100 * t_last / t_init:.3f} % "
        f"of initial  ->  {'PLATEAU' if plateau else 'STILL HARDENING'}")
    log(f"    MODE = {mode}; worst DR spend over the final tenth = "
        f"{tail_it:.0f}/{args.maxit} "
        f"({'free' if free else 'AT THE ALLOWANCE'})")
    log(f"    q_num / q_exact = {qmax / q_exact:.4f}   -->  "
        f"{'CAPACITY' if capacity else 'CONTROLLER ALLOWANCE (' + mode + ')'}")
    return dict(leg=tag, elem=args.elem, form=form, h0=args.h0,
                ndof=3 * len(nodes), nel=len(cells), qmax=qmax,
                q_exact=q_exact, ratio=qmax / q_exact, plateau=plateau,
                t_frac=100 * t_last / t_init, s_end=float(s[-1]) / H.B_FOOT,
                ninc=len(rows), iters=int(its.sum()), tail_it=tail_it,
                mode=mode, free=free, capacity=capacity, wall=wall,
                verdict=verdict, csv=path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--elem", default="h8bbar", choices=sorted(ELEMS))
    ap.add_argument("--h0", type=float, default=0.5)
    ap.add_argument("--sfrac", type=float, default=0.05,
                    help="s/B target (the linear control plateaus well before "
                         "0.05; note 82's h8bbar_h5 ran to 0.0744)")
    ap.add_argument("--ds", type=float, default=1.0e-3,
                    help="held displacement increment, m")
    ap.add_argument("--dt", type=float, default=1.0,
                    help="the integrator's -dt, i.e. the dt the fictitious "
                         "Gershgorin mass is SIZED for.  With -mass gershgorin "
                         "only the RATIO dtfac = dt_march/dt matters (see the "
                         "module docstring)")
    ap.add_argument("--dtfac", type=float, default=0.5,
                    help="dt_march / dt.  Gershgorin sizes M* so that "
                         "omega_max*dt <= 2 -- the central-difference stability "
                         "BOUNDARY, where the amplification matrix is defective "
                         "and the march grows linearly.  dtfac < 1 is the "
                         "safety factor: it makes omega_max*dt_march <= 2*dtfac "
                         "and scales the per-step relaxation by dtfac^2")
    ap.add_argument("--mass", default="gershgorin",
                    choices=["gershgorin", "lumped", "unity"])
    ap.add_argument("--mass-scale", type=float, default=1.0)
    ap.add_argument("--damping", default="kinetic",
                    choices=["kinetic", "viscous"])
    ap.add_argument("--zeta", type=float, default=1.0)
    ap.add_argument("--no-auto-refresh", action="store_true")
    ap.add_argument("--recompute", type=int, default=0)
    ap.add_argument("--chunk", type=int, default=500)
    ap.add_argument("--maxit", type=int, default=40000)
    ap.add_argument("--tol-r", type=float, default=1.0e-5,
                    help="REPORTED only (see relax()); does not gate settling")
    ap.add_argument("--tol-res", type=float, default=1.0e-6,
                    help="||f_ext-f_int||_inf as a fraction of the total "
                         "surcharge (300 kN), matching the static ladder's "
                         "absolute NormUnbalance tolerance")
    ap.add_argument("--tmax", type=float, default=1.0e9)
    ap.add_argument("--ramp", type=int, default=50,
                    help="sub-jumps the held displacement increment is applied "
                         "in.  1 = step it (measured to leave 29 % of spurious "
                         "plastic softening behind -- see ramp_push())")
    ap.add_argument("--rampevery", type=int, default=40,
                    help="DR steps between ramp sub-jumps")
    ap.add_argument("--mode", default="held", choices=["held", "continuous"])
    ap.add_argument("--rate", type=float, default=5.0e-7,
                    help="continuous mode: imposed footing displacement per DR "
                         "step, m.  This is THE controlled quantity for "
                         "quasi-staticness (see ramp_push / sec 4.2)")
    ap.add_argument("--holdevery", type=int, default=0,
                    help="continuous mode: hold and relax every N samples, to "
                         "measure the inertial lag directly")
    ap.add_argument("--report", type=int, default=5)
    ap.add_argument("--patch", dest="patch", action="store_true", default=True)
    ap.add_argument("--no-patch", dest="patch", action="store_false")
    ap.add_argument("--modes", action="store_true")
    ap.add_argument("--trace", action="store_true")
    ap.add_argument("--dr0-steps", type=int, default=3000,
                    help="length of the zero-push hold in control DR-0.  Long "
                         "on purpose: the marginal-stability failure it exists "
                         "to catch takes thousands of steps to surface")
    ap.add_argument("--suffix", default="")
    args = ap.parse_args()
    args.dt_march = args.dt * args.dtfac

    try:
        stamp = ops.ladrunoBuild()
    except Exception as exc:                                  # pragma: no cover
        raise SystemExit(f"ladrunoBuild() unavailable ({exc}) -- refusing to "
                         f"measure")
    log(f"[control 0] engine {ops.__file__}")
    log(f"[control 0] ladrunoBuild() = {stamp}")

    r = run_leg(args)
    print()
    print(f"{'leg':>26} {'DOF':>6} {'q_num':>9} {'num/exact':>10} "
          f"{'tail %':>8} {'plateau':>8} {'s_end/B':>8} {'MODE':>9} "
          f"{'free':>5} {'CAPACITY':>9} {'DRsteps':>9} {'wall_s':>7}")
    print(f"{r['leg']:>26} {r['ndof']:>6} {r['qmax']:9.2f} {r['ratio']:10.4f} "
          f"{r['t_frac']:8.3f} {'yes' if r['plateau'] else 'NO':>8} "
          f"{r['s_end']:8.4f} {r['mode']:>9} "
          f"{'yes' if r['free'] else 'NO':>5} "
          f"{'YES' if r['capacity'] else 'no':>9} {r['iters']:9d} "
          f"{r['wall']:7.0f}")


if __name__ == "__main__":
    main()
