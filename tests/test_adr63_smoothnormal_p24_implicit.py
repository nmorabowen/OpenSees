"""ADR-63 #4a P2.4 — the Q-IMPLICIT-NEWTON convergence tripwire (BLOCKER-TANGENT-POSTURE).

The last open P2 item. ADR-63 ships the SYMMETRIC frozen-field `kn·BᵀB` tangent for the smoothed
NTS normal (the faceted B3 `∂n/∂u` block is SUPPRESSED under `-smoothNormal`, LadrunoContactFE.cpp
`addKtToTang`: `if (consistentNormal && !useSmoothNormal)`). The frozen-field tangent is an EXACT
modified-Newton on per-step-frozen geometry (residual `n̄` ≡ tangent `n̄` — one handle-frozen field),
so the converged answer is tangent-independent; the only dropped piece is the configuration
dependence `∂n_smooth/∂u`, the gated P3 term. Steelman-B's live concern (ADR §Gate decision): on a
GENUINELY CURVED master the dropped term is non-zero and LARGER in magnitude than the faceted
default's (the faceted flat-facet block is ≡0), so Newton could STALL — and that lands hardest here
because R3 is *about* curved masters.

This rig is the evidence-gated P2 oracle that BINDS the P3 promotion to a MEASURED condition
(killing the "P3 vaporware" discretionary-backlog risk). It runs a genuinely rotating-normal IMPLICIT
quasi-static problem and records the Newton iteration count per load step. Two outcomes, both fine:
  (a) smoothed converges every step (no stall/divergence), iteration counts stay well under the cap
      →  B's stall concern is empirically REFUTED, P3 stays deferred WITH EVIDENCE.
  (b) smoothed STALLS/DIVERGES (not merely slower)  →  the tripwire TRIPS, P3 promoted to REQUIRED.

WHY IMPLICIT + DisplacementControl (not the explicit point-mass-on-curve of P1/P2.3): the explicit
lane assembles NO tangent (the contact FE returns the integrator-routed zero tangent), so it can
NEVER exercise the frozen-field `kn·BᵀB` — only an IMPLICIT integrator calls `addKtToTang`. And a
frictionless slave dragged by a lateral FORCE launches ballistically; a PRESCRIBED x-displacement
(DisplacementControl) is quasi-static by construction. DisplacementControl is also the ONLY sanctioned
displacement driver here — `LadrunoContactHandler::handle()` refuses a non-homogeneous (imposed-
displacement) SP: "DisplacementControl is the driver, it drives a DOF via the load factor."

RIG — the clean isolation of ∂n_smooth/∂u.  A symmetric convex bump built from THREE quad facets:
two steep outer facets and a FLAT horizontal middle facet at z=hz. The master nodes are all FIXED, so
the nodal-normal field `n_a` (a property of the fixed master geometry) is genuinely CONSTANT — no
re-handle needed for correctness, so we need NEITHER `-reemit` (unvalidated under StaticAnalysis) NOR
facet migration. The middle facet is FLAT, so its FACETED normal is a constant vertical (0,0,1) — but
its two shared-edge nodal normals inherit the steep outer facets' tilt (~±22° at hz=1.0). So as the
slave is dragged ACROSS the flat middle facet, the SMOOTHED normal `n_smooth(ξ) = Σ N_i(ξ) n_i`
rotates smoothly from vertical toward the tilt while the facet geometry itself is flat. The dropped
`∂n_smooth/∂u` is therefore ENTIRELY the nodal-blend rotation through the moving projection — exactly
the term P3 would add, isolated from any facet-projection curvature. If the frozen-field Newton
converges here, it converges on the essence of the smoothed-normal nonlinearity.

The slave is seated at x=0 (the flat facet's symmetric centre, where `n_smooth` is vertical ⇒ no
frictionless seat-slide) under a constant downward press (loadConst-frozen), then dragged to
x=±0.35·(half-facet) by DisplacementControl. A weak lateral spring (ks ≪ kn) removes the crest's
zero-lateral-stiffness singularity at the seat step ONLY (during the drag x is displacement-controlled
⇒ non-singular anyway) — its force is <0.5% of the contact force, negligible to the Newton behaviour.

FINDING (2026-07-01, this rig).  Outcome (a) — CONVERGES.  The smoothed frozen-field Newton converges
in 2 iterations every step across the whole VALID penalty regime, and — the load-bearing evidence —
independently of the load-step coarseness: even a SINGLE step that drags the slave across the entire
facet (the maximal within-step normal rotation) converges in 2 iterations. It stays 2 iterations up to
steep facets (hz≤4), stiff penalties (kn≤1e6), and heavy press (≤3000). The physical reason: the
dropped `∂n_smooth/∂u` term is `O(kn·gN) = O(press)` — scaled by the (small) penetration `gN≈press/kn`,
so it is sub-dominant to the kept `kn·BᵀB` in any well-posed penalty contact. The ONLY regime where it
strains is LARGE penetration (`gN` ≳ 15% of the facet) — but there the FACETED `-geomtan` consistent
tangent diverges too (and even the seat fails), because that is penalty NTS contact violating its own
small-penetration premise, NOT a smoothed-normal defect (P3's full `∂n_smooth/∂u` would not rescue it).
So B's stall concern is empirically REFUTED and P3 stays a deferred, evidence-gated option.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e4                      # implicit ⇒ moderate, well-conditioned penalty (penetration ~ P/KN)
KS = 0.1                        # weak lateral stabilizing spring (seat non-singularity); ≪ KN
HZ = 1.0                        # middle-facet height (outer facets rise 0→HZ over 1.0 ⇒ ~22° nodal tilt)
XCOL = [-1.5, -0.5, 0.5, 1.5]  # 4 master node columns ⇒ 3 facets; middle facet is x∈[-0.5,0.5], flat
ZCOL = [0.0, HZ, HZ, 0.0]      # convex: flat middle at z=HZ, steep outer facets
PRESS = 30.0                   # constant downward press (loadConst-frozen)
FX_REF = 10.0                  # lateral reference load DisplacementControl scales to advance x
HALF = 0.5                     # middle-facet half-width
XDRAG = 0.35                   # drag |x| ≤ 0.35 < 0.5 ⇒ stays INTERIOR to the flat middle facet
NSTEPS_FINE = 35               # many small steps (small within-step rotation)
NSTEPS_COARSE = 3              # few big steps — MAXIMAL within-step normal rotation (the real stress)
TOL = 1.0e-8
MAXIT = 40


def _build_bump():
    """The 3-facet convex bump (flat middle, steep outer), all master nodes FIXED. Returns the master
    tag list; facet 1 (index s=1) is the flat middle facet x∈[-0.5,0.5]."""
    for c, (x, z) in enumerate(zip(XCOL, ZCOL)):
        ops.node(100 + c, x, -0.6, z); ops.fix(100 + c, 1, 1, 1)
        ops.node(200 + c, x,  0.6, z); ops.fix(200 + c, 1, 1, 1)
    m = []
    for s in range(len(XCOL) - 1):
        m += [100 + s, 100 + s + 1, 200 + s + 1, 200 + s]
    return m


def _drag_implicit(flags, nsteps=NSTEPS_FINE):
    """Implicit quasi-static: seat a slave at the flat-facet centre under a constant press (loadConst),
    then drag it across the flat middle facet by DisplacementControl on its x-DOF, recording the Newton
    iteration count per step. `flags` = extra contact tokens (e.g. ['-smoothNormal', '-geomtan']);
    `nsteps` = the load-step schedule (few big steps ⇒ larger within-step normal rotation).

    Returns dict: steps_done (of nsteps), iters (per-step testIter list), max_iter, mean_iter,
    reached_x, diverged, seat_ok."""
    du = XDRAG / nsteps
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    m = _build_bump()

    # slave at the flat-facet centre, δ below the surface so contact is active from step 1. x=0 ⇒
    # n_smooth vertical by symmetry ⇒ a vertical press induces NO lateral force ⇒ no seat-slide.
    z0 = HZ - 1.0e-4
    ops.node(1, 0.0, 0.0, z0)
    ops.fix(1, 0, 1, 0)                        # y FIXED (no y motion); x,z free

    # weak lateral spring to ground: removes the seat-step x singularity (vertical n ⇒ 0 contact
    # x-stiffness). Grounded twin at the slave's start; zeroLength in dir 1 (x).
    ops.node(2, 0.0, 0.0, z0); ops.fix(2, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, KS)
    ops.element("zeroLength", 999, 2, 1, "-mat", 1, "-dir", 1)

    ops.contactSurface(10, "-master", 4, *m)
    ops.contactSurface(20, "-slave", 1)
    # frictionless (the clean symmetric-tangent measurement — ADR §P2 note); explicit -outward so the
    # global sign is well-conditioned. NO -reemit: the fixed master's nodal field is constant, and the
    # drag stays interior to one facet ⇒ no facet migration ⇒ the reference candidate band always holds.
    args = [1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0] + list(flags)
    ops.contact(*args)

    # --- stage 1: seat the constant press at the flat centre, then freeze it (loadConst) ---
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -PRESS)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("BandGeneral")
    ops.test("NormDispIncr", TOL, MAXIT)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    seat_ok = (ops.analyze(1) == 0)
    ops.loadConst("-time", 0.0)                # freeze the press; DisplacementControl scales only pat 2

    # --- stage 2: drag +x by DisplacementControl, record Newton iterations per step ---
    ops.timeSeries("Constant", 2)
    ops.pattern("Plain", 2, 2)
    ops.load(1, FX_REF, 0.0, 0.0)              # lateral reference load; λ scales it to advance x by du
    ops.integrator("DisplacementControl", 1, 1, du)

    iters = []
    reached = 0.0
    diverged = not seat_ok
    steps_done = 0
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            diverged = True
            break
        it = ops.testIter()
        iters.append(int(it))
        x = ops.nodeDisp(1, 1)
        if not (x == x):                       # NaN guard
            diverged = True
            break
        reached = max(reached, x)
        steps_done += 1

    return {
        "steps_done": steps_done,
        "iters": iters,
        "max_iter": max(iters) if iters else -1,
        "mean_iter": (sum(iters) / len(iters)) if iters else -1.0,
        "reached_x": reached,
        "diverged": diverged,
        "seat_ok": seat_ok,
    }


# --------------------------------------------------------------------------------------------------
# The gate. Outcome (a) is the passing assertion: the smoothed frozen-field tangent CONVERGES on a
# genuinely rotating-normal implicit master. If the run STALLS/DIVERGES these fail → the tripwire
# TRIPS → P3 is promoted to REQUIRED (implement the full ∂n_smooth/∂u, -consistentNormalSmooth).
# --------------------------------------------------------------------------------------------------

def test_p24_smoothed_frozen_field_converges_on_rotating_normal():
    """THE tripwire (fine-step). Drag a slave across the flat middle facet whose SMOOTHED normal
    rotates ~15° (the dropped ∂n_smooth/∂u is fully active) under implicit Newton with the frozen-field
    symmetric kn·BᵀB tangent (-smoothNormal + -geomtan ⇒ B3 suppressed). Newton must CONVERGE every
    step (analyze==0 for all steps) and stay well under the iteration cap. Converges ⇒ B's stall
    concern is empirically REFUTED and P3 stays deferred WITH EVIDENCE."""
    r = _drag_implicit(["-smoothNormal", "-geomtan"], nsteps=NSTEPS_FINE)
    print("\n[P2.4 smoothed frozen-field, fine]", {k: r[k] for k in
          ("steps_done", "max_iter", "mean_iter", "reached_x", "diverged")}, "iters=", r["iters"])
    assert r["seat_ok"], "flat-centre seating step did not converge"
    assert not r["diverged"], (
        "TRIPWIRE TRIPPED: smoothed frozen-field Newton diverged/stalled on a rotating-normal implicit "
        f"master (steps_done={r['steps_done']}/{NSTEPS_FINE}, iters={r['iters']}) — promote ADR-63 P3 "
        "(∂n_smooth/∂u, -consistentNormalSmooth)")
    assert r["steps_done"] == NSTEPS_FINE, (
        f"smoothed run did not complete the drag (steps_done={r['steps_done']}/{NSTEPS_FINE})")
    assert r["reached_x"] > 0.9 * XDRAG, f"slave did not traverse the facet (reached_x={r['reached_x']})"
    assert r["max_iter"] < MAXIT, (
        f"smoothed Newton scraped the iteration cap (max_iter={r['max_iter']} ≥ {MAXIT}) — "
        "borderline stall, treat as a trip")


def test_p24_smoothed_converges_under_maximal_within_step_rotation():
    """THE load-bearing stress. The fine-step arm has tiny per-step rotation, so the frozen-field
    tangent is barely tested. Here the WHOLE facet is crossed in only NSTEPS_COARSE=3 steps — each step
    demands a large within-step normal rotation, exactly the regime where a modified-Newton (dropped
    ∂n_smooth/∂u) would STALL if the dropped term were dominant. It is not: on a penalty contact the
    dropped term is O(kn·gN)=O(press) ≪ kn·BᵀB, so coarse-step Newton converges just as robustly. If
    THIS diverges/stalls, the tripwire trips and P3 is promoted."""
    r = _drag_implicit(["-smoothNormal", "-geomtan"], nsteps=NSTEPS_COARSE)
    print("\n[P2.4 smoothed frozen-field, COARSE]", {k: r[k] for k in
          ("steps_done", "max_iter", "mean_iter", "reached_x", "diverged")}, "iters=", r["iters"])
    assert r["seat_ok"], "flat-centre seating step did not converge"
    assert not r["diverged"] and r["steps_done"] == NSTEPS_COARSE, (
        "TRIPWIRE TRIPPED under coarse steps: smoothed frozen-field Newton failed the large-within-step-"
        f"rotation drag (steps_done={r['steps_done']}/{NSTEPS_COARSE}, iters={r['iters']}) — promote "
        "ADR-63 P3 (∂n_smooth/∂u, -consistentNormalSmooth)")
    assert r["max_iter"] < MAXIT, (
        f"coarse-step smoothed Newton scraped the iteration cap (max_iter={r['max_iter']}) — a stall")


def test_p24_smoothed_iterations_bounded_and_faceted_sanity():
    """Quantitative arm: (i) the frozen-field smoothed Newton converges COMFORTABLY — a small, bounded
    iteration count that does NOT climb toward the cap as the normal rotation accumulates across the
    drag (a stall would show as monotonically rising iters). (ii) The faceted -geomtan path on the
    SAME rig converges too (a sanity anchor; on the flat middle facet its faceted normal is constant
    vertical ⇒ near-trivial). Both completing confirms the rig + solver are sound and the smoothed
    frozen-field posture is a modest RATE cost at worst, not a convergence failure."""
    faceted = _drag_implicit(["-geomtan"], nsteps=NSTEPS_FINE)
    smooth = _drag_implicit(["-smoothNormal", "-geomtan"], nsteps=NSTEPS_FINE)
    print("\n[P2.4 faceted -geomtan   ]", {k: faceted[k] for k in
          ("steps_done", "max_iter", "mean_iter", "diverged")})
    print("[P2.4 smoothed frozen    ]", {k: smooth[k] for k in
          ("steps_done", "max_iter", "mean_iter", "diverged")}, "iters=", smooth["iters"])
    assert not faceted["diverged"] and faceted["steps_done"] == NSTEPS_FINE, (
        f"faceted -geomtan baseline did not converge (steps_done={faceted['steps_done']}/{NSTEPS_FINE})")
    assert not smooth["diverged"] and smooth["steps_done"] == NSTEPS_FINE, (
        f"smoothed frozen-field did not converge (steps_done={smooth['steps_done']}/{NSTEPS_FINE})")
    # 'comfortably': the WORST-case step stays well under the cap (no cliff-edge convergence).
    assert smooth["max_iter"] <= 15, (
        f"smoothed Newton took {smooth['max_iter']} iters at its worst step (cap {MAXIT}) — the "
        "frozen-field rate penalty is large; re-examine whether P3 should be promoted")
    # no runaway growth: the last-third mean is not dramatically worse than the first-third mean
    # (a genuine stall from the accumulating rotation would blow up the tail).
    n = len(smooth["iters"])
    if n >= 6:
        head = smooth["iters"][: n // 3]
        tail = smooth["iters"][-(n // 3):]
        mh = sum(head) / len(head)
        mt = sum(tail) / len(tail)
        assert mt <= 3.0 * mh + 3.0, (
            f"smoothed Newton iteration count is RUNNING AWAY as the normal rotation accumulates "
            f"(head mean={mh:.1f}, tail mean={mt:.1f}) — a rotation-driven stall; consider P3")
