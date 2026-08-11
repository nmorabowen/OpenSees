"""R3 — Prandtl–Reissner strip footing: the fork's exact-answer collapse gate.

SLOW TIER (`@pytest.mark.slow`, opt-in via `--runslow` / `LADRUNO_RUN_SLOW=1`).
MEASURED WALL TIME: see `_WALL_TIME_NOTE` below — this file is minutes, not
seconds, and must never be pulled into the default Zone-A battery.

WHAT THIS GATES
---------------
A strip footing on WEIGHTLESS soil (gamma = 0) under a surcharge q0 is the
Prandtl–Reissner problem.  Its limit-analysis upper and lower bounds COINCIDE at

    q_u = q0 * N_q,   N_q = e^(pi tan phi) tan^2(45 + phi/2)

so there is no shape factor, no N_gamma, no empirical fit and nothing to
calibrate: an EXACT oracle, owed to nothing in this repository.  The gate drives
it with `LadrunoBrick -formulation bbar` — the one configuration with an
established fork baseline (`Ladruno_files/testbed/hypo_bearing/COLLAPSE.md`,
leg `nq20_nonassoc`, 1.0020 of exact) — at three resolutions, and asks the three
questions a collapse-load measurement actually has to answer.

Credit: the deck, the guards and both calibration constraints below come from
the TIMs / APE campaign (`ADR 63` R3; vault note 71, "The Exact Answer the
Tetrahedron Cannot Reach").  Reference drivers: `r3_prandtl_tet10.py` and
`r3_hex_control.py` in that campaign's fork bundle.  This file reproduces the
hex b-bar control leg of `r3_hex_control.py` on the same domain, cone,
surcharge and push as the fork's own `hypo_bearing/dp_strip.py`.

CONSTRAINT 1 — PER-RESOLUTION BANDS, NEVER A FLAT PERCENTAGE
------------------------------------------------------------
The measured hex-bbar sequence is

    h0 = 1.00  ->  1.0842      (200 hexes,   1 386 DOF)
    h0 = 0.50  ->  0.9938      (390 hexes,   2 592 DOF)
    h0 = 0.25  ->  0.9513      (782 hexes,   5 040 DOF)

each one a GENUINE PLATEAU.  That is textbook convergence FROM ABOVE for
displacement-based elements at a limit load: too few DOF to form the true
velocity discontinuity, so the cheapest available mechanism is stiffer and
stronger, and refinement relaxes it.  A flat 5 % accuracy gate would FAIL the
coarse mesh for being correct.  So each resolution carries its own band, and
the real assertion is the SHAPE of the sequence — monotone decrease, every leg
on a plateau — not proximity to the exact answer at any one rung.

CONSTRAINT 2 — THE SUBDIVISION BUDGET IS PINNED, NOT INHERITED
--------------------------------------------------------------
`SUBDIV_BUDGET` below is load-bearing.  Note 71 records that at ADR-63 D16's
default of 24 the h0 = 0.25 leg read **0.8988, "still hardening"**, terminal
tangent 5.75 % of initial — which reads as a sequence turning over, i.e. a real
loss of reach.  Raising the budget to 80 and changing nothing else, the same
mesh runs on to a clean plateau at **0.9513**, terminal tangent 0.02 %.  The
0.899 was the stepping controller, not the mesh.  A convergence study whose
termination criterion binds before the mechanics do measures the criterion, so
the budget is pinned here as a named constant with that history attached; a
gate that inherited a default would be silently measuring its own guard.

WHAT MAKES A NUMBER A CAPACITY (three clauses, not one)
--------------------------------------------------------
    A leg is a CAPACITY only if it reached the target with a flat tangent.
    Otherwise it is a CONTROLLER ALLOWANCE, and the allowance must be named.

A flat tail tangent alone is NOT enough, because a curve can flatten for two
completely different reasons: a mechanism formed, or **the run seized**.  A step
floor reached just as the curve rolls over, or a stall firing on what is really
numerical seizure, both present as a flat tail.  So:

    A flat tangent is a capacity only if the leg was still ADVANCING FREELY
    when it flattened.

Every leg is therefore classified on three clauses, and `capacity` requires all
three (`_classify`):

  1. **plateau** — tail dq/ds < 2 % of the initial tangent;
  2. **termination mode** — `TARGET` (reached s/B) or `BUDGET` (spent the pinned
     subdivision budget).  `FLOOR` (step collapsed to DS_MIN), `WALL`
     (wall-clock) and `TRUNCATED` are SEIZURE modes and are never a capacity;
  3. **free advance** — the step size over the final tenth of the run stayed at
     least `FREE_ADVANCE_FLOOR_FACTOR` x DS_MIN above the floor.

`BUDGET` is a capacity WITH A NAMED ALLOWANCE, not a silent pass: the reference
records the h0 = 0.5 control as having "plateaued first and then been halted by
D16's subdivision budget ... the load had been flat for the last tenth of the run
... the guards simply declined to spend an hour confirming it".  Every failure
message here prints the mode, the terminal and tail-minimum step size against the
floor, and the subdivisions used against the pinned budget — so a future reader of
a red gate can tell "the element got worse" from "the controller ran out of room",
which is the exact confusion this rule exists to prevent.

The clauses separate by three orders of magnitude on this problem (measured):
the three gate legs sit at tail 0.001 / 0.003 / 0.016 % with a tail step 2500 /
1250 / 800 x the floor; the associated control sits at tail 39 % with a tail step
1.6 x the floor — it is seizing, not plateauing.

THE FALSIFICATION CONTROL
-------------------------
The gate leg is NON-ASSOCIATED (psi = 0).  The ASSOCIATED leg is run as the
control and must NOT be a CAPACITY: dilating at psi = phi has to push the
surrounding soil aside, a bounded mesh resists it, and the leg hardens straight
past its own exact answer (COLLAPSE.md: `nq20_assoc` = 1.1171, no plateau).  If
the associated leg ALSO produced a capacity, the non-associated agreement would
be a coincidence rather than a measurement, and this file says so.

THE SURCHARGE STEP IS CHECKED TWICE, ON PURPOSE
-----------------------------------------------
The reference driver validates the surcharge with the resultant identity
`sum(R_z) == q0 * sum(A_trib)`.  That identity is kept here — and it is
**structurally blind to load-DISTRIBUTION error**: a tributary lumping can
reproduce its applied total to +0.0000000 % while putting O(1) error into
sigma_zz.  So it is paired with a **1-D elastic stress patch check**: this
domain under a uniform surcharge on a weightless, laterally-rollered,
base-fixed box has the exact closed-form answer sigma_zz = -q0 and
sigma_xx = sigma_yy = -K0 q0 at EVERY Gauss point, with K0 = nu/(1-nu).  The two
checks discriminate totally — the identity cannot fail on a distribution error,
the patch cannot pass one.

MEASURED WALL TIME
------------------
**61.1 min** (3668 s) for the whole module, 11 tests, one process, on the
reference box (Windows 11, Intel oneAPI build, `system Pardiso`, MKL default
threading): 537 s (h0 = 1.0) + 517 s (h0 = 0.5) + 1566 s (h0 = 0.25) for the
non-associated sequence, plus 1047 s for the associated control.  That is why
this file is slow-tier and opt-in; do not let it into the default battery.

MEASURED, on the build stamped in the output (`ladrunoBuild`
1db3394b8c4d4c23e12f2629d5e7ef92a106acfa):

    leg               DOF     q_num   ratio   tail %  mode    ds/floor  CAP
    h1.0_nonassoc    1386    150.71  1.0849    0.001  BUDGET      2500  yes
    h0.5_nonassoc    2592    138.58  0.9977    0.003  BUDGET      1250  yes
    h0.25_nonassoc   5040    132.16  0.9514    0.017  BUDGET       800  yes
    h0.5_assoc       2592    154.73  1.1139   38.933  BUDGET       1.6   NO   <- control

against the reference 1.0842 / 0.9938 / 0.9513 (deltas +0.06 / +0.39 / +0.01 %)
and the reference's associated 1.1171.  Surcharge: resultant identity 2.7e-15 to
8.9e-15, 1-D stress patch 1.0e-14 to 1.4e-13.

Note every gate leg terminates on `BUDGET`, not `TARGET` — a NAMED controller
allowance, and the reason clause 3 carries the weight here: all three stop with
a dead-flat tail while their step is still 800-2500x the floor, whereas the
control stops at 1.6x with a 39 % tail.  Reaching s/B = 0.15 is not affordable
at this budget and MUST NOT be bought by lowering `SUBDIV_BUDGET` (that is the
note-71 trap) or by moving the target down to meet the run (that is the same
trap wearing a different hat).
"""

import csv
import math
import os
import re
import time

import numpy as np
import pytest

from _testbed import ops


pytestmark = [pytest.mark.zone_a, pytest.mark.slow, pytest.mark.t2a]


# --- the problem (identical to hypo_bearing/dp_strip.py `nq20_*`) -----------
PHI_TXC = 20.0        # the MILD cone: mesh-resolvable, and where the
                      # Prandtl/Vesic coefficients were fitted
Q0 = 10.0             # kPa surcharge
K_EL = 1.5e5          # kPa
NU = 0.45             # free parameter: a collapse load of an elastic-perfectly-
                      # plastic body does not depend on the elastic constants,
                      # but at a LOW nu this cone starts at 95 % of yield (the
                      # 1-D state m0 = (1-K0)/(sqrt3 alpha (1+2 K0))) and the leg
                      # would measure nothing but its own initial condition.
                      # That is how the fork's first `verify_phi20` leg was void.
SY = 0.2              # kPa apex regularizer, NOT a soil property
B_FOOT, THICK = 2.0, 0.5
XLIM, ZBOT = 30.0, -20.0     # 14.5 B of clearance, 10 B of depth
R_GRADE = 1.35
SFRAC = 0.15          # push to s/B; the plateau arrives well inside this

# --- ADR-63 D16 adaptive-stepping guards ------------------------------------
# The plain halve-on-failure / double-on-success controller livelocked the
# reference campaign for 17 hours pinned at one settlement.  These are the
# guards that made it terminate WITH A VERDICT.
DS_BASE, DS_MIN, DS_MAX = 2.0e-5, 2.0e-7, 1.0e-3
GROW_AFTER = 6        # recovery latch: 6 clean steps before the step grows

# *** CONSTRAINT 2 — PINNED, see the module docstring. ***
# D16 ships 24.  At 24 the h0 = 0.25 leg reads 0.8988 STILL HARDENING (tail
# 5.75 % of initial); at 80 the same mesh plateaus at 0.9513 (tail 0.02 %).
# Do not "clean this up" into a default, and do not lower it to make the gate
# faster: the number below is the difference between measuring the element and
# measuring the guard.
SUBDIV_BUDGET = 80

DROP_TRUNCATE = 0.02  # cut the curve where q first falls >2 % below its max —
                      # a bottomed-out run writes its collapse as a final
                      # converged step, which then reads as the answer
PLATEAU_FRAC = 0.02   # tail dq/ds < 2 % of the initial tangent == plateau

# Clause 3: how far above the step floor the controller must still be over the
# final tenth of the run for its flat tail to count as a mechanism rather than a
# seizure.  100x is chosen from the measured separation, not from taste: the
# three gate legs end their tails at 2500 / 1250 / 800 x DS_MIN, the associated
# control at 1.6 x.  100 sits ~8x below the worst genuine leg and ~60x above the
# seizing one — the widest gap available on this problem.
FREE_ADVANCE_FLOOR_FACTOR = 100.0

# The 1-D elastic stress patch (see the module docstring).  Exact to round-off
# on a tensor grid, where the tributary lumping IS the consistent load vector
# for a rectangular bilinear face; 1e-6 is a distribution-error detector with
# six decades of headroom, not a fitted tolerance.
PATCH_RTOL = 1.0e-6

# Termination modes.  The first two can be a capacity; the rest are SEIZURE.
_CAPACITY_MODES = ("TARGET", "BUDGET")
_SEIZURE_MODES = ("FLOOR", "WALL", "TRUNCATED")

# A hard stop so a pathological run cannot hang a battery.  It is NOT a
# termination criterion for the physics: a leg that ends here is mode `WALL`,
# reported INADMISSIBLE, and FAILS the gate rather than quietly contributing a
# q_max.  Sized ~2x the slowest measured leg (1807 s, h0 = 0.25) so an ordinary
# busy box cannot turn the gate red — if it ever binds, RAISE IT, and never
# lower SUBDIV_BUDGET to make a leg fit inside it.
WALL_BUDGET_S = float(os.environ.get("LADRUNO_R3_WALL_BUDGET", 3600.0))

# --- CONSTRAINT 1 — the per-resolution bands --------------------------------
# Centres are the TIMs/APE measured hex-bbar sequence (note 71).  The half-width
# is +/- 3 %, justified by measurement rather than taste:
#   * the largest observed spread on this problem between two INDEPENDENT
#     drivers of the same leg is 0.8 % (fork `dp_strip.py` 1.0020 vs campaign
#     `r3_hex_control.py` 0.9938 at h0 = 0.5; at h0 = 0.25 they agree to 0.04 %,
#     0.9517 vs 0.9513).  3 % is ~4x that, which is the room a threaded PARDISO
#     factorisation plus an adaptive controller (one step converging differently
#     re-sequences everything after it) needs.
#   * it is still tight enough to catch every failure signature this gate exists
#     for: turning b-bar off on the same mesh (`-formulation std`) reads 0.3999,
#     the budget-truncated fine leg 0.8988 (5.5 % low — outside the h0 = 0.25
#     band, which is the point), the associated leg 1.1171.  NOTE those are
#     PATH-FAILURE readings — where a run stopped — not measured capacities, and
#     must not be quoted as ceilings; that is what the three-clause capacity rule
#     above is for.
#   * the three bands are DISJOINT, so they cannot be satisfied by a sequence
#     that is not monotone.
_BAND_HALFWIDTH = 0.03
_MEASURED = {1.0: 1.0842, 0.5: 0.9938, 0.25: 0.9513}
BANDS = {h: (c * (1.0 - _BAND_HALFWIDTH), c * (1.0 + _BAND_HALFWIDTH))
         for h, c in _MEASURED.items()}
RESOLUTIONS = [1.0, 0.5, 0.25]        # coarse -> fine, the refinement sequence
CONTROL_H0 = 0.5                      # the associated falsification control

_HEX40 = re.compile(r"^[0-9a-f]{40}$")


# --- oracle -----------------------------------------------------------------
def _alpha_from_phi_txc(phi_deg):
    s = math.sin(math.radians(phi_deg))
    return 2.0 * s / (math.sqrt(3.0) * (3.0 - s))


def _phi_ps_from_cone(alpha):
    """Mohr-Coulomb angle matching a Drucker-Prager cone of slope alpha in
    PLANE STRAIN (Chen & Han).  Matching at a triaxial corner and then solving
    a plane-strain problem is the error this whole deck exists to avoid."""
    d = 1.0 - 12.0 * alpha ** 2
    return math.degrees(math.atan(math.sqrt(9.0 * alpha ** 2 / d)))


def _n_q(phi_deg):
    ph = math.radians(phi_deg)
    return math.exp(math.pi * math.tan(ph)) * math.tan(math.pi / 4 + ph / 2) ** 2


# --- the structured graded strip mesh (no gmsh — Zone A stays dep-free) -----
def _graded(lo, hi, n, h0, fine_at_lo):
    if abs(h0 * n - (hi - lo)) < 1e-12:
        return np.linspace(lo, hi, n + 1)
    a, b = 1.0, 4.0
    for _ in range(200):
        r = 0.5 * (a + b)
        tot = h0 * n if abs(r - 1) < 1e-12 else h0 * (r ** n - 1) / (r - 1)
        a, b = (r, b) if tot < hi - lo else (a, r)
    r = 0.5 * (a + b)
    e = np.concatenate([[0.0], np.cumsum(h0 * r ** np.arange(n))])
    e *= (hi - lo) / e[-1]
    return lo + e if fine_at_lo else hi - e[::-1]


def _n_graded(length, h0, r=R_GRADE):
    """Elements needed to span `length` starting at h0 and growing at ~r, so
    that refining h0 does NOT also shrink the domain or steepen the grading —
    otherwise the refinement sequence would confound three changes at once."""
    return max(1, int(np.ceil(np.log(1.0 + length * (r - 1.0) / h0) / np.log(r))))


def _strip_mesh(h0):
    nf = int(round(B_FOOT / h0))
    no, nzo = _n_graded(XLIM - 1.0, h0), _n_graded(-ZBOT - 3.0, h0)
    x = np.unique(np.concatenate([_graded(-XLIM, -1.0, no, h0, False),
                                  np.linspace(-1.0, 1.0, nf + 1),
                                  _graded(1.0, XLIM, no, h0, True)]))
    y = np.array([0.0, THICK])
    z = np.unique(np.concatenate([
        _graded(ZBOT, -3.0, nzo, h0, False),
        np.linspace(-3.0, 0.0, int(round(3.0 / h0)) + 1)]))
    nx, ny, nz = len(x), len(y), len(z)
    idx = np.arange(nx * ny * nz).reshape(nx, ny, nz)
    nodes = np.stack(np.meshgrid(x, y, z, indexing="ij"), axis=-1).reshape(-1, 3)
    hexes = np.array([[idx[i, j, k], idx[i+1, j, k], idx[i+1, j+1, k], idx[i, j+1, k],
                       idx[i, j, k+1], idx[i+1, j, k+1], idx[i+1, j+1, k+1],
                       idx[i, j+1, k+1]]
                      for i in range(nx - 1) for j in range(ny - 1)
                      for k in range(nz - 1)], dtype=np.int32)
    p = nodes[hexes]
    vol = np.einsum("ij,ij->i", p[:, 6] - p[:, 0],
                    np.cross(p[:, 1] - p[:, 0], p[:, 3] - p[:, 0]))
    assert (vol > 0).all(), "inverted hexes in the strip mesh"

    tol = 1e-9
    top = np.abs(nodes[:, 2]) < tol
    sets = dict(top=np.where(top)[0],
                bottom=np.where(np.abs(nodes[:, 2] - ZBOT) < tol)[0],
                xface=np.where(np.abs(np.abs(nodes[:, 0]) - XLIM) < tol)[0],
                footing=np.where(top & (np.abs(nodes[:, 0]) <= 1.0 + tol))[0])
    assert len(sets["footing"]) == 2 * (nf + 1), len(sets["footing"])

    trib = np.zeros(len(nodes))
    for h in hexes:                       # linear quad face -> A/4 per node
        face = h[[4, 5, 6, 7]]
        q = nodes[face]
        if np.any(np.abs(q[:, 2]) > tol):
            continue
        # shoelace on the diagonals; np.cross lost its 2-D form in numpy 2
        d1, d2 = q[2, :2] - q[0, :2], q[3, :2] - q[1, :2]
        trib[face] += 0.5 * abs(d1[0] * d2[1] - d1[1] * d2[0]) / 4.0
    assert abs(trib.sum() - 2 * XLIM * THICK) < 1e-6, trib.sum()
    return nodes, hexes, trib, sets


def _pick_system():
    # The tangent is UNSYMMETRIC on the non-associated leg (rho_bar != rho).
    try:
        ops.system("Pardiso")
        return "Pardiso"
    except Exception:
        ops.system("UmfPack")
        return "UmfPack"


# --- one leg ----------------------------------------------------------------
def _run_leg(h0, assoc, out_dir, wall_budget=None):
    """Push the strip footing to collapse and return the measurement.

    Never raises on a physics outcome: a leg that fails to reach a plateau, or
    that is stopped by a guard, comes back with `verdict` set and the tests
    decide.  Only a broken MODEL (bad mesh, failed surcharge step, violated
    equilibrium identity) raises here.
    """
    wall_budget = WALL_BUDGET_S if wall_budget is None else wall_budget
    nodes, hexes, trib, sets = _strip_mesh(h0)
    alpha = _alpha_from_phi_txc(PHI_TXC)
    rho = math.sqrt(2.0) * alpha
    rho_bar = rho if assoc else 0.0
    g_el = 3.0 * K_EL * (1.0 - 2.0 * NU) / (2.0 * (1.0 + NU))
    phi_ps = _phi_ps_from_cone(alpha)
    q_exact = Q0 * _n_q(phi_ps)

    # The initial 1-D state must be far from yield or the leg is void.
    k0 = NU / (1.0 - NU)
    m0 = (1.0 - k0) / (math.sqrt(3.0) * alpha * (1.0 + 2.0 * k0))
    assert m0 < 0.8, (
        f"initial 1-D state at m0 = {m0:.4f} of yield — the leg would measure "
        "its own initial condition, not a collapse load; raise nu"
    )

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(nodes, start=1):
        ops.node(i, float(x), float(y), float(z))
    # K G sigma_y rho rho_bar Kinf Ko delta1 delta2 H theta -> perfectly plastic
    ops.nDMaterial("DruckerPrager", 1, K_EL, g_el, SY, rho, rho_bar,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    for e, conn in enumerate(hexes, start=1):
        ops.element("LadrunoBrick", e, *[int(c) + 1 for c in conn], 1,
                    "-geom", "linear", "-b", 0.0, 0.0, 0.0,
                    "-formulation", "bbar")
    # ONE fix() per node — a second fix() on the same node adds a DUPLICATE
    # SP_Constraint rather than replacing it.  u_y = 0 on EVERY node is what
    # makes this plane strain, and it is exact for any element type.
    bt, xf = set(sets["bottom"].tolist()), set(sets["xface"].tolist())
    for n in range(len(nodes)):
        if n in bt:
            ops.fix(n + 1, 1, 1, 1)
        else:
            ops.fix(n + 1, 1 if n in xf else 0, 1, 0)

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for n in sets["top"]:
        if trib[n] > 0:
            ops.load(int(n) + 1, 0.0, 0.0, -Q0 * float(trib[n]))

    want = Q0 * float(trib.sum())
    tol = 1.0e-5 * want
    ops.constraints("Transformation")
    ops.numberer("RCM")
    solver = _pick_system()
    ops.test("NormUnbalance", tol, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, f"h0={h0} assoc={assoc}: surcharge step failed"
    ops.reactions()
    rz = sum(ops.nodeReaction(int(n) + 1, 3) for n in sets["bottom"])
    resultant_err = abs(rz / want - 1)
    assert resultant_err < 1e-6, (
        f"initial equilibrium identity violated: {rz} vs {want} kN — the "
        "tributary/body-force convention is wrong, nothing downstream is valid"
    )
    # ...and the check the identity above CANNOT make (see the module docstring):
    # the exact 1-D elastic state.  A distribution error reproduces the total to
    # machine zero and still puts O(1) into sigma_zz.
    k0_el = NU / (1.0 - NU)
    e_zz, e_xx = 0.0, 0.0
    for e in range(1, len(hexes) + 1):
        st = ops.eleResponse(e, "material", 1, "stress")
        if not st or len(st) < 6:
            continue
        e_zz = max(e_zz, abs(st[2] + Q0) / Q0)
        e_xx = max(e_xx, abs(st[0] + k0_el * Q0) / (k0_el * Q0))
    patch_err = max(e_zz, e_xx)

    # ---- the push ----------------------------------------------------------
    foot = [int(n) + 1 for n in sets["footing"]]
    r_corr = Q0 * float(trib[sets["footing"]].sum())   # add back what we applied
    area = B_FOOT * THICK
    uz0 = ops.nodeDisp(foot[0], 3)
    smax = SFRAC * B_FOOT
    # Pseudo-time IS the imposed footing displacement: a Linear series makes the
    # pattern factor equal to lambda, the sp value is 1.0, and starting the clock
    # at uz0 means the first increment carries no jump.  (dp_strip.py uses a Path
    # series over [0,1]; that form does not survive the trip into openseespy.)
    ops.loadConst("-time", uz0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for t in foot:
        ops.sp(t, 3, 1.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    _pick_system()
    ops.analysis("Static")

    # A perfectly plastic collapse is where Newton is worst: the tangent goes
    # singular in the mechanism mode exactly when the answer arrives.  Each
    # increment gets this ladder before the step is shrunk, and any step that
    # needed the RELAXED tolerance is flagged so it cannot pass unnoticed.
    ladder = [("Newton", tol, 25, 0),
              ("NewtonLineSearch", tol, 40, 0),
              ("KrylovNewton", 10.0 * tol, 60, 1)]

    tag = f"h{h0}_{'assoc' if assoc else 'nonassoc'}"
    csv_path = os.path.join(out_dir, f"r3_{tag}.csv")
    fh = open(csv_path, "w", newline="")
    w = csv.writer(fh)
    w.writerow(["s_m", "s_over_B", "q_kPa", "ds_mm", "relaxed", "wall_s"])
    rows, ds, good, nfail, nrelax, nsub = [], DS_BASE, 0, 0, 0, 0
    mode = "TARGET"
    verdict, t0 = "reached the target settlement", time.time()
    while True:
        s_now = uz0 - ops.getTime()
        if s_now >= smax - 1e-12:
            break
        if time.time() - t0 > wall_budget:
            mode = "WALL"
            verdict = f"WALL-CLOCK budget spent at s/B = {s_now / B_FOOT:.5f}"
            break
        ds = min(ds, smax - s_now)
        ops.integrator("LoadControl", -ds)
        ok, relaxed = False, 0
        for algo, tl, it, rl in ladder:
            ops.test("NormUnbalance", tl, it, 0)
            ops.algorithm(algo)
            if ops.analyze(1) == 0:
                ok, relaxed = True, rl
                break
            nfail += 1
        if not ok:
            good = 0
            nsub += 1
            ds *= 0.5
            if nsub > SUBDIV_BUDGET:
                mode = "BUDGET"
                verdict = (f"subdivision budget of {SUBDIV_BUDGET} spent at "
                           f"s/B = {s_now / B_FOOT:.5f}")
                break
            if ds < DS_MIN:
                mode = "FLOOR"
                verdict = (f"step collapsed to the DS_MIN floor at s/B = "
                           f"{s_now / B_FOOT:.5f} (every ladder rung failed at "
                           f"ds = {ds*1000:.6f} mm)")
                break
            continue
        nrelax += relaxed
        good += 1
        if good >= GROW_AFTER and ds < DS_MAX:          # the recovery latch
            ds, good = min(2 * ds, DS_MAX), 0
        ops.reactions()
        q = (-sum(ops.nodeReaction(t, 3) for t in foot) + r_corr) / area
        s = uz0 - ops.getTime()
        rows.append((s, s / B_FOOT, q, ds * 1000, relaxed, time.time() - t0))
        w.writerow([f"{v:.9g}" for v in rows[-1]])
        fh.flush()
    fh.close()
    wall = time.time() - t0

    assert len(rows) > 8, (
        f"h0={h0} assoc={assoc}: only {len(rows)} converged steps — the push "
        f"never started ({verdict})"
    )
    a = np.array([r[:4] for r in rows])
    s, q, dsm = a[:, 0], a[:, 2], a[:, 3]      # dsm is the step in MILLIMETRES
    # D16's truncation rule (see DROP_TRUNCATE).  A truncated curve is a SEIZURE
    # mode: the run bottomed out and wrote its collapse as a converged step.
    peak = np.maximum.accumulate(q)
    bad = np.where(q < (1.0 - DROP_TRUNCATE) * peak)[0]
    if len(bad) and bad[0] > 8:
        mode = "TRUNCATED"
        verdict += f"; curve TRUNCATED ({len(q) - bad[0]} steps cut)"
        s, q, dsm = s[:bad[0]], q[:bad[0]], dsm[:bad[0]]

    msk = s >= 0.9 * s[-1]                     # "the final tenth of the run"
    t_last = float(np.polyfit(s[msk], q[msk], 1)[0]) if msk.sum() > 2 else float("nan")
    n0 = max(4, len(s) // 50)
    t_init = float(np.polyfit(s[:n0], q[:n0], 1)[0])
    plateau = bool(abs(t_last) < PLATEAU_FRAC * abs(t_init))
    qmax = float(q.max())

    # Clause 3 — was the controller still advancing freely when it flattened?
    floor_mm = DS_MIN * 1000.0
    ds_end_mm = float(dsm[-1])
    ds_tail_min_mm = float(dsm[msk].min()) if msk.sum() else ds_end_mm
    headroom = ds_tail_min_mm / floor_mm
    free = bool(headroom >= FREE_ADVANCE_FLOOR_FACTOR)
    budget_used_frac = nsub / float(SUBDIV_BUDGET)
    capacity = bool(plateau and free and mode in _CAPACITY_MODES)

    return dict(h0=h0, assoc=assoc, tag=tag, solver=solver,
                nodes=len(nodes), hexes=len(hexes), dof=3 * len(nodes),
                qmax=qmax, q_exact=q_exact, ratio=qmax / q_exact,
                plateau=plateau, free=free, capacity=capacity, mode=mode,
                t_init=t_init, t_last=t_last,
                tail_pct=100.0 * t_last / t_init, s_end_over_B=float(s[-1]) / B_FOOT,
                s_target_over_B=SFRAC, ds_end_mm=ds_end_mm,
                ds_tail_min_mm=ds_tail_min_mm, ds_floor_mm=floor_mm,
                headroom=headroom, budget_used_frac=budget_used_frac,
                resultant_err=resultant_err, patch_err=patch_err,
                steps=len(rows), nfail=nfail, nrelax=nrelax, nsub=nsub,
                verdict=verdict, wall_s=wall, m0=m0, phi_ps=phi_ps,
                csv=csv_path)


def _why(r):
    """The whole termination story in one string.  Every failure message below
    ends with this, because a gate that fails with only a number reproduces the
    exact confusion the capacity rule was written to prevent."""
    return (
        f"[{r['tag']}] mode={r['mode']} ({r['verdict']}); "
        f"ended s/B={r['s_end_over_B']:.5f} of target {r['s_target_over_B']:.3f}; "
        f"plateau={'yes' if r['plateau'] else 'NO'} (tail {r['tail_pct']:.3f} % of "
        f"initial, needs < {100*PLATEAU_FRAC:.0f} %); "
        f"free-advance={'yes' if r['free'] else 'NO'} (tail-min step "
        f"{r['ds_tail_min_mm']:.6g} mm = {r['headroom']:.1f}x the "
        f"{r['ds_floor_mm']:.6g} mm floor, needs >= "
        f"{FREE_ADVANCE_FLOOR_FACTOR:.0f}x; terminal step {r['ds_end_mm']:.6g} mm); "
        f"subdivisions {r['nsub']}/{SUBDIV_BUDGET} pinned "
        f"({100*r['budget_used_frac']:.0f} % of budget); "
        f"{r['steps']} steps, {r['nfail']} failed attempts, {r['nrelax']} relaxed, "
        f"{r['wall_s']:.0f} s"
    )


# --- the run (once per session; every test below reads the same measurement) -
@pytest.fixture(scope="module")
def campaign(tmp_path_factory):
    out = tmp_path_factory.mktemp("r3_prandtl")
    build = ops.ladrunoBuild() if hasattr(ops, "ladrunoBuild") else "unknown"
    t0 = time.time()
    print(f"\n=== R3 Prandtl-Reissner collapse gate =========================")
    print(f"    engine build (ladrunoBuild)  : {build}")
    print(f"    openseespy module            : {os.path.abspath(ops.__file__)}")
    print(f"    SUBDIV_BUDGET (PINNED)       : {SUBDIV_BUDGET}   "
          f"(D16 default 24 reads 0.8988 still hardening at h0=0.25)")
    print(f"    wall budget per leg          : {WALL_BUDGET_S:.0f} s")
    legs = {}
    for h0 in RESOLUTIONS:
        legs[(h0, False)] = _run_leg(h0, False, str(out))
    legs[(CONTROL_H0, True)] = _run_leg(CONTROL_H0, True, str(out))
    total = time.time() - t0

    hdr = (f"{'leg':>18} {'DOF':>7} {'q_num':>9} {'ratio':>8} {'tail %':>9} "
           f"{'plat':>5} {'mode':>10} {'ds/floor':>9} {'free':>5} {'CAP':>5} "
           f"{'s_end/B':>8} {'wall s':>8}")
    print(f"\n    exact q_u = q0*N_q = {next(iter(legs.values()))['q_exact']:.3f} kPa")
    print(f"{hdr}\n{'-' * len(hdr)}")
    for r in legs.values():
        print(f"{r['tag']:>18} {r['dof']:>7} {r['qmax']:9.2f} {r['ratio']:8.4f} "
              f"{r['tail_pct']:9.3f} {'yes' if r['plateau'] else 'NO':>5} "
              f"{r['mode']:>10} {r['headroom']:9.1f} "
              f"{'yes' if r['free'] else 'NO':>5} "
              f"{'yes' if r['capacity'] else 'NO':>5} {r['s_end_over_B']:8.4f} "
              f"{r['wall_s']:8.1f}")
    print("\n    termination stories (a red gate must be readable):")
    for r in legs.values():
        print(f"      {_why(r)}")
    print("\n    surcharge step, checked twice (the second is the discriminating one):")
    for r in legs.values():
        print(f"      [{r['tag']}] resultant identity err {r['resultant_err']:.3e}"
              f"  |  1-D stress patch err {r['patch_err']:.3e}")
    print(f"\n    TOTAL WALL TIME: {total:.1f} s ({total/60:.1f} min), "
          f"solver = {next(iter(legs.values()))['solver']}")
    print("=" * 63)
    return dict(build=build, legs=legs, total_wall_s=total)


# --- provenance -------------------------------------------------------------
def test_r3_records_the_engine_build_stamp(campaign):
    """Provenance is part of the measurement in this campaign: a collapse load
    with no engine hash attached cannot be re-attributed later.  `ladrunoBuild`
    (#718) is the machine-readable form of the banner's build line."""
    build = campaign["build"]
    if not hasattr(ops, "ladrunoBuild"):
        pytest.skip("openseespy wheel fallback in use — fork build required")
    assert _HEX40.match(build), (
        f"ladrunoBuild() returned {build!r}: this run has no provenance and its "
        "numbers must not be quoted"
    )


# --- CONSTRAINT 1a: per-resolution bands ------------------------------------
@pytest.mark.parametrize("h0", RESOLUTIONS)
def test_r3_ratio_within_its_own_resolution_band(campaign, h0):
    """Each resolution has its OWN band.  Do not replace these with a flat
    percentage: the coarse mesh is 8.4 % above exact and is CORRECT there."""
    r = campaign["legs"][(h0, False)]
    lo, hi = BANDS[h0]
    assert lo <= r["ratio"] <= hi, (
        f"h0 = {h0}: q_num/q_exact = {r['ratio']:.4f} outside its band "
        f"[{lo:.4f}, {hi:.4f}] (centre {_MEASURED[h0]:.4f} +/- "
        f"{100*_BAND_HALFWIDTH:.0f} %). q_max = {r['qmax']:.2f} vs exact "
        f"{r['q_exact']:.2f} kPa; ended at s/B = {r['s_end_over_B']:.4f} "
        f"[{r['verdict']}]"
    )


# --- CONSTRAINT 1b: the real assertion --------------------------------------
def test_r3_sequence_decreases_monotonically_under_refinement(campaign):
    """Convergence FROM ABOVE.  Displacement-based elements over-predict a limit
    load on a coarse mesh — too few DOF to form the true velocity discontinuity,
    so the cheapest available mechanism is stiffer and stronger — and refinement
    relaxes it.  A sequence that turns over is the signature of a guard binding
    before the mechanics (see SUBDIV_BUDGET), not of convergence."""
    seq = [campaign["legs"][(h, False)]["ratio"] for h in RESOLUTIONS]
    pretty = " -> ".join(f"{h}: {v:.4f}" for h, v in zip(RESOLUTIONS, seq))
    for i in range(len(seq) - 1):
        assert seq[i + 1] < seq[i], (
            f"refinement sequence is not monotone decreasing ({pretty}); "
            f"h0 = {RESOLUTIONS[i+1]} did not fall below h0 = {RESOLUTIONS[i]}"
        )
    assert seq[0] > 1.0, (
        f"the coarse mesh should over-predict the limit load ({pretty}); a "
        "coarse rung already below exact means the sequence is not converging "
        "from above and the reading is not the classical one"
    )


def test_r3_every_leg_of_the_sequence_is_on_a_plateau(campaign):
    """Clause 1.  A q_max without a plateau is NOT a collapse load: a leg still
    hardening reports where its solver stopped, not where the soil failed."""
    bad = [r for r in (campaign["legs"][(h, False)] for h in RESOLUTIONS)
           if not r["plateau"]]
    assert not bad, "\n".join(
        f"h0 = {r['h0']}: STILL HARDENING at q = {r['qmax']:.2f} kPa — tail "
        f"dq/ds = {r['t_last']:.1f} kPa/m = {r['tail_pct']:.3f} % of initial. "
        + _why(r) for r in bad
    )


def test_r3_every_leg_was_still_advancing_freely_when_it_flattened(campaign):
    """Clause 3 — the companion to the plateau test, and the one that makes it
    mean anything.  A curve can flatten because a mechanism formed OR because
    the run SEIZED: a step floor reached just as the curve rolls over presents
    as a flat tail too.  So a flat tangent counts only if the controller still
    had room when it flattened."""
    bad = [r for r in (campaign["legs"][(h, False)] for h in RESOLUTIONS)
           if not r["free"]]
    assert not bad, "\n".join(
        f"h0 = {r['h0']}: the tail is flat but the controller was AT THE FLOOR — "
        f"this is SEIZURE, not a mechanism, and q = {r['qmax']:.2f} kPa is not a "
        f"capacity. " + _why(r) for r in bad
    )


def test_r3_termination_mode_is_admissible_and_named(campaign):
    """Clause 2, plus the naming rule.  `TARGET` is a capacity outright;
    `BUDGET` is a capacity WITH A NAMED CONTROLLER ALLOWANCE (the load was
    already flat and the guards simply declined to spend another hour confirming
    it — the reference records exactly this for its own control leg).  `FLOOR`,
    `WALL` and `TRUNCATED` are seizure and are never a capacity, whatever the
    tail tangent says.

    Note what this test does NOT do: it does not silently accept a budget halt.
    Every leg's mode, step headroom and budget usage are printed by the fixture
    and repeated in this message, so the allowance is on the record."""
    legs = [campaign["legs"][(h, False)] for h in RESOLUTIONS]
    seized = [r for r in legs if r["mode"] in _SEIZURE_MODES]
    assert not seized, "\n".join(
        f"h0 = {r['h0']}: terminated in SEIZURE mode {r['mode']} — its q_max is "
        f"INADMISSIBLE as a capacity"
        + ("; raise LADRUNO_R3_WALL_BUDGET and re-run, and do NOT lower "
           "SUBDIV_BUDGET to make it fit" if r["mode"] == "WALL" else "")
        + ". " + _why(r) for r in seized
    )
    allowance = [r for r in legs if r["mode"] == "BUDGET"]
    if allowance:                       # named, not silent
        print("\n    CONTROLLER ALLOWANCE named for "
              + ", ".join(f"h0={r['h0']} (stopped at s/B={r['s_end_over_B']:.4f} "
                          f"of {r['s_target_over_B']:.3f}, tail "
                          f"{r['tail_pct']:.3f} %, step {r['headroom']:.0f}x floor)"
                          for r in allowance))


def test_r3_every_leg_of_the_sequence_is_a_capacity(campaign):
    """The three clauses together — this is the assertion the ratios above are
    only meaningful under."""
    bad = [r for r in (campaign["legs"][(h, False)] for h in RESOLUTIONS)
           if not r["capacity"]]
    assert not bad, "\n".join(
        f"h0 = {r['h0']}: q = {r['qmax']:.2f} kPa is NOT a capacity "
        f"(plateau={r['plateau']}, free-advance={r['free']}, mode={r['mode']}). "
        + _why(r) for r in bad
    )


# --- the surcharge step, checked twice ---------------------------------------
def test_r3_surcharge_load_distribution_not_only_its_resultant(campaign):
    """The resultant identity `sum(R_z) == q0 * sum(A_trib)` is STRUCTURALLY
    BLIND to load-distribution error — a tributary lumping can reproduce its
    applied total to +0.0000000 % while putting O(1) error into sigma_zz.  So it
    is paired with the exact 1-D elastic state this domain has in closed form:
    sigma_zz = -q0 and sigma_xx = sigma_yy = -K0 q0 at every Gauss point.  The
    two discriminate totally: only the second can see a distribution error."""
    bad = [r for r in campaign["legs"].values() if r["patch_err"] > PATCH_RTOL]
    assert not bad, "\n".join(
        f"[{r['tag']}] 1-D elastic stress patch FAILED: max relative error "
        f"{r['patch_err']:.3e} > {PATCH_RTOL:.0e}, while the resultant identity "
        f"passed at {r['resultant_err']:.3e}. That combination IS the signature "
        "of a load-DISTRIBUTION error (right total, wrong distribution) — the "
        "surcharge is not the one the oracle assumes and no q_max below is valid."
        for r in bad
    )


# --- the falsification control ----------------------------------------------
def test_r3_associated_control_is_not_a_capacity(campaign):
    """The control that makes the agreement above a measurement rather than a
    coincidence.  Associated flow (psi = phi) has to push the surrounding soil
    aside; a bounded mesh resists it and the leg hardens straight past its own
    exact answer.  If the associated leg ALSO produced a capacity, this deck
    would be plateauing on something other than the collapse mechanism and the
    non-associated agreement would carry no information.

    Note this is asserted on CAPACITY, not on `plateau` alone — precisely
    because the associated leg's characteristic ending is a step floor, and a
    seized run can present a flat tail.  Measured, it fails both clauses at
    once: tail ~39 % of initial AND a terminal step ~1.6x the floor."""
    r = campaign["legs"][(CONTROL_H0, True)]
    assert r["mode"] != "WALL", (
        f"associated control INCONCLUSIVE: stopped by the wall budget, so 'not a "
        f"capacity' would be an artefact of the clock, not of the flow rule. "
        f"Raise LADRUNO_R3_WALL_BUDGET and re-run. " + _why(r)
    )
    assert not r["capacity"], (
        f"associated control PRODUCED A CAPACITY (ratio {r['ratio']:.4f}) — the "
        "reference records associated legs hardening past their own exact answer "
        "without ever reaching one. Both flow rules yielding a capacity means the "
        "non-associated agreement is a coincidence, not a measurement of the "
        "Prandtl mechanism. " + _why(r)
    )
