"""ADR-85 T1b -- the 2D NTS lane goes LIVE (gate G-T1b b / e / f / g).

WHAT THIS FILE IS
    T0 opened the 2D DECLARATION path and refused every 2D pair by name.  T1b
    builds the machinery behind that refusal -- the FE SEGMENT 2D construction
    path, the handler's 2D branch, the new interface-level orientation vote, the
    2-component `-outward`, the ndm-generalized auto-kn -- and this file is the
    acceptance half of its gate.  The corner primitive has its own file
    (test_adr85_contact2d_t1b_corner.py) and the slice twin / plane-stress
    headline has a third (test_adr85_contact2d_t1b_slicetwin.py).

COVERAGE MAP -> ADR-85 G-T1b
    (b) 2D compression patch across a NON-MATCHING interface transmits the load
            test_2d_patch_compression             (+ an embedded reference solve)
            test_2d_patch_initial_tangent         (the same, on the `-initial` arm)
    (e) `-geomtan` gate                       test_2d_geomtan_noop
    (f) auto-kn sizes a usable penalty        test_2d_autokn
    (g) flush-interface degenerate orientation vote -> named abort
            test_2d_orientation_vote_flush        (+ the -outward control leg)
            test_2d_orientation_vote_resolves     (the vote's POSITIVE control)
            test_2d_outward_two_components        (the 2-component parse + flip)
    regression (T1b must not WIDEN T0's guards)
            test_2d_ndf3_still_refused, test_2d_mortar_still_refused

    G-T1b(a) (3D battery N-unchanged + contact_dump bit-identical), (c) (the
    slice twin) and (d) (plane stress) are NOT here: (a) is the whole shipped
    battery plus tests/_testbed/contact_dump.py, run by the PR, and (c)/(d) live
    in the slicetwin file.

SIBLING FILES IN THE T1b BATTERY
    test_adr85_contact2d_t1b_corner.py     the 2D corner primitive
    test_adr85_contact2d_t1b_slicetwin.py  G-T1b(c) + (d)
    test_adr85_contact2d_t1b_guards.py     the merge panel's guard matrix --
        the degenerate-vote magnitude gate, master-chain integrity, the SPLIT
        vote (a U-channel whose normals span 180 degrees), and the mirrored
        cross-dimension pair.  All named FATALs, all child-process.
    test_adr85_contact2d_t1b_explicit.py   `-visc` and `-soft` on the SEGMENT
        lane (the banner claims both; the dashpot is a stride-3 PORT).

THE `-geomtan` CONTRACT IS A DECLARED NO-OP -- NOT AN IDENTITY
    An earlier draft of this file said the 2D geometric term is "exactly zero
    for a straight segment, the d(n)/du and d(xi)/du contributions cancel
    identically".  That was WRONG, and the ADR How/5 correction says why: what
    the T1a oracle FD-gated is that the 2D B-operator is the exact FIRST
    VARIATION of the gap.  The consistent tangent needs the SECOND variation,
    and its curvature term is NOT zero even on a straight segment -- a master
    segment's unit normal rotates as its nodes move, so d(n)/du has content
    regardless of how straight the segment is.

    T1b therefore does not implement it.  `-geomtan` on a 2D pair is ACCEPTED
    and does nothing, by DECISION (the term is deferred), and the engine emits
    a warnOnce saying so.  The gate is that pair of facts:

      * BYTE-IDENTITY of the solved state with and without the flag -- the
        decision is "no-op", so anything else means a term was formed;
      * the WARNING is emitted -- a silent no-op flag is a user writing
        `-geomtan`, believing they bought the consistent tangent, and getting
        the main term.  That is the ADR-76 `-initial` failure mode exactly
        (a flag that silently does not do what its name says).

    Byte-identity, not the 3D iteration-count idiom: an iteration-count test on
    a flag that is currently a no-op asserts nothing, and would silently pass a
    `-geomtan` that added a small WRONG term.

WHY A CHILD PROCESS FOR ONE TEST ONLY
    Only test_2d_orientation_vote_flush asserts a FATAL.  Contact refusals are
    not ordinary Python exceptions (serial: -1 -> OpenSeesError; MPI:
    MPI_Abort; handler-time: a failed analyze after domainChanged), so they run
    in a child interpreter whose exit code AND message the parent asserts.  The
    full rationale, including both Windows traps (the PYTHONPATH 32 KB
    environment-block trap and the stdin=DEVNULL Popen trap), is in the module
    docstring of tests/test_adr85_contact2d_t0_refusals.py, from which
    ENGINE_DIR and the two regression rows are imported so the two files cannot
    drift.  Everything else here is in-process with ops.wipe().

TOLERANCE POLICY (ADR-85 How/1)
    Every bound below is RELATIVE to a deck quantity and carries the reason it
    has the value it has.  Two classes, kept apart on purpose:

      * EQUILIBRIUM IDENTITIES (summed contact force == applied load; the
        closed-form penalty equilibrium of a single pair) are exact statements
        about a converged solve, so they are gated at 1e-6 relative against a
        NormDispIncr 1e-12 solve;
      * DISCRETIZATION statements (interface pressure uniformity across a
        non-matching NTS interface) are single-digit-percent statements -- NTS
        does not pass the patch test, which is exactly why the mortar lane (T3)
        exists and is gated at machine precision instead.  Those bounds print
        their measured value so the first real run can tighten them.
"""
import math

import numpy as np
import pytest

from _testbed import ops
from test_adr85_contact2d_t0_refusals import ENGINE_DIR, _assert_refused

pytestmark = [pytest.mark.zone_a]

# ---- deck constants (every tolerance below is expressed against one of these)
E = 1.0e4          # Young's modulus of both blocks
NU = 0.0           # nu = 0 => the PlaneStrain uniaxial modulus is exactly E,
                   # so the analytic column compression is P*h/(W*E) with no
                   # (1-nu)/((1+nu)(1-2nu)) factor to argue about
THICK = 1.0        # 2D plane elements assemble ACTUAL forces (h is inside the
                   # stiffness); kn is force/length on a NODAL gap and carries
                   # no thickness at all -- ADR-85 How/7
KN = 1.0e6         # NTS penalty
SEED = 1.0e-8      # seeded penetration: the pair is active from step 1 (a slave
                   # held only by a penalty is singular while separated)
W = 2.0            # interface width == the deck's length scale
H = 1.0            # block height

# The 2D element is upstream FourNodeQuad (`element quad`), the same one the T0
# rigid-plane file uses (tests/test_adr85_contact2d_t0_rigidplane.py:84) and the
# same one the slice twin pairs against `stdBrick`.  Chosen over the fork's
# LadrunoQuad deliberately: the twin's element-formulation floor is smallest
# between two plain isoparametric siblings, and the contact lane is what is
# under test here, not the element.


# --------------------------------------------------------------------------
# builders
# --------------------------------------------------------------------------
def _tributary(n, width):
    """Consistent nodal tributary lengths of an edge meshed into n equal
    elements over `width`: end nodes own half an element, interior nodes a
    whole one.  Sums to `width` by construction, so a uniform traction maps to
    P*w_i/width and the mean nodal pressure f_i/w_i is the applied one."""
    le = width / n
    return [le / 2.0 if j in (0, n) else le for j in range(n + 1)]


def _patch_model(ntop=5, nbot=2, geomtan=False, kn=KN):
    """Two stacked PlaneStrain blocks with NON-MATCHING interface meshes.

    bottom: x in [0,W], y in [0,H],  nbot elements, base fully fixed
    top:    x in [0,W], y in [H,2H], ntop elements, riding on the bottom block

    x is fixed on every node.  That is not decoration: with nu = 0 the exact
    solution of the laterally confined column is uniform sigma_yy and u_x == 0,
    which the bilinear quad reproduces EXACTLY, so every deviation measured
    below is attributable to the interface and nothing else.  The interface
    still engages shear, because NTS hands the master a load pattern that is
    NOT the consistent one for the master mesh (that is the discretization
    error being measured).

    `-outward 0 1` is given explicitly and is the T1b 2-component parser branch:
    the slave row is seeded SEED below the master line so the pair is active
    from step 1, which leaves the interface-level centroid vote degenerate
    (ADR-85 How/2) -- the flush case that must be told, not guessed.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)

    base, mid = [], []
    for i in range(nbot + 1):
        x = W * i / nbot
        ops.node(100 + i, x, 0.0)
        ops.fix(100 + i, 1, 1)
        ops.node(110 + i, x, H)
        ops.fix(110 + i, 1, 0)             # x fixed, y free
        base.append(100 + i)
        mid.append(110 + i)
    for i in range(nbot):
        ops.element("quad", 300 + i, base[i], base[i + 1], mid[i + 1], mid[i],
                    THICK, "PlaneStrain", 1)

    slave, top = [], []
    for j in range(ntop + 1):
        x = W * j / ntop
        ops.node(200 + j, x, H - SEED)
        ops.fix(200 + j, 1, 0)
        ops.node(210 + j, x, 2.0 * H)
        ops.fix(210 + j, 1, 0)
        slave.append(200 + j)
        top.append(210 + j)
    for j in range(ntop):
        ops.element("quad", 400 + j, slave[j], slave[j + 1], top[j + 1], top[j],
                    THICK, "PlaneStrain", 1)

    mseg = []
    for i in range(nbot):                  # nps = 2 => the node list is read as
        mseg += [mid[i], mid[i + 1]]       # consecutive 2-node LINE segments
    ops.contactSurface(10, "-master", 2, *mseg)
    ops.contactSurface(20, "-slave", *slave)
    args = [1, 10, 20, kn, 0.0, 0.0, "-outward", 0.0, 1.0]
    if geomtan:
        args.append("-geomtan")
    ops.contact(*args)
    return dict(base=base, mid=mid, slave=slave, top=top, ntop=ntop, nbot=nbot)


def _load_uniform(m, P):
    """Uniform traction on the top edge, as CONSISTENT nodal loads (so the top
    block's own discretization contributes no error to the pressure profile)."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, w in zip(m["top"], _tributary(m["ntop"], W)):
        ops.load(t, 0.0, -P * w / W)


def _static(nsteps=1, tol=1.0e-12, iters=40, what="2D NTS patch", algo=("Newton",)):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", tol, iters, 0)
    ops.algorithm(*algo)
    ops.analysis("Static")
    ops.integrator("LoadControl", 1.0 / nsteps)
    for k in range(nsteps):
        assert ops.analyze(1) == 0, f"{what}: step {k} did not converge"


def _wmean(tags, n):
    """Tributary-weighted mean of u_y over an edge row of n+1 nodes.  A plain
    mean would be biased by the half-tributary end nodes whenever u_y varies
    along the row, which is precisely the situation being measured."""
    w = _tributary(n, W)
    return sum(wi * ops.nodeDisp(t, 2) for t, wi in zip(tags, w)) / W


# --------------------------------------------------------------------------
# the REFERENCE SOLVE -- an independent closed-form answer for the patch deck
#
# WHY THIS IS EMBEDDED RATHER THAN QUOTED.  The patch test is the only case in
# the battery that exercises multi-segment master + ownership + the terminal
# acceptance window TOGETHER, and its only quantitative bound used to be a
# +-8% uniformity band -- so a transfer bug anywhere inside that band passed.
# Quoting numbers from an uncommitted script in a docstring is not a fix: the
# script is not run, not reviewed and not kept in sync.  This is the script.
#
# WHY THE REDUCTION IS EXACT, NOT AN APPROXIMATION.  The patch deck fixes u_x on
# EVERY node (a BC, not a symmetry argument) and uses nu = 0, so for a bilinear
# quad the strain from the surviving y-DOFs is eps_xx = 0, eps_yy = u,y,
# gamma_xy = u,x and the energy collapses to
#
#     0.5 * t * integral[ E*(u,y)^2 + G*(u,x)^2 ],   G = mu = E/2  at nu = 0
#
# integrated on the same 2x2 Gauss rule FourNodeQuad uses -- exact for these
# quadratic integrands on a rectangle.  The NTS penalty is B = [1 | -N0 | -N1]
# per slave node with the shape functions of the segment the slave projects
# into, and the seeded gap enters as a constant on the RHS.
#
# WHAT IS *NOT* MODELLED, and why it does not matter at 1e-3: the engine
# recomputes xi and the normal on the DEFORMED master every iteration, while
# this solve uses the reference-configuration ones.  The master tilts by a slope
# of ~2e-4 here and the gaps are ~5e-5, so the projection moves by ~1e-8 in x
# and the normal tilts by ~2e-4 rad -- relative effects of ~1e-8 and ~2e-8 on
# the forces.  The assertion runs at 1e-3, four orders above that.
# --------------------------------------------------------------------------
_G1 = 1.0 / math.sqrt(3.0)
_GP = ((-_G1, -_G1), (_G1, -_G1), (_G1, _G1), (-_G1, _G1))


def _ke_yonly(a, b, t):
    """y-DOF stiffness of an a x b rectangle of thickness t at nu = 0, nodes in
    the deck's CCW order (bottom-left, bottom-right, top-right, top-left)."""
    K = np.zeros((4, 4))
    for xi, eta in _GP:
        dxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
        det = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
        dNdx, dNdy = dxi * (2.0 / a), det * (2.0 / b)
        K += (E * np.outer(dNdy, dNdy) + 0.5 * E * np.outer(dNdx, dNdx)) \
            * (a * b / 4.0) * t
    return K


def _reference_interface_forces(nbot, ntop, load_total, kn=KN, seed=SEED):
    """Per-slave-node normal contact force of the patch deck, solved directly.

    Returns a list of ntop+1 forces (positive = compression), in the same order
    as the deck's slave row.
    """
    nm, ns = nbot + 1, ntop + 1
    n = nm + 2 * ns
    K = np.zeros((n, n))
    F = np.zeros(n)

    def iM(i):
        return i

    def iS(j):
        return nm + j

    def iT(j):
        return nm + ns + j

    kb = _ke_yonly(W / nbot, H, THICK)
    for i in range(nbot):                       # base row is FIXED -> mid dofs only
        d = (-1, -1, iM(i + 1), iM(i))
        for p in range(4):
            if d[p] < 0:
                continue
            for q in range(4):
                if d[q] < 0:
                    continue
                K[d[p], d[q]] += kb[p, q]
    kt = _ke_yonly(W / ntop, H, THICK)
    for j in range(ntop):
        d = (iS(j), iS(j + 1), iT(j + 1), iT(j))
        for p in range(4):
            for q in range(4):
                K[d[p], d[q]] += kt[p, q]

    seg_len = W / nbot
    rows = []
    for j in range(ns):
        x = W * j / ntop
        seg = min(int(x / seg_len), nbot - 1)   # the last slave lands at xi = 1
        xi = (x - seg * seg_len) / seg_len
        b = np.zeros(n)
        b[iS(j)] += 1.0
        b[iM(seg)] -= 1.0 - xi
        b[iM(seg + 1)] -= xi
        rows.append(b)
        K += kn * np.outer(b, b)
        F += kn * seed * b
    for j, w in enumerate(_tributary(ntop, W)):
        F[iT(j)] -= load_total * w / W

    u = np.linalg.solve(K, F)
    return [-kn * float(b @ u - seed) for b in rows]


# --------------------------------------------------------------------------
# 1. the compression patch across a NON-MATCHING interface -- G-T1b(b)
# --------------------------------------------------------------------------
def test_2d_patch_compression():
    """A 2-element block under a 5-element block, uniform pressure.

    FOUR claims at FOUR different tolerances, because they are four different
    kinds of statement.  The two loose ones are bounded by a hand-solved
    reference (see PREDICTED VALUES below) rather than by a guess:

    1. TOTAL FORCE (equilibrium identity, 1e-6).  The summed
       `ladrunoContactForce` over the slave row equals the applied P.  The top
       block's only vertical forces are the applied loads and the contact, so
       this is exact at convergence -- and it is the assertion that a lane which
       "assembles but transmits nothing" cannot fake.  A zero here with a
       converged solve IS the ADR-78 P0 signature.

    2. BLOCK COMPRESSION (identity, 1e-6).  Each block's TRIBUTARY-WEIGHTED mean
       compression is the uniaxial sigma*H/E with sigma = P/(W*THICK) -- nu = 0
       makes the PlaneStrain uniaxial modulus exactly E.  This is tight, not
       loose, and the weighting is why: the interface's transfer error is a
       SELF-EQUILIBRATED nodal pattern, and the tributary-weighted mean is
       work-conjugate to a uniform traction, so the perturbation drops out
       exactly (the reference solve gives 1e-14 for every mesh pair tried).  A
       plain arithmetic mean would NOT have this property and would need a
       percent-scale bound.  Without this claim, 1 and 3 would still pass on an
       interface that transmitted force into thin air.

    3. INTERFACE PRESSURE UNIFORMITY (discretization, 8%).  f_i / w_i must be
       the applied traction P/W at every slave node, w_i being the slave node's
       tributary length.  This is where NTS is ALLOWED to be wrong: it transfers
       the slave's nodal forces through the MASTER shape functions, which is not
       the consistent load for the master mesh -- for this 2:5 pair the master
       row receives (0.26, 0.48, 0.26)*P where consistency wants
       (0.25, 0.50, 0.25)*P -- so the master settles non-uniformly and feeds a
       non-uniform reaction back.  T3's mortar lane is what makes this machine
       precision; asserting machine precision HERE would be asserting a bug.

    4. INTERFACE SETTLEMENT UNIFORMITY (discretization, 5%).  The same error
       seen on the smoother quantity: the master row's u_y should be constant
       for a uniform traction.  Kept alongside 3 because forces amplify the
       error and displacements smooth it, so together they bracket it.

    5. THE PER-STATION REFERENCE (1e-3).  Claims 3 and 4 are BANDS, and a band
       is exactly what a transfer bug hides in: this deck is the only case in
       the battery that exercises multi-segment master + ownership + the
       terminal acceptance window together, and anything inside +-8% used to
       pass it.  So every one of the six nodal forces is also compared against
       `_reference_interface_forces`, an independent solve of the same deck
       (see its comment block for why the reduction is exact and what the 1e-3
       is buying).  This is the assertion that actually pins the load transfer;
       3 and 4 are kept because they say, in engineering terms, WHAT a
       disagreement would mean.

    PREDICTED VALUES (the embedded reference, evaluated here):

        total force        200.000000 (rel 1e-15)
        nodal forces       19.0536, 39.7388, 41.2076, 41.2076, 39.7388, 19.0536
        nodal pressures    95.27, 99.35, 103.02, 103.02, 99.35, 95.27
        pressure deviation 4.73%       <- bound 8%
        settlement dev     2.28%       <- bound 5%
        both compressions  exact to 1e-14

    The 2:5 mesh pair was CHOSEN on that reference rather than assumed: the
    obvious 2:3 pair puts the pressure deviation at 9.50% and the settlement at
    6.89%, which leaves no headroom under a single-digit bound and would make
    this test a tolerance-tuning exercise on its first real run.  2:5 is just as
    non-matching (no slave node lands on a master node, and none lands on the
    shared master vertex at x = 1) with half the transfer error.
    """
    P = 2.0e2
    m = _patch_model(ntop=5, nbot=2)
    _load_uniform(m, P)
    _static(nsteps=2)

    f = [ops.ladrunoContactForce(t) for t in m["slave"]]
    tot = sum(f)
    assert tot > 0.0, (
        f"the 2D NTS interface transmitted NOTHING (per-node forces {f}) on a "
        "converged solve -- the ADR-78 P0 silent-drop signature")
    assert abs(tot - P) / P < 1.0e-6, (
        f"summed contact force {tot:.8e} != applied {P:.8e}: the top block's "
        f"vertical equilibrium is violated (per-node {f})")

    # 2. both blocks carry the stress (weighted-mean identity)
    comp_ref = (P / (W * THICK)) * H / E
    comp_bot = -_wmean(m["mid"], m["nbot"])                    # base is fixed
    comp_top = _wmean(m["slave"], m["ntop"]) - _wmean(m["top"], m["ntop"])
    for name, comp in (("bottom", comp_bot), ("top", comp_top)):
        assert abs(comp - comp_ref) / comp_ref < 1.0e-6, (
            f"{name} block tributary-weighted compression {comp:.8e} != "
            f"uniaxial sigma*H/E {comp_ref:.8e} -- the transmitted force is not "
            "reaching the material with the right magnitude")

    # 3. pressure uniformity, RELATIVE to the applied traction
    p_ref = P / W
    w = _tributary(m["ntop"], W)
    p = [fi / wi for fi, wi in zip(f, w)]
    dev_p = max(abs(pi / p_ref - 1.0) for pi in p)
    # 4. settlement uniformity on the master row
    u_mid = [ops.nodeDisp(t, 2) for t in m["mid"]]
    dev_u = max(abs(ui - (-comp_bot)) for ui in u_mid) / comp_bot
    print(f"\n[2D NTS patch 2:5] pressures {[float('%.6g' % pi) for pi in p]} "
          f"vs applied {p_ref:.6g}; pressure dev {100.0 * dev_p:.3f}% "
          f"(predicted 4.73%), settlement dev {100.0 * dev_u:.3f}% "
          f"(predicted 2.28%)")
    assert dev_p < 0.08, (
        f"interface pressure non-uniformity {100.0 * dev_p:.2f}% exceeds the "
        f"NTS band (pressures {p}, applied {p_ref}); the reference solve of this "
        "deck says 4.73%, so this is a transfer bug, not discretization")
    assert dev_u < 0.05, (
        f"master-row settlement non-uniformity {100.0 * dev_u:.2f}% (u_y "
        f"{u_mid}, mean {-comp_bot:.6e}); the reference solve says 2.28%")

    # 5. the per-station reference -- the assertion that actually pins the
    #    transfer.  Bands 3-4 say what a disagreement MEANS; this says whether
    #    there is one.
    ref = _reference_interface_forces(m["nbot"], m["ntop"], P)
    print(f"[2D NTS patch 2:5] measured {[float('%.6g' % v) for v in f]} vs "
          f"reference {[float('%.6g' % v) for v in ref]}")
    for j, (got, want) in enumerate(zip(f, ref)):
        assert abs(got - want) / abs(want) < 1.0e-3, (
            f"slave station {j} (x = {W * j / m['ntop']:.4f}): contact force "
            f"{got:.8e} vs the independent reference {want:.8e} "
            f"(rel {abs(got - want) / abs(want):.3e}).  The uniformity bands "
            "above are wide enough to hide a transfer bug; this is not")


# --------------------------------------------------------------------------
# 2. the orientation vote -- G-T1b(g)
# --------------------------------------------------------------------------
# The child is the T0 shape; see this module's docstring for why, and the T0
# module docstring for the two Windows traps encoded below.  Note the template
# is rendered by the percent operator against a one-key mapping, so NOTHING
# inside it -- code or comment -- may contain a bare percent conversion.
CHILD_VOTE = r'''
import os, sys

_D = %(ENGINE_DIR)r
if os.path.isdir(_D):
    os.environ["PATH"] = _D + os.pathsep + os.environ.get("PATH", "")
    _add = getattr(os, "add_dll_directory", None)      # WINDOWS-ONLY: probe, do
    if _add is not None:                               # not try/except OSError
        try:
            _add(_D)
        except OSError:
            pass
    sys.path.insert(0, _D)

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

CASE = sys.argv[1]
KN, KS, P = 1.0e6, 1.0e3, 1.0e3

# A FLUSH, ZERO-GAP interface: the slave node lies EXACTLY on the master line.
# The master surface centroid is (0.5, 0) and so is the slave surface centroid,
# so the interface-level reference vote (ADR-85 How/2) is degenerate -- the
# datum has no length and there is no honest way to pick a side.  This is the
# masonry-joint / footing-on-soil geometry the 2D lane exists for, so it must
# ABORT BY NAME rather than silently drop the pair the way the shipped 3D NTS
# per-pair datum does.
#
# The soft y-spring is what makes the CONTROL leg solvable: at gap == 0 the
# penalty contributes no stiffness on the first iteration, so a slave held only
# by the contact would present a singular tangent at step 1.  It is present in
# BOTH legs so the two differ in exactly one token, the -outward flag.
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.node(101, 0.0, 0.0)
ops.node(102, 1.0, 0.0)
ops.fix(101, 1, 1)
ops.fix(102, 1, 1)
ops.node(1, 0.5, 0.0)
ops.fix(1, 1, 0)
ops.node(2, 0.5, 0.0)
ops.fix(2, 1, 1)
ops.uniaxialMaterial("Elastic", 1, KS)
ops.element("zeroLength", 1, 2, 1, "-mat", 1, "-dir", 2)
ops.contactSurface(10, "-master", 2, 101, 102)
ops.contactSurface(20, "-slave", 1)

if CASE == "vote_flush_no_outward":
    ops.contact(1, 10, 20, KN, 0.0, 0.0)
elif CASE == "vote_flush_with_outward":
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 1.0)
else:
    print("UNKNOWN CASE " + CASE)
    sys.exit(4)

ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(1, 0.0, -P)
ops.constraints("LadrunoContact")
ops.numberer("Plain")
ops.system("FullGeneral")
ops.test("NormDispIncr", 1.0e-12, 40, 0)
ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0)
ops.analysis("Static")
rc = ops.analyze(1)
if rc != 0:
    print("VOTE_ANALYZE_REFUSED " + CASE + " rc=" + str(rc))
    sys.exit(3)
print("VOTE_RAN " + CASE + " u=" + repr(ops.nodeDisp(1, 2))
      + " f=" + repr(ops.ladrunoContactForce(1)))
sys.exit(0)
'''


def _run_vote_child(case):
    import os
    import subprocess
    import sys

    script = CHILD_VOTE % {"ENGINE_DIR": ENGINE_DIR}
    env = dict(os.environ)
    env["LADRUNO_OPENSEES_QUIET"] = "1"     # the banner's UTF-8 glyphs break a
                                            # text-mode capture on cp1252
    proc = subprocess.run(
        [sys.executable, "-c", script, case],
        stdin=subprocess.DEVNULL,           # load-bearing on Windows (T0 docstring)
        capture_output=True, text=True, timeout=300,
        encoding="utf-8", errors="replace",
        env=env,
    )
    return proc, (proc.stdout or "") + "\n" + (proc.stderr or "")


def _grab(out, marker, key):
    for line in out.splitlines():
        if not line.startswith(marker):
            continue
        for tok in line.split():
            if tok.startswith(key + "="):
                return float(tok[len(key) + 1:])
    raise AssertionError(f"no {key!r} on a {marker!r} line:\n{out}")


def test_2d_orientation_vote_flush():
    """G-T1b(g): a flush interface with no `-outward` -> NAMED ABORT; the same
    deck WITH `-outward 0 1` -> runs.

    This is NEW behavior relative to the shipped 3D default, and the ADR says so
    (How/2): the 3D per-pair datum `slave_ref - segment_ref_centroid` degenerates
    SILENTLY for a slave lying on its master segment -- normalOriented's
    fail-safe returns false and the pair is quietly inactive.  A quietly
    inactive contact on a masonry joint is the failure this lane cannot ship
    with, so 2D refuses instead of guessing.

    Both legs run in the SAME child template so they cannot drift apart; only
    the abort leg needs a child, and giving the control leg one costs nothing
    and removes the "the two decks were not really identical" objection.

    The control leg is checked against its closed form, not just "exit 0":
    starting from a zero gap, kn*|u| + ks*|u| = P, so u = -P/(KN+KS) and the
    reported contact force is KN*|u|.  A build that accepted `-outward` and then
    oriented the normal the WRONG way would exit 0 with u = -P/KS.
    """
    proc, out = _run_vote_child("vote_flush_no_outward")
    assert proc.returncode != 0, (
        "a flush 2D interface with a degenerate orientation vote and no "
        f"-outward must ABORT, not run (exit 0).\n{out}")
    assert "VOTE_RAN" not in out, (
        f"the flush deck solved instead of refusing.\n{out}")
    for needle in ("ADR-85", "outward"):
        assert needle in out, (
            f"the degenerate-vote abort fired but does not name itself -- "
            f"missing {needle!r}.  ADR-85 How/2 requires the message to name "
            f"the ADR and to point at the -outward override.\n{out}")

    proc, out = _run_vote_child("vote_flush_with_outward")
    assert proc.returncode == 0, (
        f"`-outward 0.0 1.0` must resolve the flush interface (exit "
        f"{proc.returncode}).  Note the 2-COMPONENT form is itself T1b work: "
        f"T0 left the parser reading exactly 3 doubles.\n{out}")
    assert "VOTE_RAN" in out, f"the control leg did not finish.\n{out}"
    u = _grab(out, "VOTE_RAN", "u")
    f = _grab(out, "VOTE_RAN", "f")
    KN_, KS_, P_ = 1.0e6, 1.0e3, 1.0e3
    u_ref = -P_ / (KN_ + KS_)
    assert abs(u - u_ref) / abs(u_ref) < 1.0e-6, (
        f"oriented flush contact equilibrium {u!r} != closed form {u_ref!r} "
        f"(a normal pointing the wrong way would give {-P_ / KS_!r}).\n{out}")
    assert abs(f - KN_ * abs(u)) / (KN_ * abs(u)) < 1.0e-6, (
        f"contact force {f!r} != kn*penetration {KN_ * abs(u)!r}.\n{out}")


def test_2d_orientation_vote_resolves():
    """THE POSITIVE CONTROL on the vote: when the vote is NOT degenerate it must
    be USED, and used correctly -- otherwise "abort on degenerate" could be
    implemented as "abort on everything without -outward", which would pass the
    test above while removing a capability.

    The 3D idiom this copies is test_adr39_contact_p2b.py::
    test_p2b_auto_orientation_no_outward: a slave clearly ABOVE the master,
    flung down explicitly, with NO `-outward`.  The reference-configuration
    surface centroids are then separated by gap0 along +y, the vote resolves to
    +y, and the block must rebound.  A wrong-signed vote lets it pass straight
    through; a refused vote aborts.

    Explicit, not static, on purpose: with a genuine positive initial gap the
    static problem is singular until contact engages, which is why the shipped
    3D control is explicit too.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -0.5, 0.0)
    ops.node(102, 0.5, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    gap0, m, v0 = 0.01, 1.0, -2.0
    ops.node(1, 0.0, gap0)
    ops.mass(1, m, m)
    ops.fix(1, 1, 0)                        # free y only; the pair is rank-1
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0)    # NO -outward: the vote decides
    ops.setNodeVel(1, 2, v0, "-commit")
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(m / KN)      # half the contact CFL
    max_pen, v_out = 0.0, 0.0
    # 1000 steps at dt = 1e-3 covers the whole event with room to spare: the
    # slave reaches the master at t = gap0/|v0| = 5 ms (step 5) and is clear of
    # it again well before step 50.
    for _ in range(1000):
        assert ops.analyze(1, dt) == 0, "the auto-oriented 2D drop diverged"
        y = gap0 + ops.nodeDisp(1, 2)
        vy = ops.nodeVel(1, 2)
        if y < 0.0:
            max_pen = max(max_pen, -y)
        if y > 0.5 * gap0 and vy > 0.0:
            v_out = vy
    e = v_out / abs(v0)
    pen_ref = abs(v0) * math.sqrt(m / KN)   # the deck's impact penalty depth
    assert e > 0.9, (
        f"restitution {e:.4f}: the resolved vote did not produce an upward "
        "normal (a flipped or refused vote gives pass-through, e = 0)")
    assert max_pen < 2.0 * pen_ref, (
        f"penetration {max_pen:.3e} vs the impact depth {pen_ref:.3e}")


def test_2d_outward_two_components():
    """The 2-COMPONENT `-outward ox oy` parses AND orients (ADR-85 How/2, the
    T1b parser branch -- T0 left the option reading exactly 3 doubles, so a 2D
    deck writing `-outward 0 1` got an unnamed 3D-signature error).

    Parsing is the cheap half.  The load-bearing half is that the two components
    mean what they say, so all three legs share one deck and differ only in the
    orientation:

      A  `-outward 0 1`, slave seeded just BELOW  -> penetrating -> ACTIVE
      B  `-outward 0 -1` on deck A                -> the slave is on the allowed
                                                     side -> INACTIVE (the
                                                     spring alone carries P)
      C  the mirror of A (slave above, load up, `-outward 0 -1`) -> ACTIVE

    A and B differ by three orders of magnitude (KN/KS = 1e3), so "inactive"
    cannot be confused with "slightly softer".  C == -A EXACTLY by mirror
    symmetry, which is a stronger statement than either leg alone: it says the
    sign is applied to the derived normal, not baked into a hard-coded +y.
    """
    KS, P = 1.0e3, 1.0e3

    def _leg(outward, seed_sign, load_sign):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.node(101, -0.5, 0.0)
        ops.node(102, 0.5, 0.0)
        ops.fix(101, 1, 1)
        ops.fix(102, 1, 1)
        y0 = seed_sign * SEED
        ops.node(1, 0.0, y0)
        ops.fix(1, 1, 0)
        ops.node(2, 0.0, y0)                # spring ground, coincident
        ops.fix(2, 1, 1)
        ops.uniaxialMaterial("Elastic", 1, KS)
        ops.element("zeroLength", 1, 2, 1, "-mat", 1, "-dir", 2)
        ops.contactSurface(10, "-master", 2, 101, 102)
        ops.contactSurface(20, "-slave", 1)
        ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", outward[0], outward[1])
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(1, 0.0, load_sign * P)
        _static(nsteps=1, what="-outward 2-component")
        return ops.nodeDisp(1, 2)

    u_a = _leg((0.0, 1.0), -1.0, -1.0)      # active from below
    u_b = _leg((0.0, -1.0), -1.0, -1.0)     # same deck, normal flipped
    u_c = _leg((0.0, -1.0), +1.0, +1.0)     # mirror of A

    u_a_ref = (KN * SEED - P) / (KN + KS)   # kn*(seed - u) = P + ks*u
    u_b_ref = -P / KS                       # spring only
    assert abs(u_a - u_a_ref) / abs(u_a_ref) < 1.0e-8, (
        f"`-outward 0 1` equilibrium {u_a!r} != closed form {u_a_ref!r}")
    assert abs(u_b - u_b_ref) / abs(u_b_ref) < 1.0e-9, (
        f"`-outward 0 -1` must leave the pair INACTIVE (the slave is on the "
        f"allowed side): got {u_b!r}, spring-only is {u_b_ref!r}")
    assert abs(u_a) < 1.0e-2 * abs(u_b), (
        f"flipping -outward changed nothing: {u_a!r} vs {u_b!r} -- the two "
        "components were parsed but not used")
    assert abs(u_c + u_a) <= 1.0e-12 * abs(u_a), (
        f"the mirrored deck is not the mirror image: {u_c!r} vs {-u_a!r} -- the "
        "orientation is baked into an axis instead of applied to the derived "
        "normal")


# --------------------------------------------------------------------------
# 3. -geomtan is a DECLARED no-op in 2D (deferred, not zero) -- G-T1b(e)
# --------------------------------------------------------------------------
def _patch_state(geomtan):
    """Solve the patch deck and return every nodal displacement as an EXACT
    IEEE-754 hexfloat.  float.hex() and not repr(): the claim is byte-identity,
    and hex() is the contact_dump harness's own convention for exactly that."""
    P = 2.0e2
    m = _patch_model(ntop=5, nbot=2, geomtan=geomtan)
    _load_uniform(m, P)
    _static(nsteps=2, what=f"-geomtan={geomtan} patch")
    tags = m["base"] + m["mid"] + m["slave"] + m["top"]     # explicit order
    return [(t, d, ops.nodeDisp(t, d).hex()) for t in tags for d in (1, 2)]


def test_2d_geomtan_noop(capfd):
    """`-geomtan` on a 2D pair is ACCEPTED, does nothing, and SAYS SO.

    TWO assertions, and the second is the one the merge panel added.

    (a) BYTE-IDENTITY.  The 2D geometric normal term is deferred, not
        implemented, so the flag must not move a single bit of the solved
        state.  Compared as exact IEEE-754 hexfloats over every nodal
        displacement.  A non-identical result is not ambiguous: adding a term
        that rounds to +-0.0 leaves any nonzero operand unchanged, so a diff
        here means a term was actually formed.

    (b) THE WARNING.  A no-op is only acceptable if the user is told.  Silence
        means someone writes `-geomtan`, believes they bought the consistent
        tangent, and gets kn*B'B -- which is ADR-76's `-initial` finding
        replayed on a different flag, and that one cost a whole ADR to find.
        So the engine emits a warnOnce naming the option, and this asserts the
        keyword appears.

    NOT the old claim.  This test used to say the term is "analytically ZERO on
    straight segments".  That is false and the module docstring now records the
    correction: the T1a FD gate proved B is the exact FIRST variation of the
    gap; the curvature term in the SECOND variation is nonzero even for a
    straight segment, because the segment's unit normal rotates as its nodes
    move.  The no-op is a decision, and byte-identity gates the decision.

    WARNONCE ORDERING HAZARD.  The warning is latched per process, and this
    battery runs a whole file in one interpreter.  The geomtan leg therefore
    runs FIRST here, and this is the only 2D `-geomtan` declaration in the
    battery -- if another is ever added ahead of it in collection order the
    latch will already have fired and (b) will fail spuriously.  The failure
    message says so, and the shipped precedent for exactly this trap is
    tests/test_contact_review_p5_percontact.py::test_p5_warning_latch_is_per_contact.
    """
    capfd.readouterr()                       # drop anything earlier in this test
    got = _patch_state(geomtan=True)         # FIRST: the warnOnce latch is fresh
    txt = "".join(capfd.readouterr()).lower()
    ref = _patch_state(geomtan=False)

    assert len(ref) == len(got)
    diffs = [(t, d, a, b) for (t, d, a), (_t, _d, b) in zip(ref, got) if a != b]
    assert not diffs, (
        "-geomtan changed the 2D solution.  T1b defers the 2D geometric normal "
        "term, so the flag is a declared no-op and any difference means a term "
        "was formed. First few (node, dof, without, with): "
        f"{diffs[:5]}")
    assert "geomtan" in txt, (
        "no warning naming `geomtan` was emitted for a 2D pair that declared it. "
        "A flag that silently does nothing is the ADR-76 `-initial` failure "
        "mode.  (If this fails while the warning demonstrably works, check "
        "whether an earlier test in this process already tripped the warnOnce "
        f"latch -- see this test's docstring.)  Captured tail:\n{txt[-800:]}")


# --------------------------------------------------------------------------
# 3b. the INITIAL-tangent arm (addKiToTang) -- ADR-76's lesson, applied
# --------------------------------------------------------------------------
def test_2d_patch_initial_tangent():
    """The patch deck under `ModifiedNewton -initial` reaches the SAME answer.

    WHY THIS ROW EXISTS.  The 2D arm of `addKiToTang` is not reached by any
    other test in this battery: every other case here uses full Newton, which
    goes through `addKtToTang`.  ADR-76 is the fork's recorded lesson about
    exactly this shape of gap -- `NDMaterial::getInitialTangent()` DEFAULTS to
    `getTangent()`, so `-initial` was silently full Newton on most solid models
    and nobody noticed, because a silently-wrong tangent still converges to the
    right answer whenever it happens to be the right tangent.  An UNPOPULATED
    2D initial tangent does not have that luxury: the contact stiffness here is
    kn = 1e6 against a block stiffness of ~1e4, so an initial tangent missing
    the contact block has a modified-Newton iteration matrix with spectral
    radius far above 1 and the run does not converge at all.  This test is
    therefore loud in both directions.

    THE ASSERTIONS.  The deck is linear once the active set is fixed (penalty
    contact, seeded so every pair is active from step 1), so `-initial` and full
    Newton must reach the SAME state, not merely a similar one:

      * convergence itself -- the iteration count is raised to 200 because
        modified Newton on a stiff penalty is linearly convergent where full
        Newton is exact, not because the answer is expected to drift;
      * the transmitted-load identity, sum(f) == P at 1e-6 (equilibrium);
      * every per-station force against the same independent reference the
        Newton row uses, at 1e-3 -- so an initial tangent that converges to a
        DIFFERENT load split (the plausible-wrong outcome, e.g. a contact block
        assembled with a stale or zeroed B) is caught rather than absorbed.
    """
    P = 2.0e2
    m = _patch_model(ntop=5, nbot=2)
    _load_uniform(m, P)
    _static(nsteps=2, tol=1.0e-10, iters=200, what="2D NTS patch (-initial)",
            algo=("ModifiedNewton", "-initial"))

    f = [ops.ladrunoContactForce(t) for t in m["slave"]]
    tot = sum(f)
    ref = _reference_interface_forces(m["nbot"], m["ntop"], P)
    print(f"\n[2D NTS patch -initial] forces {[float('%.6g' % v) for v in f]} "
          f"vs reference {[float('%.6g' % v) for v in ref]}")

    assert tot > 0.0, (
        f"-initial: the interface transmitted NOTHING (per-node {f}) on a "
        "converged solve")
    assert abs(tot - P) / P < 1.0e-6, (
        f"-initial: summed contact force {tot:.8e} != applied {P:.8e}")
    for j, (got, want) in enumerate(zip(f, ref)):
        assert abs(got - want) / abs(want) < 1.0e-3, (
            f"-initial: slave station {j} force {got:.8e} vs reference "
            f"{want:.8e} -- the initial-tangent arm converges to a different "
            "load split than the tangent arm")


# --------------------------------------------------------------------------
# 4. auto-kn on 2-DOF/node plane elements -- G-T1b(f)
# --------------------------------------------------------------------------
def _autokn_block(Emod, P, kn="auto"):
    """The p2b2b fixture in 2D: one fixed-base PlaneStrain quad whose TOP EDGE
    is a deformable master, with a slave node just penetrated at mid-span.

    `auto` is a POSITIONAL literal in the kn slot (`contact tag ms ss auto ...`,
    OpenSeesOutputCommands.cpp:527), not a `-kn auto` flag -- the same spelling
    tests/_testbed/contact_dump.py's deck 3 uses.

    Returns (penetration, contact force).  Penetration is measured against the
    DEFORMED master segment interpolated at the slave's x, so the deformable
    master's own settlement is not counted as gap.
    """
    L = 1.0
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, Emod, NU)
    for t, (x, y) in zip((1, 2, 3, 4), ((0.0, 0.0), (L, 0.0), (L, L), (0.0, L))):
        ops.node(t, x, y)
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    ops.fix(3, 1, 0)                        # uniaxial: x fixed on the top edge
    ops.fix(4, 1, 0)
    ops.element("quad", 1, 1, 2, 3, 4, THICK, "PlaneStrain", 1)
    ops.node(99, L / 2.0, L - SEED)
    ops.fix(99, 1, 0)                       # free y only
    ops.contactSurface(10, "-master", 2, 4, 3)
    ops.contactSurface(20, "-slave", 99)
    ops.contact(1, 10, 20, kn, "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(99, 0.0, -P)
    _static(nsteps=1, tol=1.0e-10, iters=50, what="2D auto-kn")
    y_master = L + 0.5 * (ops.nodeDisp(3, 2) + ops.nodeDisp(4, 2))
    y_slave = (L - SEED) + ops.nodeDisp(99, 2)
    return y_master - y_slave, ops.ladrunoContactForce(99)


def test_2d_autokn():
    """`auto` sizes a positive, usable kn from a 2-DOF/node plane element, and
    the pair still transmits.

    ADR-85 Where flags this as more than the advertised `3*nn -> ndm*nn`
    generalization: `ladrunoResolveAutoKn`'s coordinate gather is an unguarded
    3-loop (the same getCrds()(2) out-of-bounds class as the historical
    incident) and it calls the tri/quad-only
    `LadrunoContactProjection::normalOriented`.  Both need 2D siblings, and the
    failure path is FATAL -- so this is the first thing a 2D deck that omits an
    explicit kn hits.

    Three claims, the middle one non-circular:

    1. it converges and lands a WORKING penalty: the penetration is a small
       fraction of the element size (not O(L) = inert, not locked).  Loaded at
       the same 0.1%-of-modulus contact pressure the 3D gate uses, scaled to 2D
       as a force per unit thickness over a length: P = 1e-3 * E * L * THICK;
    2. kn is proportional to E, so the penetration scales as 1/E.  This is
       formulation-INDEPENDENT physics -- it holds for any reduction of the
       element stiffness through the normal, so it cannot be satisfied by a
       constant fallback that merely happens to be positive;
    3. the transmitted force equals the applied load (equilibrium identity),
       which is also the 2D check that `ladrunoContactForce` is live on this
       lane.

    NOT gated here: the ABSOLUTE kn value against the f_si*mean(n'Kn) oracle.
    The 3D battery can do that because proto_p2b2b_autokn.py supplies the same
    hex stiffness the C++ reads; there is no 2D sibling oracle, and hand-rolling
    a quad stiffness inside a test would gate the test's algebra, not the
    engine's.  Recorded as a gap rather than approximated.
    """
    L, Emod = 1.0, 2.0e4
    P = 1.0e-3 * Emod * L * THICK
    pen, force = _autokn_block(Emod, P)
    assert pen > 0.0, f"no penetration: `auto` produced an inert pair (pen={pen:.3e})"
    kn_auto = P / pen
    print(f"\n[2D auto-kn] kn_auto = P/pen = {kn_auto:.6e} (pen/L = {pen / L:.3e})")
    assert pen / L < 0.10, (
        f"penetration {pen / L:.3e} of the element size -- `auto` sized a "
        "penalty that is far too soft (or the pair is barely engaged)")
    assert abs(force - P) / P < 1.0e-6, (
        f"transmitted force {force:.8e} != applied {P:.8e}: auto-kn assembled "
        "but the pair does not carry the load")

    P2 = 1.0e2                              # fixed load, two moduli
    pen1, _f1 = _autokn_block(1.0e4, P2)
    pen3, _f3 = _autokn_block(3.0e4, P2)
    ratio = pen1 / pen3
    assert abs(ratio - 3.0) / 3.0 < 0.01, (
        f"pen(E)/pen(3E) = {ratio:.4f}, expected 3 -- kn is not proportional to "
        "the element stiffness, so `auto` is not reading it")


# --------------------------------------------------------------------------
# 5. T1b must not WIDEN T0's guards
# --------------------------------------------------------------------------
# Both rows are re-run through the T0 file's own child + contract table, not
# re-implemented here: a second copy of the deck could drift from the one the T0
# matrix pins, and then "still refused" would be asserting something else.
def test_2d_ndf3_still_refused():
    """A `-ndm 2 -ndf 3` node on an NTS surface is STILL refused after T1b.

    The temptation this guards is specific.  T1b's job is to make ndm = 2 work
    everywhere the handler previously hard-coded 3, and the six shipped
    `getNumberDOF() != 3` guards are exactly such sites -- ADR-85 How/8 says
    they become `!= ndm`, an EQUALITY, and a wiring pass that relaxed them to
    `>= ndm` or dropped them while chasing the stride-3 hardcodes would look
    like progress.  It is the historical out-of-bounds reproducer:
    FE_Element::setID() packs each node's full DOF_Group ID sequentially, so a
    third DOF pushes rotation equations into the contact's myID slots and
    shifts every later slot.  Silent mis-assembly, on a deck that now RUNS.
    """
    _assert_refused("ndf3_frame_node_on_nts")


def test_2d_mortar_still_refused():
    """A well-formed all-2D MORTAR pair is STILL refused by name, pointing at T3.

    The NTS and mortar lanes share the handler, the surfaces, the orientation
    machinery and (after T1b) the 2D branch of nearly everything except the
    interval integrator.  That is precisely why this row is re-asserted from the
    T1b side: opening NTS must not let a `-mortar` declaration fall through the
    same widened gate and reach a lane whose 2D integrator does not exist yet.
    """
    _assert_refused("mortar_2d_not_yet_supported")
