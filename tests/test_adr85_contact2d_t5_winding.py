"""ADR-85 F1 -- `-outward winding`: the 2D NTS chain declares its own side.

WHAT THIS FILE IS
    T1b gave the 2D lanes ONE sigma per interface, resolved by an
    interface-level reference vote over the surface centroids
    (`ladruno2DOrientationVote`, SRC/analysis/handler/LadrunoContactHandler.cpp).
    That vote refuses two whole classes of deck BY NAME, and both refusals are
    correct as far as they go:

      * FLUSH interfaces -- coincident master/slave centroids.  The masonry
        joint, the footing seated on soil, every zero-gap deck the 2D lane was
        built for.  Before F1 they had to pass `-outward ox oy` or abort
        (LadrunoContact2D_guide.md, the "Orientation" section).
      * STRONGLY CURVED masters -- the vote SPLITS, because no single direction
        vector agrees with every segment's own perpendicular.  A CLOSED-LOOP
        master (a ring, a full indenter profile) can therefore not be oriented
        at all, even WITH `-outward`.

    F1 adds a third form of the option: `-outward winding`.  It fixes
    sigma = +1, i.e. every master segment's normal is `perp(t) = (-t_y, t_x)`
    of ITS OWN travel direction, so the slave lies to the LEFT of the chain's
    traversal.  The centroid vote is bypassed entirely.  Nothing else changes:
    sigma was already one scalar and the engine already evaluates
    `n = sigma*perp(t_k)/L_k` per segment from that segment's own tangent
    (LadrunoContact2DKernel.h::projectSegment2D), so the kernel, the FE
    adapter, the concave-vertex classifier, the D4 end-cap and
    `ladrunoResolveAutoKn2D` are untouched.

THE GUARD THIS FILE EXISTS TO PIN (read before touching any row here)
    The handle()-time chain-integrity scan is **vacuous for disjoint
    segments**: it skips any node used only once (`if (v.size() <= 1)
    continue;`), and a DISJOINT master -- two separate runs with a hole between
    them -- is explicitly legal and documented as such
    (LadrunoContact2D_guide.md, the stride-2 pair-list section).  Today the
    SPLIT-VOTE refusal is what catches a mis-wound second run.  Bypassing the
    vote would convert that named FATAL into a SILENT WRONG-SIDE NORMAL -- the
    ADR-78 P0 shape, converged and balanced and wrong.  So `-outward winding`
    carries its OWN connectivity requirement: the master must form ONE
    connected chain, refused by name otherwise
    (`winding_disjoint_master` below).  A winding mode without that row is a
    silent-wrong feature; do not retire it.

SCOPE: NTS ONLY (deliberate asymmetry, documented)
    The chain-integrity scan is NTS-only.  The 2D MORTAR lane has no chain scan
    at all -- it goes Lref -> vote -> the double loop -- so the connectivity
    invariant winding leans on is not established there, and hoisting the scan
    would newly FATAL permuted or reversed mortar masters that ship ACCEPTED
    today.  That is a behavior change on a shipped lane and needs its own ADR
    decision, so: winding is refused by name on the mortar lane
    (`winding_on_mortar`), flush MORTAR interfaces still abort, and mortar
    users still pass `-outward ox oy`.

ONE ROW HERE ASSERTS BEHAVIOUR THAT IS WRONG
    `test_coarse_closed_loop_transmits_nothing` pins a DISCLOSED LIMITATION --
    a closed loop whose facets are coarse relative to its own diameter is
    silently inert.  It is a pin, not a blessing; read its docstring (and the
    LEDGER_quirks row it names) before touching it, and delete it if the
    behaviour is ever fixed.

STRUCTURE (the tests/test_adr85_contact2d_t1b_guards.py pattern)
    Acceptance rows run IN-PROCESS (`from _testbed import ops`, `ops.wipe()`
    between decks).  EVERY FATAL row runs in a CHILD interpreter via
    `_run_child`/`_assert_refused`, because a contact refusal is not an
    ordinary Python exception: a declaration refusal returns -1 (an
    `OpenSeesError` in serial, an `MPI_Abort` under np > 1) and a handler-time
    one surfaces later still, as a nonzero `analyze()` after `domainChanged()`
    fails.  `ENGINE_DIR` is imported from
    tests/test_adr85_contact2d_t0_refusals.py so the child binds the SAME
    binary the parent did; that file's module docstring documents both Windows
    traps (the PYTHONPATH environment-block trap and the stdin=DEVNULL Popen
    trap).

MESSAGE CONTRACT
    Exit code alone is too weak -- a deck that merely fails to converge also
    exits nonzero -- so every falsifier asserts SUBSTRINGS on the child's
    combined output.  `ADR-85` is required on every row (the shipped precedent
    that a refusal names the ADR that owns it); the second substring names the
    specific rule.
"""
import math
import os
import subprocess
import sys

import pytest

from _testbed import ops
from test_adr85_contact2d_t0_refusals import ENGINE_DIR

pytestmark = [pytest.mark.zone_a]


# --------------------------------------------------------------------------
# rig constants -- shared by the in-process rows and the child template.
# --------------------------------------------------------------------------
KN = 1.0e6
KS = 1.0e3
KT = 1.0e5
P = 1.0e3
SEED = 1.0e-8
EPSN = 1.0e6


# ==========================================================================
# IN-PROCESS ACCEPTANCE ROWS
# ==========================================================================
def _solve():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 40, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    return ops.analyze(1)


def _flush_leg(with_contact, master_order=(101, 102), seed_sign=-1.0, load_sign=-1.0):
    """The T1b flush deck (tests/test_adr85_contact2d_t1b_guards.py, case
    `vote_seed_noise_*`) with the orientation supplied by the MASTER LISTING
    ORDER instead of by an `-outward` vector.

    Master line y = 0 from x = -0.5 to +0.5; the slave sits `seed_sign*SEED`
    off it and is held by a soft y-spring so the deck stays solvable with or
    without a live pair (that is what makes leg 2 a quantitative control
    rather than a singular matrix).  Returns (rc, u_y, contactForce).
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -0.5, 0.0)
    ops.node(102, 0.5, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.node(1, 0.0, seed_sign * SEED)
    ops.fix(1, 1, 0)                       # x fixed, y free
    ops.node(2, 0.0, seed_sign * SEED)     # spring ground, coincident
    ops.fix(2, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, KS)
    ops.element("zeroLength", 1, 2, 1, "-mat", 1, "-dir", 2)
    if with_contact:
        ops.contactSurface(10, "-master", 2, *master_order)
        ops.contactSurface(20, "-slave", 1)
        ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", "winding")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, load_sign * P)
    rc = _solve()
    if rc != 0:
        return rc, 0.0, 0.0
    f = ops.ladrunoContactForce(1) if with_contact else 0.0
    return 0, ops.nodeDisp(1, 2), f


def test_flush_interface_runs_and_transmits_under_winding():
    """THE HEADLINE ROW.  The flush deck that FATALs today runs under winding
    -- and does not merely run, it TRANSMITS the closed-form force.

    The identical deck WITHOUT `-outward winding` draws the named
    degenerate-vote abort; that verdict is pinned in the child row
    `flush_no_winding_refused` below, so the pair of rows differs by exactly
    one token and the acceptance here is attributable to winding and nothing
    else.

    "It ran" is deliberately NOT the assertion: a 2D pair that assembles and
    transmits nothing converges and balances its reactions just as
    convincingly as a live one (the ADR-78 P0 incident).  So the deck is
    solved TWICE and both legs are checked against closed forms, exactly as
    tests/test_adr85_contact2d_t0_refusals.py::test_pair_2d_now_live does:

        leg 2 (no pair)  u = -P/KS                            (pure spring)
        leg 1 (pair)     KN*(SEED - u) = P + KS*u
                         =>  u = (KN*SEED - P)/(KN + KS)      (spring || penalty)
        force            f = KN*(SEED - u)

    Master listed 101 -> 102, i.e. t = +x, so perp(t) = (0, +1): the slave
    (nominally above the line, seeded SEED into it) is to the LEFT of the
    traversal, which is precisely what sigma = +1 means.
    """
    rc_c, u_c, f_c = _flush_leg(True)
    rc_f, u_f, _ = _flush_leg(False)
    assert rc_c == 0, "the flush deck must now RUN under -outward winding"
    assert rc_f == 0

    u_f_ref = -P / KS
    u_c_ref = (KN * SEED - P) / (KN + KS)
    f_c_ref = KN * (SEED - u_c_ref)

    assert abs(u_f - u_f_ref) / abs(u_f_ref) < 1.0e-9, (
        f"the contact-free control is not the pure spring: {u_f!r} vs "
        f"{u_f_ref!r} -- the deck changed, so leg 2 no longer isolates the pair")
    assert abs(u_c - u_c_ref) / abs(u_c_ref) < 1.0e-8, (
        f"winding pair equilibrium {u_c!r} != closed form {u_c_ref!r}: the pair "
        f"assembled but does not carry kn")
    assert abs(u_c) < 1.0e-2 * abs(u_f), (
        f"the winding pair transmitted (almost) nothing: |u| with contact "
        f"{abs(u_c):.6e} vs free {abs(u_f):.6e} -- the ADR-78 P0 silent-drop "
        f"signature")
    assert f_c > 0.0, "ladrunoContactForce reports 0 on a demonstrably loaded pair"
    assert abs(f_c - f_c_ref) / f_c_ref < 1.0e-6, (
        f"ladrunoContactForce {f_c!r} != kn*penetration {f_c_ref!r}")


def test_winding_sign_flip_is_equal_and_opposite():
    """THE SIGN-FLIP CONTROL.  Reversing the master's listing order reverses
    the normal, and nothing else.

    Two checks, and the second is the one that would catch a winding mode
    wired to a hardcoded direction rather than to the segment's own tangent:

      1. MIRROR.  Deck A (master 101->102, slave seeded BELOW, pushed DOWN)
         and deck B (master 102->101, slave seeded ABOVE, pushed UP) are exact
         y-mirrors of each other.  Their answers must be exactly opposite and
         their transmitted forces equal.  Tolerance 1e-12 relative: the
         mirrored deck is the same linear system with a sign-flipped RHS, so
         the only floating-point difference is sign propagation through the
         geometry; 1e-12 is far above that and far below any real asymmetry.

      2. DISARM.  Deck C keeps deck A's geometry and load but lists the master
         BACKWARDS (102->101).  perp(t) now points -y, the slave's SEED offset
         becomes a POSITIVE gap, the pair is inactive, and the answer collapses
         onto the free-spring control with ZERO contact force.  This is the
         hazard the winding mode creates for the caller -- a wrongly wound
         master runs and gives a converged wrong answer -- and it is recorded
         here as a fact, not a defect: the engine cannot distinguish it from a
         genuinely separated interface, which is exactly why the connectivity
         guard (`winding_disjoint_master`) is mandatory and why the declaring
         layer above must own the "master faces the slave" check.
    """
    rc_a, u_a, f_a = _flush_leg(True, master_order=(101, 102),
                                seed_sign=-1.0, load_sign=-1.0)
    rc_b, u_b, f_b = _flush_leg(True, master_order=(102, 101),
                                seed_sign=+1.0, load_sign=+1.0)
    assert rc_a == 0 and rc_b == 0, "both winding orientations must run"
    assert abs(u_b + u_a) <= 1.0e-12 * abs(u_a), (
        f"the two winding sides are not mirror images: u_A={u_a!r} u_B={u_b!r} "
        f"-- sigma=+1 must give perp(t) of each segment's OWN tangent")
    assert abs(f_b - f_a) <= 1.0e-12 * f_a, (
        f"mirrored decks transmit different force: {f_a!r} vs {f_b!r}")

    rc_c, u_c, f_c = _flush_leg(True, master_order=(102, 101),
                                seed_sign=-1.0, load_sign=-1.0)
    assert rc_c == 0
    assert f_c == 0.0, (
        f"a BACKWARDS-wound master still transmitted {f_c!r}: the normal did "
        f"not follow the listing order, so `-outward winding` is not reading "
        f"the segment tangent")
    assert abs(u_c - (-P / KS)) / (P / KS) < 1.0e-9, (
        f"the backwards-wound leg {u_c!r} is not the free-spring answer "
        f"{-P / KS!r}; the pair is neither armed nor cleanly inert")


def test_negative_outward_vector_still_parses_after_the_keyword_peek():
    """THE REGRESSION GUARD ON THE PARSE ITSELF.  `-outward` with a NEGATIVE
    component still reads as two doubles.

    F1 inserts a one-token peek in front of the shipped 2-double read, and the
    peek CONSUMES its token on both interpreter paths.  If the un-read were
    wrong by one, or if a negative number were misclassified (on the Tcl path
    every token is a string, and `-1.0` starts with `-` exactly like an option
    flag -- the hazard this file's `ladrunoCountLeadingNumbers` sibling
    documents at length), this deck would die in the double read.  No shipped
    2D test passes a negative `-outward` component, so nothing else covers it.

    The deck is the flush rig mirrored: slave seeded SEED ABOVE the master and
    pulled UP, oriented with `-outward 0.0 -1.0`.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -0.5, 0.0)
    ops.node(102, 0.5, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.node(1, 0.0, SEED)
    ops.fix(1, 1, 0)
    ops.node(2, 0.0, SEED)
    ops.fix(2, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, KS)
    ops.element("zeroLength", 1, 2, 1, "-mat", 1, "-dir", 2)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", -0.0, -1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, P)
    assert _solve() == 0, "a negative -outward component must still parse and run"

    u_ref = -(KN * SEED - P) / (KN + KS)
    u = ops.nodeDisp(1, 2)
    assert abs(u - u_ref) / abs(u_ref) < 1.0e-8, (
        f"the mirrored explicit-vector deck gave {u!r}, closed form {u_ref!r}")
    assert ops.ladrunoContactForce(1) > 0.0


def test_chained_multisegment_master_runs_under_winding():
    """THE CONTROL ON THE CONNECTIVITY GUARD.  A properly chained three-segment
    master under winding must RUN -- otherwise "refuse every winding master"
    would pass the `winding_disjoint_master` row for free.

    Same four nodes as the disjoint child case; the ONLY difference is six
    tags (chained pairs) instead of four (the guide's "silently legal" holed
    declaration).  The slave sits over the middle segment, so the closed form
    is deck A's verbatim.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for i, x in enumerate((-1.5, -0.5, 0.5, 1.5)):
        ops.node(101 + i, x, 0.0)
        ops.fix(101 + i, 1, 1)
    ops.node(1, 0.0, -SEED)
    ops.fix(1, 1, 0)
    ops.node(2, 0.0, -SEED)
    ops.fix(2, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, KS)
    ops.element("zeroLength", 1, 2, 1, "-mat", 1, "-dir", 2)
    ops.contactSurface(10, "-master", 2, 101, 102, 102, 103, 103, 104)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", "winding")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, -P)
    assert _solve() == 0, "a chained 3-segment winding master must run"

    u_ref = (KN * SEED - P) / (KN + KS)
    u = ops.nodeDisp(1, 2)
    f = ops.ladrunoContactForce(1)
    assert abs(u - u_ref) / abs(u_ref) < 1.0e-8, f"{u!r} != {u_ref!r}"
    assert abs(f - KN * (SEED - u_ref)) / (KN * (SEED - u_ref)) < 1.0e-6


def _ring(nseg, R=1.0, winding=True):
    """A regular `nseg`-gon master ring, listed CLOCKWISE so that left-of-travel
    (sigma = +1) is OUTWARD on every segment, with a slave seeded SEED inside
    the midpoint of segment 0 and held by a soft x-spring.

    Vertices are placed at 180 + 180/nseg - 360*k/nseg degrees, which puts
    segment 0 on a VERTICAL left face at x = -R*cos(pi/nseg): the deck is then
    the x-rotation of the flush deck, so its closed form is the same one
    mirrored,

        u = (P - KN*SEED)/(KN + KS)        f = KN*(SEED + u)

    A closed loop has no open terminal (so no D4 end-cap) and every corner is
    convex seen from outside (turn < 0, so no concave vertex pair): the ring
    exercises the pure segment path plus the closed-loop seam.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    half = 180.0 / nseg
    pts = [(R * math.cos(math.radians(180.0 + half - 360.0 * k / nseg)),
            R * math.sin(math.radians(180.0 + half - 360.0 * k / nseg)))
           for k in range(nseg)]
    for i, (x, y) in enumerate(pts):
        ops.node(101 + i, x, y)
        ops.fix(101 + i, 1, 1)
    x0 = pts[0][0] + SEED
    ops.node(1, x0, 0.0)
    ops.fix(1, 0, 1)                       # x free, y fixed
    ops.node(2, x0, 0.0)                   # spring ground, coincident
    ops.fix(2, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, KS)
    ops.element("zeroLength", 1, 2, 1, "-mat", 1, "-dir", 1)
    tags = []
    for k in range(nseg):
        tags += [101 + k, 101 + ((k + 1) % nseg)]
    ops.contactSurface(10, "-master", 2, *tags)
    ops.contactSurface(20, "-slave", 1)
    args = [1, 10, 20, KN, 0.0, 0.0]
    if winding:
        args += ["-outward", "winding"]
    ops.contact(*args)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, P, 0.0)
    rc = _solve()
    if rc != 0:
        return rc, 0.0, 0.0
    return 0, ops.nodeDisp(1, 1), ops.ladrunoContactForce(1)


def test_closed_loop_master_runs_under_winding():
    """THE CAPABILITY ROW.  A CLOSED-LOOP master runs, and transmits.

    A closed loop already passes the chain-integrity scan -- the wrap
    `headSeg == (tailSeg + 1) % nSeg` is explicitly legal -- so the only thing
    that ever blocked it was the SPLIT-VOTE refusal, which is pinned in the
    child row `closed_loop_no_winding_refused` below (the same deck, minus the
    winding token, still aborts by name).

    TWELVE segments, not four, and the number is load-bearing -- see
    test_coarse_closed_loop_transmits_nothing below and the LEDGER_quirks row
    it names.  A ring whose facets are coarse relative to its own diameter is
    silently inert; twelve is comfortably past the measured transition (n = 9
    at R = 1 with the default `-cell`).
    """
    rc, u, f = _ring(12)
    assert rc == 0, (
        "a closed-loop 2D master must RUN under -outward winding (it already "
        "passes chain integrity; only the split vote blocked it)")

    u_ref = (P - KN * SEED) / (KN + KS)
    f_ref = KN * (SEED + u_ref)
    assert abs(u - u_ref) / abs(u_ref) < 1.0e-8, (
        f"closed-loop winding equilibrium {u!r} != closed form {u_ref!r}")
    assert f > 0.0, "the closed-loop pair transmitted nothing"
    assert abs(f - f_ref) / f_ref < 1.0e-6, f"{f!r} != {f_ref!r}"


def test_coarse_closed_loop_transmits_nothing():
    """A DISCLOSED LIMITATION, pinned so it cannot drift unnoticed.

    THIS ROW ASSERTS BEHAVIOUR THAT IS WRONG.  A four-segment ring is a
    perfectly well-formed closed-loop master: it passes chain integrity, it
    passes the winding connectivity guard, it runs, it converges -- and it
    transmits EXACTLY ZERO.  Measured mechanism (attributed, not guessed; full
    write-up in LEDGER_quirks under "closed-loop 2D master: coarse rings go
    silently inert"):

      * the broad-phase grid is capped at `NX*NY*NZ <= nSeg` cells with a
        +/-1-cell candidate search (LadrunoContactBucketSort.h), so a ring with
        <= 8 segments cannot be spatially separated from its OWN far side --
        every segment is a candidate;
      * the far-side segment then arms with a body-diameter "penetration" and
        kicks the slave toward the ring's interior on the first iterate;
      * from anywhere near the interior EVERY segment projects the slave
        in-bounds, so the ordered-ownership rule ("stand down if your
        PREDECESSOR is in-bounds") cycles all the way around the closed chain
        and no segment owns it.  Zero contact stiffness, zero force, and the
        spring-only answer is a genuine self-consistent equilibrium.

    Measured transition for a regular n-gon of circumradius 1 at the default
    `-cell 1.0`: n <= 8 inert, n >= 9 exact.  Forcing one bucket (`-cell` huge)
    reproduces the inert result at ANY n, which is what identifies the broad
    phase rather than the ring itself as the discriminator.

    IF THIS ROW EVER FAILS BECAUSE THE RING NOW TRANSMITS, that is an
    IMPROVEMENT: update LadrunoContact2D_guide.md's closed-loop paragraph and
    the LEDGER_quirks row, and delete this test.  It is here so that outcome is
    noticed rather than assumed.
    """
    rc, u, f = _ring(4)
    assert rc == 0, "the coarse ring converges (that is the whole problem)"
    assert f == 0.0, (
        f"the 4-segment ring transmitted {f!r}: the documented coarse-ring "
        f"stand-down no longer reproduces -- see this test's docstring")
    assert abs(u - P / KS) / (P / KS) < 1.0e-9, (
        f"the coarse ring settled at {u!r}, which is neither the free-spring "
        f"answer {P / KS!r} nor the contact closed form -- the failure mode "
        f"changed shape and the documentation is now wrong")


def test_winding_friction_saturates_at_mu_n():
    """`-mu > 0` keys correctly on a winding interface.

    The T2 flat slip rig (tests/test_adr85_contact2d_t2_friction.py::
    _flat_slip_anchor) with `-outward 0.0 1.0` swapped for
    `-outward winding` -- the master is already listed left-to-right, so the
    declared side is unchanged and the T2 closed forms carry over verbatim:
    driven above the cone, the slave slips until the anchor takes the excess,

        x = (Q - mu*P)/ka        friction = mu*P

    This matters because friction is the one consumer that reads the normal
    TWICE (the cone is keyed off the normal traction, the slip direction off
    the tangent): a winding mode that got sigma into the gap but not into the
    friction cone would pass every frictionless row above and fail here.
    """
    mu, Q, ka = 0.5, 8.0e2, 1.0e4
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -2.0, 0.0)
    ops.node(102, 2.0, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.node(1, 0.0, -SEED)
    ops.node(900, 0.0, -SEED)
    ops.fix(900, 1, 1)
    ops.uniaxialMaterial("Elastic", 900, ka)
    ops.element("zeroLength", 900, 900, 1, "-mat", 900, "-dir", 1)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, mu, "-outward", "winding")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-9, 60, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "a winding pair with -mu must converge"

    x = ops.nodeDisp(1, 1)
    expect = (Q - mu * P) / ka
    assert abs(x - expect) / expect < 1.0e-6, (
        f"winding + friction slip {x!r} != closed form {expect!r}")
    friction = Q - ka * x
    assert abs(friction - mu * P) / (mu * P) < 1.0e-6, (
        f"the cone did not saturate at mu*N on a winding interface: "
        f"{friction!r} vs {mu * P!r}")


# ==========================================================================
# THE REFUSAL CONTRACT.  case -> substrings that MUST appear in the output.
# ==========================================================================
REFUSALS = {
    # `-outward winding` AND `-outward ox oy` on the same contact.  Two
    # declarations of the same thing; the engine never picks one.
    "winding_plus_outward": ("ADR-85", "winding"),

    # `-outward winding` on a 3D pair.  There is no 3D analogue: a facet
    # surface has no traversal.  It must be refused BY NAME rather than
    # surfacing as the generic "need ox oy oz" double-read failure, which
    # would name a layout the user never wrote.
    "winding_on_3d_pair": ("ADR-85", "winding"),

    # THE MANDATORY GUARD.  Four tags = two DISJOINT segments with a hole
    # between them (the guide's "silently legal" holed declaration).  The
    # chain-integrity scan is VACUOUS here -- every node is used exactly once,
    # so every node is skipped -- and with the vote bypassed a mis-wound
    # second run would be a silent wrong-side normal.  Refused by name.
    "winding_disjoint_master": ("ADR-85", "disjoint"),

    # Winding on the MORTAR lane.  The mortar lane has no chain-integrity scan,
    # so the one-connected-chain invariant is not established there; hoisting
    # the scan would newly FATAL mortar masters that ship accepted today.
    # Refused by name, with the asymmetry stated.
    "winding_on_mortar": ("ADR-85", "mortar"),

    # CONTROL on the headline row: the SAME flush deck as
    # test_flush_interface_runs_and_transmits_under_winding, minus the winding
    # token, still draws the degenerate-vote abort.  One token, opposite
    # verdict.
    "flush_no_winding_refused": ("ADR-85", "outward"),

    # CONTROL on the capability row: the SAME closed-loop deck, minus the
    # winding token, still draws the SPLIT-vote abort -- which is the only
    # thing that ever blocked a closed loop (chain integrity already accepts
    # the wrap).
    "closed_loop_no_winding_refused": ("ADR-85", "SPLIT"),
}


# --------------------------------------------------------------------------
# the child.  NOTE: this template is rendered by the percent operator against
# a one-key mapping, so NOTHING inside it -- code or comment -- may contain a
# bare percent conversion; use repr() and concatenation.
# --------------------------------------------------------------------------
CHILD = r'''
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

KN, KS, KT, P, SEED, EPSN = 1.0e6, 1.0e3, 1.0e5, 1.0e3, 1.0e-8, 1.0e6


def _solve():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-10, 40, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.integrator("LoadControl", 1.0)
    return ops.analyze(1)


def _tether(slave, x, y, direction, k=KS):
    ops.node(900, x, y)
    ops.fix(900, 1, 1)
    ops.uniaxialMaterial("Elastic", 900, k)
    ops.element("zeroLength", 900, 900, slave, "-mat", 900, "-dir", direction)


rc = 0

if CASE in ("winding_plus_outward", "flush_no_winding_refused"):
    # The flush deck.  `winding_plus_outward` declares BOTH forms of -outward
    # (a parse-time refusal); `flush_no_winding_refused` declares NEITHER and
    # must still draw the T1b degenerate-vote abort (the control proving the
    # in-process acceptance row is attributable to winding alone).
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -0.5, 0.0)
    ops.node(102, 0.5, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.node(1, 0.0, -SEED)
    ops.fix(1, 1, 0)
    _tether(1, 0.0, -SEED, 2)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.contactSurface(20, "-slave", 1)
    if CASE == "winding_plus_outward":
        ops.contact(1, 10, 20, KN, 0.0, 0.0,
                    "-outward", "winding", "-outward", 0.0, 1.0)
    else:
        ops.contact(1, 10, 20, KN, 0.0, 0.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, -P)
    rc = _solve()

elif CASE == "winding_on_3d_pair":
    # A well-formed 3D NTS pair.  `winding` has no 3D meaning (a facet surface
    # has no traversal), and the refusal must NAME it rather than fall through
    # to the generic 3-double read failure.
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y) in enumerate([(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]):
        ops.node(101 + i, x, y, 0.0)
        ops.fix(101 + i, 1, 1, 1)
    ops.node(1, 0.0, 0.0, -SEED)
    ops.fix(1, 1, 1, 0)
    ops.contactSurface(10, "-master", 4, 101, 102, 103, 104)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", "winding")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -P)
    rc = _solve()

elif CASE == "winding_disjoint_master":
    # THE MANDATORY GUARD.  Four master nodes, FOUR tags -- two disjoint
    # segments (101,102) and (103,104) with a hole where (102,103) should be.
    # Every node is used exactly ONCE, so the chain-integrity scan skips them
    # all and sees nothing.  With the vote bypassed, a mis-wound second run
    # would be a silent wrong-side normal, so winding must refuse this by name.
    # (The properly chained SIX-tag twin is the in-process control
    # test_chained_multisegment_master_runs_under_winding.)
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for i, x in enumerate([-1.5, -0.5, 0.5, 1.5]):
        ops.node(101 + i, x, 0.0)
        ops.fix(101 + i, 1, 1)
    ops.node(1, 0.0, -SEED)
    ops.fix(1, 1, 0)
    _tether(1, 0.0, -SEED, 2)
    ops.contactSurface(10, "-master", 2, 101, 102, 103, 104)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", "winding")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, -P)
    rc = _solve()

elif CASE == "winding_on_mortar":
    # The T0/T3 matched 2D mortar pair with `-outward winding`.  NTS-only:
    # the mortar lane establishes no chain invariant, so winding is refused by
    # name and flush mortar interfaces still need `-outward ox oy`.
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -0.5, 0.0)
    ops.node(102, 0.5, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.node(1, -0.5, -SEED)
    ops.node(2, 0.5, -SEED)
    ops.fix(1, 1, 0)
    ops.fix(2, 1, 0)
    ops.node(3, -0.5, -SEED)
    ops.node(4, 0.5, -SEED)
    ops.fix(3, 1, 1)
    ops.fix(4, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, KS)
    ops.element("zeroLength", 1, 3, 1, "-mat", 1, "-dir", 2)
    ops.element("zeroLength", 2, 4, 2, "-mat", 1, "-dir", 2)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.contactSurface(20, "-slave-segments", 2, 1, 2)
    ops.contact(1, 10, 20, "-mortar", "-epsN", EPSN, "-outward", "winding")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, -P / 2.0)
    ops.load(2, 0.0, -P / 2.0)
    rc = _solve()

elif CASE == "closed_loop_no_winding_refused":
    # The in-process closed-loop deck (the SAME 12-gon ring) MINUS the winding
    # token.  The centroid datum points -x; the left face votes +1, the right
    # face votes -1, so the vote SPLITS.  This is the ONLY thing that ever
    # blocked a closed-loop master -- chain integrity already accepts the wrap.
    import math
    NSEG, R = 12, 1.0
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    half = 180.0 / NSEG
    pts = [(R * math.cos(math.radians(180.0 + half - 360.0 * k / NSEG)),
            R * math.sin(math.radians(180.0 + half - 360.0 * k / NSEG)))
           for k in range(NSEG)]
    for i, xy in enumerate(pts):
        ops.node(101 + i, xy[0], xy[1])
        ops.fix(101 + i, 1, 1)
    x0 = pts[0][0] + SEED
    ops.node(1, x0, 0.0)
    ops.fix(1, 0, 1)
    _tether(1, x0, 0.0, 1)
    tags = []
    for k in range(NSEG):
        nxt = k + 1 if k + 1 < NSEG else 0     # NB: no bare modulo in this
        tags += [101 + k, 101 + nxt]           # template -- see _run_child
    ops.contactSurface(10, "-master", 2, *tags)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, P, 0.0)
    rc = _solve()

else:
    print("UNKNOWN CASE " + CASE)
    sys.exit(4)

if rc == 0:
    print("RAN " + CASE)
    sys.exit(0)

print("ANALYZE_REFUSED " + CASE + " rc=" + str(rc))
sys.exit(3)
'''


def _run_child(case):
    script = CHILD % {"ENGINE_DIR": ENGINE_DIR}
    env = dict(os.environ)
    env["LADRUNO_OPENSEES_QUIET"] = "1"     # the banner's UTF-8 glyphs break a
                                            # text-mode capture on cp1252
    proc = subprocess.run(
        [sys.executable, "-c", script, case],
        stdin=subprocess.DEVNULL,           # load-bearing on Windows
        capture_output=True, text=True, timeout=300,
        encoding="utf-8", errors="replace",
        env=env,
    )
    return proc, (proc.stdout or "") + "\n" + (proc.stderr or "")


def _assert_refused(case):
    proc, out = _run_child(case)
    assert proc.returncode != 0, (
        f"{case}: the child exited 0 -- the deck was ACCEPTED and ran. "
        f"ADR-85 F1 must refuse it.\n{out}")
    assert "RAN " not in out, (
        f"{case}: the deck assembled and solved without a refusal.\n{out}")
    for needle in REFUSALS[case]:
        assert needle in out, (
            f"{case}: the refusal fired but does not name itself -- missing "
            f"{needle!r}. See the REFUSALS contract in this file.\n{out}")


# --------------------------------------------------------------------------
# the falsifier matrix
# --------------------------------------------------------------------------
def test_refuse_winding_plus_outward():
    """`-outward winding` together with `-outward ox oy` -> named refusal.

    Two declarations of the same fact.  Accepting the pair would mean picking
    a precedence rule silently, and a user who wrote both has a wrong mental
    model of at least one of them.
    """
    _assert_refused("winding_plus_outward")


def test_refuse_winding_on_3d_pair():
    """`-outward winding` on a 3D pair -> named refusal, not a token error.

    The 3D arm of the `-outward` parse is otherwise textually untouched by F1;
    the keyword is intercepted BEFORE the 3-double read so the message names
    winding instead of demanding `ox oy oz`, a layout the user never wrote.
    """
    _assert_refused("winding_on_3d_pair")


def test_refuse_winding_on_disjoint_master():
    """THE MANDATORY GUARD -- a DISJOINT master under winding -> named refusal.

    This is the row the whole feature hangs on.  The chain-integrity scan skips
    any node used only once, so a disjoint master (explicitly legal, and
    documented as such) passes it vacuously.  Today the SPLIT-VOTE refusal is
    what catches a mis-wound second run; bypassing the vote without adding this
    guard converts that named FATAL into a silent wrong-side normal.

    Its falsifier is test_chained_multisegment_master_runs_under_winding: the
    same four nodes, chained, must still run.
    """
    _assert_refused("winding_disjoint_master")


def test_refuse_winding_on_mortar():
    """`-outward winding` with `-mortar` -> named refusal.

    Deliberate scope fence, not an oversight: the mortar lane runs no
    chain-integrity scan, so the one-connected-chain invariant winding relies
    on is not established there, and hoisting the scan would newly FATAL
    permuted or reversed mortar masters that ship ACCEPTED today.  Consequence
    the message must carry: flush 2D MORTAR interfaces still require
    `-outward ox oy`.
    """
    _assert_refused("winding_on_mortar")


def test_flush_deck_without_winding_is_still_refused():
    """CONTROL: the acceptance row's deck, minus the winding token, still
    aborts on the degenerate centroid vote.  Without this the headline row
    could be passing because the vote quietly started resolving flush decks.
    """
    _assert_refused("flush_no_winding_refused")


def test_closed_loop_without_winding_is_still_split_refused():
    """CONTROL: the closed-loop deck, minus the winding token, still aborts on
    the SPLIT vote.  This is what makes
    test_closed_loop_master_runs_under_winding a capability claim rather than a
    deck that would have run anyway.
    """
    _assert_refused("closed_loop_no_winding_refused")
