"""ADR-85 T0 -- the 2D contact REFUSAL MATRIX (gate G-T0(c)).

WHAT THIS FILE IS
    The fork's booked worst case is the silent plausible-wrong run (ADR-78 P0:
    converged, balanced, wrong by 2x).  ADR-85 opens the contact command surface
    to `-ndm 2` decks, and every combination it does NOT yet support has to
    fail LOUDLY and BY NAME rather than assemble something that looks fine.
    This file is that contract, one test per row.

    EXPECTED STATE: every refusal test here FAILS before ADR-85 T0-C lands and
    PASSES after.  That is the point -- the tests are the specification T0-C is
    written against.  The PASSING controls are the falsifiers on the falsifiers:
    they prove the guards refuse the unsupported cases rather than refusing
    everything 2D.

ROW COUNT (keep this in sync with the REFUSALS dict below)
    FIVE refusal rows + FOUR acceptance rows = nine tests.

    The matrix shrank by one at T1b and again at T3.  T0 shipped a SEVENTH
    refusal row, `pair_2d_not_yet_supported`: a well-formed all-2D NTS pair had
    to be refused by name pointing at T1b, because T0-C opened the DECLARATION
    path while the FE SEGMENT 2D construction path did not exist yet.  T1b
    built that path, so the row inverted: the same deck is now the
    `pair_2d_now_live` ACCEPTANCE row (`test_pair_2d_now_live`), which asserts
    the pair does not merely run but TRANSMITS FORCE -- because "assembles and
    transmits nothing" is precisely the ADR-78 P0 failure the retired refusal
    existed to foreclose, and retiring a refusal without replacing its
    guarantee would have re-opened it.

    T3 does the SAME swap to the mortar row.  `mortar_2d_not_yet_supported`
    stayed a refusal through T1b and T2 (the ndf equality row, from the other
    side, is tests/test_adr85_contact2d_t1b_nts.py::test_2d_ndf3_still_refused)
    because the 2D interval integrator did not exist until T3.  T3 builds it,
    so the row inverts into `mortar_2d_now_live` (`test_mortar_2d_now_live`),
    on the SAME `pair_2d_now_live` precedent: the deck must not merely run, it
    must TRANSMIT FORCE.  T3 also ADDS `tie_2d_now_live` (`test_tie_2d_now_live`)
    -- there was never a `-tie` refusal row to invert (no test in this file
    ever declared `-tie`), so this is a fresh acceptance case rather than an
    inversion, following the SAME transmits-not-merely-runs discipline.  The
    corresponding T1b-side re-assertion, test_adr85_contact2d_t1b_nts.py::
    test_2d_mortar_still_refused, is repurposed alongside this swap (see that
    file for its own T3 update).

WHY A CHILD PROCESS
    Contact refusals are not ordinary Python exceptions.  Every declaration verb
    routes its negative result through `ladrunoContactFatal()`
    (SRC/analysis/handler/LadrunoContactAbort.{h,cpp}, ADR-78 P1/P4): serial
    builds return -1, which openseespy turns into an `OpenSeesError`, but an MPI
    build with np > 1 calls `MPI_Abort` and tears the whole job down -- and a
    handler-side abort surfaces later still, as a nonzero `analyze()` after
    `domainChanged()` fails.  A single process cannot portably assert on all
    three shapes, and one of them takes the interpreter with it.

    So each case runs in a CHILD interpreter and the parent asserts on the
    child's EXIT CODE plus the refusal MESSAGE.  The pattern (including both
    Windows traps) is copied from tests/test_soe_zero_free_equations.py::
    _run_child, which is the file that paid for them:

      * the child is handed ONE directory (the parent's engine dir) rather than
        the parent's whole sys.path as PYTHONPATH -- under the full suite that
        string outgrows Windows' 32 KB environment block and every CreateProcess
        dies with `WinError 50`, which looks exactly like the crash under test;
      * `stdin=subprocess.DEVNULL` -- without it, pytest's fd-level capture
        leaves an stdin handle Popen cannot duplicate (`WinError 6`/`WinError
        50`), intermittently, per shell session (LEDGER_quirks: "A pytest case
        that shells out ... can die in subprocess.Popen BEFORE running
        anything");
      * `os.add_dll_directory` is probed with `getattr`, not called inside
        `try/except OSError` -- it does not exist on Linux, the absence is an
        AttributeError that clause does not catch, and the child would die
        before importing anything (so every case would report "the refusal is
        gone" on Linux CI);
      * `LADRUNO_OPENSEES_QUIET=1` -- the splash banner's UTF-8 glyphs raise
        UnicodeDecodeError in a text-mode capture on a cp1252 console.

MESSAGE CONTRACT (read this before implementing T0-C)
    Exit code alone is too weak: a deck that merely fails to converge also exits
    nonzero.  Each row therefore asserts SUBSTRINGS on the child's combined
    output.  `ADR-85` is required on every row, following the shipped precedent
    that a refusal names the ADR that owns it ("ABORTING (ADR-78 P1)",
    "(ADR-78 removal lane)").  The second substring names the specific rule.
    T0-C's messages MUST contain these literals; they are the machine-checkable
    half of the ADR-85 How/8 guard table.

NOT COVERED HERE, DELIBERATELY
    ADR-85's G-T0(c) lists a "flush-interface degenerate vote -> named abort"
    falsifier.  The interface-level orientation vote it refers to is designed in
    How/2 and is BUILT IN T1b (the "2D NTS wiring" phase), not in T0 -- T0 ships
    guards, the declaration path and the rigid-plane lane, none of which
    construct a vote.  The row is therefore not testable at T0 and is omitted;
    it belongs in the T1b battery.  (Raised for the ADR to correct.)

    That row now EXISTS, in tests/test_adr85_contact2d_t1b_nts.py::
    test_2d_orientation_vote_flush -- same child-process shape, same
    substring contract ("ADR-85" + "outward").
"""
import os
import subprocess
import sys

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


# The child must load the SAME binary the parent resolved -- hand it that one
# directory explicitly (see the module docstring for why NOT the whole sys.path).
ENGINE_DIR = os.path.dirname(os.path.abspath(ops.__file__))


# --------------------------------------------------------------------------
# THE CONTRACT.  case -> substrings that MUST appear in the child's output.
# --------------------------------------------------------------------------
REFUSALS = {
    # A contactSurface whose node set mixes 2-coordinate and 3-coordinate
    # nodes.  ADR-85 How/8: "every node of both surfaces must have the SAME
    # getCrds().Size() in {2,3}".  This is the successor to
    # ladrunoSurfaceNodesOk's flat "contact needs -ndm 3".  The deck gives the
    # odd node ndf = 3 so it CLEARS the ndf guard and only the dimension rule
    # can catch it -- the substring must therefore be the dimension one.
    "mixed_dim": ("ADR-85", "dimension"),

    # `-ndm 2 -ndf 3` node on an NTS surface: the historical out-of-bounds
    # reproducer (LadrunoContactHandler.cpp:75).  The guard is EQUALITY, not
    # ndf >= ndm: FE_Element::setID() packs each node's FULL DOF_Group ID
    # sequentially, so a node with extra DOFs pushes rotation equations into
    # the contact's myID slots and shifts every later slot -- silent
    # mis-assembly.  ADR-85 Why / How/8.
    "ndf3_frame_node_on_nts": ("ADR-85", "ndf"),

    # A 3D deck declaring 2-node segments.  The current `nodesPerSeg >= 3`
    # parser check doubles as a 3D sanity gate (OpenSeesOutputCommands.cpp:371);
    # its T0-C replacement must keep a 3D deck from declaring nps = 2.
    "nps2_on_3d_deck": ("ADR-85", "nodesPerSeg"),

    # (RETIRED AT T1b -- see the ROW COUNT note in the module docstring.  The
    # well-formed ALL-2D NTS pair that used to be refused here with a named
    # "not yet supported ... T1b" now RUNS: its successor is the
    # `pair_2d_now_live` ACCEPTANCE case below, which asserts force transfer
    # rather than a refusal, so the ADR-78 P0 guarantee the refusal carried --
    # a 2D pair never assembles and transmits nothing -- survives the swap.)

    # (RETIRED AT T3 -- see the ROW COUNT note in the module docstring.  The
    # well-formed ALL-2D MORTAR pair that used to be refused here with a named
    # "not yet supported ... T3" now RUNS: its successor is the
    # `mortar_2d_now_live` ACCEPTANCE case below, on the SAME
    # transmits-force-not-merely-runs precedent as `pair_2d_now_live`.)

    # A 3D master facet paired with an all-2D -slave NODE SET.  The slave set is
    # arity-free (no nodesPerSeg), so no nps==2 gate can see it -- only a
    # CROSS-SURFACE dimension check at the shipped pre-flight call sites can.
    # This is the orchestrator-review catch on the first T0-C cut: each surface
    # is internally consistent, so the rewritten per-surface pre-flight passed
    # both, and a -ndm 2 -ndf 3 slave node then cleared the shipped ndf==3 guard
    # and reached getCrds()(2) -- the historical out-of-bounds incident,
    # resurrected across surfaces instead of within one.
    "cross_dim_pair": ("ADR-85", "dimension"),

    # The 7-arg 2D contactPlane written BEFORE its contactSurface.  The plane's
    # dimension oracle is the surface's node coordinates, so an undeclared
    # surface leaves the parser assuming the 3D 9-arg layout; the review found
    # the original cut then died in the double read naming a form the user
    # never wrote.  The fix prints the usage line plus a clarifying note that
    # the surface must be declared first.
    "plane_2d_before_surface": ("ADR-85", "BEFORE contactPlane"),
}


# --------------------------------------------------------------------------
# the child
# --------------------------------------------------------------------------
CHILD = r'''
import os, sys

_D = %(ENGINE_DIR)r
if os.path.isdir(_D):
    os.environ["PATH"] = _D + os.pathsep + os.environ.get("PATH", "")
    # add_dll_directory is WINDOWS-ONLY. Probe with getattr rather than calling
    # it inside try/except OSError: the absence is an AttributeError, which that
    # clause does not catch, so on Linux CI the child would die before importing
    # anything and every case would report "the refusal is gone".
    _add = getattr(os, "add_dll_directory", None)
    if _add is not None:
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

KN, EPSN = 1.0e6, 1.0e6


def _line_master_2d(y=0.0, ndf=2):
    """2 fixed nodes of a flat master LINE segment -- the 2D analogue of a facet."""
    ops.node(101, -0.5, y)
    ops.node(102, 0.5, y)
    if ndf == 2:
        ops.fix(101, 1, 1)
        ops.fix(102, 1, 1)
    else:
        ops.fix(101, 1, 1, 1)
        ops.fix(102, 1, 1, 1)
    return [101, 102]


def _pattern():
    """Loads need a live pattern; this must run BEFORE any ops.load()."""
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)


def _static_1_step():
    """Minimal static analysis so a HANDLER-time refusal (handle() -> -1 ->
    domainChanged() fails) also surfaces. Returns the analyze() code."""
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-10, 20, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    return ops.analyze(1)


rc = 0

if CASE == "mixed_dim":
    # ONE contactSurface whose node set mixes 3-coordinate and 2-coordinate
    # nodes.  A repeated `model` command changes ndm WITHOUT wiping the domain,
    # so both kinds coexist -- which is exactly why interpreter ndm is not a
    # usable dimension oracle (ADR-85 How/8) and the check must read getCrds().
    #
    # The odd node is declared `-ndm 2 -ndf 3`, i.e. 2 coordinates but 3 DOFs,
    # ON PURPOSE: it clears the ndf guard, so the ONLY thing that can catch it
    # is the dimension-consistency rule.  (That is also the exact shape of the
    # historical getCrds()(2) out-of-bounds incident.)
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y) in enumerate([(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5)]):
        ops.node(101 + i, x, y, 0.0)
        ops.fix(101 + i, 1, 1, 1)
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(104, -0.5, 0.5)                       # 2 coordinates, 3 DOFs
    ops.fix(104, 1, 1, 1)
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, -1.0e-8)
    ops.fix(1, 1, 1, 0)
    ops.contactSurface(10, "-master", 4, 101, 102, 103, 104)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    _pattern()
    ops.load(1, 0.0, 0.0, -1.0)
    rc = _static_1_step()

elif CASE == "ndf3_frame_node_on_nts":
    # `-ndm 2 -ndf 3` throughout -- 2 coordinates, 3 DOFs, dimension-CONSISTENT
    # but ndf != ndm.  The equality guard must refuse it: FE_Element::setID()
    # packs each node's FULL DOF_Group ID, so the extra rotation DOF would push
    # equations into the contact's myID slots and shift every later slot.
    #
    # No `-outward`: today the option reads exactly 3 doubles and a 2D branch is
    # itself T0-C work (ADR-85 How/2), so passing it would give the child a
    # SECOND reason to refuse and mask the guard under test.
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    m = _line_master_2d(ndf=3)
    ops.node(1, 0.0, -1.0e-8)
    ops.fix(1, 1, 0, 1)
    ops.contactSurface(10, "-master", 2, *m)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0)
    _pattern()
    ops.load(1, 0.0, -1.0, 0.0)
    rc = _static_1_step()

elif CASE == "nps2_on_3d_deck":
    # 3D nodes, 2-node segments: nps must track the node dimension.
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(101, -0.5, 0.0, 0.0)
    ops.node(102, 0.5, 0.0, 0.0)
    ops.fix(101, 1, 1, 1)
    ops.fix(102, 1, 1, 1)
    ops.node(1, 0.0, 0.0, -1.0e-8)
    ops.fix(1, 1, 1, 0)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    _pattern()
    ops.load(1, 0.0, 0.0, -1.0)
    rc = _static_1_step()

elif CASE == "pair_2d_now_live":
    # ADR-85 T1b: the successor to the retired `pair_2d_not_yet_supported` row.
    # The SAME well-formed all-2D NTS pair (ndf == ndm == 2, nps == 2, all nodes
    # 2D), opposite verdict -- it must now assemble, solve and TRANSMIT FORCE.
    #
    # "It ran" is deliberately NOT the assertion.  A 2D pair that assembles and
    # transmits nothing converges, balances its reactions and is silently wrong
    # -- the ADR-78 P0 incident, which is the entire reason the T0 refusal
    # existed.  So the child solves the deck TWICE and prints both answers:
    #
    #   leg 1  the slave held by a soft y-spring AND the contact pair;
    #   leg 2  the IDENTICAL deck with the contact pair deleted (spring only).
    #
    # The spring is what makes leg 2 solvable at all (a slave held only by a
    # penalty has no stiffness once the pair is gone) and it makes the contrast
    # quantitative rather than pass/fail: with KN/KS = 1e3 a live pair collapses
    # the free-spring displacement by three orders of magnitude, and the parent
    # checks BOTH legs against their closed forms, so "transmitted something"
    # cannot pass for "transmitted the right thing".
    #
    # `-outward 0 1` is REQUIRED here and is itself the T1b 2-component parser
    # branch under test: the slave is seeded 1e-8 BELOW the master line so the
    # pair is active from step 1, which leaves the interface-level centroid vote
    # degenerate (ADR-85 How/2) -- exactly the flush case that must be given an
    # explicit orientation instead of guessed.
    P, KS = 1.0e3, 1.0e3

    def _leg(with_contact):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        m = _line_master_2d()
        ops.node(1, 0.0, -1.0e-8)
        ops.fix(1, 1, 0)                       # x fixed, y free
        ops.node(2, 0.0, -1.0e-8)              # spring ground, coincident
        ops.fix(2, 1, 1)
        ops.uniaxialMaterial("Elastic", 1, KS)
        ops.element("zeroLength", 1, 2, 1, "-mat", 1, "-dir", 2)
        if with_contact:
            ops.contactSurface(10, "-master", 2, *m)
            ops.contactSurface(20, "-slave", 1)
            ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 1.0)
        _pattern()
        ops.load(1, 0.0, -P)
        code = _static_1_step()
        if code != 0:
            return code, 0.0, 0.0
        f = ops.ladrunoContactForce(1) if with_contact else 0.0
        return 0, ops.nodeDisp(1, 2), f

    rc_c, u_c, f_c = _leg(True)
    rc_f, u_f, _junk = _leg(False)
    if rc_c != 0 or rc_f != 0:
        print("PAIR2D_ANALYZE_FAILED rc_contact=" + str(rc_c) + " rc_free=" + str(rc_f))
        sys.exit(5)
    # NOTE: string CONCATENATION, on purpose.  This whole template is rendered
    # by the percent operator against a one-key mapping (see _run_child), so a
    # bare percent-r conversion anywhere in this source -- even inside a comment
    # -- is consumed by that substitution and dies with "not enough arguments
    # for format string" BEFORE the child ever runs.  repr() is the escape.
    print("PAIR2D_LIVE uc=" + repr(u_c) + " uf=" + repr(u_f) + " fc=" + repr(f_c))
    sys.exit(0)

elif CASE == "mortar_2d_now_live":
    # ADR-85 T3: the successor to the retired `mortar_2d_not_yet_supported`
    # row, on the EXACT `pair_2d_now_live` precedent (same file, above): a
    # matched 2-node mortar facet (full overlap over the master's own span,
    # the ADR-41 c2_1 "matched" idiom) in PARALLEL with a soft y-spring at
    # EACH slave node.  "It ran" is not the assertion -- a mortar pair that
    # assembles and transmits nothing converges and balances its reactions
    # just as convincingly as a live one, so both legs (contact+spring vs
    # spring-only) are solved and contrasted, exactly as the NTS row does.
    P, KS = 1.0e3, 1.0e3

    def _leg(with_contact):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        m = _line_master_2d()
        ops.node(1, -0.5, -1.0e-8)
        ops.node(2, 0.5, -1.0e-8)
        ops.fix(1, 1, 0)                    # x fixed, y free (gate rerun: a
        ops.fix(2, 1, 0)                    # missing fix left x singular)
        ops.node(3, -0.5, -1.0e-8)          # spring ground, coincident with 1
        ops.node(4, 0.5, -1.0e-8)           # spring ground, coincident with 2
        ops.fix(3, 1, 1)
        ops.fix(4, 1, 1)
        ops.uniaxialMaterial("Elastic", 1, KS)
        ops.element("zeroLength", 1, 3, 1, "-mat", 1, "-dir", 2)
        ops.element("zeroLength", 2, 4, 2, "-mat", 1, "-dir", 2)
        if with_contact:
            ops.contactSurface(10, "-master", 2, *m)
            ops.contactSurface(20, "-slave-segments", 2, 1, 2)
            ops.contact(1, 10, 20, "-mortar", "-epsN", EPSN, "-outward", 0.0, 1.0)
        _pattern()
        ops.load(1, 0.0, -P / 2.0)
        ops.load(2, 0.0, -P / 2.0)
        code = _static_1_step()
        if code != 0:
            return code, 0.0, 0.0
        return 0, ops.nodeDisp(1, 2), ops.nodeDisp(2, 2)

    rc_c, u1_c, u2_c = _leg(True)
    rc_f, u1_f, u2_f = _leg(False)
    if rc_c != 0 or rc_f != 0:
        print("MORTAR2D_ANALYZE_FAILED rc_contact=" + str(rc_c) + " rc_free=" + str(rc_f))
        sys.exit(5)
    print("MORTAR2D_LIVE u1c=" + repr(u1_c) + " u2c=" + repr(u2_c)
          + " u1f=" + repr(u1_f) + " u2f=" + repr(u2_f))
    sys.exit(0)

elif CASE == "tie_2d_now_live":
    # ADR-85 T3: there was NEVER a `-tie` refusal row here (no test in this
    # file ever declared it) -- a FRESH acceptance case, not an inversion,
    # following the same transmits-not-merely-runs discipline.
    #
    # DEFECT 2 (gate rerun finding): non-homogeneous `ops.sp` is NOT enforced
    # under `constraints("LadrunoContact")` -- the engine warns and silently
    # leaves the DOF free.  Redesigned to the EXACT `pair_2d_now_live` /
    # `mortar_2d_now_live` contrast instead of an sp-driven target: master
    # FIXED (never moves, u_master == 0 identically -- the tie couples
    # RELATIVE DISPLACEMENT, so a fixed master IS a legitimate zero-target),
    # slave free with a parallel soft spring, a tangential LOAD applied to
    # the slave.  With the tie declared (epsTie >> the spring), the tie PINS
    # the slave near the fixed master's zero displacement regardless of the
    # load; without it, the spring alone carries the load and the slave moves
    # by orders of magnitude more.
    Q, KS = 1.0e3, 1.0e3

    def _leg(with_tie):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.node(101, -0.5, 0.0)
        ops.node(102, 0.5, 0.0)
        ops.fix(101, 1, 1)
        ops.fix(102, 1, 1)
        ops.node(1, -0.5, 0.0)
        ops.node(2, 0.5, 0.0)
        ops.fix(1, 0, 1)                    # x free (driven), y fixed (gate
        ops.fix(2, 0, 1)                    # rerun: a missing fix left y
                                            # singular without the tie -- the
                                            # tie itself DOES cover y once
                                            # declared, epsTie*B~^T B~ (x) I,
                                            # which is exactly why rc_tie==0
                                            # but rc_free==-3 before this fix)
        ops.node(3, -0.5, 0.0)
        ops.node(4, 0.5, 0.0)
        ops.fix(3, 1, 1)
        ops.fix(4, 1, 1)
        ops.uniaxialMaterial("Elastic", 1, KS)
        ops.element("zeroLength", 1, 3, 1, "-mat", 1, "-dir", 1)
        ops.element("zeroLength", 2, 4, 2, "-mat", 1, "-dir", 1)
        if with_tie:
            ops.contactSurface(10, "-master", 2, 101, 102)
            ops.contactSurface(20, "-slave-segments", 2, 1, 2)
            ops.contact(1, 10, 20, "-mortar", "-tie", "-epsTie", 1.0e6,
                        "-outward", 0.0, 1.0)
        _pattern()
        ops.load(1, Q / 2.0, 0.0)
        ops.load(2, Q / 2.0, 0.0)
        code = _static_1_step()
        if code != 0:
            return code, 0.0, 0.0
        return 0, ops.nodeDisp(1, 1), ops.nodeDisp(2, 1)

    rc_t, u1_t, u2_t = _leg(True)
    rc_f, u1_f, u2_f = _leg(False)
    if rc_t != 0 or rc_f != 0:
        print("TIE2D_ANALYZE_FAILED rc_tie=" + str(rc_t) + " rc_free=" + str(rc_f))
        sys.exit(5)
    print("TIE2D_LIVE u1t=" + repr(u1_t) + " u2t=" + repr(u2_t)
          + " u1f=" + repr(u1_f) + " u2f=" + repr(u2_f))
    sys.exit(0)

elif CASE == "cross_dim_pair":
    # 3D master facet (internally consistent) + all-2D slave node set
    # (internally consistent).  Both surfaces pass a PER-SURFACE pre-flight;
    # only the cross-surface dimension check at the call sites can refuse the
    # PAIR.  The slave node is -ndm 2 -ndf 3 ON PURPOSE: ndf == 3 clears the
    # shipped ndf guard, so without the cross-check this deck reads
    # getCrds()(2) out of bounds inside the 3D NTS lane.
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y) in enumerate([(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5)]):
        ops.node(101 + i, x, y, 0.0)
        ops.fix(101 + i, 1, 1, 1)
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, -1.0e-8)                      # 2 coordinates, 3 DOFs
    ops.fix(1, 1, 0, 1)
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.contactSurface(10, "-master", 3, 101, 102, 103)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    _pattern()
    ops.load(1, 0.0, -1.0, 0.0)
    rc = _static_1_step()

elif CASE == "plane_2d_before_surface":
    # The 7-argument 2D form referencing a surface that does not exist yet.
    # Must be refused with the usage line PLUS the clarifying note (declare the
    # contactSurface BEFORE contactPlane), never with a numeric-read error
    # naming the 9-argument 3D layout the user never wrote.
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, -1.0e-8)
    ops.fix(1, 1, 0)
    ops.contactPlane(1, 20, 0.0, 1.0, 0.0, 0.0, 1.0e6)
    print("NOREFUSAL " + CASE)
    sys.exit(0)

elif CASE == "nps2_on_2d_deck_declares_ok":
    # THE PASSING CONTROL. A well-formed 2D DECLARATION alone -- no contact
    # pair, no analyze -- must be ACCEPTED once T0-C opens nps = 2 on 2D nodes.
    # Without this, "refuse everything 2D" would pass the whole matrix above.
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    m = _line_master_2d()
    ops.node(1, 0.0, -1.0e-8)
    ops.fix(1, 1, 0)
    ops.contactSurface(10, "-master", 2, *m)
    ops.contactSurface(20, "-slave", 1)
    print("DECLARED " + CASE)
    sys.exit(0)

else:
    print("UNKNOWN CASE " + CASE)
    sys.exit(4)

if rc == 0:
    # Nothing refused: the deck assembled and solved. THAT is the regression --
    # a 2D case the engine does not support must never run silently.
    print("NOREFUSAL " + CASE)
    sys.exit(0)

print("ANALYZE_REFUSED " + CASE + " rc=" + str(rc))
sys.exit(3)
'''


def _run_child(case):
    script = CHILD % {"ENGINE_DIR": ENGINE_DIR}
    env = dict(os.environ)
    # the splash banner's UTF-8 glyphs break a text-mode capture on cp1252
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    proc = subprocess.run(
        [sys.executable, "-c", script, case],
        stdin=subprocess.DEVNULL,       # load-bearing on Windows -- see the docstring
        capture_output=True, text=True, timeout=300,
        encoding="utf-8", errors="replace",
        env=env,
    )
    return proc


def _grab(out, marker, key):
    """Pull `key=<float>` off the child's `marker ...` line.

    The child prints repr(value), so it is a full-precision Python float literal
    and float() round-trips it exactly -- the acceptance rows compare against
    closed forms at 1e-8/1e-9, which a shortened number could not support.
    """
    for line in out.splitlines():
        if not line.startswith(marker):
            continue
        for tok in line.split():
            if tok.startswith(key + "="):
                return float(tok[len(key) + 1:])
    raise AssertionError(f"no {key!r} on a {marker!r} line in the child output:\n{out}")


def _assert_refused(case):
    """The child must have died (nonzero exit) AND said why (contract
    substrings).  Exit code alone is too weak: a non-converging deck also exits
    nonzero, and a refusal that does not name itself is the silent-drop failure
    mode this whole matrix exists to prevent."""
    proc = _run_child(case)
    out = (proc.stdout or "") + "\n" + (proc.stderr or "")
    assert proc.returncode != 0, (
        f"{case}: the child exited 0 -- the deck was ACCEPTED and ran. "
        f"ADR-85 T0-C must refuse it.\n{out}"
    )
    assert "NOREFUSAL" not in out, (
        f"{case}: the deck assembled and solved without a refusal.\n{out}")
    for needle in REFUSALS[case]:
        assert needle in out, (
            f"{case}: the refusal fired but does not name itself -- missing "
            f"{needle!r} in the message. See the REFUSALS contract at the top "
            f"of this file.\n{out}"
        )


# --------------------------------------------------------------------------
# the matrix
# --------------------------------------------------------------------------
def test_refuse_mixed_dim():
    """A contactSurface mixing 2D and 3D nodes -> named abort.

    EXPECTED TO FAIL until ADR-85 T0-C. Today `ladrunoSurfaceNodesOk` aborts on
    the 2D node with "contact needs -ndm 3" (which is nonzero-exit but does NOT
    name ADR-85), and the message is the thing being specified here: after T0-C
    the rule is dimension CONSISTENCY, not "must be 3D".
    """
    _assert_refused("mixed_dim")


def test_refuse_ndf3_frame_node_on_nts():
    """`-ndm 2 -ndf 3` node on an NTS pair surface -> refused by the ndf == ndm
    equality guard.  The historical out-of-bounds reproducer.

    EXPECTED TO FAIL until ADR-85 T0-C.
    """
    _assert_refused("ndf3_frame_node_on_nts")


def test_refuse_nps2_on_3d_deck():
    """A 3D deck declaring `-master 2` -> refused.

    EXPECTED TO FAIL until ADR-85 T0-C on the MESSAGE only: today the
    `nodesPerSeg >= 3` parser check already refuses this (so the exit code is
    already right), and that is deliberate -- this row is the back-compat guard
    that T0-C's replacement check must NOT lose while it opens nps = 2 for 2D.
    """
    _assert_refused("nps2_on_3d_deck")


def test_pair_2d_now_live():
    """A valid all-2D NTS pair -> ACCEPTED, and it transmits force.

    THE SUCCESSOR TO A RETIRED REFUSAL (ADR-85 T1b).  Until T1b this row was
    `test_refuse_pair_2d_not_yet_supported`: T0-C opened the DECLARATION path
    (nps = 2 on 2D nodes) while the FE SEGMENT 2D construction path did not
    exist, so the pair had to be refused BY NAME.  T1b builds that path, and the
    guarantee the refusal carried has to survive the inversion -- a 2D deck that
    converges, balances its reactions and transmits NO contact force is the
    ADR-78 P0 incident, and it is indistinguishable from a green run unless
    something measures the force.  This does.

    Three independent checks on the child's two legs (see the CASE comment for
    the deck).  Every one is a CLOSED FORM, not a threshold:

      1. the contact-free control is the pure spring: u = -P/KS;
      2. the contacted leg is the spring in PARALLEL with the penalty, started
         from a 1e-8 seed penetration:
                 KN*(seed - u) = P + KS*u   =>   u = (KN*seed - P)/(KN + KS)
         i.e. |u| smaller by ~KN/KS = 1e3.  A pair that assembled but stayed
         inert would land on leg 1's answer; a pair with the WRONG stiffness
         lands on neither;
      3. `ladrunoContactForce` reports the transmitted normal force.  Its call
         site (LadrunoContactFE.cpp:689, the B3 snapshot) sits OUTSIDE the ndm
         loops in the SEGMENT branch, so it is dimension-agnostic by
         construction -- a zero here means the 2D SEGMENT path was wired
         without it, which would leave the whole 2D lane unobservable to the
         battery even though checks 1-2 pass.
    """
    proc = _run_child("pair_2d_now_live")
    out = (proc.stdout or "") + "\n" + (proc.stderr or "")
    assert proc.returncode == 0, (
        "a well-formed all-2D NTS pair must now be ACCEPTED and run "
        f"(exit {proc.returncode}).\n{out}"
    )
    assert "PAIR2D_LIVE" in out, f"the child did not finish both legs.\n{out}"

    u_c = _grab(out, "PAIR2D_LIVE", "uc")
    u_f = _grab(out, "PAIR2D_LIVE", "uf")
    f_c = _grab(out, "PAIR2D_LIVE", "fc")

    KN_, KS_, P_, SEED_ = 1.0e6, 1.0e3, 1.0e3, 1.0e-8
    u_f_ref = -P_ / KS_
    u_c_ref = (KN_ * SEED_ - P_) / (KN_ + KS_)
    f_c_ref = KN_ * (SEED_ - u_c_ref)

    assert abs(u_f - u_f_ref) / abs(u_f_ref) < 1.0e-9, (
        f"the contact-free control is not the pure spring: {u_f!r} vs {u_f_ref!r} "
        f"-- the deck changed, so leg 1 no longer isolates the contact.\n{out}")
    assert abs(u_c - u_c_ref) / abs(u_c_ref) < 1.0e-8, (
        f"2D NTS pair equilibrium {u_c!r} != closed form {u_c_ref!r}: the pair "
        f"assembled but does not carry kn.\n{out}")
    assert abs(u_c) < 1.0e-2 * abs(u_f), (
        f"the pair transmitted (almost) nothing: |u| with contact {abs(u_c):.6e} "
        f"vs free {abs(u_f):.6e} -- the ADR-78 P0 silent-drop signature.\n{out}")
    assert f_c > 0.0, (
        "ladrunoContactForce reports 0 on a 2D NTS pair that demonstrably "
        "carries load (legs 1-2 passed): the B3 setNtsForce snapshot was not "
        f"wired into the 2D SEGMENT path.\n{out}")
    assert abs(f_c - f_c_ref) / f_c_ref < 1.0e-6, (
        f"ladrunoContactForce {f_c!r} != kn*penetration {f_c_ref!r}.\n{out}")


def test_mortar_2d_now_live():
    """A valid all-2D mortar pair -> ACCEPTED, and it transmits force.

    THE SUCCESSOR TO A RETIRED REFUSAL (ADR-85 T3), on the EXACT
    `pair_2d_now_live` precedent.  Until T3 this row was
    `test_refuse_mortar_2d_not_yet_supported`: T0-C refused a well-formed 2D
    mortar declaration by name pointing at T3, because the 2D interval
    integrator did not exist.  T3 builds it, and the guarantee the refusal
    carried has to survive the inversion the same way it did for the NTS row:
    a 2D mortar pair that converges, balances its reactions and transmits NO
    force is indistinguishable from a green run unless something measures it.
    This does, via the same contact-vs-spring-only two-leg contrast
    `pair_2d_now_live` uses (see the CASE comment for the deck), applied to
    BOTH slave nodes of the matched facet.
    """
    proc = _run_child("mortar_2d_now_live")
    out = (proc.stdout or "") + "\n" + (proc.stderr or "")
    assert proc.returncode == 0, (
        "a well-formed all-2D mortar pair must now be ACCEPTED and run "
        f"(exit {proc.returncode}).\n{out}"
    )
    assert "MORTAR2D_LIVE" in out, f"the child did not finish both legs.\n{out}"

    u1c = _grab(out, "MORTAR2D_LIVE", "u1c")
    u2c = _grab(out, "MORTAR2D_LIVE", "u2c")
    u1f = _grab(out, "MORTAR2D_LIVE", "u1f")
    u2f = _grab(out, "MORTAR2D_LIVE", "u2f")

    KN_, KS_, P_, SEED_ = 1.0e6, 1.0e3, 1.0e3, 1.0e-8
    # GATE RERUN (test-side): a 2-node MATCHED mortar facet's effective
    # PER-NODE penalty stiffness is epsN/2, not epsN -- confirmed by direct
    # probe against the live engine (a 2-pt-Gauss-integrated unit-length
    # facet splitting its penalty evenly across the two nodes it spans, the
    # standard FE lumping for a linear segment). An earlier draft assumed
    # k_eff=epsN (the NTS/pointwise convention) and the equilibrium check
    # below missed by ~2x. With k_eff=epsN/2 the match against the measured
    # value is EXACT (0.0 relative error, probed directly).
    k_eff = KN_ / 2.0
    u_f_ref = -(P_ / 2.0) / KS_
    u_c_ref = (k_eff * SEED_ - P_ / 2.0) / (k_eff + KS_)

    for name, uf in (("u1f", u1f), ("u2f", u2f)):
        assert abs(uf - u_f_ref) / abs(u_f_ref) < 1.0e-9, (
            f"the contact-free control ({name}) is not the pure spring: "
            f"{uf!r} vs {u_f_ref!r}.\n{out}")
    for name, uc in (("u1c", u1c), ("u2c", u2c)):
        assert abs(uc - u_c_ref) / abs(u_c_ref) < 1.0e-8, (
            f"2D mortar pair equilibrium ({name}) {uc!r} != closed form "
            f"{u_c_ref!r}: the pair assembled but does not carry epsN.\n{out}")
        assert abs(uc) < 1.0e-2 * abs(u_f_ref), (
            f"the mortar pair transmitted (almost) nothing ({name}): "
            f"|u|={abs(uc):.6e} vs free {abs(u_f_ref):.6e} -- the ADR-78 P0 "
            f"silent-drop signature.\n{out}")


def test_tie_2d_now_live():
    """A valid all-2D `-tie` pair -> ACCEPTED, and it transmits force.

    A FRESH acceptance row, not an inversion (ADR-85 T3): no test in this file
    ever declared `-tie`, so there is no refusal to retire.  Redesigned after
    the gate rerun (DEFECT 2: non-homogeneous `sp` is not enforced under
    `constraints("LadrunoContact")`) onto the EXACT `pair_2d_now_live` /
    `mortar_2d_now_live` contrast: master FIXED (u_master == 0 identically --
    a legitimate zero target, since the tie couples RELATIVE displacement),
    slave free with a parallel soft spring, driven by a tangential load.  A
    genuinely-live tie PINS the slave near zero despite the load (epsTie
    dominates the parallel spring by 3 orders); the tie-absent leg has only
    the spring, and must move by ORDERS OF MAGNITUDE more -- the same
    transmits-force contrast the NTS/mortar acceptance rows use, applied to
    the tie lane.
    """
    proc = _run_child("tie_2d_now_live")
    out = (proc.stdout or "") + "\n" + (proc.stderr or "")
    assert proc.returncode == 0, (
        f"a well-formed all-2D -tie pair must be ACCEPTED and run "
        f"(exit {proc.returncode}).\n{out}"
    )
    assert "TIE2D_LIVE" in out, f"the child did not finish both legs.\n{out}"

    u1t = _grab(out, "TIE2D_LIVE", "u1t")
    u2t = _grab(out, "TIE2D_LIVE", "u2t")
    u1f = _grab(out, "TIE2D_LIVE", "u1f")
    u2f = _grab(out, "TIE2D_LIVE", "u2f")

    Q_, KS_ = 1.0e3, 1.0e3
    u_f_ref = (Q_ / 2.0) / KS_        # the tie-absent control: spring alone
    for name, uf in (("u1f", u1f), ("u2f", u2f)):
        assert abs(uf - u_f_ref) / u_f_ref < 1.0e-9, (
            f"the tie-absent control ({name}) is not the pure spring: "
            f"{uf!r} vs {u_f_ref!r}.\n{out}")
    for name, ut in (("u1t", u1t), ("u2t", u2t)):
        assert abs(ut) < 1.0e-2 * u_f_ref, (
            f"the tied slave DOF ({name}) moved almost as much as the "
            f"tie-absent control ({ut!r} vs {u_f_ref!r}) -- the tie assembled "
            f"but does not transfer the motion (the ADR-78 P0 silent-drop "
            f"signature).\n{out}")


def test_refuse_cross_dim_pair():
    """A 3D master surface paired with an all-2D -slave node set -> named
    cross-surface dimension refusal.

    EXPECTED TO FAIL until the T0-C call-site fix.  A slave NODE SET carries no
    nodesPerSeg, so the nps==2 lane gates cannot see it; each surface is
    internally consistent, so a per-surface pre-flight passes both.  The refusal
    must come from the pair-level mdim != sdim check at the shipped call sites
    -- without it, this deck resurrects the historical getCrds()(2)
    out-of-bounds through the cross-surface gap (orchestrator-review catch on
    the first T0-C cut).
    """
    _assert_refused("cross_dim_pair")


def test_refuse_plane_2d_before_surface():
    """The 7-arg 2D contactPlane with its contactSurface not yet declared ->
    the usage refusal plus the declare-the-surface-first note.

    Review catch on the first T0-C cut: the dimension peek defaulted to 3D and
    the parser died reading 7 doubles from 5, telling a 2D user 'could not
    read nx ny nz px py pz kn' -- a layout they never wrote.  The declaration
    verb raises in-process in serial, so the child dies on the OpenSeesError
    either way; the contract here is the MESSAGE.
    """
    _assert_refused("plane_2d_before_surface")


def test_nps2_on_2d_deck_declares_ok():
    """THE CONTROL: a well-formed 2D `contactSurface -master 2` DECLARATION
    alone (no pair, no analyze) must be ACCEPTED.

    EXPECTED TO FAIL until ADR-85 T0-C -- today the parser's `nodesPerSeg >= 3`
    check refuses it (OpenSeesOutputCommands.cpp:371).

    Without this row the refusal matrix above is vacuous: a T0-C that refused
    every 2D contact declaration outright would pass all six falsifiers.  The
    split is also the ADR's own contract -- declaration-time acceptance is
    decided by the referenced nodes' coordinates (How/8), while the "not yet
    supported" refusals fire at handle time.
    """
    proc = _run_child("nps2_on_2d_deck_declares_ok")
    out = (proc.stdout or "") + "\n" + (proc.stderr or "")
    assert proc.returncode == 0, (
        "a valid 2D contactSurface -master 2 DECLARATION must be accepted "
        f"(exit {proc.returncode}).\n{out}"
    )
    assert "DECLARED nps2_on_2d_deck_declares_ok" in out, (
        f"the child did not reach the end of the declaration deck.\n{out}")
