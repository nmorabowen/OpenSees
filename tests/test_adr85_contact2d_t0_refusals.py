"""ADR-85 T0 -- the 2D contact REFUSAL MATRIX (gate G-T0(c)).

WHAT THIS FILE IS
    The fork's booked worst case is the silent plausible-wrong run (ADR-78 P0:
    converged, balanced, wrong by 2x).  ADR-85 opens the contact command surface
    to `-ndm 2` decks, and every combination it does NOT yet support has to
    fail LOUDLY and BY NAME rather than assemble something that looks fine.
    This file is that contract, one test per row.

    EXPECTED STATE: every refusal test here FAILS before ADR-85 T0-C lands and
    PASSES after.  That is the point -- the tests are the specification T0-C is
    written against.  The one PASSING control (the last test) is the falsifier
    on the falsifiers: it proves the guards refuse the unsupported cases rather
    than refusing everything 2D.

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

    # A well-formed ALL-2D NTS pair.  The kernel and the FE SEGMENT 2D path do
    # not exist until T1b, so the handler must REFUSE and say where the
    # capability lands -- never silently accept a 2D pair and assemble nothing
    # (the ADR-78 P0 failure mode).
    "pair_2d_not_yet_supported": ("ADR-85", "T1b"),

    # Same for -mortar: the 2D interval integrator lands in T3.
    "mortar_2d_not_yet_supported": ("ADR-85", "T3"),

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

elif CASE == "pair_2d_not_yet_supported":
    # A well-formed all-2D NTS pair (ndf == ndm == 2, nps == 2, all nodes 2D).
    # DECLARATION may be accepted -- T0-C opens nps = 2 on 2D nodes -- but
    # ASSEMBLY must be refused by name: the FE SEGMENT 2D path lands in T1b.
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    m = _line_master_2d()
    ops.node(1, 0.0, -1.0e-8)
    ops.fix(1, 1, 0)
    ops.contactSurface(10, "-master", 2, *m)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0)
    _pattern()
    ops.load(1, 0.0, -1.0)
    rc = _static_1_step()

elif CASE == "mortar_2d_not_yet_supported":
    # Same, for the mortar lane: the 2D interval integrator lands in T3.
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    m = _line_master_2d()
    ops.node(1, -0.5, -1.0e-4)
    ops.node(2, 0.5, -1.0e-4)
    ops.fix(1, 1, 0)
    ops.fix(2, 1, 0)
    ops.contactSurface(10, "-master", 2, *m)
    ops.contactSurface(20, "-slave-segments", 2, 1, 2)
    ops.contact(1, 10, 20, "-mortar", "-epsN", EPSN)
    _pattern()
    ops.load(1, 0.0, -1.0)
    ops.load(2, 0.0, -1.0)
    rc = _static_1_step()

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


def test_refuse_pair_2d_not_yet_supported():
    """A valid all-2D NTS pair -> refused by name, pointing at T1b.

    EXPECTED TO FAIL until ADR-85 T0-C.  The failure mode this forecloses is
    the expensive one: T0-C opens the DECLARATION path (nps = 2 on 2D nodes)
    while the FE SEGMENT 2D construction path does not exist until T1b.  If the
    handler simply skipped such a pair, a 2D deck would converge, balance its
    reactions and transmit no contact force -- the ADR-78 P0 incident, again.
    """
    _assert_refused("pair_2d_not_yet_supported")


def test_refuse_mortar_2d_not_yet_supported():
    """A valid all-2D mortar pair -> refused by name, pointing at T3.

    EXPECTED TO FAIL until ADR-85 T0-C.  Same argument as the NTS row; the 2D
    interval-overlap integrator lands in T3.
    """
    _assert_refused("mortar_2d_not_yet_supported")


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
    every 2D contact declaration outright would pass all five falsifiers.  The
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
