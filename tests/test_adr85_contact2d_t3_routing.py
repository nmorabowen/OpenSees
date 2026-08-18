"""ADR-85 T3 -- routing (NTS ≻ mortar stand-down) + the T3 scope fences (gate G-T3(d)).

WHAT THIS FILE IS
    Written against `_adr85_t3_design.md` §2b.4 (routing) and §2b.6/§2b.8
    (fences), and the ADR's How/3 "Routing" paragraph, for C++ that does not
    exist yet (see test_adr85_contact2d_t3_mortar.py's module docstring for
    the general disclaimer).

    A 2D mortar pair that is geometrically near-PERPENDICULAR can never
    integrate a meaningful overlap (two straight lines at a large angle meet
    at, at most, a point -- the exact 2D analogue of the shipped 3D edge-edge
    ADR's "exactly-perpendicular headline case where the face-mortar clip is
    genuinely zero", test_adr57_edge_edge_1.py). The 2D lane has no edge-edge
    tier to catch that geometry, so ADR-85 adopts the ADR-57 E7 ownership
    protocol minus the edge tier: NTS >= mortar, with a NEW named latch,
    `WARN_2D_PERP_NO_NTS`, firing once whenever a perpendicular-refused mortar
    pair is NOT covered by an NTS contact declared on the same nodes.

COVERAGE MAP -> ADR-85 G-T3(d)
    test_2d_perp_pair_nts_codeclared_transmits_no_warn   NTS + mortar declared
        on the same interface nodes: NTS owns the transfer, no double-count,
        WARN_2D_PERP_NO_NTS does NOT fire.
    test_2d_perp_pair_mortar_only_warns_once             mortar-only: the pair
        is inert AND the latch fires -- exactly once, despite two analyze()
        calls (the warnOnce contract, ADR-85 t1b_nts.py::test_2d_geomtan_noop's
        own latch-per-process precedent).
    test_2d_soft_on_2d_mortar_refused             -soft on a 2D mortar deck
    test_2d_edgeedge_on_2d_mortar_refused         -edgeedge on a 2D mortar deck
    test_2d_thickness_without_mortar_refused      -thickness needs -mortar
    test_2d_thickness_nonpositive_refused         -thickness must be > 0
    test_2d_thickness_on_3d_mortar_refused        -thickness is a 2D-only knob

MESSAGE CONTRACT (the exact substrings asserted below -- report these to the
C++ implementer for reconciliation, per the phase brief)
    WARN_2D_PERP_NO_NTS  ..  "WARN" + "NTS" + the mortar contact's own tag ("1"
                             -- mortar is declared with a FIXED contact tag of
                             1 in `_build()`, whether or not NTS is co-declared)
    -soft on 2D mortar   ..  "SOFT=2" + "ADR-85"
    -edgeedge on 2D       ..  "no 2D analogue" + "scope fence"
    -thickness (generic)  ..  "-thickness"
    -thickness on 3D      ..  "-thickness" + "2D plane-model option"

WHY THE "PERPENDICULAR PAIR" DECK IS SHAPED THE WAY IT IS
    Master: a fixed horizontal segment spanning x in [-1,1]. Mortar slave: a
    2-node VERTICAL segment (node 1 at x=0, seeded a hair below the master;
    node 2 one unit further down) -- exactly perpendicular to the master
    (t_s . t_m == 0 by construction: both slave endpoints share x=0, so the
    closed-form projection ξ_M(point) = (point.x - X0.x)/L depends ONLY on
    the point's x-coordinate here, and both endpoints map to the SAME ξ,
    giving a genuinely empty/degenerate overlap -- there is no finite-width
    "tilted-but-still-perpendicular-ish" case to construct in pure 2D, since
    two non-parallel straight lines meet at one point). Node 1 alone is also
    declared on an NTS point-set surface, seeded PENETRATING the master, so
    NTS treats it as an ordinary active pair independent of the mortar
    segment's orientation (NTS never looks at slave-segment tangents at all).
    Node 2 is fully fixed and carries no load -- it exists only to give the
    mortar slave surface its second (perpendicular-defining) node.

    A soft y-spring is present on node 1 in EVERY leg (the `pair_2d_now_live`
    idiom, test_adr85_contact2d_t0_refusals.py) so the mortar-ONLY leg is
    still well-posed even though its contact pair is inert.
"""
import os
import subprocess
import sys

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

ENGINE_DIR = os.path.dirname(os.path.abspath(ops.__file__))

# --------------------------------------------------------------------------
# the routing legs -- child process (stderr capture is the point)
# --------------------------------------------------------------------------
CHILD_ROUTE = r'''
import os, sys

_D = %(ENGINE_DIR)r
if os.path.isdir(_D):
    os.environ["PATH"] = _D + os.pathsep + os.environ.get("PATH", "")
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
KN, KS, EPSN, P, SEED = 1.0e6, 1.0e0, 1.0e6, 1.0e0, 1.0e-4


def _build(with_nts, with_mortar):
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -1.0, 0.0)
    ops.node(102, 1.0, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.contactSurface(10, "-master", 2, 101, 102)

    ops.node(1, 0.0, -SEED)
    ops.node(2, 0.0, -SEED - 1.0)
    ops.fix(1, 1, 0)                     # x fixed, y free
    ops.fix(2, 1, 1)                     # fully fixed -- carries no load

    ops.node(901, 0.0, -SEED)
    ops.fix(901, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, KS)
    ops.element("zeroLength", 1, 901, 1, "-mat", 1, "-dir", 2)

    # ORCHESTRATOR UPDATE (post-review): the implemented ownership rule is
    # CONJUNCTIVE -- an NTS contact suppresses the mortar pair's ALIGN warn
    # iff BOTH hold: (1) the NTS contact's MASTER surface tag == the mortar
    # contact's MASTER surface tag (both 10, here), AND (2) the UNION of
    # that NTS contact's slave NODE tags contains BOTH node tags of the
    # specific mortar SLAVE segment that refuses ALIGN (nodes 1 AND 2, not
    # just the touching node 1). Mortar is declared with a FIXED contact tag
    # (1) and slave surface tag (20) regardless of whether NTS is also
    # present, so the WARN's own contact-tag naming stays consistent between
    # the "mortar_only" and "both_codeclared" cases.
    if with_mortar:
        ops.contactSurface(20, "-slave-segments", 2, 1, 2)
        ops.contact(1, 10, 20, "-mortar", "-epsN", EPSN, "-outward", 0.0, 1.0)
    if with_nts:
        ops.contactSurface(21, "-slave", 1, 2)
        ops.contact(2, 10, 21, KN, 0.0, 0.0, "-outward", 0.0, 1.0)

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


if CASE == "nts_only":
    _build(True, False)
    rc = ops.analyze(1)
    assert rc == 0, "nts_only did not converge rc=" + str(rc)
    print("ROUTE_RESULT nts_only u=" + repr(ops.nodeDisp(1, 2)))
elif CASE == "mortar_only":
    _build(False, True)
    rc1 = ops.analyze(1)
    assert rc1 == 0, "mortar_only step1 did not converge rc=" + str(rc1)
    # SECOND analyze() call, at a HELD load (LoadControl 0.0 -- the
    # analyze_augmented idiom): a plain LoadControl(1.0) here would apply a
    # SECOND increment of P, doubling the total load and making "u" reflect
    # 2x the intended equilibrium (gate rerun finding, test-side -- an
    # earlier draft called analyze(1) twice under LoadControl(1.0) and the
    # resulting "u=-2.0" looked like extra transmitted force but was really
    # just P applied twice). The point of the second call is only to
    # re-exercise the warnOnce latch at the SAME converged state.
    ops.integrator("LoadControl", 0.0)
    rc2 = ops.analyze(1)
    assert rc2 == 0, "mortar_only step2 did not converge rc=" + str(rc2)
    print("ROUTE_RESULT mortar_only u=" + repr(ops.nodeDisp(1, 2)))
elif CASE == "both_codeclared":
    _build(True, True)
    rc = ops.analyze(1)
    assert rc == 0, "both_codeclared did not converge rc=" + str(rc)
    print("ROUTE_RESULT both_codeclared u=" + repr(ops.nodeDisp(1, 2)))
elif CASE == "tie_align_always_warns":
    # NOT YET VALIDATED against dist/bin -- see the owning test's docstring.
    # The SAME conjunctive-satisfying shape as both_codeclared (master tag
    # 10 shared; NTS `-slave` node set == exactly {1, 2}), except the
    # mortar-side declaration is `-tie` instead of plain penalty: an
    # equality bond cannot be owned by a unilateral NTS contact, so the
    # ALIGN warn must fire regardless of the co-declared NTS contact.
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -1.0, 0.0)
    ops.node(102, 1.0, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.contactSurface(10, "-master", 2, 101, 102)

    ops.node(1, 0.0, -SEED)
    ops.node(2, 0.0, -SEED - 1.0)
    ops.fix(1, 1, 0)
    ops.fix(2, 1, 1)

    ops.node(901, 0.0, -SEED)
    ops.fix(901, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, KS)
    ops.element("zeroLength", 1, 901, 1, "-mat", 1, "-dir", 2)

    ops.contactSurface(20, "-slave-segments", 2, 1, 2)
    ops.contact(1, 10, 20, "-mortar", "-tie", "-epsTie", EPSN, "-outward", 0.0, 1.0)
    ops.contactSurface(21, "-slave", 1, 2)
    ops.contact(2, 10, 21, KN, 0.0, 0.0, "-outward", 0.0, 1.0)

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
    assert rc == 0, "tie_align_always_warns did not converge rc=" + str(rc)
    print("ROUTE_RESULT tie_align_always_warns u=" + repr(ops.nodeDisp(1, 2)))
else:
    print("UNKNOWN CASE " + CASE)
    sys.exit(4)
sys.exit(0)
'''


def _run_route_child(case):
    script = CHILD_ROUTE % {"ENGINE_DIR": ENGINE_DIR}
    env = dict(os.environ)
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    proc = subprocess.run(
        [sys.executable, "-c", script, case],
        stdin=subprocess.DEVNULL,
        capture_output=True, text=True, timeout=300,
        encoding="utf-8", errors="replace",
        env=env,
    )
    return proc, (proc.stdout or "") + "\n" + (proc.stderr or "")


def _grab_u(out, marker):
    for line in out.splitlines():
        if line.startswith(marker):
            for tok in line.split():
                if tok.startswith("u="):
                    return float(tok[2:])
    raise AssertionError(f"no u= on a {marker!r} line:\n{out}")


def test_2d_perp_pair_nts_codeclared_transmits_no_warn():
    """G-T3(d): NTS and mortar declared on the SAME interface node -> NTS owns
    the transfer (its closed form, kn dominant over the soft spring), no
    double-count (mortar contributes nothing on the perpendicular-refused
    pair), and WARN_2D_PERP_NO_NTS must NOT fire.

    THE OWNERSHIP RULE IS CONJUNCTIVE (orchestrator update, post-review round
    -- this correction landed after the first gate rerun, which tried three
    "same surfaces" readings that each satisfied only one half and reported
    the mismatch rather than guessing further). An NTS contact suppresses the
    mortar pair's ALIGN warn iff BOTH hold simultaneously:
      1. the NTS contact's MASTER surface tag == the mortar contact's MASTER
         surface tag (declaration-tag identity on the master side), AND
      2. the UNION of that NTS contact's slave NODE tags contains BOTH node
         tags of the specific mortar SLAVE segment that refuses ALIGN.
    `_build()` above declares mortar with a FIXED master/slave-surface/
    contact-tag triple (10/20/1) regardless of whether NTS is also present,
    and -- when present -- NTS with master tag 10 (same as mortar) and a
    `-slave` NODE set containing EXACTLY {1, 2} (both of mortar's slave
    segment's nodes, not just the touching one) -- confirmed directly against
    the live dist/bin build to suppress the warning.
    """
    proc_n, out_n = _run_route_child("nts_only")
    proc_b, out_b = _run_route_child("both_codeclared")
    assert proc_n.returncode == 0, f"nts_only failed:\n{out_n}"
    assert proc_b.returncode == 0, f"both_codeclared failed:\n{out_b}"

    u_n = _grab_u(out_n, "ROUTE_RESULT nts_only")
    u_b = _grab_u(out_b, "ROUTE_RESULT both_codeclared")
    assert abs(u_n - u_b) <= 1.0e-9 * max(abs(u_n), 1.0), (
        f"co-declaring a (refused) mortar pair alongside NTS changed the "
        f"equilibrium: nts_only u={u_n!r} vs both_codeclared u={u_b!r} -- "
        "either double-counting or the NTS stand-down is wrong")

    # WARN_2D_PERP_NO_NTS's latch text must be ABSENT: the pair is covered by
    # the co-declared NTS contact, so nothing should warn.
    assert not ("WARN" in out_b and "NTS" in out_b and " 1" in out_b), (
        "WARN_2D_PERP_NO_NTS fired even though an NTS contact was declared on "
        f"the same interface node -- it should be silently covered:\n{out_b}")


def test_2d_tie_align_always_warns_regardless_of_nts():
    """G-T3(d) post-review addendum: a `-tie` pair (mesh-tying, an EQUALITY
    bond) that refuses ALIGN (near-perpendicular) warns even when an NTS
    contact is co-declared on the same interface in the EXACT shape that
    silences the ordinary (penalty) mortar case above -- because a unilateral
    NTS contact cannot substitute for an equality bond the way it substitutes
    for a one-sided penalty contact. The message additionally contains
    "mesh-tie bond cannot be owned by an NTS contact" (tie-specific text,
    distinct from the plain WARN_2D_PERP_NO_NTS wording).
    """
    proc, out = _run_route_child("tie_align_always_warns")
    assert proc.returncode == 0, f"tie_align_always_warns failed:\n{out}"
    assert "WARN" in out and "NTS" in out, (
        f"the tie-ALIGN warn did not fire even though the pair is "
        f"near-perpendicular:\n{out}")
    assert "mesh-tie bond cannot be owned by an NTS contact" in out, (
        f"the tie-ALIGN warn fired but is missing the tie-specific "
        f"substring:\n{out}")


def test_2d_perp_pair_mortar_only_warns_once():
    """G-T3(d): mortar-only on a perpendicular pair -> the pair is inert (the
    soft spring alone equilibrates the load) AND `WARN_2D_PERP_NO_NTS` fires
    -- exactly ONCE across two analyze() calls (the warnOnce contract).
    """
    proc, out = _run_route_child("mortar_only")
    assert proc.returncode == 0, f"mortar_only failed:\n{out}"

    u = _grab_u(out, "ROUTE_RESULT mortar_only")
    u_ref = -1.0 / 1.0e0     # -P/KS: the spring alone carries the load
    assert abs(u - u_ref) / abs(u_ref) < 1.0e-6, (
        f"mortar-only perpendicular pair transmitted something: u={u!r} vs "
        f"pure-spring {u_ref!r} -- the perpendicular refusal should leave it "
        "inert")

    assert "WARN" in out and "NTS" in out, (
        f"WARN_2D_PERP_NO_NTS did not fire for an uncovered perpendicular "
        f"mortar pair:\n{out}")
    assert "1" in out, (
        f"the WARN_2D_PERP_NO_NTS message does not name the contact tag "
        f"(1, mortar's own fixed contact tag in _build()):\n{out}")
    n_fires = sum(1 for line in out.splitlines() if "WARN" in line and "NTS" in line)
    assert n_fires == 1, (
        f"WARN_2D_PERP_NO_NTS fired {n_fires} times across two analyze() "
        f"calls on the SAME pair -- the latch must be warnOnce:\n{out}")


# --------------------------------------------------------------------------
# the T3 scope fences -- child process (refusal pattern, T0's shape)
# --------------------------------------------------------------------------
FENCES = {
    "soft_on_2d_mortar": ("ADR-85", "SOFT=2"),
    "edgeedge_on_2d_mortar": ("no 2D analogue", "scope fence"),
    "thickness_without_mortar": ("-thickness",),
    "thickness_nonpositive": ("-thickness",),
    "thickness_on_3d_mortar": ("-thickness", "2D plane-model option"),
}

CHILD_FENCE = r'''
import os, sys

_D = %(ENGINE_DIR)r
if os.path.isdir(_D):
    os.environ["PATH"] = _D + os.pathsep + os.environ.get("PATH", "")
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
EPSN = 1.0e6


def _static():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-10, 20, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    return ops.analyze(1)


rc = 0

if CASE == "soft_on_2d_mortar":
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -1.0, 0.0)
    ops.node(102, 1.0, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.node(1, 0.0, -1.0e-4)
    ops.node(2, 1.0, -1.0e-4)
    ops.fix(1, 1, 0)
    ops.fix(2, 1, 0)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.contactSurface(20, "-slave-segments", 2, 1, 2)
    ops.contact(1, 10, 20, "-mortar", "-epsN", EPSN, "-soft", 0.1,
                "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, -1.0)
    ops.load(2, 0.0, -1.0)
    rc = _static()

elif CASE == "edgeedge_on_2d_mortar":
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -1.0, 0.0)
    ops.node(102, 1.0, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.node(1, 0.0, -1.0e-4)
    ops.node(2, 1.0, -1.0e-4)
    ops.fix(1, 1, 0)
    ops.fix(2, 1, 0)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.contactSurface(20, "-slave-segments", 2, 1, 2)
    ops.contact(1, 10, 20, "-mortar", "-epsN", EPSN, "-edgeedge",
                "-edgeBand", 0.15, "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, -1.0)
    ops.load(2, 0.0, -1.0)
    rc = _static()

elif CASE == "thickness_without_mortar":
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -1.0, 0.0)
    ops.node(102, 1.0, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.node(1, 0.0, -1.0e-4)
    ops.fix(1, 1, 0)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, EPSN, 0.0, 0.0, "-thickness", 2.0,
                "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, -1.0)
    rc = _static()

elif CASE == "thickness_nonpositive":
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -1.0, 0.0)
    ops.node(102, 1.0, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.node(1, 0.0, -1.0e-4)
    ops.node(2, 1.0, -1.0e-4)
    ops.fix(1, 1, 0)
    ops.fix(2, 1, 0)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.contactSurface(20, "-slave-segments", 2, 1, 2)
    ops.contact(1, 10, 20, "-mortar", "-epsN", EPSN, "-thickness", 0.0,
                "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, -1.0)
    ops.load(2, 0.0, -1.0)
    rc = _static()

elif CASE == "thickness_on_3d_mortar":
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, x, y in [(1, 0.0, 0.0), (2, 1.0, 0.0), (3, 1.0, 1.0), (4, 0.0, 1.0)]:
        ops.node(t, x, y, 0.0)
        ops.fix(t, 1, 1, 1)
    for t, x, y in [(5, 0.0, 0.0), (6, 1.0, 0.0), (7, 1.0, 1.0), (8, 0.0, 1.0)]:
        ops.node(t, x, y, -1.0e-4)
        ops.fix(t, 1, 1, 0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPSN, "-thickness", 2.0,
                "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -1.0)
    rc = _static()

else:
    print("UNKNOWN CASE " + CASE)
    sys.exit(4)

if rc == 0:
    print("NOREFUSAL " + CASE)
    sys.exit(0)

print("ANALYZE_REFUSED " + CASE + " rc=" + str(rc))
sys.exit(3)
'''


def _run_fence_child(case):
    script = CHILD_FENCE % {"ENGINE_DIR": ENGINE_DIR}
    env = dict(os.environ)
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    proc = subprocess.run(
        [sys.executable, "-c", script, case],
        stdin=subprocess.DEVNULL,
        capture_output=True, text=True, timeout=300,
        encoding="utf-8", errors="replace",
        env=env,
    )
    return proc


def _assert_fenced(case):
    proc = _run_fence_child(case)
    out = (proc.stdout or "") + "\n" + (proc.stderr or "")
    assert proc.returncode != 0, (
        f"{case}: the child exited 0 -- ADR-85 T3 must refuse it.\n{out}")
    assert "NOREFUSAL" not in out, f"{case}: the deck ran without a refusal.\n{out}"
    for needle in FENCES[case]:
        assert needle in out, (
            f"{case}: the refusal fired but is missing the substring "
            f"{needle!r} (see this file's message-contract table).\n{out}")


def test_2d_soft_on_2d_mortar_refused():
    """`-soft` (SOFT=2) on a 2D mortar deck is refused by name.
    ADR-85 design note §2b.6: "2D mortar + -soft: SOFT=2 is fenced out of 2D
    ... named not-supported refusal, never silent."
    """
    _assert_fenced("soft_on_2d_mortar")


def test_2d_edgeedge_on_2d_mortar_refused():
    """`-edgeedge` on a 2D mortar deck is refused with the PERMANENT scope-
    fence message (edge-edge has no 2D analogue, ADR-85 design note §2 recon).
    """
    _assert_fenced("edgeedge_on_2d_mortar")


def test_2d_thickness_without_mortar_refused():
    """`-thickness` requires `-mortar` (ADR-85 design note §2b.5)."""
    _assert_fenced("thickness_without_mortar")


def test_2d_thickness_nonpositive_refused():
    """`-thickness h` with h <= 0 is refused (parse-time h > 0 validation)."""
    _assert_fenced("thickness_nonpositive")


def test_2d_thickness_on_3d_mortar_refused():
    """`-thickness` on a well-formed 3D mortar deck is refused by name -- it
    is a 2D plane-model-only option (ADR-85 design note §2b.5: "handle-time
    named FATAL if the pair is 3D and h != 1.0").
    """
    _assert_fenced("thickness_on_3d_mortar")
