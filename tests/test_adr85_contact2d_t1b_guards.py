"""ADR-85 T1b -- the guard matrix the merge panel sustained.

WHAT THIS FILE IS
    T0's refusal matrix (tests/test_adr85_contact2d_t0_refusals.py) gated the
    DECLARATION surface: dimensions, arities, not-yet-supported lanes.  T1b built
    the machinery behind those refusals -- the interface-level orientation vote,
    the master chain, the cross-surface dimension peek -- and the merge panel
    found that each of them had a silent-wrong mode with no test on it.  This
    file is those rows.  Every one of them was a PROBE-CONFIRMED silent failure,
    not a hypothetical; the probe result is recorded in each test's docstring so
    the next reader knows what the row is buying.

    Each row is a FALSIFIER (a named abort) paired, where one exists, with a
    CONTROL that must still run.  Without the controls a guard implemented as
    "refuse everything interesting" would pass the whole matrix -- the same
    argument the T0 file makes for its two acceptance rows.

MESSAGE CONTRACT
    Exit code alone is too weak (a deck that merely fails to converge also exits
    nonzero), so every falsifier asserts SUBSTRINGS on the child's combined
    output.  `ADR-85` is required on every row, following the shipped precedent
    that a refusal names the ADR that owns it; the second substring names the
    specific rule.  The substrings below are KEYWORD contracts fixed by the
    orchestrator, not guesses at prose -- the engine may word the sentence any
    way it likes as long as these tokens appear.

WHY A CHILD PROCESS
    Every row here asserts a FATAL, and contact refusals are not ordinary Python
    exceptions: a declaration refusal returns -1 (an OpenSeesError in serial, an
    MPI_Abort under np > 1) and a handler-time refusal surfaces later still, as a
    nonzero `analyze()` after `domainChanged()` fails.  A single process cannot
    portably assert on all three shapes and one of them takes the interpreter
    with it.  The pattern, including both Windows traps (the PYTHONPATH
    environment-block trap and the stdin=DEVNULL Popen trap), is documented in
    the module docstring of tests/test_adr85_contact2d_t0_refusals.py, from which
    ENGINE_DIR is imported so the child loads the SAME binary the parent bound.
"""
import os
import subprocess
import sys

import pytest

from test_adr85_contact2d_t0_refusals import ENGINE_DIR

pytestmark = [pytest.mark.zone_a]


# --------------------------------------------------------------------------
# THE CONTRACT.  case -> substrings that MUST appear in the child's output.
# --------------------------------------------------------------------------
REFUSALS = {
    # THE BLOCKER.  The T0 acceptance deck (pair_2d_now_live) with `-outward`
    # REMOVED.  Its slave is seeded 1e-8 below a flush master, so the
    # interface-level centroid datum is 1e-8 long against a segment of length 1
    # -- pure seed noise, not a direction.  The magnitude gate
    # (|refDir| > 1e-6*Lref) is what turns that into the named degenerate-vote
    # abort instead of a coin flip.
    "vote_seed_noise_refused": ("ADR-85", "outward"),

    # Master segments listed OUT OF CHAIN ORDER: (0-1)(2-3)(1-2).
    "chain_permuted": ("ADR-85", "chain"),
    # One segment listed with its endpoints REVERSED: (0-1)(2-1)(2-3).
    "chain_reversed": ("ADR-85", "chain"),

    # A correctly-chained, consistently-wound master whose outward normals span
    # 180 degrees (a U-channel): no single centroid datum can orient all three,
    # so the vote SPLITS and must refuse rather than orient part of the surface
    # backwards.
    "split_vote_u_channel": ("ADR-85", "outward"),

    # The mirror of T0's cross_dim_pair: 2D master (nps = 2) + 3D slave node
    # set.  T0 covered 3D-master + 2D-slave; this is the other order, and it is
    # the one that reaches the `-outward` arity peek first.
    "cross_dim_pair_mirror": ("ADR-85", "dimension"),
}

# cases that must RUN (exit 0) -- the falsifiers on the falsifiers
CONTROLS = (
    "vote_seed_noise_control",
    "chain_ordered_control",
    "split_vote_straight_control",
)


# --------------------------------------------------------------------------
# the child.  NOTE: this template is rendered by the percent operator against a
# one-key mapping, so NOTHING inside it -- code or comment -- may contain a bare
# percent conversion; use repr() and concatenation.
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

KN, KS, P = 1.0e6, 1.0e3, 1.0e3
SEED = 1.0e-8
E, NU, THICK, W, H, PATCH_P = 1.0e4, 0.0, 1.0, 2.0, 1.0, 2.0e2


def _solve():
    """One static step.  A HANDLER-time refusal surfaces here as a nonzero
    analyze() after domainChanged() fails; a declaration-time one has already
    raised before we get here."""
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-10, 40, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.integrator("LoadControl", 1.0)
    return ops.analyze(1)


def _tributary(n, width):
    le = width / n
    return [le / 2.0 if j in (0, n) else le for j in range(n + 1)]


def _patch(nbot, ntop, mseg):
    """The T1b patch deck (tests/test_adr85_contact2d_t1b_nts.py::_patch_model)
    with the master SEGMENT LIST supplied by the caller, so the only thing that
    varies between the chain legs is the order the segments are declared in."""
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    base, mid = [], []
    for i in range(nbot + 1):
        x = W * i / nbot
        ops.node(100 + i, x, 0.0)
        ops.fix(100 + i, 1, 1)
        ops.node(110 + i, x, H)
        ops.fix(110 + i, 1, 0)
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
    tags = []
    for a, b in mseg:
        tags += [mid[a], mid[b]]
    ops.contactSurface(10, "-master", 2, *tags)
    ops.contactSurface(20, "-slave", *slave)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, w in zip(top, _tributary(ntop, W)):
        ops.load(t, 0.0, -PATCH_P * w / W)
    return slave


def _tether(slave, x, y):
    ops.node(900, x, y)
    ops.fix(900, 1, 1)
    ops.uniaxialMaterial("Elastic", 900, KS)
    ops.element("zeroLength", 900, 900, slave, "-mat", 900, 900, "-dir", 1, 2)


rc = 0
report = ""

if CASE in ("vote_seed_noise_refused", "vote_seed_noise_control"):
    # The T0 acceptance deck, byte-for-byte, MINUS `-outward` on the refused leg.
    # Master line y = 0 from -0.5 to 0.5; slave seeded SEED below it, held by a
    # soft y-spring so the deck is solvable with or without a live pair.  The
    # centroid datum is therefore (0, -SEED): 1e-8 long against Lref = 1.
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -0.5, 0.0)
    ops.node(102, 0.5, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.node(1, 0.0, -SEED)
    ops.fix(1, 1, 0)
    ops.node(2, 0.0, -SEED)
    ops.fix(2, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, KS)
    ops.element("zeroLength", 1, 2, 1, "-mat", 1, "-dir", 2)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.contactSurface(20, "-slave", 1)
    if CASE.endswith("control"):
        ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 1.0)
    else:
        ops.contact(1, 10, 20, KN, 0.0, 0.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, -P)
    rc = _solve()
    if rc == 0:
        report = (" u=" + repr(ops.nodeDisp(1, 2))
                  + " f=" + repr(ops.ladrunoContactForce(1)))

elif CASE.startswith("chain_"):
    ORDERS = {
        "chain_ordered_control": [(0, 1), (1, 2), (2, 3)],
        "chain_permuted":        [(0, 1), (2, 3), (1, 2)],
        "chain_reversed":        [(0, 1), (2, 1), (2, 3)],
    }
    slave = _patch(3, 5, ORDERS[CASE])
    rc = _solve()
    if rc == 0:
        tot = 0.0
        for t in slave:
            tot += ops.ladrunoContactForce(t)
        report = " total=" + repr(tot)

elif CASE.startswith("split_vote"):
    # U-CHANNEL: a consistently-wound 3-segment chain whose outward normals are
    # +x (left wall), +y (floor), -x (right wall) -- 180 degrees of spread.  The
    # STRAIGHT control is the same chain flattened, all normals +y.
    if CASE == "split_vote_u_channel":
        pts = [(-1.0, 1.0), (-1.0, 0.0), (1.0, 0.0), (1.0, 1.0)]
        sx, sy = 0.3, 0.2                      # inside the channel
    else:
        pts = [(-1.0, 0.0), (-1.0 / 3.0, 0.0), (1.0 / 3.0, 0.0), (1.0, 0.0)]
        sx, sy = 0.0, 0.2                      # clearly above -> the vote resolves
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for i, (x, y) in enumerate(pts):
        ops.node(101 + i, x, y)
        ops.fix(101 + i, 1, 1)
    ops.node(1, sx, sy)
    _tether(1, sx, sy)
    ops.contactSurface(10, "-master", 2, 101, 102, 102, 103, 103, 104)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0)       # NO -outward: the vote decides
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    rc = _solve()
    if rc == 0:
        report = (" ux=" + repr(ops.nodeDisp(1, 1))
                  + " uy=" + repr(ops.nodeDisp(1, 2)))

elif CASE == "cross_dim_pair_mirror":
    # 2D master (nps = 2, 2-coordinate nodes) + 3D slave node set.  Each surface
    # is internally consistent, so a PER-SURFACE pre-flight passes both; only a
    # cross-surface check can refuse the pair.  `-outward` is given with THREE
    # components because the slave surface is 3D -- which is exactly the peek
    # that must not decide the arity from one surface alone.  Without the
    # cross-check the deck dies reading an unexpected token, naming a layout
    # nobody wrote, instead of naming the dimension mismatch.
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(101, -0.5, 0.0)
    ops.node(102, 0.5, 0.0)
    ops.fix(101, 1, 1)
    ops.fix(102, 1, 1)
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, -SEED)
    ops.fix(1, 1, 1, 0)
    ops.contactSurface(10, "-master", 2, 101, 102)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -1.0)
    rc = _solve()

else:
    print("UNKNOWN CASE " + CASE)
    sys.exit(4)

if rc == 0:
    print("RAN " + CASE + report)
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


def _grab(out, key):
    for line in out.splitlines():
        if not line.startswith("RAN "):
            continue
        for tok in line.split():
            if tok.startswith(key + "="):
                return float(tok[len(key) + 1:])
    raise AssertionError(f"no {key!r} on the RAN line:\n{out}")


def _assert_refused(case):
    proc, out = _run_child(case)
    assert proc.returncode != 0, (
        f"{case}: the child exited 0 -- the deck was ACCEPTED and ran. "
        f"ADR-85 T1b must refuse it.\n{out}")
    assert "RAN " not in out, (
        f"{case}: the deck assembled and solved without a refusal.\n{out}")
    for needle in REFUSALS[case]:
        assert needle in out, (
            f"{case}: the refusal fired but does not name itself -- missing "
            f"{needle!r}.  See the REFUSALS contract at the top of this file; "
            f"these are KEYWORD contracts, the prose around them is free.\n{out}")


def _assert_ran(case):
    proc, out = _run_child(case)
    assert proc.returncode == 0, (
        f"{case}: this CONTROL must still run (exit {proc.returncode}).  A guard "
        "implemented as 'refuse everything' would pass every falsifier in this "
        f"file without it.\n{out}")
    assert "RAN " + case in out, f"{case}: the child did not finish.\n{out}"
    return out


# --------------------------------------------------------------------------
# 1. THE BLOCKER -- the degenerate-vote magnitude gate
# --------------------------------------------------------------------------
def test_2d_vote_seed_noise_refused():
    """A flush interface whose slave is SEEDED 1e-8 into the master, with no
    `-outward` -> named degenerate-vote FATAL.

    THIS IS THE ROW BOTH FIX-FIRST PANELISTS DEMANDED, and the probe is why.
    Without a magnitude gate on the centroid datum, the vote does not see a
    degenerate interface -- it sees a direction of length 1e-8 and happily
    normalises it.  The probe ran exactly this deck and it went GREEN: it
    converged, it balanced, and `ladrunoContactForce` was 0.0.  Silent nothing,
    on the very deck the T0 matrix uses as its acceptance row.

    That is what makes this the blocker rather than a nicety.  Seeding the slave
    a hair inside the master is not an exotic input -- it is the idiom the whole
    battery uses (and must use: a slave held only by a penalty is singular while
    separated), so the ONE deck shape every 2D contact user will write first is
    the one that silently transmitted nothing.

    The fix is the magnitude gate |refDir| > 1e-6*Lref: 1e-8 against Lref = 1
    fails it, the interface is declared degenerate, and the abort names
    `-outward` as the override.  The CONTROL below proves the gate did not
    simply swallow the deck.
    """
    _assert_refused("vote_seed_noise_refused")


def test_2d_vote_seed_noise_control():
    """THE CONTROL on the gate: the identical deck WITH `-outward 0 1` must run,
    and must run to the CONTACT answer rather than the spring-only one.

    Checked against closed forms, not just exit 0, because the failure the
    blocker row describes is precisely a green run with the wrong number:

        contact live   u = (KN*SEED - P)/(KN + KS),  f = KN*(SEED - u)
        silent nothing u = -P/KS,                    f = 0.0

    The two differ by three orders of magnitude, so there is no way to satisfy
    this assertion with an inert pair.
    """
    out = _assert_ran("vote_seed_noise_control")
    u, f = _grab(out, "u"), _grab(out, "f")
    KN_, KS_, P_, SEED_ = 1.0e6, 1.0e3, 1.0e3, 1.0e-8
    u_ref = (KN_ * SEED_ - P_) / (KN_ + KS_)
    f_ref = KN_ * (SEED_ - u_ref)
    assert abs(u - u_ref) / abs(u_ref) < 1.0e-8, (
        f"oriented flush deck settled at {u!r}, closed form {u_ref!r} "
        f"(the spring-only silent-nothing answer is {-P_ / KS_!r})\n{out}")
    assert f > 0.0 and abs(f - f_ref) / f_ref < 1.0e-6, (
        f"contact force {f!r} != kn*penetration {f_ref!r}\n{out}")


# --------------------------------------------------------------------------
# 2. CHAIN INTEGRITY
# --------------------------------------------------------------------------
def test_2d_chain_permuted_refused():
    """Master segments declared OUT OF CHAIN ORDER -> named FATAL naming `chain`.

    Deck: the T1b patch geometry with a 3-segment master, declared (0-1)(2-3)
    (1-2) instead of (0-1)(1-2)(2-3).  Every segment is present, every node is
    right, only the ORDER is wrong -- which is what a hand-written deck or a
    mesher emitting an unsorted edge set produces.

    THE PROBE RESULT IS WHY THIS IS A FATAL AND NOT A WARNING.  A permuted
    listing did not fail, and it did not produce noise: it produced the exact
    DOUBLE-OWNER closed form -- twice the interface stiffness -- because the
    ownership scan walks the surface in declaration order and a discontinuous
    walk lets two segments claim the same slave.  A result that is exactly 2x
    and perfectly converged is the worst shape a wrong answer can take; nothing
    downstream flags it.
    """
    _assert_refused("chain_permuted")


def test_2d_chain_reversed_refused():
    """One master segment declared with its endpoints REVERSED -> the same named
    `chain` FATAL.

    Same deck, listing (0-1)(2-1)(2-3): the segment is geometrically identical
    and the chain is even still traversable if you are willing to flip it, which
    is the trap.  A reversed segment inverts that segment's derived normal
    relative to its neighbours, so on a curved or stepped master it produces the
    same class of half-backwards interface the split-vote row covers -- but
    silently, because each individual segment is well-formed.  The endpoint
    continuity rule (`seg[i].end == seg[i+1].start`) catches both this and the
    permutation, which is why they share a keyword.
    """
    _assert_refused("chain_reversed")


def test_2d_chain_ordered_control():
    """THE CONTROL: the same 3-segment master declared IN ORDER must run and
    transmit the full applied load.

    Without it, a chain check implemented as "refuse any multi-segment master"
    would pass both falsifiers above while deleting the capability the patch
    test depends on.  The transmitted total is asserted (not just exit 0) so a
    chain check that accepts the declaration and then drops the pairs is caught
    here rather than in the patch test's tolerance.
    """
    out = _assert_ran("chain_ordered_control")
    tot = _grab(out, "total")
    assert abs(tot - 2.0e2) / 2.0e2 < 1.0e-6, (
        f"the in-order 3-segment chain transmitted {tot!r}, applied 200.0 "
        f"-- accepted but not assembled\n{out}")


# --------------------------------------------------------------------------
# 3. SPLIT VOTE
# --------------------------------------------------------------------------
def test_2d_split_vote_refused():
    """A correctly-chained master whose normals span 180 degrees -> named FATAL
    naming `outward`.

    Deck: a U-channel, master nodes (-1,1) (-1,0) (1,0) (1,1) chained
    (0-1)(1-2)(2-3).  The winding is CONSISTENT -- this deck passes every chain
    check -- and the outward normals are +x, +y, -x.  A single reference datum
    cannot have a positive dot with all three: with the slave inside the channel
    at (0.3, 0.2) the datum is (0.3, -0.133), which votes IN for the left wall
    and OUT for the floor and the right wall.

    This is the mode the magnitude gate cannot catch (the datum here is 0.33
    long, four orders above the gate) and the chain check cannot catch (the
    chain is perfect).  Left unrefused it orients part of the surface backwards
    -- the segments that lost the vote push slaves INTO the master -- which is a
    converging, energy-plausible, completely wrong answer.

    A channel, a notch and a re-entrant corner are the geometries the 2D lane
    exists for (masonry joints, footings, keyed shear connections), so this is
    not a pathological input; `-outward` cannot fix it either, which is why the
    refusal is the right answer and per-segment orientation is a T3+ question.
    """
    _assert_refused("split_vote_u_channel")


def test_2d_split_vote_straight_control():
    """THE CONTROL: the same 3-segment chain FLATTENED (all normals +y), slave
    clearly above, still no `-outward` -> the vote resolves and the deck runs.

    This isolates SPLIT from MULTI-SEGMENT.  A vote that refused every chain of
    more than one segment, or refused whenever `-outward` is absent, would pass
    the row above; here the datum (0, 0.2) votes IN for all three segments and
    the deck must be accepted.  The slave has a positive gap, so the claim is
    exactly "no spurious refusal" -- the engaged-contact claims live in the NTS
    file's vote tests.
    """
    _assert_ran("split_vote_straight_control")


# --------------------------------------------------------------------------
# 5. CROSS-DIMENSION, THE OTHER ORDER
# --------------------------------------------------------------------------
def test_2d_cross_dim_pair_mirror_refused():
    """A 2D master (nps = 2) paired with a 3D `-slave` node set -> named
    cross-surface `dimension` FATAL.

    T0's cross_dim_pair covers 3D-master + 2D-slave.  This is the mirror, and it
    is not redundant: it is the order that reaches the `-outward` ARITY PEEK
    first.  The parser has to decide whether `-outward` carries 2 or 3 doubles
    before it can finish reading the command, and if it decides that from the
    master surface alone it reads 2, leaves `1.0` on the stack, and dies
    complaining about an unexpected token -- naming an argument layout the user
    never wrote, for a deck whose real problem is that the two surfaces do not
    live in the same space.

    So the assertion is deliberately on the MESSAGE, not merely on the nonzero
    exit: the pre-fix build already fails this deck, just uselessly.  The peek
    must consult BOTH surfaces and refuse the mismatch by name.
    """
    _assert_refused("cross_dim_pair_mirror")
