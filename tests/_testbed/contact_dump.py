"""ADR-85 T0 -- the contact byte-identity harness (the G-T0 observable).

WHY THIS EXISTS
    The shipped contact battery (ADR-39 / 41 / 57) is pass/fail WITH TOLERANCE.
    It cannot certify that a change left the 3D lanes BIT-FOR-BIT alone: a
    rewrite that perturbs the last two digits of every displacement passes every
    one of those asserts.  ADR-85's rev-2 adversarial panel called that out (see
    the ADR's Review record, finding 6): "byte-identical battery" was a prose
    claim with no observable behind it.  This script IS the observable.

    It runs ~7 small, deterministic, SERIAL 3D decks -- one per shipped contact
    lane -- and writes their full-precision final-state observables to a text
    artifact.  Two builds whose contact physics is unchanged must produce
    artifacts that are IDENTICAL AS BYTES (`fc /b`, `cmp`, `git diff --no-index`).

USAGE
    python tests/_testbed/contact_dump.py --out artifact.txt

    exit 0 : every deck ran, artifact written
    exit 1 : at least one deck FAILED (its failure is still recorded in the
             artifact, deterministically, so the diff stays meaningful)
    exit 2 : the harness could not bind this worktree's engine (never a physics
             signal -- fix the build/bootstrap and re-run)

THE ARTIFACT CONTAINS PHYSICS AND NOTHING ELSE
    No build stamp, no timestamp, no path, no version string, no wall-clock, no
    step timings.  Anything that legitimately differs between two builds of the
    SAME physics would make every diff red and the gate worthless.  The
    `ladrunoBuild()` stamp, the engine path and the per-deck progress go to
    STDERR, where the caller records them alongside the artifact (the ADR's
    "both ladrunoBuild stamps recorded" clause).

    Line format, one value per line, ASCII only:

        <deck> <kind> <tag> <component> <float.hex()>

    `float.hex()` is the exact IEEE-754 double -- no rounding, no locale, no
    repr-shortening.  `-` fills a slot a kind does not use (an engine-wide
    scalar has no tag and no component).  Kinds:

        disp          nodal displacement component (1-based dof), final
                      COMMITTED state
        contactForce  ladrunoContactForce(slave) -- NTS lanes only; the query is
                      fed from LadrunoContactFE's SEGMENT branch
                      (LadrunoContactFE.cpp:689), so the rigid-plane / mortar /
                      tie / edge decks would only ever emit a constant 0.0 and
                      it is deliberately NOT emitted for them
        mortarPenPenalty / mortarPenAugmented
                      ladrunoMortarPenetration(), before and after the fixed
                      held-load Uzawa sweep (two separate kinds so a diff names
                      which one moved)
        tieRes        ladrunoMortarTieResidual()
        edgePen       ladrunoEdgePenetration()

    Ordering is fixed: decks in DECKS order, then each deck's records in the
    order its builder appends them (nodes in an explicit sorted tag list, dofs
    ascending).  Nothing iterates a dict/set whose order could drift.

    A deck that fails emits exactly one line

        DECK <name> FAILED <code>

    and the run exits 1 AFTER attempting every remaining deck.  A deterministic
    failure is still a diffable artifact; a half-written one is not.

DECK PROVENANCE -- every deck is copied from a SHIPPED, PASSING test, so its
syntax and its convergence are already known-good.  Only the step count and the
load-step size were changed (5..20 fixed steps, no tolerance-driven loops, no
randomness, no wall-clock):

    rigid_plane_static   tests/test_adr39_contact_p2a.py
                         ::test_p2a_axis_aligned_penetration_static
    nts_penalty_static   tests/test_adr39_contact_p2b.py
                         ::test_p2b_static_penetration_quad
    nts_autokn_static    tests/test_adr39_contact_p2b2b.py::_run_autokn
    nts_friction_expl    tests/test_adr39_contact_p3_friction.py::_driven_block
    mortar_alm           tests/test_adr41_mortar_c2_2.py
                         ::_rigid_master_slave_facet
    mortar_tie           tests/test_adr41_mortar_c4_1.py::_split_column
    edge_edge            tests/test_adr57_edge_edge_1.py
                         ::test_e2_crossed_bars_restrained  (+ -edgeAlm, so the
                         ladrunoEdgePenetration query is live -- the E6 idiom)

    mortar_alm's held-load augmentation is a FIXED count of Uzawa commits, not
    Ladruno_scripts/analyze_augmented's tolerance-driven loop: a loop that stops
    on a tolerance can change its ITERATION COUNT across builds, which changes
    the SHAPE of the artifact and buries the actual physics diff in a line-count
    diff.  The bracket is still the ADR-41 D1 one (ladrunoBeginAugment /
    ladrunoEndAugment around zero-increment LoadControl re-solves).
"""
import argparse
import os
import sys
from pathlib import Path

# --------------------------------------------------------------------------
# bootstrap: THIS worktree's engine, ahead of everything.
#
# LADRUNO_OPENSEES_QUIET is set BEFORE the import: PyInit prints the splash
# banner otherwise, and the banner's UTF-8 box glyphs break a text-mode
# subprocess capture on a cp1252 Windows console (LEDGER_quirks: "The splash
# banner kills naive subprocess capture on Windows").  A harness whose whole
# job is to be captured must not emit it.
# --------------------------------------------------------------------------
_HERE = Path(__file__).resolve()
_TESTS = _HERE.parents[1]                       # <worktree>/tests
_ROOT = _HERE.parents[2]                        # <worktree>
_DIST = _ROOT / "dist" / "bin"

os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

if not ((_DIST / "opensees.pyd").is_file() or (_DIST / "opensees.so").is_file()):
    sys.stderr.write(
        "contact_dump: FATAL - no engine at %s (build `OpenSees OpenSeesPy` "
        "in THIS worktree first)\n" % _DIST)
    raise SystemExit(2)

if str(_TESTS) not in sys.path:
    sys.path.insert(0, str(_TESTS))

try:
    from _engine import bind_worktree_engine
    ops = bind_worktree_engine(str(_DIST))
except Exception as exc:                        # noqa: BLE001 - report and die loudly
    sys.stderr.write("contact_dump: FATAL - could not bind the worktree engine: "
                     "%s: %s\n" % (type(exc).__name__, exc))
    raise SystemExit(2)

if not hasattr(ops, "ladrunoBuild"):
    sys.stderr.write(
        "contact_dump: FATAL - the bound module has no ladrunoBuild(): this is "
        "not a Ladruno fork build (an upstream openseespy wheel cannot run the "
        "contact lanes at all). Engine: %s\n" % getattr(ops, "__file__", "?"))
    raise SystemExit(2)


# --------------------------------------------------------------------------
# record plumbing
# --------------------------------------------------------------------------
class DeckFailure(Exception):
    """A deck did not run to completion.  `code` is the analyze() return code
    (or -999 for a raised exception), and it lands in the artifact so that even
    a failure is diffable."""

    def __init__(self, code):
        super().__init__("deck failed with code %s" % code)
        self.code = code


def _step(rc):
    """analyze() return-code check -- every deck routes through it."""
    if rc != 0:
        raise DeckFailure(rc)


class Recorder(object):
    """Accumulates (kind, tag, component, value) in append order."""

    def __init__(self):
        self.rows = []

    def add(self, kind, tag, comp, value):
        self.rows.append((kind, str(tag), str(comp), float(value)))

    def disps(self, node_tags, ndf):
        """Every listed node's displacement components at the final COMMITTED
        state.  `node_tags` is an explicit, already-ordered list (never a set /
        dict view) so the ordering cannot drift between runs."""
        for t in node_tags:
            for d in range(1, ndf + 1):
                self.add("disp", t, d, ops.nodeDisp(t, d))

    def contact_forces(self, slave_tags):
        for t in slave_tags:
            self.add("contactForce", t, "-", ops.ladrunoContactForce(t))

    def scalar(self, kind, value):
        self.add(kind, "-", "-", value)


def _fresh(ndm=3, ndf=3):
    ops.wipe()
    ops.model("basic", "-ndm", ndm, "-ndf", ndf)


def _static(tol=1.0e-12, iters=40):
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", tol, iters, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


# --------------------------------------------------------------------------
# deck 1 -- RIGID_PLANE lane (analytical plane, ADR-39 P2a)
# --------------------------------------------------------------------------
KN = 1.0e6


def deck_rigid_plane_static(rec):
    """test_adr39_contact_p2a::test_p2a_axis_aligned_penetration_static, ramped
    over 5 LoadControl steps instead of 1.  Plane x=0 with normal +x; the slave
    starts just penetrated so the pair is active from step 1 (a slave held ONLY
    by the penalty is singular while separated).  The tiny z-anchor truss is the
    shipped idiom: it gives the model an element-backed FE_Element and has zero
    effect on the contact physics (no case has z-motion)."""
    P = 1.0e3
    _fresh()
    ops.node(1, -1.0e-8, 0.0, 0.0)
    ops.fix(1, 0, 1, 1)                       # free x only
    ops.node(9991, -1.0e-8, 0.0, 1.0)
    ops.fix(9991, 1, 1, 1)
    ops.uniaxialMaterial("Elastic", 9999, 1.0)
    ops.element("Truss", 9991, 1, 9991, 1.0e-9, 9999)
    ops.contactSurface(1, "-slave", 1)
    ops.contactPlane(1, 1, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, KN)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, -P, 0.0, 0.0)
    _static(tol=1.0e-12, iters=30)
    ops.integrator("LoadControl", 0.2)
    for _ in range(5):
        _step(ops.analyze(1))
    rec.disps([1, 9991], 3)


# --------------------------------------------------------------------------
# deck 2 -- SEGMENT (NTS) lane, explicit kn
# --------------------------------------------------------------------------
def deck_nts_penalty_static(rec):
    """test_adr39_contact_p2b::test_p2b_static_penetration_quad, ramped over 5
    LoadControl steps.  ONE fixed quad master facet, one slave node just
    penetrated, `-outward` orientation."""
    P = 1.0e3
    _fresh()
    mtags = []
    for i, (x, y) in enumerate([(-0.5, -0.5), (0.5, -0.5), (0.5, 0.5), (-0.5, 0.5)]):
        t = 101 + i
        ops.node(t, x, y, 0.0)
        ops.fix(t, 1, 1, 1)
        mtags.append(t)
    ops.node(1, 0.0, 0.0, -1.0e-8)
    ops.fix(1, 1, 1, 0)                       # free z only (normal is +z)
    ops.contactSurface(10, "-master", 4, *mtags)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -P)
    _static(tol=1.0e-12, iters=30)
    ops.integrator("LoadControl", 0.2)
    for _ in range(5):
        _step(ops.analyze(1))
    rec.disps([1] + mtags, 3)
    rec.contact_forces([1])


# --------------------------------------------------------------------------
# deck 3 -- SEGMENT (NTS) lane, `-kn auto`
# --------------------------------------------------------------------------
def deck_nts_autokn_static(rec):
    """test_adr39_contact_p2b2b::_run_autokn, ramped over 5 LoadControl steps.
    A DEFORMABLE LadrunoBrick supplies the owning solid stiffness `-kn auto`
    reduces (kn = f_si * mean(n^T K_block n), f_si = 0.10)."""
    L, E, nu, P = 1.0, 2.0e4, 0.0, 1.0e2
    _fresh()
    pts = [(0, 0, 0), (L, 0, 0), (L, L, 0), (0, L, 0),
           (0, 0, L), (L, 0, L), (L, L, L), (0, L, L)]
    tags = []
    for i, (x, y, z) in enumerate(pts):
        t = 1 + i
        ops.node(t, float(x), float(y), float(z))
        tags.append(t)
    for t in tags[:4]:
        ops.fix(t, 1, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    ops.element("LadrunoBrick", 1, *tags, 1)
    top = tags[4:]
    ops.node(99, L / 2.0, L / 2.0, L - 1.0e-8)
    ops.fix(99, 1, 1, 0)                      # free z
    ops.contactSurface(10, "-master", 4, *top)
    ops.contactSurface(20, "-slave", 99)
    ops.contact(1, 10, 20, "auto", "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(99, 0.0, 0.0, -P)
    _static(tol=1.0e-10, iters=50)
    ops.integrator("LoadControl", 0.2)
    for _ in range(5):
        _step(ops.analyze(1))
    rec.disps(tags + [99], 3)
    rec.contact_forces([99])


# --------------------------------------------------------------------------
# deck 4 -- SEGMENT (NTS) lane + Coulomb friction, EXPLICIT
# --------------------------------------------------------------------------
def deck_nts_friction_expl(rec):
    """test_adr39_contact_p3_friction::_driven_block, 20 central-difference
    steps at a fixed dt.  The slave is pre-penetrated to the normal equilibrium
    (z = -P/kn, so N = P with no vertical ringing) and driven by a constant +x
    force above the cone -- it slips, which exercises the return map, the
    committed friction state and the commit hook every step."""
    KT, mu, P, Q, m, dt = 1.0e5, 0.3, 1000.0, 600.0, 1.0, 1.0e-3
    _fresh()
    mtags = []
    for i, (x, y) in enumerate([(-10.0, -10.0), (10.0, -10.0), (10.0, 10.0), (-10.0, 10.0)]):
        t = 101 + i
        ops.node(t, x, y, 0.0)
        ops.fix(t, 1, 1, 1)
        mtags.append(t)
    ops.node(1, 0.0, 0.0, -P / KN)
    ops.mass(1, m, m, m)
    ops.contactSurface(10, "-master", 4, *mtags)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KT, mu, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, Q, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    for _ in range(20):
        _step(ops.analyze(1, dt))
    rec.disps([1] + mtags, 3)
    rec.contact_forces([1])


# --------------------------------------------------------------------------
# deck 5 -- MORTAR lane + commit-cycle Uzawa ALM
# --------------------------------------------------------------------------
def deck_mortar_alm(rec):
    """test_adr41_mortar_c2_2::_rigid_master_slave_facet: one slave facet on a
    coincident FIXED master facet under a uniform pressure, then a FIXED count
    of held-load Uzawa augmentations.

    The ADR-41 D1 bracket is what makes the held re-commits legal: between
    ladrunoBeginAugment and ladrunoEndAugment, Domain::commit() performs the
    contact lambda update but fires NO recorders and bumps NO commitTag.  The
    count is fixed (3) rather than tolerance-driven ON PURPOSE -- see the module
    docstring."""
    epsN, PTOT = 1.0e5, 400.0
    _fresh()
    for t, (x, y) in zip((1, 2, 3, 4), [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]):
        ops.node(t, x, y, 0.0)
        ops.fix(t, 1, 1, 1)
    for t, (x, y) in zip((5, 6, 7, 8), [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]):
        ops.node(t, x, y, -1.0e-4)            # seeded penetrating (c2 idiom)
        ops.fix(t, 1, 1, 0)                   # x,y fixed; z free
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    # coincident facets: the per-pair centroid datum degenerates, so -outward is
    # given explicitly (the shipped C2/C3/C4 idiom).
    ops.contact(1, 1, 2, "-mortar", "-epsN", epsN, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):                    # uniform pressure: a_I = 1/4 each, A = 1
        ops.load(t, 0.0, 0.0, -PTOT / 4.0)
    _static(tol=1.0e-11, iters=50)
    ops.integrator("LoadControl", 1.0)
    _step(ops.analyze(1))                     # penalty equilibrium (1 augmentation)
    rec.scalar("mortarPenPenalty", ops.ladrunoMortarPenetration())
    ops.integrator("LoadControl", 0.0)        # hold the load
    try:
        ops.ladrunoBeginAugment()
        for _ in range(3):
            _step(ops.analyze(1))
    finally:
        ops.ladrunoEndAugment()               # never leave recorders muted
    rec.disps([1, 2, 3, 4, 5, 6, 7, 8], 3)
    rec.scalar("mortarPenAugmented", ops.ladrunoMortarPenetration())


# --------------------------------------------------------------------------
# deck 6 -- MORTAR mesh-tying lane
# --------------------------------------------------------------------------
def deck_mortar_tie(rec):
    """test_adr41_mortar_c4_1::_split_column (tension, epsTie = 1e8), ramped
    over 5 LoadControl steps.  A column split at z=1 into a 1-brick bottom
    (master facet) and a 2-brick top (NON-MATCHING slave facets sharing an edge
    node at x=0.5) -- the shared-node global-accumulator path."""
    E, NU, PTOT, epsTie = 2.0e4, 0.0, 20.0, 1.0e8
    _fresh()
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    bot = [(1, 0, 0, 0), (2, 1, 0, 0), (3, 1, 1, 0), (4, 0, 1, 0),
           (5, 0, 0, 1), (6, 1, 0, 1), (7, 1, 1, 1), (8, 0, 1, 1)]
    for t, x, y, z in bot:
        ops.node(t, float(x), float(y), float(z))
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    top = [(9, 0, 0, 1), (10, 0.5, 0, 1), (11, 1, 0, 1),
           (12, 0, 1, 1), (13, 0.5, 1, 1), (14, 1, 1, 1),
           (15, 0, 0, 2), (16, 0.5, 0, 2), (17, 1, 0, 2),
           (18, 0, 1, 2), (19, 0.5, 1, 2), (20, 1, 1, 2)]
    for t, x, y, z in top:
        ops.node(t, float(x), float(y), float(z))
    ops.element("LadrunoBrick", 2, 9, 10, 13, 12, 15, 16, 19, 18, 1)
    ops.element("LadrunoBrick", 3, 10, 11, 14, 13, 16, 17, 20, 19, 1)
    for t in range(5, 21):
        ops.fix(t, 1, 1, 0)                   # uniaxial column
    ops.contactSurface(1, "-master", 4, 5, 6, 7, 8)
    ops.contactSurface(2, "-slave-segments", 4, 9, 10, 13, 12, 10, 11, 14, 13)
    ops.contact(1, 1, 2, "-mortar", "-tie", "-epsTie", epsTie,
                "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (15, 17, 18, 20):                # consistent nodal loads, uniform end traction
        ops.load(t, 0.0, 0.0, PTOT / 8.0)
    for t in (16, 19):
        ops.load(t, 0.0, 0.0, PTOT / 4.0)
    _static(tol=1.0e-12, iters=60)
    ops.integrator("LoadControl", 0.2)
    for _ in range(5):
        _step(ops.analyze(1))
    rec.disps(list(range(1, 21)), 3)
    rec.scalar("tieRes", ops.ladrunoMortarTieResidual())


# --------------------------------------------------------------------------
# deck 7 -- EDGE_EDGE lane (ADR-57)
# --------------------------------------------------------------------------
def deck_edge_edge(rec):
    """test_adr57_edge_edge_1::test_e2_crossed_bars_restrained, ramped over 5
    LoadControl steps, with `-edgeAlm` added (the E6 idiom) so the
    ladrunoEdgePenetration query is live -- it reports only over ALM pairs.

    Two perpendicular bars touching edge-on: no slave NODE projects onto the
    master FACE (the NTS blind spot) and the perpendicular facets have no clip
    overlap (the face-mortar blind spot).  Soft z-springs regularize the
    point-like pair's differential bottom mode."""
    EPS, KSPRING, PTOT, ZC = 1.0e6, 500.0, 400.0, -1.0e-4
    _fresh()
    for t, x, y, z in [(1, 0.0, -1.0, 0.0), (2, 0.0, 1.0, 0.0),
                       (3, 0.0, 1.0, -0.3), (4, 0.0, -1.0, -0.3)]:
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1, 1)
    for t, x, y, z in [(5, -1.0, 0.0, ZC), (6, 1.0, 0.0, ZC),
                       (7, 1.0, 0.0, ZC + 0.3), (8, -1.0, 0.0, ZC + 0.3)]:
        ops.node(t, x, y, z)
    for t in (7, 8):
        ops.fix(t, 1, 1, 1)                   # anchor the slave TOP edge
    for t in (5, 6):
        ops.fix(t, 1, 1, 0)                   # bottom edge: x,y fixed; z free
    ops.uniaxialMaterial("Elastic", 1, KSPRING)
    for gnd, sn, x in [(305, 5, -1.0), (306, 6, 1.0)]:
        ops.node(gnd, x, 0.0, ZC)
        ops.fix(gnd, 1, 1, 1)
        ops.element("zeroLength", 900 + sn, gnd, sn, "-mat", 1, "-dir", 3)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-edgeedge", "-edgeBand", 0.15,
                "-edgeAlm", "-edgeAugTol", 1.0e-9, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6):
        ops.load(t, 0.0, 0.0, -PTOT / 2.0)
    _static(tol=1.0e-11, iters=40)
    ops.integrator("LoadControl", 0.2)
    for _ in range(5):
        _step(ops.analyze(1))
    rec.disps([1, 2, 3, 4, 5, 6, 7, 8, 305, 306], 3)
    rec.scalar("edgePen", ops.ladrunoEdgePenetration())


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------
DECKS = [
    ("rigid_plane_static", deck_rigid_plane_static),
    ("nts_penalty_static", deck_nts_penalty_static),
    ("nts_autokn_static", deck_nts_autokn_static),
    ("nts_friction_expl", deck_nts_friction_expl),
    ("mortar_alm", deck_mortar_alm),
    ("mortar_tie", deck_mortar_tie),
    ("edge_edge", deck_edge_edge),
]


def run_all():
    """Run every deck.  Returns (lines, n_failed).  A deck failure never stops
    the sweep: a partial artifact is not diffable, a complete one with a
    recorded failure is."""
    lines = []
    n_failed = 0
    for name, fn in DECKS:
        sys.stderr.write("contact_dump: deck %s\n" % name)
        sys.stderr.flush()
        rec = Recorder()
        try:
            fn(rec)
        except DeckFailure as fail:
            lines.append("DECK %s FAILED %d" % (name, fail.code))
            n_failed += 1
            sys.stderr.write("contact_dump:   FAILED (analyze rc %d)\n" % fail.code)
            continue
        except Exception as exc:               # noqa: BLE001 - a raised refusal is a failure too
            lines.append("DECK %s FAILED -999" % name)
            n_failed += 1
            sys.stderr.write("contact_dump:   FAILED (%s: %s)\n"
                             % (type(exc).__name__, exc))
            continue
        for kind, tag, comp, value in rec.rows:
            lines.append("%s %s %s %s %s" % (name, kind, tag, comp, value.hex()))
    ops.wipe()
    return lines, n_failed


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="ADR-85 T0 contact byte-identity harness: run the shipped 3D "
                    "contact lanes and dump full-precision observables.")
    ap.add_argument("--out", required=True,
                    help="artifact path (overwritten; physics only, ASCII, LF)")
    args = ap.parse_args(argv)

    # provenance to STDERR ONLY -- it must never reach the artifact.
    sys.stderr.write("contact_dump: engine  = %s\n" % os.path.abspath(ops.__file__))
    sys.stderr.write("contact_dump: ladrunoBuild = %s\n" % ops.ladrunoBuild())
    sys.stderr.flush()

    lines, n_failed = run_all()

    # newline="\n" is load-bearing: the artifact is compared AS BYTES, and the
    # default Windows translation would make every line differ from a Linux run.
    with open(args.out, "w", encoding="ascii", newline="\n") as fh:
        for line in lines:
            fh.write(line + "\n")

    sys.stderr.write("contact_dump: wrote %d lines, %d deck(s) failed\n"
                     % (len(lines), n_failed))
    return 1 if n_failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
