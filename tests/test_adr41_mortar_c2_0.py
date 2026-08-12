"""ADR-41 C2.0 — mortar command surface: SLAVE_SEGMENTS + the `mortar` contact selector.

C2.0 is the inert PLUMBING phase (the P1a "zero adapter" philosophy): the faceted slave
surface (`contactSurface -slave-segments`) and the mortar contact definition
(`contact ... -mortar -epsN -augTol -maxAug -ngp`) are PARSED, VALIDATED, and STORED on the
Domain — but in a SEPARATE list the NTS handler loop never reads, so a defined mortar contact
is byte-identical/inert until the narrow-phase adapter is wired in C2.1. This battery gates:
  - the definition round-trips (ladrunoContactInfo() reports the mortar count, the 4th element);
  - a model WITH a mortar contact defined is BITWISE identical to no-contact (graph-neutral);
  - the parser/Domain guards reject malformed input (nps<3, wrong surface kinds).

See Ladruno_implementation/_adr41_c2_design.md (C2.0) + 48_ladruno_contact_capstone_adr.md.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


MASTER_FACET = ((101, 0.0, 0.0, 10.0), (102, 1.0, 0.0, 10.0), (103, 0.0, 1.0, 10.0))
# 10 length units above the master: a real positive gap. Drop this to 9.5 and the
# facets interpenetrate, the penalty engages, and the byte-identity below breaks --
# that is the mutation which proves this test can fail.
SLAVE_FACET = ((201, 0.0, 0.0, 20.0), (202, 1.0, 0.0, 20.0), (203, 0.0, 1.0, 20.0))


def _chain(handler, define_mortar=False, nsteps=40):
    """A small 3D truss chain (SP constraints only) plus a free-floating pair of
    facets, optionally with a mortar contact DECLARED on those facets — used to prove
    the declaration is graph-neutral.

    The facets are built in BOTH runs and only the contactSurface/contact declarations
    are conditional, so the comparison isolates the declaration itself rather than the
    presence of six extra nodes. Their z DOFs are FREE and carry mass, and their
    displacements are part of the compared tuple: that is what gives the test teeth. An
    earlier repair pinned these nodes fully fixed and compared only the truss, which
    made the assertion unfalsifiable — no mortar behaviour, correct or broken, could
    have moved the compared quantity. Caught by mutation, not by review.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    n = 5
    for i in range(n):
        ops.node(i + 1, float(i), 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    for i in range(2, n + 1):
        ops.fix(i, 0, 1, 1)
        ops.mass(i, 1.0, 1.0, 1.0)
    ops.uniaxialMaterial("Elastic", 1, 1000.0)
    for i in range(1, n):
        ops.element("Truss", i, i, i + 1, 1.0, 1)
    ops.setNodeVel(n, 1, 0.5, "-commit")

    # the facet nodes: x/y held, z free + massive (so the transient tangent
    # M/(beta dt^2) is nonsingular without needing springs), slave drifting gently
    # downward so the compared z displacements are non-trivial and demonstrably live.
    for t, x, y, z in MASTER_FACET + SLAVE_FACET:
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1, 0)
        ops.mass(t, 1.0, 1.0, 1.0)
    for t, _x, _y, _z in SLAVE_FACET:
        ops.setNodeVel(t, 3, -0.25, "-commit")

    ops.constraints(handler)
    if define_mortar:
        # A legitimate, fully-specified mortar contact that contributes nothing because
        # the two surfaces are 10 length units apart.
        #
        # This deck used to put the facets on the truss nodes 1..5 with `-epsN auto`,
        # and was green for the WRONG REASON: collinear nodes give degenerate facets,
        # auto-sizing failed, and the contact was inert because it had been silently
        # DROPPED -- the exact degradation ADR-78 P1 abolished. The test was measuring
        # the bug that made it pass. Two things changed: the penalty is explicit, so
        # there is nothing left to fail to size; and the inertness now rests on a
        # durable physical reason (a positive gap) instead of on degenerate geometry
        # that a future zero-area guard would rightly reject.
        ops.contactSurface(1, "-master", 3, *[t for t, *_ in MASTER_FACET])
        ops.contactSurface(2, "-slave-segments", 3, *[t for t, *_ in SLAVE_FACET])
        ops.contact(1, 1, 2, "-mortar", "-epsN", 1.0e6, "-outward", 0.0, 0.0, 1.0,
                    "-augTol", 1e-8, "-maxAug", 20, "-ngp", 2)

    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.test("NormDispIncr", 1e-12, 25, 0)
    ops.algorithm("Newton")
    ops.analysis("Transient")

    out = []
    facets = [t for t, *_ in MASTER_FACET + SLAVE_FACET]
    for _ in range(nsteps):
        assert ops.analyze(1, 1.0e-3) == 0
        # BOTH the truss x-displacements and the facet z-displacements. The facet
        # terms are what make this falsifiable: an engaged mortar contact moves those
        # nodes and nothing else, so omitting them leaves an assertion that no mortar
        # behaviour could ever break.
        out.append(tuple(ops.nodeDisp(i + 1, 1) for i in range(n))
                   + tuple(ops.nodeDisp(t, 3) for t in facets))
    return out


def test_c2_0_mortar_definition_round_trips():
    """A slave-segments surface + a mortar contact are stored; ladrunoContactInfo()
    reports the mortar count in the 4th slot (NTS count stays 0)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i in range(1, 6):
        ops.node(i, float(i), 0.0, 0.0)
    ops.contactSurface(1, "-master", 3, 1, 2, 3)
    ops.contactSurface(2, "-slave-segments", 3, 3, 4, 5)
    ops.contact(1, 1, 2, "-mortar", "-epsN", 1.0e6, "-maxAug", 12, "-ngp", 3)  # raises on fail
    info = ops.ladrunoContactInfo()      # [numContacts, numCommits, numReverts, numMortar]
    assert len(info) == 4
    assert info[0] == 0, "a mortar contact must NOT count as an NTS contact"
    assert info[3] == 1, "the mortar contact was not stored"
    ops.wipe()
    assert ops.ladrunoContactInfo() == [0, 0, 0, 0], "engine leaked across wipe"


def test_c2_0_mortar_out_of_contact_byte_identical():
    """A model with an OUT-OF-CONTACT mortar contact defined is BITWISE identical to
    no-contact: declaring surfaces and a mortar pair adds no DOF-graph edges and
    perturbs no result.

    Renamed from `..._mortar_contact_is_inert_byte_identical`. That name dates from
    C2.0, when the mortar lane was inert PLUMBING by construction; the lane has been
    live since C2.1 (ADR-41), so "is inert" no longer describes anything true about
    the code -- only about this particular deck, whose surfaces do not touch. The
    graph-neutrality property the test actually gates survives the change intact.
    See _chain() for why the old fixture was green for the wrong reason.
    """
    ref = _chain("Plain", define_mortar=False)
    mor = _chain("LadrunoContact", define_mortar=True)
    assert mor == ref


def test_c2_0_slave_segments_requires_nps_ge_3():
    """contactSurface -slave-segments needs nodesPerSeg >= 3 (a facet)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i in range(1, 5):
        ops.node(i, float(i), 0.0, 0.0)
    with pytest.raises(Exception):
        ops.contactSurface(9, "-slave-segments", 2, 1, 2)   # nps=2 < 3


def test_c2_0_mortar_requires_faceted_surfaces():
    """A mortar contact needs a MASTER_SEGMENTS master and a SLAVE_SEGMENTS slave —
    a node-set (-slave) slave carries no facet, so it cannot give the mortar matrix D."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i in range(1, 7):
        ops.node(i, float(i), 0.0, 0.0)
    ops.contactSurface(1, "-master", 3, 1, 2, 3)
    ops.contactSurface(2, "-slave", 4, 5, 6)            # SLAVE_NODES (wrong for mortar)
    ops.contactSurface(3, "-slave-segments", 3, 4, 5, 6)
    with pytest.raises(Exception):
        ops.contact(1, 1, 2, "-mortar")                 # slave is a node-set -> reject
    with pytest.raises(Exception):
        ops.contact(2, 3, 3, "-mortar")                 # master is slave-segments -> reject
    # the correct kinds are accepted (a successful openseespy command returns None, not 0)
    ops.contact(3, 1, 3, "-mortar")
    assert ops.ladrunoContactInfo()[3] == 1
