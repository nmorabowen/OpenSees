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


def _chain(handler, define_mortar=False, nsteps=40):
    """A small 3D truss chain (SP constraints only), optionally with an inert mortar
    contact defined on top — used to prove the mortar definition is graph-neutral."""
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

    ops.constraints(handler)
    if define_mortar:
        # faceted surfaces over existing nodes (geometry irrelevant — C2.0 never reads it);
        # the mortar contact lands in theMortarContacts, which the handler does not iterate.
        ops.contactSurface(1, "-master", 3, 1, 2, 3)
        ops.contactSurface(2, "-slave-segments", 3, 3, 4, 5)
        ops.contact(1, 1, 2, "-mortar", "-epsN", "auto",
                    "-augTol", 1e-8, "-maxAug", 20, "-ngp", 2)

    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.test("NormDispIncr", 1e-12, 25, 0)
    ops.algorithm("Newton")
    ops.analysis("Transient")

    out = []
    for _ in range(nsteps):
        assert ops.analyze(1, 1.0e-3) == 0
        out.append(tuple(ops.nodeDisp(i + 1, 1) for i in range(n)))
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


def test_c2_0_mortar_contact_is_inert_byte_identical():
    """A model with a mortar contact DEFINED is BITWISE identical to no-contact — the
    definition adds no DOF-graph edges (the handler never reads theMortarContacts in C2.0)."""
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
