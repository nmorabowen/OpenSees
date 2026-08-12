"""Contact review fix P3 (2026-07) — friction-cap semantics + choke-point validation.

1. τmax-only refusal: `-tauMax` without `-mu`/`-cohesion` is the unified cone
   min(0, τmax) = 0 (free slip) — the OLD kernel silently turned it into an UNBOUNDED
   elastic bond; the command surface now refuses it loudly (the kernel itself now
   free-slips as defense-in-depth: contact_prototypes/proto_friction_validation.cpp,
   6/12 checks fail on the pre-fix kernels).
2. Choke-point validation (all previously silent): duplicate contact tags (the tag is
   the leading key of every path-state store — a duplicate ALIASES friction/re-emit
   state across physical contacts), kt < 0, a flat connectivity that is not a whole
   number of segments (the trailing partial segment was silently dropped), a missing
   node / a 2D node in a contact surface (previously skipped in silence; a 2D node
   read getCrds()(2) out of bounds).
3. Mortar large-facet gate: inverseIsomap2D converged on an ABSOLUTE 1e-13 residual
   with LENGTH units (aux-plane UV ~ h), so facets ≳ 2000 units silently lost their
   Gauss points (378/400 dead at h=5e4) — mortar contact quietly vanished at exactly
   the scales the PR-1 projection fix rescued for NTS. The settled-penetration test
   below FAILS on the pre-fix build (free fall through the master facet).
"""
import math
import pytest

from _testbed import ops


pytestmark = [pytest.mark.zone_a]


def _two_facets(h, z0=0.0):
    """Matched slave-over-master quad facets of edge h; slave plane at z0 (seed z0 < 0 =
    slightly penetrating, the shipped c2_1 idiom — a slave held ONLY by the penalty is
    singular/inert while separated)."""
    for t, (x, y) in zip((1, 2, 3, 4), [(0, 0), (h, 0), (h, h), (0, h)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)
    for t, (x, y) in zip((5, 6, 7, 8), [(0, 0), (h, 0), (h, h), (0, h)]):
        ops.node(t, float(x), float(y), z0)
        ops.mass(t, 1.0, 1.0, 1.0)
        ops.fix(t, 1, 1, 0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)


def _fresh():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)


# ------------------------------------------------------------------ 1. τmax-only refusals
def test_taumax_only_refused():
    _fresh()
    _two_facets(1.0)
    with pytest.raises(Exception):
        ops.contact(10, 1, 2, "-mortar", "-epsN", 1e6, "-tauMax", 1e6)


def test_taumax_with_mu_or_cohesion_accepted():
    _fresh()
    _two_facets(1.0)
    ops.contact(10, 1, 2, "-mortar", "-epsN", 1e6, "-mu", 0.3, "-tauMax", 1e6)
    _fresh()
    _two_facets(1.0)
    ops.contact(10, 1, 2, "-mortar", "-epsN", 1e6, "-cohesion", 100.0, "-tauMax", 1e6)


def test_nts_taumax_refused_as_mortar_only():
    """-tauMax on the NTS lane was silently INERT (the positional-mu NTS cone has no
    shear cap); now refused outright as a mortar-only option (fail-loud, gate MINOR-1)."""
    _fresh()
    ops.node(1, -1.0, -1.0, 0.0); ops.node(2, 1.0, -1.0, 0.0)
    ops.node(3, 1.0, 1.0, 0.0); ops.node(4, -1.0, 1.0, 0.0)
    ops.node(5, 0.0, 0.0, 1e-4)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave", 5)
    with pytest.raises(Exception):
        ops.contact(7, 1, 2, 1e6, 1e5, 0.4, "-outward", 0.0, 0.0, 1.0, "-tauMax", 50.0)


def test_edge_taumax_only_refused():
    _fresh()
    _two_facets(1.0)
    with pytest.raises(Exception):
        ops.contact(10, 1, 2, "-mortar", "-epsN", 1e6,
                    "-edgeedge", "-edgeKn", 1e6, "-edgeTauMax", 1e6)


# ------------------------------------------------------------- 2. choke-point validation
def test_duplicate_contact_tag_refused():
    _fresh()
    g0 = 1.0e-4
    ops.node(1, -1.0, -1.0, 0.0); ops.fix(1, 1, 1, 1)
    ops.node(2,  1.0, -1.0, 0.0); ops.fix(2, 1, 1, 1)
    ops.node(3,  1.0,  1.0, 0.0); ops.fix(3, 1, 1, 1)
    ops.node(4, -1.0,  1.0, 0.0); ops.fix(4, 1, 1, 1)
    ops.node(5, 0.0, 0.0, g0); ops.mass(5, 1.0, 1.0, 1.0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave", 5)
    ops.contact(7, 1, 2, 1e6, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    with pytest.raises(Exception):          # same tag, NTS lane
        ops.contact(7, 1, 2, 1e6, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    with pytest.raises(Exception):          # same tag, rigid-plane lane
        ops.contactPlane(7, 2, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1e6)
    ops.contact(8, 1, 2, 1e6, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)   # fresh tag fine


def test_negative_kt_refused():
    _fresh()
    ops.node(1, -1.0, -1.0, 0.0); ops.node(2, 1.0, -1.0, 0.0)
    ops.node(3, 1.0, 1.0, 0.0); ops.node(4, -1.0, 1.0, 0.0)
    ops.node(5, 0.0, 0.0, 1e-4)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave", 5)
    with pytest.raises(Exception):
        ops.contact(7, 1, 2, 1e6, -1.0, 0.4, "-outward", 0.0, 0.0, 1.0)


def test_partial_segment_connectivity_refused():
    _fresh()
    for t in range(1, 8):
        ops.node(t, float(t), 0.0, 0.0)
    with pytest.raises(Exception):          # 7 tags is not a whole number of quads
        ops.contactSurface(1, "-master", 4, 1, 2, 3, 4, 5, 6, 7)


def test_missing_node_contact_aborts():
    """A typo'd node tag is FATAL -- the analysis refuses to run.

    This test has always gated one thing: a bad node tag must never be quietly
    tolerated. What "not tolerated" means has been ratcheted twice, and the test
    tracks the ratchet rather than pinning one rung of it:

        drop the affected pairs, silently   (pre contact-review P3)
      -> skip the whole contact, warn       (contact-review P3 -- what this used
                                             to assert, hence the old name)
      -> abort                              (ADR-78 P1 -- asserted here)

    P1's rule is that a DECLARED contact never silently fails to happen, and a
    skipped-but-warned contact is precisely that failure under MPI, where the
    warning is one line in one rank's log among many.

    Inverted rather than repaired: fixing the typo would have deleted the gate
    outright, since the missing node IS the subject. Contrast
    test_c2_0_mortar_out_of_contact_byte_identical in the ADR-41 battery, whose
    subject was byte-identity and whose broken precondition was incidental.
    """
    _fresh()
    ops.node(1, -1.0, -1.0, 0.0); ops.fix(1, 1, 1, 1)
    ops.node(2,  1.0, -1.0, 0.0); ops.fix(2, 1, 1, 1)
    ops.node(3,  1.0,  1.0, 0.0); ops.fix(3, 1, 1, 1)
    ops.node(4, -1.0,  1.0, 0.0); ops.fix(4, 1, 1, 1)
    ops.node(5, 0.0, 0.0, 1e-4); ops.mass(5, 1.0, 1.0, 1.0)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 999)   # 999 does not exist
    ops.contactSurface(2, "-slave", 5)
    ops.contact(7, 1, 2, 1e6, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno"); ops.algorithm("Linear")
    ops.analysis("Transient")
    # handle() FATALs -> domainChanged() fails -> analyze() reports -1. Asserting
    # the exact code, not merely != 0: a bare "it failed" would also pass if the
    # deck stopped converging for an unrelated reason.
    assert ops.analyze(1, 1e-4) == -1


# ------------------------------------------------- 3. mortar large-facet (isomap) gate
def _mortar_press(h):
    """The shipped c2_1 STATIC rig at facet scale h: slave quad seeded slightly
    penetrating a coincident fixed master quad, uniform load, LoadControl + Newton.
    Returns (static penetration, the penalty prediction P_tot/(epsN*A))."""
    _fresh()
    z0 = -1.0e-4 * max(1.0, h)               # seeded slightly penetrating (c2_1 idiom)
    _two_facets(h, z0)
    A = h * h
    epsN = 1.0e-3
    Ptot = 0.1 * epsN * A                    # target penetration = 0.1 length units
    ops.contact(10, 1, 2, "-mortar", "-epsN", epsN,
                "-outward", 0.0, 0.0, 1.0)   # coincident facets: explicit orientation
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -Ptot / 4.0)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-10 * max(1.0, h), 40, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "static mortar press did not converge"
    pen = -(z0 + sum(ops.nodeDisp(t, 3) for t in (5, 6, 7, 8)) / 4.0)
    return pen, Ptot / (epsN * A)


def test_mortar_contact_alive_at_large_facet_scale():
    """h=5e4 facets (mm-unit building scale): pre-fix the inverse isomap's absolute
    tolerance killed the Gauss points (contact stiffness vanished — Newton fails or the
    slave shoots through); post-fix the static penetration equals the penalty prediction
    δ = P/(epsN·A) at h=1 AND h=5e4."""
    pen1, ref1 = _mortar_press(1.0)
    assert abs(pen1 - ref1) < 0.05 * ref1, f"h=1 baseline off: pen={pen1} vs {ref1}"
    penL, refL = _mortar_press(5.0e4)
    assert abs(penL - refL) < 0.05 * refL, (
        f"mortar contact dead/wrong at h=5e4: pen={penL} vs prediction {refL} "
        "(the isomap absolute-tolerance regression is back)")
