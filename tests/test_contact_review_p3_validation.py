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


def _two_facets(h):
    """Matched slave-over-master quad facets of edge h, initially touching (zero gap)."""
    for t, (x, y) in zip((1, 2, 3, 4), [(0, 0), (h, 0), (h, h), (0, h)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)
    for t, (x, y) in zip((5, 6, 7, 8), [(0, 0), (h, 0), (h, h), (0, h)]):
        ops.node(t, float(x), float(y), 0.0)
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


def test_missing_node_contact_skipped_loudly():
    """A typo'd node tag no longer silently drops pairs: the whole contact is skipped
    with a warning, and the analysis still runs (no crash, no partial surface)."""
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
    assert ops.analyze(1, 1e-4) == 0        # skipped contact, healthy analysis


# ------------------------------------------------- 3. mortar large-facet (isomap) gate
def _mortar_press(h):
    """Slave facet pressed onto a matched master facet (edge h) under -mortar penalty +
    -visc settling; returns (settled mean penetration, the penalty prediction P_tot/(epsN*A))."""
    _fresh()
    _two_facets(h)
    A = h * h
    epsN = 1.0e-3
    Ptot = 0.1 * epsN * A                    # target penetration = 0.1 length units
    # Kelvin-Voigt settling (D2.2): the viscous scatter rides the mortar area weights,
    # so the per-node dashpot ~ muc*A/4 — size muc for zeta ~ 0.5 on the contact mode
    # (k_node ~ epsN*A/4, m = 1) so the press SETTLES at every scale.
    k_node = epsN * A / 4.0
    muc = math.sqrt(k_node) / (A / 4.0)
    ops.contact(10, 1, 2, "-mortar", "-epsN", epsN, "-visc", muc)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -Ptot / 4.0)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno"); ops.algorithm("Linear")
    ops.analysis("Transient")
    dt = 0.25 * 2.0 / math.sqrt(k_node)
    for _ in range(4000):
        assert ops.analyze(1, dt) == 0, "mortar press step failed"
    pen = -(sum(ops.nodeDisp(t, 3) for t in (5, 6, 7, 8)) / 4.0)
    return pen, Ptot / (epsN * A)


def test_mortar_contact_alive_at_large_facet_scale():
    """h=5e4 facets (mm-unit building scale): pre-fix the inverse isomap's absolute
    tolerance killed the Gauss points and the slave fell straight through; post-fix the
    settled penetration matches the penalty prediction at h=1 AND h=5e4."""
    pen1, ref1 = _mortar_press(1.0)
    assert abs(pen1 - ref1) < 0.25 * ref1, f"h=1 baseline off: pen={pen1} vs {ref1}"
    penL, refL = _mortar_press(5.0e4)
    assert abs(penL - refL) < 0.25 * refL, (
        f"mortar contact dead/wrong at h=5e4: pen={penL} vs prediction {refL} "
        "(the isomap absolute-tolerance regression is back)")
