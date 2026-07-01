"""ADR-63 #4a P1 — averaged nodal-normal smoothing: the convex-ridge no-flip gate.

A convex ridge (a tent of two quad facets meeting at a sharp angle) is the geometry that
breaks the shipped FACETED NTS normal (ADR-60 R3): the per-segment normal's outward SIGN is
chosen from the slave's REFERENCE side (auto `orientDir`). For a slave sitting on one facet
the auto direction is ~perpendicular to that facet's true normal, so the derived normal FLIPS
-> the gap sign flips -> the pair reads "not penetrating" -> SILENT pass-through (the slave is
driven straight through the master with no contact reaction).

`-smoothNormal` (ADR-63) replaces the per-segment normal with a smooth nodal-normal field whose
outward sense is a GLOBAL surface datum (the coherent winding + the slave-centroid / `-outward`
seed), captured once -> no per-slave sign to flip -> the repulsive normal is maintained.

  * `-smoothNormal` ON : the inward-pushed slave is RESISTED (bounded penetration along the facet
    normal).
  * default (faceted, auto orientDir) : the ridge-facet normal flips -> no contact -> the slave is
    driven through the master (penetration -> large). This companion case reproduces R3 and proves
    the gate above is meaningful.

A flat-master companion proves `-smoothNormal` does NOT perturb the flat case (D5 consistency:
on a flat master every facet normal is equal, so n_smooth == n_facet to round-off).

Build-free oracle: contact_prototypes/proto_nodal_normal_selfcheck.cpp (winding / nodal normals /
C0 continuity / refusal). Here it runs in compiled OpenSees.
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6
DELTA = 1.0e-3          # initial penetration of the slave below the right facet
INV_SQRT2 = 1.0 / math.sqrt(2.0)


def _build_tent():
    """A convex ridge: ridge along y at (x=0,z=1); left eave at x=-1,z=0; right eave at x=+1,z=0.
    Two quad facets ([L0,R0,R1,L1], [R0,Rt0,Rt1,R1]); all FIXED. Returns the flat master tags.
    The right facet's outward normal is (1,0,1)/sqrt2; the left facet's is (-1,0,1)/sqrt2."""
    #            tag        x     y     z
    ops.node(100, -1.0, -0.5, 0.0); ops.fix(100, 1, 1, 1)   # L0
    ops.node(101, -1.0,  0.5, 0.0); ops.fix(101, 1, 1, 1)   # L1
    ops.node(110,  0.0, -0.5, 1.0); ops.fix(110, 1, 1, 1)   # R0 (ridge)
    ops.node(111,  0.0,  0.5, 1.0); ops.fix(111, 1, 1, 1)   # R1 (ridge)
    ops.node(120,  1.0, -0.5, 0.0); ops.fix(120, 1, 1, 1)   # Rt0
    ops.node(121,  1.0,  0.5, 0.0); ops.fix(121, 1, 1, 1)   # Rt1
    left  = [100, 110, 111, 101]   # CCW: L0 -> R0 -> R1 -> L1
    right = [110, 120, 121, 111]   # CCW: R0 -> Rt0 -> Rt1 -> R1
    return left + right


def _right_facet_signed_dist(px, pz):
    """Signed distance of (px, y, pz) from the right-facet plane (through R0=(0,*,1),
    normal (1,0,1)/sqrt2): (P - R0).n = (px + pz - 1)/sqrt2. < 0 ⇒ penetrating."""
    return (px + pz - 1.0) * INV_SQRT2


def _press_into_right_facet(smooth):
    """Slave on the right facet at (0.5, 0, 0.5-DELTA) (penetrating by DELTA), mass 1, pushed
    INTO the tent along -(1,0,1)/sqrt2 by a constant load. Run explicit; return the min signed
    distance from the facet over the run (how far the slave is driven in).

    The slave's x,y DOFs are FIXED (only z free) so this gate isolates the R3 SIGN behaviour — the
    thing -smoothNormal fixes — from frictionless lateral motion. Rationale: a frictionless slave
    under a load with a lateral component on a CONVEX ridge has NO equilibrium (it slides up-slope and
    launches off), so a laterally-free rig can't cleanly test "held". Historically this gate ran
    laterally-free and only 'passed' because the P1 spurious adjacent-facet ejection happened to pin
    min_d ≥ 0 early (the bug masquerading as a hold); once P2.1 removes that ejection the free slave
    correctly slides off. Fixing x,y makes the well-posed statement: with the global-sign field the
    normal stays REPULSIVE (slave held, min_d≈-DELTA); with the faceted auto sign it FLIPS (slave
    driven straight down through the master, min_d≪-1). See ADR-63 P2.1."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mtags = _build_tent()
    x0, z0 = 0.5, 0.5 - DELTA
    ops.node(1, x0, 0.0, z0)
    ops.fix(1, 1, 1, 0)                 # x,y FIXED, z FREE — isolate the R3 sign (no frictionless slide)
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *mtags)
    ops.contactSurface(20, "-slave", 1)
    # NO -outward: the auto orientDir is the R3 trigger. -smoothNormal opts into the global-sign field.
    args = [1, 10, 20, KN, 0.0, 0.0]
    if smooth:
        args.append("-smoothNormal")
    ops.contact(*args)
    # constant load pushing the slave INTO the right facet (-x, -z), magnitude ~100
    f = 100.0 * INV_SQRT2
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, -f, 0.0, -f)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)
    min_d = _right_facet_signed_dist(x0, z0)
    for _ in range(3000):
        if ops.analyze(1, dt) != 0:
            break
        px = x0 + ops.nodeDisp(1, 1)
        pz = z0 + ops.nodeDisp(1, 3)
        d = _right_facet_signed_dist(px, pz)
        if d < min_d:
            min_d = d
    return min_d


def test_p1_smoothnormal_holds_the_ridge_facet():
    """-smoothNormal: the global-sign field keeps the ridge-facet normal repulsive, so the
    inward-pushed slave is held near the surface (bounded penetration)."""
    min_d = _press_into_right_facet(smooth=True)
    assert min_d > -0.05, f"slave was driven through despite -smoothNormal (min_d={min_d})"


def test_p1_default_faceted_flips_and_passes_through():
    """Default (faceted, auto orientDir) reproduces R3: the ridge-facet normal flips, the pair
    reads non-penetrating, and the inward load drives the slave through the master."""
    min_d = _press_into_right_facet(smooth=False)
    assert min_d < -1.0, f"expected R3 pass-through without -smoothNormal, but min_d={min_d}"


# --------------------------------------------------------- tri-3 chain (GAP-1 coverage)
# A FLAT square split into two triangles (sharing the diagonal). This exercises the full
# handler→engine→adapter→evalSegmentSmooth chain at nps=3 (the stride / shape-fn / install path that
# the quad tests never touch — review GAP-1) WITHOUT the sharp-ridge facet-ownership pathology (a
# slave over one triangle's interior projects cleanly out of bounds on the other). The R3 sign fix
# itself is covered by the quad convex-ridge gate above.
def _build_tri_patch():
    ops.node(300, 0.0, 0.0, 0.0); ops.fix(300, 1, 1, 1)
    ops.node(310, 1.0, 0.0, 0.0); ops.fix(310, 1, 1, 1)
    ops.node(311, 1.0, 1.0, 0.0); ops.fix(311, 1, 1, 1)
    ops.node(301, 0.0, 1.0, 0.0); ops.fix(301, 1, 1, 1)
    t0 = [300, 310, 311]   # lower-right tri (CCW from +z): contains (0.5,0.3)
    t1 = [300, 311, 301]   # upper-left tri  (shared diagonal 300-311, traversed opposite ⇒ coherent)
    return t0 + t1


def _press_tri(smooth):
    """Slave over the lower-right triangle's interior, pressed straight down; settled depth."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mtags = _build_tri_patch()
    gap0 = 1.0e-3
    ops.node(1, 0.6, 0.3, gap0); ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 3, *mtags)        # nps = 3
    ops.contactSurface(20, "-slave", 1)
    args = [1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0]
    if smooth:
        args.append("-smoothNormal")
    ops.contact(*args)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(1, 0.0, 0.0, -100.0)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno"); ops.algorithm("Linear"); ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)
    last_z = gap0
    for _ in range(2000):
        if ops.analyze(1, dt) != 0:
            break
        last_z = gap0 + ops.nodeDisp(1, 3)
    return last_z


def test_p1_tri3_chain_engages_and_holds():
    """tri-3 (nps=3) full chain: -smoothNormal builds + installs the field and the slave is held near
    the surface (bounded penetration) — the previously solver-unverified triangular plumbing."""
    z = _press_tri(smooth=True)
    assert z > -0.05, f"tri-3 smooth contact did not hold (z={z})"


def test_p1_tri3_consistency_smooth_equals_faceted():
    """On a FLAT tri patch n_smooth == n_facet, so -smoothNormal matches the faceted settled response."""
    z_faceted = _press_tri(smooth=False)
    z_smooth = _press_tri(smooth=True)
    assert abs(z_smooth - z_faceted) < 1.0e-9, (
        f"tri-3 -smoothNormal perturbed the flat case (faceted z={z_faceted}, smooth z={z_smooth})")


# --------------------------------------------------------- flat-master D5 consistency
NSEG = 4


def _build_flat_strip():
    for i in range(NSEG + 1):
        ops.node(1000 + i, float(i), -0.5, 0.0); ops.fix(1000 + i, 1, 1, 1)
        ops.node(2000 + i, float(i),  0.5, 0.0); ops.fix(2000 + i, 1, 1, 1)
    mtags = []
    for s in range(NSEG):
        mtags += [1000 + s, 1000 + s + 1, 2000 + s + 1, 2000 + s]
    return mtags


def _press_flat(smooth):
    """Slave pressed straight down onto a flat strip; return the settled penetration depth."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mtags = _build_flat_strip()
    gap0 = 1.0e-3
    ops.node(1, 1.5, 0.0, gap0)
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *mtags)
    ops.contactSurface(20, "-slave", 1)
    args = [1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0]
    if smooth:
        args.append("-smoothNormal")
    ops.contact(*args)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -100.0)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)
    last_z = gap0
    for _ in range(2000):
        if ops.analyze(1, dt) != 0:
            break
        last_z = gap0 + ops.nodeDisp(1, 3)
    return last_z


def test_p1_flat_master_consistency_smooth_equals_faceted():
    """On a FLAT master n_smooth == n_facet to round-off, so -smoothNormal must give the SAME
    settled response as the faceted path (D5 consistency — smoothing must not perturb the flat case)."""
    z_faceted = _press_flat(smooth=False)
    z_smooth  = _press_flat(smooth=True)
    assert abs(z_smooth - z_faceted) < 1.0e-9, (
        f"-smoothNormal perturbed the flat case (faceted z={z_faceted}, smooth z={z_smooth})")
