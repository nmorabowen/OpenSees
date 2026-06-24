"""ADR-41 D2.2 follow-up — viscous stabilization on the SOFT=2 segment-based explicit penalty.

B2 (#406) shipped the SOFT=2 segment-based explicit penalty (`contact … -mortar -soft <SOFSCL>`) as a
FRICTIONLESS, single-pass MVP that refused `-visc`. This adds the D2.2 velocity-proportional NORMAL
damper to the SOFT=2 active set: per IN-CONTACT slave node (the soft-penalty KKT mask p_I<0) the
weighted normal gap RATE drives `p_visc_I = μ_c·ḡ̇_I`, scattered via D/−M·n exactly like the normal
penalty — the SAME oracle-validated D2.2 operator (proto_d2_viscous.py T6), now on the SOFT=2 path. It
bleeds contact chatter / controls restitution for explicit impact at the structural dt (the
progressive-collapse / element-removal regime SOFT=2 targets, where damping the recontact bounce matters).

Gates:
  1. EXPLICIT RESTITUTION — a SOFT=2 impact with `-visc` rebounds with restitution e < 1 (energy
     dissipated), monotone-decreasing in μ_c; `-visc 0` (absent) rebounds e ≈ 1 (the B2 conservative
     penalty). The damper is genuinely active on the SOFT=2 path.
  2. STATIC BYTE-IDENTITY — under an implicit/static solve `-visc` is inert (v≡0 + StaticIntegrator
     never forms the C tangent; SOFT=2 falls through to the regular mortar penalty), so adding `-visc`
     does not perturb a static solution.
  3. COMMAND — `-mortar -soft -visc` is accepted (covered also by test_adr39_contact_p5_soft2).

The D2.2 viscous operator was validated oracle-first in contact_prototypes/proto_d2_viscous.py (T6);
SOFT=2 reuses it unchanged on the soft-penalty active set.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

M = 1.0
DT = 1.0e-3
KN = 1.0e6


def _fixed_master():
    for t, (x, y) in zip((1, 2, 3, 4), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), 0.0)
        ops.fix(t, 1, 1, 1)
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)


def _matched_slave(z):
    for t, (x, y) in zip((5, 6, 7, 8), [(0, 0), (1, 0), (1, 1), (0, 1)]):
        ops.node(t, float(x), float(y), z)
        ops.mass(t, M, M, M)
        ops.fix(t, 1, 1, 0)


def _soft2_impact(muc, soft=0.1, v_mag=2.0, gap0=0.02, nsteps=3000):
    """A free matched slave facet gap0 above the master, flung DOWN at v_mag, SOFT=2 explicit. With a
    viscous damper (μ_c>0) the rebound loses energy ⇒ restitution e<1. Returns |v_out/v_mag|."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _fixed_master()
    _matched_slave(gap0)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    if muc > 0.0:
        ops.contact(1, 1, 2, "-mortar", "-epsN", KN, "-soft", soft, "-visc", muc,
                    "-outward", 0.0, 0.0, 1.0)
    else:
        ops.contact(1, 1, 2, "-mortar", "-epsN", KN, "-soft", soft, "-outward", 0.0, 0.0, 1.0)
    for t in (5, 6, 7, 8):
        ops.setNodeVel(t, 3, -v_mag, "-commit")
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    v_out = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, DT) != 0:
            break
        z = gap0 + ops.nodeDisp(5, 3)
        if z > 0.5 * gap0 and ops.nodeVel(5, 3) > 0.0:
            v_out = ops.nodeVel(5, 3)
    return abs(v_out / v_mag)


def test_soft2_visc_dissipates_restitution():
    """GATE 1 — the headline. A SOFT=2 explicit impact with the viscous damper rebounds with reduced
    restitution, and e DECREASES monotonically with μ_c (energy bled through the dashpot); the
    damper-free SOFT=2 impact rebounds near-conservatively. Same geometry + dt — only μ_c differs,
    proving the viscous term is genuinely active on the SOFT=2 active set."""
    e0 = _soft2_impact(0.0)         # no viscous: near-conservative (the B2 penalty, e ≈ 1)
    e1 = _soft2_impact(200.0)       # moderate damping
    e2 = _soft2_impact(500.0)       # heavier damping
    assert e0 > 0.95, f"SOFT=2 (no -visc) should be near-conservative (e={e0:.4f})"
    assert e1 < e0 - 0.1, f"-visc must dissipate on the SOFT=2 path (e={e1:.4f} !< e0={e0:.4f})"
    assert e2 < e1 - 0.1, f"restitution must decrease with μ_c (e(500)={e2:.4f} !< e(200)={e1:.4f})"
    assert e2 < 0.8, f"heavy -visc should clearly dissipate the rebound (e={e2:.4f})"


def _soft2_static(muc, soft=0.1, P=100.0):
    """The matched slave facet pressed onto the fixed master by a STATIC load (implicit Newton). SOFT=2
    is explicit-only ⇒ falls through to the regular mortar penalty; -visc is inert in statics (v≡0,
    StaticIntegrator forms no C tangent). Returns node-5 penetration."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _fixed_master()
    _matched_slave(0.0)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    if muc > 0.0:
        ops.contact(1, 1, 2, "-mortar", "-epsN", KN, "-soft", soft, "-visc", muc,
                    "-outward", 0.0, 0.0, 1.0)
    else:
        ops.contact(1, 1, 2, "-mortar", "-epsN", KN, "-soft", soft, "-outward", 0.0, 0.0, 1.0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6, 7, 8):
        ops.load(t, 0.0, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ops.integrator("LoadControl", 1.0)
    assert ops.analyze(1) == 0, "static SOFT=2 contact did not converge"
    return ops.nodeDisp(5, 3)


def _big_slave(z):
    """4 corner nodes of a [-0.5,1.5]² slave facet (all OUTSIDE master [0,1]²) at height z — the
    partial-overlap geometry NTS misses but SOFT=2 catches by integrating the clipped overlap."""
    for t, (x, y) in zip((5, 6, 7, 8), [(-0.5, -0.5), (1.5, -0.5), (1.5, 1.5), (-0.5, 1.5)]):
        ops.node(t, float(x), float(y), z)
        ops.mass(t, M, M, M)
        ops.fix(t, 1, 1, 0)


def _soft2_edge_impact(muc, soft=0.1, v_mag=1.0, z0=-0.005, nsteps=2000):
    """The big-slave EDGE-OVERLAP geometry (no slave node in-bounds — only the clipped interior
    overlaps the master) flung DOWN, SOFT=2 explicit. Returns the max rebound height of node 5 above
    the deepest point (a damped contact rebounds less). Exercises viscous on a PARTIAL clip."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _fixed_master()
    _big_slave(z0)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    if muc > 0.0:
        ops.contact(1, 1, 2, "-mortar", "-epsN", KN, "-soft", soft, "-visc", muc,
                    "-outward", 0.0, 0.0, 1.0)
    else:
        ops.contact(1, 1, 2, "-mortar", "-epsN", KN, "-soft", soft, "-outward", 0.0, 0.0, 1.0)
    for t in (5, 6, 7, 8):
        ops.setNodeVel(t, 3, -v_mag, "-commit")
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    min_z, max_after = z0, -1e9
    for _ in range(nsteps):
        if ops.analyze(1, DT) != 0:
            break
        z = z0 + ops.nodeDisp(5, 3)
        min_z = min(min_z, z)
        if z > min_z:
            max_after = max(max_after, z)
    return min_z, max_after


def test_soft2_visc_on_partial_clip_edge_overlap():
    """The SOFT=2 headline geometry (edge overlap, NO slave node in-bounds) WITH viscous. The viscous
    mask keys on the soft-penalty active set p_I<0, so on a PARTIAL clip only the clipped-in nodes are
    damped — a distinct code state from the matched facet. Both runs stay restrained (SOFT=2 catches
    the overlap); the heavily-damped run rebounds LESS (energy bled on the partial overlap)."""
    min_off, reb_off = _soft2_edge_impact(0.0)
    min_on, reb_on = _soft2_edge_impact(500.0)
    # SOFT=2 restrains the edge overlap in BOTH cases (no fall-through — the NTS-missed contact is caught)
    assert min_off > -0.1, f"SOFT=2 (no -visc) should restrain the edge overlap (min z={min_off:.3e})"
    assert min_on > -0.1, f"SOFT=2 (+ -visc) should restrain the edge overlap (min z={min_on:.3e})"
    # the damped run dissipates ⇒ a lower rebound than the conservative one (viscous active on the clip)
    assert reb_on < reb_off - 1e-3, (
        f"-visc must damp the partial-clip rebound (reb_on={reb_on:.3e} !< reb_off={reb_off:.3e})")


def test_soft2_visc_static_byte_identical():
    """GATE 2 — `-visc` is inert in statics: the velocity is identically zero and the StaticIntegrator
    never forms the C (damping) tangent, so a static SOFT=2 solve WITH `-visc` is BIT-IDENTICAL to one
    without it (the D2 statics-inert contract, on the SOFT=2 path)."""
    z_off = _soft2_static(0.0)
    z_on = _soft2_static(20.0)
    assert z_off == z_on, f"-visc perturbed the static (v≡0) SOFT=2 solution: {z_off!r} vs {z_on!r}"
    assert z_off < 0.0, f"the static SOFT=2 contact should penetrate under load (z={z_off:.3e})"
