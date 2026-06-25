"""ADR-57 E2 — perpendicular edge-edge penalty contact (the cos_t→0 narrow phase).

The `EDGE_EDGE` LadrunoContactFE mode (a 4th mode of the one contact adapter) runs the
LadrunoEdgeKernel segment-segment closest point at the trial config → the common-
perpendicular normal n=(ê_s×ê_m)/‖·‖ + the signed gap gN with a BODY-FIXED committed sign,
and assembles the penalty force f = tN·B (tN = εN⟨−gN⟩, B = [(1−s)n|s n|−(1−t)n|−t n]) +
the main tangent K = εN·BᵀB. The handler routes a -mortar -edgeedge facet pair to it when
the facets are perpendicular (align_fn→0, the clip degenerates) — exactly where face-mortar
AND NTS both go blind. The mechanics were validated oracle-first in
contact_prototypes/proto_e2_penalty.py (the geometry in proto_e0/proto_e1).

Two crossed bars, edge-on: a slave bar in the x-z plane whose BOTTOM edge runs along x,
crossing perpendicular over a master bar in the y-z plane whose TOP edge runs along y.
Only the two crossing edges touch — no slave NODE projects onto the master FACE interior
(the node-to-segment lane's blind spot), and the perpendicular facets have no clip overlap
(the face-mortar lane's blind spot). Soft springs hold each bottom node (the point-like edge
pair shares ONE gap ⇒ the differential bottom mode needs regularizing) and double as the NTS
baseline.

Gates:
  - the slave bottom edge, pressed down, settles at the penalty penetration δ ≈ P/εN
    (the springs barely shift it) and static Newton CONVERGES (the K_c-correctness check);
  - THE HEADLINE FALSIFIER: the identical model with an NTS contact (slave NODES vs the
    master facet) transmits ZERO contact force — no slave node penetrates the master face —
    so the bottom edge sinks to the pure-spring depth δ_nts = P/(2k), ~10³× deeper. NTS
    passes through where edge-edge restrains;
  - byte-identity: a -mortar -edgeedge model whose bars are far apart (no pair routes
    edge-edge) evolves BITWISE identically to the same model without -edgeedge.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

EPS = 1.0e6      # edge-edge normal penalty εN
KSPRING = 500.0  # soft regularizing spring (also the NTS pure-spring baseline)
PTOT = 400.0     # total downward load on the slave bottom edge
ZC = -1.0e-4     # slave bottom edge seeded slightly penetrating (contact active from step 1)


def _build_bars(xoff=0.0):
    """Two perpendicular bars. Master (fixed) in the y-z plane at x=0, its TOP edge (m1-m2)
    along y at z=0. Slave in the x-z plane at y=0 (shifted +xoff for the far-apart case), its
    BOTTOM edge (s1-s2) along x at z=ZC. Returns (master tags, slave tags)."""
    # master bar: quad in y-z plane, top edge m1-m2 at z=0
    m = [(1, 0.0, -1.0, 0.0), (2, 0.0, 1.0, 0.0), (3, 0.0, 1.0, -0.3), (4, 0.0, -1.0, -0.3)]
    for t, x, y, z in m:
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1, 1)
    # slave bar: quad in x-z plane, bottom edge s1-s2 at z=ZC
    s = [(5, -1.0 + xoff, 0.0, ZC), (6, 1.0 + xoff, 0.0, ZC),
         (7, 1.0 + xoff, 0.0, ZC + 0.3), (8, -1.0 + xoff, 0.0, ZC + 0.3)]
    for t, x, y, z in s:
        ops.node(t, x, y, z)
    return [1, 2, 3, 4], [5, 6, 7, 8]


def _add_springs():
    """Soft z-springs on the bottom slave nodes 5,6 (regularize the point-pair differential
    mode AND serve as the NTS pure-spring baseline). Ground nodes 305/306 fixed, coincident."""
    ops.uniaxialMaterial("Elastic", 1, KSPRING)
    for gnd, sn, x in [(305, 5, -1.0), (306, 6, 1.0)]:
        ops.node(gnd, x, 0.0, ZC)
        ops.fix(gnd, 1, 1, 1)
        ops.element("zeroLength", 900 + sn, gnd, sn, "-mat", 1, "-dir", 3)


def _solve_static():
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1.0e-11, 40, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    return ops.analyze(1)


def _load_bottom():
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t in (5, 6):
        ops.load(t, 0.0, 0.0, -PTOT / 2.0)


def test_e2_crossed_bars_restrained():
    """Edge-on crossing: the slave bottom edge, pressed down, is restrained by the edge-edge
    penalty at δ ≈ P/εN. Newton convergence ⇒ the K = εN BᵀB tangent is right."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _build_bars()
    for t in (7, 8):
        ops.fix(t, 1, 1, 1)               # anchor the slave TOP edge
    for t in (5, 6):
        ops.fix(t, 1, 1, 0)               # bottom edge: x,y fixed; z free
    _add_springs()
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-edgeedge",
                "-edgeBand", 0.15, "-outward", 0.0, 0.0, 1.0)
    _load_bottom()

    assert _solve_static() == 0, "edge-edge static Newton did not converge (K_c?)"
    # the bottom edge settles at the penalty penetration: gap-mode εN + 2 springs ⇒
    # δ = P/(εN + 2k) ≈ P/εN. Both bottom nodes equal (symmetry).
    delta_pred = PTOT / (EPS + 2.0 * KSPRING)
    z5 = ZC + ops.nodeDisp(5, 3)
    z6 = ZC + ops.nodeDisp(6, 3)
    assert abs(z5 - z6) < 1.0e-9, "bottom edge not symmetric"
    assert abs(z5 - (-delta_pred)) < 1.0e-7, f"penetration {z5} != {-delta_pred}"
    # the restraint is εN-dominated: ~10³× stiffer than the springs alone (P/(2k))
    assert abs(z5) < 1.0e-3, f"edge-edge did not restrain the bottom edge ({z5})"


def test_e2_nts_passes_through():
    """THE FALSIFIER. Identical geometry + springs, but an NTS contact (slave NODES 5,6 vs the
    master facet): the nodes are 1 unit off the master FACE in x, so NONE penetrates ⇒ the NTS
    net force is ZERO and the bottom edge sinks to the pure-spring depth δ = P/(2k) — ~10³×
    deeper than the edge-edge restraint. NTS passes through where edge-edge holds."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    _build_bars()
    for t in (7, 8):
        ops.fix(t, 1, 1, 1)
    for t in (5, 6):
        ops.fix(t, 1, 1, 0)
    _add_springs()
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(3, "-slave", 5, 6)         # NTS slave NODE-set (the bottom edge nodes)
    ops.contact(1, 1, 3, 1.0e6, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)  # NTS penalty
    _load_bottom()

    assert _solve_static() == 0, "NTS static Newton did not converge"
    z5 = ZC + ops.nodeDisp(5, 3)
    # contact force ZERO ⇒ pure-spring response: the spring deforms by P/(2k) from ZC, so the
    # bottom edge reaches z = ZC − P/(2k) (the springs took the WHOLE load — NTS transmitted none).
    delta_spring = PTOT / (2.0 * KSPRING)
    assert abs(z5 - (ZC - delta_spring)) < 1.0e-6, f"NTS not at the pure-spring depth ({z5})"
    # the headline: NTS sinks ~10³× deeper than the edge-edge restraint (≈P/εN) — it passed through
    delta_ee = PTOT / (EPS + 2.0 * KSPRING)
    assert abs(z5) > 100.0 * delta_ee, "NTS did NOT pass through (it restrained the edge?)"


def test_e2_oblique_no_double_count():
    """EE-1 REGRESSION (the 3-reviewer gate MAJOR). An OBLIQUE pair: both facets tilted 45° so their
    normals are ~60° apart (align_fn≈0.5 < cos50° ⇒ non-facing) AND the face-mortar clip is
    NON-degenerate (unlike the exactly-perpendicular headline case where it is genuinely zero). The
    crossing edges stay axis-aligned (slave bottom ∥x, master top ∥y ⇒ n=ẑ) so δ stays clean. Face-
    mortar must STAND DOWN where edge-edge owns the pair ⇒ a SINGLE primitive: the bottom edge settles
    at the edge-edge depth δ≈P/εN, NOT ~2× stiffer from a superposed mortar clip pressure."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    # master: top edge m1-m2 ∥ y at z=0, facet tilted 45° (offset +x,−z) ⇒ normal ~(−1,0,−1)
    for t, x, y, z in [(1, 0.0, -1.0, 0.0), (2, 0.0, 1.0, 0.0),
                       (3, 0.3, 1.0, -0.3), (4, 0.3, -1.0, -0.3)]:
        ops.node(t, x, y, z)
        ops.fix(t, 1, 1, 1)
    # slave: bottom edge s1-s2 ∥ x at z=ZC, facet tilted 45° (offset +y,+z) ⇒ normal ~(0,−1,1)
    for t, x, y, z in [(5, -1.0, 0.0, ZC), (6, 1.0, 0.0, ZC),
                       (7, 1.0, 0.3, ZC + 0.3), (8, -1.0, 0.3, ZC + 0.3)]:
        ops.node(t, x, y, z)
    for t in (7, 8):
        ops.fix(t, 1, 1, 1)
    for t in (5, 6):
        ops.fix(t, 1, 1, 0)
    _add_springs()
    ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
    ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
    ops.contact(1, 1, 2, "-mortar", "-epsN", EPS, "-edgeedge",
                "-edgeBand", 0.15, "-outward", 0.0, 0.0, 1.0)
    _load_bottom()

    assert _solve_static() == 0, "oblique edge-edge static Newton did not converge"
    delta_pred = PTOT / (EPS + 2.0 * KSPRING)        # the SINGLE-primitive edge-edge depth
    z5 = ZC + ops.nodeDisp(5, 3)
    assert abs(z5 - (-delta_pred)) < 1.0e-6, \
        f"oblique pair not single-primitive ({z5}); mortar did NOT stand down (double-count)?"


def test_e2_edgeedge_absent_byte_identical():
    """Byte-identity: bars far apart (no edge pair routes) ⇒ a -mortar -edgeedge model evolves
    BITWISE identically to the same model with -edgeedge ABSENT (zero edge adapters either way;
    the mortar pair is inert too). A transient on the free slave bar with mass + velocity."""
    def _run(edgeedge):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        _build_bars(xoff=10.0)            # slave bar 10 units away in +x ⇒ no proximity
        for t in (5, 6, 7, 8):
            ops.fix(t, 1, 1, 0)           # x,y fixed; z free with mass + a z-velocity
            ops.mass(t, 1.0, 1.0, 1.0)
            ops.setNodeVel(t, 3, 0.2, "-commit")
        ops.contactSurface(1, "-master", 4, 1, 2, 3, 4)
        ops.contactSurface(2, "-slave-segments", 4, 5, 6, 7, 8)
        args = [1, 1, 2, "-mortar", "-epsN", EPS, "-outward", 0.0, 0.0, 1.0]
        if edgeedge:
            args += ["-edgeedge", "-edgeBand", 0.15]
        ops.contact(*args)
        ops.constraints("LadrunoContact")
        ops.numberer("Plain")
        ops.system("FullGeneral")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.test("NormDispIncr", 1.0e-12, 25, 0)
        ops.algorithm("Newton")
        ops.analysis("Transient")
        out = []
        for _ in range(8):
            assert ops.analyze(1, 1.0e-3) == 0
            out.append(tuple(ops.nodeDisp(t, 3) for t in (5, 6, 7, 8)))
        return out
    assert _run(True) == _run(False)
