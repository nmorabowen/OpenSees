"""ADR-63 #4a P2.2 + P2.3 — friction under the smoothed normal, and composition with ADR-60
finite-sliding re-emit over a CURVED master (the ADR-60 "exposed combo" closed).

P2.2 (friction). `segmentActive` already threads the smoothed normal `n` into BOTH the gap operator
`B` and the friction slip `gTvec` (`tangentPart(drel, n, gTvec)`), so friction composes with
`-smoothNormal` for free — no new machinery. On a FLAT master n_smooth == n_facet, so a frictional
slide with `-smoothNormal` is identical to the faceted one (D5 consistency, now for friction).

P2.3 (compose with -reemit on a curved master). A frictional block dragged ACROSS a convex curved
master needs `-reemit` to re-pair onto the next facet as it slides off its reference-config
candidate segments (the ADR-60 pass-through bug). With `-smoothNormal` the outward sense is a GLOBAL
surface datum captured once (never a per-slave sign that flips at a ridge crossing — ADR-60 R3), so
`-reemit + -mu + -smoothNormal` sustains contact across the crest. Without `-reemit` the block slides
off its frozen candidates and the press drives it through.

CAVEATS (documented, not gated here): (a) the AUTO global sign is ill-conditioned for a slave cloud
that grazes the master edge-on (seed ~⟂ the field — the F2/F3/F5 warning) → pass an explicit
`-outward`, as these tests do; (b) at a SHARP crest the P2.1 gap-aware ownership guard still keeps the
small-gap non-owner ⇒ mild near-apex over-stiffness (not a pass-through, not destabilizing on
realistic arcs) — full single-owner selection is deferred to ADR-57 #4b; (c) traction-continuity
(the C0-normal chatter fix) is a quality property, asserted only qualitatively via sustained contact.
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e5
H, L = 0.4, 2.0                 # convex arc z = H(1-(x/L)^2), crest z=H at x=0
XS = [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]


def _zarc(x):
    return H * (1.0 - (x / L) ** 2) if abs(x) <= L else -9.0


def _build_arc():
    """A faceted convex arc (8 quads) — a curved master with real inter-facet normal jumps."""
    for i, x in enumerate(XS):
        z = H * (1.0 - (x / L) ** 2)
        ops.node(100 + i, x, -0.6, z); ops.fix(100 + i, 1, 1, 1)
        ops.node(200 + i, x,  0.6, z); ops.fix(200 + i, 1, 1, 1)
    m = []
    for s in range(len(XS) - 1):
        m += [100 + s, 100 + s + 1, 200 + s + 1, 200 + s]
    return m


def _drag_over_arc(smooth, reemit, mu, vx=5.0):
    """Slave starts on the left of the arc, coasts across the crest (initial +x velocity) under a
    downward press with friction. Returns (reached_x, maxpen_on_arc, fell) — maxpen is the deepest
    penetration BELOW the arc surface while over it; fell = it dropped well below the master."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    m = _build_arc()
    x0 = -1.5
    z0 = _zarc(x0) + 1.0e-3
    ops.node(1, x0, 0.0, z0); ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *m)
    ops.contactSurface(20, "-slave", 1)
    kt = KN if mu > 0.0 else 0.0
    args = [1, 10, 20, KN, kt, mu, "-outward", 0.0, 0.0, 1.0]   # explicit global sign (caveat a)
    if reemit:
        args.append("-reemit")
    if smooth:
        args.append("-smoothNormal")
    ops.contact(*args)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(1, 0.0, 0.0, -15.0)
    ops.setNodeVel(1, 1, vx, "-commit")
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno"); ops.algorithm("Linear"); ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)
    maxpen = 0.0; reached = -9.0; z = z0
    for _ in range(40000):
        if ops.analyze(1, dt) != 0:
            break
        x = x0 + ops.nodeDisp(1, 1); z = z0 + ops.nodeDisp(1, 3)
        if not (math.isfinite(x) and math.isfinite(z)):
            break
        reached = max(reached, x)
        if -2.0 <= x <= 2.0:
            maxpen = max(maxpen, _zarc(x) - z)
        if reached > 1.7:
            break
    return reached, maxpen, (z < -1.0)


def test_p23_smooth_reemit_friction_sustains_curved_crossing():
    """-smoothNormal + -reemit + -mu over a convex curved master: the frictional block crosses the
    crest and stays in contact the whole way (bounded penetration, never driven through). This is
    the ADR-60 'exposed combo' (`-reemit + -mu` on a curved master) closed — the global-sign field
    holds the outward normal across the ridge, re-emit re-pairs across facets, friction acts."""
    reached, maxpen, fell = _drag_over_arc(smooth=True, reemit=True, mu=0.3)
    assert not fell, "frictional block fell through the curved master"
    assert reached > 1.0, f"block did not cross the crest (reached={reached})"
    assert maxpen < 0.05, f"contact not sustained across the crossing (maxpen={maxpen})"


def test_p23_curved_crossing_needs_reemit():
    """Composition is load-bearing: WITHOUT -reemit the smoothed frictional block slides off its
    reference-config candidate segments past the crest and the press drives it through (the ADR-60
    pass-through bug is unaffected by smoothing — smoothing fixes the SIGN, re-emit fixes the SEARCH)."""
    reached, maxpen, fell = _drag_over_arc(smooth=True, reemit=False, mu=0.3)
    assert fell or maxpen > 0.5, f"expected pass-through without -reemit (maxpen={maxpen}, fell={fell})"


# ------------------------------------------------------------------ P2.2 flat-friction consistency
NSEG = 14


def _build_flat_strip():
    for i in range(NSEG + 1):
        ops.node(1000 + i, float(i), -0.5, 0.0); ops.fix(1000 + i, 1, 1, 1)
        ops.node(2000 + i, float(i),  0.5, 0.0); ops.fix(2000 + i, 1, 1, 1)
    m = []
    for s in range(NSEG):
        m += [1000 + s, 1000 + s + 1, 2000 + s + 1, 2000 + s]
    return m


def _slide_flat(smooth, mu=0.1, vx=3.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    m = _build_flat_strip()
    ops.node(1, 0.5, 0.0, 1.0e-3); ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *m); ops.contactSurface(20, "-slave", 1)
    kt = KN if mu > 0.0 else 0.0
    args = [1, 10, 20, KN, kt, mu, "-outward", 0.0, 0.0, 1.0, "-reemit"]
    if smooth:
        args.append("-smoothNormal")
    ops.contact(*args)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(1, 0.0, 0.0, -10.0)
    ops.setNodeVel(1, 1, vx, "-commit")
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno"); ops.algorithm("Linear"); ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)
    for _ in range(4000):
        if ops.analyze(1, dt) != 0:
            break
    return 0.5 + ops.nodeDisp(1, 1)


def test_p22_friction_flat_smooth_equals_faceted():
    """P2.2 D5 consistency for FRICTION: on a flat master n_smooth == n_facet, so a frictional slide
    with -smoothNormal gives the SAME travel as the faceted one — friction composes with smoothing
    without perturbing the flat result (the slip uses the smoothed n, which equals the facet n here)."""
    x_faceted = _slide_flat(smooth=False)
    x_smooth = _slide_flat(smooth=True)
    assert abs(x_smooth - x_faceted) < 1.0e-9, (
        f"-smoothNormal perturbed the flat frictional slide (faceted x={x_faceted}, smooth x={x_smooth})")
