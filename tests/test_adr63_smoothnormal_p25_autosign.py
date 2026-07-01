"""ADR-63 #4a P2.5 — auto-sign robustness: `-smoothNormal` lifts the `-outward` requirement for
edge-grazing / side-approaching slaves on a curved master.

THE RECURRING CAVEAT (F2/F3/F5, and the P2.3 note). The frozen global outward sign was decided by a
single AGGREGATE vote  sign( Σ_a σ_a n_a · (slaveCentroid − masterCentroid) ). That dots an aggregate
normal (which nearly cancels on a curved/domed master) against a GLOBAL chord (which goes ~tangent to
the field when the slave cloud grazes the master edge-on) ⇒ vote·seed ≈ 0 is a coin-flip; a tiny
wrong-signed component flips the WHOLE field inward ⇒ silent pass-through even WITH `-smoothNormal`
(gap>0, "not penetrating"). So `-smoothNormal` was NOT a blanket `-outward` lift — the P2.3 gates all
pass an explicit `-outward` for exactly this reason.

THE FIX (this phase). On the AUTO path the sign now comes from a per-slave MAJORITY vote
(`LadrunoContactProjection::voteSignRobust`): each slave votes its LOCAL nearest-facet coherent normal
against its LOCAL separation vector. A closest-point separation is PARALLEL to the local normal by
construction, so each dot is well-conditioned (|cos|≈1) no matter where the cloud sits around the
master. The surface takes the weighted majority; a genuinely two-sided cloud still yields a low margin
⇒ the handler warns and recommends `-outward` (the ambiguity is DETECTED, not silently mis-signed).
`-outward` remains the explicit override + the refuse-fallback; `-smoothNormal` OFF stays byte-identical.

These gates run WITHOUT `-outward` — they PASS on the robust vote and would FAIL on the old aggregate
vote (the field would point inward for the edge-grazing slave ⇒ pass-through from step 1).
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


# --------------------------------------------------------------- (1) SIGN isolation: held, x,y fixed
def _press_edge_grazing(smooth, outward, x_slave=1.5):
    """A slave placed above the arc OFF TO THE SIDE at x=x_slave (edge-grazing — where the aggregate
    seed slave−masterCentroid is ~horizontal with a NEGATIVE z (slave below the master centroid height)
    ⇒ the old auto sign flips inward). x,y FIXED (only z free, the corrected-P1 discipline: a
    curved-master free body has no lateral equilibrium and confounds a 'held' assertion), pressed
    straight down. Returns maxpen = max(zarc(x) − z) over the run: a correct outward sign HOLDS it
    (small penetration ≈ press/kn); a flipped sign passes it through (maxpen → huge as it falls)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    m = _build_arc()
    z0 = _zarc(x_slave) + 1.0e-3
    ops.node(1, x_slave, 0.0, z0); ops.mass(1, 1.0, 1.0, 1.0)
    ops.fix(1, 1, 1, 0)                                    # x,y FIXED, z free — isolate the sign
    ops.contactSurface(10, "-master", 4, *m)
    ops.contactSurface(20, "-slave", 1)
    args = [1, 10, 20, KN, 0.0, 0.0]
    if outward:
        args += ["-outward", 0.0, 0.0, 1.0]
    if smooth:
        args.append("-smoothNormal")
    ops.contact(*args)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(1, 0.0, 0.0, -20.0)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno"); ops.algorithm("Linear"); ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)
    maxpen = 0.0
    for _ in range(3000):
        if ops.analyze(1, dt) != 0:
            break
        z = z0 + ops.nodeDisp(1, 3)
        if not math.isfinite(z):
            break
        maxpen = max(maxpen, _zarc(x_slave) - z)   # depth below the arc surface (small if held)
    return maxpen


def test_p25_autosign_holds_edge_grazing_slave():
    """`-smoothNormal` with the AUTO sign (NO -outward) holds an edge-grazing slave on the curved arc:
    the per-slave robust vote signs the field OUTWARD for it, so the downward press is resisted
    (bounded penetration), not passed through. The old aggregate vote would flip inward here."""
    maxpen = _press_edge_grazing(smooth=True, outward=False)
    assert maxpen < 0.05, f"auto-sign smoothed field did not hold the edge-grazing slave (maxpen={maxpen})"


def test_p25_autosign_matches_explicit_outward():
    """Sanity: for this over-the-master slave the AUTO robust sign gives the SAME held result as an
    explicit -outward — the robust vote recovers the direction the user would have had to supply."""
    p_auto = _press_edge_grazing(smooth=True, outward=False)
    p_out = _press_edge_grazing(smooth=True, outward=True)
    assert p_out < 0.05 and p_auto < 0.05, f"held? auto_pen={p_auto} outward_pen={p_out}"
    assert abs(p_auto - p_out) < 1.0e-3, f"auto sign diverged from -outward (auto={p_auto}, outward={p_out})"


# --------------------------------------------------- (2) curved crossing sustained WITHOUT -outward
def _drag_over_arc(reemit, mu, outward, vx=5.0):
    """Mirror of the P2.3 curved-arc drag but the sign is chosen by the AUTO robust vote (no -outward
    unless `outward`). Slave coasts across the crest under a downward press with friction. Returns
    (reached_x, maxpen_on_arc, fell)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    m = _build_arc()
    x0 = -1.5
    z0 = _zarc(x0) + 1.0e-3
    ops.node(1, x0, 0.0, z0); ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *m)
    ops.contactSurface(20, "-slave", 1)
    kt = KN if mu > 0.0 else 0.0
    args = [1, 10, 20, KN, kt, mu]
    if outward:
        args += ["-outward", 0.0, 0.0, 1.0]
    if reemit:
        args.append("-reemit")
    args.append("-smoothNormal")
    ops.contact(*args)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(1, 0.0, 0.0, -15.0)
    ops.setNodeVel(1, 1, vx, "-commit")
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno"); ops.algorithm("Linear"); ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)
    maxpen = 0.0; reached = -9.0
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


def test_p25_autosign_sustains_curved_crossing_no_outward():
    """The P2.3 gate re-run WITHOUT -outward: `-smoothNormal + -reemit + -mu` on a convex curved master
    now sustains the crossing on the AUTO robust sign alone (block crosses the crest, bounded
    penetration). This is the caveat lifted — smoothing fixes the SIGN robustly, re-emit fixes the
    SEARCH, and no explicit -outward is needed for an over-the-master cloud."""
    reached, maxpen, fell = _drag_over_arc(reemit=True, mu=0.3, outward=False)
    assert not fell, "frictional block fell through the curved master (auto sign flipped inward?)"
    assert reached > 1.0, f"block did not cross the crest on the auto sign (reached={reached})"
    assert maxpen < 0.05, f"contact not sustained across the crossing (maxpen={maxpen})"
