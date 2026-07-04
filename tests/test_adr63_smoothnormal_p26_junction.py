"""ADR-63 #4a P2.6 — QUANTITATIVE junction-traction-continuity gate (the C0-normal chatter fix).

The SECOND headline benefit of `-smoothNormal` (besides the R3 sign fix) is that the master normal is
C0-CONTINUOUS across facet junctions: adjacent facets share the same nodal normals at their shared edge,
so `n_smooth` has NO jump at a junction (the Abaqus TG §5.1.2 chatter driver — a jumping normal makes the
contact point oscillate under sliding). ADR-60/P1 asserted this only QUALITATIVELY (sustained contact).
This gate measures it directly.

THE OBSERVABLE — the contact normal's TILT as a function of position. A frictionless slave pressed onto
the master with a constant normal load P, held laterally at position x by a stiff spring, needs a lateral
hold force F_lat = P·(n_x/n_z) (the tangential component of the normal reaction). So F_lat/P = n_x/n_z is
the contact normal's tilt at x, and the stiff spring's force reads it out directly (a clean static solve —
no DisplacementControl, and reactions() does not expose contact forces under the LadrunoContact handler).
Sampling tilt(x) across a junction:
  * FACETED normal: piecewise-CONSTANT per facet ⇒ tilt(x) is a STAIRCASE that JUMPS at the junction.
  * SMOOTHED normal: `n_smooth` rotates CONTINUOUSLY ⇒ tilt(x) has no jump at the junction (C0).

This is a direct snapshot of the normal FIELD's continuity (each x is an independent static probe, no path
dependence), which is exactly the C0 property. Co-activation near the junction is transparent: two
co-active pairs split the normal force but share the same continuous `n_smooth`, so the tilt (a force
ratio) is unchanged.

RIG — a FLAT middle facet (z=HZ) flanked by two STEEP outer facets, so the junction at x=HALF is a
flat→steep transition where the FACETED normal jumps from vertical to a large tilt (the sharpest faceted
jump). All master nodes FIXED (constant nodal field ⇒ no -reemit). At each probe x: slave free in x,z (y
fixed) at the surface with a small penetration, a stiff lateral spring to a grounded twin holds x, a
constant press seats z; the spring force is the lateral hold force. Frictionless, `-outward` explicit
(this gate is about the normal DIRECTION continuity, not the sign — keep the sign trivially correct).

Observed (2026-07-02): faceted tilt 0 → −0.60 in ONE step at x=0.5; smoothed ramps 0 → −0.30 → −0.49
continuously (junction step −0.008). The gate locks that contrast.
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e4
KLAT = 1.0e5                    # stiff lateral hold spring (>> contact lateral force ⇒ probe x accurate)
HZ = 0.6                        # flat-middle height; outer facets rise HZ over 1.0 ⇒ steep-facet tilt 0.6
XCOL = [-1.5, -0.5, 0.5, 1.5]  # middle facet [-0.5,0.5] flat; junctions at x=±0.5
ZCOL = [0.0, HZ, HZ, 0.0]
HALF = 0.5                     # the +x junction (flat→steep) is at x=HALF
PRESS = 30.0
XEND = 1.15                    # probe span: flat facet, across the x=0.5 junction, onto the steep facet
NPROBE = 30


def _zsurf(x):
    for i in range(len(XCOL) - 1):
        if XCOL[i] - 1e-9 <= x <= XCOL[i + 1] + 1e-9:
            t = (x - XCOL[i]) / (XCOL[i + 1] - XCOL[i])
            return ZCOL[i] + t * (ZCOL[i + 1] - ZCOL[i])
    return -9.0


def _build_bump():
    for c, (x, z) in enumerate(zip(XCOL, ZCOL)):
        ops.node(100 + c, x, -0.6, z); ops.fix(100 + c, 1, 1, 1)
        ops.node(200 + c, x,  0.6, z); ops.fix(200 + c, 1, 1, 1)
    m = []
    for s in range(len(XCOL) - 1):
        m += [100 + s, 100 + s + 1, 200 + s + 1, 200 + s]
    return m


def _tilt_at(x, smooth):
    """One static probe: slave at (x, 0, surface−δ), free x,z (y fixed), stiff lateral spring to a
    grounded twin, constant press. Returns (rc, tilt) with tilt = F_spring/P = n_x/n_z (the normal tilt)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    m = _build_bump()
    z0 = _zsurf(x) - 1.0e-3
    ops.node(1, x, 0.0, z0); ops.fix(1, 0, 1, 0)          # y fixed; x,z free
    ops.node(2, x, 0.0, z0); ops.fix(2, 1, 1, 1)          # grounded twin (spring anchor)
    ops.uniaxialMaterial("Elastic", 1, KLAT)
    ops.element("zeroLength", 999, 2, 1, "-mat", 1, "-dir", 1)   # stiff lateral hold spring
    ops.contactSurface(10, "-master", 4, *m)
    ops.contactSurface(20, "-slave", 1)
    args = [1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0]
    if smooth:
        args.append("-smoothNormal")
    ops.contact(*args)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(1, 0.0, 0.0, -PRESS)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-9, 60); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0); ops.analysis("Static")
    rc = ops.analyze(1)
    f = ops.eleResponse(999, "force")
    fx = f[0] if isinstance(f, (list, tuple)) and f else float("nan")
    return rc, fx / PRESS


def _trace(smooth):
    xs = [XEND * i / (NPROBE - 1) for i in range(NPROBE)]
    tilts = []
    for x in xs:
        rc, t = _tilt_at(x, smooth)
        assert rc == 0, f"probe did not converge at x={x} (smooth={smooth})"
        tilts.append(t)
    return xs, tilts


def _jumps(xs, tilts):
    max_jump = 0.0
    junction_jump = 0.0
    for i in range(1, len(tilts)):
        d = abs(tilts[i] - tilts[i - 1])
        max_jump = max(max_jump, d)
        if xs[i - 1] < HALF <= xs[i]:
            junction_jump = max(junction_jump, d)
    return max_jump, junction_jump


def test_p26_smoothed_normal_tilt_continuous_across_junction():
    """Sample the contact normal's tilt across a flat→steep junction. FACETED tilt is a staircase that
    jumps ~0.6 at the junction (the chatter driver); SMOOTHED tilt is C0 — it varies continuously with a
    tiny junction step. Assert the smoothed junction jump and max per-step jump are far below the faceted
    staircase, quantifying the C0-normal chatter fix that P1 could only assert qualitatively."""
    xf, tf = _trace(smooth=False)
    xs, ts = _trace(smooth=True)

    fac_max, fac_jx = _jumps(xf, tf)
    smo_max, smo_jx = _jumps(xs, ts)
    print(f"\n[P2.6] faceted: max_jump={fac_max:.4f} junction_jump={fac_jx:.4f} "
          f"tilt[0]={tf[0]:.3f} tilt[-1]={tf[-1]:.3f}")
    print(f"[P2.6] smooth : max_jump={smo_max:.4f} junction_jump={smo_jx:.4f} "
          f"tilt[0]={ts[0]:.3f} tilt[-1]={ts[-1]:.3f}")

    # the faceted normal genuinely JUMPS at the junction (the discontinuity we remove)
    assert fac_jx > 0.3, (
        f"faceted junction jump unexpectedly small ({fac_jx:.4f}) — rig not exercising a real normal jump")
    # the smoothed normal is C0: the junction-crossing step is tiny in absolute terms ...
    assert smo_jx < 0.05, (
        f"smoothed normal jumped {smo_jx:.4f} at the junction — not C0-continuous")
    # ... and its whole trace has no step approaching the faceted staircase
    assert smo_max < 0.2 * fac_max, (
        f"smoothed max step-jump ({smo_max:.4f}) not clearly below the faceted staircase ({fac_max:.4f})")
    # sanity: both reach a comparable steep-facet tilt (continuity, not a different answer)
    assert abs(ts[-1] - tf[-1]) < 0.2, (
        f"smoothed vs faceted steep-facet tilt diverge (smooth {ts[-1]:.3f}, faceted {tf[-1]:.3f})")
