"""ADR-60 R5 — slide-off-the-surface is already handled cleanly under -reemit.

The R5 remediation item (BLOCKER-SLIDE-OFF) asked whether a slave that departs the
master ENTIRELY can hold spurious traction (a stale adapter / edge clamp). Verified
EMPIRICALLY (Ladruno_scripts/_probe_r5_slideoff.py): under -reemit the existing
machinery already covers it, so R5 needs no code change — this test is the regression
guard for that conclusion:

  * the narrow phase gates on penetrating-AND-in-bounds (LadrunoContactProjection::
    evalSegment), so an out-of-bounds (slid-off) slave gets ZERO force;
  * the migration trigger + friction-slot GC drop the stale adapter within a bounded
    window;
  * D4 fresh-slot re-engagement keeps any later re-pairing traction-continuous.

A frictional slave is flung off the +x END of a FINITE strip. Once past the end it must
be in free flight: contact force 0, falling under the press, and NOT decelerated in x by
any phantom friction. (The frozen path is dominated by the pre-existing pass-through bug
+ a stale force snapshot, so only the -reemit path is asserted here.)
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6
NSEG = 6                 # FINITE strip: x in [0, 6]
MU = 0.1
PRESS = 10.0
VX0 = 8.0


def test_r5_slideoff_end_no_spurious_traction():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i in range(NSEG + 1):
        ops.node(1000 + i, float(i), -0.5, 0.0); ops.fix(1000 + i, 1, 1, 1)
        ops.node(2000 + i, float(i),  0.5, 0.0); ops.fix(2000 + i, 1, 1, 1)
    mtags = []
    for s in range(NSEG):
        mtags += [1000 + s, 1000 + s + 1, 2000 + s + 1, 2000 + s]
    gap0 = 1.0e-3
    x0 = 0.5
    ops.node(1, x0, 0.0, gap0)
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *mtags)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, KN, MU, "-outward", 0.0, 0.0, 1.0, "-reemit")
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -PRESS)
    ops.setNodeVel(1, 1, VX0, "-commit")
    ops.constraints("LadrunoContact")
    ops.numberer("Plain"); ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno"); ops.algorithm("Linear")
    ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)

    vx_on_strip = None
    for _ in range(4000):
        if ops.analyze(1, dt) != 0:
            break
        x = x0 + ops.nodeDisp(1, 1)
        if x < float(NSEG) - 0.2:        # last reading while still solidly on the strip
            vx_on_strip = ops.nodeVel(1, 1)

    x_final = x0 + ops.nodeDisp(1, 1)
    z_final = gap0 + ops.nodeDisp(1, 3)
    vx_final = ops.nodeVel(1, 1)
    fn_final = ops.ladrunoContactForce(1)

    assert x_final > float(NSEG) + 0.5, f"slave did not slide off the end (x={x_final})"
    # departed the surface ⇒ no contact force is held (stale-snapshot is cleared by re-handle)
    assert fn_final == 0.0, f"spurious contact force after slide-off (Fn={fn_final})"
    # free fall under the press ⇒ no phantom NORMAL traction propping it up
    assert z_final < -5.0, f"slave not in free fall after departure (z={z_final}); normal force held it"
    # no phantom TANGENTIAL traction ⇒ x-velocity is retained off-surface (only on-strip friction acted)
    assert vx_on_strip is not None
    assert vx_final > 0.95 * vx_on_strip, (
        f"x-velocity decayed off-surface (vx_final={vx_final} vs on-strip {vx_on_strip}); phantom friction")
