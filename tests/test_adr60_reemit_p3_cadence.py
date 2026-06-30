"""ADR-60 R7 — `-resortEvery` forced cadence (LS-DYNA BSORT) runtime smoke.

R7 split `-resortEvery N` from the migration floor: it is now a FORCED re-sort every N
commits (the Trigger's forceEvery / BSORT cycle override), independent of drift. The
Trigger logic itself is proven build-free in proto_reemit_selfcheck.cpp; this test only
confirms the runtime path — parser -> engine persistent Trigger -> forced re-handle —
runs end to end and keeps the finite-sliding fix intact (the forced cadence must not
break the no-pass-through guarantee).

Same sliding-block rig as test_adr60_reemit_p1 but with `-reemit -resortEvery 5`.
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6
NSEG = 6


def _build_strip():
    for i in range(NSEG + 1):
        ops.node(1000 + i, float(i), -0.5, 0.0); ops.fix(1000 + i, 1, 1, 1)
        ops.node(2000 + i, float(i),  0.5, 0.0); ops.fix(2000 + i, 1, 1, 1)
    mtags = []
    for s in range(NSEG):
        mtags += [1000 + s, 1000 + s + 1, 2000 + s + 1, 2000 + s]
    return mtags


def test_r7_resortEvery_forced_cadence_sustains_contact():
    """`-reemit -resortEvery 5`: a forced re-sort every 5 commits keeps the slave in
    contact across the whole strip (the BSORT cadence drives re-emit even where the
    migration metric alone might wait out its floor). No pass-through, slides across."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mtags = _build_strip()
    gap0 = 1.0e-3
    ops.node(1, 0.5, 0.0, gap0)
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *mtags)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0, "-reemit", "-resortEvery", 5)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -100.0)
    ops.setNodeVel(1, 1, 1.0, "-commit")
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)
    min_z = gap0
    for _ in range(4000):
        if ops.analyze(1, dt) != 0:
            break
        z = gap0 + ops.nodeDisp(1, 3)
        if z < min_z:
            min_z = z
    x_final = 0.5 + ops.nodeDisp(1, 1)
    assert x_final > 3.5, f"slave did not slide across with -resortEvery (x_final={x_final})"
    assert min_z > -0.05, f"slave fell through with -resortEvery (min_z={min_z})"
