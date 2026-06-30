"""ADR-60 P2 — Coulomb friction under finite-sliding re-emission.

The 5-lens gate (decision D4) established that NTS friction needs NO transfer code at a
segment crossing: when a slave re-pairs onto a new segment, getOrCreateFrictionState
returns a FRESH slot (gT0 = current, gpT = 0 -> zero stick force at the crossing instant
-> already traction-continuous); a SLIPPING block is in kinetic friction (tT = mu*N)
regardless of the reset. This test validates that empirically — friction adapters are
correctly re-emitted across crossings and the friction state survives the re-handle.

A slave slides (+x velocity, frictional) across a long flat master strip under a gentle
downward press, with `-reemit`:
  * frictional vs frictionless travel: friction measurably shortens the slide AND the
    block stays in contact the whole way -> friction acted across the re-emitted crossings
    (not just the initial candidate set).
  * default (no -reemit) with friction still pass-through (the bug is friction-agnostic).
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6
NSEG = 14                # unit quads, x in [0, NSEG], y in [-0.5, 0.5], z = 0
PRESS = 10.0             # gentle downward load so friction lets the block slide far
VX0 = 3.0


def _build_strip():
    for i in range(NSEG + 1):
        ops.node(1000 + i, float(i), -0.5, 0.0); ops.fix(1000 + i, 1, 1, 1)
        ops.node(2000 + i, float(i),  0.5, 0.0); ops.fix(2000 + i, 1, 1, 1)
    mtags = []
    for s in range(NSEG):
        mtags += [1000 + s, 1000 + s + 1, 2000 + s + 1, 2000 + s]
    return mtags


def _slide(reemit, mu):
    """Frictional slide; return (x_final, min_z)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mtags = _build_strip()
    gap0 = 1.0e-3
    ops.node(1, 0.5, 0.0, gap0)
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *mtags)
    ops.contactSurface(20, "-slave", 1)
    kt = KN if mu > 0.0 else 0.0
    args = [1, 10, 20, KN, kt, mu, "-outward", 0.0, 0.0, 1.0]
    if reemit:
        args.append("-reemit")
    ops.contact(*args)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -PRESS)
    ops.setNodeVel(1, 1, VX0, "-commit")
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
    return 0.5 + ops.nodeDisp(1, 1), min_z


def test_p2_friction_decelerates_across_crossings():
    """-reemit: friction stays engaged across every segment crossing -> the frictional
    block stays in contact and slides measurably LESS far than the frictionless one (the
    energy it loses to friction is dissipated past the initial candidate set, so re-emit
    must have kept the friction adapters live across the crossings)."""
    x_free, zmin_free = _slide(reemit=True, mu=0.0)
    x_fric, zmin_fric = _slide(reemit=True, mu=0.1)
    assert zmin_free > -0.05, f"frictionless fell through despite -reemit (min_z={zmin_free})"
    assert zmin_fric > -0.05, f"frictional fell through despite -reemit (min_z={zmin_fric})"
    assert x_fric > 3.5, f"frictional block did not cross past the initial band (x={x_fric})"
    assert x_fric < 0.85 * x_free, f"friction had no sustained effect: x_fric={x_fric} x_free={x_free}"


def test_p2_friction_no_reemit_passes_through():
    """The pass-through bug is friction-agnostic: without -reemit the frictional block
    still loses contact past the frozen candidate set and the press drives it through."""
    x_fric, zmin_fric = _slide(reemit=False, mu=0.1)
    assert zmin_fric < -1.0, f"expected pass-through without -reemit (min_z={zmin_fric})"
