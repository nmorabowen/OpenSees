"""ADR-60 P1 — finite-sliding NTS re-emission: the no-pass-through gate.

A slave node slides (initial +x velocity, frictionless) across a multi-segment flat
master strip while a constant downward load presses it into contact. The NTS broad
phase builds its candidate set ONCE per handle() from the slave's position; as the
slave slides off that set's segments it must re-pair onto the segment it is now over.

  * `-reemit` ON  (ADR-60): the migration trigger forces a re-handle from the
    DEFORMED config, so the active segment is always a candidate -> the slave stays
    in contact across every segment crossing (bounded |z|).
  * default (no `-reemit`): the candidate set is frozen at the start cell, so once the
    slave slides past it there is no adapter -> SILENT pass-through -> the downward
    load drives the slave through the master (z -> large negative). This companion
    test reproduces the bug and proves the gate above is meaningful.

Build-free oracles: contact_prototypes/proto_reemit.py (no-pass-through) +
proto_reemit_selfcheck.cpp (the metric/trigger). Here it runs in compiled OpenSees.
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6
NSEG = 6                 # unit quads tiling x in [0, NSEG], y in [-0.5, 0.5], z = 0


def _build_strip():
    """A flat master strip of NSEG unit quads sharing edge nodes; all FIXED. Returns
    the flat master node-tag list (4 per segment, CCW: bot-i, bot-i+1, top-i+1, top-i)."""
    # node tags: bottom row 1000+i (y=-0.5), top row 2000+i (y=+0.5), i = 0..NSEG
    for i in range(NSEG + 1):
        ops.node(1000 + i, float(i), -0.5, 0.0); ops.fix(1000 + i, 1, 1, 1)
        ops.node(2000 + i, float(i),  0.5, 0.0); ops.fix(2000 + i, 1, 1, 1)
    mtags = []
    for s in range(NSEG):
        mtags += [1000 + s, 1000 + s + 1, 2000 + s + 1, 2000 + s]
    return mtags


def _slide(reemit):
    """Slave at x=0.5 just above the strip, mass 1, flung in +x, pressed down by a
    constant load. Run explicit; return (x_final, min_z) over the slide."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mtags = _build_strip()
    gap0 = 1.0e-3
    ops.node(1, 0.5, 0.0, gap0)
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *mtags)
    ops.contactSurface(20, "-slave", 1)
    args = [1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0]
    if reemit:
        args.append("-reemit")
    ops.contact(*args)
    # constant downward press + initial slide velocity
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -100.0)
    ops.setNodeVel(1, 1, 1.0, "-commit")          # +x slide at 1.0
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)
    min_z = gap0
    for _ in range(4000):                          # ~4 s -> x: 0.5 -> ~4.5
        if ops.analyze(1, dt) != 0:
            break
        z = gap0 + ops.nodeDisp(1, 3)
        if z < min_z:
            min_z = z
    x_final = 0.5 + ops.nodeDisp(1, 1)
    return x_final, min_z


def test_p1_reemit_sustains_contact_across_segments():
    """-reemit: the slave slides across several segments and STAYS in contact
    (bounded penetration); no silent pass-through."""
    x_final, min_z = _slide(reemit=True)
    assert x_final > 3.5, f"slave did not slide across (x_final={x_final})"
    assert min_z > -0.05, f"slave fell through despite -reemit (min_z={min_z})"


def test_p1_default_no_reemit_passes_through():
    """Default (no -reemit) reproduces the bug: once the slave slides off the frozen
    candidate set it loses contact and the load drives it through the master."""
    x_final, min_z = _slide(reemit=False)
    assert x_final > 3.5, f"slave did not slide (x_final={x_final})"
    assert min_z < -1.0, f"expected pass-through without -reemit, but min_z={min_z}"
