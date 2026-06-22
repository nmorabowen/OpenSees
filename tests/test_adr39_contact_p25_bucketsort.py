"""ADR-39 P2.5 — bucket-sort broad phase: verify == brute force.

The handler now bins the master segments into a spatial grid (LS-DYNA §26.11) and,
per slave node, emits adapters only for the segments in its 27-neighbour cell — a
SUPERSET of every geometrically-near pair — instead of all nSeg (brute force). The
narrow phase (projection/gap/penalty) is unchanged, so the contact RESULT is
identical; only dead far-pair adapters are dropped.

THE GATE = `verify == brute force`: run the SAME multi-segment model twice — once on
the default bucket grid, once with `-cell 1e12` (one giant bucket ⇒ every segment is
a candidate = brute force) — and assert the slave displacements are bitwise identical.
If the broad phase ever dropped a truly-near pair, the default-cell result would
DIFFER from brute force and this test would fail. (Algorithm + 33× pruning validated
oracle-first in contact_prototypes/proto_bucket_sort.py, P0a 4/4.)
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6
P = 1.0e2


def _ntag(ix, iy, N):
    return 1 + ix * (N + 1) + iy


def _build(N, cell=None):
    """A rigid N×N master grid of unit quad segments at z=0 (all master nodes fixed),
    with one slave just-penetrated at the centre of each segment, free-z, loaded −P.
    Returns the list of slave tags. `cell` (if given) sets the `-cell` lmax frac."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    # master grid nodes (fixed)
    for ix in range(N + 1):
        for iy in range(N + 1):
            t = _ntag(ix, iy, N)
            ops.node(t, float(ix), float(iy), 0.0)
            ops.fix(t, 1, 1, 1)
    # flat master connectivity: N*N quads, CCW from +z
    master = []
    for cx in range(N):
        for cy in range(N):
            master += [_ntag(cx, cy, N), _ntag(cx + 1, cy, N),
                       _ntag(cx + 1, cy + 1, N), _ntag(cx, cy + 1, N)]
    # one slave at each segment centre, just penetrated (−1e-8 z), free z, loaded −P
    slaves = []
    stag = 10000
    for cx in range(N):
        for cy in range(N):
            stag += 1
            ops.node(stag, cx + 0.5, cy + 0.5, -1.0e-8)
            ops.fix(stag, 1, 1, 0)
            slaves.append(stag)
    ops.contactSurface(10, "-master", 4, *master)
    ops.contactSurface(20, "-slave", *slaves)
    if cell is None:
        ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)
    else:
        ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0, "-cell", cell)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for s in slaves:
        ops.load(s, 0.0, 0.0, -P)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-10, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "multi-segment contact did not converge"
    return slaves


def test_bucketsort_equals_bruteforce():
    """THE gate: default bucket grid vs `-cell 1e12` (= brute force) → identical."""
    N = 5  # 25 segments → the grid meaningfully prunes (each slave sees ≤9, not 25)
    slaves = _build(N, cell=None)
    disp_bucket = [ops.nodeDisp(s, 3) for s in slaves]
    slaves2 = _build(N, cell=1.0e12)   # one bucket → every segment a candidate
    disp_brute = [ops.nodeDisp(s, 3) for s in slaves2]
    assert slaves == slaves2
    for s, db, bf in zip(slaves, disp_bucket, disp_brute):
        assert db == bf, f"slave {s}: bucket dz={db} != brute dz={bf} (broad phase dropped a pair?)"
    # and every slave actually contacted (dz ≈ −P/kn, not free-fall)
    for db in disp_bucket:
        assert abs(db - (-P / KN)) / (P / KN) < 0.02, f"dz={db} != −P/kn={-P/KN}"


def test_bucketsort_interior_segment():
    """A slave over an INTERIOR segment of a 4×4 grid penetrates with the correct
    P/kn — the broad phase found the right (non-boundary) segment as a candidate."""
    N = 4
    slaves = _build(N, cell=None)
    # the centre-most segment's slave is an interior one; check them all = P/kn
    for s in slaves:
        dz = ops.nodeDisp(s, 3)
        assert abs(dz - (-P / KN)) / (P / KN) < 0.02, f"interior slave {s} dz={dz}"


def test_bucketsort_single_segment_regression():
    """A 1-segment master (the P2b base case) still works through the broad phase
    (grid trivially returns the one segment)."""
    slaves = _build(1, cell=None)
    assert len(slaves) == 1
    dz = ops.nodeDisp(slaves[0], 3)
    assert abs(dz - (-P / KN)) / (P / KN) < 0.02, f"single-segment dz={dz}"
