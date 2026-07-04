"""Contact review fix (2026-07) — offset-coordinate projection gate (the BLOCKER).

`LadrunoContactProjection::project()` converged on an ABSOLUTE residual tolerance
(1e-12) whose operand R = d.g has units length^2: its floating-point noise floor is
~eps*|X|*|g|, so for a model whose coordinates sit far from the origin — e.g. any
mm-unit building (h~500 mm facets at x~5e4 mm) — the tolerance was UNREACHABLE,
every projection returned "no valid projection", every pair evaluated inactive, and
contact SILENTLY vanished (no warning; the slave free-falls through the master).
Every shipped contact gate ran at unit scale near the origin, where the residual
products are exact in binary — which is exactly why this was never caught.

The fix adds a SCALE-FREE parametric-step escape (|dxi|+|deta| < 1e-8, parent
coordinates are O(1)) checked after the Newton update. No previously-converging
input is lost; on warped facets the escape's footpoint sits within ~tolP parametric
(~1e-8*h in position) of the residual-converged answer — physically nil (the
adversarial gate measured drift <= ~1e-9 over 5M trials) — and contact is alive at
any offset.

  * mm-building rig (rotated strip at (5e4,5e4,5e4), h=500): pre-fix the slave is
    driven straight through (min_d ~ -200); post-fix it is held at ~F/kn.
  * origin-vs-offset consistency: the settled penetration must be offset-invariant.

Build-free oracle: contact_prototypes/proto_projection_offset.cpp (13/21 checks fail
on the pre-fix header). Here the same physics runs in compiled OpenSees.
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6
NSEG = 4
TH_Z = 0.37          # rotation about z: irrational node coords (kills exact-binary luck)
TH_X = 0.21          # small tilt about x: the facet normal leaves the z axis too


def _rot(x, y, z):
    """Rotate (x,y,z) by TH_Z about z then TH_X about x (the oracle's facet transform)."""
    x1 = x * math.cos(TH_Z) - y * math.sin(TH_Z)
    y1 = x * math.sin(TH_Z) + y * math.cos(TH_Z)
    y2 = y1 * math.cos(TH_X) - z * math.sin(TH_X)
    z2 = y1 * math.sin(TH_X) + z * math.cos(TH_X)
    return x1, y2, z2


def _normal():
    """Outward unit normal of the rotated strip plane (the rotated +z)."""
    return _rot(0.0, 0.0, 1.0)


def _press_strip(h, L):
    """A flat NSEG-facet strip of facet size h, rotated by (TH_Z, TH_X), translated to
    (L, L, L); one slave over the strip middle, pressed along -n by a constant load.
    Explicit CDL, frictionless. Returns the min signed distance from the strip plane
    over the run (post-fix: held at ~ -F/KN; pre-fix at offset: free-fall, ~ -200*h/500)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    nx, ny, nz = _normal()
    mtags = []
    for i in range(NSEG + 1):
        xa, ya, za = _rot(i * h, -0.5 * h, 0.0)
        xb, yb, zb = _rot(i * h,  0.5 * h, 0.0)
        ops.node(1000 + i, xa + L, ya + L, za + L); ops.fix(1000 + i, 1, 1, 1)
        ops.node(2000 + i, xb + L, yb + L, zb + L); ops.fix(2000 + i, 1, 1, 1)
    for s in range(NSEG):
        mtags += [1000 + s, 1000 + s + 1, 2000 + s + 1, 2000 + s]
    # slave over the strip middle, lifted gap0 along +n
    gap0 = 1.0e-3 * h
    cx, cy, cz = _rot(1.5 * h, 0.0, 0.0)
    sx, sy, sz = cx + L + gap0 * nx, cy + L + gap0 * ny, cz + L + gap0 * nz
    ops.node(1, sx, sy, sz)
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", 4, *mtags)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", nx, ny, nz)
    # constant load along -n (no tangential component on the flat frictionless strip)
    F = 100.0
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, -F * nx, -F * ny, -F * nz)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)
    min_d = gap0
    nsteps = 3000
    done = 0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        done += 1
        ux, uy, uz = ops.nodeDisp(1, 1), ops.nodeDisp(1, 2), ops.nodeDisp(1, 3)
        # signed distance from the strip plane = ((P - plane point).n); the slave started
        # at gap0 along n from an on-plane point, so d = gap0 + u.n
        d = gap0 + ux * nx + uy * ny + uz * nz
        if d < min_d:
            min_d = d
    # anti-vacuous guards (adversarial gate F4): a run that dies at step 0 must FAIL,
    # and contact must actually ENGAGE (the slave crosses the plane: settled ~ -F/KN)
    assert done == nsteps, f"analysis died at step {done}/{nsteps}"
    assert min_d < 0.0, f"contact never engaged (min_d={min_d} stayed above the plane)"
    return min_d


def test_offset_mm_building_contact_is_alive():
    """h=500 facets at (5e4,5e4,5e4) — the canonical mm-unit building. Pre-fix the
    absolute residual tolerance is unreachable, project() dies, and the slave is driven
    straight through the strip (min_d ~ -200); post-fix it is held at ~ -F/KN."""
    min_d = _press_strip(h=500.0, L=5.0e4)
    assert min_d > -1.0, (
        f"contact is DEAD at mm-building coordinates (min_d={min_d}): "
        "the projection tolerance regression is back")


def test_offset_penetration_matches_origin():
    """The settled penetration must be offset-invariant: the same rig at the origin and
    at (5e4,5e4,5e4) resolves the same contact physics (footpoint error <= 1e-8*h)."""
    d_origin = _press_strip(h=500.0, L=0.0)
    d_offset = _press_strip(h=500.0, L=5.0e4)
    assert abs(d_offset - d_origin) < 1.0e-6, (
        f"offset changed the contact answer (origin min_d={d_origin}, "
        f"offset min_d={d_offset})")


def test_far_offset_unit_facets_alive():
    """Unit facets a million units out (the harshest ratio the parametric escape must
    cover): contact still resolves."""
    min_d = _press_strip(h=1.0, L=1.0e6)
    assert min_d > -1.0e-2, (
        f"contact is DEAD at L=1e6 with unit facets (min_d={min_d})")
