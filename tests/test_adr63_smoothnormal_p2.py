"""ADR-63 #4a P2.1 — facet OWNERSHIP at a sharp convex ridge (the one real P1 limitation).

P1 shipped the smooth nodal-normal field (the R3 sign fix + C0 continuity). Its documented
bring-up limitation: a slave clearly on ONE facet near a sharp convex ridge has its closest point
on the ADJACENT (non-owning) facet land ON their SHARED ridge edge, where it reads a penetration
∝ the slave's LATERAL distance from the ridge — a large spurious EJECTING force → instability. The
faceted path only *incidentally* prunes it (at a 90° ridge the neighbour's per-pair orientDir makes
its faceted gap read non-penetrating); ADR-63 D2 replaced that per-pair sign with a global one, so
the smoothed path has no such prune.

Test rig: the slave's x,y DOFs are FIXED (only z free) so there is NO frictionless sliding to
confound the signal — a spurious adjacent-facet force shows purely as an upward z EJECTION. The
slave sits on the RIGHT facet near the ridge and is pushed straight DOWN. With the P2.1 ownership
guard the owning facet holds it near the surface (max upward excursion ≈ 0); WITHOUT the guard the
non-owning facet's large shared-edge gap throws it UP toward the ridge height (max_up ~ the lateral
offset). The guard keeps the OWNER and the at-apex contact (small gap) and drops only the spurious
LARGE-gap non-owner, so it does NOT reintroduce the pass-through that the blunt near-edge reject did.

Build-free oracle: contact_prototypes/proto_nodal_normal_selfcheck.cpp (shared-edge mask + the
gap-aware reject vs keep). Here it runs in compiled OpenSees.
"""
import math
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

KN = 1.0e6
DELTA = 1.0e-3
INV_SQRT2 = 1.0 / math.sqrt(2.0)


def _build_quad_tent():
    """Convex ridge along y at (x=0,z=1); eaves at x=±1,z=0. Two quad facets (planes x∓z=… ;
    the right facet is the plane x+z=1, outward (1,0,1)/√2). All FIXED. Returns (mtags, nps)."""
    ops.node(100, -1.0, -0.5, 0.0); ops.fix(100, 1, 1, 1)   # L0
    ops.node(101, -1.0,  0.5, 0.0); ops.fix(101, 1, 1, 1)   # L1
    ops.node(110,  0.0, -0.5, 1.0); ops.fix(110, 1, 1, 1)   # R0 (ridge)
    ops.node(111,  0.0,  0.5, 1.0); ops.fix(111, 1, 1, 1)   # R1 (ridge)
    ops.node(120,  1.0, -0.5, 0.0); ops.fix(120, 1, 1, 1)   # Rt0
    ops.node(121,  1.0,  0.5, 0.0); ops.fix(121, 1, 1, 1)   # Rt1
    left  = [100, 110, 111, 101]    # CCW: L0 → R0 → R1 → L1
    right = [110, 120, 121, 111]    # CCW: R0 → Rt0 → Rt1 → R1
    return left + right, 4


def _build_tri_tent():
    """The same convex ridge, two TRIANGLES sharing the ridge edge R0-R1 (nps=3). Right tri on the
    plane x+z=1 (outward (1,0,1)/√2)."""
    ops.node(110,  0.0, -0.5, 1.0); ops.fix(110, 1, 1, 1)   # R0 (ridge)
    ops.node(111,  0.0,  0.5, 1.0); ops.fix(111, 1, 1, 1)   # R1 (ridge)
    ops.node(120,  1.0,  0.0, 0.0); ops.fix(120, 1, 1, 1)   # Rt (right eave apex)
    ops.node(130, -1.0,  0.0, 0.0); ops.fix(130, 1, 1, 1)   # Lt (left eave apex)
    right = [110, 120, 111]         # right tri, shared edge R1-R0
    left  = [110, 111, 130]         # left tri,  shared edge R0-R1 (opposite ⇒ coherent)
    return right + left, 3


def _press_near_ridge_zonly(build, nps, xs_x):
    """Slave near the ridge over the RIGHT facet (plane x+z=1), x,y DOFs FIXED so it can only move in
    z. Start δ below the right facet, pushed straight DOWN. Returns (max_up, final_pen, diverged):
    max_up = the largest z RISE above the start (the ejection signature — the spurious non-owner
    force at the ridge pushes the slave UP toward the ridge height); final_pen = the settled
    penetration below the right facet surface (>0 ⇒ inside); diverged = the run went non-finite."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    mtags, _nps = build()
    assert _nps == nps
    surf_z = 1.0 - xs_x                       # right-facet surface height at x=xs_x (plane x+z=1)
    z0 = surf_z - DELTA                        # start δ below the surface (straight down)
    ops.node(1, xs_x, 0.0, z0)
    ops.fix(1, 1, 1, 0)                        # x,y FIXED, z FREE — no sliding
    ops.mass(1, 1.0, 1.0, 1.0)
    ops.contactSurface(10, "-master", nps, *mtags)
    ops.contactSurface(20, "-slave", 1)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0, "-smoothNormal")
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, 0.0, 0.0, -100.0)             # straight down (only z responds)
    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.integrator("CentralDifferenceLadruno")
    ops.algorithm("Linear")
    ops.analysis("Transient")
    dt = 0.5 * 2.0 * math.sqrt(1.0 / KN)
    max_up = 0.0
    z = z0
    diverged = False
    for _ in range(3000):
        if ops.analyze(1, dt) != 0:
            diverged = True
            break
        uz = ops.nodeDisp(1, 3)
        if not math.isfinite(uz):
            diverged = True
            break
        z = z0 + uz
        if uz > max_up:
            max_up = uz
    final_pen = surf_z - z                    # >0 ⇒ still below the right surface (in contact)
    return max_up, final_pen, diverged


def test_p2_quad_ridge_no_spurious_ejection():
    """Quad convex ridge, slave near the ridge (x,y fixed, z free): the ownership guard drops the
    spurious adjacent-facet force, so the slave stays near the surface (no upward ejection) and is
    held (not driven through). Without the guard the non-owner's shared-edge gap (≈0.15 here) throws
    the slave UP toward the ridge (max_up ~ 0.15)."""
    max_up, final_pen, diverged = _press_near_ridge_zonly(_build_quad_tent, 4, 0.15)
    assert not diverged, "run diverged — spurious adjacent-facet ejection at the ridge"
    assert max_up < 0.02, f"slave was EJECTED upward by the spurious adjacent facet (max_up={max_up})"
    assert final_pen > -0.05, f"slave left the surface / driven through (final_pen={final_pen})"


def test_p2_tri_ridge_no_spurious_ejection():
    """tri-3 convex ridge (nps=3) — the case the P1 tri coverage deliberately avoided (flat patch).
    Same ownership guard: no spurious upward ejection, slave held on the owning triangle."""
    max_up, final_pen, diverged = _press_near_ridge_zonly(_build_tri_tent, 3, 0.20)
    assert not diverged, "run diverged — spurious adjacent-triangle ejection at the ridge"
    assert max_up < 0.02, f"slave was EJECTED upward by the spurious adjacent triangle (max_up={max_up})"
    assert final_pen > -0.05, f"slave left the surface / driven through (final_pen={final_pen})"
