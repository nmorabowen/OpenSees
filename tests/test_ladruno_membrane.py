"""WP-3 Zone-A verification battery for LadrunoMembrane (ELE_TAG_LadrunoMembrane = 33022).

Covers:
  1. Parse-time & geometry refusals (< 3 nodes, unknown flag/method, collinear/zero-area,
     non-planar warped polygon, invalid -orthoRatio).
  2. Two-way rectangular bay ('tributary') onto ElasticBeam3d and ForceBeamColumn3d:
     yield-line tributary areas, partial trapezoidal frameLoads, global reaction
     equilibrium (sum Fz = q*A, sum Mx, sum My), and ForceBeamColumn3d internal section
     bending moment.
  3. One-way strip projection ('oneway') at 0 deg, 90 deg, oblique 30 deg, and
     single-sided cantilever overhang.
  4. Poisson membrane analogy ('fea_poisson') isotropic (eta=1) and orthotropic
     (eta=1e-4 -> oneway limit).
  5. DKT plate bending ('fea_plate') on 4-edge simple support (rank >= 2, signed
     corner uplift preserved) and single-edge cantilever (rank < 2, clamped rotation).
  6. Complex geometries: non-convex L-shape, interior '-hole' cutout, interior
     T-junction beam, partial-span edge beam, and isolated '-supportNodes' with mixed
     3-DOF / 6-DOF nodes.
  7. Multi-pattern loading ('Dead' + loadConst + 'Live') on ElasticBeam3d and
     ForceBeamColumn3d verifying zero double-counting via m_activeBeamLoads.
  8. Lazy frame resolution: declaring LadrunoMembrane BEFORE supporting beam elements.
  9. Guyan-condensed in-plane CST diaphragm ('-stiff E nu t' and '-thick t -mat tag'),
     consistent vs '-lumped' mass, ground-motion rigid-body acceleration probe
     (a_rel = -a_g), and Rayleigh damping (betaK > 0, betaKc > 0).
  10. Database sendSelf/recvSelf round-trip and live-domain state restoration.
"""
import math

import pytest
from _testbed import ops
from _testbed.roundtrip import database_roundtrip

pytestmark = [pytest.mark.zone_a]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _build_rect_bay_3d(
    Lx: float = 6.0,
    Ly: float = 4.0,
    q: float = 10.0,
    method: str = "tributary",
    beam_type: str = "elastic",
    declare_membrane_first: bool = False,
    extra_membrane_args: list | None = None,
):
    """Build a 3D rectangular slab (0..Lx) x (0..Ly) at z=0 with 4 perimeter beams."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)

    coords = {
        1: (0.0, 0.0, 0.0),
        2: (Lx, 0.0, 0.0),
        3: (Lx, Ly, 0.0),
        4: (0.0, Ly, 0.0),
    }
    for tag, (x, y, z) in coords.items():
        ops.node(tag, x, y, z)
        ops.fix(tag, 1, 1, 1, 1, 1, 1)

    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    E, G, A, J, Iy, Iz = 30.0e6, 12.5e6, 0.18, 0.002, 0.003, 0.003

    if beam_type == "force":
        ops.section("Elastic", 1, E, A, Iz, Iy, G, J)
        ops.beamIntegration("Lobatto", 1, 1, 5)

    def _add_beams():
        edges = [(1, 1, 2), (2, 2, 3), (3, 3, 4), (4, 4, 1)]
        for etag, ni, nj in edges:
            if beam_type == "elastic":
                ops.element("elasticBeamColumn", etag, ni, nj, A, E, G, J, Iy, Iz, 1)
            else:
                ops.element("forceBeamColumn", etag, ni, nj, 1, 1)

    def _add_membrane():
        args = [
            "LadrunoMembrane", 100,
            "-nodes", 1, 2, 3, 4,
            "-frames", 1, 2, 3, 4,
            "-method", method,
            "-pressure", q,
        ]
        if extra_membrane_args:
            args.extend(extra_membrane_args)
        ops.element(*args)

    if declare_membrane_first:
        _add_membrane()
        _add_beams()
    else:
        _add_beams()
        _add_membrane()

    return coords


def _solve_static(load_factor: float = 1.0, pattern_tag: int = 1, ele_tag: int = 100):
    ops.timeSeries("Linear", pattern_tag)
    ops.pattern("Plain", pattern_tag, pattern_tag)
    ops.eleLoad("-ele", ele_tag, "-type", "-surfaceLoad", load_factor)
    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    rc = ops.analyze(1)
    assert rc == 0
    ops.reactions()


def _check_global_equilibrium_3d(coords: dict[int, tuple[float, float, float]], q_total: float, area: float, x_c: float, y_c: float):
    """Verify sum(R_z) == -q_total * area and overturning reaction moments about origin."""
    expected_fz = -q_total * area
    expected_mx = expected_fz * y_c
    expected_my = -expected_fz * x_c

    sum_fx = 0.0
    sum_fy = 0.0
    sum_fz = 0.0
    sum_mx = 0.0
    sum_my = 0.0
    for n, (x, y, _z) in coords.items():
        rx = ops.nodeReaction(n, 1)
        ry = ops.nodeReaction(n, 2)
        rz = ops.nodeReaction(n, 3)
        mx = ops.nodeReaction(n, 4)
        my = ops.nodeReaction(n, 5)
        sum_fx += rx
        sum_fy += ry
        sum_fz += rz
        sum_mx += mx + y * rz
        sum_my += my - x * rz

    assert math.isclose(sum_fx, 0.0, abs_tol=1e-8)
    assert math.isclose(sum_fy, 0.0, abs_tol=1e-8)
    assert math.isclose(sum_fz, expected_fz, rel_tol=1e-6, abs_tol=1e-8)
    assert math.isclose(sum_mx, expected_mx, rel_tol=1e-5, abs_tol=1e-7)
    assert math.isclose(sum_my, expected_my, rel_tol=1e-5, abs_tol=1e-7)


# ---------------------------------------------------------------------------
# 1. Parse-time & geometry refusals
# ---------------------------------------------------------------------------

_OPENSEES_ERR = (RuntimeError, ValueError, getattr(ops, "OpenSeesError", RuntimeError))


def test_refusal_fewer_than_three_nodes():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 4.0, 0.0, 0.0)
    with pytest.raises(_OPENSEES_ERR):
        ops.element("LadrunoMembrane", 1, "-nodes", 1, 2, "-supportNodes", 1, 2)


def test_refusal_unknown_method():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 4.0, 0.0, 0.0)
    ops.node(3, 0.0, 4.0, 0.0)
    with pytest.raises(_OPENSEES_ERR):
        ops.element(
            "LadrunoMembrane", 1,
            "-nodes", 1, 2, 3,
            "-supportNodes", 1, 2, 3,
            "-method", "bogus_engine",
        )


def test_refusal_unknown_flag():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 4.0, 0.0, 0.0)
    ops.node(3, 0.0, 4.0, 0.0)
    with pytest.raises(_OPENSEES_ERR):
        ops.element(
            "LadrunoMembrane", 1,
            "-nodes", 1, 2, 3,
            "-supportNodes", 1, 2, 3,
            "-notAValidFlag", 123,
        )


def test_refusal_invalid_ortho_ratio():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 4.0, 0.0, 0.0)
    ops.node(3, 0.0, 4.0, 0.0)
    with pytest.raises(_OPENSEES_ERR):
        ops.element(
            "LadrunoMembrane", 1,
            "-nodes", 1, 2, 3,
            "-supportNodes", 1, 2, 3,
            "-orthoRatio", 0.0,
        )


def test_refusal_collinear_zero_area_polygon():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 2.0, 0.0, 0.0)
    ops.node(3, 4.0, 0.0, 0.0)
    with pytest.raises(_OPENSEES_ERR):
        ops.element(
            "LadrunoMembrane", 1,
            "-nodes", 1, 2, 3,
            "-supportNodes", 1, 2, 3,
        )


def test_refusal_warped_nonplanar_polygon():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 4.0, 0.0, 0.0)
    ops.node(3, 4.0, 4.0, 0.25)  # 25 cm warp on 4m panel >> 1e-3 * L_char
    ops.node(4, 0.0, 4.0, 0.0)
    with pytest.raises(_OPENSEES_ERR):
        ops.element(
            "LadrunoMembrane", 1,
            "-nodes", 1, 2, 3, 4,
            "-supportNodes", 1, 2, 3, 4,
        )


# ---------------------------------------------------------------------------
# 2. Two-way rectangular bay ('tributary')
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("beam_type", ["elastic", "force"])
def test_rectangular_tributary_two_way_yield_line(beam_type: str):
    Lx, Ly, q = 6.0, 4.0, 12.0
    coords = _build_rect_bay_3d(Lx=Lx, Ly=Ly, q=q, method="tributary", beam_type=beam_type)
    _solve_static(1.0)

    area_resp = ops.eleResponse(100, "area")
    assert len(area_resp) == 1
    assert math.isclose(area_resp[0], Lx * Ly, rel_tol=1e-12)

    tot_resp = ops.eleResponse(100, "totalLoad")
    assert len(tot_resp) == 3
    assert math.isclose(tot_resp[2], q * Lx * Ly, rel_tol=1e-6)

    # Short beams (2 and 4, length Ly=4) get triangle area Ly^2 / 4 = 4.0 m^2 each
    # Long beams (1 and 3, length Lx=6) get trapezoid area (Lx*Ly - 2*4.0)/2 = 8.0 m^2 each
    trib = ops.eleResponse(100, "tributaryAreas")
    assert len(trib) == 4
    expected_short = (Ly ** 2) / 4.0
    expected_long = 0.5 * (Lx * Ly - 2.0 * expected_short)
    assert math.isclose(trib[0], expected_long, rel_tol=2e-3)
    assert math.isclose(trib[1], expected_short, rel_tol=2e-3)
    assert math.isclose(trib[2], expected_long, rel_tol=2e-3)
    assert math.isclose(trib[3], expected_short, rel_tol=2e-3)

    # frameLoads returns [aL, bL, wa, wb] per supporting frame (4 * 4 = 16 floats)
    floads = ops.eleResponse(100, "frameLoads")
    assert len(floads) == 16
    beam_lengths = [Lx, Ly, Lx, Ly]
    total_f_seg = 0.0
    for k in range(4):
        aL, bL, wa, wb = floads[4 * k : 4 * k + 4]
        assert 0.0 <= aL < bL <= 1.0 + 1e-12
        total_f_seg += 0.5 * (wa + wb) * (bL - aL) * beam_lengths[k]
    assert math.isclose(total_f_seg, q * Lx * Ly, rel_tol=1e-6)

    _check_global_equilibrium_3d(coords, q, Lx * Ly, 0.5 * Lx, 0.5 * Ly)


# ---------------------------------------------------------------------------
# 3. One-way strip projection ('oneway')
# ---------------------------------------------------------------------------

def test_oneway_strip_axis_aligned_and_oblique():
    # Axis-aligned spanning along Y (0, 1, 0): beams 1 and 3 carry 50% each, beams 2 and 4 carry 0%
    Lx, Ly, q = 6.0, 4.0, 8.0
    coords = _build_rect_bay_3d(
        Lx=Lx,
        Ly=Ly,
        q=q,
        method="oneway",
        extra_membrane_args=["-spanAxis", 0.0, 1.0, 0.0],
    )
    _solve_static(1.0)
    trib = ops.eleResponse(100, "tributaryAreas")
    assert math.isclose(trib[0], 0.5 * Lx * Ly, rel_tol=1e-4)
    assert math.isclose(trib[1], 0.0, abs_tol=1e-8)
    assert math.isclose(trib[2], 0.5 * Lx * Ly, rel_tol=1e-4)
    assert math.isclose(trib[3], 0.0, abs_tol=1e-8)
    _check_global_equilibrium_3d(coords, q, Lx * Ly, 0.5 * Lx, 0.5 * Ly)

    # Oblique 30-degree span direction: total load and overturning moments still conserved
    theta = math.radians(30.0)
    coords = _build_rect_bay_3d(
        Lx=Lx,
        Ly=Ly,
        q=q,
        method="oneway",
        extra_membrane_args=["-spanAxis", math.cos(theta), math.sin(theta), 0.0],
    )
    _solve_static(1.0)
    assert math.isclose(sum(ops.eleResponse(100, "tributaryAreas")), Lx * Ly, rel_tol=1e-5)
    _check_global_equilibrium_3d(coords, q, Lx * Ly, 0.5 * Lx, 0.5 * Ly)


# ---------------------------------------------------------------------------
# 4. Poisson membrane analogy ('fea_poisson') & orthotropic limit
# ---------------------------------------------------------------------------

def test_fea_poisson_isotropic_and_orthotropic():
    Lx, Ly, q = 6.0, 4.0, 10.0
    # Isotropic (eta = 1.0)
    coords = _build_rect_bay_3d(Lx=Lx, Ly=Ly, q=q, method="fea_poisson")
    _solve_static(1.0)
    trib_iso = ops.eleResponse(100, "tributaryAreas")
    assert math.isclose(sum(trib_iso), Lx * Ly, rel_tol=1e-6)
    # Symmetry: opposite beams match, long beams carry more than short beams
    assert math.isclose(trib_iso[0], trib_iso[2], rel_tol=1e-5)
    assert math.isclose(trib_iso[1], trib_iso[3], rel_tol=1e-5)
    assert trib_iso[0] > trib_iso[1]
    _check_global_equilibrium_3d(coords, q, Lx * Ly, 0.5 * Lx, 0.5 * Ly)

    # Orthotropic (eta = 1e-4 along Y): approaches one-way spanning to beams 1 & 3
    coords = _build_rect_bay_3d(
        Lx=Lx,
        Ly=Ly,
        q=q,
        method="fea_poisson",
        extra_membrane_args=["-spanAxis", 0.0, 1.0, 0.0, "-orthoRatio", 1e-4],
    )
    _solve_static(1.0)
    trib_ortho = ops.eleResponse(100, "tributaryAreas")
    assert math.isclose(sum(trib_ortho), Lx * Ly, rel_tol=1e-6)
    # > 96% of total load goes to the two beams perpendicular to the stiff axis (1 and 3)
    assert (trib_ortho[0] + trib_ortho[2]) / (Lx * Ly) > 0.96
    _check_global_equilibrium_3d(coords, q, Lx * Ly, 0.5 * Lx, 0.5 * Ly)


# ---------------------------------------------------------------------------
# 5. DKT plate bending ('fea_plate') — 2D support rank & 1D cantilever rank < 2
# ---------------------------------------------------------------------------

def test_fea_plate_two_way_and_single_girder_cantilever():
    Lx, Ly, q = 5.0, 5.0, 10.0
    coords = _build_rect_bay_3d(Lx=Lx, Ly=Ly, q=q, method="fea_plate")
    _solve_static(1.0)
    trib = ops.eleResponse(100, "tributaryAreas")
    assert math.isclose(sum(trib), Lx * Ly, rel_tol=1e-6)
    # Fourfold symmetry on square bay
    for a in trib:
        assert math.isclose(a, 0.25 * Lx * Ly, rel_tol=1e-3)
    _check_global_equilibrium_3d(coords, q, Lx * Ly, 0.5 * Lx, 0.5 * Ly)

    # Single-girder cantilever (support rank < 2): only beam 1 (along y=0) supports the slab
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 4.0, 0.0, 0.0)
    ops.node(3, 4.0, 2.0, 0.0)
    ops.node(4, 0.0, 2.0, 0.0)
    for n in (1, 2):
        ops.fix(n, 1, 1, 1, 1, 1, 1)
    for n in (3, 4):
        ops.fix(n, 1, 1, 1, 1, 1, 1)

    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.element("elasticBeamColumn", 1, 1, 2, 0.18, 30e6, 12.5e6, 0.002, 0.003, 0.003, 1)
    ops.element(
        "LadrunoMembrane", 100,
        "-nodes", 1, 2, 3, 4,
        "-frames", 1,
        "-method", "fea_plate",
        "-pressure", q,
    )
    _solve_static(1.0)
    trib_cant = ops.eleResponse(100, "tributaryAreas")
    assert len(trib_cant) == 1
    assert math.isclose(trib_cant[0], 4.0 * 2.0, rel_tol=1e-5)
    rz_beam = ops.nodeReaction(1, 3) + ops.nodeReaction(2, 3)
    assert math.isclose(rz_beam, -q * 8.0, rel_tol=1e-5)


# ---------------------------------------------------------------------------
# 6. Complex topologies: L-shape, -hole, T-junction, partial beam, -supportNodes
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("method", ["tributary", "fea_poisson"])
def test_l_shaped_nonconvex_polygon_and_hole(method: str):
    """Non-convex L-shaped slab (area = 36 - 9 = 27 m^2) with an interior 1x1 m hole (net = 26 m^2)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)

    outer = [
        (1, 0.0, 0.0),
        (2, 6.0, 0.0),
        (3, 6.0, 3.0),
        (4, 3.0, 3.0),
        (5, 3.0, 6.0),
        (6, 0.0, 6.0),
    ]
    hole = [
        (11, 1.0, 1.0),
        (12, 2.0, 1.0),
        (13, 2.0, 2.0),
        (14, 1.0, 2.0),
    ]
    coords = {}
    for n, x, y in outer + hole:
        ops.node(n, x, y, 0.0)
        ops.fix(n, 1, 1, 1, 1, 1, 1)
        coords[n] = (x, y, 0.0)

    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    E, G, A, J, Iy, Iz = 30e6, 12.5e6, 0.18, 0.002, 0.003, 0.003
    edges = [(1, 1, 2), (2, 2, 3), (3, 3, 4), (4, 4, 5), (5, 5, 6), (6, 6, 1)]
    for etag, ni, nj in edges:
        ops.element("elasticBeamColumn", etag, ni, nj, A, E, G, J, Iy, Iz, 1)

    q = 5.0
    ops.element(
        "LadrunoMembrane", 100,
        "-nodes", 1, 2, 3, 4, 5, 6,
        "-hole", 11, 12, 13, 14,
        "-frames", 1, 2, 3, 4, 5, 6,
        "-method", method,
        "-pressure", q,
    )
    _solve_static(1.0)

    net_area = 26.0
    assert math.isclose(ops.eleResponse(100, "area")[0], net_area, rel_tol=1e-10)
    assert math.isclose(sum(ops.eleResponse(100, "tributaryAreas")), net_area, rel_tol=1e-5)

    # Hole nodes have no beams attached and must receive zero vertical reaction
    for hn, _, _ in hole:
        assert math.isclose(ops.nodeReaction(hn, 3), 0.0, abs_tol=1e-12)

    # Total vertical reaction on outer perimeter equals -q * net_area
    rz_outer = sum(ops.nodeReaction(n, 3) for n, _, _ in outer)
    assert math.isclose(rz_outer, -q * net_area, rel_tol=1e-5)


def test_interior_t_junction_beam_and_partial_edge_beam():
    """6m x 4m slab split by an interior beam (5) from (3,0) to (3,4), plus partial edge."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)

    pts = {
        1: (0.0, 0.0, 0.0),
        2: (6.0, 0.0, 0.0),
        3: (6.0, 4.0, 0.0),
        4: (0.0, 4.0, 0.0),
        5: (3.0, 0.0, 0.0),
        6: (3.0, 4.0, 0.0),
    }
    for n, (x, y, z) in pts.items():
        ops.node(n, x, y, z)
        ops.fix(n, 1, 1, 1, 1, 1, 1)

    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    E, G, A, J, Iy, Iz = 30e6, 12.5e6, 0.18, 0.002, 0.003, 0.003
    # Perimeter beams 1..4 + interior beam 5 connecting midpoints (5 -> 6)
    ops.element("elasticBeamColumn", 1, 1, 2, A, E, G, J, Iy, Iz, 1)
    ops.element("elasticBeamColumn", 2, 2, 3, A, E, G, J, Iy, Iz, 1)
    ops.element("elasticBeamColumn", 3, 3, 4, A, E, G, J, Iy, Iz, 1)
    ops.element("elasticBeamColumn", 4, 4, 1, A, E, G, J, Iy, Iz, 1)
    ops.element("elasticBeamColumn", 5, 5, 6, A, E, G, J, Iy, Iz, 1)

    q = 10.0
    ops.element(
        "LadrunoMembrane", 100,
        "-nodes", 1, 2, 3, 4,
        "-frames", 1, 2, 3, 4, 5,
        "-method", "tributary",
        "-pressure", q,
    )
    _solve_static(1.0)

    trib = ops.eleResponse(100, "tributaryAreas")
    assert len(trib) == 5
    assert math.isclose(sum(trib), 24.0, rel_tol=1e-5)
    # Interior beam 5 receives two-sided tributary load > each exterior short beam (2 or 4)
    assert trib[4] > trib[1] * 1.5
    _check_global_equilibrium_3d(pts, q, 24.0, 3.0, 2.0)


def test_support_nodes_mixed_3dof_and_6dof():
    """Flat slab supported on 4 corner columns (-supportNodes) with mixed 3-DOF and 6-DOF nodes."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 4.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 1, 1, 1, 1, 1, 1)

    # Switch to 3-DOF builder for nodes 3 and 4 to verify mixed-NDF handling
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(3, 4.0, 4.0, 0.0)
    ops.node(4, 0.0, 4.0, 0.0)
    ops.fix(3, 1, 1, 1)
    ops.fix(4, 1, 1, 1)

    q = 12.0
    ops.element(
        "LadrunoMembrane", 100,
        "-nodes", 1, 2, 3, 4,
        "-supportNodes", 1, 2, 3, 4,
        "-method", "tributary",
        "-pressure", q,
    )
    _solve_static(1.0)

    # Each corner point support in a 4x4 square receives 1/4 of -q * 16 = -48.0 kN
    for n in (1, 2, 3, 4):
        assert math.isclose(ops.nodeReaction(n, 3), -0.25 * q * 16.0, rel_tol=2e-3)

    nodal_resp = ops.eleResponse(100, "supportReactions")
    assert len(nodal_resp) == 4
    for r in nodal_resp:
        assert math.isclose(r, 0.25 * q * 16.0, rel_tol=2e-3)


# ---------------------------------------------------------------------------
# 7. Multi-pattern loading (Dead + loadConst + Live) & lazy frame resolution
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("beam_type", ["elastic", "force"])
def test_multi_pattern_dead_plus_live_no_double_counting(beam_type: str):
    Lx, Ly = 6.0, 4.0
    q_unit = 1.0
    coords = _build_rect_bay_3d(Lx=Lx, Ly=Ly, q=q_unit, method="tributary", beam_type=beam_type)

    # Pattern 1: Dead load = 10 kPa
    q_dead = 10.0
    ops.timeSeries("Linear", 1, "-factor", q_dead)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 100, "-type", "-surfaceLoad")

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.reactions()
    _check_global_equilibrium_3d(coords, q_dead, Lx * Ly, 0.5 * Lx, 0.5 * Ly)

    # Lock Dead load with loadConst and apply Pattern 2: Live load = 5 kPa
    ops.loadConst("-time", 0.0)
    q_live = 5.0
    ops.timeSeries("Linear", 2, "-factor", q_live)
    ops.pattern("Plain", 2, 2)
    ops.eleLoad("-ele", 100, "-type", "-surfaceLoad")
    assert ops.analyze(1) == 0
    ops.reactions()

    # Total reaction must equal -(q_dead + q_live) * A = -15 * 24 = -360 kN (never 2*q_dead + q_live)
    _check_global_equilibrium_3d(coords, q_dead + q_live, Lx * Ly, 0.5 * Lx, 0.5 * Ly)


def test_lazy_frame_resolution_membrane_declared_before_beams():
    Lx, Ly, q = 6.0, 4.0, 10.0
    coords = _build_rect_bay_3d(
        Lx=Lx,
        Ly=Ly,
        q=q,
        method="tributary",
        declare_membrane_first=True,
    )
    _solve_static(1.0)
    _check_global_equilibrium_3d(coords, q, Lx * Ly, 0.5 * Lx, 0.5 * Ly)


# ---------------------------------------------------------------------------
# 8. In-plane Guyan-condensed CST diaphragm, mass, rigid-body motion & Rayleigh
# ---------------------------------------------------------------------------

def test_in_plane_diaphragm_stiff_and_ndmaterial_equivalence():
    """In-plane shear displacement of a 4x4 diaphragm matches between -stiff and -thick/-mat."""
    E, nu, t = 25.0e6, 0.20, 0.15

    def _run_shear(use_ndmat: bool) -> tuple[float, float]:
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, 4.0, 0.0, 0.0)
        ops.node(3, 4.0, 4.0, 0.0)
        ops.node(4, 0.0, 4.0, 0.0)
        # Fix bottom edge (nodes 1, 2) completely; fix out-of-plane & rotations on (3, 4)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.fix(2, 1, 1, 1, 1, 1, 1)
        ops.fix(3, 0, 0, 1, 1, 1, 1)
        ops.fix(4, 0, 0, 1, 1, 1, 1)

        if use_ndmat:
            ops.nDMaterial("ElasticIsotropic", 10, E, nu)
            stiff_args = ["-thick", t, "-mat", 10]
        else:
            stiff_args = ["-stiff", E, nu, t]

        ops.element(
            "LadrunoMembrane", 100,
            "-nodes", 1, 2, 3, 4,
            "-supportNodes", 1, 2,
            *stiff_args,
        )
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(3, 500.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        ops.load(4, 500.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Linear")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        ops.reactions()
        return ops.nodeDisp(3, 1), ops.nodeReaction(1, 1) + ops.nodeReaction(2, 1)

    u_stiff, rx_stiff = _run_shear(use_ndmat=False)
    u_ndmat, rx_ndmat = _run_shear(use_ndmat=True)
    assert u_stiff > 0.0
    assert math.isclose(u_stiff, u_ndmat, rel_tol=1e-12)
    assert math.isclose(rx_stiff, -1000.0, rel_tol=1e-10)
    assert math.isclose(rx_ndmat, -1000.0, rel_tol=1e-10)


@pytest.mark.parametrize("lumped", [False, True])
def test_dynamic_ground_motion_rigid_body_probe_and_rayleigh(lumped: bool):
    """Free-sliding diaphragm under uniform ground acceleration a_g has a_rel == -a_g."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    Lx, Ly = 4.0, 3.0
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, Lx, 0.0, 0.0)
    ops.node(3, Lx, Ly, 0.0)
    ops.node(4, 0.0, Ly, 0.0)
    # Free in X (DOF 1), fixed in Y, Z and rotations so the diaphragm translates as a 1D rigid body
    for n in (1, 2, 3, 4):
        ops.fix(n, 0, 1, 1, 1, 1, 1)

    extra = ["-lumped"] if lumped else []
    ops.element(
        "LadrunoMembrane", 100,
        "-nodes", 1, 2, 3, 4,
        "-supportNodes", 1, 2,
        "-stiff", 25.0e6, 0.2, 0.20,
        "-rho", 2.5,
        *extra,
    )

    ag = 9.81
    ops.timeSeries("Constant", 1, "-factor", ag)
    ops.pattern("UniformExcitation", 1, 1, "-accel", 1)

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.algorithm("Linear")
    ops.analysis("Transient")
    assert ops.analyze(1, 0.01) == 0

    for n in (1, 2, 3, 4):
        assert math.isclose(ops.nodeAccel(n, 1), -ag, rel_tol=1e-6)


def test_rayleigh_damping_betaK_and_betaKc_damps_diaphragm():
    """Verify Element::getResistingForceIncInertia() applies Rayleigh damping (betaK, betaKc)."""
    def _run_free_vibration(beta_k: float) -> float:
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 6)
        ops.node(1, 0.0, 0.0, 0.0)
        ops.node(2, 4.0, 0.0, 0.0)
        ops.node(3, 4.0, 4.0, 0.0)
        ops.node(4, 0.0, 4.0, 0.0)
        ops.fix(1, 1, 1, 1, 1, 1, 1)
        ops.fix(2, 1, 1, 1, 1, 1, 1)
        ops.fix(3, 0, 1, 1, 1, 1, 1)
        ops.fix(4, 0, 1, 1, 1, 1, 1)

        ops.element(
            "LadrunoMembrane", 100,
            "-nodes", 1, 2, 3, 4,
            "-supportNodes", 1, 2,
            "-stiff", 20.0e4, 0.2, 0.10,
            "-rho", 2.5,
            "-lumped",
        )
        if beta_k > 0.0:
            ops.rayleigh(0.0, beta_k, 0.0, beta_k)

        ops.setNodeVel(3, 1, 1.0, "-commit")
        ops.setNodeVel(4, 1, 1.0, "-commit")

        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("Newmark", 0.5, 0.25)
        ops.algorithm("Linear")
        ops.analysis("Transient")
        peak = 0.0
        for step in range(25):
            assert ops.analyze(1, 0.005) == 0
            if step >= 15:
                peak = max(peak, abs(ops.nodeDisp(3, 1)))
        return peak

    undamped_tail = _run_free_vibration(0.0)
    damped_tail = _run_free_vibration(0.02)
    assert undamped_tail > 0.0
    assert damped_tail < 0.65 * undamped_tail


# ---------------------------------------------------------------------------
# 9. Database sendSelf / recvSelf round-trip
# ---------------------------------------------------------------------------

def test_send_recv_database_roundtrip():
    Lx, Ly, q = 6.0, 4.0, 10.0

    def _build():
        _build_rect_bay_3d(
            Lx=Lx,
            Ly=Ly,
            q=q,
            method="tributary",
            extra_membrane_args=["-stiff", 25.0e6, 0.2, 0.15, "-rho", 2.4],
        )
        # Release node 3 in X so we also probe non-zero in-plane nodal displacement
        ops.remove("sp", 3, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(3, 250.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        ops.system("BandGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Linear")
        ops.analysis("Static")
        assert ops.analyze(1) == 0

    def _probe():
        return list(ops.eleResponse(100, "tributaryAreas")) + list(ops.eleResponse(100, "stiffness"))[:12]

    database_roundtrip(
        _build,
        probe_nodes=[1, 2, 3, 4],
        ndf=6,
        probe_fn=_probe,
        rounds=2,
    )


def test_t2_apegmsh_subsegment_loads_and_diagram_reconstruction():
    """
    Verify all post-processor eleResponse queries (subSegmentLoads / beamSubLoads,
    subSegments, overhangTorques, meshStats, meshVertices, meshTriangles,
    nodalReactions, tributaryWidths) and verify that reconstructing internal
    force diagrams (N, Vy, Vz, My, Mz) from localForce + subSegmentLoads matches
    ForceBeamColumn3d internal Gauss-Lobatto section forces to < 1e-10.
    """
    import numpy as np

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)

    Lx, Ly, q = 6.0, 4.0, -12.0
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, Lx, 0.0, 0.0)
    ops.node(3, Lx, Ly, 0.0)
    ops.node(4, 0.0, Ly, 0.0)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1, 1, 1, 1)

    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    E, G, A, J, Iy, Iz = 30e6, 12.5e6, 0.18, 0.002, 0.003, 0.003
    ops.section("Elastic", 1, E, A, Iz, Iy, G, J)
    ops.beamIntegration("Lobatto", 1, 1, 5)

    ops.element("forceBeamColumn", 1, 1, 2, 1, 1)
    ops.element("forceBeamColumn", 2, 2, 3, 1, 1)
    ops.element("forceBeamColumn", 3, 3, 4, 1, 1)
    ops.element("forceBeamColumn", 4, 4, 1, 1, 1)

    ops.element(
        "LadrunoMembrane", 100,
        "-nodes", 1, 2, 3, 4,
        "-frames", 1, 2, 3, 4,
        "-method", "tributary",
        "-pressure", q,
    )

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 100, "-type", "-surfaceLoad")

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    # 1. Verify mesh telemetry and aliases
    stats = ops.eleResponse(100, "meshStats")
    assert len(stats) == 4
    n_verts, n_tris, n_sub, rank = (round(v) for v in stats)
    assert n_verts >= 5 and n_tris >= 4 and n_sub >= 4 and rank == 2

    sub_count = round(ops.eleResponse(100, "subSegments")[0])
    assert sub_count == n_sub

    verts3d = ops.eleResponse(100, "meshVertices")
    tris = ops.eleResponse(100, "meshTriangles")
    assert len(verts3d) == 3 * n_verts
    assert len(tris) == 3 * n_tris

    torques = ops.eleResponse(100, "overhangTorques")
    assert len(torques) == 4
    assert np.allclose(torques, 0.0, atol=1e-12)

    assert np.allclose(
        ops.eleResponse(100, "tributaryAreas"),
        ops.eleResponse(100, "tributaryWidths"),
        atol=1e-12,
    )
    assert np.allclose(
        ops.eleResponse(100, "supportReactions"),
        ops.eleResponse(100, "nodalReactions"),
        atol=1e-12,
    )

    # 2. Verify subSegmentLoads == beamSubLoads and reconstruct exact internal Vz(x), My(x)
    raw_sub = ops.eleResponse(100, "subSegmentLoads")
    assert len(raw_sub) == 9 * n_sub
    assert np.allclose(raw_sub, ops.eleResponse(100, "beamSubLoads"), atol=1e-12)

    # Lobatto 5-point stations on [0, Lx] for Beam 1
    lobatto_xi = np.array([0.0, (1.0 - np.sqrt(3.0 / 7.0)) / 2.0, 0.5, (1.0 + np.sqrt(3.0 / 7.0)) / 2.0, 1.0])
    x_pts = lobatto_xi * Lx

    lf = ops.eleResponse(1, "localForce")
    Vz_I, My_I = lf[2], lf[4]

    Vz_rec = np.full_like(x_pts, Vz_I)
    My_rec = -My_I - Vz_I * x_pts

    for s in range(n_sub):
        base = 9 * s
        f_tag = round(raw_sub[base + 0])
        if f_tag != 1:
            continue
        a = raw_sub[base + 1] * Lx
        b = raw_sub[base + 2] * Lx
        wza = raw_sub[base + 5]
        wzb = raw_sub[base + 6]
        c = b - a
        if c <= 1e-14:
            continue
        k = (wzb - wza) / c
        for idx, x in enumerate(x_pts):
            if x <= a:
                continue
            elif x < b:
                xi = x - a
                dF = wza * xi + 0.5 * k * (xi ** 2)
                dM = 0.5 * wza * (xi ** 2) + (1.0 / 6.0) * k * (xi ** 3)
            else:
                F_seg = 0.5 * (wza + wzb) * c
                Ma_seg = (c ** 2) * (wza + 2.0 * wzb) / 6.0
                dF = F_seg
                dM = F_seg * (x - a) - Ma_seg
            Vz_rec[idx] += dF
            My_rec[idx] -= dM

    # Check end consistency at x = Lx against Node J local forces
    assert abs(Vz_rec[-1] - (-lf[8])) < 1e-10
    assert abs(My_rec[-1] - lf[10]) < 1e-10

    # Check against ForceBeamColumn3d internal section forces [P, Mz, My, T] at all 5 Lobatto points
    for ip in range(1, 6):
        sec_f = ops.eleResponse(1, "section", ip, "force")
        assert abs(My_rec[ip - 1] - sec_f[2]) < 1e-9


@pytest.mark.parametrize("method", ["tributary", "oneway", "fea_poisson", "fea_plate"])
def test_t2_apegmsh_fem_layout_and_surface_load_flow(method: str):
    """
    Verify all 7 new apeGmsh internal FEM layout & surface load-distribution
    eleResponse queries: frameTags, supportNodeTags, meshSubSegments,
    meshVertexFlags, meshNodalField, meshTriangleLoadFlow, and yieldLines.
    """
    import numpy as np

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)

    Lx, Ly, q = 6.0, 4.0, -10.0
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, Lx, 0.0, 0.0)
    ops.node(3, Lx, Ly, 0.0)
    ops.node(4, 0.0, Ly, 0.0)
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1, 1, 1, 1)

    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    E, G, A, J, Iy, Iz = 30e6, 12.5e6, 0.18, 0.002, 0.003, 0.003
    ops.element("elasticBeamColumn", 10, 1, 2, A, E, G, J, Iy, Iz, 1)
    ops.element("elasticBeamColumn", 20, 2, 3, A, E, G, J, Iy, Iz, 1)
    ops.element("elasticBeamColumn", 30, 3, 4, A, E, G, J, Iy, Iz, 1)
    ops.element("elasticBeamColumn", 40, 4, 1, A, E, G, J, Iy, Iz, 1)

    extra = ["-spanAxis", 1.0, 0.0, 0.0] if method == "oneway" else []
    ops.element(
        "LadrunoMembrane", 100,
        "-nodes", 1, 2, 3, 4,
        "-frames", 10, 20, 30, 40,
        "-method", method,
        "-pressure", q,
        *extra,
    )

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.eleLoad("-ele", 100, "-type", "-surfaceLoad")

    ops.system("BandGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    # 1. frameTags and supportNodeTags (empty Vector returns [0.0] sentinel)
    frame_tags = [round(v) for v in ops.eleResponse(100, "frameTags") if round(v) > 0]
    assert frame_tags == [10, 20, 30, 40]
    supp_tags = [round(v) for v in ops.eleResponse(100, "supportNodeTags") if round(v) > 0]
    assert supp_tags == []

    # 2. meshVertices, meshTriangles, meshVertexFlags, meshNodalField
    verts = np.array(ops.eleResponse(100, "meshVertices"), dtype=float).reshape(-1, 3)
    tris = np.array(ops.eleResponse(100, "meshTriangles"), dtype=float).round().astype(int).reshape(-1, 3)
    vflags = np.array(ops.eleResponse(100, "meshVertexFlags"), dtype=float).round().astype(int)
    field = np.array(ops.eleResponse(100, "meshNodalField"), dtype=float)
    assert len(vflags) == len(verts)
    assert len(field) == len(verts)
    assert np.allclose(field, ops.eleResponse(100, "meshField"), atol=1e-12)
    # Every vertex flag must be in {0, 1, 2, 3, 4}, with constrained frame vertices (3) and interior Steiner vertices (0)
    assert set(vflags.tolist()).issubset({0, 1, 2, 3, 4})
    assert np.any(vflags == 3)
    assert np.any(vflags == 0)
    # Constrained frame vertices have zero deflection/distance on active supports
    if method == "oneway":
        on_active_supports = (vflags == 3) & (np.isclose(verts[:, 0], 0.0) | np.isclose(verts[:, 0], Lx))
        assert np.allclose(field[on_active_supports], 0.0, atol=1e-10)
    else:
        assert np.allclose(field[vflags == 3], 0.0, atol=1e-10)
    assert np.max(np.abs(field[vflags == 0])) > 0.0

    # 3. meshSubSegments: [uIdx, vIdx, frameTag, aOverL, bOverL, w_avg]
    sub_segs = np.array(ops.eleResponse(100, "meshSubSegments"), dtype=float).reshape(-1, 6)
    assert len(sub_segs) >= 4
    total_sub_force = 0.0
    for row in sub_segs:
        u_idx, v_idx, f_tag = round(row[0]), round(row[1]), round(row[2])
        a_l, b_l, w_avg = row[3], row[4], row[5]
        assert 0 <= u_idx < len(verts)
        assert 0 <= v_idx < len(verts)
        assert f_tag in (10, 20, 30, 40)
        assert 0.0 <= a_l <= b_l <= 1.0 + 1e-12
        seg_len = float(np.linalg.norm(verts[v_idx] - verts[u_idx]))
        total_sub_force += w_avg * seg_len
    # Integrated w_avg * L_seg across all constrained frame sub-segments equals |q| * Lx * Ly
    assert math.isclose(abs(total_sub_force), abs(q) * Lx * Ly, rel_tol=1e-6)

    # 4. meshTriangleLoadFlow: [targetTag, qx, qy, qz] per triangle
    tri_flow = np.array(ops.eleResponse(100, "meshTriangleLoadFlow"), dtype=float).reshape(-1, 4)
    assert len(tri_flow) == len(tris)
    assert np.allclose(
        tri_flow.ravel(),
        ops.eleResponse(100, "meshLoadVectors"),
        atol=1e-12,
    )
    target_tags = np.round(tri_flow[:, 0]).astype(int)
    assert set(target_tags.tolist()).issubset({10, 20, 30, 40})
    if method == "oneway":
        # One-way along X transfers load exclusively to the Y-parallel beams (20 and 40)
        assert set(target_tags.tolist()) == {20, 40}
    else:
        assert set(target_tags.tolist()) == {10, 20, 30, 40}
    # In-plane load-flow vectors have non-zero magnitude
    flow_norms = np.linalg.norm(tri_flow[:, 1:4], axis=1)
    assert np.max(flow_norms) > 0.0

    # 5. yieldLines: [X1, Y1, Z1, X2, Y2, Z2, tagA, tagB]
    yl = np.array(ops.eleResponse(100, "yieldLines"), dtype=float).reshape(-1, 8)
    assert len(yl) >= 1
    for row in yl:
        p1, p2 = row[0:3], row[3:6]
        tag_a, tag_b = round(row[6]), round(row[7])
        assert tag_a != tag_b
        assert tag_a in (10, 20, 30, 40)
        assert tag_b in (10, 20, 30, 40)
        assert float(np.linalg.norm(p2 - p1)) > 0.0


@pytest.mark.parametrize("transf_type", ["Linear", "PDelta", "Corotational"])
def test_wp209_timoshenko_local_axes_and_etabs_horizontal_membrane_projection(transf_type: str):
    """WP-209: ElasticTimoshenkoBeam3d exposes 'localAxes' (case 30) matching elasticBeamColumn,
    so a horizontal LadrunoMembrane with ETABS-oriented horizontal beams (vecxz = (0, -1, 0) along X,
    (1, 0, 0) along Y) projects vertical gravity onto local +y (ETABS local 2 = +Z) instead of
    falling back to vxz = (0, 0, 1) and loading the beams sideways."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)

    Lx, Ly, q = 6.0, 4.0, 12.0
    coords = {
        1: (0.0, 0.0, 0.0),
        2: (Lx, 0.0, 0.0),
        3: (Lx, Ly, 0.0),
        4: (0.0, Ly, 0.0),
    }
    for tag, (x, y, z) in coords.items():
        ops.node(tag, x, y, z)
        ops.fix(tag, 1, 1, 1, 1, 1, 1)

    # Standard ETABS horizontal beam local axes: local 2 = +Z => OpenSees local y = +Z
    # X-beams (+X): vecxz = (0, -1, 0) => local z = -Y, local y = (-Y) x (+X) = +Z
    # Y-beams (+Y): vecxz = (1, 0, 0)  => local z = +X, local y = (+X) x (+Y) = +Z
    ops.geomTransf(transf_type, 1, 0.0, -1.0, 0.0)
    ops.geomTransf(transf_type, 2, 1.0, 0.0, 0.0)
    ops.geomTransf(transf_type, 3, 0.0, 1.0, 0.0)
    ops.geomTransf(transf_type, 4, -1.0, 0.0, 0.0)

    E, G, A, J, Iy, Iz, Avy, Avz = 30.0e6, 12.5e6, 0.18, 0.002, 0.001, 0.004, 0.15, 0.15
    edges = [(1, 1, 2, 1), (2, 2, 3, 2), (3, 3, 4, 3), (4, 4, 1, 4)]
    for etag, ni, nj, ttag in edges:
        ops.element("ElasticTimoshenkoBeam", etag, ni, nj, E, G, A, J, Iy, Iz, Avy, Avz, ttag)
        ops.element("elasticBeamColumn", etag + 10, ni, nj, A, E, G, J, Iy, Iz, ttag)
        la_timo = ops.eleResponse(etag, "localAxes")
        la_eb = ops.eleResponse(etag + 10, "localAxes")
        assert len(la_timo) == 9
        assert la_timo == pytest.approx(la_eb, rel=1e-12, abs=1e-12)
        # Local y must point in +Z (0, 0, 1) for every horizontal beam
        assert la_timo[3:6] == pytest.approx([0.0, 0.0, 1.0], abs=1e-12)
        ops.remove("element", etag + 10)

    ops.element(
        "LadrunoMembrane", 100,
        "-nodes", 1, 2, 3, 4,
        "-frames", 1, 2, 3, 4,
        "-method", "tributary",
        "-pressure", q,
    )
    _solve_static(1.0)
    _check_global_equilibrium_3d(coords, q, Lx * Ly, 0.5 * Lx, 0.5 * Ly)

    # All bending must be in major axis Mz (localForce[5], localForce[11]) and shear in Vy (localForce[1], localForce[7]);
    # minor-axis Vz (2, 8) and My (4, 10) must be identically zero.
    for etag, _, _, _ in edges:
        lf = ops.eleResponse(etag, "localForce")
        assert abs(lf[1]) > 1.0 and abs(lf[7]) > 1.0
        assert abs(lf[5]) > 1.0 and abs(lf[11]) > 1.0
        assert lf[2] == pytest.approx(0.0, abs=1e-10)
        assert lf[8] == pytest.approx(0.0, abs=1e-10)
        assert lf[4] == pytest.approx(0.0, abs=1e-10)
        assert lf[10] == pytest.approx(0.0, abs=1e-10)


