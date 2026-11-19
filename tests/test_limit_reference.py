"""WP-237b: tests of the Python reference assembler ``tests/limit_reference.py``.

L0 operator tests need numpy and scipy only (Zone A). They check every operator against an independent
route: the contract's helpers, Gauss quadrature, closed-form fields, or the case derivations of interface
16, and they check that each family of program rows is live (a field breaking it leaves a residual).

Tests that run the open solver are Zone B (``pytest.importorskip("clarabel")``); their numbers are
estimates of the discrete optima, never certified values. The ``slow`` tier holds the larger estimates:
measured 2026-11-19 on a 16-thread desktop, 52 s in total (Clarabel 0.11.1); the hard stop is
``limit_reference.CLARABEL_TIME_LIMIT_S``, and a leg ending on it carries no number and fails.
"""
import math

import numpy as np
import pytest
import scipy.sparse as sp

import limit_reference as lr
import test_limit_contract as contract


S2, S3 = math.sqrt(2.0), math.sqrt(3.0)
CASES = contract.parse_cases()
EXACT = sorted(name for name, case in CASES.items()
               if case["status"] == "active" and case["class"] == "exact_on_mesh")
EXACT_LOWER = sorted(name for name, case in CASES.items()
                     if case["status"] == "active" and case["class"] == "exact_lower_on_mesh")
BUILDABLE = sorted(name for name, case in CASES.items() if case["status"] == "active" and case["boundaries"])
FORMULATIONS = (lr.LOWER_ID, lr.CONTROL_ID, lr.JUMP_ID)

_CACHE = {}


def problem_of(name, mesh=0):
    key = (name, mesh)
    if key not in _CACHE:
        problem, arrays = lr.case_problem(name, mesh)
        _CACHE[key] = (problem, arrays, lr.Geometry(problem, arrays))
    return _CACHE[key]


def random_tet(rng):
    while True:
        x = rng.normal(size=(4, 3))
        if np.linalg.det(x[1:] - x[0]) > 0.3:
            return x


def single_tet_problem(x):
    """A one-element problem around four vertices, with every facet free."""
    nodes, tets = contract.tet10_from_tet4(np.asarray(x, float), np.array([[0, 1, 2, 3]]))
    problem = {"regions": [{"name": "r", "criterion": "von_mises", "k": 1.0, "b0": [0.0] * 3, "b1": [0.0] * 3}],
               "footings": [], "boundary_sets": [{"name": "all", "kind": "free"}], "queries": []}
    arrays = {"nodes": nodes, "tets": tets, "element_region": np.zeros(1, dtype=np.int64),
              "facets": np.array([[0, f] for f in range(4)], dtype=np.int64), "facet_set": np.zeros(4, np.int64),
              "t0": np.zeros((4, 3, 3)), "t1": np.zeros((4, 3, 3))}
    return problem, arrays


def gauss_triangle():
    return contract.gauss_triangle()


def gauss_tet():
    return contract.gauss_tet()


# ---------------------------------------------------------------------------------------------------
# Tables and conventions
# ---------------------------------------------------------------------------------------------------

@pytest.mark.zone_a
def test_tables_and_constants_match_the_contract():
    np.testing.assert_array_equal(lr.Q, contract.Q)
    np.testing.assert_array_equal(lr.ONE, contract.ONE)
    assert lr.EDGES == contract.EDGES
    assert lr.FACET_VERTICES == contract.FACET_VERTICES
    assert lr.FACET_NODES == contract.FACET_NODES
    np.testing.assert_allclose(lr.FACET_WORK, contract.facet_work_table(1.0), rtol=0, atol=1e-16)
    rng = np.random.default_rng(1)
    nodal = rng.normal(size=(6, 3))
    np.testing.assert_allclose(lr.BERNSTEIN @ nodal, contract.bernstein_from_nodal(nodal), atol=1e-15)
    for _ in range(20):
        n = rng.normal(size=3)
        n /= np.linalg.norm(n)
        for mine, theirs in zip(lr.tangent_basis(n), contract.tangent_basis(n)):
            np.testing.assert_allclose(mine, theirs, atol=1e-15)
        np.testing.assert_allclose(lr.tangent_bases(n[None])[0], np.vstack(lr.tangent_basis(n)), atol=1e-15)
    # The tie rule: first index of the smallest |n_m|.
    t1, t2 = lr.tangent_basis(np.array([0.0, 0.0, 1.0]))
    np.testing.assert_array_equal(t1, [0.0, -1.0, 0.0])     # e_x x e_z
    np.testing.assert_array_equal(t2, [1.0, 0.0, 0.0])      # e_z x T1
    for _ in range(5):
        n = rng.normal(size=3)
        np.testing.assert_array_equal(lr.voigt_n(n), contract.voigt_n(n))


@pytest.mark.zone_a
def test_volume_and_facet_weights_are_the_exact_integrals():
    weights = np.zeros(10)
    for bary, weight in gauss_tet():
        weights += weight * contract.p2_values(bary)
    np.testing.assert_allclose(lr.BODY_WEIGHTS, weights, atol=1e-13)
    numeric = np.zeros((3, 6))
    for bary, weight in gauss_triangle():
        shape = [bary[j] * (2.0 * bary[j] - 1.0) for j in range(3)]
        shape += [4.0 * bary[a] * bary[b] for a, b in lr.FACET_EDGES]
        numeric += weight * np.outer(bary, shape)
    np.testing.assert_allclose(lr.FACET_WORK, numeric, atol=1e-13)


@pytest.mark.zone_a
def test_support_value_and_jump_constants():
    rng = np.random.default_rng(2)
    for alpha in (0.0, 0.1, 0.25):
        for _ in range(10):
            d = rng.normal(size=6)
            if alpha == 0.0:
                d[:3] -= d[:3].mean()
            else:
                theta, enorm = lr.strain_split(d)
                d[:3] += (3.0 * S2 * alpha * enorm - theta + 0.3) / 3.0
            assert lr.support_value(d, 1.3, alpha) == pytest.approx(contract.support_function(d, 1.3, alpha), rel=1e-12)
    tau, kj = lr.jump_constants(np.array([0.0, 0.2]))
    ctau, croot = contract.jump_constants(0.2)
    assert tau[0] == 0.0 and kj[0] == 1.0
    assert tau[1] == pytest.approx(ctau, rel=1e-15) and kj[1] == pytest.approx(1.0 / croot, rel=1e-15)
    with pytest.raises(contract.ContractError) as error:
        lr.jump_constants(1.0 / math.sqrt(12.0))
    assert error.value.code == "formulation_not_applicable"


# ---------------------------------------------------------------------------------------------------
# Geometry and kinematic operators (L0)
# ---------------------------------------------------------------------------------------------------

@pytest.mark.zone_a
def test_gradients_agree_with_both_routes_and_normals_point_out():
    rng = np.random.default_rng(3)
    for _ in range(20):
        x = random_tet(rng)
        geo = lr.Geometry(*single_tet_problem(x))
        grads, volume = contract.tet_gradients(x)
        np.testing.assert_allclose(geo.grad[0], grads, atol=1e-12)
        assert geo.volume[0] == pytest.approx(volume, rel=1e-13)
        for i in range(4):
            n, area = geo.facet_normal[0, i], geo.facet_area[0, i]
            np.testing.assert_allclose(geo.grad[0, i], -(area / (3.0 * volume)) * n, atol=1e-12)
            assert n @ (x[i] - x[lr.FACET_VERTICES[i][0]]) < 0.0


@pytest.mark.zone_a
def test_rigid_velocity_gives_zero_strain_rate_and_linear_velocity_its_symmetric_gradient():
    rng = np.random.default_rng(4)
    for _ in range(10):
        x = random_tet(rng)
        geo = lr.Geometry(*single_tet_problem(x))
        xn = geo.xn[0]
        u, w = rng.normal(size=3), rng.normal(size=3)
        rigid = u + np.cross(w, xn)
        assert np.abs(geo.strain_rates(rigid[None])).max() < 1e-12
        g = rng.normal(size=(3, 3))
        linear = xn @ g.T + u
        sym = 0.5 * (g + g.T)
        expected = np.array([sym[0, 0], sym[1, 1], sym[2, 2], 2 * sym[0, 1], 2 * sym[1, 2], 2 * sym[2, 0]])
        np.testing.assert_allclose(geo.strain_rates(linear[None])[0], np.tile(expected, (4, 1)), atol=1e-12)


@pytest.mark.zone_a
def test_quadratic_velocity_strain_matches_the_analytic_gradient_and_the_contract():
    rng = np.random.default_rng(5)
    for _ in range(10):
        x = random_tet(rng)
        geo = lr.Geometry(*single_tet_problem(x))
        a, g, h = rng.normal(size=3), rng.normal(size=(3, 3)), rng.normal(size=(3, 3, 3))
        h = 0.5 * (h + np.transpose(h, (0, 2, 1)))

        def field(p):
            return a + g @ p + np.einsum("cij,i,j->c", h, p, p)

        def gradient(p):
            return g + 2.0 * np.einsum("cij,j->ci", h, p)

        v = np.array([field(p) for p in geo.xn[0]])
        d = geo.strain_rates(v[None])[0]
        for m in range(4):
            grad = gradient(x[m])
            sym = 0.5 * (grad + grad.T)
            expected = [sym[0, 0], sym[1, 1], sym[2, 2], 2 * sym[0, 1], 2 * sym[1, 2], 2 * sym[2, 0]]
            np.testing.assert_allclose(d[m], expected, atol=1e-10)
            np.testing.assert_allclose(d[m], contract.strain_rate_at_vertex(geo.xn[0], v, m), atol=1e-10)


@pytest.mark.zone_a
def test_mesh_facets_match_the_contract():
    for name, mesh in (("vm_plate_shear", 1), ("tresca_strip_vm_thin_slab", 0), ("upper3d_fixed_phi00", 0)):
        problem, arrays, geo = problem_of(name, mesh)
        interior, exterior = contract.mesh_facets(arrays["tets"])
        got = [((int(a), int(b)), (int(c), int(d))) for a, b, c, d in
               zip(geo.e_minus, geo.f_minus, geo.e_plus, geo.f_plus)]
        assert got == interior
        assert [tuple(int(v) for v in row) for row in geo.facets] == exterior
        tets = arrays["tets"]
        for fi in range(0, geo.F, max(1, geo.F // 50)):
            np.testing.assert_array_equal(tets[geo.e_minus[fi], geo.minus_nodes[fi]],
                                          tets[geo.e_plus[fi], geo.plus_nodes[fi]])


@pytest.mark.zone_a
def test_geometry_refuses_curved_inverted_and_reserved_input():
    problem, arrays, _ = problem_of("vm_plate_shear")
    bent = {key: value.copy() for key, value in arrays.items()}
    bent["nodes"][bent["tets"][0, 4]] += 1e-3
    with pytest.raises(contract.ContractError) as error:
        lr.Geometry(problem, bent)
    assert error.value.code == "curved_element"
    flipped = {key: value.copy() for key, value in arrays.items()}
    flipped["tets"][0, [1, 2]] = flipped["tets"][0, [2, 1]]
    with pytest.raises(contract.ContractError) as error:
        lr.Geometry(problem, flipped)
    assert error.value.code == "element_orientation"
    for path, value, code in (
        (("regions", 0, "criterion"), "tresca", "criterion_reserved"),
        (("footings", 0, "interface_law"), "coulomb", "interface_law_reserved"),
        (("boundary_sets", 0, "kind"), "interface", "boundary_kind_reserved"),
    ):
        changed = contract.copy.deepcopy(problem)
        changed[path[0]][path[1]][path[2]] = value
        with pytest.raises(contract.ContractError) as error:
            lr.Geometry(changed, arrays)
        assert error.value.code == code


@pytest.mark.zone_a
def test_formulations_are_named_and_the_jump_refuses_steep_cones():
    problem, arrays, geo = problem_of("vm_plate_shear")
    with pytest.raises(contract.ContractError) as error:
        lr.assemble(problem, arrays, "upper.p2d-jump-other.v0", geometry=geo)
    assert error.value.code == "formulation_not_admitted"
    steep = contract.copy.deepcopy(problem)
    steep["regions"][0] = {"name": "soil", "criterion": "drucker_prager", "k": 1.0, "alpha": 0.3,
                           "b0": [0.0] * 3, "b1": [0.0] * 3}
    lr.assemble(steep, arrays, lr.CONTROL_ID)
    with pytest.raises(contract.ContractError) as error:
        lr.assemble(steep, arrays, lr.JUMP_ID)
    assert error.value.code == "formulation_not_applicable"


# ---------------------------------------------------------------------------------------------------
# Loads and resultants (L0)
# ---------------------------------------------------------------------------------------------------

@pytest.mark.zone_a
def test_rigid_motion_gravity_work_is_exact():
    case = CASES["gravity_work_rigid_motion"]
    probe, region = case["probe"], case["regions"][0]
    u, w, about = (np.array(probe[key]) for key in ("u", "w", "about"))
    for mesh in range(len(case["mesh"])):
        vertices, tet4 = contract.case_tet4(case, mesh)
        nodes, tets = contract.tet10_from_tet4(vertices, tet4)
        _, exterior = contract.mesh_facets(tets)
        problem = {"regions": [region], "footings": [], "boundary_sets": [{"name": "all", "kind": "free"}],
                   "queries": []}
        arrays = {"nodes": nodes, "tets": tets, "element_region": np.zeros(len(tets), np.int64),
                  "facets": np.array(exterior, np.int64), "facet_set": np.zeros(len(exterior), np.int64),
                  "t0": np.zeros((len(exterior), 3, 3)), "t1": np.zeros((len(exterior), 3, 3))}
        geo = lr.Geometry(problem, arrays)
        field = u + np.cross(w, geo.xn - about)
        work = float(np.sum(lr.element_node_work(geo, 0) * field))
        assert work == pytest.approx(case["target"], rel=1e-12)
        # A wrong (lumped vertex) rule misses a quadratic field.
        z = geo.xn[..., 2:3] ** 2 * np.array([0.0, 0.0, 1.0])
        exact = float(np.sum(lr.element_node_work(geo, 0) * z))
        lumped = float(np.sum((geo.volume[:, None] * np.array([0.25] * 4 + [0.0] * 6))[:, :, None]
                              * geo.b0[:, None, :] * z))
        box = case["box"]
        analytic = region["b0"][2] * (box[3] - box[0]) * (box[4] - box[1]) * (box[5] ** 3 - box[2] ** 3) / 3.0
        assert exact == pytest.approx(analytic, rel=1e-12)
        assert abs(lumped - analytic) > 1e-3 * abs(analytic)


@pytest.mark.zone_a
def test_traction_work_is_exact_for_linear_traction_and_quadratic_velocity():
    rng = np.random.default_rng(6)
    x = random_tet(rng)
    problem, arrays = single_tet_problem(x)
    problem["boundary_sets"] = [{"name": "all", "kind": "traction"}]
    arrays["t0"] = rng.normal(size=(4, 3, 3))
    geo = lr.Geometry(problem, arrays)
    v = rng.normal(size=(10, 3))
    work = float(np.sum(lr.element_node_work(geo, 0)[0] * v))
    numeric = 0.0
    for f in range(4):
        verts, nodes6 = lr.FACET_VERTICES[f], lr.FACET_NODES[f]
        area = geo.facet_area[0, f]
        for bary, weight in gauss_triangle():
            t = bary @ arrays["t0"][f]
            shape = [bary[j] * (2.0 * bary[j] - 1.0) for j in range(3)]
            shape += [4.0 * bary[a] * bary[b] for a, b in lr.FACET_EDGES]
            numeric += weight * area * t @ (np.asarray(shape) @ v[list(nodes6)])
        assert verts == tuple(nodes6[:3])
    assert work == pytest.approx(numeric, rel=1e-12)


@pytest.mark.zone_a
@pytest.mark.parametrize("name", ["vm_plate_shear", "vm_plate_shear_offset_ref", "tresca_strip_vm_thin_slab"])
def test_footing_resultant_matches_the_contract_and_exact_integration(name):
    problem, arrays, geo = problem_of(name)
    rng = np.random.default_rng(7)
    sigma = rng.normal(size=(geo.E, 4, 6))
    np.testing.assert_allclose(lr.lower_resultant(problem, arrays, sigma, geo),
                               contract.fixture_resultant(problem, arrays, sigma), atol=1e-12)
    # Exact for a linear traction: compare one footing facet with Gauss quadrature.
    row = next(r for r in range(len(geo.facets)) if geo.kinds[geo.facet_set[r]] == "footing")
    e, f = geo.facets[row]
    x_ref = np.array(problem["footings"][0]["x_ref"])
    verts = list(lr.FACET_VERTICES[f])
    rmat = lr.footing_resultant_matrix(geo).toarray()
    n = geo.boundary_normal[row]
    area = geo.boundary_area[row]
    traction = np.array([lr.voigt_n(n) @ sigma[e, a] for a in verts])
    force, moment = np.zeros(3), np.zeros(3)
    for bary, weight in gauss_triangle():
        point, t = bary @ geo.xv[e, verts], bary @ traction
        force += weight * area * t
        moment += weight * area * np.cross(point - x_ref, t)
    columns = np.concatenate([24 * e + 6 * a + np.arange(6) for a in verts])
    # Only this facet's contribution: the columns of element e's facet vertices, minus other footing
    # facets of the same element (none on these meshes: a footing facet's element has one footing facet).
    others = [r for r in range(len(geo.facets)) if r != row and geo.facets[r, 0] == e
              and geo.kinds[geo.facet_set[r]] == "footing"]
    assert not others
    got = rmat[:6, columns] @ np.concatenate([sigma[e, a] for a in verts])
    np.testing.assert_allclose(got, np.concatenate([force, moment]), atol=1e-12)


# ---------------------------------------------------------------------------------------------------
# Exact fields in the programs (L0): rows hold, cones hold, the objective is the case value
# ---------------------------------------------------------------------------------------------------

def exact_lower_field(name, geo):
    """(lambda, sigma) of interface 16 for the exact and sandwich cases."""
    case = CASES[name]
    x = geo.xv
    sigma = np.zeros((geo.E, 4, 6))
    target = case.get("target")
    if name in ("vm_uniaxial_tension", "dp_uniaxial_tension"):
        sigma[..., 0] = target
    elif name == "dp_uniaxial_compression":
        sigma[..., 0] = -target
    elif name == "dp_hydrostatic_tension":
        sigma[..., :3] = target
    elif name == "vm_tension_dead_load":
        sigma[..., 0] = S3
    elif name in ("vm_plate_shear", "vm_plate_shear_offset_ref"):
        sigma[..., 3] = 1.0
    elif name == "vm_body_force_shear_block":
        sigma[..., 3] = target * (1.0 - x[..., 1])
    elif name == "vm_dead_gravity_shear_layer":
        sigma[..., 5] = target + 0.3 * (1.0 - x[..., 2])
    elif name == "dp_weighted_apex_column":
        sigma[..., :3] = (target + 0.2 * x[..., 2])[..., None]
    elif name == "vm_weighted_column_sandwich":
        target = contract.TARGET_FORMULAS["column_static"]({"k": 1.0, "gamma": 0.25, "H": 2.0})
        sigma[..., 2] = -target - 0.25 * (2.0 - x[..., 2])
    else:
        raise KeyError(name)
    return target, sigma


def lower_vector(program, lam, sigma):
    x = np.zeros(program.A.shape[1])
    if program.layout["lambda"] is not None:
        x[program.layout["lambda"]] = lam
    x[program.layout["sigma"]:] = sigma.reshape(-1)
    return x


def in_cone(slack, tol):
    return bool(np.all(slack[:, 0] >= np.linalg.norm(slack[:, 1:], axis=1) - tol))


LOWER_FIELD_CASES = EXACT + EXACT_LOWER + ["vm_weighted_column_sandwich"]


@pytest.mark.zone_a
@pytest.mark.parametrize("name", LOWER_FIELD_CASES)
@pytest.mark.parametrize("mesh", [0, 1])
def test_exact_lower_field_satisfies_the_lower_program(name, mesh):
    problem, arrays, geo = problem_of(name, mesh)
    program = lr.assemble(problem, arrays, lr.LOWER_ID, geometry=geo)
    lam, sigma = exact_lower_field(name, geo)
    x = lower_vector(program, lam, sigma)
    assert np.abs(program.zero_residual(x)).max() < 1e-12
    (cones,) = program.cone_slacks(x)
    assert in_cone(cones, 1e-12)
    assert program.value(x) == pytest.approx(lam, rel=1e-15)
    if name != "vm_weighted_column_sandwich":
        # At the exact optimum some vertex is on the yield surface.
        assert np.max(lr.yield_margin(problem, arrays, sigma, geo)) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.zone_a
def test_lower_rows_are_live():
    """Zeroed gravity, a facet traction jump, a sign-flipped resultant and a weakened cone each show."""
    problem, arrays, geo = problem_of("vm_dead_gravity_shear_layer", 1)
    program = lr.assemble(problem, arrays, lr.LOWER_ID, geometry=geo)
    lam, sigma = exact_lower_field("vm_dead_gravity_shear_layer", geo)
    zero_gravity = sigma.copy()
    zero_gravity[..., 5] = lam
    assert np.abs(program.zero_residual(lower_vector(program, lam, zero_gravity))[: 3 * geo.E]).max() > 0.1
    jump = sigma.copy()
    jump[0, :, 0] += 1e-2
    residual = program.zero_residual(lower_vector(program, lam, jump))
    assert np.abs(residual[3 * geo.E: 3 * geo.E + 9 * geo.F]).max() > 1e-3
    problem, arrays, geo = problem_of("vm_plate_shear_offset_ref", 0)
    program = lr.assemble(problem, arrays, lr.LOWER_ID, geometry=geo)
    lam, sigma = exact_lower_field("vm_plate_shear_offset_ref", geo)
    flipped = contract.copy.deepcopy(problem)
    flipped["queries"][0]["F1"][0][5] *= -1.0
    flipped["problem_sha256"] = contract.document_hash(flipped, arrays, "problem_sha256")
    wrong = lr.assemble(flipped, arrays, lr.LOWER_ID, geometry=geo)
    assert np.abs(program.zero_residual(lower_vector(program, lam, sigma))).max() < 1e-12
    assert np.abs(wrong.zero_residual(lower_vector(wrong, lam, sigma))).max() > 0.5
    (cones,) = program.cone_slacks(lower_vector(program, lam, 1.01 * sigma))
    assert not in_cone(cones, 1e-6)


def upper_vector(program, velocity, footing_velocity):
    """A program vector for a continuous P2 field ``velocity`` (N, 3), with cone variables at their optimum."""
    layout, geo = program.layout, program.layout["geometry"]
    x = np.zeros(program.A.shape[1])
    field = velocity[geo.tets]
    if layout["broken"]:
        x[: 30 * geo.E] = field.reshape(-1)
    else:
        x[: 3 * geo.N] = velocity.reshape(-1)
    x[layout["V"]: layout["V"] + footing_velocity.size] = footing_velocity.reshape(-1)
    d = geo.strain_rates(field)
    z = np.linalg.solve(lr.Q @ lr.Q.T, np.einsum("ab,eib->eia", lr.Q, d)[..., None])[..., 0]
    theta = d[..., 0] + d[..., 1] + d[..., 2]
    alpha = geo.alpha[:, None]
    r = np.where(alpha > 0.0, theta / (3.0 * S2 * np.where(alpha > 0.0, alpha, 1.0)), np.linalg.norm(z, axis=-1))
    x[layout["rz"]: layout["rz"] + 24 * geo.E] = np.concatenate([r[..., None], z], axis=-1).reshape(-1)
    if layout["broken"] and geo.F:
        coeff = lr.jump_coefficients(geo, field)
        basis = np.concatenate([geo.interior_normal[:, None], lr.tangent_bases(geo.interior_normal)], axis=1)
        local = np.einsum("fqc,fbc->fbq", basis, coeff)
        tau, _ = lr.jump_constants(geo.alpha[geo.e_minus])
        rho = np.where(tau[:, None] > 0, local[..., 0] / np.where(tau > 0, tau, 1.0)[:, None],
                       np.linalg.norm(local[..., 1:], axis=-1))
        x[layout["jump"]:] = np.concatenate([rho[..., None], local[..., 1:]], axis=-1).reshape(-1)
    return x


def exact_upper_field(name, geo, problem):
    """(value, velocity (N, 3), footing velocity) of interface 16 for the exact cases, with W1 = 1."""
    case = CASES[name]
    xg = geo.global_node_coordinates()
    n_f = len(problem["footings"])
    footing = np.zeros((n_f, 6))
    target = case.get("target")
    k = problem["regions"][0]["k"]
    alpha = problem["regions"][0].get("alpha", 0.0)
    if name in ("vm_uniaxial_tension", "vm_tension_dead_load"):
        v = xg * np.array([1.0, -0.5, -0.5])
    elif name == "dp_uniaxial_tension":
        a, c = 1.0 / S3 + alpha, -0.5 / S3 + alpha
        v = xg * np.array([1.0, c / a, c / a])
    elif name == "dp_uniaxial_compression":
        a, c = -1.0 / S3 + alpha, 0.5 / S3 + alpha
        v = -xg * np.array([1.0, c / a, c / a])
    elif name == "dp_hydrostatic_tension":
        v = xg / 3.0
    elif name in ("vm_plate_shear", "vm_plate_shear_offset_ref"):
        v = np.zeros_like(xg)
        v[:, 0] = xg[:, 1]
        footing[0, 0] = 1.0
    else:
        raise KeyError(name)
    assert target is not None and k > 0
    return target, v, footing


def band_upper_field(name, geo, case, mesh):
    """The one-layer band fields of interface 16.3 and 16.5, unnormalised."""
    assertion = next(tokens for tokens in case["assert"] if tokens[0] == "upper_le_formula")
    params = contract.parse_key_values(assertion[2:])
    axis = int(params["axis"])
    h = 1.0 / case["mesh"][mesh][1][axis]
    xg = geo.global_node_coordinates()
    v = np.zeros_like(xg)
    if name == "vm_body_force_shear_block":
        v[:, 0] = np.minimum(xg[:, 1] / h, 1.0)
    elif name == "vm_dead_gravity_shear_layer":
        v[:, 0] = np.minimum(xg[:, 2] / h, 1.0)
    elif name == "dp_weighted_apex_column":
        v[:, 2] = -np.minimum((1.0 - xg[:, 2]) / h, 1.0)
    return contract.band_formula(assertion[1], params, h), v


@pytest.mark.zone_a
@pytest.mark.parametrize("formulation", [lr.CONTROL_ID, lr.JUMP_ID])
@pytest.mark.parametrize("name", EXACT)
@pytest.mark.parametrize("mesh", [0, 1])
def test_exact_upper_field_satisfies_the_upper_program(name, mesh, formulation):
    problem, arrays, geo = problem_of(name, mesh)
    program = lr.assemble(problem, arrays, formulation, geometry=geo)
    value, v, footing = exact_upper_field(name, geo, problem)
    x = upper_vector(program, v, footing)
    assert np.abs(program.zero_residual(x)).max() < 1e-12
    for cones in program.cone_slacks(x):
        assert in_cone(cones, 1e-12)
    assert program.value(x) == pytest.approx(value, rel=1e-12)
    field_value = lr.upper_field_value(problem, arrays, "q", v, footing, geo)[0]
    assert field_value == pytest.approx(value, rel=1e-12)


@pytest.mark.zone_a
@pytest.mark.parametrize("formulation", [lr.CONTROL_ID, lr.JUMP_ID])
@pytest.mark.parametrize("name", EXACT_LOWER)
@pytest.mark.parametrize("mesh", [0, 1])
def test_band_field_attains_its_formula_in_the_upper_program(name, mesh, formulation):
    """The gravity work enters with the right sign and quadrature: the band value is exact."""
    problem, arrays, geo = problem_of(name, mesh)
    program = lr.assemble(problem, arrays, formulation, geometry=geo)
    formula, v = band_upper_field(name, geo, CASES[name], mesh)
    footing = np.zeros((len(problem["footings"]), 6))
    w1 = program.layout["w1"]
    raw = upper_vector(program, v, footing)
    v = v / float(w1 @ raw[: len(w1)])
    x = upper_vector(program, v, footing)
    assert np.abs(program.zero_residual(x)).max() < 1e-12
    for cones in program.cone_slacks(x):
        assert in_cone(cones, 1e-12)
    assert program.value(x) == pytest.approx(formula, rel=1e-12)
    # Zeroing gravity breaks the value: the dead-load cases are sensitive.
    if any(problem["regions"][0]["b0"]):
        weightless = contract.copy.deepcopy(problem)
        weightless["regions"][0]["b0"] = [0.0, 0.0, 0.0]
        weightless["problem_sha256"] = contract.document_hash(weightless, arrays, "problem_sha256")
        other = lr.assemble(weightless, arrays, formulation)
        assert abs(other.value(upper_vector(other, v, footing)) - formula) > 1e-3


@pytest.mark.zone_a
def test_upper_rows_are_live():
    """A footprint off the rigid motion, a moving wall, an isochory break and a normal jump each show."""
    problem, arrays, geo = problem_of("vm_plate_shear", 1)
    value, v, footing = exact_upper_field("vm_plate_shear", geo, problem)
    control = lr.assemble(problem, arrays, lr.CONTROL_ID, geometry=geo)
    jump = lr.assemble(problem, arrays, lr.JUMP_ID, geometry=geo)
    xg = geo.global_node_coordinates()
    top = int(np.nonzero(np.abs(xg[:, 1] - 1.0) < 1e-12)[0][0])
    base = int(np.nonzero(np.abs(xg[:, 1]) < 1e-12)[0][0])
    centre = int(np.argmin(np.linalg.norm(xg - 0.5, axis=1)))
    for node, component in ((top, 1), (base, 0), (centre, 0)):
        bad = v.copy()
        bad[node, component] += 1e-3
        x = upper_vector(control, bad, footing)
        x[control.layout["rz"]:] = upper_vector(control, v, footing)[control.layout["rz"]:]
        assert np.abs(control.zero_residual(x)).max() > 1e-6
    # A normal jump at one interior facet node, in the broken space.
    x = upper_vector(jump, v, footing)
    fi = next(f for f in range(geo.F) if 0.0 < xg[arrays["tets"][geo.e_plus[f], geo.plus_nodes[f, 3]], 1] < 1.0)
    e, a = geo.e_plus[fi], geo.plus_nodes[fi, 3]
    x[30 * e + 3 * a: 30 * e + 3 * a + 3] += 1e-3 * geo.interior_normal[fi]
    rows = jump.zero_residual(x)[24 * geo.E: 24 * geo.E + 18 * geo.F].reshape(geo.F, 6, 3)
    assert np.abs(rows[fi, :, 0]).max() > 1e-4


@pytest.mark.zone_a
def test_jump_rows_are_the_bernstein_coefficients_of_the_contract_jump():
    problem, arrays, geo = problem_of("dp_uniaxial_tension", 1)
    program = lr.assemble(problem, arrays, lr.JUMP_ID, geometry=geo)
    rng = np.random.default_rng(8)
    field = rng.normal(size=(geo.E, 10, 3))
    coeff = lr.jump_coefficients(geo, field)
    for fi, ((e1, f1), (e2, _)) in enumerate(contract.mesh_facets(arrays["tets"])[0]):
        nodal = contract.facet_nodal_jump(arrays["tets"], field, e1, f1, e2)
        np.testing.assert_allclose(coeff[fi], contract.bernstein_from_nodal(nodal), atol=1e-12)
    x = np.zeros(program.A.shape[1])
    x[: 30 * geo.E] = field.reshape(-1)
    rows = (program.A @ x)[24 * geo.E: 24 * geo.E + 18 * geo.F].reshape(geo.F, 6, 3)
    basis = np.concatenate([geo.interior_normal[:, None], lr.tangent_bases(geo.interior_normal)], axis=1)
    np.testing.assert_allclose(rows, np.einsum("fqc,fbc->fbq", basis, coeff), atol=1e-12)
    # The ADR 10.2 pattern: opening at the facet vertices only is admissible at the nodes, not on the facet.
    n = geo.interior_normal[0]
    nodal = np.zeros((6, 3))
    nodal[:3] = 1e-3 * n
    bern = lr.BERNSTEIN @ nodal
    assert np.all(nodal @ n >= 0.0) and np.min(bern @ n) < 0.0


@pytest.mark.zone_a
def test_upper_value_matches_the_contract_fixture_tooling():
    """The work and dissipation formulas agree with the contract's fixture values, jumps included."""
    for stem, formulation in (("vm_dead_gravity_shear_layer.m0", lr.CONTROL_ID),
                              ("vm_ground_patch_halfspace.m0", lr.JUMP_ID),
                              ("dp_free_hydrostatic.m0", lr.JUMP_ID)):
        problem, arrays = contract.load_document(contract.DATA / f"{stem}.problem.json")
        result, fields = contract.load_document(contract.DATA / f"{stem}.q.{formulation}.result.json")
        velocity = fields["element_velocity"] if formulation == lr.JUMP_ID else fields["velocity"]
        mine = lr.upper_field_value(problem, arrays, "q", velocity, fields["footing_velocity"])[0]
        assert mine == pytest.approx(result["value"], rel=1e-12)
        assert mine == pytest.approx(contract.fixture_upper_value(problem, arrays, velocity, fields["footing_velocity"]),
                                     rel=1e-12)


# ---------------------------------------------------------------------------------------------------
# Programs: shape, export, queries
# ---------------------------------------------------------------------------------------------------

@pytest.mark.zone_a
@pytest.mark.parametrize("name", BUILDABLE)
def test_program_shape_matches_the_contract_on_every_case(name):
    problem, arrays, geo = problem_of(name, 0)
    case = CASES[name]
    for formulation in case["formulation"] + case["control"]:
        program = lr.assemble(problem, arrays, formulation, geometry=geo)
        contract.validate_program(program.document("p.npz"), program.arrays(), problem, arrays)
        a = program.A
        assert a.has_canonical_format and not np.any(a.data == 0.0)


@pytest.mark.zone_a
def test_program_and_result_export_round_trip(tmp_path):
    problem, arrays, geo = problem_of("vm_plate_shear", 0)
    path = lr.export_problem(tmp_path, problem, arrays)
    assert path.name == "vm_plate_shear.m0.problem.json"
    loaded, loaded_arrays = contract.load_document(path)
    contract.validate_problem(loaded, loaded_arrays)
    value, v, footing = exact_upper_field("vm_plate_shear", geo, problem)
    for formulation in FORMULATIONS:
        program = lr.assemble(loaded, loaded_arrays, formulation)
        path = lr.export_program(tmp_path, program)
        assert path.name == f"vm_plate_shear.m0.q.{formulation}.program.json"
        doc, parts = contract.load_document(path)
        contract.validate_program(doc, parts, loaded, loaded_arrays)
        a = sp.csc_matrix((parts["A_data"], parts["A_indices"], parts["A_indptr"]), shape=(doc["n_rows"], doc["n_variables"]))
        assert (a != program.A).nnz == 0
        np.testing.assert_array_equal(parts["b"], program.b)
        np.testing.assert_array_equal(parts["c"], program.c)
        if program.bound == "lower":
            x = lower_vector(program, *exact_lower_field("vm_plate_shear", program.layout["geometry"]))
        else:
            x = upper_vector(program, v, footing)
        record = {"name": "none", "version": "0", "settings": {}}
        rdoc, rarrays = lr.export_result(tmp_path, program, "solved", x, record, 0.0)
        reloaded, rfields = contract.load_document(tmp_path / f"vm_plate_shear.m0.q.{formulation}.result.json")
        contract.validate_result(reloaded, rfields, loaded, loaded_arrays)
        assert reloaded["value"] == pytest.approx(1.0, rel=1e-12)
        assert reloaded["result_sha256"] == rdoc["result_sha256"]
    nothing, empty = lr.result_document(program, "time_limit")
    assert nothing["value"] is None and nothing["objective"] is None and not empty
    contract.validate_result(nothing, empty, loaded, loaded_arrays)


PLATE_FREE = """
case plate_free_sides
units       length m; force kN; stress kPa
box         0 0 0 1 1 1
mesh        1 1 1
region      soil von_mises k 1.0
footing     plate xref 0.5 1 0.5 law bonded
boundary    base fixed plane y 0
boundary    top footing plane y 1 footing plate
boundary    x0 free plane x 0
boundary    x1 free plane x 1
boundary    z0 free plane z 0
boundary    z1 free plane z 1
query       q ray
F1          plate 1 0 0 0 0 0
end
"""


def plate_with_queries(queries):
    case = contract.parse_case_text(PLATE_FREE)["plate_free_sides"]
    problem, arrays = contract.build_case_problem(case, 0)
    problem["queries"] = queries
    problem["problem_sha256"] = contract.document_hash(problem, arrays, "problem_sha256")
    contract.validate_problem(problem, arrays)
    return problem, arrays


@pytest.mark.zone_a
def test_support_query_program():
    d = [[1.0, 0.0, 0.0, 0.0, 0.5, 0.0]]
    problem, arrays = plate_with_queries([{"id": "s", "form": "support", "restrained": [[0, 1]], "direction": d}])
    geo = lr.Geometry(problem, arrays)
    lower = lr.assemble(problem, arrays, lr.LOWER_ID, "s", geometry=geo)
    contract.validate_program(lower.document("p.npz"), lower.arrays(), problem, arrays)
    assert lower.layout["lambda"] is None
    rmat = lr.footing_resultant_matrix(geo)
    np.testing.assert_allclose(lower.c, -(np.asarray(d).reshape(-1) @ rmat), atol=1e-15)
    rng = np.random.default_rng(9)
    sigma = rng.normal(size=(geo.E, 4, 6))
    x = lower_vector(lower, None, sigma)
    assert lower.value(x) == pytest.approx(float(np.asarray(d).reshape(-1) @ (rmat @ sigma.reshape(-1))), rel=1e-12)
    for formulation in (lr.CONTROL_ID, lr.JUMP_ID):
        upper = lr.assemble(problem, arrays, formulation, "s", geometry=geo)
        contract.validate_program(upper.document("p.npz"), upper.arrays(), problem, arrays)
        # V_fc = d_fc over unrestrained components, in order; the restrained one is V = 0.
        tail = upper.b[upper.n_zero - 5: upper.n_zero]
        np.testing.assert_array_equal(tail, [1.0, 0.0, 0.0, 0.5, 0.0])
        # A shear field with plate velocity d's translation: U = D - W0 = k * volume, W0 = 0.
        xg = geo.global_node_coordinates()
        v = np.zeros_like(xg)
        v[:, 0] = xg[:, 1]
        footing = np.array([[1.0, 0, 0, 0, 0, 0]])
        x = upper_vector(upper, v, footing)
        residual = upper.zero_residual(x)
        assert np.abs(residual[:-5]).max() < 1e-12      # everything but V = d holds (ry differs)
        assert upper.value(x) == pytest.approx(1.0, rel=1e-12)
        assert lr.upper_field_value(problem, arrays, "s", v, footing, geo)[0] == pytest.approx(1.0, rel=1e-12)


@pytest.mark.zone_a
def test_fixed_components_and_restrained_rows():
    query = {"id": "f", "form": "fixed_components", "restrained": [[0, 5]], "fixed": [[0, 1, -0.5]],
             "direction": [[1.0, 0.0, 0.0, 0.0, 0.0, 0.0]]}
    problem, arrays = plate_with_queries([query])
    geo = lr.Geometry(problem, arrays)
    lower = lr.assemble(problem, arrays, lr.LOWER_ID, "f", geometry=geo)
    contract.validate_program(lower.document("p.npz"), lower.arrays(), problem, arrays)
    # Ray rows for the five unrestrained components: R - lambda F1 = F0, F0 = the fixed values.
    start = lower.n_zero - 5
    np.testing.assert_array_equal(lower.b[start: lower.n_zero], [0.0, -0.5, 0.0, 0.0, 0.0])
    lam_col = lower.A[:, 0].toarray().ravel()[start: lower.n_zero]
    np.testing.assert_array_equal(lam_col, [-1.0, 0.0, 0.0, 0.0, 0.0])
    for formulation in (lr.CONTROL_ID, lr.JUMP_ID):
        upper = lr.assemble(problem, arrays, formulation, "f", geometry=geo)
        contract.validate_program(upper.document("p.npz"), upper.arrays(), problem, arrays)
        layout = upper.layout
        np.testing.assert_array_equal(layout["w0"][layout["V"]:], [0.0, -0.5, 0.0, 0.0, 0.0, 0.0])
        np.testing.assert_array_equal(layout["w1"][layout["V"]:], [1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        restrained_row = upper.n_zero - 2
        row = upper.A.getrow(restrained_row)
        assert list(row.indices) == [layout["V"] + 5] and upper.b[restrained_row] == 0.0
        assert upper.b[upper.n_zero - 1] == 1.0


# ---------------------------------------------------------------------------------------------------
# Open solver (Zone B): estimates only
# ---------------------------------------------------------------------------------------------------

def solved(name, mesh, formulation, interior=False):
    pytest.importorskip("clarabel")
    key = ("solved", name, mesh, formulation, interior)
    if key not in _CACHE:
        problem, arrays, geo = problem_of(name, mesh)
        program = lr.assemble(problem, arrays, formulation, geometry=geo)
        status, x, wall, record = lr.solve(program)
        extra = None
        if interior:
            found = lr.lower_interior_field(program)
            extra = None if found is None else found[:2]
        doc, fields = lr.result_document(program, status, x, record, wall, extra)
        contract.validate_result(doc, fields, problem, arrays)
        _CACHE[key] = (program, x, doc, fields)
    return _CACHE[key]


SOLVER_SLACK = 1e-6   # relative; Clarabel 0.11.1 default tolerances are 1e-8, see the module docstring


@pytest.mark.zone_b
@pytest.mark.parametrize("formulation", FORMULATIONS)
@pytest.mark.parametrize("name", EXACT)
def test_estimates_of_exact_cases_reach_the_target(name, formulation):
    program, x, doc, _ = solved(name, 0, formulation)
    target = CASES[name]["target"]
    assert doc["solver_status"] == "solved"
    value = doc["value"]
    if formulation == lr.LOWER_ID:
        assert target - SOLVER_SLACK * target <= value <= target + 1e-9 * target
    else:
        assert target - 1e-9 * target <= value <= target + SOLVER_SLACK * target
    assert np.abs(program.zero_residual(x)).max() < 1e-7


@pytest.mark.zone_b
@pytest.mark.parametrize("name", EXACT_LOWER)
def test_estimates_of_exact_lower_cases_respect_the_band_formula(name):
    case = CASES[name]
    formula, _ = band_upper_field(name, problem_of(name, 0)[2], case, 0)
    target = case["target"]
    lower = solved(name, 0, lr.LOWER_ID)[2]["value"]
    assert target - SOLVER_SLACK * target <= lower <= target + 1e-9 * target
    for formulation in (lr.CONTROL_ID, lr.JUMP_ID):
        upper = solved(name, 0, formulation)[2]["value"]
        assert target <= upper <= formula + 1e-9 * formula


@pytest.mark.zone_b
def test_estimates_of_the_column_sit_inside_the_sandwich():
    params = {"k": 1.0, "gamma": 0.25, "H": 2.0}
    low, high = (contract.TARGET_FORMULAS[name](params) for name in ("column_static", "column_uniform"))
    lower = solved("vm_weighted_column_sandwich", 0, lr.LOWER_ID)[2]["value"]
    jump = solved("vm_weighted_column_sandwich", 0, lr.JUMP_ID)[2]["value"]
    assert low - 1e-9 <= lower <= jump <= high + 1e-9


@pytest.mark.zone_b
def test_zero_gravity_mutant_is_killed_by_the_dead_shear_layer():
    pytest.importorskip("clarabel")
    problem, arrays, geo = problem_of("vm_dead_gravity_shear_layer", 1)
    case = CASES["vm_dead_gravity_shear_layer"]
    weightless = contract.copy.deepcopy(problem)
    weightless["regions"][0]["b0"] = [0.0, 0.0, 0.0]
    weightless["problem_sha256"] = contract.document_hash(weightless, arrays, "problem_sha256")
    formula, _ = band_upper_field("vm_dead_gravity_shear_layer", geo, case, 1)
    lower = lr.assemble(weightless, arrays, lr.LOWER_ID)
    upper = lr.assemble(weightless, arrays, lr.JUMP_ID)
    low = lower.value(lr.solve(lower)[1])
    high = lr.result_document(upper, "solved", lr.solve(upper)[1])[0]["value"]
    assert low < case["target"] - 0.1 and high > formula + 0.1


@pytest.mark.zone_b
def test_lower_interior_field_is_strictly_admissible():
    program, x, doc, fields = solved("vm_dead_gravity_shear_layer", 1, lr.LOWER_ID, interior=True)
    problem, arrays, geo = problem_of("vm_dead_gravity_shear_layer", 1)
    sigma_int = fields["interior_sigma"]
    margin = lr.yield_margin(problem, arrays, sigma_int, geo)
    assert margin.max() < -1e-3
    xi = lower_vector(program, doc["interior_value"], sigma_int)
    assert np.abs(program.zero_residual(xi)).max() < 1e-7


@pytest.mark.zone_b
def test_strip_estimates_on_mesh_0_straddle_the_target():
    target = CASES["tresca_strip_vm_thin_slab"]["target"]
    lower = solved("tresca_strip_vm_thin_slab", 0, lr.LOWER_ID)[2]["value"]
    jump = solved("tresca_strip_vm_thin_slab", 0, lr.JUMP_ID)[2]["value"]
    control = solved("tresca_strip_vm_thin_slab", 0, lr.CONTROL_ID)[2]["value"]
    assert lower < target < jump < control


@pytest.mark.zone_b
@pytest.mark.slow
@pytest.mark.parametrize("name", ["upper3d_fixed_phi00", "dp_prandtl_cohesion_thin_slab"])
def test_slow_estimates_are_ordered(name):
    lower = solved(name, 0, lr.LOWER_ID)[2]
    upper = solved(name, 0, lr.JUMP_ID)[2]
    assert lower["solver_status"] in ("solved", "almost_solved")
    assert upper["solver_status"] in ("solved", "almost_solved")
    assert lower["value"] < upper["value"]
