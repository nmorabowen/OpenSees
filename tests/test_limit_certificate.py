"""WP-237c: tests of the independent certificate checker, ``tests/limit_certificate.py``.

The checker is written from the frozen contract alone. These tests pin its tables to the contract test, compare its
operators with the contract's convention helpers (which use other routes), run every committed valid fixture and
every negative fixture through the full procedure, and read the checker's certificates back through the contract's
validators and case evaluator. They also hold the review probes: a nearly dependent row block, a von Mises opening
with the normal-jump mutant, final checks that must refuse, a restrained component, capacity refusals, and the method
report without which no certificate reaches a gate. ``tests/limit_conditioning_report.json`` is the committed output
of ``python tests/limit_certificate.py conditioning``.

Uses the standard library, numpy, scipy and pytest.
"""
import copy
import math
from pathlib import Path

import numpy as np
import pytest

import limit_certificate as lc
import test_limit_contract as ct

pytestmark = pytest.mark.zone_a

DATA = ct.DATA
CASES = ct.parse_cases()
T = ct.TOLERANCES
NO_FOOTING = np.zeros((0, 6))


# ---------------------------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------------------------

def stored(name, formulation, mesh=0):
    problem, arrays = ct.load_document(DATA / f"{ct.case_stem(name, mesh)}.problem.json")
    result, fields = ct.load_document(DATA / f"{ct.record_stem(name, formulation, mesh)}.result.json")
    return problem, arrays, result, fields


def make_result(name, problem, arrays, bound, value, fields, formulation, mesh=0, **extra):
    """A result document of an analytic producer, valid by the contract's validator. The hash is the contract's."""
    fields = {key: np.asarray(item, dtype="<f8") for key, item in fields.items()}
    doc = ct._result(name, problem, bound, value, fields, formulation, mesh)
    doc.update(extra)
    doc["result_sha256"] = ct.document_hash(doc, fields, "result_sha256")
    ct.validate_result(doc, fields, problem, arrays)
    return doc, fields


def upper_result(name, problem, arrays, velocity, footing, formulation, mesh=0, **more):
    """An Upper result whose value is that of its field, computed by the contract's fixture tooling."""
    broken = formulation == ct.JUMP_ID
    field = velocity[arrays["tets"]] if broken and velocity.ndim == 2 else velocity
    fields = {"element_velocity" if broken else "velocity": field, "footing_velocity": footing, **more}
    value = ct.fixture_upper_value(problem, arrays, field, footing)
    return make_result(name, problem, arrays, "upper", value, fields, formulation, mesh)


def lower_result(name, problem, arrays, value, sigma, mesh=0, extra_fields=None, **extra):
    fields = {"sigma": sigma, "resultant": ct.fixture_resultant(problem, arrays, sigma), **(extra_fields or {})}
    return make_result(name, problem, arrays, "lower", value, fields, ct.LOWER_ID, mesh, **extra)


def certified(problem, arrays, result, fields, **options):
    """Run the checker and read its certificate back through the contract: valid, real and gate-admissible."""
    outcome = lc.check(problem, arrays, result, fields, **options)
    certificate = outcome.certificate
    assert certificate["reasons"] == [], certificate["reasons"]
    ct.validate_certificate(certificate, problem, arrays, result)
    form = ct.query_by_id(problem, result["query_id"])["form"]
    assert ct.gate_admissible(certificate, form) and lc.gate_admissible(certificate, outcome.method)
    assert not lc.gate_admissible(certificate), "without its method report a certificate satisfies no gate"
    assert outcome.method["query_form"] == form and lc.method_failures(certificate, outcome.method) == []
    assert certificate["checker"]["name"] == lc.CHECKER_NAME and outcome.passed is None
    return outcome


def estimate(problem, arrays, result, fields, code, **options):
    """Run the checker and require an estimate with exactly one reason code."""
    outcome = lc.check(problem, arrays, result, fields, **options)
    return is_estimate(outcome, code)


def is_estimate(outcome, code):
    certificate = outcome.certificate
    assert certificate["classification"] == "estimate" and certificate["reasons"] == [code], certificate["reasons"]
    assert certificate["value_certified"] is None and certificate["floor"] is None
    assert not ct.gate_admissible(certificate) and not lc.gate_admissible(certificate, outcome.method)
    return outcome


def refused_final(problem, arrays, result, fields, code):
    """The final check alone must refuse the stored field with exactly one code."""
    outcome = lc.final_check(problem, arrays, result, fields)
    assert outcome.passed is False
    return is_estimate(outcome, code)


def passes_final(problem, arrays, result, fields):
    """A final check that finds nothing gives a verdict and no certificate of any kind."""
    outcome = lc.final_check(problem, arrays, result, fields)
    assert outcome.passed is True and outcome.certificate is None and outcome.method is None and not outcome.certified
    return outcome


def case_problem(name, mesh=0):
    return ct.build_case_problem(CASES[name], mesh)


def exact_fields(name, arrays):
    """Analytic fields of interface 16 on any mesh of a case: (value, sigma, velocity, footing velocity)."""
    x = arrays["nodes"]
    sigma = np.zeros((len(arrays["tets"]), 4, 6))
    velocity = np.zeros_like(x)
    if name == "vm_uniaxial_tension":
        sigma[..., 0] = ct.S3
        velocity[:, 0], velocity[:, 1], velocity[:, 2] = x[:, 0], -0.5 * x[:, 1], -0.5 * x[:, 2]
        return ct.S3, sigma, velocity, NO_FOOTING
    if name in ("vm_plate_shear", "vm_plate_shear_offset_ref"):
        sigma[..., 3] = 1.0
        velocity[:, 0] = x[:, 1]
        return 1.0, sigma, velocity, np.array([[1.0, 0, 0, 0, 0, 0]])
    if name == "dp_hydrostatic_tension":
        sigma[..., :3] = 1.0 / 0.6
        return 1.0 / 0.6, sigma, x.copy(), NO_FOOTING
    raise KeyError(name)


# ---------------------------------------------------------------------------------------------
# The checker's tables are the frozen ones
# ---------------------------------------------------------------------------------------------

def test_tables_match_the_frozen_contract():
    assert lc.TOLERANCES == ct.TOLERANCES
    assert lc.TAU_EQ == ct.tau_eq() == 64.0 * 2.0 ** -53
    assert set(lc.CHECK_KEYS) == ct.CHECK_KEYS and set(lc.REASON_CODES) == ct.REASON_CODES
    assert lc.CONVENTIONS == ct.CONVENTIONS
    assert (lc.EDGES, lc.FACET_VERTICES, lc.FACET_NODES) == (ct.EDGES, ct.FACET_VERTICES, ct.FACET_NODES)
    assert (lc.LOWER_ID, lc.CONTROL_ID, lc.JUMP_ID) == (ct.LOWER_ID, ct.CONTROL_ID, ct.JUMP_ID)
    assert {lc.SCHEMA_FIXTURE_CHECKER, lc.PROCEDURE_DEMO_CHECKER} == ct.NON_GATE_CHECKERS
    assert lc.CHECKER_NAME not in ct.NON_GATE_CHECKERS


def test_tolerance_table_is_the_one_printed_in_the_interface():
    text = ct.INTERFACE.read_text(encoding="utf-8")
    for name, value in lc.TOLERANCES.items():
        row = [line for line in text.splitlines() if line.startswith(f"| `{name}` |")]
        assert len(row) == 1 and float(row[0].split("|")[2]) == value, name


def test_hash_is_the_frozen_one():
    doc = {"b": [1, None, True], "a": "x"}
    arrays = {"v": np.array([0.5, -0.0], dtype="<f8")}
    assert lc.document_hash(doc, arrays, "h") == "bba0e8c06001808c4f198c1956e7d81423a19daf45369b633cfc14b0c9f12a46"
    for path in sorted(DATA.glob("*.problem.json")) + sorted(DATA.glob("*.result.json")):
        key = "problem_sha256" if path.name.endswith("problem.json") else "result_sha256"
        loaded, arrays = lc.load_document(path)
        assert lc.document_hash(loaded, arrays, key) == loaded[key] == ct.document_hash(loaded, arrays, key), path.name


def test_quadrature_tables_come_out_of_the_monomial_formula():
    np.testing.assert_allclose(lc.volume_weights(), [-1 / 20] * 4 + [1 / 5] * 6, rtol=1e-15)
    np.testing.assert_allclose(lc.facet_work_table(), ct.facet_work_table(1.0), rtol=1e-14)
    np.testing.assert_allclose(lc.traction_moment_weights(), (np.ones((3, 3)) + np.eye(3)) / 12.0, rtol=1e-15)
    nodal = np.random.default_rng(0).normal(size=(6, 3))
    np.testing.assert_allclose(lc.bernstein_map() @ nodal, ct.bernstein_from_nodal(nodal), rtol=1e-15)


# ---------------------------------------------------------------------------------------------
# L0: operators rebuilt from the problem, against the contract's helpers on an unstructured mesh
# ---------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def gate_model():
    problem, arrays = case_problem("upper3d_fixed_phi10", 0)
    return problem, arrays, lc.Model(problem, arrays)


def test_geometry_by_the_facet_route_matches_the_inverse_route(gate_model):
    _, arrays, model = gate_model
    nodes, tets = arrays["nodes"], arrays["tets"]
    for e in range(0, model.n_elements, 7):
        grads, volume = ct.tet_gradients(nodes[tets[e][:4]])
        np.testing.assert_allclose(model.grad[e], grads, rtol=1e-11, atol=1e-13)
        assert model.volume[e] == pytest.approx(volume, rel=1e-13)
        for f in range(4):
            normal, area = ct.facet_geometry(nodes[tets[e][:4]], f)
            np.testing.assert_allclose(model.normal[e, f], normal, atol=1e-14)
            assert model.area[e, f] == pytest.approx(area, rel=1e-13)
        for m in range(4):
            np.testing.assert_allclose(model.dn[e, m], ct.p2_gradients_at_vertex(grads, m), rtol=1e-11, atol=1e-12)
    interior, exterior = ct.mesh_facets(tets)
    assert [tuple(row) for row in model.ext.tolist()] == exterior
    mine = [((int(a), int(b)), (int(c), int(d))) for (a, b), (c, d) in zip(model.int_minus, model.int_plus)]
    assert mine == interior


def test_interior_facet_nodes_are_matched_by_node_id(gate_model):
    _, arrays, model = gate_model
    tets = arrays["tets"]
    e1, e2 = model.int_minus[:, 0], model.int_plus[:, 0]
    np.testing.assert_array_equal(tets[e1[:, None], model.int_local_minus], tets[e2[:, None], model.int_local_plus])


def test_principal_routes_match_the_invariant_routes(gate_model):
    _, arrays, model = gate_model
    rng = np.random.default_rng(1)
    sigma = rng.normal(size=(40, 6))
    for alpha, k in ((0.0, 1.0), (0.2, 0.3)):
        np.testing.assert_allclose(lc.principal_yield(sigma, k, alpha), [ct.yield_value(s, k, alpha) for s in sigma], atol=1e-13)
    field = rng.normal(size=(model.n_elements, 10, 3))
    theta, enorm, dnorm = lc.principal_rates(model, field)
    nodes, tets = arrays["nodes"], arrays["tets"]
    for e in range(0, model.n_elements, 11):
        for m in range(4):
            d = ct.strain_rate_at_vertex(nodes[tets[e]], field[e], m)
            expected = ct.strain_split(d)
            assert theta[e, m] == pytest.approx(expected[0], abs=1e-11)
            assert enorm[e, m] == pytest.approx(expected[1], rel=1e-11)
            assert dnorm[e, m] == pytest.approx(math.sqrt(d[:3] @ d[:3] + 0.5 * d[3:] @ d[3:]), rel=1e-11)


def test_lower_rows_vanish_on_a_linear_equilibrated_stress_field(gate_model):
    problem, arrays, model = gate_model
    query = lc.Query(model, problem["queries"][0])
    system = lc.LowerSystem(model, query)
    # sigma = G z with constant divergence d: sigma_zz = c z, sigma_zx = s y. div = (0, 0, c).
    x = model.x
    sigma = np.zeros((model.n_elements, 4, 6))
    sigma[..., 2] = 0.7 * x[..., 2]
    sigma[..., 5] = -0.3 * x[..., 1]
    residual = system.matrix @ system.pack(0.0, sigma) - system.rhs
    equilibrium = residual[system.kinds == "equilibrium"].reshape(-1, 3)
    np.testing.assert_allclose(equilibrium, np.tile([0.0, 0.0, 0.7], (model.n_elements, 1)), atol=1e-12)
    np.testing.assert_allclose(residual[system.kinds == "continuity"], 0.0, atol=1e-12)


def test_footing_resultant_operator_is_exact_for_linear_traction(gate_model):
    problem, arrays, model = gate_model
    query = lc.Query(model, problem["queries"][0])
    system = lc.LowerSystem(model, query)
    # sigma_zz = a + b x + c y on the footing z = 0, [-0.5, 0.5]^2, normal +z, x_ref at the origin.
    a, b, c = -2.0, 0.6, -0.4
    x = model.x
    sigma = np.zeros((model.n_elements, 4, 6))
    sigma[..., 2] = a + b * x[..., 0] + c * x[..., 1]
    resultant = (system.resultant @ system.pack(0.0, sigma)).reshape(1, 6)
    # F_z = a; M = int r x t = (c/12, -b/12, 0) for the unit square.
    np.testing.assert_allclose(resultant[0], [0, 0, a, c / 12.0, -b / 12.0, 0], atol=1e-13)
    np.testing.assert_allclose(resultant, ct.fixture_resultant(problem, arrays, sigma), atol=1e-13)


def test_rigid_motion_has_no_strain_and_exact_gravity_work():
    case = CASES["gravity_work_rigid_motion"]
    probe = case["probe"]
    for mesh in range(len(case["mesh"])):
        vertices, tet4 = ct.case_tet4(case, mesh)
        nodes, tets = ct.tet10_from_tet4(vertices, tet4)
        problem, arrays = case_problem("vm_uniaxial_tension", 0)
        problem = copy.deepcopy(problem)
        problem["regions"][0]["b1"] = case["regions"][0]["b0"]
        _, exterior = ct.mesh_facets(tets)
        arrays = {"nodes": nodes, "tets": tets, "element_region": np.zeros(len(tets), dtype="<i8"),
                  "facets": np.array(exterior, dtype="<i8"), "facet_set": np.arange(len(exterior), dtype="<i8") % 6,
                  "t0": np.zeros((len(exterior), 3, 3)), "t1": np.zeros((len(exterior), 3, 3))}
        for item in problem["boundary_sets"]:
            item["kind"] = "free"
        model = lc.Model(problem, arrays)
        velocity = np.array(probe["u"]) + np.cross(np.array(probe["w"]), nodes - np.array(probe["about"]))
        field = velocity[tets]
        theta, enorm, dnorm = lc.principal_rates(model, field)
        assert max(np.abs(theta).max(), enorm.max(), dnorm.max()) < 1e-14
        work, _ = lc.work_rate(model, field, NO_FOOTING, model.b1, model.t1, None)
        assert work == pytest.approx(case["target"], rel=T["TOL_VALUE_MATCH"] * 1e-6)


def test_rigid_footing_rows_and_the_power_identity(gate_model):
    """Interface 6: v = u + w x (x - x_ref) meets the footing rows, has no jump, and int t . v dA = R^T V."""
    problem, arrays, model = gate_model
    query = lc.Query(model, problem["queries"][0])
    motion = np.array([[0.3, -0.2, 0.5, 0.11, -0.07, 0.23]])
    for formulation in (ct.CONTROL_ID, ct.JUMP_ID):
        system = lc.UpperSystem(model, query, formulation)
        velocity = motion[0, :3] + np.cross(motion[0, 3:], system.key_x - model.x_ref[0])
        x = system.pack(velocity, motion)
        residual = system.matrix @ x - system.rhs
        np.testing.assert_allclose(residual[system.kinds == "footing"], 0.0, atol=1e-14)
        assert np.sum(system.kinds == "footing") == 3 * len(system.footing_keys[0])
        flipped = system.pack(velocity, motion * np.array([1, 1, 1, -1, -1, -1]))
        assert np.max(np.abs((system.matrix @ flipped - system.rhs)[system.kinds == "footing"])) > 0.05, "the sign of w is seen"
        field = system.element_field(x)
        theta, enorm, dnorm = lc.principal_rates(model, field)
        assert dnorm.max() < 1e-13
        if formulation == ct.JUMP_ID:
            assert np.max(np.abs(lc.jump_coefficients(model, system, field))) < 1e-14
    # A linear traction on the footing: its power on the rigid motion, by the P2 facet table, against R^T V.
    lower = lc.LowerSystem(model, query)
    sigma = np.zeros((model.n_elements, 4, 6))
    xv = model.x
    sigma[..., 2] = -2.0 + 0.6 * xv[..., 0] - 0.4 * xv[..., 1]
    sigma[..., 5] = 0.3 + 0.2 * xv[..., 1]
    sigma[..., 4] = -0.1 * xv[..., 0]
    resultant = (lower.resultant @ lower.pack(0.0, sigma)).reshape(1, 6)
    power = 0.0
    for row in np.nonzero(model.facet_footing == 0)[0]:
        e, f = model.ext[row]
        traction = np.einsum("cv,iv->ic", lc.traction_matrix(model.ext_normal[row]), sigma[e, list(lc.FACET_VERTICES[f])])
        nodal = field[e, list(lc.FACET_NODES[f])]
        power += model.ext_area[row] * np.einsum("ic,ij,jc->", traction, lc.facet_work_table(), nodal)
    assert power == pytest.approx(float(np.sum(resultant * motion)), rel=1e-12)


def test_malformed_documents_are_refused_not_certified():
    problem, arrays, result, fields = stored("vm_plate_shear", ct.LOWER_ID)
    broken = copy.deepcopy(problem)
    del broken["queries"][0]["F1"]
    broken["problem_sha256"] = ct.document_hash(broken, arrays, "problem_sha256")
    moved = dict(result, problem_sha256=broken["problem_sha256"])
    moved["result_sha256"] = ct.document_hash(moved, fields, "result_sha256")
    with pytest.raises(lc.Refusal) as error:
        lc.check(broken, arrays, moved, fields)
    assert error.value.code == "malformed_document"
    short = dict(arrays, facets=arrays["facets"][:-1], facet_set=arrays["facet_set"][:-1], t0=arrays["t0"][:-1], t1=arrays["t1"][:-1])
    with pytest.raises(lc.Refusal) as error:
        lc.check(problem, short, result, fields)
    assert error.value.code == "facets_mismatch", "an exterior facet left out of the list would drop its boundary condition"
    wrong = dict(fields, sigma=fields["sigma"][:, :3])
    estimate(problem, arrays, result, wrong, "hash_mismatch")
    rehashed = dict(result)
    rehashed["result_sha256"] = ct.document_hash(rehashed, wrong, "result_sha256")
    with pytest.raises(lc.Refusal) as error:
        lc.check(problem, arrays, rehashed, wrong)
    assert error.value.code == "array_shape"


def test_tangent_basis_and_traction_operator_follow_the_interface():
    rng = np.random.default_rng(2)
    for normal in list(np.eye(3)) + [-np.eye(3)[0]] + [v / np.linalg.norm(v) for v in rng.normal(size=(5, 3))]:
        for mine, theirs in zip(lc.tangent_pair(normal), ct.tangent_basis(normal)):
            np.testing.assert_allclose(mine, theirs, atol=1e-15)
        np.testing.assert_array_equal(lc.traction_matrix(normal), ct.voigt_n(normal))


# ---------------------------------------------------------------------------------------------
# Equality machinery: literal zeros, row-normwise omega, conditioning
# ---------------------------------------------------------------------------------------------

def test_literal_zeros_run_to_a_fixed_point_and_skip_inhomogeneous_rows():
    import scipy.sparse as sp
    matrix = sp.csr_matrix(np.array([
        [1.0, 0, 0, 0, 0],      # x0 = 0
        [2.0, 3.0, 0, 0, 0],    # then x1 = 0
        [0, 1.0, 1.0, 1.0, 0],  # two entries remain: not literal
        [0, 0, 0, 0, 4.0],      # x4 = 5: inhomogeneous, not a literal zero
    ]))
    rhs = np.array([0.0, 0.0, 0.0, 5.0])
    np.testing.assert_array_equal(lc.literal_zeros(matrix, rhs), [True, True, False, False, False])


def test_omega_is_row_normwise_with_scales():
    import scipy.sparse as sp
    rng = np.random.default_rng(3)
    dense = rng.normal(size=(6, 9))
    rhs, x, scale = rng.normal(size=6), rng.normal(size=9), np.abs(rng.normal(size=9)) + 1.0
    assert lc.backward_error(sp.csr_matrix(dense), rhs, x, scale) == pytest.approx(ct.row_backward_error(dense, rhs, x, scale), rel=1e-13)
    # A homogeneous row of rounding-small entries is measured against its own size, with the scales (interface 17.1):
    # a componentwise form would put this row at 1 for every field.
    tiny = sp.csr_matrix(np.array([[1e-17, -1e-17, 0.0], [0.0, 0.0, 0.0]]))
    assert lc.backward_error(tiny, np.zeros(2), np.array([1.0, 1.0 - 1e-15, 0.0]), np.ones(3)) < lc.TAU_EQ
    # A row with nothing in it and a right side holds for no field; a row with zero denominator must be exact.
    assert lc.backward_error(sp.csr_matrix((1, 3)), np.array([1e-30]), np.zeros(3), np.ones(3)) == 1.0
    single = sp.csr_matrix(np.array([[1.0, 0.0, 0.0]]))
    assert lc.backward_error(single, np.zeros(1), np.array([1e-30, 0.0, 0.0]), np.zeros(3)) == math.inf
    assert lc.backward_error(single, np.zeros(1), np.zeros(3), np.zeros(3)) == 0.0


@pytest.mark.parametrize("name,formulation", [
    ("vm_plate_shear", ct.LOWER_ID), ("vm_plate_shear", ct.CONTROL_ID), ("vm_plate_shear", ct.JUMP_ID),
    ("tresca_strip_vm_thin_slab", ct.CONTROL_ID),
])
def test_sparse_conditioning_proves_a_bound_above_the_dense_value(name, formulation):
    problem, arrays = case_problem(name, 1 if name == "vm_plate_shear" else 0)
    dense = lc.conditioning_report(problem, arrays, "q", formulation)
    sparse = lc.conditioning_report(problem, arrays, "q", formulation, dense_entries_max=0)
    assert dense["method"].startswith("dense") and sparse["method"].startswith("shift-invert")
    assert dense["capacity_refusal"] is None and sparse["capacity_refusal"] is None
    # The value reported is a bound the route proves: never below the true ratio. Its slack is that of the norm bound
    # on sigma_max and of the second shift; the Ritz estimate beside it is the ratio itself.
    assert dense["kappa_eq"] <= sparse["kappa_eq"] <= 4.0 * dense["kappa_eq"]
    assert sparse["kappa_estimate"] == pytest.approx(dense["kappa_eq"], rel=1e-6)
    dependent = dense["reduced_rows"] - dense["rank"]
    assert dependent > 0, "these systems have dependent rows, so the count has to be proved"
    assert sparse["rank"] == dense["rank"] and sparse["dependent_rows_skipped"] == dependent == dense["dependent_rows_skipped"]
    assert f"{dependent} dependent rows proved by dense singular value decompositions" in sparse["proof"]
    assert sparse["ritz_residual"] < 1e-8
    assert 1.0 <= dense["kappa_eq"] <= T["TOL_KAPPA_CAP"]


def near_dependent_block(perturbation):
    """A 60 by 60 block whose last row is a combination of the others plus a perturbation, among easy blocks."""
    import scipy.sparse as sp
    rng = np.random.default_rng(42)
    block = rng.normal(size=(60, 60))
    block[59] = rng.normal(size=59) @ block[:59] / 8.0 + perturbation * rng.normal(size=60)
    easy = [np.eye(60) + 0.1 * rng.normal(size=(60, 60)) for _ in range(5)]
    matrix = sp.block_diag([block] + easy, format="csr")
    n = matrix.shape[1]
    return matrix, np.zeros(matrix.shape[0]), np.ones(n), np.zeros(n, dtype=bool)


def test_near_dependent_block_is_never_counted_as_dependent_rows():
    """The reviewer's probe: a singular value of 1e-8 sigma_max sits below the first shift of the normal matrix,
    where a count alone would call it a dependent row and report the next singular value. The row blocks show it is
    not dependent at TOL_RANK, so the sparse route gives no kappa_eq at all."""
    matrix, rhs, scale, zero = near_dependent_block(1e-7)
    dense = lc.Conditioning(matrix, rhs, scale, zero)
    assert dense.method.startswith("dense") and dense.rank == 360 and dense.skipped == 0
    assert 1e6 < dense.kappa < 1e10, "the block is nearly dependent, not dependent"
    sparse = lc.Conditioning(matrix, rhs, scale, zero, dense_entries_max=0)
    assert sparse.kappa is None and sparse.method == "unavailable" and sparse.rank is None and sparse.skipped is None
    assert "1 eigenvalues of the normal matrix lie below the first shift" in sparse.refusal
    assert "prove 0 dependent rows" in sparse.refusal
    # The same block with the row exactly dependent: one dependent row is proved, and the bound holds.
    matrix, rhs, scale, zero = near_dependent_block(0.0)
    dense = lc.Conditioning(matrix, rhs, scale, zero)
    sparse = lc.Conditioning(matrix, rhs, scale, zero, dense_entries_max=0)
    assert dense.rank == 359 == sparse.rank and sparse.skipped == 1
    assert dense.kappa <= sparse.kappa <= 4.0 * dense.kappa
    assert "1 dependent rows proved" in sparse.proof and "the largest of 60 rows" in sparse.proof
    # When the block is too large to decompose, an explicit orthonormal basis of the null space proves the count.
    basis = lc.Conditioning(matrix, rhs, scale, zero, dense_entries_max=0, block_entries_max=100)
    assert basis.rank == 359 and "orthonormal basis Y of 1 vectors" in basis.proof and dense.kappa <= basis.kappa <= 4.0 * dense.kappa
    # That basis does not rescue the nearly dependent block either: its vector misses TOL_RANK sigma_max by decades.
    matrix, rhs, scale, zero = near_dependent_block(1e-7)
    near = lc.Conditioning(matrix, rhs, scale, zero, dense_entries_max=0, block_entries_max=100)
    assert near.kappa is None and "an explicit null-space basis reaches only" in near.refusal
    # With neither a decomposable block nor room for a basis, nothing is proved.
    monkey = lc.NULL_BASIS_ENTRIES_MAX
    try:
        lc.NULL_BASIS_ENTRIES_MAX = 10.0
        matrix, rhs, scale, zero = near_dependent_block(0.0)
        none = lc.Conditioning(matrix, rhs, scale, zero, dense_entries_max=0, block_entries_max=100)
        assert none.kappa is None and "is too large as well" in none.refusal
    finally:
        lc.NULL_BASIS_ENTRIES_MAX = monkey


def test_kappa_eq_is_reported_on_a_gate_mesh_before_any_result_exists():
    """Interface 17.1: 237c reports kappa_eq on the gate meshes. The first one is run here, by both routes."""
    problem, arrays = case_problem("upper3d_fixed_phi00", 0)
    report = lc.conditioning_report(problem, arrays, "q", ct.JUMP_ID, dense_entries_max=0)
    assert report["elements"] == 243 and report["dependent_rows_skipped"] == 0 and "forced by the shape" in report["proof"]
    assert 10.0 < report["kappa_estimate"] <= report["kappa_eq"] < 1e3
    frictional = lc.conditioning_report(*case_problem("upper3d_fixed_phi30", 0), "q", ct.JUMP_ID)
    assert frictional["method"].startswith("dense") and frictional["kappa_eq"] < 1e3, "only footprint rows remain"
    # Those rows all share the six footing unknowns: the sparse route takes the other normal matrix, whose zero
    # eigenvalues are the six forced by its shape.
    arrow = lc.conditioning_report(*case_problem("upper3d_fixed_phi30", 0), "q", ct.JUMP_ID, dense_entries_max=0)
    assert "A_hat^T A_hat" in arrow["method"] and "the 6 eigenvalues" in arrow["proof"]
    assert frictional["kappa_eq"] <= arrow["kappa_eq"] <= 4.0 * frictional["kappa_eq"]
    assert arrow["rank"] == frictional["rank"] and arrow["dependent_rows_skipped"] == 0


def test_conditioning_above_the_cap_is_an_estimate_and_a_capacity_limit_is_a_refusal():
    problem, arrays, result, fields = stored("vm_plate_shear", ct.CONTROL_ID)
    model = lc.Model(problem, arrays)
    query = lc.Query(model, problem["queries"][0])
    system = lc.UpperSystem(model, query, ct.CONTROL_ID)
    x = system.pack(fields["velocity"], fields["footing_velocity"])
    zero = lc.literal_zeros(system.matrix, system.rhs)

    class Stub:
        kappa = None

    for kappa, expected in ((2.0 * T["TOL_KAPPA_CAP"], "equality_ill_conditioned"), (50.0, None)):
        Stub.kappa = kappa
        run = lc._Run(problem, result, model, "upper", None, {})
        assert run.recheck_equalities(system, x, system.scales(x), zero, Stub) == expected
    Stub.kappa = None
    with pytest.raises(lc.CapacityRefusal):
        lc._Run(problem, result, model, "upper", None, {}).recheck_equalities(system, x, system.scales(x), zero, Stub)
    # A system the checker cannot factorise is refused outright: no estimate, no certificate, no feasibility record.
    for call in (lc.check, lc.final_check, lc.feasibility):
        with pytest.raises(lc.CapacityRefusal) as error:
            call(problem, arrays, result, fields, dense_entries_max=0, sparse_rows_max=10)
        assert error.value.code == "checker_capacity" and "exceed the sparse factorisation limit of 10" in str(error.value)
    report = lc.conditioning_report(problem, arrays, "q", ct.CONTROL_ID, dense_entries_max=0, sparse_rows_max=10)
    assert report["kappa_eq"] is None and report["rank"] is None and "exceed the sparse factorisation limit" in report["capacity_refusal"]


# ---------------------------------------------------------------------------------------------
# Valid fixtures
# ---------------------------------------------------------------------------------------------

VALID = [
    ("vm_uniaxial_tension", ct.LOWER_ID, ct.S3), ("vm_uniaxial_tension", ct.CONTROL_ID, ct.S3),
    ("vm_uniaxial_tension", ct.JUMP_ID, ct.S3), ("vm_plate_shear", ct.LOWER_ID, 1.0), ("vm_plate_shear", ct.CONTROL_ID, 1.0),
    ("vm_plate_shear", ct.JUMP_ID, 1.0), ("dp_uniaxial_tension", ct.CONTROL_ID, 1.0 / (1.0 / ct.S3 + 0.2)),
    ("dp_uniaxial_tension", ct.JUMP_ID, 1.0 / (1.0 / ct.S3 + 0.2)), ("dp_free_hydrostatic", ct.JUMP_ID, 1.0 / 0.6),
    ("vm_ground_patch_halfspace", ct.JUMP_ID, None), ("vm_dead_gravity_shear_layer", ct.CONTROL_ID, None),
]


@pytest.mark.parametrize("name,formulation,target", VALID, ids=lambda item: str(item))
def test_committed_valid_fixture_is_certified(name, formulation, target):
    problem, arrays, result, fields = stored(name, formulation)
    outcome = certified(problem, arrays, result, fields)
    certificate, checks = outcome.certificate, outcome.certificate["checks"]
    value = certificate["value_certified"]
    assert value == pytest.approx(result["value"], rel=1e-11)
    if target is not None:
        assert value == pytest.approx(target, rel=1e-11)
    assert checks["omega_eq"] <= lc.TAU_EQ and checks["projection_move"] <= 1e-12 and checks["repair_mu"] <= 1e-11
    assert certificate["floor"] == max(T["FLOOR_MIN_REL"], checks["kappa_eq"] * checks["omega_eq"]) * abs(value)
    assert certificate["domain_claim"] == ct.supported_domain_claim(problem, arrays, result["bound"])
    assert outcome.report["kappa_method"]
    # The certified field itself passes the final check when it is handed back, and that check certifies nothing.
    again, again_fields = make_result(name, problem, arrays, result["bound"], value, outcome.fields, formulation)
    passes_final(problem, arrays, again, again_fields)


def test_half_space_fixture_carries_the_zero_extension_claim():
    problem, arrays, result, fields = stored("vm_ground_patch_halfspace", ct.JUMP_ID)
    assert lc.check(problem, arrays, result, fields).certificate["domain_claim"] == "halfspace_by_zero_extension"


def test_exact_lower_at_yield_under_dead_load_needs_a_supplied_interior_field():
    """The committed Lower fixture of the dead-gravity layer is exact and sits on the yield surface at the base.

    Step 4 must repair it, zero stress is not admissible under dead loads, and the fixture supplies no
    ``interior_sigma``: it stays an estimate with ``no_interior_field``. With the field sigma_zx = gamma (H - z) at
    multiplier 0 supplied, the same candidate is certified.
    """
    problem, arrays, result, fields = stored("vm_dead_gravity_shear_layer", ct.LOWER_ID)
    estimate(problem, arrays, result, fields, "no_interior_field")
    z = arrays["nodes"][arrays["tets"][:, :4], 2]
    interior = np.zeros_like(fields["sigma"])
    interior[..., 5] = 0.3 * (1.0 - z)
    doc, more = lower_result("vm_dead_gravity_shear_layer", problem, arrays, result["value"], fields["sigma"],
                             extra_fields={"interior_sigma": interior}, interior_value=0.0)
    outcome = certified(problem, arrays, doc, more)
    value = outcome.certificate["value_certified"]
    assert 0.7 * (1.0 - 1e-11) <= value <= 0.7 and 0.0 < outcome.certificate["checks"]["repair_mu"] < 1e-11
    assert outcome.report["interior_field"] == "supplied"
    # An interior field that is not strictly inside, or does not meet the equalities, is refused.
    doc, more = lower_result("vm_dead_gravity_shear_layer", problem, arrays, result["value"], fields["sigma"],
                             extra_fields={"interior_sigma": fields["sigma"]}, interior_value=result["value"])
    estimate(problem, arrays, doc, more, "interior_field_invalid")
    doc, more = lower_result("vm_dead_gravity_shear_layer", problem, arrays, result["value"], fields["sigma"],
                             extra_fields={"interior_sigma": np.zeros_like(interior)}, interior_value=0.0)
    estimate(problem, arrays, doc, more, "interior_field_invalid")


def test_zero_gravity_mutant_of_the_checker_is_killed(monkeypatch):
    problem, arrays, result, fields = stored("vm_dead_gravity_shear_layer", ct.LOWER_ID)
    original = lc.Model._regions

    def weightless(self, arrays_):
        original(self, arrays_)
        self.b0 = np.zeros_like(self.b0)

    monkeypatch.setattr(lc.Model, "_regions", weightless)
    estimate(problem, arrays, result, fields, "projection_move_too_large")


def test_moment_sign_mutant_of_the_checker_is_killed(monkeypatch):
    problem, arrays = case_problem("vm_plate_shear_offset_ref", 0)
    value, sigma, _, _ = exact_fields("vm_plate_shear_offset_ref", arrays)
    doc, fields = lower_result("vm_plate_shear_offset_ref", problem, arrays, value, sigma)
    assert fields["resultant"][0, 5] == pytest.approx(0.3), "Mz = 0.3 Fx with x_ref above the plate"
    certified(problem, arrays, doc, fields)
    skew = lc.skew
    monkeypatch.setattr(lc, "skew", lambda vector: -skew(vector))
    estimate(problem, arrays, doc, fields, "value_mismatch")


# ---------------------------------------------------------------------------------------------
# Negative fixtures: one expected code each
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize("entry", [e for e in ct.NEGATIVE if e["stage"] == "checker"], ids=lambda e: e["id"])
def test_checker_negative_fixture_ends_with_its_one_code(entry):
    problem, arrays, result, fields, _ = ct.apply_negative(entry)
    outcome = estimate(problem, arrays, result, fields, entry["expect"])
    ct.validate_certificate(outcome.certificate, problem, arrays, result)


@pytest.mark.parametrize("entry", [e for e in ct.NEGATIVE if e["stage"] == "final_check"], ids=lambda e: e["id"])
def test_final_check_negative_fixture_is_refused_without_projection(entry):
    problem, arrays, result, fields, _ = ct.apply_negative(entry)
    outcome = refused_final(problem, arrays, result, fields, entry["expect"])
    omega = outcome.certificate["checks"]["omega_eq"]
    assert lc.TAU_EQ < omega < 6.0 * lc.TAU_EQ, "just above the floor: a few rounding units, far below any solver tolerance"
    # The same candidate through the whole procedure is projected back and certified: the defect is a rounding-level one.
    certified(problem, arrays, result, fields)


def test_schema_negative_fixtures_never_become_certificates():
    seen = set()
    for entry in [e for e in ct.NEGATIVE if e["stage"] == "schema"]:
        problem, arrays, result, fields, forged = ct.apply_negative(entry)
        kind = entry["mutation"]["kind"]
        seen.add(kind)
        if kind == "relabel_upper":
            with pytest.raises(lc.Refusal) as error:
                lc.check(problem, arrays, result, fields)
            assert error.value.code == entry["expect"] == "result_wrong_bound_arrays"
        elif kind == "move_mid_node":
            good, good_fields = ct.load_document(DATA / f"{ct.record_stem('vm_plate_shear', ct.LOWER_ID)}.result.json")
            good["problem_sha256"] = problem["problem_sha256"]
            good["result_sha256"] = ct.document_hash(good, good_fields, "result_sha256")
            outcome = estimate(problem, arrays, good, good_fields, "curved_element")
            assert outcome.certificate["checks"]["straightness_max"] > T["TOL_GEOM_STRAIGHT"]
        elif kind == "problem_hash":
            estimate(problem, arrays, result, fields, "hash_mismatch")
        else:
            # A forged claim: the checker derives the claim itself and gives an estimate for the one asked for.
            outcome = estimate(problem, arrays, result, fields, "unsupported_domain_claim", requested_claim=forged["domain_claim"])
            claim = outcome.certificate["domain_claim"]
            assert claim == ct.supported_domain_claim(problem, arrays, "upper") != forged["domain_claim"]
            assert lc.check(problem, arrays, result, fields).certificate["domain_claim"] == claim
    assert seen == {"relabel_upper", "move_mid_node", "problem_hash", "certificate_claim", "drop_truncation_set"}


def test_every_negative_fixture_is_exercised():
    stages = {entry["stage"] for entry in ct.NEGATIVE}
    assert stages == {"schema", "checker", "final_check"} and len(ct.NEGATIVE) == 23


def test_other_identity_failures_have_their_codes():
    problem, arrays, result, fields = stored("vm_plate_shear", ct.CONTROL_ID)
    for other_id in ("upper.p2c-vertex.v1", "p2c-vertex", "lower.p1d-vertex.v0", ""):       # prefixed or not (interface 8)
        late = dict(result, formulation=other_id)
        late["result_sha256"] = ct.document_hash(late, fields, "result_sha256")
        estimate(problem, arrays, late, fields, "formulation_not_admitted")
    other = dict(result, query_id="other")
    other["result_sha256"] = ct.document_hash(other, fields, "result_sha256")
    estimate(problem, arrays, other, fields, "query_mismatch")
    stale = dict(result, wall_time_s=1.0)
    estimate(problem, arrays, stale, fields, "hash_mismatch")
    empty = dict(result, solver_status="time_limit", value=None, objective=None)
    empty["result_sha256"] = ct.document_hash(empty, {}, "result_sha256")
    outcome = estimate(problem, arrays, empty, {}, "no_candidate")
    assert outcome.certificate["value_candidate"] is None
    # alpha >= 1/sqrt(12): the jump formulation does not apply.
    wide, wide_arrays = case_problem("dp_uniaxial_tension", 0)
    wide = copy.deepcopy(wide)
    wide["regions"][0]["alpha"] = 0.3
    wide["problem_sha256"] = ct.document_hash(wide, wide_arrays, "problem_sha256")
    _, _, jump, jump_fields = stored("dp_uniaxial_tension", ct.JUMP_ID)
    jump = dict(jump, problem_sha256=wide["problem_sha256"])
    jump["result_sha256"] = ct.document_hash(jump, jump_fields, "result_sha256")
    estimate(wide, wide_arrays, jump, jump_fields, "formulation_not_applicable")


VM_FREE_PULL = """
case vm_free_pull
units       length m; force kN; stress kPa
box         0 0 0 1 1 1
mesh        1 1 1
region      soil von_mises k 1.0
boundary    x0 traction plane x 0 t1 -1 0 0
boundary    x1 traction plane x 1 t1 1 0 0
boundary    y0 free plane y 0
boundary    y1 free plane y 1
boundary    z0 free plane z 0
boundary    z1 free plane z 1
query       q ray
end
"""


def opened_cube(amount):
    """A von Mises cube with no velocity condition, and its exact field with one side of the plane x = y moved
    rigidly along the plane's normal: strain-free, so the only defect is an opening on whole facets."""
    problem, arrays = ct.build_case_problem(ct.parse_case_text(VM_FREE_PULL)["vm_free_pull"], 0)
    ct.validate_problem(problem, arrays)
    nodes, tets = arrays["nodes"], arrays["tets"]
    exact = np.column_stack([nodes[:, 0], -0.5 * nodes[:, 1], -0.5 * nodes[:, 2]])[tets]
    side = np.nonzero(nodes[tets[:, :4]].mean(axis=1)[:, 1] > nodes[tets[:, :4]].mean(axis=1)[:, 0])[0]
    opened = exact.copy()
    opened[side] += amount * np.array([-1.0, 1.0, 0.0]) / ct.S2
    return problem, arrays, exact, opened


def test_von_mises_jump_opening_is_refused_by_the_normal_jump_rows():
    problem, arrays, exact, opened = opened_cube(1e-3)
    doc, fields = upper_result("vm_free_pull", problem, arrays, exact, NO_FOOTING, ct.JUMP_ID)
    assert certified(problem, arrays, doc, fields).certificate["value_certified"] == pytest.approx(ct.S3, rel=1e-12)
    doc, fields = upper_result("vm_free_pull", problem, arrays, opened, NO_FOOTING, ct.JUMP_ID)
    # Nothing but the normal jump is wrong: every vertex is isochoric and no boundary row exists.
    model = lc.Model(problem, arrays)
    system = lc.UpperSystem(model, lc.Query(model, problem["queries"][0]), ct.JUMP_ID)
    residual = np.abs(system.matrix @ system.pack(opened, NO_FOOTING) - system.rhs)
    assert set(system.kinds.tolist()) == {"isochory", "normal_jump"}
    assert residual[system.kinds == "isochory"].max() < 1e-14 and residual[system.kinds == "normal_jump"].max() == pytest.approx(1e-3)
    assert np.sum(residual[system.kinds == "normal_jump"] > 9e-4) == 12, "all six coefficients of both facets of the plane"
    outcome = estimate(problem, arrays, doc, fields, "projection_move_too_large")
    assert outcome.certificate["checks"]["projection_move"] > 1e-4
    refused_final(problem, arrays, doc, fields, "equality_residual_above_floor")


def test_normal_jump_mutant_of_the_checker_is_killed(monkeypatch):
    """Dropping the normal-jump rows must turn a test red. The committed fixture ``normal_jump_one_node`` does not
    do it, because its one moved node also strains its element; the opening above does."""
    problem, arrays, _, opened = opened_cube(1e-3)
    doc, fields = upper_result("vm_free_pull", problem, arrays, opened, NO_FOOTING, ct.JUMP_ID)
    entry = [e for e in ct.NEGATIVE if e["id"] == "normal_jump_one_node"][0]
    fixture = ct.apply_negative(entry)[:4]
    monkeypatch.setattr(lc.UpperSystem, "_normal_jump_rows", lambda self, rows: None)
    estimate(*fixture, "projection_move_too_large")           # the mutant survives the committed fixture
    mutant = lc.check(problem, arrays, doc, fields)
    assert mutant.certified, "without the rows an open von Mises facet is certified: the opening test is what kills the mutant"
    monkeypatch.undo()
    estimate(problem, arrays, doc, fields, "projection_move_too_large")


def test_final_check_never_certifies():
    """A final check gives a verdict. A field that meets every final condition gets no certificate from it."""
    for name, formulation in (("vm_plate_shear", ct.CONTROL_ID), ("vm_plate_shear", ct.JUMP_ID), ("dp_free_hydrostatic", ct.JUMP_ID),
                              ("vm_ground_patch_halfspace", ct.JUMP_ID)):
        outcome = passes_final(*stored(name, formulation))
        assert outcome.report["kappa_method"]
    problem, arrays, result, fields = stored("vm_plate_shear", ct.LOWER_ID)
    doc, more = lower_result("vm_plate_shear", problem, arrays, 0.5, 0.5 * fields["sigma"])       # strictly inside yield
    passes_final(problem, arrays, doc, more)
    with pytest.raises(TypeError):
        lc.check(problem, arrays, doc, more, final_check=True)      # the old switch is gone: check always runs the whole procedure


def test_final_check_refuses_a_stress_outside_the_yield_surface():
    problem, arrays, result, fields = stored("vm_plate_shear", ct.LOWER_ID)
    for factor in (1.01, 1.0 + 1e-9):       # far outside, and inside what the whole procedure would repair
        doc, more = lower_result("vm_plate_shear", problem, arrays, factor, factor * fields["sigma"])
        outcome = refused_final(problem, arrays, doc, more, "cone_violation_too_large")
        assert outcome.certificate["checks"]["rho_cone"] > 0.0 and outcome.certificate["checks"]["omega_eq"] <= lc.TAU_EQ
    certified(problem, arrays, doc, more)   # the 1e-9 one is repaired by the whole procedure, which the final check never does
    entry = [e for e in ct.NEGATIVE if e["id"] == "stress_scaled_1p01"][0]
    refused_final(*ct.apply_negative(entry)[:4], "cone_violation_too_large")


def test_final_check_refuses_a_velocity_outside_the_flow_and_jump_cones():
    for fixture in ("dp_flow_cone_one_vertex", "dissipation_charge_as_repair", "jump_between_nodes"):
        entry = [e for e in ct.NEGATIVE if e["id"] == fixture][0]
        outcome = refused_final(*ct.apply_negative(entry)[:4], "flow_violation_too_large")
        checks = outcome.certificate["checks"]
        assert max(checks["rho_flow"], checks["rho_jump_cone"] or 0.0) > T["TOL_CONE_MARGIN"]
    # A slip without opening of 1e-9 on a Drucker-Prager facet: inside TOL_FLOW_ACCEPT, outside the final margin.
    problem, arrays, result, fields = stored("dp_free_hydrostatic", ct.JUMP_ID)
    nodes, tets = arrays["nodes"], arrays["tets"]
    side = np.nonzero(nodes[tets[:, :4]].mean(axis=1)[:, 1] > nodes[tets[:, :4]].mean(axis=1)[:, 0])[0]
    slipped = nodes[tets].copy()
    slipped[side] += 1e-9 * np.array([0.0, 0.0, 1.0])
    doc, more = upper_result("dp_free_hydrostatic", problem, arrays, slipped, NO_FOOTING, ct.JUMP_ID)
    outcome = refused_final(problem, arrays, doc, more, "flow_violation_too_large")
    assert outcome.certificate["checks"]["rho_flow"] == 0.0 and T["TOL_CONE_MARGIN"] < outcome.certificate["checks"]["rho_jump_cone"] < 1e-8


def test_nonzero_velocity_on_a_restrained_component_is_refused():
    """v = (y, 0, a y) under a plate moving with (1, 0, a): every row holds except the restrained uz."""
    base, arrays = case_problem("vm_plate_shear", 0)
    problem = variant(base, arrays, dict(base["queries"][0], restrained=[[0, 2], [0, 4]]))
    nodes = arrays["nodes"]
    for formulation in (ct.CONTROL_ID, ct.JUMP_ID):
        clean = np.column_stack([nodes[:, 1], 0.0 * nodes[:, 1], 0.0 * nodes[:, 1]])
        doc, fields = upper_result("vm_plate_shear", problem, arrays, clean, np.array([[1.0, 0, 0, 0, 0, 0]]), formulation)
        outcome = certified(problem, arrays, doc, fields)
        assert outcome.fields["footing_velocity"][0, 2] == 0.0 and outcome.fields["footing_velocity"][0, 4] == 0.0
        lifted = np.column_stack([nodes[:, 1], 0.0 * nodes[:, 1], 1e-3 * nodes[:, 1]])
        footing = np.array([[1.0, 0, 1e-3, 0, 0, 0]])
        doc, fields = upper_result("vm_plate_shear", problem, arrays, lifted, footing, formulation)
        model = lc.Model(problem, arrays)
        system = lc.UpperSystem(model, lc.Query(model, problem["queries"][0]), formulation)
        key = "element_velocity" if formulation == ct.JUMP_ID else "velocity"
        residual = np.abs(system.matrix @ system.pack(fields[key], footing) - system.rhs)
        assert residual[system.kinds != "restrained"].max() < 1e-14, "only the restrained row is violated"
        assert sorted(residual[system.kinds == "restrained"].tolist()) == [0.0, 1e-3]
        outcome = estimate(problem, arrays, doc, fields, "projection_move_too_large")
        assert outcome.certificate["checks"]["projection_move"] == pytest.approx(1e-3)
        refused_final(problem, arrays, doc, fields, "literal_zero_violated")
        # The same field is fine when uz is not restrained.
        free = variant(base, arrays, dict(base["queries"][0], restrained=[[0, 4]]))
        doc, fields = upper_result("vm_plate_shear", free, arrays, lifted, footing, formulation)
        certified(free, arrays, doc, fields)


def test_a_literal_zero_must_be_exactly_zero_in_the_final_field():
    problem, arrays, result, fields = stored("vm_plate_shear", ct.CONTROL_ID)
    base = int(np.nonzero(np.abs(arrays["nodes"][:, 1]) < 1e-12)[0][0])
    velocity = fields["velocity"].copy()
    velocity[base, 2] = 1e-30      # far below every residual threshold, and still not zero
    doc, more = make_result("vm_plate_shear", problem, arrays, "upper", result["value"],
                            {"velocity": velocity, "footing_velocity": fields["footing_velocity"]}, ct.CONTROL_ID)
    refused_final(problem, arrays, doc, more, "literal_zero_violated")
    outcome = certified(problem, arrays, doc, more)
    assert outcome.fields["velocity"][base, 2] == 0.0 and outcome.report["literal_zeros"] > 0


# ---------------------------------------------------------------------------------------------
# Projection and repair on fields that are not exact
# ---------------------------------------------------------------------------------------------

def test_projection_brings_a_solver_like_lower_field_to_the_rounding_floor():
    problem, arrays, result, fields = stored("vm_plate_shear", ct.LOWER_ID)
    rng = np.random.default_rng(7)
    sigma = fields["sigma"] + 1e-9 * rng.normal(size=fields["sigma"].shape)
    doc, more = lower_result("vm_plate_shear", problem, arrays, 1.0, sigma)
    outcome = certified(problem, arrays, doc, more)
    checks = outcome.certificate["checks"]
    assert max(checks["rho_div"], checks["rho_jump"], checks["rho_bc"]) > 1e-10, "the candidate misses the equalities"
    assert 1e-11 < checks["projection_move"] < 1e-7 and checks["omega_eq"] <= lc.TAU_EQ
    assert checks["rho_cone"] <= 0.0 and 0.0 < checks["repair_mu"] < 1e-7
    assert outcome.certificate["value_certified"] <= 1.0, "a Lower bound of an exact case never exceeds the collapse load"
    assert outcome.certificate["value_certified"] == pytest.approx(1.0, rel=1e-7)


@pytest.mark.parametrize("formulation", [ct.CONTROL_ID, ct.JUMP_ID])
def test_projection_brings_a_solver_like_upper_field_to_the_rounding_floor(formulation):
    problem, arrays, result, fields = stored("vm_plate_shear", formulation)
    rng = np.random.default_rng(8)
    key = "element_velocity" if formulation == ct.JUMP_ID else "velocity"
    velocity = fields[key] + 1e-9 * rng.normal(size=fields[key].shape)
    footing = fields["footing_velocity"] + 1e-9 * rng.normal(size=(1, 6))
    doc, more = upper_result("vm_plate_shear", problem, arrays, velocity, footing, formulation)
    outcome = certified(problem, arrays, doc, more)
    checks = outcome.certificate["checks"]
    assert 1e-11 < checks["projection_move"] < 1e-7 and checks["omega_eq"] <= lc.TAU_EQ
    value = outcome.certificate["value_certified"]
    assert value >= 1.0 - outcome.certificate["floor"], "an Upper bound of an exact case never undercuts the collapse load"
    assert value == pytest.approx(1.0, rel=1e-7)


@pytest.mark.parametrize("bound,formulation", [("lower", ct.LOWER_ID), ("upper", ct.CONTROL_ID), ("upper", ct.JUMP_ID)])
def test_sparse_route_certifies_the_same_field_with_the_same_numbers(bound, formulation):
    """The whole procedure with the iterative conditioning and the factorised projection, against the dense one."""
    problem, arrays = case_problem("vm_plate_shear", 1)
    value, sigma, velocity, footing = exact_fields("vm_plate_shear", arrays)
    rng = np.random.default_rng(11)
    if bound == "lower":
        doc, fields = lower_result("vm_plate_shear", problem, arrays, value, sigma + 1e-9 * rng.normal(size=sigma.shape), 1)
    else:
        noisy = velocity[arrays["tets"]] if formulation == ct.JUMP_ID else velocity
        noisy = noisy + 1e-9 * rng.normal(size=noisy.shape)
        doc, fields = upper_result("vm_plate_shear", problem, arrays, noisy, footing, formulation, 1)
    dense = certified(problem, arrays, doc, fields)
    sparse = certified(problem, arrays, doc, fields, dense_entries_max=0)
    assert dense.method["kappa_method"].startswith("dense") and sparse.method["kappa_method"].startswith("shift-invert")
    one, two = dense.certificate, sparse.certificate
    assert two["value_certified"] == pytest.approx(one["value_certified"], rel=1e-12)
    assert one["checks"]["kappa_eq"] <= two["checks"]["kappa_eq"] <= 4.0 * one["checks"]["kappa_eq"]
    assert sparse.method["kappa_estimate"] == pytest.approx(one["checks"]["kappa_eq"], rel=1e-6)
    assert sparse.method["dependent_rows_skipped"] == dense.method["dependent_rows_skipped"] > 0
    assert "dependent rows proved" in sparse.method["kappa_proof"] and sparse.method["ritz_residual"] < 1e-8
    assert two["checks"]["projection_move"] == pytest.approx(one["checks"]["projection_move"], rel=1e-6)
    assert max(one["checks"]["omega_eq"], two["checks"]["omega_eq"]) <= lc.TAU_EQ
    assert sparse.report["rank"] == dense.report["rank"]


def test_drucker_prager_volume_cone_is_repaired_with_a_relative_interior_field():
    problem, arrays, result, fields = stored("dp_uniaxial_tension", ct.CONTROL_ID)
    target = CASES["dp_uniaxial_tension"]["target"]
    mutation = [e for e in ct.NEGATIVE if e["id"] == "dissipation_charge_as_repair"][0]["mutation"]
    velocity = fields["velocity"].copy()
    velocity[mutation["node"], mutation["component"]] += mutation["amount"]
    signed, size = ct.signed_flow_conditions(problem, arrays, velocity)
    assert T["TOL_CONE_MARGIN"] < signed.max() / size < T["TOL_FLOW_ACCEPT"], "outside the cone, within reach of a repair"
    interior = {"interior_velocity": arrays["nodes"].copy(), "interior_footing_velocity": NO_FOOTING}     # v = x: strictly dilating
    doc, more = upper_result("dp_uniaxial_tension", problem, arrays, velocity, NO_FOOTING, ct.CONTROL_ID, **interior)
    outcome = certified(problem, arrays, doc, more)
    checks = outcome.certificate["checks"]
    assert 0.0 < checks["repair_mu"] < 1e-7 and checks["rho_flow"] <= T["TOL_CONE_MARGIN"]
    final, final_size = ct.signed_flow_conditions(problem, arrays, outcome.fields["velocity"])
    assert final.max() <= T["TOL_CONE_MARGIN"] * final_size, "the contract's own evaluation agrees"
    value = outcome.certificate["value_certified"]
    assert target - outcome.certificate["floor"] <= value <= target * (1.0 + 1e-6)
    # Without the interior field, and with a field that is not strictly inside, there is no certificate.
    doc, more = upper_result("dp_uniaxial_tension", problem, arrays, velocity, NO_FOOTING, ct.CONTROL_ID)
    estimate(problem, arrays, doc, more, "no_admissible_correction")
    flat = {"interior_velocity": fields["velocity"], "interior_footing_velocity": NO_FOOTING}             # on the cone boundary
    doc, more = upper_result("dp_uniaxial_tension", problem, arrays, velocity, NO_FOOTING, ct.CONTROL_ID, **flat)
    estimate(problem, arrays, doc, more, "interior_field_invalid")
    still = {"interior_velocity": np.zeros_like(velocity), "interior_footing_velocity": NO_FOOTING}       # D_scale_int = 0
    doc, more = upper_result("dp_uniaxial_tension", problem, arrays, velocity, NO_FOOTING, ct.CONTROL_ID, **still)
    estimate(problem, arrays, doc, more, "interior_field_invalid")


def test_drucker_prager_jump_cone_is_repaired_on_bernstein_coefficients():
    """A slip of 1e-9 with no opening across one plane leaves the jump cone; an opening field restores it."""
    problem, arrays, result, fields = stored("dp_free_hydrostatic", ct.JUMP_ID)
    nodes, tets = arrays["nodes"], arrays["tets"]
    centroids = nodes[tets[:, :4]].mean(axis=1)
    side = np.nonzero(centroids[:, 1] > centroids[:, 0])[0]
    normal = np.array([-1.0, 1.0, 0.0]) / ct.S2
    candidate, opening = nodes[tets].copy(), nodes[tets].copy()
    candidate[side] += 1e-9 * np.array([0.0, 0.0, 1.0])
    opening[side] += 0.1 * normal
    rows = ct.signed_jump_conditions(problem, arrays, candidate)
    assert max(float(coefficient.max()) for _, coefficient in rows) > 1e-10
    interior = {"interior_element_velocity": opening, "interior_footing_velocity": NO_FOOTING}
    doc, more = upper_result("dp_free_hydrostatic", problem, arrays, candidate, NO_FOOTING, ct.JUMP_ID, **interior)
    outcome = certified(problem, arrays, doc, more)
    checks = outcome.certificate["checks"]
    assert 0.0 < checks["repair_mu"] < 1e-7 and checks["rho_jump_cone"] <= T["TOL_CONE_MARGIN"]
    rows = ct.signed_jump_conditions(problem, arrays, outcome.fields["element_velocity"])
    assert max(float(coefficient.max()) for _, coefficient in rows) <= T["TOL_CONE_MARGIN"] * math.sqrt(3.0)
    value = outcome.certificate["value_certified"]
    assert 5.0 / 3.0 - 2.0 * outcome.certificate["floor"] <= value <= (5.0 / 3.0) * (1.0 + 1e-6)
    doc, more = upper_result("dp_free_hydrostatic", problem, arrays, candidate, NO_FOOTING, ct.JUMP_ID)
    estimate(problem, arrays, doc, more, "no_admissible_correction")
    # A repair that would move the field by more than TOL_REPAIR_MAX is refused.
    weak = nodes[tets].copy()
    weak[side] += 2e-6 * normal
    weak_interior = {"interior_element_velocity": weak, "interior_footing_velocity": NO_FOOTING}
    far = nodes[tets].copy()
    far[side] += 1e-7 * np.array([0.0, 0.0, 1.0])
    doc, more = upper_result("dp_free_hydrostatic", problem, arrays, far, NO_FOOTING, ct.JUMP_ID, **weak_interior)
    estimate(problem, arrays, doc, more, "repair_too_large")


def test_lower_repair_larger_than_the_limit_is_refused():
    problem, arrays, result, fields = stored("vm_plate_shear", ct.LOWER_ID)
    factor = 1.0 + 5e-7           # inside TOL_CONE_ACCEPT, but zero stress repairs it only with mu of about 5e-7
    doc, more = lower_result("vm_plate_shear", problem, arrays, factor, fields["sigma"] * factor)
    outcome = certified(problem, arrays, doc, more)
    assert outcome.certificate["checks"]["repair_mu"] == pytest.approx(5e-7, rel=1e-3)
    assert outcome.certificate["value_certified"] <= 1.0
    # Under a dead load zero is not available; an interior field 1e-4 inside cannot repair 5e-7 within TOL_REPAIR_MAX.
    dead, dead_arrays = case_problem("vm_tension_dead_load", 0)
    sigma = np.zeros((len(dead_arrays["tets"]), 4, 6))
    sigma[..., 0] = ct.S3 * factor
    value = 0.6 * ct.S3 + ct.S3 * (factor - 1.0)
    near = np.zeros_like(sigma)
    near[..., 0] = ct.S3 * (1.0 - 1e-4)
    doc, more = lower_result("vm_tension_dead_load", dead, dead_arrays, value, sigma, extra_fields={"interior_sigma": near},
                             interior_value=0.6 * ct.S3 - 1e-4 * ct.S3)
    estimate(dead, dead_arrays, doc, more, "repair_too_large")


# ---------------------------------------------------------------------------------------------
# The producer cannot certify itself, and stored producer claims are ignored
# ---------------------------------------------------------------------------------------------

def test_a_field_produced_by_the_checker_itself_gets_no_gate_admissible_certificate():
    problem, arrays, result, fields = stored("vm_plate_shear", ct.CONTROL_ID)
    own = dict(result, producer={"name": lc.CHECKER_NAME, "version": lc.CHECKER_VERSION})
    own["result_sha256"] = ct.document_hash(own, fields, "result_sha256")
    certificate = lc.check(problem, arrays, own, fields).certificate
    ct.validate_certificate(certificate, problem, arrays, own)
    assert certificate["checker"]["name"] == ct.PROCEDURE_DEMO_CHECKER
    outcome = lc.check(problem, arrays, own, fields)
    assert not ct.gate_admissible(certificate) and not lc.gate_admissible(certificate, outcome.method)
    record, _ = lc.feasibility(problem, arrays, own, fields)
    assert not ct.gate_admissible(record)


def test_stored_strain_rate_and_dissipation_are_never_read():
    problem, arrays, result, fields = stored("vm_plate_shear", ct.CONTROL_ID)
    n_e = len(arrays["tets"])
    lies = dict(fields, strain_rate=np.full((n_e, 4, 6), 7.0), dissipation=np.full((n_e, 4), -3.0))
    doc, more = make_result("vm_plate_shear", problem, arrays, "upper", result["value"], lies, ct.CONTROL_ID)
    assert certified(problem, arrays, doc, more).certificate["value_certified"] == pytest.approx(1.0, rel=1e-12)


def test_reserved_checker_names_stay_out_of_the_checker():
    for path in DATA.glob("*.certificate.json"):
        committed = ct.load_document(path)[0]
        assert committed["checker"]["name"] == ct.SCHEMA_FIXTURE_CHECKER and not lc.gate_admissible(committed)


# ---------------------------------------------------------------------------------------------
# Domain claims
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize("target,kind,weight,free,bound,expected", [
    ("box", "fixed", False, True, "upper", "box_stated"),
    ("box", "free", False, True, "lower", "box_stated"),
    ("halfspace", "fixed", False, True, "upper", "halfspace_by_zero_extension"),
    ("halfspace", "symmetry", False, True, "upper", "unverified"),
    ("halfspace", "fixed", False, True, "lower", "unverified"),
    ("halfspace", "free", False, True, "lower", "halfspace_by_zero_extension"),
    ("halfspace", "free", True, True, "lower", "unverified"),
    ("halfspace", "free", False, False, "lower", "unverified"),
    ("halfspace", "free", True, True, "upper", "unverified"),
])
def test_domain_claim_is_derived_per_bound(target, kind, weight, free, bound, expected):
    problem, arrays = ct.halfspace_problem()
    for item in problem["boundary_sets"]:
        if item["name"] in problem["domain"]["truncation_sets"]:
            item["kind"] = kind
    problem["regions"][0]["b0"] = [0.0, -1.0 if weight else 0.0, 0.0]
    problem["domain"]["exterior_load_free"] = free
    if target == "box":
        problem["domain"].update(target="box", truncation_sets=[], ground=None)
    model = lc.Model(problem, arrays)
    assert lc.supported_domain_claim(model, bound) == expected == ct.supported_domain_claim(problem, arrays, bound)


@pytest.mark.parametrize("change,code", [
    ({"truncation_sets": ["xmin", "nowhere"]}, "domain_truncation"),
    ({"truncation_sets": []}, "domain_truncation"),
    ({"symmetry_sets": ["nowhere"]}, "domain_symmetry"),
    ({"symmetry_sets": ["xmax"]}, "domain_symmetry"),            # a fixed set declared as a symmetry plane
    ({"ground": None}, "domain_ground"),
    ({"ground": {"axis": "w", "coordinate": 0.0}}, "domain_ground"),
    ({"target": "strip"}, "domain_target"),
    ({"exterior_load_free": "yes"}, "domain_target"),
    ({"target": "box"}, "domain_truncation"),                    # a box with truncation sets
])
def test_invalid_domain_declarations_are_refused_by_name(change, code):
    problem, arrays = ct.halfspace_problem()
    problem["domain"].update(change)
    problem["problem_sha256"] = ct.document_hash(problem, arrays, "problem_sha256")
    with pytest.raises(lc.Refusal) as error:
        lc.Model(problem, arrays)
    assert error.value.code == code
    _, _, result, fields = stored("vm_ground_patch_halfspace", ct.JUMP_ID)
    moved = dict(result, problem_sha256=problem["problem_sha256"])
    moved["result_sha256"] = ct.document_hash(moved, fields, "result_sha256")
    for call in (lc.check, lc.final_check, lc.feasibility):
        with pytest.raises(lc.Refusal) as error:
            call(problem, arrays, moved, fields)
        assert error.value.code == code


def test_half_space_claim_needs_coverage_of_every_face_off_the_ground():
    problem, arrays = ct.halfspace_problem()
    problem["domain"]["truncation_sets"].remove("xmin")
    assert lc.supported_domain_claim(lc.Model(problem, arrays), "upper") == "unverified"
    problem["domain"]["truncation_sets"].append("xmin")
    problem["domain"]["ground"]["coordinate"] = -1.0          # the loaded face is now off the ground plane
    assert lc.supported_domain_claim(lc.Model(problem, arrays), "upper") == "unverified"
    plate, plate_arrays = case_problem("vm_plate_shear", 0)
    plate = copy.deepcopy(plate)
    plate["domain"].update(target="halfspace", truncation_sets=["base"], ground={"axis": "y", "coordinate": 1.0})
    assert lc.supported_domain_claim(lc.Model(plate, plate_arrays), "upper") == "unverified"


# ---------------------------------------------------------------------------------------------
# Brackets, files and the contract's case evaluator
# ---------------------------------------------------------------------------------------------

def case_records(name):
    """Results, certificates and feasibility records of the exact fields on every listed mesh of a case."""
    case = CASES[name]
    results, certificates, feasibility = [], [], []
    for mesh in range(len(case["mesh"])):
        problem, arrays = case_problem(name, mesh)
        value, sigma, velocity, footing = exact_fields(name, arrays)
        pairs = [lower_result(name, problem, arrays, value, sigma, mesh)]
        pairs += [upper_result(name, problem, arrays, velocity, footing, form, mesh) for form in (ct.CONTROL_ID, ct.JUMP_ID)]
        for doc, fields in pairs:
            results.append(doc)
            certificates.append(lc.check(problem, arrays, doc, fields).certificate)
            if doc["bound"] == "upper":
                feasibility.append(lc.feasibility(problem, arrays, doc, fields)[0])
    return case, results, certificates, feasibility


@pytest.mark.parametrize("name", ["vm_uniaxial_tension", "vm_plate_shear", "vm_plate_shear_offset_ref", "dp_hydrostatic_tension"])
def test_exact_case_passes_the_contract_evaluator_on_the_checker_certificates(name):
    case, results, certificates, feasibility = case_records(name)
    assert all(item["classification"] == "certified_bound" for item in certificates)
    verdict = ct.evaluate_case(case, certificates, results, feasibility)
    assert verdict["passed"] and not verdict["invalid_records"] and not verdict["stop_program"], ct.failed(verdict)
    assert {item["kind"] for item in verdict["items"]} == {"lower_le_target", "upper_ge_target", "gap_le"}
    # An estimate in place of one certificate fails the case.
    spoiled = [dict(item) for item in certificates]
    spoiled[0].update(classification="estimate", reasons=["no_interior_field"], value_certified=None, floor=None)
    assert not ct.evaluate_case(case, spoiled, results, feasibility)["passed"]


def test_bracket_record_and_files_are_in_the_frozen_format(tmp_path):
    name = "vm_plate_shear"
    problem, arrays, low, low_fields = stored(name, ct.LOWER_ID)
    lower = lc.check(problem, arrays, low, low_fields)
    for formulation in (ct.CONTROL_ID, ct.JUMP_ID):
        _, _, high, high_fields = stored(name, formulation)
        upper = lc.check(problem, arrays, high, high_fields)
        record = lc.bracket(lower.certificate, upper.certificate, lower.method, upper.method)
        ct.validate_bracket(record, lower.certificate, upper.certificate)
        assert record["status"] == "bracket" and 0.0 <= record["gap"] <= 1e-11 and record["domain_claim"] == "box_stated"
        paths = [lc.write_certificate(tmp_path, low, lower), lc.write_certificate(tmp_path, high, upper),
                 lc.write_bracket(tmp_path, low, high, record)]
        assert [path.name for path in paths] == [
            f"{ct.record_stem(name, ct.LOWER_ID)}.certificate.json", f"{ct.record_stem(name, formulation)}.certificate.json",
            f"{ct.bracket_stem(name, ct.LOWER_ID, formulation)}.bracket.json"]
        for path, document in zip(paths, (lower.certificate, upper.certificate, record)):
            assert ct.load_document(path)[0] == document, "strict JSON, no non-finite number, nothing lost"
        feasibility, _ = lc.feasibility(problem, arrays, high, high_fields)
        path = lc.write_feasibility(tmp_path, high, feasibility)
        assert path.name == f"{ct.record_stem(name, formulation)}.feasibility.json"
        ct.validate_feasibility(ct.load_document(path)[0], problem)


def test_method_report_is_kept_beside_the_certificate_and_gates_need_it(tmp_path):
    """The frozen certificate cannot say how kappa_eq was obtained. The checker's method report can, and a gate
    reader takes no certificate without it."""
    problem, arrays, result, fields = stored("vm_plate_shear", ct.LOWER_ID)
    outcome = certified(problem, arrays, result, fields)
    path = lc.write_certificate(tmp_path, result, outcome)
    beside = path.with_name(f"{ct.record_stem('vm_plate_shear', ct.LOWER_ID)}.certificate.method.json")
    assert beside.exists()
    certificate, method = lc.read_gate_certificate(path)
    assert certificate == outcome.certificate and method == outcome.method == ct.load_document(beside)[0]
    assert method["format"] == lc.FORMAT_METHOD and method["certificate_sha256"] == lc.document_hash(certificate, {}, "")
    assert method["kappa_method"] == "dense singular value decomposition" and method["kappa_eq"] == certificate["checks"]["kappa_eq"]
    assert method["dependent_rows_skipped"] == method["reduced_rows"] - method["rank"] >= 0
    assert method["ritz_residual"] is None and method["interior_field"] == "zero", "this exact field is repaired with zero stress"
    assert method["kappa_proof"] and method["query_form"] == "ray" and method["checker"] == certificate["checker"]

    # Every way of losing or bending the report closes the gate.
    for key, value, failure in (("kappa_method", "unavailable", "kappa_method"), ("kappa_method", None, "kappa_method"),
                                ("kappa_proof", "", "kappa_proof"), ("kappa_eq", 1.0, "kappa_eq"),
                                ("dependent_rows_skipped", method["dependent_rows_skipped"] + 1, "dependent_rows_skipped"),
                                ("dependent_rows_skipped", None, "dependent_rows_skipped"), ("ritz_residual", "small", "ritz_residual"),
                                ("interior_field", "none", "interior_field"), ("interior_field", "magic", "interior_field"),
                                ("query_form", "other", "query_form"), ("certificate_sha256", "0" * 64, "certificate_sha256"),
                                ("checker", {"name": "someone", "version": "1"}, "checker"), ("format", "other", "certificate_method_missing")):
        bent = dict(method, **{key: value})
        assert failure in lc.method_failures(certificate, bent), key
        assert not lc.gate_admissible(certificate, bent), key
    dropped = {key: value for key, value in method.items() if key != "ritz_residual"}
    assert not lc.gate_admissible(certificate, dropped)
    iterative = dict(method, kappa_method="shift-invert Lanczos on A_hat A_hat^T")
    assert lc.method_failures(certificate, iterative) == ["ritz_residual"], "an iterative kappa_eq reports its Ritz residual"
    # A report cannot be moved to another certificate, and a certificate that changed no longer matches its report.
    other = certified(*stored("vm_plate_shear", ct.CONTROL_ID))
    assert not lc.gate_admissible(other.certificate, method) and not lc.gate_admissible(certificate, other.method)
    assert not lc.gate_admissible(dict(certificate, value_certified=certificate["value_certified"] * 2.0), method)

    beside.unlink()
    assert lc.read_certificate(path) == (certificate, None)
    with pytest.raises(lc.Refusal) as error:
        lc.read_gate_certificate(path)
    assert error.value.code == "certificate_method_missing"
    lc.write_document(beside, dict(method, kappa_eq=1.0))
    with pytest.raises(lc.Refusal) as error:
        lc.read_gate_certificate(path)
    assert error.value.code == "certificate_method_invalid"

    # An estimate has a report too, and stays out of every gate. A final check writes nothing.
    low = estimate(*stored("vm_dead_gravity_shear_layer", ct.LOWER_ID), "no_interior_field")
    estimate_path = lc.write_certificate(tmp_path, stored("vm_dead_gravity_shear_layer", ct.LOWER_ID)[2], low)
    with pytest.raises(lc.Refusal) as error:
        lc.read_gate_certificate(estimate_path)
    assert error.value.code in ("certificate_not_gate_admissible", "certificate_method_invalid")
    with pytest.raises(lc.Refusal) as error:
        lc.write_certificate(tmp_path, result, lc.final_check(problem, arrays, result, fields))
    assert error.value.code == "no_certificate"


def test_interior_field_of_a_repair_is_recorded():
    problem, arrays, result, fields = stored("dp_uniaxial_tension", ct.CONTROL_ID)
    assert certified(problem, arrays, result, fields).method["interior_field"] == "none", "no repair, no interior field"
    mutation = [e for e in ct.NEGATIVE if e["id"] == "dissipation_charge_as_repair"][0]["mutation"]
    velocity = fields["velocity"].copy()
    velocity[mutation["node"], mutation["component"]] += mutation["amount"]
    interior = {"interior_velocity": arrays["nodes"].copy(), "interior_footing_velocity": NO_FOOTING}
    doc, more = upper_result("dp_uniaxial_tension", problem, arrays, velocity, NO_FOOTING, ct.CONTROL_ID, **interior)
    outcome = certified(problem, arrays, doc, more)
    assert outcome.method["interior_field"] == "supplied" and outcome.certificate["checks"]["repair_mu"] > 0.0


def rebound(certificate, method, **changes):
    """A certificate with changed numbers and a method report bound to it again: test tooling for bracket statuses."""
    changed = dict(certificate, **changes)
    checks = changed["checks"]
    changed["floor"] = max(T["FLOOR_MIN_REL"], checks["kappa_eq"] * checks["omega_eq"]) * abs(changed["value_certified"])
    return changed, dict(method, certificate_sha256=lc.document_hash(changed, {}, ""), problem_sha256=changed["problem_sha256"],
                         result_sha256=changed["result_sha256"])


def test_bracket_statuses_and_refusals():
    problem, arrays, low, low_fields = stored("vm_plate_shear", ct.LOWER_ID)
    _, _, high, high_fields = stored("vm_plate_shear", ct.CONTROL_ID)
    one, two = lc.check(problem, arrays, low, low_fields), lc.check(problem, arrays, high, high_fields)
    lower, upper, lower_method, upper_method = one.certificate, two.certificate, one.method, two.method
    within, within_method = rebound(lower, lower_method, value_certified=upper["value_certified"] + 0.5 * upper["floor"])
    record = lc.bracket(within, upper, within_method, upper_method)
    ct.validate_bracket(record, within, upper)
    assert record["status"] == "order_within_floor" and record["gap"] is None
    above, above_method = rebound(lower, lower_method, value_certified=1.01)
    record = lc.bracket(above, upper, above_method, upper_method)
    ct.validate_bracket(record, above, upper)
    assert record["status"] == "contradiction" and record["gap"] is None, "a certified Lower above a certified Upper"

    def refused(code, low_pair, high_pair):
        with pytest.raises(lc.Refusal) as error:
            lc.bracket(low_pair[0], high_pair[0], low_pair[1], high_pair[1])
        assert error.value.code == code

    good_low, good_high = (lower, lower_method), (upper, upper_method)
    refused("bracket_labels", good_high, good_low)
    refused("bracket_uncertified", (dict(lower, classification="estimate"), lower_method), good_high)
    refused("hash_mismatch", rebound(lower, lower_method, problem_sha256="0" * 64), good_high)
    refused("query_mismatch", rebound(lower, lower_method, query_id="other"), good_high)
    # Only gate-admissible inputs form a bracket: no missing report, no reserved checker name, no bent certificate.
    refused("bracket_not_gate_admissible", (lower, None), good_high)
    refused("bracket_not_gate_admissible", good_low, (upper, None))
    refused("bracket_not_gate_admissible", good_low, (upper, lower_method))
    refused("bracket_not_gate_admissible", (dict(lower, value_certified=0.5), lower_method), good_high)
    own = dict(low, producer={"name": lc.CHECKER_NAME, "version": lc.CHECKER_VERSION})
    own["result_sha256"] = ct.document_hash(own, low_fields, "result_sha256")
    demo = lc.check(problem, arrays, own, low_fields)
    assert demo.certified and demo.certificate["checker"]["name"] == ct.PROCEDURE_DEMO_CHECKER
    refused("bracket_not_gate_admissible", (demo.certificate, demo.method), good_high)
    for committed in ("vm_plate_shear", "vm_uniaxial_tension"):
        fixture = ct.load_document(DATA / f"{ct.record_stem(committed, ct.LOWER_ID)}.certificate.json")[0]
        refused("bracket_not_gate_admissible", (fixture, lower_method), good_high)

    zero = "halfspace_by_zero_extension"
    assert lc.bracket(*[item for pair in zip(rebound(lower, lower_method, domain_claim=zero), rebound(upper, upper_method, domain_claim=zero))
                        for item in pair])["domain_claim"] == zero
    wide, wide_method = rebound(upper, upper_method, domain_claim=zero)
    assert lc.bracket(lower, wide, lower_method, wide_method)["domain_claim"] == "box_stated"


# ---------------------------------------------------------------------------------------------
# Kinematic feasibility
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize("name,formulation", [
    ("vm_plate_shear", ct.CONTROL_ID), ("vm_plate_shear", ct.JUMP_ID), ("vm_uniaxial_tension", ct.JUMP_ID),
    ("vm_ground_patch_halfspace", ct.JUMP_ID), ("dp_uniaxial_tension", ct.CONTROL_ID), ("dp_free_hydrostatic", ct.JUMP_ID),
])
def test_certified_upper_is_a_feasibility_witness(name, formulation):
    problem, arrays, result, fields = stored(name, formulation)
    record, report = lc.feasibility(problem, arrays, result, fields)
    ct.validate_feasibility(record, problem)
    assert record["status"] == "formulation_feasible" and ct.gate_admissible(record) and lc.gate_admissible(record)
    assert record["witness_result_sha256"] == result["result_sha256"] and record["witness_residual"] <= T["TOL_FEAS_OK"]
    assert record["rho_feas"] <= T["TOL_FEAS_OK"] and record["denominator"] > 0.0
    assert record["nullity"] == record["n_unknowns"] - record["rank"]
    if problem["footings"]:
        assert record["footing_modes"] >= 1
    else:
        assert record["footing_modes"] is None, "no footing: the count does not apply"
    if any(region.get("alpha", 0.0) > 0.0 for region in problem["regions"]):
        assert record["rho_feas"] == record["witness_residual"], "frictional: the witness alone"
    assert report["witness_reasons"] == []


LOCKED_PULL = """
case vm_locked_pull
units       length m; force kN; stress kPa
box         0 0 0 1 1 1
mesh        1 1 1
mesh        2 2 2
region      soil von_mises k 1.0
boundary    x0 fixed plane x 0
boundary    x1 fixed plane x 1
boundary    z0 fixed plane z 0
boundary    z1 fixed plane z 1
boundary    y0 fixed plane y 0
boundary    y1 traction plane y 1 t1 0 1 0
query       q ray
end
"""


@pytest.mark.parametrize("mesh", [0, 1])
@pytest.mark.parametrize("formulation", [ct.CONTROL_ID, ct.JUMP_ID])
def test_locked_space_is_reported_as_a_formulation_result(mesh, formulation):
    """Five fixed faces and a normal pull on the sixth: the divergence theorem gives W1 = 0 for every isochoric
    field with no normal jump, so the space has no mechanism for this load, whatever its nullity."""
    case = ct.parse_case_text(LOCKED_PULL)["vm_locked_pull"]
    problem, arrays = ct.build_case_problem(case, mesh)
    nodes = arrays["nodes"]
    velocity = np.zeros_like(nodes)
    velocity[:, 1] = np.where(np.abs(nodes[:, 1] - 1.0) < 1e-12, 1.0, 0.0)       # does work, is not isochoric
    doc, fields = upper_result("vm_locked_pull", problem, arrays, velocity, NO_FOOTING, formulation, mesh)
    record, report = lc.feasibility(problem, arrays, doc, fields)
    ct.validate_feasibility(record, problem)
    assert record["status"] == "formulation_infeasible" and record["rho_feas"] >= 0.1 > T["TOL_FEAS_LOCKED"]
    assert record["witness_result_sha256"] is None and record["witness_residual"] is None
    assert report["witness_reasons"] == ["projection_move_too_large"] and not ct.gate_admissible(record)
    if formulation == ct.JUMP_ID and mesh == 1:
        assert record["nullity"] > 100, "a large kernel and still locked: nullity does not decide"


def test_locked_verdict_shows_its_least_squares_status_and_bounds():
    case = ct.parse_case_text(LOCKED_PULL)["vm_locked_pull"]
    problem, arrays = ct.build_case_problem(case, 1)
    nodes = arrays["nodes"]
    velocity = np.zeros_like(nodes)
    velocity[:, 1] = np.where(np.abs(nodes[:, 1] - 1.0) < 1e-12, 1.0, 0.0)
    doc, fields = upper_result("vm_locked_pull", problem, arrays, velocity, NO_FOOTING, ct.JUMP_ID, 1)
    record, report = lc.feasibility(problem, arrays, doc, fields)
    for part in ("numerator", "denominator"):
        item = report["least_squares"][part]
        assert all(run["istop"] in (0, 1, 2, 4, 5) for run in item["lsmr"]) and len(item["lsmr"]) == 2
        assert item["lsmr"][0]["relative_residual"] <= lc.LSMR_COMPATIBLE, "the reference work lies in the row space: no mechanism"
        assert 0.0 < item["lower"] <= item["upper"] and item["lower"] == pytest.approx(item["upper"], rel=1e-6)
    assert "stop statuses" in record["method"] and "relative residual" in record["method"]
    assert record["numerator"] == report["least_squares"]["numerator"]["upper"]
    assert record["denominator"] == report["least_squares"]["denominator"]["lower"]


def test_least_squares_that_did_not_converge_gives_no_feasibility_record(monkeypatch):
    case = ct.parse_case_text(LOCKED_PULL)["vm_locked_pull"]
    problem, arrays = ct.build_case_problem(case, 0)
    velocity = np.zeros_like(arrays["nodes"])
    velocity[:, 1] = np.where(np.abs(arrays["nodes"][:, 1] - 1.0) < 1e-12, 1.0, 0.0)
    doc, fields = upper_result("vm_locked_pull", problem, arrays, velocity, NO_FOOTING, ct.JUMP_ID)
    real = lc.spla.lsmr

    def stalled(*args, **kwargs):
        out = list(real(*args, **kwargs))
        out[1] = 7                                       # the iteration limit
        return tuple(out)

    monkeypatch.setattr(lc.spla, "lsmr", stalled)
    with pytest.raises(lc.CapacityRefusal) as error:
        lc.feasibility(problem, arrays, doc, fields)
    assert "LSMR stopped with status 7" in str(error.value)
    monkeypatch.undo()

    # Explicit points alone bound the minimum from above only: without a bound from below there is no locked verdict.
    seminorm = lc._least_seminorm

    def unbounded(*args, **kwargs):
        return dict(seminorm(*args, **kwargs), lower=0.0)

    monkeypatch.setattr(lc, "_least_seminorm", unbounded)
    with pytest.raises(lc.CapacityRefusal):
        lc.feasibility(problem, arrays, doc, fields)


def test_feasibility_rank_leaves_out_the_query_normalisation_rows():
    case = copy.deepcopy(CASES["vm_plate_shear"])
    for item in case["boundaries"]:
        if item["kind"] == "traction":
            item["kind"] = "free"
            item.pop("t1")
    base, arrays = ct.build_case_problem(case, 0)
    problem = variant(base, arrays, {"id": "q", "form": "support", "restrained": [[0, 2]], "direction": [[1.0, 0, 0, 0, 0, 0]]})
    velocity = np.zeros_like(arrays["nodes"])
    velocity[:, 0] = arrays["nodes"][:, 1]
    made = {"velocity": velocity, "footing_velocity": np.array([[1.0, 0, 0, 0, 0, 0]])}
    doc, fields = make_result("vm_plate_shear", problem, arrays, "upper", 1.0, made, ct.CONTROL_ID)
    record, report = lc.feasibility(problem, arrays, doc, fields)
    ct.validate_feasibility(record, problem)
    model = lc.Model(problem, arrays)
    system = lc.UpperSystem(model, lc.Query(model, problem["queries"][0]), ct.CONTROL_ID)
    assert np.sum(system.kinds == "support") == 5 and np.sum(system.kinds == "restrained") == 1
    assert record["n_rows"] == system.matrix.shape[0] - 5, "V = d is the normalisation, not a row of E"
    ray, _ = lc.feasibility(*stored("vm_plate_shear", ct.CONTROL_ID))
    assert record["rank"] == ray["rank"] + 1 and record["footing_modes"] == ray["footing_modes"] - 1, "one restrained component more"
    assert record["status"] == "formulation_feasible" and record["rho_feas"] <= T["TOL_FEAS_OK"]


def test_feasible_space_without_a_valid_witness_is_inconclusive():
    problem, arrays, result, fields = stored("vm_plate_shear", ct.JUMP_ID)
    backwards = {"element_velocity": -fields["element_velocity"], "footing_velocity": -fields["footing_velocity"]}
    doc, more = make_result("vm_plate_shear", problem, arrays, "upper", -result["value"], backwards, ct.JUMP_ID)
    record, report = lc.feasibility(problem, arrays, doc, more)
    ct.validate_feasibility(record, problem)
    assert report["witness_reasons"] == ["work_nonpositive"]
    assert record["status"] == "inconclusive" and record["rho_feas"] <= T["TOL_FEAS_OK"] and record["witness_residual"] is None
    assert not ct.gate_admissible(record)


def test_frictional_witness_outside_the_cone_is_inconclusive_never_infeasible():
    entry = [e for e in ct.NEGATIVE if e["id"] == "dp_flow_cone_one_vertex"][0]
    problem, arrays, result, fields, _ = ct.apply_negative(entry)
    record, report = lc.feasibility(problem, arrays, result, fields)
    ct.validate_feasibility(record, problem)
    assert record["status"] == "inconclusive" and record["witness_residual"] is None
    assert record["rho_feas"] > T["TOL_FEAS_OK"] and report["witness_reasons"] == ["flow_violation_too_large"]


# ---------------------------------------------------------------------------------------------
# Other query forms
# ---------------------------------------------------------------------------------------------

def variant(problem, arrays, query):
    doc = copy.deepcopy(problem)
    doc["queries"] = [query]
    doc["problem_sha256"] = ct.document_hash(doc, arrays, "problem_sha256")
    ct.validate_problem(doc, arrays)
    return doc


def test_fixed_components_query_is_a_ray_with_its_fixed_values():
    base, arrays = case_problem("vm_plate_shear", 0)
    direction = base["queries"][0]["F1"]
    problem = variant(base, arrays, {"id": "q", "form": "fixed_components", "restrained": [], "fixed": [[0, 1, 0.0], [0, 2, 0.0]],
                                     "direction": direction})
    value, sigma, velocity, footing = exact_fields("vm_plate_shear", arrays)
    doc, fields = lower_result("vm_plate_shear", problem, arrays, value, sigma)
    low = certified(problem, arrays, doc, fields)
    lower = low.certificate
    # The contract's fixture tooling reads F0 and F1 of a ray; the value of this field is k = 1 (interface 16.2).
    made = {"element_velocity": velocity[arrays["tets"]], "footing_velocity": footing}
    doc, fields = make_result("vm_plate_shear", problem, arrays, "upper", 1.0, made, ct.JUMP_ID)
    high = certified(problem, arrays, doc, fields)
    upper = high.certificate
    assert lower["value_certified"] == pytest.approx(1.0, rel=1e-11) and upper["value_certified"] == pytest.approx(1.0, rel=1e-11)
    assert upper["checks"]["work_reference"] == pytest.approx(1.0)
    ct.validate_bracket(lc.bracket(lower, upper, low.method, high.method), lower, upper)


def test_support_query_gives_a_point_and_a_half_space():
    """Plate on a block with free sides and no reference load: d^T R with d the plate's x translation."""
    case = copy.deepcopy(CASES["vm_plate_shear"])
    for item in case["boundaries"]:
        if item["kind"] == "traction":
            item["kind"] = "free"
            item.pop("t1")
    base, arrays = ct.build_case_problem(case, 0)
    problem = variant(base, arrays, {"id": "q", "form": "support", "restrained": [], "direction": [[1.0, 0, 0, 0, 0, 0]]})
    n_e = len(arrays["tets"])
    doc, fields = lower_result("vm_plate_shear", problem, arrays, 0.0, np.zeros((n_e, 4, 6)))
    low = certified(problem, arrays, doc, fields)
    lower = low.certificate
    assert lower["value_certified"] == 0.0 and lower["floor"] == 0.0 and lower["checks"]["repair_mu"] == 0.0
    velocity = np.zeros_like(arrays["nodes"])
    velocity[:, 0] = arrays["nodes"][:, 1]
    footing = np.array([[1.0, 0, 0, 0, 0, 0]])
    # The contract's fixture tooling divides by W1, which a support query does not have: the value is D - W0 = k V.
    made = {"velocity": velocity, "footing_velocity": footing}
    doc, fields = make_result("vm_plate_shear", problem, arrays, "upper", 1.0, made, ct.CONTROL_ID)
    high = certified(problem, arrays, doc, fields)
    upper = high.certificate
    assert upper["value_certified"] == pytest.approx(1.0, rel=1e-12) and upper["checks"]["work_reference"] is None
    record = lc.bracket(lower, upper, low.method, high.method)
    ct.validate_bracket(record, lower, upper)
    assert record["status"] == "bracket" and record["gap"] is None
    # A plate velocity other than d is not this query's field.
    made = {"velocity": 2.0 * velocity, "footing_velocity": 2.0 * footing}
    doc, fields = make_result("vm_plate_shear", problem, arrays, "upper", 2.0, made, ct.CONTROL_ID)
    estimate(problem, arrays, doc, fields, "projection_move_too_large")


# ---------------------------------------------------------------------------------------------
# Public files
# ---------------------------------------------------------------------------------------------

def test_conditioning_command_writes_the_table(tmp_path):
    out = tmp_path / "table.json"
    assert lc.main(["conditioning", "--out", str(out), "--cases", "vm_plate_shear", "--meshes", "0"]) == 0
    document = ct.load_document(out)[0]
    assert document["format"] == lc.FORMAT_CONDITIONING and len(document["rows"]) == 3
    assert [row["formulation"] for row in document["rows"]] == CASES["vm_plate_shear"]["formulation"]
    assert all(row["capacity_refusal"] is None and row["kappa_eq"] >= 1.0 for row in document["rows"])
    refusing = lc.conditioning_table(["vm_plate_shear"], [0], dense_entries_max=0, sparse_rows_max=10)
    assert all(row["kappa_eq"] is None and "sparse factorisation limit" in row["capacity_refusal"] for row in refusing["rows"])


CONDITIONING_RECORD = Path(__file__).with_name("limit_conditioning_report.json")


def test_committed_conditioning_record_covers_every_gate_mesh_and_hides_no_refusal():
    """kappa_eq on meshes 0 to 5 of the eight true-3D cases and on the strips, as the command wrote it. A system the
    checker could not evaluate is in the record as a capacity refusal, with no number."""
    document = ct.load_document(CONDITIONING_RECORD)[0]
    assert document["format"] == lc.FORMAT_CONDITIONING and document["checker"]["name"] == lc.CHECKER_NAME
    hashes = ct.load_problem_hashes()
    gate_cases = [name for name in ct.ACTIVE if CASES[name]["class"] in lc.GATE_CLASSES]
    assert len(gate_cases) == 12 and len([name for name in gate_cases if CASES[name]["class"] == "upper_gate_3d"]) == 8
    expected = {(name, mesh, formulation) for name in gate_cases for mesh in range(len(CASES[name]["mesh"]))
                for formulation in CASES[name]["formulation"]}
    rows = {(row["case"], row["mesh"], row["formulation"]): row for row in document["rows"]}
    assert set(rows) == expected and len(document["rows"]) == len(expected) == 144
    evaluated, refused = 0, []
    for (name, mesh, formulation), row in rows.items():
        assert row["problem_sha256"] == hashes[name][f"m{mesh}"], "the record is of this case on this mesh"
        if row["kappa_eq"] is None:
            assert isinstance(row["capacity_refusal"], str) and row["capacity_refusal"] and row["rank"] is None
            assert row["method"] == "unavailable"
            refused.append((name, mesh, formulation))
        else:
            assert row["capacity_refusal"] is None and row["proof"] and row["method"] != "unavailable"
            assert 1.0 <= row["kappa_eq"] <= T["TOL_KAPPA_CAP"] and row["kappa_eq"] * lc.TAU_EQ <= T["TOL_FLOOR_MAX"]
            assert row["dependent_rows_skipped"] == row["reduced_rows"] - row["rank"] >= 0
            evaluated += 1
    assert evaluated + len(refused) == 144
    # What this implementation cannot reach is known and stated: the finest true-3D meshes.
    assert all(CASES[name]["class"] == "upper_gate_3d" and mesh == 5 for name, mesh, _ in refused), refused
    # One row is recomputed here, with the capacity the record states.
    key = ("upper3d_fixed_phi30", 0, ct.JUMP_ID)
    again = lc.conditioning_report(*case_problem(key[0], key[1]), "q", key[2], **document["capacity"])
    assert again["kappa_eq"] == pytest.approx(rows[key]["kappa_eq"], rel=1e-9) and again["rank"] == rows[key]["rank"]


def test_checker_files_carry_no_private_material():
    forbidden = ["TI" + "Ms", "Op" + "tum", "Fa" + "ble", "C:" + "\\Users", "nmo" + "ra"]
    for path in (Path(lc.__file__), Path(__file__), CONDITIONING_RECORD):
        text = path.read_text(encoding="utf-8")
        assert text.isascii(), path.name
        for word in forbidden:
            assert word not in text, f"{word!r} in {path.name}"
