"""Verification battery for the OpenSees-free CMS oracle.

These tests validate the mathematics and the data contracts required by a
future OpenSees implementation.  They do not claim to validate MPI transport or
the C++ parser; those are integration gates in ADR 1000.
"""

import numpy as np
import pytest
from scipy import linalg as sla

from _testbed.ladruno_cms_reference import (
    assembly_signature,
    assign_contribution_owner,
    assemble_global,
    build_cms_basis,
    build_spring_mass_elements,
    classify_interface,
    classify_effective_ownership,
    collective_compatibility_plan,
    compare_assembly_signatures,
    condense_coordinate_massless_dofs,
    compute_mac,
    contiguous_labels,
    estimate_dense_root_memory,
    refine_subspace_original_pencil,
    run_single_level_cms,
    run_single_level_matrices,
    run_paper_hierarchical_cms,
    run_two_level_cms,
    solve_coordinate_massless_pencil,
    subspace_mac,
)


pytestmark = [pytest.mark.zone_a]


N = 120
N_COARSE = 4
N_FINE_PER_COARSE = 3
TOTAL_FINE = N_COARSE * N_FINE_PER_COARSE


def reference(model):
    K, M = assemble_global(model)
    values, vectors = sla.eigh(K, M)
    return K, M, values, vectors


def test_lumped_chain_matches_closed_form():
    n = 30
    k = 2.5
    m = 1.7
    model = build_spring_mass_elements(n, k, m, "lumped")
    _, _, values, _ = reference(model)
    analytical = np.array(
        [4.0 * k / m * np.sin(j * np.pi / (2.0 * (n + 1))) ** 2 for j in range(1, n + 1)]
    )
    assert np.max(np.abs(values - analytical)) < 2.0e-14


def test_boundary_only_partition_is_exact():
    model = build_spring_mass_elements(3)
    K, M, values_ref, _ = reference(model)
    result = run_single_level_cms(model, n_sub=3, k_modes=1)
    assert all(len(i) == 0 for i in result["basis"].interiors.values())
    assert result["r_dim"] == 3
    assert np.max(np.abs(result["eigenvalues"] - values_ref)) < 1.0e-12
    assert np.max(result["residuals"]) < 1.0e-13
    assert np.allclose(result["T"] @ result["T"].T, np.eye(3), atol=1.0e-14)


def test_global_interface_blocks_are_present_once():
    model = build_spring_mass_elements(48, mass_type="consistent")
    K, M = assemble_global(model)
    labels = contiguous_labels(48, 4)
    basis = build_cms_basis(K, M, labels, k_modes=3)

    b = basis.boundary
    n_modal = sum(sub["k_actual"] for sub in basis.subdomains)
    expected_kbb = K[np.ix_(b, b)].copy()
    expected_mbb = M[np.ix_(b, b)].copy()
    has_mib = False
    for sub in basis.subdomains:
        interior = sub["interior"]
        if not len(interior):
            continue
        psi = sub["Psi"]
        K_BI = K[np.ix_(b, interior)]
        M_BI = M[np.ix_(b, interior)]
        M_IB = M[np.ix_(interior, b)]
        M_II = M[np.ix_(interior, interior)]
        has_mib = has_mib or np.linalg.norm(M_IB) > 0.0
        expected_kbb += K_BI @ psi
        expected_mbb += M_BI @ psi + psi.T @ M_IB + psi.T @ M_II @ psi

    assert has_mib, "the consistent-mass gate must exercise a nonzero M_IB"
    assert np.allclose(basis.K_star[n_modal:, n_modal:], expected_kbb, atol=2.0e-13)
    assert np.allclose(basis.M_star[n_modal:, n_modal:], expected_mbb, atol=2.0e-13)
    assert np.linalg.norm(basis.K_star[:n_modal, n_modal:], ord=np.inf) < 2.0e-12


@pytest.mark.parametrize("mass_type", ["lumped", "consistent"])
def test_single_level_full_basis_is_exact(mass_type):
    model = build_spring_mass_elements(72, mass_type=mass_type)
    K, M, values_ref, vectors_ref = reference(model)
    result = run_single_level_cms(model, n_sub=6, k_modes=None)

    assert result["r_dim"] == 72
    assert np.max(np.abs(result["eigenvalues"] - values_ref)) < 2.0e-11
    assert np.max(result["residuals"]) < 2.0e-12
    assert np.min(compute_mac(result["eigenvectors"], vectors_ref, M)) > 1.0 - 2.0e-10
    assert np.allclose(result["T"].T @ K @ result["T"], result["K_star"], atol=2.0e-12)
    assert np.allclose(result["T"].T @ M @ result["T"], result["M_star"], atol=2.0e-12)


def test_rayleigh_ritz_upper_bound_and_monotone_enrichment():
    model = build_spring_mass_elements(N, mass_type="consistent")
    _, _, values_ref, _ = reference(model)
    retained = [1, 2, 3, 5, 8]
    approximations = []
    for keep in retained:
        result = run_single_level_cms(model, TOTAL_FINE, keep)
        approximations.append(result["eigenvalues"][:10])

    approximations = np.asarray(approximations)
    # Ritz eigenvalues are upper bounds and decrease as the nested basis grows.
    assert np.all(approximations >= values_ref[:10] - 2.0e-12)
    assert np.all(np.diff(approximations, axis=0) <= 2.0e-11)
    rel_first_five = np.abs(approximations[-1, :5] - values_ref[:5]) / values_ref[:5]
    assert np.max(rel_first_five) < 2.0e-6


def test_two_level_performs_a_real_second_condensation():
    model = build_spring_mass_elements(N, mass_type="consistent")
    result = run_two_level_cms(
        model,
        n_coarse=N_COARSE,
        n_fine_per_coarse=N_FINE_PER_COARSE,
        k_per_fine=5,
        k_per_coarse=6,
    )
    assert result["r_dim"] < result["r_level2"] < result["n_original"]
    assert np.allclose(
        result["T"],
        result["fine_basis"].T @ result["coarse_basis"].T,
        atol=1.0e-14,
    )

    _, _, values_ref, _ = reference(model)
    rel = np.abs(result["eigenvalues"][:5] - values_ref[:5]) / values_ref[:5]
    assert np.max(rel) < 1.0e-3
    # This deliberately aggressive 120 -> 30 reduction is a structural gate,
    # not the production accuracy target.  The ADR requires adaptive enrichment
    # until a much tighter residual is reached.
    assert np.max(result["residuals"][:5]) < 1.5e-1


def test_two_level_full_bases_are_exact():
    model = build_spring_mass_elements(60, mass_type="consistent")
    K, M, values_ref, vectors_ref = reference(model)
    result = run_two_level_cms(
        model,
        n_coarse=3,
        n_fine_per_coarse=3,
        k_per_fine=None,
        k_per_coarse=None,
    )
    assert result["r_level2"] == 60
    assert result["r_dim"] == 60
    assert np.max(np.abs(result["eigenvalues"] - values_ref)) < 3.0e-11
    assert np.max(result["residuals"]) < 3.0e-12
    assert np.min(compute_mac(result["eigenvectors"], vectors_ref, M)) > 1.0 - 3.0e-10


def test_noncontiguous_partition_and_dense_spd_mass():
    rng = np.random.default_rng(20260721)
    n = 24
    A = rng.normal(size=(n, n))
    B = rng.normal(size=(n, n))
    K = A.T @ A + 3.0 * np.eye(n)
    M = B.T @ B + 2.0 * np.eye(n)
    labels = np.arange(n) % 4

    # Dense matrices make every equation an interface equation.  The reduction
    # must then become the identity and remain exact despite noncontiguous labels.
    result = run_single_level_matrices(K, M, labels, k_modes=2)
    values_ref, _ = sla.eigh(K, M)
    assert result["r_dim"] == n
    assert np.max(np.abs(result["eigenvalues"] - values_ref)) < 2.0e-12


def test_permuted_chain_has_genuinely_noncontiguous_interiors():
    model = build_spring_mass_elements(48, mass_type="consistent")
    K, M = assemble_global(model)
    labels = contiguous_labels(48, 4)
    rng = np.random.default_rng(7401)
    permutation = rng.permutation(48)
    Kp = K[np.ix_(permutation, permutation)]
    Mp = M[np.ix_(permutation, permutation)]
    labels_p = labels[permutation]

    result = run_single_level_matrices(Kp, Mp, labels_p, k_modes=None)
    values_ref, vectors_ref = sla.eigh(Kp, Mp)
    assert any(
        len(interior) > 2 and np.any(np.diff(interior) > 1)
        for interior in result["basis"].interiors.values()
    )
    assert np.max(np.abs(result["eigenvalues"] - values_ref)) < 2.0e-11
    assert np.max(result["residuals"]) < 2.0e-12
    assert np.min(compute_mac(result["eigenvectors"], vectors_ref, Mp)) > 1.0 - 2.0e-10


def test_mass_only_cross_partition_coupling_creates_interface():
    K = np.diag([2.0, 3.0, 4.0, 5.0])
    M = np.eye(4)
    M[1, 2] = M[2, 1] = 0.2
    labels = np.array([0, 0, 1, 1])
    boundary, interiors = classify_interface(K, M, labels)
    assert np.array_equal(boundary, [1, 2])
    assert np.array_equal(interiors[0], [0])
    assert np.array_equal(interiors[1], [3])

    result = run_single_level_matrices(K, M, labels, k_modes=None)
    values_ref, _ = sla.eigh(K, M)
    assert np.max(np.abs(result["eigenvalues"] - values_ref)) < 2.0e-12


def test_cross_coarse_contribution_is_owned_once_and_reaches_s1():
    # Preliminary equation labels put equation 2 in coarse group 1, while the
    # whole three-equation contribution is owned by fine subdomain 0 by
    # majority.  A second contribution owned by fine 1 also touches equation
    # 2, so effective ownership must promote it to the S1 interface.
    fine_labels = np.array([0, 0, 1, 1])
    contribution_ids = [[0, 1, 2], [2, 3], [-1, 0]]
    owners = [assign_contribution_owner(ids, fine_labels) for ids in contribution_ids]
    assert owners == [0, 1, 0]

    effective = classify_effective_ownership(contribution_ids, owners, coarse_of_fine=[0, 1])
    assert effective["level1_interface"] == {2}
    assert effective["level2_interface"] == set()
    assert effective["equation_owners"][2] == {0, 1}
    assert effective["interiors"][0] == {0, 1}
    assert effective["interiors"][1] == {3}

    # The tie rule is deterministic and therefore identical for addA/addM.
    assert assign_contribution_owner([1, 2], fine_labels) == 0


def test_clustered_modes_are_compared_as_a_subspace():
    theta = np.pi / 4.0
    A = np.eye(4)[:, :2]
    B = A @ np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    diagonal_mac = compute_mac(A, B, np.eye(4))
    cluster_mac = subspace_mac(A, B, np.eye(4))
    assert np.allclose(diagonal_mac, 0.5)
    assert np.allclose(cluster_mac, 1.0, atol=2.0e-14)


def test_base_cb_oracle_rejects_indefinite_mass_and_singular_interior():
    with pytest.raises(ValueError, match="positive definite"):
        run_single_level_matrices(np.eye(3), np.diag([1.0, -0.1, 1.0]), [0, 0, 0], None)

    with pytest.raises(ValueError, match="K_II is not positive definite"):
        run_single_level_matrices(np.diag([1.0, 0.0, 2.0]), np.eye(3), [0, 0, 0], None)


def test_coordinate_massless_condensation_recovers_finite_modes_and_full_residual():
    K = np.array(
        [
            [6.0, -2.0, 1.0],
            [-2.0, 5.0, -1.5],
            [1.0, -1.5, 4.0],
        ]
    )
    M = np.diag([2.0, 1.0, 0.0])
    result = solve_coordinate_massless_pencil(K, M)
    general_values = sla.eigvals(K, M)
    finite = np.sort(np.real(general_values[np.isfinite(general_values)]))

    assert np.allclose(result["eigenvalues"], finite, rtol=2.0e-13, atol=2.0e-13)
    assert np.max(result["residuals"]) < 2.0e-15
    assert np.linalg.norm(result["eigenvectors"][2, :]) > 0.0
    assert np.allclose(result["G"].T @ M @ result["G"], result["M"])


def test_multiple_coordinate_massless_dofs_are_reconstructed():
    K = np.array(
        [
            [8.0, -1.0, 2.0, 0.5],
            [-1.0, 7.0, -0.5, 1.5],
            [2.0, -0.5, 6.0, -1.0],
            [0.5, 1.5, -1.0, 5.0],
        ]
    )
    M = np.diag([3.0, 2.0, 0.0, 0.0])
    result = solve_coordinate_massless_pencil(K, M)
    assert np.array_equal(result["dynamic"], [0, 1])
    assert np.array_equal(result["massless"], [2, 3])
    assert np.max(result["residuals"]) < 3.0e-15


def test_massless_contract_rejects_general_nullspace_and_singular_kzz():
    with pytest.raises(ValueError, match="not coordinate aligned"):
        condense_coordinate_massless_dofs(
            np.eye(2), np.array([[0.0, 0.1], [0.1, 1.0]])
        )
    with pytest.raises(ValueError, match="dynamic mass block is not positive definite"):
        condense_coordinate_massless_dofs(
            np.eye(2), np.array([[1.0, -1.0], [-1.0, 1.0]])
        )
    with pytest.raises(ValueError, match="K_ZZ is not positive definite"):
        condense_coordinate_massless_dofs(
            np.array([[2.0, 0.0], [0.0, 0.0]]), np.diag([1.0, 0.0])
        )


def test_collective_compatibility_plan_handles_duplicates_and_empty_leader():
    groups = [
        [("mode", 0, 0), ("eq", 12)],
        [],
        [("eq", 12), ("mode", 2, 0), ("eq", 31)],
    ]
    plan = collective_compatibility_plan(groups)
    assert np.array_equal(plan["counts"], [2, 0, 3])
    assert np.array_equal(plan["displacements"], [0, 2, 2])
    assert plan["r_dup"] == 5
    assert plan["r"] == 4
    assert len(plan["local_maps"][1]) == 0
    assert plan["local_maps"][0][1] == plan["local_maps"][2][0]
    assert plan["unique_keys"] == sorted(plan["unique_keys"])


def test_dense_memory_preflight_matches_adr_formula():
    # p=8, k1=200 and |union B_g|=400 imply r=2000.  Use a duplicated
    # boundary count of 100 per leader for this explicit gather estimate.
    estimate = estimate_dense_root_memory([300] * 8, 2000, 50)
    assert estimate["matrix_and_vectors_bytes"] == 64_800_000
    assert estimate["static_tile_bytes"] == 0
    assert estimate["packed_entries"] == 361_200
    assert estimate["packed_gather_bytes"] == 8_668_800
    assert estimate["root_lower_bound_bytes"] == 73_468_800

    with_massless = estimate_dense_root_memory(
        [300] * 8, 2000, 50, massless_dimension=400, static_block_size=32
    )
    assert with_massless["static_tile_bytes"] == 102_400
    assert with_massless["root_lower_bound_bytes"] == 73_571_200


def test_assembly_signature_is_exact_in_structure_and_tolerant_in_values():
    values = np.array([[2.0, -1.0], [-1.0, 3.0]])
    reference = assembly_signature(7, "A", [4, -1], values, factor=0.5)
    close = assembly_signature(7, "A", [4, -1], values + 1.0e-13, factor=0.5)
    far = assembly_signature(7, "A", [4, -1], values + 1.0e-7, factor=0.5)
    wrong_ids = assembly_signature(7, "A", [5, -1], values, factor=0.5)

    assert compare_assembly_signatures(reference, close)
    assert not compare_assembly_signatures(reference, far)
    assert not compare_assembly_signatures(reference, wrong_ids)


@pytest.mark.parametrize("mass_type", ["lumped", "consistent"])
def test_paper_four_transformations_are_exact_with_full_local_bases(mass_type):
    model = build_spring_mass_elements(60, mass_type=mass_type)
    K, M, values_ref, vectors_ref = reference(model)
    result = run_paper_hierarchical_cms(model, 3, 3, None, None)

    assert result["dimensions"]["final"] == 60
    assert np.max(np.abs(result["eigenvalues"] - values_ref)) < 2.0e-11
    assert np.max(result["residuals"]) < 3.0e-12
    assert np.min(compute_mac(result["eigenvectors"], vectors_ref, M)) > 1.0 - 3.0e-10
    assert np.allclose(result["T"].T @ K @ result["T"], result["K_star"], atol=3.0e-12)
    assert np.allclose(result["T"].T @ M @ result["T"], result["M_star"], atol=3.0e-12)


def test_paper_compatibility_maps_merge_duplicates_exactly_once():
    result = run_paper_hierarchical_cms(
        build_spring_mass_elements(48, mass_type="consistent"), 3, 2, None, None
    )

    # Every duplicated coordinate selects exactly one independent coordinate.
    assert np.all(np.sum(result["S_level1"], axis=1) == 1.0)
    assert np.any(result["compatibility_counts_level1"] > 1)
    for assembled in result["level1_assembled"]:
        assert np.all(np.sum(assembled["S"], axis=1) == 1.0)
        assert np.any(assembled["counts"] > 1)

    # Reverse S1 -> level-1 CB -> S2 -> level-2 CB gives identical values on
    # every copy of a shared physical equation.
    assert result["max_duplicate_jump"] < 2.0e-14


def test_paper_hierarchy_executes_both_truncations_and_ritz_bounds():
    model = build_spring_mass_elements(120, mass_type="consistent")
    _, _, values_ref, _ = reference(model)
    result = run_paper_hierarchical_cms(model, 4, 3, 8, 10, num_modes=10)
    dims = result["dimensions"]

    assert dims["after_level1_cb_before_compatibility"] < dims["after_level2_compatibility"]
    assert dims["final"] < dims["original"]
    assert len(result["eigenvalues"]) == 10
    assert np.all(result["eigenvalues"] >= values_ref[:10] - 2.0e-12)
    assert np.max(np.abs(result["eigenvalues"][:5] - values_ref[:5]) / values_ref[:5]) < 2.0e-4


def test_global_subspace_refinement_reduces_original_residual():
    model = build_spring_mass_elements(120, mass_type="consistent")
    K, M = assemble_global(model)
    cms = run_paper_hierarchical_cms(model, 4, 3, 8, 10, num_modes=13)
    refined = refine_subspace_original_pencil(
        K, M, cms["eigenvectors"], number_of_modes=5,
        tolerance=1.0e-8, maximum_iterations=8,
    )

    assert cms["residuals"][:5].max() > 7.0e-2
    assert refined["converged"]
    assert refined["iterations"] == 5
    assert refined["residuals"].max() < 2.0e-9
    assert refined["history"][1]["maximum_residual"] < 5.0e-4
    assert max(item["m_orthogonality_loss"] for item in refined["history"]) < 2.0e-12


def test_reduced_only_iteration_is_a_negative_control():
    model = build_spring_mass_elements(120, mass_type="consistent")
    K, M = assemble_global(model)
    cms = run_paper_hierarchical_cms(model, 4, 3, 8, 10, num_modes=13)
    reduced = refine_subspace_original_pencil(
        cms["K_star"], cms["M_star"], cms["reduced_eigenvectors"],
        number_of_modes=5, tolerance=1.0e-30, maximum_iterations=3,
    )
    reconstructed = cms["T"] @ reduced["eigenvectors"]
    residual = np.max(
        np.linalg.norm(
            K @ reconstructed - (M @ reconstructed) * reduced["eigenvalues"], axis=0
        )
        / np.maximum(
            np.linalg.norm(K @ reconstructed, axis=0)
            + np.abs(reduced["eigenvalues"]) * np.linalg.norm(M @ reconstructed, axis=0),
            1.0e-30,
        )
    )

    assert reduced["iterations"] == 3
    assert abs(residual - cms["residuals"][:5].max()) < 2.0e-12


def test_bathe_iteration_vector_rule_accelerates_requested_modes():
    model = build_spring_mass_elements(120, mass_type="consistent")
    K, M = assemble_global(model)
    residuals = {}
    for q in (8, 16):
        cms = run_paper_hierarchical_cms(model, 4, 3, 8, 10, num_modes=q)
        refined = refine_subspace_original_pencil(
            K, M, cms["eigenvectors"], number_of_modes=8,
            tolerance=1.0e-14, maximum_iterations=8,
        )
        residuals[q] = refined["residuals"].max()

    assert residuals[8] > 1.0e-6
    assert residuals[16] < 1.0e-8
    assert residuals[16] < residuals[8] * 2.0e-4


def test_subspace_refinement_supports_coordinate_semidefinite_mass():
    K = np.array(
        [
            [8.0, -1.0, 2.0, 0.5],
            [-1.0, 7.0, -0.5, 1.5],
            [2.0, -0.5, 6.0, -1.0],
            [0.5, 1.5, -1.0, 5.0],
        ]
    )
    M = np.diag([3.0, 2.0, 0.0, 0.0])
    initial = np.random.default_rng(1000).normal(size=(4, 2))
    reference_result = solve_coordinate_massless_pencil(K, M)
    refined = refine_subspace_original_pencil(
        K, M, initial, number_of_modes=2,
        tolerance=1.0e-12, maximum_iterations=4,
    )

    assert refined["converged"]
    assert np.allclose(
        refined["eigenvalues"], reference_result["eigenvalues"],
        rtol=2.0e-12, atol=2.0e-12,
    )
    assert refined["residuals"].max() < 2.0e-12


def test_subspace_refinement_reports_mass_rank_loss_and_rigid_pencil():
    with pytest.raises(ValueError, match="lost rank in the M product"):
        refine_subspace_original_pencil(
            np.eye(2), np.diag([1.0, 0.0]), np.array([[0.0], [1.0]]),
            number_of_modes=1,
        )
    with pytest.raises(ValueError, match="K is not positive definite"):
        refine_subspace_original_pencil(
            np.diag([0.0, 2.0]), np.eye(2), np.eye(2), number_of_modes=1,
        )


def test_paper_hierarchy_rejects_impossible_partition_and_mode_request():
    model = build_spring_mass_elements(8)
    with pytest.raises(ValueError, match=r"level1\*level2"):
        run_paper_hierarchical_cms(model, 5, 2, 1, 1)
    with pytest.raises(ValueError, match="num_modes"):
        run_paper_hierarchical_cms(model, 2, 2, 0, 0, num_modes=2)
