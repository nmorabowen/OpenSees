"""OpenSees-free Craig-Bampton oracles for the CMS feasibility study.

The generic matrix oracle checks algebraic identities on an assembled,
constrained equation-space pencil.  The paper oracle additionally represents
the four transformations of Yu et al. explicitly: local level-2
Craig-Bampton, level-2 compatibility, local level-1 Craig-Bampton, and global
level-1 compatibility.  Neither oracle imports or emulates an OpenSees
eigensolver; they are independent numerical references for the proposed new
``LadrunoCMS`` library.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple, Union

import numpy as np
from scipy import linalg as sla


ModeSpec = Union[int, Sequence[int], Mapping[int, int], None]


@dataclass
class CMSBasis:
    """One Craig-Bampton transformation and its reduced pencil."""

    T: np.ndarray
    K_star: np.ndarray
    M_star: np.ndarray
    labels: np.ndarray
    boundary: np.ndarray
    interiors: Dict[int, np.ndarray]
    metadata: List[dict]
    subdomains: List[dict]


def condense_coordinate_massless_dofs(
    K: np.ndarray,
    M: np.ndarray,
    mass_rtol: float = 1.0e-12,
    mass_atol: float = 1.0e-14,
) -> dict:
    """Condense coordinate-aligned zero-mass DOFs from a symmetric pencil.

    This is the deliberately narrow v1 contract needed by frame/shell models:
    a massless coordinate must have a numerically zero *row and column* in M.
    A general nullspace is rejected.  The stiffness block on massless
    coordinates must be positive definite so that static reconstruction is
    unique.
    """

    K = np.asarray(K, dtype=float)
    M = np.asarray(M, dtype=float)
    if K.ndim != 2 or K.shape[0] != K.shape[1] or M.shape != K.shape:
        raise ValueError("K and M must be square matrices of equal size")
    if not np.allclose(K, K.T, rtol=mass_rtol, atol=mass_atol):
        raise ValueError("K must be symmetric")
    if not np.allclose(M, M.T, rtol=mass_rtol, atol=mass_atol):
        raise ValueError("M must be symmetric")

    scale = max(1.0, float(np.max(np.abs(M))))
    threshold = mass_atol + mass_rtol * scale
    diagonal = np.diag(M)
    if np.any(diagonal < -threshold):
        raise ValueError("M is not positive semidefinite")

    massless = np.flatnonzero(np.abs(diagonal) <= threshold)
    dynamic = np.flatnonzero(np.abs(diagonal) > threshold)
    if not len(dynamic):
        raise ValueError("the pencil has no dynamic coordinates")
    if len(massless) and np.max(np.abs(M[massless, :])) > threshold:
        raise ValueError("M nullspace is not coordinate aligned")

    M_dd = M[np.ix_(dynamic, dynamic)]
    try:
        sla.cholesky(M_dd, lower=True, check_finite=True)
    except sla.LinAlgError as exc:
        raise ValueError("dynamic mass block is not positive definite") from exc

    G = np.zeros((K.shape[0], len(dynamic)))
    G[dynamic, :] = np.eye(len(dynamic))
    if len(massless):
        K_zz = K[np.ix_(massless, massless)]
        K_zd = K[np.ix_(massless, dynamic)]
        try:
            factor = sla.cho_factor(K_zz, lower=True, check_finite=True)
        except sla.LinAlgError as exc:
            raise ValueError("massless stiffness block K_ZZ is not positive definite") from exc
        G[massless, :] = sla.cho_solve(factor, -K_zd, check_finite=True)

    K_hat = 0.5 * (G.T @ K @ G + (G.T @ K @ G).T)
    M_hat = 0.5 * (G.T @ M @ G + (G.T @ M @ G).T)
    return {
        "G": G,
        "K": K_hat,
        "M": M_hat,
        "dynamic": dynamic,
        "massless": massless,
    }


def solve_coordinate_massless_pencil(
    K: np.ndarray,
    M: np.ndarray,
    num_modes: int | None = None,
) -> dict:
    """Solve finite modes after static condensation and reconstruct all DOFs."""

    reduced = condense_coordinate_massless_dofs(K, M)
    n_dynamic = reduced["K"].shape[0]
    wanted = n_dynamic if num_modes is None else int(num_modes)
    if wanted < 1 or wanted > n_dynamic:
        raise ValueError("num_modes must satisfy 1 <= num_modes <= dynamic dimension")
    values, q = sla.eigh(
        reduced["K"], reduced["M"], subset_by_index=[0, wanted - 1]
    )
    vectors = reduced["G"] @ q
    return {
        **reduced,
        "eigenvalues": values,
        "eigenvectors": vectors,
        "residuals": compute_residual(np.asarray(K), np.asarray(M), values, vectors),
    }


def collective_compatibility_plan(key_groups: Sequence[Sequence[tuple]]) -> dict:
    """Pure oracle for the counts, maps and stable key union used by S2/S1."""

    counts = np.asarray([len(group) for group in key_groups], dtype=np.int64)
    displacements = np.zeros(len(counts), dtype=np.int64)
    if len(counts) > 1:
        displacements[1:] = np.cumsum(counts[:-1])
    flat = [key for group in key_groups for key in group]
    unique = sorted(set(flat))
    index = {key: position for position, key in enumerate(unique)}
    local_maps = [np.asarray([index[key] for key in group], dtype=np.int64) for group in key_groups]
    return {
        "counts": counts,
        "displacements": displacements,
        "flat_keys": flat,
        "unique_keys": unique,
        "local_maps": local_maps,
        "r_dup": len(flat),
        "r": len(unique),
    }


def estimate_dense_root_memory(
    local_dimensions: Sequence[int],
    dynamic_dimension: int,
    num_modes: int,
    massless_dimension: int = 0,
    static_block_size: int = 32,
    int_bytes: int = 4,
    double_bytes: int = 8,
) -> dict:
    """Return the ADR's deterministic lower-bound memory preflight."""

    r = int(dynamic_dimension)
    r_z = int(massless_dimension)
    block_size = int(static_block_size)
    n_modes = int(num_modes)
    dims = [int(value) for value in local_dimensions]
    if r < 0 or r_z < 0 or block_size < 1 or n_modes < 0 or any(value < 0 for value in dims):
        raise ValueError("dimensions must be nonnegative")
    matrix_and_vectors = 2 * double_bytes * r * r + double_bytes * r * n_modes
    static_tile = double_bytes * r_z * min(block_size, max(1, r))
    packed_entries = sum(value * (value + 1) // 2 for value in dims)
    packed_gather = (2 * int_bytes + 2 * double_bytes) * packed_entries
    return {
        "matrix_and_vectors_bytes": matrix_and_vectors,
        "static_tile_bytes": static_tile,
        "packed_gather_bytes": packed_gather,
        "root_lower_bound_bytes": matrix_and_vectors + static_tile + packed_gather,
        "packed_entries": packed_entries,
    }


def assembly_signature(
    ordinal: int,
    kind: str,
    ids: Sequence[int],
    values: np.ndarray,
    factor: float = 1.0,
) -> dict:
    """Compact per-call assembly signature computed before owner filtering."""

    matrix = np.asarray(values, dtype=float)
    active = tuple(int(value) for value in ids if int(value) >= 0)
    structural_text = repr((int(ordinal), str(kind), tuple(int(v) for v in ids), matrix.shape))
    structural_fields = (int(ordinal), str(kind), tuple(int(v) for v in ids), matrix.shape)
    structural = hashlib.sha256(structural_text.encode("utf-8")).hexdigest()
    scaled = (float(factor) * matrix).ravel()
    weights = ((np.arange(scaled.size, dtype=float) * 131.0 + 17.0) % 257.0 + 1.0) / 257.0
    metrics = np.array(
        [
            np.max(np.abs(scaled)) if scaled.size else 0.0,
            np.sum(scaled),
            np.sum(np.abs(scaled)),
            np.dot(scaled, scaled),
            np.dot(weights, scaled),
        ],
        dtype=float,
    )
    return {
        "ordinal": int(ordinal),
        "active_ids": active,
        "structural_fields": structural_fields,
        "structural": structural,
        "metrics": metrics,
        "n_values": int(scaled.size),
    }


def compare_assembly_signatures(
    reference: dict,
    candidate: dict,
    rtol: float = 1.0e-12,
    atol: float = 1.0e-14,
) -> bool:
    """Compare exact structure and tolerant numerical invariants."""

    if reference["structural_fields"] != candidate["structural_fields"]:
        return False
    if reference["n_values"] != candidate["n_values"]:
        return False
    scale = np.maximum(np.abs(reference["metrics"]), np.abs(candidate["metrics"]))
    tolerance = atol * max(1, reference["n_values"]) + rtol * np.maximum(1.0, scale)
    return bool(np.all(np.abs(reference["metrics"] - candidate["metrics"]) <= tolerance))


# ---------------------------------------------------------------------------
# Small deterministic models
# ---------------------------------------------------------------------------


def build_spring_mass_elements(
    n: int,
    k: float = 1.0,
    m: float = 1.0,
    mass_type: str = "lumped",
) -> dict:
    """Describe a fixed-fixed 1D spring chain with ``n`` free equations.

    ``mass_type='consistent'`` uses the standard assembled two-node consistent
    mass stencil.  It deliberately creates nonzero M_IB blocks, which the old
    lumped-only oracle did not exercise.
    """

    if n < 1:
        raise ValueError("n must be positive")
    if k <= 0.0 or m <= 0.0:
        raise ValueError("k and m must be positive")
    if mass_type not in {"lumped", "consistent"}:
        raise ValueError("mass_type must be 'lumped' or 'consistent'")
    return {"global_n": n, "k": float(k), "m": float(m), "mass_type": mass_type}


def assemble_global(model: dict) -> Tuple[np.ndarray, np.ndarray]:
    """Assemble K and M for :func:`build_spring_mass_elements`."""

    n = int(model["global_n"])
    k = float(model["k"])
    m = float(model["m"])

    K = np.diag(np.full(n, 2.0 * k))
    if n > 1:
        K += np.diag(np.full(n - 1, -k), 1)
        K += np.diag(np.full(n - 1, -k), -1)

    if model.get("mass_type", "lumped") == "lumped":
        M = np.eye(n) * m
    else:
        # Fixed-fixed uniform two-node chain: each free node receives two
        # element diagonal terms 2m/6 and each adjacent pair one m/6 term.
        M = np.diag(np.full(n, 4.0 * m / 6.0))
        if n > 1:
            M += np.diag(np.full(n - 1, m / 6.0), 1)
            M += np.diag(np.full(n - 1, m / 6.0), -1)
    return K, M


def contiguous_labels(n: int, n_sub: int) -> np.ndarray:
    """Balanced contiguous partition labels in ``[0, n_sub)``."""

    if n < 1 or n_sub < 1 or n_sub > n:
        raise ValueError("require 1 <= n_sub <= n")
    return np.minimum((np.arange(n, dtype=int) * n_sub) // n, n_sub - 1)


# ---------------------------------------------------------------------------
# Matrix and partition contracts
# ---------------------------------------------------------------------------


def _validate_pencil(K: np.ndarray, M: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    K = np.asarray(K, dtype=float)
    M = np.asarray(M, dtype=float)
    if K.ndim != 2 or K.shape[0] != K.shape[1] or M.shape != K.shape:
        raise ValueError("K and M must be square matrices with the same shape")
    scale_k = max(1.0, float(np.linalg.norm(K, ord=np.inf)))
    scale_m = max(1.0, float(np.linalg.norm(M, ord=np.inf)))
    if not np.allclose(K, K.T, rtol=0.0, atol=1.0e-12 * scale_k):
        raise ValueError("K must be symmetric")
    if not np.allclose(M, M.T, rtol=0.0, atol=1.0e-12 * scale_m):
        raise ValueError("M must be symmetric")
    try:
        sla.cholesky(M, lower=True, check_finite=True)
    except sla.LinAlgError as exc:
        raise ValueError("M must be positive definite for fixed-interface CMS") from exc
    return K, M


def _normalise_labels(labels: Iterable[int], n: int) -> np.ndarray:
    labels = np.asarray(list(labels), dtype=int)
    if labels.shape != (n,):
        raise ValueError("one partition label is required per equation")
    unique = np.unique(labels)
    remap = {int(old): new for new, old in enumerate(unique)}
    return np.array([remap[int(x)] for x in labels], dtype=int)


def assign_contribution_owner(active_ids: Iterable[int], fine_labels: Iterable[int]) -> int:
    """Apply the ADR's deterministic whole-contribution ownership rule.

    Negative OpenSees equation IDs are ignored.  The contribution goes to the
    fine subdomain containing the most active IDs; ties go to the lowest label.
    Keeping the contribution whole is what prevents duplicate K_BB/M_BB terms.
    """

    labels = np.asarray(list(fine_labels), dtype=int)
    ids = np.asarray([int(i) for i in active_ids if int(i) >= 0], dtype=int)
    if not len(ids):
        raise ValueError("a contribution needs at least one active equation ID")
    if np.any(ids >= len(labels)):
        raise ValueError("active equation ID is outside the partition map")
    owners, counts = np.unique(labels[ids], return_counts=True)
    return int(np.min(owners[counts == np.max(counts)]))


def classify_effective_ownership(
    contribution_ids: Sequence[Iterable[int]],
    contribution_owners: Sequence[int],
    coarse_of_fine: Sequence[int],
) -> dict:
    """Classify interfaces from actual contribution owners, not METIS labels.

    An equation seen by multiple fine owners in one coarse group belongs to
    S2.  If those owners span coarse groups it belongs to S1.  An equation seen
    by only one owner migrates to that owner and remains interior there.
    """

    if len(contribution_ids) != len(contribution_owners):
        raise ValueError("one owner is required per contribution")
    coarse = np.asarray(coarse_of_fine, dtype=int)
    equation_owners: Dict[int, set] = {}
    for ids, owner in zip(contribution_ids, contribution_owners):
        owner = int(owner)
        if owner < 0 or owner >= len(coarse):
            raise ValueError("contribution owner is outside coarse_of_fine")
        for equation in {int(i) for i in ids if int(i) >= 0}:
            equation_owners.setdefault(equation, set()).add(owner)

    level2 = set()
    level1 = set()
    interiors: Dict[int, set] = {fine: set() for fine in range(len(coarse))}
    for equation, owners in equation_owners.items():
        if len(owners) == 1:
            interiors[next(iter(owners))].add(equation)
        elif len({int(coarse[owner]) for owner in owners}) > 1:
            level1.add(equation)
        else:
            level2.add(equation)
    return {
        "equation_owners": equation_owners,
        "interiors": interiors,
        "level2_interface": level2,
        "level1_interface": level1,
    }


def classify_interface(
    K: np.ndarray,
    M: np.ndarray,
    labels: Iterable[int],
    pattern_tol: float = 0.0,
) -> Tuple[np.ndarray, Dict[int, np.ndarray]]:
    """Classify equation-space interiors and the global vertex interface.

    An equation is on the interface when it participates in a nonzero K or M
    coupling whose other endpoint has a different partition label.  Both
    endpoints become boundary equations.  This is a vertex-separator contract
    over the union pattern of K and M.
    """

    K, M = _validate_pencil(K, M)
    labels = _normalise_labels(labels, K.shape[0])
    pattern = (np.abs(K) > pattern_tol) | (np.abs(M) > pattern_tol)
    np.fill_diagonal(pattern, False)
    cross = pattern & (labels[:, None] != labels[None, :])
    boundary = np.flatnonzero(np.any(cross, axis=1))
    bset = set(int(i) for i in boundary)
    interiors = {
        s: np.array([i for i in np.flatnonzero(labels == s) if int(i) not in bset], dtype=int)
        for s in range(int(labels.max()) + 1)
    }
    return boundary, interiors


def _mode_count(spec: ModeSpec, subdomain: int, n_interior: int) -> int:
    if spec is None:
        return n_interior
    if isinstance(spec, Mapping):
        requested = int(spec.get(subdomain, 0))
    elif np.isscalar(spec):
        requested = int(spec)
    else:
        requested = int(spec[subdomain])
    if requested < 0:
        raise ValueError("retained mode counts must be nonnegative")
    return min(requested, n_interior)


def build_cms_basis(
    K: np.ndarray,
    M: np.ndarray,
    labels: Iterable[int],
    k_modes: ModeSpec,
    pattern_tol: float = 0.0,
) -> CMSBasis:
    """Build a fixed-interface Craig-Bampton basis in equation space.

    The global transformation is

    ``T = [[blockdiag(Phi_s), Psi], [0, I_B]]``.

    Reduced matrices are formed by a global congruence.  Therefore global
    K_BB/M_BB occur once, independently of how interface contributions were
    originally assembled by OpenSees.
    """

    K, M = _validate_pencil(K, M)
    n = K.shape[0]
    labels = _normalise_labels(labels, n)
    boundary, interiors = classify_interface(K, M, labels, pattern_tol)
    n_sub = int(labels.max()) + 1

    columns: List[np.ndarray] = []
    metadata: List[dict] = []
    subdomains: List[dict] = []

    for s in range(n_sub):
        interior = interiors[s]
        n_i = len(interior)
        keep = _mode_count(k_modes, s, n_i)
        K_II = K[np.ix_(interior, interior)]
        M_II = M[np.ix_(interior, interior)]
        K_IB = K[np.ix_(interior, boundary)]

        if n_i:
            try:
                chol_k = sla.cho_factor(K_II, lower=True, check_finite=True)
            except sla.LinAlgError as exc:
                raise ValueError(
                    f"K_II is not positive definite in subdomain {s}; "
                    "fixed-interface CMS needs an invertible interior block"
                ) from exc
            Psi = sla.cho_solve(chol_k, -K_IB, check_finite=True)
        else:
            Psi = np.zeros((0, len(boundary)))

        if keep:
            evals, Phi = sla.eigh(K_II, M_II, subset_by_index=[0, keep - 1])
            for local_mode in range(keep):
                col = np.zeros(n)
                col[interior] = Phi[:, local_mode]
                columns.append(col)
                metadata.append(
                    {
                        "kind": "modal",
                        "subdomain": s,
                        "local_mode": local_mode,
                        "interior_eigenvalue": float(evals[local_mode]),
                    }
                )
        else:
            evals = np.empty(0)
            Phi = np.zeros((n_i, 0))

        subdomains.append(
            {
                "subdomain": s,
                "interior": interior,
                "boundary": boundary,
                "Phi": Phi,
                "Psi": Psi,
                "Lambda": evals,
                "k_actual": keep,
            }
        )

    for ib, dof in enumerate(boundary):
        col = np.zeros(n)
        col[dof] = 1.0
        for sub in subdomains:
            if len(sub["interior"]):
                col[sub["interior"]] = sub["Psi"][:, ib]
        columns.append(col)
        metadata.append({"kind": "boundary", "equation": int(dof)})

    T = np.column_stack(columns) if columns else np.zeros((n, 0))
    K_star = T.T @ K @ T
    M_star = T.T @ M @ T
    # Remove harmless round-off asymmetry before passing the pencil to LAPACK.
    K_star = 0.5 * (K_star + K_star.T)
    M_star = 0.5 * (M_star + M_star.T)
    return CMSBasis(T, K_star, M_star, labels, boundary, interiors, metadata, subdomains)


def solve_reduced(
    basis: CMSBasis,
    num_modes: int | None = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Solve a reduced pencil and back-substitute full eigenvectors."""

    r = basis.K_star.shape[0]
    if r == 0:
        return np.empty(0), np.empty((0, 0)), np.empty((basis.T.shape[0], 0))
    wanted = r if num_modes is None else min(int(num_modes), r)
    if wanted < 1:
        raise ValueError("num_modes must be positive")
    values, reduced = sla.eigh(
        basis.K_star,
        basis.M_star,
        subset_by_index=[0, wanted - 1],
    )
    full = basis.T @ reduced
    return values, reduced, full


# ---------------------------------------------------------------------------
# Accuracy metrics
# ---------------------------------------------------------------------------


def compute_residual(
    K: np.ndarray,
    M: np.ndarray,
    eigenvalues: np.ndarray,
    eigenvectors: np.ndarray,
) -> np.ndarray:
    """Scale-invariant generalized-eigenpair backward residual."""

    residuals = np.zeros(len(eigenvalues))
    for j, value in enumerate(eigenvalues):
        v = eigenvectors[:, j]
        Kv = K @ v
        Mv = M @ v
        residuals[j] = np.linalg.norm(Kv - value * Mv) / max(
            np.linalg.norm(Kv) + abs(value) * np.linalg.norm(Mv), 1.0e-30
        )
    return residuals


def refine_subspace_original_pencil(
    K: np.ndarray,
    M: np.ndarray,
    initial_vectors: np.ndarray,
    number_of_modes: int,
    tolerance: float = 1.0e-8,
    maximum_iterations: int = 20,
) -> dict:
    """Refine a CMS starting space by Bathe subspace iteration on ``(K, M)``.

    Rayleigh--Ritz extraction and mass orthonormalization are performed at
    every iteration.  A semidefinite mass is admitted only when the supplied
    block has a positive-definite Gram matrix.  The stiffness must be positive
    definite; shifted/free-free pencils are outside the v1 contract.
    """

    stiffness = np.asarray(K, dtype=float)
    mass = np.asarray(M, dtype=float)
    vectors = np.asarray(initial_vectors, dtype=float).copy()
    if stiffness.ndim != 2 or stiffness.shape[0] != stiffness.shape[1]:
        raise ValueError("K must be square")
    if mass.shape != stiffness.shape:
        raise ValueError("M must have the same shape as K")
    if vectors.ndim != 2 or vectors.shape[0] != stiffness.shape[0]:
        raise ValueError("initial vectors have incompatible dimensions")
    if not np.allclose(stiffness, stiffness.T, rtol=1.0e-12, atol=1.0e-14):
        raise ValueError("K must be symmetric")
    if not np.allclose(mass, mass.T, rtol=1.0e-12, atol=1.0e-14):
        raise ValueError("M must be symmetric")
    requested = int(number_of_modes)
    iterations = int(maximum_iterations)
    if requested < 1 or requested > vectors.shape[1]:
        raise ValueError("number_of_modes must satisfy 1 <= p <= q")
    if not np.isfinite(tolerance) or tolerance <= 0.0:
        raise ValueError("tolerance must be positive and finite")
    if iterations < 0:
        raise ValueError("maximum_iterations must be nonnegative")

    try:
        stiffness_factor = sla.cho_factor(stiffness, lower=True, check_finite=True)
    except sla.LinAlgError as error:
        raise ValueError("K is not positive definite") from error

    def rayleigh_ritz(block: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
        gram = 0.5 * (block.T @ mass @ block + (block.T @ mass @ block).T)
        try:
            lower = sla.cholesky(gram, lower=True, check_finite=True)
        except sla.LinAlgError as error:
            raise ValueError("iteration block lost rank in the M product") from error
        inverse_transpose = sla.solve_triangular(
            lower.T,
            np.eye(lower.shape[0]),
            lower=False,
            check_finite=True,
        )
        orthonormal = block @ inverse_transpose
        projected = orthonormal.T @ stiffness @ orthonormal
        projected = 0.5 * (projected + projected.T)
        values, coordinates = sla.eigh(projected, check_finite=True)
        ritz = orthonormal @ coordinates
        orthogonality = np.linalg.norm(
            ritz.T @ mass @ ritz - np.eye(ritz.shape[1]), ord=np.inf
        )
        return values, ritz, float(orthogonality)

    values, vectors, orthogonality = rayleigh_ritz(vectors)
    residuals = compute_residual(stiffness, mass, values[:requested], vectors[:, :requested])
    history = [
        {
            "iteration": 0,
            "maximum_residual": float(np.max(residuals)),
            "m_orthogonality_loss": orthogonality,
        }
    ]
    converged = bool(np.max(residuals) <= tolerance)
    completed = 0
    while not converged and completed < iterations:
        right_hand_sides = mass @ vectors
        inverse_block = sla.cho_solve(
            stiffness_factor, right_hand_sides, check_finite=True
        )
        values, vectors, orthogonality = rayleigh_ritz(inverse_block)
        residuals = compute_residual(
            stiffness, mass, values[:requested], vectors[:, :requested]
        )
        completed += 1
        history.append(
            {
                "iteration": completed,
                "maximum_residual": float(np.max(residuals)),
                "m_orthogonality_loss": orthogonality,
            }
        )
        converged = bool(np.max(residuals) <= tolerance)

    return {
        "eigenvalues": values[:requested],
        "eigenvectors": vectors[:, :requested],
        "iteration_values": values,
        "iteration_vectors": vectors,
        "residuals": residuals,
        "history": history,
        "iterations": completed,
        "converged": converged,
        "number_of_iteration_vectors": vectors.shape[1],
    }


def compute_mac(v_cms: np.ndarray, v_ref: np.ndarray, M: np.ndarray | None = None) -> np.ndarray:
    """Diagonal Euclidean or mass-weighted MAC for simple eigenvalues."""

    n_modes = min(v_cms.shape[1], v_ref.shape[1])
    metric = np.eye(v_cms.shape[0]) if M is None else np.asarray(M, dtype=float)
    mac = np.zeros(n_modes)
    for k in range(n_modes):
        a = v_cms[:, k]
        b = v_ref[:, k]
        num = abs(a.T @ metric @ b) ** 2
        den = (a.T @ metric @ a) * (b.T @ metric @ b)
        mac[k] = num / max(float(den), 1.0e-30)
    return mac


def subspace_mac(A: np.ndarray, B: np.ndarray, M: np.ndarray) -> np.ndarray:
    """Squared principal cosines between two modal subspaces.

    This is the appropriate comparison for repeated or tightly clustered modes,
    whose individual vectors may rotate while spanning the same eigenspace.
    """

    gram_a = A.T @ M @ A
    gram_b = B.T @ M @ B
    Qa = A @ sla.inv(sla.sqrtm(gram_a))
    Qb = B @ sla.inv(sla.sqrtm(gram_b))
    singular = sla.svdvals(Qa.T @ M @ Qb)
    return np.clip(singular**2, 0.0, 1.0)


# ---------------------------------------------------------------------------
# High-level single- and two-level drivers
# ---------------------------------------------------------------------------


def _result(
    K: np.ndarray,
    M: np.ndarray,
    basis: CMSBasis,
    num_modes: int | None,
    **extra,
) -> dict:
    values, reduced, vectors = solve_reduced(basis, num_modes)
    data = {
        "eigenvalues": values,
        "eigenvectors": vectors,
        "reduced_eigenvectors": reduced,
        "residuals": compute_residual(K, M, values, vectors),
        "K_star": basis.K_star,
        "M_star": basis.M_star,
        "T": basis.T,
        "basis": basis,
        "r_dim": basis.K_star.shape[0],
        "n_original": K.shape[0],
    }
    data.update(extra)
    return data


def run_single_level_matrices(
    K: np.ndarray,
    M: np.ndarray,
    labels: Iterable[int],
    k_modes: ModeSpec,
    num_modes: int | None = None,
) -> dict:
    """Single-level CMS directly on an assembled equation-space pencil."""

    K, M = _validate_pencil(K, M)
    basis = build_cms_basis(K, M, labels, k_modes)
    return _result(K, M, basis, num_modes, level=1)


def run_single_level_cms(model: dict, n_sub: int, k_modes: ModeSpec) -> dict:
    """Compatibility driver for the deterministic spring-chain tests."""

    K, M = assemble_global(model)
    labels = contiguous_labels(K.shape[0], n_sub)
    return run_single_level_matrices(K, M, labels, k_modes)


def _coarse_coordinate_labels(
    fine_basis: CMSBasis,
    coarse_of_fine: Sequence[int],
) -> np.ndarray:
    coarse_of_fine = np.asarray(coarse_of_fine, dtype=int)
    n_fine = int(fine_basis.labels.max()) + 1
    if coarse_of_fine.shape != (n_fine,):
        raise ValueError("coarse_of_fine must map every fine subdomain")
    result = []
    for item in fine_basis.metadata:
        if item["kind"] == "modal":
            fine = int(item["subdomain"])
        else:
            # A global boundary coordinate is assigned to its primary fine
            # partition.  Couplings to another coarse label make it a coarse
            # boundary automatically when the reduced graph is classified.
            fine = int(fine_basis.labels[int(item["equation"])])
        result.append(int(coarse_of_fine[fine]))
    return np.asarray(result, dtype=int)


def run_two_level_matrices(
    K: np.ndarray,
    M: np.ndarray,
    fine_labels: Iterable[int],
    coarse_of_fine: Sequence[int],
    k_fine: ModeSpec,
    k_coarse: ModeSpec,
    num_modes: int | None = None,
) -> dict:
    """Genuine two-level CMS with ``T_total = T_fine @ T_coarse``."""

    K, M = _validate_pencil(K, M)
    fine = build_cms_basis(K, M, fine_labels, k_fine)
    coarse_labels = _coarse_coordinate_labels(fine, coarse_of_fine)
    coarse = build_cms_basis(fine.K_star, fine.M_star, coarse_labels, k_coarse)

    T_total = fine.T @ coarse.T
    final = CMSBasis(
        T=T_total,
        K_star=coarse.K_star,
        M_star=coarse.M_star,
        labels=coarse.labels,
        boundary=coarse.boundary,
        interiors=coarse.interiors,
        metadata=coarse.metadata,
        subdomains=coarse.subdomains,
    )
    return _result(
        K,
        M,
        final,
        num_modes,
        level=2,
        fine_basis=fine,
        coarse_basis=coarse,
        r_level2=fine.K_star.shape[0],
    )


def run_two_level_cms(
    model: dict,
    n_coarse: int,
    n_fine_per_coarse: int,
    k_per_fine: ModeSpec,
    k_per_coarse: ModeSpec | None = None,
) -> dict:
    """Two-level spring-chain driver.

    ``k_per_coarse=None`` means a full coarse interior basis.  Use an integer to
    exercise a genuine second truncation and obtain ``r_dim < r_level2``.
    """

    K, M = assemble_global(model)
    total_fine = n_coarse * n_fine_per_coarse
    fine_labels = contiguous_labels(K.shape[0], total_fine)
    coarse_of_fine = np.repeat(np.arange(n_coarse, dtype=int), n_fine_per_coarse)
    return run_two_level_matrices(
        K,
        M,
        fine_labels,
        coarse_of_fine,
        k_per_fine,
        k_per_coarse,
    )


# ---------------------------------------------------------------------------
# Paper-faithful four-transformation hierarchy
# ---------------------------------------------------------------------------


def _chain_elements(model: dict) -> List[Tuple[Tuple[int | None, int | None], np.ndarray, np.ndarray]]:
    """Return fixed-fixed chain elements before assembly.

    ``None`` denotes an eliminated support equation.  Keeping element matrices
    separate is essential here: compatibility matrices must merge duplicated
    physical interface coordinates exactly once.
    """

    n = int(model["global_n"])
    k = float(model["k"])
    m = float(model["m"])
    ke = k * np.array([[1.0, -1.0], [-1.0, 1.0]])
    if model.get("mass_type", "lumped") == "lumped":
        me = 0.5 * m * np.eye(2)
    else:
        me = (m / 6.0) * np.array([[2.0, 1.0], [1.0, 2.0]])

    nodes: List[Tuple[int | None, int | None]] = [(None, 0)]
    nodes.extend((i - 1, i) for i in range(1, n))
    nodes.append((n - 1, None))
    return [(pair, ke.copy(), me.copy()) for pair in nodes]


def _balanced_element_labels(n_elements: int, n_parts: int) -> np.ndarray:
    if n_parts < 1 or n_parts > n_elements:
        raise ValueError("hierarchy requires 1 <= level1*level2 <= number of elements")
    return np.minimum((np.arange(n_elements, dtype=int) * n_parts) // n_elements, n_parts - 1)


def _assemble_element_subset(elements, element_ids: Sequence[int]) -> Tuple[List[int], np.ndarray, np.ndarray]:
    nodes = sorted(
        {
            int(node)
            for eid in element_ids
            for node in elements[int(eid)][0]
            if node is not None
        }
    )
    lookup = {node: i for i, node in enumerate(nodes)}
    K = np.zeros((len(nodes), len(nodes)))
    M = np.zeros_like(K)
    for eid in element_ids:
        endpoints, ke, me = elements[int(eid)]
        for a, node_a in enumerate(endpoints):
            if node_a is None:
                continue
            ia = lookup[int(node_a)]
            for b, node_b in enumerate(endpoints):
                if node_b is None:
                    continue
                ib = lookup[int(node_b)]
                K[ia, ib] += ke[a, b]
                M[ia, ib] += me[a, b]
    return nodes, K, M


def _local_cb(
    K: np.ndarray,
    M: np.ndarray,
    interior: Sequence[int],
    boundary: Sequence[int],
    keep_spec: ModeSpec,
    subdomain: int,
) -> dict:
    """Local fixed-interface Craig-Bampton transform in ``[modes, B]`` order."""

    interior = np.asarray(interior, dtype=int)
    boundary = np.asarray(boundary, dtype=int)
    n_i = len(interior)
    keep = _mode_count(keep_spec, subdomain, n_i)
    K_II = K[np.ix_(interior, interior)]
    M_II = M[np.ix_(interior, interior)]
    K_IB = K[np.ix_(interior, boundary)]

    if n_i:
        try:
            factor = sla.cho_factor(K_II, lower=True, check_finite=True)
        except sla.LinAlgError as exc:
            raise ValueError(
                f"K_II is not positive definite in hierarchical subdomain {subdomain}"
            ) from exc
        psi = sla.cho_solve(factor, -K_IB, check_finite=True)
    else:
        psi = np.zeros((0, len(boundary)))

    if keep:
        values, phi = sla.eigh(K_II, M_II, subset_by_index=[0, keep - 1])
    else:
        values = np.empty(0)
        phi = np.zeros((n_i, 0))

    T = np.zeros((K.shape[0], keep + len(boundary)))
    if keep:
        T[np.ix_(interior, np.arange(keep))] = phi
    if len(boundary):
        bcols = np.arange(keep, keep + len(boundary))
        if n_i:
            T[np.ix_(interior, bcols)] = psi
        T[np.ix_(boundary, bcols)] = np.eye(len(boundary))

    K_star = 0.5 * (T.T @ K @ T + (T.T @ K @ T).T)
    M_star = 0.5 * (T.T @ M @ T + (T.T @ M @ T).T)
    return {
        "T": T,
        "K_star": K_star,
        "M_star": M_star,
        "interior": interior,
        "boundary": boundary,
        "Phi": phi,
        "Psi": psi,
        "Lambda": values,
        "k_actual": keep,
    }


def _compatibility(key_groups: Sequence[Sequence[tuple]]) -> Tuple[np.ndarray, List[tuple], np.ndarray]:
    """Build the Boolean compatibility map from unique to duplicated coordinates."""

    unique: List[tuple] = []
    index: Dict[tuple, int] = {}
    flat = [key for group in key_groups for key in group]
    for key in flat:
        if key not in index:
            index[key] = len(unique)
            unique.append(key)
    S = np.zeros((len(flat), len(unique)))
    for row, key in enumerate(flat):
        S[row, index[key]] = 1.0
    counts = np.asarray([flat.count(key) for key in unique], dtype=int)
    return S, unique, counts


def _block_diag(items: Sequence[np.ndarray]) -> np.ndarray:
    if not items:
        return np.zeros((0, 0))
    return sla.block_diag(*items)


def run_paper_hierarchical_cms(
    model: dict,
    n_level1: int,
    n_level2_per_level1: int,
    k_level2: ModeSpec,
    k_level1: ModeSpec,
    num_modes: int | None = None,
) -> dict:
    """Execute the four transformations specified by Yu et al.

    Element ownership is disjoint.  Local interface coordinates are duplicated
    only until the corresponding Boolean compatibility transformation merges
    them.  The reverse path constructs full equation-space eigenvectors, which
    mirrors the output contract required by ``AnalysisModel::setEigenvector``.
    """

    K_global, M_global = assemble_global(model)
    elements = _chain_elements(model)
    n_fine = int(n_level1) * int(n_level2_per_level1)
    fine_labels = _balanced_element_labels(len(elements), n_fine)
    coarse_of_fine = np.repeat(np.arange(n_level1, dtype=int), n_level2_per_level1)
    if len(coarse_of_fine) != n_fine:
        raise ValueError("level hierarchy is inconsistent")

    node_fine_owners: Dict[int, set] = {}
    node_coarse_owners: Dict[int, set] = {}
    for eid, fine in enumerate(fine_labels):
        coarse = int(coarse_of_fine[int(fine)])
        for node in elements[eid][0]:
            if node is None:
                continue
            node_fine_owners.setdefault(int(node), set()).add(int(fine))
            node_coarse_owners.setdefault(int(node), set()).add(coarse)
    fine_boundary_nodes = {node for node, owners in node_fine_owners.items() if len(owners) > 1}
    coarse_boundary_nodes = {node for node, owners in node_coarse_owners.items() if len(owners) > 1}

    # Transformation 1: independent Craig-Bampton reductions on level 2.
    fine_data: List[dict] = []
    for fine in range(n_fine):
        element_ids = np.flatnonzero(fine_labels == fine)
        nodes, K_local, M_local = _assemble_element_subset(elements, element_ids)
        boundary = [i for i, node in enumerate(nodes) if node in fine_boundary_nodes]
        bset = set(boundary)
        interior = [i for i in range(len(nodes)) if i not in bset]
        cb = _local_cb(K_local, M_local, interior, boundary, k_level2, fine)
        keys = [("level2_mode", fine, j) for j in range(cb["k_actual"])]
        keys.extend(("node", nodes[i]) for i in boundary)
        cb.update(
            {
                "fine": fine,
                "coarse": int(coarse_of_fine[fine]),
                "element_ids": element_ids,
                "nodes": nodes,
                "keys": keys,
            }
        )
        fine_data.append(cb)

    # Transformation 2: compatibility/assembly of level-2 reductions in each
    # level-1 subdomain.
    level1_assembled: List[dict] = []
    for coarse in range(n_level1):
        children = [item for item in fine_data if item["coarse"] == coarse]
        S2, keys2, counts2 = _compatibility([item["keys"] for item in children])
        Kdup = _block_diag([item["K_star"] for item in children])
        Mdup = _block_diag([item["M_star"] for item in children])
        K2 = S2.T @ Kdup @ S2
        M2 = S2.T @ Mdup @ S2
        level1_assembled.append(
            {
                "coarse": coarse,
                "children": children,
                "S": S2,
                "keys": keys2,
                "counts": counts2,
                "K": 0.5 * (K2 + K2.T),
                "M": 0.5 * (M2 + M2.T),
            }
        )

    # Transformation 3: cooperative Craig-Bampton reductions on level 1.
    coarse_data: List[dict] = []
    for assembled in level1_assembled:
        keys = assembled["keys"]
        boundary = [
            i for i, key in enumerate(keys) if key[0] == "node" and key[1] in coarse_boundary_nodes
        ]
        bset = set(boundary)
        interior = [i for i in range(len(keys)) if i not in bset]
        coarse = int(assembled["coarse"])
        cb = _local_cb(assembled["K"], assembled["M"], interior, boundary, k_level1, coarse)
        reduced_keys = [("level1_mode", coarse, j) for j in range(cb["k_actual"])]
        reduced_keys.extend(keys[i] for i in boundary)
        cb.update({"coarse": coarse, "assembled": assembled, "keys": reduced_keys})
        coarse_data.append(cb)

    # Transformation 4: leader-level compatibility/assembly and reduced solve.
    S1, final_keys, counts1 = _compatibility([item["keys"] for item in coarse_data])
    K1dup = _block_diag([item["K_star"] for item in coarse_data])
    M1dup = _block_diag([item["M_star"] for item in coarse_data])
    K_final = 0.5 * (S1.T @ K1dup @ S1 + (S1.T @ K1dup @ S1).T)
    M_final = 0.5 * (S1.T @ M1dup @ S1 + (S1.T @ M1dup @ S1).T)

    r_final = K_final.shape[0]
    if r_final < 1:
        raise ValueError("hierarchical reduction produced an empty pencil")
    wanted = r_final if num_modes is None else int(num_modes)
    if wanted < 1 or wanted > r_final:
        raise ValueError("num_modes must satisfy 1 <= num_modes <= reduced dimension")
    values, reduced_vectors = sla.eigh(K_final, M_final, subset_by_index=[0, wanted - 1])

    # Reverse the four transforms.  Equality of duplicated interface rows is a
    # direct test of S2/S1 and not an averaging assumption.
    T_global = np.zeros((K_global.shape[0], r_final))
    contributions: Dict[int, List[np.ndarray]] = {i: [] for i in range(K_global.shape[0])}
    coarse_offset = 0
    for coarse_cb in coarse_data:
        n_reduced = coarse_cb["K_star"].shape[0]
        y_coarse = S1[coarse_offset : coarse_offset + n_reduced, :]
        coarse_offset += n_reduced
        z_level1 = coarse_cb["T"] @ y_coarse
        assembled = coarse_cb["assembled"]
        p_level2 = assembled["S"] @ z_level1
        fine_offset = 0
        for fine_cb in assembled["children"]:
            n_fine_reduced = fine_cb["K_star"].shape[0]
            p_fine = p_level2[fine_offset : fine_offset + n_fine_reduced, :]
            fine_offset += n_fine_reduced
            local_physical = fine_cb["T"] @ p_fine
            for row, node in enumerate(fine_cb["nodes"]):
                contributions[int(node)].append(local_physical[row, :])

    max_duplicate_jump = 0.0
    for node, rows in contributions.items():
        if not rows:
            raise AssertionError(f"equation {node} was not reconstructed")
        reference_row = rows[0]
        for row in rows[1:]:
            max_duplicate_jump = max(
                max_duplicate_jump,
                float(np.linalg.norm(row - reference_row, ord=np.inf)),
            )
        T_global[node, :] = np.mean(rows, axis=0)

    vectors = T_global @ reduced_vectors
    return {
        "eigenvalues": values,
        "eigenvectors": vectors,
        "reduced_eigenvectors": reduced_vectors,
        "residuals": compute_residual(K_global, M_global, values, vectors),
        "K_star": K_final,
        "M_star": M_final,
        "T": T_global,
        "fine_subdomains": fine_data,
        "level1_assembled": level1_assembled,
        "coarse_subdomains": coarse_data,
        "S_level1": S1,
        "final_keys": final_keys,
        "compatibility_counts_level1": counts1,
        "max_duplicate_jump": max_duplicate_jump,
        "dimensions": {
            "original": K_global.shape[0],
            "after_level2_cb_before_compatibility": sum(x["K_star"].shape[0] for x in fine_data),
            "after_level2_compatibility": sum(x["K"].shape[0] for x in level1_assembled),
            "after_level1_cb_before_compatibility": sum(x["K_star"].shape[0] for x in coarse_data),
            "final": r_final,
        },
    }
