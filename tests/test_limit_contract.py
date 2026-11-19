"""WP-237a: the limit-analysis contract, executable part.

Schema validators for the version-0 JSON + NPZ formats, the parser and mesh builder for
``tests/limit_cases.txt``, convention checks, fixtures and the case evaluator. The mathematics is
decided in ``tests/limit_interface.md``. Nothing here assembles a conic program or certifies a bound.

Uses only the standard library, numpy and pytest. Regenerate the fixtures with
``python tests/test_limit_contract.py --write-fixtures``.
"""
import copy
import hashlib
import itertools
import json
import math
import re
import struct
import sys
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.zone_a

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "limit"
INTERFACE = HERE / "limit_interface.md"
CASES = HERE / "limit_cases.txt"
PROBE_OUTPUT = HERE.parent / "Ladruno_files" / "testbed" / "limit_analysis" / "p2_isochory_probe_output.txt"
PROBE_SCRIPT = PROBE_OUTPUT.with_name("p2_isochory_probe.py")
DEMO_SCRIPT = PROBE_OUTPUT.with_name("certificate_demo.py")
DEMO_RECORD = PROBE_OUTPUT.with_name("certificate_demo_record.json")
MESH_SCRIPT = PROBE_OUTPUT.with_name("make_gate_meshes.py")

FORMAT_VERSION = 0
FORMAT_PROBLEM = "ladruno-limit-problem"
FORMAT_RESULT = "ladruno-limit-result"
FORMAT_PROGRAM = "ladruno-limit-program"
FORMAT_FEASIBILITY = "ladruno-limit-feasibility"
FORMAT_CERTIFICATE = "ladruno-limit-certificate"
FORMAT_BRACKET = "ladruno-limit-bracket"
FORMAT_METHOD = "ladruno-limit-certificate-method"  # amendment 1, interface 15.9
FORMAT_REFUSAL = "ladruno-limit-refusal"  # amendment 1, interface 15.8

CONVENTIONS = {
    "sign": "tension_positive",
    "voigt": ["xx", "yy", "zz", "xy", "yz", "zx"],
    "stress_shear": "tensor",
    "strain_rate_shear": "engineering",
    "tet10": "vtk",
    "facets": "opposite_vertex_outward",
    "footing_resultant": "footing_on_soil",
}

UNIT_ROUNDOFF = 2.0 ** -53

TOLERANCES = {
    "TOL_GEOM_STRAIGHT": 1e-9,
    "TOL_GEOM_VOLUME": 1e-10,
    "TOL_GEOM_PLANE": 1e-9,
    "TOL_RANK": 1e-10,
    "TOL_EQ_ROUNDING_FACTOR": 64.0,
    "TOL_FLOOR_MAX": 1e-6,
    "TOL_KAPPA_CAP": 1e8,
    "TOL_PROJECTION_MOVE": 1e-6,
    "TOL_VALUE_MATCH": 1e-6,
    "TOL_CONE_ACCEPT": 1e-6,
    "TOL_CONE_MARGIN": 1e-12,
    "TOL_INTERIOR_MIN": 1e-6,
    "TOL_REPAIR_MAX": 1e-3,
    "TOL_FLOW_ACCEPT": 1e-6,
    "TOL_WORK_MIN": 1e-12,
    "TOL_FEAS_OK": 1e-9,
    "TOL_FEAS_LOCKED": 1e-6,
    "TOL_MONOTONE": 1e-6,
    "FLOOR_MIN_REL": 1e-12,
    "TOL_CASE_ONE_SIDED": 1e-6,
    "TOL_CASE_EXACT_GAP": 1e-6,
    "MAX_GAP_STRIP": 0.05,
    "MAX_GAP_PRANDTL": 0.10,
    "MAX_GAP_3D_VM": 0.30,
    "MAX_GAP_3D_DP": 0.50,
}

CRITERIA = {"von_mises", "drucker_prager"}
RESERVED_CRITERIA = {"tresca", "mohr_coulomb"}
BOUNDARY_KINDS = {"free", "fixed", "symmetry", "traction", "footing"}
RESERVED_BOUNDARY_KINDS = {"interface"}
INTERFACE_LAWS = {"bonded"}
RESERVED_INTERFACE_LAWS = {"no_tension", "coulomb"}
QUERY_FORMS = {"ray", "fixed_components", "support"}
LOWER_ID = "lower.p1d-vertex.v0"
CONTROL_ID = "upper.p2c-vertex.v0"
JUMP_ID = "upper.p2d-jump-bernstein.v0"
FORMULATIONS = {"lower": {LOWER_ID}, "upper": {CONTROL_ID, JUMP_ID}}
SOLVER_STATUS_WITH_VALUE = {"solved", "almost_solved"}
SOLVER_STATUS_NO_VALUE = {
    "primal_infeasible", "dual_infeasible", "time_limit", "iteration_limit", "numerical_error", "not_run",
}
DOMAIN_CLAIMS = {
    "box_stated", "unverified", "halfspace_by_zero_extension", "halfspace_by_extension_elements",
}
REASON_CODES = {
    "no_candidate", "hash_mismatch", "query_mismatch", "formulation_not_admitted",
    "formulation_not_applicable", "curved_element", "degenerate_element", "boundary_violation",
    "projection_move_too_large", "equality_residual_above_floor", "equality_ill_conditioned",
    "cone_violation_too_large", "no_interior_field", "interior_field_invalid", "repair_too_large",
    "flow_violation_too_large", "no_admissible_correction", "work_nonpositive", "value_mismatch",
    "unsupported_domain_claim", "literal_zero_violated",
}
CHECK_KEYS = {
    "rho_div", "rho_jump", "rho_bc", "rho_resultant", "projection_move", "omega_eq", "kappa_eq",
    "tau_eq", "rho_cone", "repair_mu", "straightness_max", "rho_flow", "rho_jump_cone",
    "work_reference", "rho_apex",  # rho_apex: amendment 2, M3
}
FEASIBILITY_STATUS = {"formulation_feasible", "formulation_infeasible", "inconclusive"}
SCHEMA_FIXTURE_CHECKER = "schema-fixture"
PROCEDURE_DEMO_CHECKER = "procedure-demo"
NON_GATE_CHECKERS = {SCHEMA_FIXTURE_CHECKER, PROCEDURE_DEMO_CHECKER}
# Amendment 1: method reports (15.9) and capacity refusals (15.8).
KAPPA_SVD_METHODS = {"dense_svd", "dense_svd_blocks"}
KAPPA_EIGEN_METHODS = {"eigen_normal_rows", "eigen_normal_columns"}
COMPLETENESS_ROUTES = {"full_spectrum", "inertia"}
DEPENDENT_ROUTES = {"shape_forced_zeros", "dense_svd_blocks", "null_space_basis"}
INTERIOR_SOURCES = {"none", "zero", "result"}
REFUSAL_KINDS = {"capacity", "conditioning"}
REFUSAL_STAGES = {"projection", "dependent_rows", "kappa_bound"}
APEX_KINDS = {0: "flow", 1: "jump"}  # amendment 2, interface 19.2.1: declared apex face Z
APEX_ROWS = {0: 6, 1: 3}  # exact rows added per condition: the strain rate d_ei, or the coefficient c_B
Z_FORCED = {"proved", "unproved", "absent"}
BUDGET_EXCEEDED = {"memory", "time"}  # amendment 2, m4

EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
FACET_VERTICES = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))
FACET_NODES = ((1, 2, 3, 5, 9, 8), (0, 3, 2, 7, 9, 6), (0, 1, 3, 4, 8, 7), (0, 2, 1, 6, 5, 4))
COMPONENTS = ("ux", "uy", "uz", "rx", "ry", "rz")
AXES = {"x": 0, "y": 1, "z": 2}

S2, S3, S6 = math.sqrt(2.0), math.sqrt(3.0), math.sqrt(6.0)
Q = np.array([
    [1 / S2, -1 / S2, 0, 0, 0, 0],
    [-1 / S6, -1 / S6, 2 / S6, 0, 0, 0],
    [0, 0, 0, S2, 0, 0],
    [0, 0, 0, 0, S2, 0],
    [0, 0, 0, 0, 0, S2],
])
ONE = np.array([1.0, 1, 1, 0, 0, 0])


class ContractError(ValueError):
    """A document or array set breaks the contract. ``code`` names the rule."""

    def __init__(self, code, detail=""):
        super().__init__(f"{code}: {detail}" if detail else code)
        self.code = code


def fail(code, detail=""):
    raise ContractError(code, detail)


def require(condition, code, detail=""):
    if not condition:
        fail(code, detail)


# ---------------------------------------------------------------------------------------------
# Canonical hash (interface 15.1)
# ---------------------------------------------------------------------------------------------

def _u64(value):
    return struct.pack("<Q", value)


def canonical_encode(obj):
    if obj is None:
        return b"n"
    if obj is True:
        return b"t"
    if obj is False:
        return b"f"
    if isinstance(obj, (int, float)):
        return b"d" + struct.pack("<d", float(obj) + 0.0)  # -0.0 hashes as 0.0
    if isinstance(obj, str):
        raw = obj.encode("utf-8")
        return b"s" + _u64(len(raw)) + raw
    if isinstance(obj, (list, tuple)):
        return b"l" + _u64(len(obj)) + b"".join(canonical_encode(item) for item in obj)
    if isinstance(obj, dict):
        keys = sorted(obj, key=lambda key: key.encode("utf-8"))
        body = b"".join(canonical_encode(key) + canonical_encode(obj[key]) for key in keys)
        return b"m" + _u64(len(keys)) + body
    raise TypeError(type(obj))


def document_hash(doc, arrays, hash_key):
    digest = hashlib.sha256()
    digest.update(canonical_encode({key: value for key, value in doc.items() if key != hash_key}))
    for name in sorted(arrays):
        array = np.asarray(arrays[name])
        if array.dtype.kind == "f":
            array = array + 0.0  # -0.0 hashes as 0.0
        digest.update(b"a" + canonical_encode(name) + canonical_encode(array.dtype.str))
        digest.update(_u64(array.ndim) + b"".join(_u64(dim) for dim in array.shape))
        digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def strict_json_loads(text):
    def refuse(token):
        fail("json_non_finite", token)

    def finite_float(token):
        value = float(token)
        require(math.isfinite(value), "json_non_finite", token)
        return value

    def finite_int(token):
        require(math.isfinite(float(token)), "json_non_finite", token)
        return int(token)
    return json.loads(text, parse_constant=refuse, parse_float=finite_float, parse_int=finite_int)


def load_document(path):
    path = Path(path)
    doc = strict_json_loads(path.read_text(encoding="utf-8"))
    arrays = {}
    if doc.get("arrays"):
        with np.load(path.parent / doc["arrays"], allow_pickle=False) as npz:
            arrays = {name: npz[name] for name in npz.files}
    return doc, arrays


# ---------------------------------------------------------------------------------------------
# Convention helpers (interface 2 and 3)
# ---------------------------------------------------------------------------------------------

def voigt_n(n):
    return np.array([
        [n[0], 0, 0, n[1], 0, n[2]],
        [0, n[1], 0, n[0], n[2], 0],
        [0, 0, n[2], 0, n[1], n[0]],
    ], dtype=float)


def tangent_basis(n):
    n = np.asarray(n, dtype=float)
    axis = int(np.argmin(np.abs(n)))
    e = np.zeros(3)
    e[axis] = 1.0
    t1 = np.cross(e, n)
    t1 /= np.linalg.norm(t1)
    return t1, np.cross(n, t1)


def symmetry_set_normal(problem, arrays, set_index):
    """n_S of interface 5, exactly as amendment 2 (m6, M5) states it: the first facet of the set in `facets` order,
    its vertices a, b, c in 3.2 order, u = x_b - x_a, w = x_c - x_a, q = u x w by components, and
    n_S = q / sqrt(q_x q_x + q_y q_y + q_z q_z), summed left to right, all in IEEE double. Every set-normal row
    builder here and in the demonstration calls it."""
    row = int(np.nonzero(arrays["facet_set"] == set_index)[0][0])
    e, f = (int(v) for v in arrays["facets"][row])
    xa, xb, xc = (arrays["nodes"][arrays["tets"][e][v]].tolist() for v in FACET_VERTICES[f])
    u = [xb[0] - xa[0], xb[1] - xa[1], xb[2] - xa[2]]
    w = [xc[0] - xa[0], xc[1] - xa[1], xc[2] - xa[2]]
    q = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
    norm = math.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2])
    return np.array([q[0] / norm, q[1] / norm, q[2] / norm])


def tet_gradients(x):
    jac = (x[1:] - x[0]).T
    grads = np.zeros((4, 3))
    grads[1:] = np.linalg.inv(jac)
    grads[0] = -grads[1:].sum(axis=0)
    return grads, np.linalg.det(jac) / 6.0


def facet_geometry(x, facet):
    a, b, c = x[list(FACET_VERTICES[facet])]
    cross = np.cross(b - a, c - a)
    area = 0.5 * np.linalg.norm(cross)
    return cross / (2.0 * area), area


def p2_gradients_at_vertex(grads, vertex):
    out = np.zeros((10, 3))
    for i in range(4):
        out[i] = (3.0 if i == vertex else -1.0) * grads[i]
    for index, (a, b) in enumerate(EDGES):
        if a == vertex:
            out[4 + index] = 4.0 * grads[b]
        elif b == vertex:
            out[4 + index] = 4.0 * grads[a]
    return out


def p2_values(bary):
    values = np.zeros(10)
    for i in range(4):
        values[i] = bary[i] * (2.0 * bary[i] - 1.0)
    for index, (a, b) in enumerate(EDGES):
        values[4 + index] = 4.0 * bary[a] * bary[b]
    return values


def strain_rate_at_vertex(x_tet10, velocity, vertex):
    grads, _ = tet_gradients(x_tet10[:4])
    dn = p2_gradients_at_vertex(grads, vertex)
    return sum(voigt_n(dn[a]).T @ velocity[a] for a in range(10))


def invariants(sigma):
    sigma = np.asarray(sigma, dtype=float)
    i1 = sigma[0] + sigma[1] + sigma[2]
    dev = sigma[:3] - i1 / 3.0
    j2 = 0.5 * float(dev @ dev) + float(sigma[3:] @ sigma[3:])
    return i1, j2


def yield_value(sigma, k, alpha=0.0):
    i1, j2 = invariants(sigma)
    return math.sqrt(j2) + alpha * i1 - k


def yield_from_principal(sigma, k, alpha=0.0):
    s = np.asarray(sigma, dtype=float)
    tensor = np.array([[s[0], s[3], s[5]], [s[3], s[1], s[4]], [s[5], s[4], s[2]]])
    p = np.linalg.eigvalsh(tensor)
    j2 = ((p[0] - p[1]) ** 2 + (p[1] - p[2]) ** 2 + (p[2] - p[0]) ** 2) / 6.0
    return math.sqrt(j2) + alpha * p.sum() - k


def strain_split(d):
    d = np.asarray(d, dtype=float)
    theta = d[0] + d[1] + d[2]
    dev = d[:3] - theta / 3.0
    return theta, math.sqrt(float(dev @ dev) + 0.5 * float(d[3:] @ d[3:]))


def support_function(d, k, alpha=0.0, tol=1e-12):
    theta, enorm = strain_split(d)
    scale = max(1.0, float(np.max(np.abs(d))))
    if alpha == 0.0:
        return S2 * k * enorm if abs(theta) <= tol * scale else math.inf
    if theta >= 3.0 * S2 * alpha * enorm - tol * scale:
        return k * theta / (3.0 * alpha)
    return math.inf


def tau_eq():
    """The equality acceptance level of interface 17.2: a rounding level, with no conditioning factor."""
    return TOLERANCES["TOL_EQ_ROUNDING_FACTOR"] * UNIT_ROUNDOFF


def row_backward_error(matrix, rhs, x, scale):
    """omega of interface 17.1: max_i |a_i x - g_i| / (sum_j |a_ij| scale_j + |g_i|)."""
    matrix, rhs, x, scale = (np.asarray(item, dtype=float) for item in (matrix, rhs, x, scale))
    residual = np.abs(matrix @ x - rhs)
    denominator = np.abs(matrix) @ scale + np.abs(rhs)
    require(not np.any((denominator == 0.0) & (residual != 0.0)), "equality_residual_above_floor")
    active = denominator > 0.0
    return float(np.max(residual[active] / denominator[active], initial=0.0))


def jump_constants(alpha):
    """tau (dilation per unit slip) and root = sqrt(1 - 12 alpha^2) of interface 12.2; k_J = k / root."""
    require(alpha * alpha * 12.0 < 1.0, "formulation_not_applicable", "alpha >= 1/sqrt(12)")
    root = math.sqrt(1.0 - 12.0 * alpha * alpha)
    return 3.0 * alpha / root, root


def bernstein_from_nodal(nodal):
    nodal = np.asarray(nodal, dtype=float)
    coefficients = [nodal[0], nodal[1], nodal[2]]
    for q, (a, b) in enumerate(((0, 1), (1, 2), (2, 0))):
        coefficients.append(2.0 * nodal[3 + q] - 0.5 * (nodal[a] + nodal[b]))
    return np.array(coefficients)


def plane_strain_match(c, phi_deg):
    t = math.tan(math.radians(phi_deg))
    root = math.sqrt(9.0 + 12.0 * t * t)
    return 3.0 * c / root, t / root


def prandtl_nq(phi_deg):
    phi = math.radians(phi_deg)
    return math.exp(math.pi * math.tan(phi)) * math.tan(math.pi / 4 + phi / 2) ** 2


def prandtl_nc(phi_deg):
    if phi_deg == 0:
        return 2.0 + math.pi
    return (prandtl_nq(phi_deg) - 1.0) / math.tan(math.radians(phi_deg))


def prandtl_extent(b, phi_deg):
    phi = math.radians(phi_deg)
    r1 = b / math.cos(math.pi / 4 + phi / 2) * math.exp(math.pi / 2 * math.tan(phi))
    return b + 2.0 * r1 * math.cos(math.pi / 4 - phi / 2), r1


TARGET_FORMULAS = {
    "vm_uniaxial": lambda p: S3 * p["k"],
    "vm_uniaxial_dead": lambda p: (1.0 - p["dead_fraction"]) * S3 * p["k"],
    "dp_uniaxial_tension": lambda p: p["k"] / (1.0 / S3 + p["alpha"]),
    "dp_uniaxial_compression": lambda p: p["k"] / (1.0 / S3 - p["alpha"]),
    "dp_hydrostatic": lambda p: p["k"] / (3.0 * p["alpha"]),
    "k_direct": lambda p: p["k"],
    "shear_block": lambda p: p["k"] / (p["gamma"] * p["H"]),
    "dead_shear": lambda p: p["k"] - p["gamma"] * p["H"],
    "dp_apex_column": lambda p: p["k"] / (3.0 * p["alpha"]) - p["gamma"] * p["H"],
    "column_static": lambda p: S3 * p["k"] - p["gamma"] * p["H"],
    "column_uniform": lambda p: S3 * p["k"] - 0.5 * p["gamma"] * p["H"],
    "prandtl_tresca": lambda p: (2.0 + math.pi) * p["cu"],
    "prandtl_nc": lambda p: p["c"] * prandtl_nc(p["phi_deg"]),
    "prandtl_nq": lambda p: p["q0"] * prandtl_nq(p["phi_deg"]),
}


def band_formula(name, params, first_layer):
    """Upper fields with one plastic layer of thickness h: interface 16.3 and 16.5."""
    h = first_layer
    if name == "shear_band":
        return params["k"] / (params["gamma"] * (params["H"] - 0.5 * h))
    if name == "dead_shear_band":
        return params["k"] - params["gamma"] * (params["H"] - 0.5 * h)
    if name == "apex_band":
        return params["k"] / (3.0 * params["alpha"]) - params["gamma"] * (params["H"] - 0.5 * h)
    fail("case_format", name)


# ---------------------------------------------------------------------------------------------
# Meshes (interface 3.4, 3.5)
# ---------------------------------------------------------------------------------------------

def tet10_from_tet4(vertices, tets):
    edges = sorted({tuple(sorted((int(t[a]), int(t[b])))) for t in tets for a, b in EDGES})
    index = {edge: len(vertices) + i for i, edge in enumerate(edges)}
    mids = np.array([0.5 * (vertices[a] + vertices[b]) for a, b in edges]).reshape(-1, 3)
    rows = [
        [int(v) for v in t] + [index[tuple(sorted((int(t[a]), int(t[b]))))] for a, b in EDGES]
        for t in tets
    ]
    return np.vstack([vertices, mids]), np.array(rows, dtype=np.int64)


def kuhn_tet4(lo, hi, counts):
    lo, hi = np.asarray(lo, dtype=float), np.asarray(hi, dtype=float)
    nx, ny, nz = counts

    def vid(i, j, k):
        return i + (nx + 1) * (j + (ny + 1) * k)

    vertices = np.zeros(((nx + 1) * (ny + 1) * (nz + 1), 3))
    for k, j, i in itertools.product(range(nz + 1), range(ny + 1), range(nx + 1)):
        vertices[vid(i, j, k)] = lo + (hi - lo) * np.array([i / nx, j / ny, k / nz])
    tets = []
    for k, j, i in itertools.product(range(nz), range(ny), range(nx)):
        for permutation in itertools.permutations(range(3)):
            corner = [i, j, k]
            tet = [vid(*corner)]
            for axis in permutation:
                corner = list(corner)
                corner[axis] += 1
                tet.append(vid(*corner))
            if np.linalg.det(vertices[tet[1:]] - vertices[tet[0]]) < 0.0:
                tet[2], tet[3] = tet[3], tet[2]
            tets.append(tet)
    return vertices, np.array(tets, dtype=np.int64)


def refine_tet4(vertices, tets):
    """One uniform refinement into eight children (interface 3.6)."""
    pairs = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    edges = sorted({tuple(sorted((int(t[a]), int(t[b])))) for t in tets for a, b in pairs})
    index = {edge: len(vertices) + i for i, edge in enumerate(edges)}
    points = np.vstack([vertices, [0.5 * (vertices[a] + vertices[b]) for a, b in edges]])
    children = []
    for t in tets:
        v = [int(item) for item in t]
        mid = {pair: index[tuple(sorted((v[pair[0]], v[pair[1]])))] for pair in pairs}
        m01, m02, m03, m12, m13, m23 = (mid[pair] for pair in pairs)
        for child in (
            [v[0], m01, m02, m03], [m01, v[1], m12, m13], [m02, m12, v[2], m23], [m03, m13, m23, v[3]],
            [m01, m02, m03, m13], [m01, m02, m12, m13], [m02, m03, m13, m23], [m02, m12, m13, m23],
        ):
            if np.linalg.det(points[child[1:]] - points[child[0]]) < 0.0:
                child[2], child[3] = child[3], child[2]
            children.append(child)
    return points, np.array(children, dtype=np.int64)


def case_tet4(case, mesh_index):
    entry = case["mesh"][mesh_index]
    if entry[0] == "kuhn":
        return kuhn_tet4(case["box"][:3], case["box"][3:], entry[1])
    with np.load(DATA / entry[1], allow_pickle=False) as npz:
        vertices, tets = npz["nodes"], npz["tets"]
    for _ in range(entry[2]):
        vertices, tets = refine_tet4(vertices, tets)
    return vertices, tets


def mesh_facets(tets):
    owners = {}
    for e, tet in enumerate(tets):
        for f, local in enumerate(FACET_VERTICES):
            owners.setdefault(tuple(sorted(int(tet[a]) for a in local)), []).append((e, f))
    interior, exterior = [], []
    for entries in owners.values():
        require(len(entries) <= 2, "mesh_nonconforming", "facet with more than two owners")
        (interior if len(entries) == 2 else exterior).append(tuple(sorted(entries)) if len(entries) == 2
                                                             else entries[0])
    return sorted(interior), sorted(exterior)


# ---------------------------------------------------------------------------------------------
# Case file (interface 16 and the header of limit_cases.txt)
# ---------------------------------------------------------------------------------------------

def parse_key_values(tokens):
    return {key: float(value) for key, value in (token.split("=", 1) for token in tokens)}


def parse_cases(path=CASES):
    return parse_case_text(Path(path).read_text(encoding="ascii"))


def parse_case_text(text):
    cases, current = {}, None
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        key, _, rest = line.partition(" ")
        rest = rest.strip()
        if key == "case":
            require(current is None and rest not in cases, "case_format", line)
            current = {"name": rest, "lines": []}
            cases[rest] = current
        elif key == "end":
            require(current is not None, "case_format", line)
            current = None
        else:
            require(current is not None, "case_format", line)
            current["lines"].append((key, rest))
    require(current is None, "case_format", "missing end")
    return {name: _digest_case(case) for name, case in cases.items()}


def _vector(tokens, count):
    values = [float(token) for token in tokens[:count]]
    require(len(values) == count, "case_format", str(tokens))
    return values


def _digest_region(tokens):
    region = {"name": tokens[0], "criterion": tokens[1], "b0": [0.0] * 3, "b1": [0.0] * 3}
    i = 2
    while i < len(tokens):
        key = tokens[i]
        if key in ("k", "alpha"):
            region[key] = float(tokens[i + 1])
            i += 2
        elif key in ("b0", "b1"):
            region[key] = _vector(tokens[i + 1:], 3)
            i += 4
        else:
            fail("case_format", key)
    return region


def _digest_boundary(tokens):
    item = {"name": tokens[0], "kind": tokens[1], "ranges": []}
    i = 2
    while i < len(tokens):
        key = tokens[i]
        if key == "plane":
            item["plane"] = (AXES[tokens[i + 1]], float(tokens[i + 2]))
            i += 3
        elif key == "range":
            item["ranges"].append((AXES[tokens[i + 1]], float(tokens[i + 2]), float(tokens[i + 3])))
            i += 4
        elif key in ("t0", "t1"):
            c = _vector(tokens[i + 1:], 3)
            item[key] = (np.array(c), np.zeros((3, 3)))
            i += 4
        elif key in ("t0_affine", "t1_affine"):
            v = _vector(tokens[i + 1:], 12)
            item[key[:2]] = (np.array(v[:3]), np.array(v[3:]).reshape(3, 3))
            i += 13
        elif key == "footing":
            item["footing"] = tokens[i + 1]
            i += 2
        else:
            fail("case_format", key)
    require("plane" in item, "case_format", "boundary without plane")
    return item


def _digest_case(raw):
    case = {
        "name": raw["name"], "source": [], "checked": [], "mesh": [], "formulation": [], "control": [],
        "regions": [], "footings": [], "boundaries": [], "queries": [], "F0": {}, "F1": {}, "restrained": [],
        "assert": [], "gate": [], "record": [], "nested": [],
    }
    for key, rest in raw["lines"]:
        tokens = rest.split()
        if key in ("status", "class", "desc", "units", "gravity_sensitive"):
            case[key] = rest
        elif key in ("source", "checked", "formulation", "control"):
            case[key].append(rest)
        elif key == "box":
            case["box"] = _vector(tokens, 6)
        elif key == "mesh":
            case["mesh"].append(("kuhn", tuple(int(token) for token in tokens)))
        elif key == "mesh_file":
            refine = 0
            if len(tokens) > 1:
                require(len(tokens) == 3 and tokens[1] == "refine", "case_format", rest)
                refine = int(tokens[2])
            case["mesh"].append(("file", tokens[0], refine))
        elif key == "nested":
            case["nested"].append((int(tokens[0]), int(tokens[1])))
        elif key == "phi_deg":
            case["phi_deg"] = float(rest)
        elif key == "region":
            case["regions"].append(_digest_region(tokens))
        elif key == "footing":
            require(tokens[1] == "xref" and tokens[5] == "law", "case_format", rest)
            case["footings"].append({"name": tokens[0], "x_ref": _vector(tokens[2:], 3), "interface_law": tokens[6]})
        elif key == "boundary":
            case["boundaries"].append(_digest_boundary(tokens))
        elif key == "query":
            case["queries"].append({"id": tokens[0], "form": tokens[1]})
        elif key in ("F0", "F1"):
            case[key][tokens[0]] = _vector(tokens[1:], 6)
        elif key == "restrained":
            case["restrained"].extend((tokens[0], name) for name in tokens[1:])
        elif key == "target":
            if tokens[0] == "none":
                case["target"] = None
            else:
                require(tokens[1] == "formula", "case_format", rest)
                case["target"] = float(tokens[0])
                case["formula"] = tokens[2]
                case["formula_params"] = parse_key_values(tokens[3:])
        elif key == "probe":
            require(tokens[0] == "rigid_motion" and tokens[1] == "u" and tokens[5] == "w" and tokens[9] == "about",
                    "case_format", rest)
            case["probe"] = {"u": _vector(tokens[2:], 3), "w": _vector(tokens[6:], 3), "about": _vector(tokens[10:], 3)}
        elif key == "assert":
            case["assert"].append(tokens)
        elif key == "gate":
            case["gate"].append(tokens)
        elif key == "record":
            case["record"].append(tokens)
        else:
            fail("case_format", f"unknown key {key}")
    return case


def case_stem(name, mesh_index=0):
    """File stem of one case on one of its meshes (interface 15.2)."""
    return f"{name}.m{mesh_index}"


def record_stem(name, formulation, mesh_index=0, query="q"):
    return f"{case_stem(name, mesh_index)}.{query}.{formulation}"


def bracket_stem(name, lower_formulation, upper_formulation, mesh_index=0, query="q"):
    return f"{case_stem(name, mesh_index)}.{query}.{lower_formulation}__{upper_formulation}"


def build_case_problem(case, mesh_index=0):
    """Build the problem document and arrays of one case on one of its meshes."""
    mesh_entry = case["mesh"][mesh_index]
    lo, hi = np.array(case["box"][:3]), np.array(case["box"][3:])
    vertices, tet4 = case_tet4(case, mesh_index)
    nodes, tets = tet10_from_tet4(vertices, tet4)
    _, exterior = mesh_facets(tets)
    tol = 1e-10 * float(np.linalg.norm(hi - lo))
    footing_index = {item["name"]: i for i, item in enumerate(case["footings"])}
    facets = np.array(exterior, dtype=np.int64)
    facet_set = np.zeros(len(facets), dtype=np.int64)
    t0, t1 = np.zeros((len(facets), 3, 3)), np.zeros((len(facets), 3, 3))
    for row, (e, f) in enumerate(exterior):
        x = nodes[tets[e][list(FACET_VERTICES[f])]]
        for index, item in enumerate(case["boundaries"]):
            axis, coordinate = item["plane"]
            if np.any(np.abs(x[:, axis] - coordinate) > tol):
                continue
            if any(np.any(x[:, a] < low - tol) or np.any(x[:, a] > high + tol) for a, low, high in item["ranges"]):
                continue
            facet_set[row] = index
            for name, target in (("t0", t0), ("t1", t1)):
                if name in item:
                    c, g = item[name]
                    target[row] = c + x @ g.T
            break
        else:
            fail("case_unmatched_facet", f"{case['name']}: facet {(e, f)}")
    sets = []
    for item in case["boundaries"]:
        entry = {"name": item["name"], "kind": item["kind"]}
        if "footing" in item:
            entry["footing"] = footing_index[item["footing"]]
        sets.append(entry)
    queries = []
    n_f = len(case["footings"])
    names = [item["name"] for item in case["footings"]]
    for query in case["queries"]:
        queries.append({
            "id": query["id"], "form": query["form"],
            "F0": [case["F0"].get(name, [0.0] * 6) for name in names],
            "F1": [case["F1"].get(name, [0.0] * 6) for name in names],
            "restrained": [[footing_index[f], COMPONENTS.index(c)] for f, c in case["restrained"]],
        })
    units = dict(part.strip().split(" ", 1) for part in case["units"].split(";")[:2])
    if mesh_entry[0] == "kuhn":
        mesh_label = "x".join(map(str, mesh_entry[1]))
    else:
        mesh_label = f"{mesh_entry[1].split('.')[0]}+{mesh_entry[2]}"
    doc = {
        "format": FORMAT_PROBLEM, "format_version": FORMAT_VERSION,
        "arrays": f"{case_stem(case['name'], mesh_index)}.problem.npz", "name": f"{case['name']}:{mesh_label}",
        "units": {"length": units["length"], "force": units["force"]},
        "conventions": copy.deepcopy(CONVENTIONS),
        "domain": {"target": "box", "truncation_sets": [], "symmetry_sets": [], "ground": None,
                   "exterior_load_free": True},
        "length_scale": float(np.max(hi - lo)),
        "regions": copy.deepcopy(case["regions"]), "boundary_sets": sets,
        "footings": copy.deepcopy(case["footings"]), "queries": queries,
    }
    arrays = {
        "nodes": nodes.astype("<f8"), "tets": tets.astype("<i8"),
        "element_region": np.zeros(len(tets), dtype="<i8"), "facets": facets.astype("<i8"),
        "facet_set": facet_set.astype("<i8"), "t0": t0.astype("<f8"), "t1": t1.astype("<f8"),
    }
    assert n_f == len(names)
    doc["problem_sha256"] = document_hash(doc, arrays, "problem_sha256")
    return doc, arrays


# ---------------------------------------------------------------------------------------------
# Validators
# ---------------------------------------------------------------------------------------------

def _keys(doc, required, optional=(), code="unknown_key"):
    keys = set(doc)
    missing = set(required) - keys
    require(not missing, "missing_key", str(sorted(missing)))
    extra = keys - set(required) - set(optional)
    require(not extra, code, str(sorted(extra)))


def _finite_number(value, code):
    require(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value), code, repr(value))


def _vector3(value, code):
    require(isinstance(value, list) and len(value) == 3, code, repr(value))
    for item in value:
        _finite_number(item, code)


def _array(arrays, name, dtype, shape):
    require(name in arrays, "missing_array", name)
    array = arrays[name]
    require(array.dtype.str == dtype, "array_dtype", f"{name}: {array.dtype.str}")
    require(array.shape == tuple(shape), "array_shape", f"{name}: {array.shape} != {tuple(shape)}")
    if dtype == "<f8":
        require(np.all(np.isfinite(array)), "array_non_finite", name)
    return array


def _load_table(query, key, n_f):
    table = query.get(key)
    require(isinstance(table, list) and len(table) == n_f, "query_shape", key)
    for row in table:
        require(isinstance(row, list) and len(row) == 6, "query_shape", key)
        for value in row:
            _finite_number(value, "query_shape")
    return np.array(table, dtype=float).reshape(n_f, 6)


def _validate_query(query, doc, arrays):
    n_f = len(doc["footings"])
    form = query.get("form")
    require(form in QUERY_FORMS, "query_form", repr(form))
    require(isinstance(query.get("id"), str) and query["id"], "query_id")
    restrained = query.get("restrained")
    require(isinstance(restrained, list), "query_restrained")
    mask = np.zeros((n_f, 6), dtype=bool)
    for item in restrained:
        require(isinstance(item, list) and len(item) == 2, "query_restrained", repr(item))
        f, c = item
        require(isinstance(f, int) and 0 <= f < n_f and isinstance(c, int) and 0 <= c < 6,
                "query_restrained", repr(item))
        require(not mask[f, c], "query_restrained", "duplicate")
        mask[f, c] = True
    has_t1 = bool(np.any(arrays["t1"] != 0.0))
    has_b1 = any(any(value != 0.0 for value in region["b1"]) for region in doc["regions"])
    if form == "ray":
        _keys(query, {"id", "form", "restrained", "F0", "F1"})
        f0, f1 = _load_table(query, "F0", n_f), _load_table(query, "F1", n_f)
        require(not np.any(f0[mask]) and not np.any(f1[mask]), "query_restrained_load")
        require(np.any(f1) or has_t1 or has_b1, "query_no_reference_load")
    elif form == "fixed_components":
        _keys(query, {"id", "form", "restrained", "fixed", "direction"})
        direction = _load_table(query, "direction", n_f)
        require(not np.any(direction[mask]), "query_restrained_load")
        require(isinstance(query["fixed"], list), "query_fixed")
        for item in query["fixed"]:
            require(isinstance(item, list) and len(item) == 3, "query_fixed", repr(item))
            f, c, value = item
            require(isinstance(f, int) and 0 <= f < n_f and isinstance(c, int) and 0 <= c < 6, "query_fixed")
            _finite_number(value, "query_fixed")
            require(direction[f, c] == 0.0, "query_fixed_direction")
            require(not mask[f, c] or value == 0.0, "query_restrained_load")
        require(np.any(direction) or has_t1 or has_b1, "query_no_reference_load")
    else:
        _keys(query, {"id", "form", "restrained", "direction"})
        require(n_f >= 1, "query_support_needs_footing")
        direction = _load_table(query, "direction", n_f)
        require(not np.any(direction[mask]), "query_restrained_load")
        require(np.any(direction), "query_zero_direction")
        require(not has_t1 and not has_b1, "query_support_reference_load")


def validate_problem(doc, arrays):
    require(doc.get("format") == FORMAT_PROBLEM, "format", repr(doc.get("format")))
    require(doc.get("format_version") == FORMAT_VERSION, "format_version")
    _keys(doc, {
        "format", "format_version", "arrays", "name", "units", "conventions", "domain", "length_scale",
        "regions", "boundary_sets", "footings", "queries", "problem_sha256",
    })
    require(doc["conventions"] == CONVENTIONS, "conventions")
    require(set(doc["units"]) == {"length", "force"}, "units")
    _finite_number(doc["length_scale"], "length_scale")
    require(doc["length_scale"] > 0.0, "length_scale")

    regions = doc["regions"]
    require(isinstance(regions, list) and regions, "regions")
    for region in regions:
        criterion = region.get("criterion")
        require(criterion not in RESERVED_CRITERIA, "criterion_reserved", repr(criterion))
        require(criterion in CRITERIA, "criterion", repr(criterion))
        k = region.get("k")
        _finite_number(k, "strength")
        _vector3(region.get("b0"), "body_force")
        _vector3(region.get("b1"), "body_force")
        if criterion == "von_mises":
            _keys(region, {"name", "criterion", "k", "b0", "b1"})
            require(k > 0.0, "strength", "von_mises needs k > 0")
        else:
            _keys(region, {"name", "criterion", "k", "alpha", "b0", "b1"})
            alpha = region["alpha"]
            _finite_number(alpha, "strength")
            require(0.0 < alpha < 1.0 / S3, "strength", "alpha out of range")
            require(k >= 0.0, "strength", "k < 0")

    footings = doc["footings"]
    require(isinstance(footings, list), "footings")
    for footing in footings:
        _keys(footing, {"name", "x_ref", "interface_law"})
        _vector3(footing["x_ref"], "footing_x_ref")
        law = footing["interface_law"]
        require(law not in RESERVED_INTERFACE_LAWS, "interface_law_reserved", repr(law))
        require(law in INTERFACE_LAWS, "interface_law", repr(law))

    sets = doc["boundary_sets"]
    require(isinstance(sets, list) and sets, "boundary_sets")
    names = [item.get("name") for item in sets]
    require(len(set(names)) == len(names), "boundary_set_name", "duplicate")
    footing_seen = set()
    for item in sets:
        kind = item.get("kind")
        require(kind not in RESERVED_BOUNDARY_KINDS, "boundary_kind_reserved", repr(kind))
        require(kind in BOUNDARY_KINDS, "boundary_kind", repr(kind))
        if kind == "footing":
            _keys(item, {"name", "kind", "footing"})
            index = item["footing"]
            require(isinstance(index, int) and 0 <= index < len(footings), "boundary_footing", repr(index))
            footing_seen.add(index)
        else:
            _keys(item, {"name", "kind"})
    require(footing_seen == set(range(len(footings))), "footing_without_set")

    domain = doc["domain"]
    _keys(domain, {"target", "truncation_sets", "symmetry_sets", "ground", "exterior_load_free"})
    require(domain["target"] in ("box", "halfspace"), "domain_target", repr(domain["target"]))
    require(isinstance(domain["exterior_load_free"], bool), "domain_target")
    require(set(domain["truncation_sets"]) <= set(names), "domain_truncation", "unknown set")
    require(set(domain["symmetry_sets"]) <= set(names), "domain_symmetry", "unknown set")
    kind_of = {item["name"]: item["kind"] for item in sets}
    require(all(kind_of[name] == "symmetry" for name in domain["symmetry_sets"]), "domain_symmetry", "not a symmetry set")
    if domain["target"] == "box":
        require(not domain["truncation_sets"] and not domain["symmetry_sets"] and domain["ground"] is None,
                "domain_truncation", "a box has no truncation sets, symmetry sets or ground plane")
    else:
        require(domain["truncation_sets"], "domain_truncation", "halfspace needs truncation sets")
        ground = domain["ground"]
        require(isinstance(ground, dict), "domain_ground", "halfspace needs a ground plane")
        _keys(ground, {"axis", "coordinate"})
        require(ground["axis"] in AXES, "domain_ground", repr(ground["axis"]))
        _finite_number(ground["coordinate"], "domain_ground")

    _validate_problem_arrays(doc, arrays)

    queries = doc["queries"]
    require(isinstance(queries, list) and queries, "queries")
    ids = [query.get("id") for query in queries]
    require(len(set(ids)) == len(ids), "query_id", "duplicate")
    for query in queries:
        _validate_query(query, doc, arrays)
    require(doc["problem_sha256"] == document_hash(doc, arrays, "problem_sha256"), "hash_mismatch")


def _validate_problem_arrays(doc, arrays):
    expected = {"nodes", "tets", "element_region", "facets", "facet_set", "t0", "t1"}
    require(set(arrays) == expected, "array_names", str(sorted(set(arrays) ^ expected)))
    require(arrays["nodes"].ndim == 2 and arrays["tets"].ndim == 2, "array_shape")
    n_nodes, n_tets = len(arrays["nodes"]), len(arrays["tets"])
    nodes = _array(arrays, "nodes", "<f8", (n_nodes, 3))
    tets = _array(arrays, "tets", "<i8", (n_tets, 10))
    regions = _array(arrays, "element_region", "<i8", (n_tets,))
    require(n_tets > 0, "array_shape", "no elements")
    require(tets.min() >= 0 and tets.max() < n_nodes, "tet_index")
    require(len(np.unique(tets)) == n_nodes, "orphan_node")
    require(regions.min() >= 0 and regions.max() < len(doc["regions"]), "element_region")

    edge_mid = {}
    for e, tet in enumerate(tets):
        require(len(set(tet.tolist())) == 10, "tet_index", f"repeated node in element {e}")
        x = nodes[tet]
        h = max(np.linalg.norm(x[a] - x[b]) for a, b in itertools.combinations(range(4), 2))
        volume = np.linalg.det(x[1:4] - x[0]) / 6.0
        require(volume > 0.0, "element_orientation", f"element {e}")
        require(volume >= TOLERANCES["TOL_GEOM_VOLUME"] * h ** 3, "element_degenerate", f"element {e}")
        for index, (a, b) in enumerate(EDGES):
            offset = np.linalg.norm(x[4 + index] - 0.5 * (x[a] + x[b]))
            require(offset <= TOLERANCES["TOL_GEOM_STRAIGHT"] * h, "element_curved", f"element {e}")
            key = tuple(sorted((int(tet[a]), int(tet[b]))))
            require(edge_mid.setdefault(key, int(tet[4 + index])) == int(tet[4 + index]),
                    "mesh_nonconforming", f"edge {key}")

    _, exterior = mesh_facets(tets)
    require(arrays["facets"].ndim == 2, "array_shape")
    n_facets = len(arrays["facets"])
    facets = _array(arrays, "facets", "<i8", (n_facets, 2))
    facet_set = _array(arrays, "facet_set", "<i8", (n_facets,))
    t0 = _array(arrays, "t0", "<f8", (n_facets, 3, 3))
    t1 = _array(arrays, "t1", "<f8", (n_facets, 3, 3))
    listed = [tuple(int(v) for v in row) for row in facets]
    require(listed == sorted(listed), "facets_order")
    require(listed == exterior, "facets_mismatch")
    require(facet_set.min() >= 0 and facet_set.max() < len(doc["boundary_sets"]), "facet_set")
    require(set(facet_set.tolist()) == set(range(len(doc["boundary_sets"]))), "boundary_set_empty")
    kinds = np.array([item["kind"] for item in doc["boundary_sets"]])[facet_set]
    other = kinds != "traction"
    require(not np.any(t0[other]) and not np.any(t1[other]), "traction_on_other_kind")

    tol = TOLERANCES["TOL_GEOM_PLANE"]
    span = float(np.linalg.norm(nodes.max(axis=0) - nodes.min(axis=0)))
    for index, item in enumerate(doc["boundary_sets"]):
        if item["kind"] != "symmetry":
            continue
        rows = np.nonzero(facet_set == index)[0]
        reference, origin = None, None
        for row in rows:
            e, f = facets[row]
            x = nodes[tets[e][:4]]
            normal, _ = facet_geometry(x, int(f))
            point = x[FACET_VERTICES[f][0]]
            if reference is None:
                reference, origin = normal, point
            require(np.linalg.norm(normal - reference) <= tol, "symmetry_not_planar", item["name"])
            require(abs(reference @ (point - origin)) <= tol * span, "symmetry_not_planar", item["name"])


def query_by_id(problem, query_id):
    for query in problem["queries"]:
        if query["id"] == query_id:
            return query
    fail("query_mismatch", repr(query_id))


def validate_formulation(bound, formulation, problem=None):
    if formulation not in FORMULATIONS.get(bound, ()):
        prefix = f"{bound}."
        known = isinstance(formulation, str) and formulation.startswith(prefix)
        fail("formulation_not_admitted" if known else "formulation_unknown", repr(formulation))
    if formulation == JUMP_ID and problem is not None:
        for region in problem["regions"]:
            jump_constants(region.get("alpha", 0.0))


def validate_result(doc, arrays, problem, problem_arrays):
    require(doc.get("format") == FORMAT_RESULT, "format", repr(doc.get("format")))
    require(doc.get("format_version") == FORMAT_VERSION, "format_version")
    optional = {"interior_value"}
    _keys(doc, {
        "format", "format_version", "arrays", "problem_sha256", "query_id", "bound", "formulation",
        "producer", "solver", "solver_status", "objective", "value", "wall_time_s", "result_sha256",
    }, optional)
    require(doc["problem_sha256"] == problem["problem_sha256"], "hash_mismatch", "problem")
    query_by_id(problem, doc["query_id"])
    bound = doc["bound"]
    require(bound in ("lower", "upper"), "bound", repr(bound))
    validate_formulation(bound, doc["formulation"], problem)
    _keys(doc["producer"], {"name", "version"})
    _keys(doc["solver"], {"name", "version", "settings"})
    status = doc["solver_status"]
    require(status in SOLVER_STATUS_WITH_VALUE | SOLVER_STATUS_NO_VALUE, "solver_status", repr(status))

    n_e, n_n = len(problem_arrays["tets"]), len(problem_arrays["nodes"])
    n_f = len(problem["footings"])
    broken = doc["formulation"] == JUMP_ID
    velocity_name = "element_velocity" if broken else "velocity"
    velocity_shape = (n_e, 10, 3) if broken else (n_n, 3)
    other_velocity = "velocity" if broken else "element_velocity"
    wrong_velocity = {other_velocity, f"interior_{other_velocity}"}
    shapes = {
        "lower": ({"sigma": (n_e, 4, 6), "resultant": (n_f, 6)},
                  {"interior_sigma": (n_e, 4, 6), "equilibrium_multiplier": (n_e, 3)}),
        "upper": ({velocity_name: velocity_shape, "footing_velocity": (n_f, 6)},
                  {f"interior_{velocity_name}": velocity_shape, "interior_footing_velocity": (n_f, 6),
                   "strain_rate": (n_e, 4, 6), "dissipation": (n_e, 4)}),
    }
    required, allowed = shapes[bound]
    # Amendment 2 (interface 19.2.1): the declared apex face is an Upper declaration; a Lower never carries one.
    face = {"apex_zero", "apex_zero_dual"} & set(arrays)
    require(not face or bound == "upper", "apex_zero_lower", str(sorted(face)))
    other = "upper" if bound == "lower" else "lower"
    foreign = (set(shapes[other][0]) | set(shapes[other][1]) | wrong_velocity) & set(arrays)
    require(not foreign, "result_wrong_bound_arrays", str(sorted(foreign)))
    unknown = set(arrays) - set(required) - set(allowed) - face
    require(not unknown, "array_names", str(sorted(unknown)))
    if status in SOLVER_STATUS_WITH_VALUE:
        _finite_number(doc["value"], "result_value")
        _finite_number(doc["objective"], "result_value")
        if bound == "lower":
            # Amendment 2, m5: the one schema rule on `objective` (interface 15.3, item 2 of amendment 1).
            require(doc["value"] == -doc["objective"], "result_objective", "Lower value = objective_sign * objective")
        for name, shape in required.items():
            _array(arrays, name, "<f8", shape)
    else:
        require(doc["value"] is None and doc["objective"] is None, "result_value_without_solution")
        require(not arrays, "result_arrays_without_solution")
    for name, shape in allowed.items():
        if name in arrays:
            _array(arrays, name, "<f8", shape)
    require(("interior_sigma" in arrays) == ("interior_value" in doc), "result_interior")
    require((f"interior_{velocity_name}" in arrays) == ("interior_footing_velocity" in arrays), "result_interior")
    validate_apex_zero(doc, arrays, problem, problem_arrays)
    require(doc["result_sha256"] == document_hash(doc, arrays, "result_sha256"), "hash_mismatch", "result")


def apex_zero_rows(arrays):
    """The declared apex face Z as sorted (kind, index) pairs, or [] (amendment 2, interface 19.2.1)."""
    if "apex_zero" not in arrays:
        return []
    return [(int(kind), int(index)) for kind, index in arrays["apex_zero"]]


def validate_apex_zero(doc, arrays, problem, problem_arrays):
    """Amendment 2 (interface 19.2.1): Z is an `<i8` array (n_Z, 2) of (kind, index), sorted and unique, n_Z >= 1.
    kind 0 is the flow condition at (e, i), index 4 e + i, alpha_e > 0; kind 1 is the jump condition at (F, B),
    index 6 F + B with F in the interior-facet order of 3.4, jump formulation only, tau > 0 from the region of e-.
    The optional producer dual `apex_zero_dual` (n_Z, 6) is a diagnostic; a jump row uses its first three entries."""
    if "apex_zero" not in arrays:
        require("apex_zero_dual" not in arrays, "apex_zero_invalid", "a dual without a declared face")
        return
    require(doc["solver_status"] in SOLVER_STATUS_WITH_VALUE, "apex_zero_invalid", "no field, no face")
    z = arrays["apex_zero"]
    require(z.dtype.str == "<i8" and z.ndim == 2 and z.shape[1:] == (2,) and z.shape[0] >= 1, "apex_zero_invalid",
            f"{z.dtype.str} {z.shape}")
    rows = apex_zero_rows(arrays)
    require(rows == sorted(set(rows)), "apex_zero_invalid", "rows sorted lexicographically and unique")
    tets, regions, region_of = problem_arrays["tets"], problem["regions"], problem_arrays["element_region"]
    interior = mesh_facets(tets)[0] if any(kind == 1 for kind, _ in rows) else []
    for kind, index in rows:
        if kind == 0:
            require(0 <= index < 4 * len(tets), "apex_zero_invalid", f"flow index {index}")
            element = index // 4
        else:
            require(kind == 1, "apex_zero_invalid", f"kind {kind}")
            require(doc["formulation"] == JUMP_ID, "apex_zero_invalid", "jump conditions exist in the jump formulation")
            require(0 <= index < 6 * len(interior), "apex_zero_invalid", f"jump index {index}")
            element = interior[index // 6][0][0]
        require(regions[int(region_of[element])].get("alpha", 0.0) > 0.0, "apex_zero_invalid",
                "Z holds inequality conditions only: alpha > 0 and tau > 0")
    if "apex_zero_dual" in arrays:
        dual = _array(arrays, "apex_zero_dual", "<f8", (len(rows), 6))
        for (kind, _), row in zip(rows, dual):
            require(kind == 0 or not np.any(row[3:]), "apex_zero_invalid", "a jump dual has three components")


def expected_program_shape(problem, arrays, query, bound, formulation):
    """Variables, rows, zero rows and cone list of the canonical program (interface 15.4)."""
    n_e, n_n = len(arrays["tets"]), len(arrays["nodes"])
    n_f = len(problem["footings"])
    free_components = 6 * n_f - len(query["restrained"])
    is_ray = query["form"] != "support"
    kinds = [item["kind"] for item in problem["boundary_sets"]]
    tets, facets, facet_set = arrays["tets"], arrays["facets"], arrays["facet_set"]
    interior, _ = mesh_facets(tets)
    if bound == "lower":
        zero = 3 * n_e + 9 * len(interior)
        for index in facet_set:
            zero += {"free": 9, "traction": 9, "symmetry": 6}.get(kinds[int(index)], 0)
        zero += free_components if is_ray else 0
        return (1 if is_ray else 0) + 24 * n_e, zero + 24 * n_e, zero, [{"dim": 6, "count": 4 * n_e}]
    broken = formulation == JUMP_ID
    node_sets = [set() for _ in kinds]
    for (e, f), index in zip(facets, facet_set):
        for a in FACET_NODES[int(f)]:
            node_sets[int(index)].add((int(e), a) if broken else int(tets[int(e)][a]))
    fixed = set().union(*[nodes for nodes, kind in zip(node_sets, kinds) if kind == "fixed"], set())
    zero = 24 * n_e + 3 * len(fixed)
    footing_nodes = {}
    for nodes, item in zip(node_sets, problem["boundary_sets"]):
        if item["kind"] == "symmetry":
            zero += len(nodes - fixed)
        elif item["kind"] == "footing":
            footing_nodes.setdefault(item["footing"], set()).update(nodes)
    zero += 3 * sum(len(nodes) for nodes in footing_nodes.values()) + len(query["restrained"])
    zero += 1 if is_ray else free_components
    cones = [{"dim": 6, "count": 4 * n_e}]
    if broken:
        zero += 18 * len(interior)
        cones.append({"dim": 3, "count": 6 * len(interior)})
        n_var = 30 * n_e + 6 * n_f + 24 * n_e + 18 * len(interior)
    else:
        n_var = 3 * n_n + 6 * n_f + 24 * n_e
    return n_var, zero + sum(c["dim"] * c["count"] for c in cones), zero, cones


def validate_program(doc, arrays, problem, problem_arrays):
    require(doc.get("format") == FORMAT_PROGRAM, "format", repr(doc.get("format")))
    require(doc.get("format_version") == FORMAT_VERSION, "format_version")
    _keys(doc, {
        "format", "format_version", "arrays", "problem_sha256", "query_id", "bound", "formulation",
        "objective_sign", "n_variables", "n_rows", "n_zero", "cones",
    })
    require(doc["problem_sha256"] == problem["problem_sha256"], "hash_mismatch")
    query = query_by_id(problem, doc["query_id"])
    bound = doc["bound"]
    require(bound in ("lower", "upper"), "bound", repr(bound))
    validate_formulation(bound, doc["formulation"], problem)
    require(doc["objective_sign"] == (-1 if bound == "lower" else 1), "program_objective_sign")
    n_var, n_row, n_zero, cones = expected_program_shape(problem, problem_arrays, query, bound, doc["formulation"])
    require(doc["cones"] == cones, "program_cone")
    require((doc["n_variables"], doc["n_rows"], doc["n_zero"]) == (n_var, n_row, n_zero), "program_shape")
    require(set(arrays) == {"A_indptr", "A_indices", "A_data", "b", "c"}, "array_names")
    indptr = _array(arrays, "A_indptr", "<i8", (n_var + 1,))
    require(arrays["A_indices"].ndim == 1 and arrays["A_indices"].shape == arrays["A_data"].shape, "array_shape")
    indices = _array(arrays, "A_indices", "<i8", arrays["A_indices"].shape)
    _array(arrays, "A_data", "<f8", arrays["A_data"].shape)
    _array(arrays, "b", "<f8", (n_row,))
    _array(arrays, "c", "<f8", (n_var,))
    require(indptr[0] == 0 and indptr[-1] == len(indices) and np.all(np.diff(indptr) >= 0), "program_csc")
    require(len(indices) == 0 or (indices.min() >= 0 and indices.max() < n_row), "program_csc")


def validate_feasibility(doc, problem):
    require(doc.get("format") == FORMAT_FEASIBILITY, "format", repr(doc.get("format")))
    require(doc.get("format_version") == FORMAT_VERSION, "format_version")
    _keys(doc, {
        "format", "format_version", "problem_sha256", "query_id", "formulation", "status", "rho_feas",
        "numerator", "denominator", "n_unknowns", "n_rows", "rank", "nullity", "footing_modes", "method",
        "witness_result_sha256", "witness_residual", "witness_projection_move", "witness_work_reference",
        "checker", "tolerances",
    })
    require(doc["problem_sha256"] == problem["problem_sha256"], "hash_mismatch")
    query = query_by_id(problem, doc["query_id"])
    validate_formulation("upper", doc["formulation"], problem)
    _keys(doc["checker"], {"name", "version"})
    require(doc["status"] in FEASIBILITY_STATUS, "feasibility_status", repr(doc["status"]))
    require(doc["tolerances"] == TOLERANCES, "tolerances")
    for key in ("rho_feas", "numerator", "denominator"):
        _finite_number(doc[key], "feasibility_number")
    require(doc["numerator"] >= 0.0 and doc["denominator"] > 0.0, "feasibility_number")
    rho = doc["rho_feas"]
    require(abs(rho - doc["numerator"] / doc["denominator"]) <= 1e-12 * max(1.0, abs(rho)), "feasibility_ratio")
    if doc["rank"] is not None:
        require(doc["nullity"] == doc["n_unknowns"] - doc["rank"], "feasibility_rank")
    else:
        require(doc["nullity"] is None and doc["footing_modes"] is None, "feasibility_rank", "skipped together")
    if not problem["footings"]:
        # Amendment 1, item 10: no footing, no V, so no footing modes and no condition on them (interface 11.1).
        require(doc["footing_modes"] is None, "feasibility_rank", "footing_modes without a footing")
    elif doc["rank"] is not None:
        footing_modes = doc["footing_modes"]
        require(isinstance(footing_modes, int) and not isinstance(footing_modes, bool) and footing_modes >= 0,
                "feasibility_rank", "computed together with the rank")
    witness, residual = doc["witness_result_sha256"], doc["witness_residual"]
    move, work = doc["witness_projection_move"], doc["witness_work_reference"]
    require((witness is None) == (residual is None) == (move is None), "feasibility_witness")
    is_ray = query["form"] in ("ray", "fixed_components")
    if residual is not None:
        _finite_number(residual, "feasibility_witness")
        require(isinstance(witness, str) and re.fullmatch("[0-9a-f]{64}", witness) is not None, "feasibility_witness")
        # Amendment 2, m2: the witness is a field the certificate procedure would accept as the same field.
        _finite_number(move, "feasibility_witness")
        require(move >= 0.0, "feasibility_witness")
        if is_ray:
            _finite_number(work, "feasibility_witness")
    if residual is None or not is_ray:
        require(work is None, "feasibility_witness", "W1 only for a ray witness")
    frictional = any(region.get("alpha", 0.0) > 0.0 for region in problem["regions"])
    ok, locked = TOLERANCES["TOL_FEAS_OK"], TOLERANCES["TOL_FEAS_LOCKED"]
    witness_passes = (residual is not None and residual <= ok and move <= TOLERANCES["TOL_PROJECTION_MOVE"]
                      and (not is_ray or work >= TOLERANCES["TOL_WORK_MIN"]))  # W1 / W_scale, as 17.4 step 5
    if frictional:
        # The minimisation is conic here, so the verdict rests on the witness alone (interface 11.1).
        require(residual is None or math.isclose(rho, residual, rel_tol=1e-12, abs_tol=1e-300), "feasibility_ratio")
    if doc["status"] == "formulation_feasible":
        require(rho <= ok and witness_passes, "feasibility_verdict")
        require(doc["footing_modes"] is None or doc["footing_modes"] >= 1, "feasibility_rank")
    elif doc["status"] == "formulation_infeasible":
        require(not frictional, "feasibility_verdict", "no infeasibility verdict without a linear minimisation")
        require(rho >= locked, "feasibility_verdict")
    else:
        require(not (rho <= ok and witness_passes), "feasibility_verdict")
        require(frictional or rho < locked, "feasibility_verdict")


def halfspace_boundary_is_covered(problem, arrays):
    """Every exterior facet is on the ground plane, in a truncation set or in a declared symmetry set."""
    domain = problem["domain"]
    axis, coordinate = AXES[domain["ground"]["axis"]], domain["ground"]["coordinate"]
    nodes, tets = arrays["nodes"], arrays["tets"]
    tol = TOLERANCES["TOL_GEOM_PLANE"] * float(np.linalg.norm(nodes.max(axis=0) - nodes.min(axis=0)))
    covered = set(domain["truncation_sets"]) | set(domain["symmetry_sets"])
    for (e, f), index in zip(arrays["facets"], arrays["facet_set"]):
        item = problem["boundary_sets"][int(index)]
        on_ground = bool(np.all(np.abs(nodes[tets[e][list(FACET_VERTICES[int(f)])], axis] - coordinate) <= tol))
        if on_ground:
            continue
        if item["kind"] in ("traction", "footing") or item["name"] not in covered:
            return False
    return True


def supported_domain_claim(problem, arrays, bound):
    """The strongest claim version 0 can support for a certified bound (interface 14)."""
    domain = problem["domain"]
    if domain["target"] == "box":
        return "box_stated"
    if not halfspace_boundary_is_covered(problem, arrays):
        return "unverified"
    kinds = {item["name"]: item["kind"] for item in problem["boundary_sets"]}
    walls = [kinds[name] for name in domain["truncation_sets"]]
    if bound == "upper":
        return "halfspace_by_zero_extension" if all(kind == "fixed" for kind in walls) else "unverified"
    weightless = all(not any(region["b0"]) and not any(region["b1"]) for region in problem["regions"])
    if all(kind == "free" for kind in walls) and weightless and domain["exterior_load_free"]:
        return "halfspace_by_zero_extension"
    return "unverified"


def _check_number(checks, key):
    value = checks.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        return None
    return float(value)


def certified_check_failures(doc, query_form=None):
    """Which numeric conditions of interface 17 a certified_bound does not show. Empty means all shown."""
    checks, t = doc["checks"], TOLERANCES
    numbers = {key: _check_number(checks, key) for key in CHECK_KEYS}
    required = ["omega_eq", "kappa_eq", "tau_eq", "projection_move", "repair_mu", "straightness_max"]
    if doc["bound"] == "lower":
        required.append("rho_cone")
    else:
        required.append("rho_flow")
        if doc["formulation"] == JUMP_ID:
            required.append("rho_jump_cone")
        if query_form in ("ray", "fixed_components"):
            required.append("work_reference")
    failures = [key for key in required if numbers[key] is None]
    if failures:
        return failures
    omega, kappa = numbers["omega_eq"], numbers["kappa_eq"]
    limits = [
        ("omega_eq", omega <= numbers["tau_eq"]), ("tau_eq", numbers["tau_eq"] == tau_eq()),
        ("kappa_eq", 1.0 <= kappa <= t["TOL_KAPPA_CAP"] and kappa * omega <= t["TOL_FLOOR_MAX"]),
        ("projection_move", 0.0 <= numbers["projection_move"] <= t["TOL_PROJECTION_MOVE"]),
        ("repair_mu", 0.0 <= numbers["repair_mu"] <= t["TOL_REPAIR_MAX"]),
        ("straightness_max", numbers["straightness_max"] <= t["TOL_GEOM_STRAIGHT"]),
    ]
    if doc["bound"] == "lower":
        limits.append(("rho_cone", numbers["rho_cone"] <= 0.0))
    else:
        limits.append(("rho_flow", numbers["rho_flow"] <= t["TOL_CONE_MARGIN"]))
        if doc["formulation"] == JUMP_ID:
            limits.append(("rho_jump_cone", numbers["rho_jump_cone"] <= t["TOL_CONE_MARGIN"]))
        if numbers["work_reference"] is not None or query_form in ("ray", "fixed_components"):
            limits.append(("work_reference", (numbers["work_reference"] or 0.0) > 0.0))
    if checks.get("rho_apex") is not None:
        # Amendment 2, M3: the Z arguments, kind-wise over D_scale, within the cone rounding margin (19.2.1 item 9).
        rho_apex = numbers["rho_apex"]
        limits.append(("rho_apex", doc["bound"] == "upper" and rho_apex is not None
                       and 0.0 <= rho_apex <= t["TOL_CONE_MARGIN"]))
    failures = [key for key, passed in limits if not passed]
    expected = max(t["FLOOR_MIN_REL"], kappa * omega) * abs(doc["value_certified"])
    if not math.isclose(doc["floor"], expected, rel_tol=1e-9, abs_tol=1e-300):
        failures.append("floor")
    return failures


def validate_certificate(doc, problem, problem_arrays, result):
    require(doc.get("format") == FORMAT_CERTIFICATE, "format", repr(doc.get("format")))
    require(doc.get("format_version") == FORMAT_VERSION, "format_version")
    _keys(doc, {
        "format", "format_version", "problem_sha256", "result_sha256", "query_id", "bound", "formulation",
        "classification", "reasons", "value_candidate", "value_certified", "floor", "domain_claim",
        "criteria", "flow_rule", "interface_laws", "checks", "checker", "tolerances",
    })
    require(doc["problem_sha256"] == problem["problem_sha256"], "hash_mismatch", "problem")
    require(doc["result_sha256"] == result["result_sha256"], "hash_mismatch", "result")
    require(doc["query_id"] == result["query_id"], "query_mismatch")
    require(doc["bound"] == result["bound"], "certificate_bound")
    require(doc["formulation"] == result["formulation"], "certificate_formulation")
    validate_formulation(doc["bound"], doc["formulation"], problem)
    require(doc["tolerances"] == TOLERANCES, "tolerances")
    require(set(doc["checks"]) == CHECK_KEYS, "certificate_checks")
    _keys(doc["checker"], {"name", "version"})
    require(doc["flow_rule"] == "associated", "certificate_flow_rule")
    require(doc["criteria"] == sorted({region["criterion"] for region in problem["regions"]}), "certificate_criteria")
    require(doc["interface_laws"] == sorted({item["interface_law"] for item in problem["footings"]}),
            "certificate_interface_laws")
    require(set(doc["reasons"]) <= REASON_CODES, "certificate_reason", str(doc["reasons"]))
    claim = doc["domain_claim"]
    require(claim in DOMAIN_CLAIMS, "domain_claim", repr(claim))
    require(claim == supported_domain_claim(problem, problem_arrays, doc["bound"]),
            "unsupported_domain_claim", repr(claim))
    classification = doc["classification"]
    require(classification in ("certified_bound", "estimate"), "certificate_classification", repr(classification))
    if classification == "certified_bound":
        require(not doc["reasons"], "certificate_reason", "certified with reasons")
        _finite_number(doc["value_candidate"], "certificate_value")
        _finite_number(doc["value_certified"], "certificate_value")
        _finite_number(doc["floor"], "certificate_value")
        require(doc["floor"] >= 0.0, "certificate_value")
        require(result["value"] is not None, "certificate_without_candidate")
        form = query_by_id(problem, doc["query_id"])["form"]
        failures = certified_check_failures(doc, form)
        require(not failures, "certificate_checks", str(failures))
    else:
        require(doc["reasons"], "certificate_reason", "estimate needs a reason")
        require(doc["value_certified"] is None and doc["floor"] is None, "certificate_value")


def certificate_hash(doc):
    """certificate_sha256 of amendment 1 (interface 15.1, 15.9): the whole certificate document, no arrays."""
    return document_hash(doc, {}, "certificate_sha256")


def interior_array_names(bound, formulation):
    if bound == "lower":
        return ["interior_sigma"]
    velocity = "interior_element_velocity" if formulation == JUMP_ID else "interior_velocity"
    return sorted([velocity, "interior_footing_velocity"])


def interior_field_sha256(bound, formulation, arrays):
    """Identity of the interior field (interface 15.9): the 15.1 hash of {} with the bound's interior arrays."""
    return document_hash({}, {name: arrays[name] for name in interior_array_names(bound, formulation)}, None)


def _count(value, code, minimum=0):
    require(isinstance(value, int) and not isinstance(value, bool) and value >= minimum, code, repr(value))


def validate_method_report(doc, certificate, problem=None, result=None, result_arrays=None, problem_arrays=None):
    """Method report of amendment 1 (interface 15.9), bound to its certificate by hash.

    Checks what the report and its certificate show: identity, the proved-bound arithmetic of 17.1, that no row is
    skipped without a proof, and the interior-field identity. With the problem, result and result arrays it also checks
    the interior field against them.
    """
    require(doc.get("format") == FORMAT_METHOD, "format", repr(doc.get("format")))
    require(doc.get("format_version") == FORMAT_VERSION, "format_version")
    _keys(doc, {
        "format", "format_version", "problem_sha256", "result_sha256", "certificate_sha256", "query_id", "bound",
        "formulation", "checker", "kappa", "dependent_rows", "rows_skipped", "interior_field", "apex_zero",
        "tolerances",
    })
    require(doc["certificate_sha256"] == certificate_hash(certificate), "method_certificate_hash")
    for key in ("problem_sha256", "result_sha256", "query_id", "bound", "formulation", "checker"):
        require(doc[key] == certificate[key], "method_identity", key)
    require(doc["tolerances"] == TOLERANCES, "tolerances")

    kappa = doc["kappa"]
    require(isinstance(kappa, dict), "method_kappa")
    _keys(kappa, {"method", "sigma_max_upper", "sigma_min_lower", "kappa_upper", "perturbation_bound", "residual",
                  "completeness"}, code="method_kappa")
    method = kappa["method"]
    require(method in KAPPA_SVD_METHODS | KAPPA_EIGEN_METHODS, "method_kappa", repr(method))
    for key in ("sigma_max_upper", "sigma_min_lower", "kappa_upper", "perturbation_bound"):
        _finite_number(kappa[key], "method_kappa")
    require(kappa["sigma_min_lower"] > 0.0 and kappa["perturbation_bound"] >= 0.0, "method_kappa")
    require(kappa["sigma_max_upper"] >= kappa["sigma_min_lower"], "method_kappa")
    if method in KAPPA_EIGEN_METHODS:
        _finite_number(kappa["residual"], "method_kappa")
        require(kappa["residual"] >= 0.0, "method_kappa", "an eigenvalue route reports its residual")
    else:
        require(kappa["residual"] is None and kappa["completeness"] is None, "method_kappa",
                "no residual or completeness record for a singular-value route")
    require(kappa["kappa_upper"] >= kappa["sigma_max_upper"] / kappa["sigma_min_lower"], "method_kappa_not_proved")
    require(certificate["checks"].get("kappa_eq") == kappa["kappa_upper"], "method_kappa_mismatch")

    rows = doc["dependent_rows"]
    require(isinstance(rows, dict), "method_dependent_rows")
    _keys(rows, {"n_rows", "n_unknowns", "rank", "proofs"}, code="method_dependent_rows")
    for key in ("n_rows", "n_unknowns", "rank"):
        _count(rows[key], "method_dependent_rows")
    require(isinstance(rows["proofs"], list), "method_dependent_rows")
    proved = 0
    for proof in rows["proofs"]:
        require(isinstance(proof, dict), "method_dependent_rows")
        _keys(proof, {"route", "count", "residual_bound"}, code="method_dependent_rows")
        require(proof["route"] in DEPENDENT_ROUTES, "method_dependent_rows", repr(proof["route"]))
        _count(proof["count"], "method_dependent_rows", minimum=1)
        bound = proof["residual_bound"]
        if proof["route"] == "shape_forced_zeros":
            require(bound is None, "method_dependent_rows", "shape-forced zeros need no residual")
        else:
            _finite_number(bound, "method_dependent_rows")
            require(bound >= 0.0, "method_dependent_rows")
            threshold = TOLERANCES["TOL_RANK"] * kappa["sigma_max_upper"]
            if proof["route"] == "null_space_basis":
                require(bound <= threshold, "method_dependent_rows", "null-space residual above the rank threshold")
            else:  # dense_svd_blocks: the factorisation error cannot move a value across the split (17.1, M1)
                require(bound <= threshold, "method_dependent_rows", "block perturbation above the rank threshold")
                require(kappa["sigma_min_lower"] > threshold + bound, "method_dependent_rows",
                        "kept values not separated from the skipped ones")
        proved += proof["count"]
    _count(doc["rows_skipped"], "method_rows_skipped")
    require(doc["rows_skipped"] == proved, "method_rows_skipped", "a row skipped without a proof")
    require(rows["rank"] == rows["n_rows"] - doc["rows_skipped"], "method_rows_skipped")
    require(rows["rank"] <= rows["n_unknowns"], "method_dependent_rows")
    if method in KAPPA_EIGEN_METHODS:
        # M1: a Ritz value bounds sigma_min only with a proof that no unskipped eigenvalue lies below it (17.1):
        # count_below eigenvalues of the normal matrix lie below `shift`, exactly the zeros skipped, and
        # sigma_min_lower^2 <= shift - residual.
        complete = kappa["completeness"]
        require(isinstance(complete, dict), "method_kappa_completeness", "an eigenvalue route proves completeness")
        _keys(complete, {"route", "shift", "count_below"}, code="method_kappa_completeness")
        require(complete["route"] in COMPLETENESS_ROUTES, "method_kappa_completeness", repr(complete["route"]))
        _finite_number(complete["shift"], "method_kappa_completeness")
        _count(complete["count_below"], "method_kappa_completeness")
        zeros = doc["rows_skipped"] if method == "eigen_normal_rows" else rows["n_unknowns"] - rows["rank"]
        require(complete["count_below"] == zeros, "method_kappa_completeness", "count below the shift != zeros skipped")
        require(complete["shift"] > kappa["residual"], "method_kappa_completeness", "shift not above the residual")
        require(kappa["sigma_min_lower"] ** 2 <= complete["shift"] - kappa["residual"], "method_kappa_not_proved",
                "sigma_min_lower^2 > shift - residual")

    interior = doc["interior_field"]
    require(isinstance(interior, dict), "method_interior")
    _keys(interior, {"source", "sha256"}, code="method_interior")
    source = interior["source"]
    require(source in INTERIOR_SOURCES, "method_interior", repr(source))
    repair = certificate["checks"].get("repair_mu")
    if source == "none":
        require(interior["sha256"] is None, "method_interior")
        require(repair in (None, 0, 0.0), "method_interior", "a repair needs an interior field")
    elif source == "zero":
        require(interior["sha256"] is None and doc["bound"] == "lower", "method_interior", "zero is a Lower field only")
        if problem is not None:
            query = query_by_id(problem, doc["query_id"])
            fixed_free = (all(not any(region["b0"]) for region in problem["regions"])
                          and all(region["k"] > 0.0 for region in problem["regions"])
                          and not any(any(row) for row in query.get("F0", []))
                          and (problem_arrays is None or not np.any(problem_arrays["t0"])))
            require(fixed_free, "method_interior", "zero is interior only without fixed loads (17.3 step 5)")
    else:
        sha = interior["sha256"]
        require(isinstance(sha, str) and re.fullmatch("[0-9a-f]{64}", sha) is not None, "method_interior")
        if result_arrays is not None:
            names = interior_array_names(doc["bound"], doc["formulation"])
            require(all(name in result_arrays for name in names), "method_interior", "no interior arrays in the result")
            require(sha == interior_field_sha256(doc["bound"], doc["formulation"], result_arrays), "method_interior")

    # Amendment 2 (interface 19.2.1): the declared apex face, its identity and the non-gating z_forced diagnostic.
    face = doc["apex_zero"]
    certified = certificate["classification"] == "certified_bound"
    if face is None:
        require(result_arrays is None or "apex_zero" not in result_arrays, "method_apex_zero",
                "the result declares Z and the report does not")
        require(certificate["checks"].get("rho_apex") is None, "method_apex_zero", "rho_apex without a face")
    else:
        # M3: a certificate of a Z-declaring result reports rho_apex.
        require(not certified or _check_number(certificate["checks"], "rho_apex") is not None, "method_apex_zero",
                "a certified Z-declaring result reports rho_apex")
        require(isinstance(face, dict), "method_apex_zero")
        _keys(face, {"count", "sha256", "z_forced"}, code="method_apex_zero")
        require(doc["bound"] == "upper", "method_apex_zero", "Z is refused for a Lower")
        _count(face["count"], "method_apex_zero", minimum=1)
        require(isinstance(face["sha256"], str) and re.fullmatch("[0-9a-f]{64}", face["sha256"]) is not None,
                "method_apex_zero")
        require(face["z_forced"] in Z_FORCED, "method_apex_zero", repr(face["z_forced"]))
        if result_arrays is not None:
            require("apex_zero" in result_arrays, "method_apex_zero", "a face the result does not declare")
            require(face["count"] == len(result_arrays["apex_zero"]), "method_apex_zero", "count")
            require(face["sha256"] == apex_zero_sha256(result_arrays), "method_apex_zero", "identity")
            require((face["z_forced"] == "absent") == ("apex_zero_dual" not in result_arrays), "method_apex_zero",
                    "z_forced is absent exactly when no dual is attached")
    if result is not None:
        require(result["result_sha256"] == doc["result_sha256"], "hash_mismatch", "result")


def apex_zero_sha256(arrays):
    """Identity of Z in a method report (amendment 2): the 15.1 hash of {} with the single array `apex_zero`."""
    return document_hash({}, {"apex_zero": arrays["apex_zero"]}, None)


def method_needs_result_arrays(method):
    """Amendment 2, m1: a report that names result arrays (interior field or Z) is verified only against them."""
    return method["interior_field"].get("source") == "result" or method.get("apex_zero") is not None


def validate_refusal(doc, problem, result):
    """Capacity or conditioning refusal of amendment 1 (interface 15.8): written instead of a certificate."""
    require(doc.get("format") == FORMAT_REFUSAL, "format", repr(doc.get("format")))
    require(doc.get("format_version") == FORMAT_VERSION, "format_version")
    _keys(doc, {
        "format", "format_version", "problem_sha256", "result_sha256", "query_id", "bound", "formulation", "kind",
        "stage", "detail", "resources", "budget", "checker", "tolerances",
    })
    require(doc["problem_sha256"] == problem["problem_sha256"], "hash_mismatch", "problem")
    require(doc["result_sha256"] == result["result_sha256"], "hash_mismatch", "result")
    require(doc["query_id"] == result["query_id"], "query_mismatch")
    require(doc["bound"] == result["bound"] and doc["formulation"] == result["formulation"], "refusal_identity")
    validate_formulation(doc["bound"], doc["formulation"], problem)
    require(doc["kind"] in REFUSAL_KINDS, "refusal_kind", repr(doc["kind"]))
    require(doc["stage"] in REFUSAL_STAGES, "refusal_stage", repr(doc["stage"]))
    require(isinstance(doc["detail"], str) and doc["detail"], "refusal_detail")
    require(isinstance(doc["resources"], dict), "refusal_resources")
    _keys(doc["resources"], {"n_rows", "n_unknowns", "nnz"}, code="refusal_resources")
    for value in doc["resources"].values():
        _count(value, "refusal_resources")
    # Amendment 2, m4: the budget the checker declared before it ran, and which one a capacity refusal exceeded.
    budget = doc["budget"]
    require(isinstance(budget, dict), "refusal_budget")
    _keys(budget, {"memory_bytes", "wall_s", "exceeded"}, code="refusal_budget")
    for key in ("memory_bytes", "wall_s"):
        _finite_number(budget[key], "refusal_budget")
        require(budget[key] > 0.0, "refusal_budget", key)
    if doc["kind"] == "capacity":
        require(budget["exceeded"] in BUDGET_EXCEEDED, "refusal_budget", "a capacity refusal names the budget exceeded")
    else:
        require(budget["exceeded"] is None, "refusal_budget", "a conditioning refusal ran within its budget")
    _keys(doc["checker"], {"name", "version"})
    require(doc["tolerances"] == TOLERANCES, "tolerances")


def is_gate_checker(record):
    return record["checker"]["name"] not in NON_GATE_CHECKERS


def gate_admissible(record, query_form="ray", method=None, result_arrays=None):
    """Only a real checker's certified bound or feasibility verdict can satisfy a gate. Every case query is a ray.

    Amendment 1: a certificate also needs a valid method report bound to it (interface 15.6, 15.9).
    Amendment 2, m1: a report that names result arrays (a `result` interior field, or a declared face Z) counts only
    when checked against the result's arrays, which must be supplied.
    """
    if not is_gate_checker(record):
        return False
    if record.get("format") == FORMAT_FEASIBILITY:
        return record["status"] == "formulation_feasible" and record["witness_result_sha256"] is not None
    if record["classification"] != "certified_bound" or certified_check_failures(record, query_form):
        return False
    if method is None:
        return False
    try:
        validate_method_report(method, record, result_arrays=result_arrays)
    except ContractError:
        return False
    return result_arrays is not None or not method_needs_result_arrays(method)


def bracket_claim(lower_claim, upper_claim):
    if lower_claim == upper_claim == "halfspace_by_zero_extension":
        return lower_claim
    return "box_stated"


def validate_bracket(doc, lower, upper, lower_method=None, upper_method=None):
    require(doc.get("format") == FORMAT_BRACKET, "format", repr(doc.get("format")))
    require(doc.get("format_version") == FORMAT_VERSION, "format_version")
    _keys(doc, {
        "format", "format_version", "problem_sha256", "query_id", "lower_formulation", "upper_formulation",
        "lower_result_sha256", "upper_result_sha256", "lower", "upper", "floor_lower", "floor_upper", "status",
        "gap", "domain_claim", "lower_domain_claim", "upper_domain_claim",
    })
    require(lower["bound"] == "lower" and upper["bound"] == "upper", "bracket_labels")
    require(doc["lower_formulation"] == lower["formulation"], "bracket_formulation")
    require(doc["upper_formulation"] == upper["formulation"], "bracket_formulation")
    require(lower["classification"] == upper["classification"] == "certified_bound", "bracket_uncertified")
    for certificate, method in ((lower, lower_method), (upper, upper_method)):
        if is_gate_checker(certificate):
            # Amendment 1, item 5: a gate certificate needs its method report (interface 15.6).
            require(method is not None, "method_report_missing", certificate["result_sha256"])
            validate_method_report(method, certificate)
    require(lower["problem_sha256"] == upper["problem_sha256"] == doc["problem_sha256"], "hash_mismatch")
    require(lower["query_id"] == upper["query_id"] == doc["query_id"], "query_mismatch")
    require(doc["lower_result_sha256"] == lower["result_sha256"], "hash_mismatch")
    require(doc["upper_result_sha256"] == upper["result_sha256"], "hash_mismatch")
    low, high = lower["value_certified"], upper["value_certified"]
    require(doc["lower"] == low and doc["upper"] == high, "bracket_values")
    require(doc["floor_lower"] == lower["floor"] and doc["floor_upper"] == upper["floor"], "bracket_values")
    floors = lower["floor"] + upper["floor"]
    expected = "bracket" if low <= high else ("order_within_floor" if low <= high + floors else "contradiction")
    require(doc["status"] == expected, "bracket_status", f"{doc['status']} != {expected}")
    require(doc["lower_domain_claim"] == lower["domain_claim"], "bracket_domain")
    require(doc["upper_domain_claim"] == upper["domain_claim"], "bracket_domain")
    require(doc["domain_claim"] == bracket_claim(lower["domain_claim"], upper["domain_claim"]), "bracket_domain")
    if expected == "bracket" and low > 0.0:
        _finite_number(doc["gap"], "bracket_gap")
        require(math.isclose(doc["gap"], (high - low) / (0.5 * (high + low)), rel_tol=1e-12, abs_tol=1e-15), "bracket_gap")
    else:
        require(doc["gap"] is None, "bracket_gap")


# ---------------------------------------------------------------------------------------------
# Case evaluator: the executable meaning of `assert`, `gate` and `record` lines (interface 17.9)
# ---------------------------------------------------------------------------------------------

def load_problem_hashes():
    return strict_json_loads((DATA / "problem_hashes.json").read_text(encoding="utf-8"))


def _formula_value(case, mesh, name, params):
    if name in TARGET_FORMULAS:
        return TARGET_FORMULAS[name](params)
    return band_formula(name, params, first_layer(case, mesh, params))


def _one_sided(kind, value, limit, tol):
    if kind.endswith("_le"):
        return value <= limit + tol * abs(limit)
    return value >= limit - tol * abs(limit)


_PROBLEM_CACHE = {}


def case_problem_cached(case, mesh):
    key = (case["name"], mesh, json.dumps(case["mesh"][mesh]))
    if key not in _PROBLEM_CACHE:
        _PROBLEM_CACHE[key] = build_case_problem(case, mesh)
    return _PROBLEM_CACHE[key]


WITHHOLDABLE_CAUSES = {"missing", "estimate"}
# Amendment 1: a refusal, a result nobody checked, a gate certificate without a valid method report, and an estimate
# refused for conditioning say nothing about whether a certificate could exist, so OD2 never withholds them.
NON_WITHHOLDABLE_CAUSES = {"contradiction", "refused", "unchecked", "unproved", "ill_conditioned"}


def _withholdable(case, item):
    """Owner decision OD2 (interface 11.3): only the gap limit and missing certificates of the 30-degree cases."""
    if case["class"] != "upper_gate_3d" or case.get("phi_deg") != 30.0 or item["cause"] in NON_WITHHOLDABLE_CAUSES:
        return False
    if item["kind"] == "gap_le":  # amendment 3: no P0 gate line is gap_le MAX_GAP_3D_DP now; the branch is kept
        return True
    return item["kind"] != "kinematic_status" and item["cause"] in WITHHOLDABLE_CAUSES


def evaluate_case(case, certificates=(), results=(), feasibility=(), methods=(), *, brackets=(), refusals=(),
                  result_arrays=None, problem_loader=None, hashes=None):
    """Verdict of one case from checker records (interface 17.9).

    Records are matched to meshes by problem hash and validated against the problem built from the case: a
    certificate with its result document, its method report (amendment 1), a feasibility record, a refusal, and any
    bracket record supplied. Invalid records are listed and excluded, and a case with an invalid record does not pass.
    Bracket status and gap are recomputed from two gate-admissible certificates; a bracket record is only checked,
    never trusted. A gate certificate counts only with a valid method report. A refusal from a gate checker turns a
    missing certificate into the cause `refused`, which OD2 never withholds.

    Every `formulation` line binds: a Lower assertion must hold for each Lower formulation, an Upper assertion or
    gate for each Upper formulation, and a pair item for each (Lower, Upper) pair. `control` lines bind nothing.
    Assertions hold on every listed mesh.
    """
    table = (load_problem_hashes() if hashes is None else hashes)[case["name"]]
    mesh_of = {sha: int(key[1:]) for key, sha in table.items()}
    load = problem_loader or (lambda mesh: case_problem_cached(case, mesh))
    meshes = range(len(case["mesh"]))
    lowers = [f for f in case["formulation"] if f in FORMULATIONS["lower"]]
    uppers = [f for f in case["formulation"] if f in FORMULATIONS["upper"]]
    by_result = {doc["result_sha256"]: doc for doc in results}
    offered = {(mesh_of[doc["problem_sha256"]], doc["formulation"]) for doc in results
               if doc.get("problem_sha256") in mesh_of}
    verdict = {"case": case["name"], "passed": False, "stop_program": False, "claim_withheld": False,
               "p0_blocking": True, "items": [], "records": [], "foreign_records": [], "invalid_records": []}

    def locate(doc, kind):
        mesh = mesh_of.get(doc.get("problem_sha256"))
        if mesh is None:
            verdict["foreign_records"].append({"kind": kind, "problem_sha256": doc.get("problem_sha256")})
        return mesh

    def invalid(kind, doc, code):
        verdict["invalid_records"].append({"kind": kind, "id": doc.get("result_sha256") or doc.get("problem_sha256"),
                                           "code": code})

    # Amendment 2, m1: result arrays, by result_sha256, let the gate path check the interior field and Z by hash.
    arrays_of = {}
    for sha, arrays in (result_arrays or {}).items():
        result = by_result.get(sha)
        if result is None or document_hash(result, arrays, "result_sha256") != sha:
            invalid("result_arrays", {"result_sha256": sha}, "hash_mismatch")
            continue
        mesh = mesh_of.get(result.get("problem_sha256"))
        if mesh is None:
            invalid("result_arrays", {"result_sha256": sha}, "hash_mismatch")
            continue
        try:  # Amendment 2, M2: the full result schema, the declared face included, before any gate reads it.
            validate_result(result, arrays, *load(mesh))
        except ContractError as error:
            invalid("result_arrays", {"result_sha256": sha}, error.code)
            continue
        arrays_of[sha] = arrays

    method_of = {}
    for doc in methods:
        if isinstance(doc, dict) and isinstance(doc.get("certificate_sha256"), str):
            method_of.setdefault(doc["certificate_sha256"], doc)
        else:
            invalid("method", doc, "format")
    by_key, by_sha, kinematic, refused, gate_method, gate_arrays = {}, {}, {}, set(), {}, {}
    for doc in certificates:
        mesh = locate(doc, "certificate")
        if mesh is None:
            continue
        problem, arrays = load(mesh)
        result = by_result.get(doc.get("result_sha256"))
        try:
            require(result is not None, "certificate_without_candidate", "no result document")
            validate_certificate(doc, problem, arrays, result)
        except ContractError as error:
            invalid("certificate", doc, error.code)
            continue
        method = method_of.get(certificate_hash(doc))
        found = arrays_of.get(doc["result_sha256"])
        if method is not None:
            try:
                validate_method_report(method, doc, problem, result, found, problem_arrays=arrays)
            except ContractError as error:
                invalid("method", method, error.code)
                method = None
        gate_method[id(doc)] = method
        gate_arrays[id(doc)] = found
        by_key.setdefault((mesh, doc["formulation"]), []).append(doc)
        by_sha[doc["result_sha256"]] = doc
    for doc in refusals:
        mesh = locate(doc, "refusal")
        if mesh is None:
            continue
        problem, _ = load(mesh)
        result = by_result.get(doc.get("result_sha256"))
        try:
            require(result is not None, "certificate_without_candidate", "no result document")
            validate_refusal(doc, problem, result)
            require(doc["result_sha256"] not in by_sha, "refusal_conflict", "a certificate exists for this result")
        except ContractError as error:
            invalid("refusal", doc, error.code)
            continue
        if is_gate_checker(doc):
            refused.add((mesh, doc["formulation"]))
    for doc in feasibility:
        mesh = locate(doc, "feasibility")
        if mesh is None:
            continue
        try:
            validate_feasibility(doc, load(mesh)[0])
        except ContractError as error:
            invalid("feasibility", doc, error.code)
            continue
        kinematic.setdefault((mesh, doc["formulation"]), []).append(doc)
    for doc in brackets:
        if locate(doc, "bracket") is None:
            continue
        try:
            lower, upper = by_sha.get(doc.get("lower_result_sha256")), by_sha.get(doc.get("upper_result_sha256"))
            require(lower is not None and upper is not None, "bracket_uncertified", "no valid certificate")
            validate_bracket(doc, lower, upper, gate_method.get(id(lower)), gate_method.get(id(upper)))
        except ContractError as error:
            invalid("bracket", doc, error.code)

    def add(kind, mesh, formulation, passed, cause="", detail=""):
        verdict["items"].append({"kind": kind, "mesh": mesh, "formulation": formulation,
                                 "status": "pass" if passed else "fail", "cause": "" if passed else cause or "violated",
                                 "detail": detail})

    def certified(mesh, formulation):
        """A gate-admissible certificate, or the cause of its absence (interface 11.3, 17.9; amendment 1)."""
        docs = by_key.get((mesh, formulation), [])
        good = [doc for doc in docs
                if gate_admissible(doc, method=gate_method.get(id(doc)), result_arrays=gate_arrays.get(id(doc)))]
        if good:
            return good[0], ""
        if (mesh, formulation) in refused:
            return None, "refused"
        if any(is_gate_checker(doc) and doc["classification"] == "certified_bound" for doc in docs):
            return None, "unproved"  # a gate certificate without a valid method report
        if any(doc["classification"] == "estimate" and "equality_ill_conditioned" in doc["reasons"] for doc in docs):
            return None, "ill_conditioned"
        if docs:
            return None, "estimate"
        return None, "unchecked" if (mesh, formulation) in offered else "missing"

    def pair(mesh, lower, upper):
        """Status and gap recomputed from two gate-admissible certificates (interface 17.9)."""
        low, why_low = certified(mesh, lower)
        high, why_high = certified(mesh, upper)
        if low is None or high is None:
            return None, why_low or why_high
        l_value, u_value = low["value_certified"], high["value_certified"]
        if l_value <= u_value:
            if l_value <= 0.0:
                return None, "no_gap"
            return (u_value - l_value) / (0.5 * (u_value + l_value)), ""
        if l_value <= u_value + low["floor"] + high["floor"]:
            return (0.0, "") if case["class"] == "exact_on_mesh" else (None, "order_within_floor")
        verdict["stop_program"] = True
        return None, "contradiction"

    def gap_item(kind, mesh, limit):
        for lower, upper in itertools.product(lowers, uppers):
            gap, why = pair(mesh, lower, upper)
            add(kind, mesh, f"{lower}__{upper}", gap is not None and gap <= limit, why, f"gap {gap!r}")

    side = {"lower_le_target": ("lower", "_le"), "upper_ge_target": ("upper", "_ge"), "lower_ge_target": ("lower", "_ge"),
            "lower_ge_formula": ("lower", "_ge"), "lower_le_formula": ("lower", "_le"),
            "upper_le_formula": ("upper", "_le"), "upper_ge_formula": ("upper", "_ge")}
    for tokens in case["assert"]:
        kind = tokens[0]
        if kind == "gravity_work_eq":
            continue  # an L0 operator check, tested directly
        for mesh in meshes:
            if kind == "gap_le":
                gap_item(kind, mesh, float(tokens[1]))
                continue
            bound, direction = side[kind]
            if kind.endswith("_target"):
                limit, tol = case["target"], float(tokens[1])
            else:
                limit = _formula_value(case, mesh, tokens[1], parse_key_values(tokens[2:]))
                tol = TOLERANCES["TOL_CASE_ONE_SIDED"]
            for formulation in (lowers if bound == "lower" else uppers):
                doc, why = certified(mesh, formulation)
                value = None if doc is None else doc["value_certified"]
                passed = value is not None and _one_sided(direction, value, limit, tol)
                add(kind, mesh, formulation, passed, why, f"value {value!r} limit {limit!r}")
    for tokens in case["gate"]:
        kind = tokens[0]
        if kind == "kinematic_status":
            for mesh, upper in itertools.product(meshes, uppers):
                docs = [d for d in kinematic.get((mesh, upper), []) if d["status"] == tokens[1] and gate_admissible(d)]
                add(kind, mesh, upper, bool(docs), "missing")
        elif kind == "certified_upper":
            mesh = int(tokens[2])
            for upper in uppers:
                doc, why = certified(mesh, upper)
                add(kind, mesh, upper, doc is not None, why)
        elif kind == "certified_bracket":
            mesh = int(tokens[2])
            for lower, upper in itertools.product(lowers, uppers):
                gap, why = pair(mesh, lower, upper)
                add(kind, mesh, f"{lower}__{upper}", why in ("", "no_gap"), why)
        elif kind == "gap_le":
            gap_item(kind, int(tokens[3]), TOLERANCES[tokens[1]])
        elif kind in ("upper_nonincrease", "lower_nondecrease"):
            coarse, fine = int(tokens[1]), int(tokens[2])
            for formulation in (uppers if kind == "upper_nonincrease" else lowers):
                first, why_first = certified(coarse, formulation)
                second, why_second = certified(fine, formulation)
                if first is None or second is None:
                    add(kind, fine, formulation, False, why_first or why_second)
                    continue
                a, b = first["value_certified"], second["value_certified"]
                allowance = TOLERANCES["TOL_MONOTONE"] * abs(a)
                passed = b <= a + allowance if kind == "upper_nonincrease" else b >= a - allowance
                add(kind, fine, formulation, passed, "violated", f"coarse {a!r} fine {b!r}")
        else:
            fail("case_format", kind)
    for tokens in case["record"]:
        require(tokens[0] == "gap" and tokens[1] == "mesh", "case_format", str(tokens))
        mesh = int(tokens[2])
        for lower, upper in itertools.product(lowers, uppers):
            gap, why = pair(mesh, lower, upper)
            verdict["records"].append({"kind": "gap", "mesh": mesh, "formulation": f"{lower}__{upper}",
                                       "gap": gap, "detail": why})
    failures = [item for item in verdict["items"] if item["status"] == "fail"]
    clean = not verdict["stop_program"] and not verdict["invalid_records"]
    verdict["passed"] = bool(verdict["items"]) and clean and not failures
    verdict["claim_withheld"] = clean and bool(failures) and all(_withholdable(case, item) for item in failures)
    verdict["p0_blocking"] = not verdict["passed"] and not verdict["claim_withheld"]
    return verdict


def evaluate_program(verdicts, required=()):
    """P0 reading of case verdicts: a contradiction stops the program, a blocking case keeps P0 open, and the
    30-degree Drucker-Prager claim may be withheld (owner decision OD2). A refused item keeps P0 open (amendment 1)."""
    seen = {verdict["case"] for verdict in verdicts}
    missing = sorted(set(required) - seen)
    stopped = any(verdict["stop_program"] for verdict in verdicts)
    blocking = sorted(verdict["case"] for verdict in verdicts if verdict["p0_blocking"])
    return {"stop_program": stopped, "p0_open": stopped or bool(blocking) or bool(missing),
            "blocking_cases": blocking, "missing_cases": missing,
            "withheld_claims": sorted(verdict["case"] for verdict in verdicts if verdict["claim_withheld"]),
            "refused_cases": sorted(verdict["case"] for verdict in verdicts
                                    if any(item["cause"] == "refused" for item in verdict["items"]))}


# ---------------------------------------------------------------------------------------------
# Fixtures: exact fields on the smallest meshes. No solver is involved.
# ---------------------------------------------------------------------------------------------

GROUND_PATCH_CASE = """
case vm_ground_patch_halfspace
units       length m; force kN; stress kPa
box         -1 -1 -1 1 1 0
mesh        2 2 1
region      soil von_mises k 1.0
boundary    xmin fixed plane x -1
boundary    xmax fixed plane x 1
boundary    ymin fixed plane y -1
boundary    ymax fixed plane y 1
boundary    base fixed plane z -1
boundary    ground traction plane z 0 t1 -1 0 0
query       q ray
end
"""

FREE_APEX_CASE = """
case dp_free_hydrostatic
units       length m; force kN; stress kPa
box         0 0 0 1 1 1
mesh        1 1 1
region      soil drucker_prager k 1.0 alpha 0.2
boundary    x0 traction plane x 0 t1 -1 0 0
boundary    x1 traction plane x 1 t1 1 0 0
boundary    y0 traction plane y 0 t1 0 -1 0
boundary    y1 traction plane y 1 t1 0 1 0
boundary    z0 traction plane z 0 t1 0 0 -1
boundary    z1 traction plane z 1 t1 0 0 1
query       q ray
end
"""


def facet_work_table(area):
    table = np.zeros((3, 6))
    for i in range(3):
        for j in range(3):
            table[i, j] = area / 30.0 if i == j else -area / 60.0
        for q, (a, b) in enumerate(((0, 1), (1, 2), (2, 0))):
            table[i, 3 + q] = 2.0 * area / 15.0 if i in (a, b) else area / 15.0
    return table


def dissipation_density(d, k, alpha):
    """The value formulas of interface 4.3, evaluated without testing the domain (17.4 step 1)."""
    theta, enorm = strain_split(d)
    return S2 * k * enorm if alpha == 0.0 else k * theta / (3.0 * alpha)


def element_velocities(arrays, velocity):
    return velocity if velocity.ndim == 3 else velocity[arrays["tets"]]


def facet_nodal_jump(tets, element_velocity, e1, f1, e2):
    """v(e+) - v(e-) at the six facet nodes, in the facet-node order of e-."""
    local1 = list(FACET_NODES[f1])
    order = [int(v) for v in tets[e2]]
    local2 = [order.index(int(tets[e1][a])) for a in local1]
    return element_velocity[e2][local2] - element_velocity[e1][local1]


def fixture_upper_value(problem, arrays, velocity, footing_velocity, quadrature="vertex", body="exact"):
    """Value of a P2 fixture field, continuous (N,3) or broken (E,10,3). Fixture tooling, not the checker."""
    dissipation, w0, w1 = fixture_upper_terms(problem, arrays, velocity, footing_velocity, quadrature, body)
    return (dissipation - w0) / w1


def fixture_upper_objective(problem, arrays, velocity, footing_velocity):
    """c^T x of the field with every cone variable at its minimum: D - W0 (amendment 1, item 2)."""
    dissipation, w0, _ = fixture_upper_terms(problem, arrays, velocity, footing_velocity)
    return dissipation - w0


def fixture_upper_terms(problem, arrays, velocity, footing_velocity, quadrature="vertex", body="exact"):
    """(D, W0, W1) of a P2 fixture field: vertex and Bernstein sums, exact work integrals."""
    nodes, tets = arrays["nodes"], arrays["tets"]
    query = problem["queries"][0]
    field = element_velocities(arrays, velocity)
    dissipation, w0, w1 = 0.0, 0.0, 0.0
    for e, tet in enumerate(tets):
        region = problem["regions"][int(arrays["element_region"][e])]
        k, alpha = region["k"], region.get("alpha", 0.0)
        _, volume = tet_gradients(nodes[tet[:4]])
        rates = [strain_rate_at_vertex(nodes[tet], field[e], m) for m in range(4)]
        if quadrature == "vertex":
            dissipation += sum(volume / 4.0 * dissipation_density(d, k, alpha) for d in rates)
        else:
            dissipation += volume * dissipation_density(np.mean(rates, axis=0), k, alpha)
        weights = volume * (np.array([-0.05] * 4 + [0.2] * 6) if body == "exact" else np.array([0.25] * 4 + [0.0] * 6))
        moment = weights @ field[e]
        w0 += np.dot(region["b0"], moment)
        w1 += np.dot(region["b1"], moment)
    if velocity.ndim == 3:
        for (e1, f1), (e2, _) in mesh_facets(tets)[0]:
            region = problem["regions"][int(arrays["element_region"][e1])]
            k, alpha = region["k"], region.get("alpha", 0.0)
            normal, area = facet_geometry(nodes[tets[e1][:4]], f1)
            coefficients = bernstein_from_nodal(facet_nodal_jump(tets, field, e1, f1, e2))
            opening = coefficients @ normal
            if alpha > 0.0:
                dissipation += area / 6.0 * k / (3.0 * alpha) * float(np.sum(opening))
            else:
                slip = coefficients - np.outer(opening, normal)
                dissipation += area / 6.0 * k * float(np.sum(np.linalg.norm(slip, axis=1)))
    for row, (e, f) in enumerate(arrays["facets"]):
        _, area = facet_geometry(nodes[tets[e][:4]], int(f))
        table = facet_work_table(area)
        facet_velocity = field[e][list(FACET_NODES[int(f)])]
        w0 += float(np.sum(arrays["t0"][row] * (table @ facet_velocity)))
        w1 += float(np.sum(arrays["t1"][row] * (table @ facet_velocity)))
    w0 += float(np.sum(np.array(query["F0"]).reshape(-1, 6) * footing_velocity))
    w1 += float(np.sum(np.array(query["F1"]).reshape(-1, 6) * footing_velocity))
    return dissipation, w0, w1


def fixture_resultant(problem, arrays, sigma):
    """Footing resultants of a stress field (interface 6). Fixture tooling, not the checker."""
    nodes, tets = arrays["nodes"], arrays["tets"]
    out = np.zeros((len(problem["footings"]), 6))
    for (e, f), index in zip(arrays["facets"], arrays["facet_set"]):
        item = problem["boundary_sets"][int(index)]
        if item["kind"] != "footing":
            continue
        normal, area = facet_geometry(nodes[tets[e][:4]], int(f))
        vertices = FACET_VERTICES[int(f)]
        x = nodes[tets[e][list(vertices)]] - np.array(problem["footings"][item["footing"]]["x_ref"])
        traction = [voigt_n(normal) @ sigma[e, a] for a in vertices]
        for b in range(3):
            out[item["footing"], :3] += area / 3.0 * traction[b]
            for a in range(3):
                out[item["footing"], 3:] += (area / 6.0 if a == b else area / 12.0) * np.cross(x[a], traction[b])
    return out


def signed_flow_conditions(problem, arrays, velocity):
    """c~ of interface 17.4 at every (e, i), -inf where alpha = 0, and the largest |D|_F."""
    nodes, tets = arrays["nodes"], arrays["tets"]
    field = element_velocities(arrays, velocity)
    out, size = np.full((len(tets), 4), -math.inf), 0.0
    for e, tet in enumerate(tets):
        alpha = problem["regions"][int(arrays["element_region"][e])].get("alpha", 0.0)
        for m in range(4):
            d = strain_rate_at_vertex(nodes[tet], field[e], m)
            theta, enorm = strain_split(d)
            size = max(size, math.sqrt(float(d[:3] @ d[:3]) + 0.5 * float(d[3:] @ d[3:])))
            if alpha > 0.0:
                out[e, m] = 3.0 * S2 * alpha * enorm - theta
    return out, size


def signed_jump_conditions(problem, arrays, element_velocity):
    """Per interior facet: tau |c_t| - n.c at the six nodal jumps and at the six Bernstein coefficients."""
    nodes, tets = arrays["nodes"], arrays["tets"]
    alpha = problem["regions"][0].get("alpha", 0.0)
    tau, _ = jump_constants(alpha)
    rows = []
    for (e1, f1), (e2, _) in mesh_facets(tets)[0]:
        normal, _ = facet_geometry(nodes[tets[e1][:4]], f1)
        nodal = facet_nodal_jump(tets, element_velocity, e1, f1, e2)

        def signed(values):
            opening = values @ normal
            return tau * np.linalg.norm(values - np.outer(opening, normal), axis=1) - opening
        rows.append((signed(nodal), signed(bernstein_from_nodal(nodal))))
    return rows


def facet_coefficients(problem_arrays, element_velocity):
    """Normals (F, 3) and Bernstein jump coefficients (F, 6, 3) on every interior facet, in 3.4 order."""
    nodes, tets = problem_arrays["nodes"], problem_arrays["tets"]
    normals, coefficients = [], []
    for (e1, f1), (e2, _) in mesh_facets(tets)[0]:
        normals.append(facet_geometry(nodes[tets[e1][:4]], f1)[0])
        coefficients.append(bernstein_from_nodal(facet_nodal_jump(tets, element_velocity, e1, f1, e2)))
    return np.array(normals), np.array(coefficients)


def apex_zero_row_values(problem_arrays, element_velocity, rows):
    """The exact rows Z adds (amendment 2, 19.2.1): d_ei (6) per flow entry, c_B (3, Cartesian) per jump entry."""
    nodes, tets = problem_arrays["nodes"], problem_arrays["tets"]
    _, coefficients = facet_coefficients(problem_arrays, element_velocity)
    out = []
    for kind, index in rows:
        if kind == 0:
            e, i = divmod(index, 4)
            out.extend(strain_rate_at_vertex(nodes[tets[e]], element_velocity[e], i))
        else:
            out.extend(coefficients[index // 6, index % 6])
    return np.array(out)


def project_onto_apex_rows(problem_arrays, element_velocity, rows):
    """Minimum-norm projection of a broken field onto the Z rows alone, for a fixture whose problem has no other
    equality row. Every unknown has the same scale v_ref, so this is the scaled correction of 17.1. Fixture tooling."""
    shape = element_velocity.shape
    x = element_velocity.reshape(-1)
    matrix = np.column_stack([apex_zero_row_values(problem_arrays, unit.reshape(shape), rows)
                              for unit in np.eye(x.size)])
    correction = np.linalg.lstsq(matrix, matrix @ x, rcond=None)[0]
    return (x - correction).reshape(shape)


def signed_jump_coefficients(problem, problem_arrays, element_velocity):
    """c~_FB = tau |c_t| - n . c_B at every (F, B), flattened as 6 F + B, and D_scale_jump (17.4 step 3)."""
    tau, _ = jump_constants(problem["regions"][0].get("alpha", 0.0))
    normals, coefficients = facet_coefficients(problem_arrays, element_velocity)
    opening = np.einsum("fbk,fk->fb", coefficients, normals)
    slip = np.linalg.norm(coefficients - opening[..., None] * normals[:, None, :], axis=2)
    return (tau * slip - opening).reshape(-1), float(np.max(np.linalg.norm(coefficients, axis=2)))


def _result(name, problem, bound, value, arrays, formulation=None, mesh=0, objective=None):
    """objective is c^T x: -lambda for Lower, D - W0 for Upper; it equals value only when W1 = 1 (amendment 1)."""
    formulation = formulation or (LOWER_ID if bound == "lower" else CONTROL_ID)
    if objective is None:
        objective = -value if bound == "lower" else value
    doc = {
        "format": FORMAT_RESULT, "format_version": FORMAT_VERSION,
        "arrays": f"{record_stem(name, formulation, mesh)}.result.npz", "problem_sha256": problem["problem_sha256"],
        "query_id": "q", "bound": bound, "formulation": formulation,
        "producer": {"name": "analytic-fixture", "version": "0"},
        "solver": {"name": "none", "version": "0", "settings": {}},
        "solver_status": "solved", "objective": objective,
        "value": value, "wall_time_s": None,
    }
    doc["result_sha256"] = document_hash(doc, arrays, "result_sha256")
    return doc


def fixture_checks(bound, formulation):
    """Hand-written metrics of an exact fixture. The checker name keeps them away from every gate."""
    checks = {key: None for key in sorted(CHECK_KEYS)}
    checks.update(omega_eq=0.0, kappa_eq=1.0, tau_eq=tau_eq(), projection_move=0.0, repair_mu=0.0, straightness_max=0.0)
    if bound == "lower":
        checks["rho_cone"] = 0.0
    else:
        checks.update(rho_flow=0.0, work_reference=1.0)
        if formulation == JUMP_ID:
            checks["rho_jump_cone"] = 0.0
    return checks


def _certificate(problem, arrays, result):
    return {
        "format": FORMAT_CERTIFICATE, "format_version": FORMAT_VERSION,
        "problem_sha256": problem["problem_sha256"], "result_sha256": result["result_sha256"],
        "query_id": result["query_id"], "bound": result["bound"], "formulation": result["formulation"],
        "classification": "certified_bound", "reasons": [], "value_candidate": result["value"],
        "value_certified": result["value"], "floor": TOLERANCES["FLOOR_MIN_REL"] * abs(result["value"]),
        "domain_claim": supported_domain_claim(problem, arrays, result["bound"]),
        "criteria": sorted({region["criterion"] for region in problem["regions"]}),
        "flow_rule": "associated",
        "interface_laws": sorted({item["interface_law"] for item in problem["footings"]}),
        "checks": fixture_checks(result["bound"], result["formulation"]),
        "checker": {"name": SCHEMA_FIXTURE_CHECKER, "version": "0"},
        "tolerances": dict(TOLERANCES),
    }


def _method_report(certificate, n_rows=4, n_unknowns=10):
    """Hand-written method report of a fixture certificate (amendment 1, interface 15.9). Its numbers are
    schema-only: a single dense block, no dependent row, no repair, and kappa_upper = sigma_max / sigma_min."""
    kappa = certificate["checks"]["kappa_eq"]
    return {
        "format": FORMAT_METHOD, "format_version": FORMAT_VERSION,
        "problem_sha256": certificate["problem_sha256"], "result_sha256": certificate["result_sha256"],
        "certificate_sha256": certificate_hash(certificate), "query_id": certificate["query_id"],
        "bound": certificate["bound"], "formulation": certificate["formulation"],
        "checker": copy.deepcopy(certificate["checker"]),
        "kappa": {"method": "dense_svd", "sigma_max_upper": kappa, "sigma_min_lower": 1.0, "kappa_upper": kappa,
                  "perturbation_bound": 0.0, "residual": None, "completeness": None},
        "dependent_rows": {"n_rows": n_rows, "n_unknowns": n_unknowns, "rank": n_rows, "proofs": []},
        "rows_skipped": 0,
        "interior_field": {"source": "none", "sha256": None},
        "apex_zero": None,
        "tolerances": dict(TOLERANCES),
    }


def _bracket(lower, upper):
    low, high = lower["value_certified"], upper["value_certified"]
    return {
        "format": FORMAT_BRACKET, "format_version": FORMAT_VERSION,
        "problem_sha256": lower["problem_sha256"], "query_id": lower["query_id"],
        "lower_formulation": lower["formulation"], "upper_formulation": upper["formulation"],
        "lower_result_sha256": lower["result_sha256"], "upper_result_sha256": upper["result_sha256"],
        "lower": low, "upper": high, "floor_lower": lower["floor"], "floor_upper": upper["floor"],
        "status": "bracket", "gap": (high - low) / (0.5 * (high + low)),
        "domain_claim": bracket_claim(lower["domain_claim"], upper["domain_claim"]),
        "lower_domain_claim": lower["domain_claim"], "upper_domain_claim": upper["domain_claim"],
    }


def _feasibility(problem, witness):
    return {
        "format": FORMAT_FEASIBILITY, "format_version": FORMAT_VERSION,
        "problem_sha256": problem["problem_sha256"], "query_id": "q", "formulation": witness["formulation"],
        "status": "formulation_feasible", "rho_feas": 0.0, "numerator": 0.0, "denominator": 1.0,
        "n_unknowns": 10, "n_rows": 4, "rank": 4, "nullity": 6, "footing_modes": 1 if problem["footings"] else None,
        "method": "analytic witness; counts and denominator are schema-only",
        "witness_result_sha256": witness["result_sha256"], "witness_residual": 0.0,
        "witness_projection_move": 0.0, "witness_work_reference": 1.0,
        "checker": {"name": SCHEMA_FIXTURE_CHECKER, "version": "0"}, "tolerances": dict(TOLERANCES),
    }


def _exact_fields(name, arrays):
    n_e, x = len(arrays["tets"]), arrays["nodes"]
    sigma = np.zeros((n_e, 4, 6))
    velocity = np.zeros_like(x)
    if name == "vm_uniaxial_tension":
        sigma[..., 0] = S3
        velocity[:, 0], velocity[:, 1], velocity[:, 2] = x[:, 0], -0.5 * x[:, 1], -0.5 * x[:, 2]
        return S3, sigma, np.zeros((0, 6)), velocity, np.zeros((0, 6))
    sigma[..., 3] = 1.0
    velocity[:, 0] = x[:, 1]
    plate = np.array([[1.0, 0, 0, 0, 0, 0]])
    return 1.0, sigma, plate.copy(), velocity, plate.copy()


def _float_arrays(**arrays):
    return {name: value.astype("<f8") for name, value in arrays.items()}


def _with_hash(doc, arrays):
    doc["problem_sha256"] = document_hash({k: v for k, v in doc.items() if k != "problem_sha256"}, arrays, "problem_sha256")
    return doc


def _store_result(documents, name, problem, bound, value, arrays, formulation=None, objective=None):
    result = _result(name, problem, bound, value, arrays, formulation, objective=objective)
    documents[f"{record_stem(name, result['formulation'])}.result"] = (result, arrays)
    return result


def _store_problem(documents, name, problem, arrays):
    documents[f"{case_stem(name)}.problem"] = (problem, arrays)


def _add_exact_set(documents, name, problem, arrays):
    value, sigma, resultant, velocity, footing_velocity = _exact_fields(name, arrays)
    _store_problem(documents, name, problem, arrays)
    lower = _store_result(documents, name, problem, "lower", value, _float_arrays(sigma=sigma, resultant=resultant))
    upper = _store_result(documents, name, problem, "upper", value,
                          _float_arrays(velocity=velocity, footing_velocity=footing_velocity))
    _store_result(documents, name, problem, "upper", value,
                  _float_arrays(element_velocity=velocity[arrays["tets"]], footing_velocity=footing_velocity), JUMP_ID)
    lower_cert = _certificate(problem, arrays, lower)
    upper_cert = _certificate(problem, arrays, upper)
    documents[f"{record_stem(name, LOWER_ID)}.certificate"] = (lower_cert, None)
    documents[f"{record_stem(name, CONTROL_ID)}.certificate"] = (upper_cert, None)
    documents[f"{record_stem(name, LOWER_ID)}.certificate.method"] = (_method_report(lower_cert), None)
    documents[f"{record_stem(name, CONTROL_ID)}.certificate.method"] = (_method_report(upper_cert), None)
    documents[f"{bracket_stem(name, LOWER_ID, CONTROL_ID)}.bracket"] = (_bracket(lower_cert, upper_cert), None)
    documents[f"{record_stem(name, CONTROL_ID)}.feasibility"] = (_feasibility(problem, upper), None)


def ground_patch_field(arrays):
    """v = curl(h^3 e_y) with h the hat function of the centre ground vertex.

    Quadratic per element, divergence free, with continuous normal component because the tangential trace
    of the potential is continuous. It vanishes on the link of the vertex, so on every wall and the base.
    """
    nodes, tets = arrays["nodes"], arrays["tets"]
    vertices = np.unique(tets[:, :4])
    centre = int(vertices[np.argmin(np.linalg.norm(nodes[vertices], axis=1))])
    field = np.zeros((len(tets), 10, 3))
    for e, tet in enumerate(tets):
        where = np.nonzero(tet[:4] == centre)[0]
        if len(where) == 0:
            continue
        i = int(where[0])
        grads, _ = tet_gradients(nodes[tet[:4]])
        hat = np.zeros(10)
        hat[i] = 1.0
        for q, (a, b) in enumerate(EDGES):
            if i in (a, b):
                hat[4 + q] = 0.5
        field[e] = 3.0 * (hat ** 2)[:, None] * np.cross(grads[i], [0.0, 1.0, 0.0])
    return field


def _first_violating_direction(problem, arrays, base, node, amount):
    for component in range(3):
        trial = base.copy()
        trial[node, component] += amount
        conditions, size = signed_flow_conditions(problem, arrays, trial)
        if conditions.max() / size > 0.1 * amount:
            return component
    fail("negative_fixture", "no violating direction")


APEX_FIXTURE_DELTA, APEX_FIXTURE_EPS, APEX_FIXTURE_HIDE = 0.1, 1e-8, 1e-4


def apex_fixture_fields(problem_arrays, base, side, direction, element, local_node, perturbation):
    """(candidate, interior) of the amendment-2 omission fixture: the free Drucker-Prager cube v = x with one side of
    the plane x = y translated so that plane opens, then one vertex copy of one element moved by `perturbation`.
    The interior field is the unperturbed one: its jump is exactly zero on every other interior facet."""
    interior = base.copy()
    interior[np.array(side)] += APEX_FIXTURE_DELTA * np.array(direction)
    candidate = interior.copy()
    candidate[element, local_node] += np.array(perturbation)
    return candidate, interior


def _apex_omission(problem, problem_arrays, base, side, direction):
    """Search the smallest instance of a Z that omits one needed condition (amendment 2, 19.2.1).

    Z_full is every jump condition the interior field leaves at the apex (all interior facets off the opening
    plane). With Z_full the projected candidate has no condition in P; with Z = Z_full minus one condition of P, that
    condition is still in P after projection while the interior field sits at the apex there."""
    margin = TOLERANCES["TOL_CONE_MARGIN"]
    normals, _ = facet_coefficients(problem_arrays, base)
    plane = np.abs(normals @ np.array(direction)) > 1.0 - 1e-12
    z_full = [(1, 6 * f + b) for f in range(len(normals)) if not plane[f] for b in range(6)]
    for element in range(len(problem_arrays["tets"])):
        for local_node in range(4):
            for axis in range(3):
                perturbation = [0.0, 0.0, 0.0]
                perturbation[axis] = APEX_FIXTURE_EPS
                candidate, _ = apex_fixture_fields(problem_arrays, base, side, direction, element, local_node,
                                                   perturbation)
                signed, scale = signed_jump_coefficients(problem, problem_arrays, candidate)
                in_p = [(1, int(index)) for index in np.nonzero(signed > margin * scale)[0]]
                full = project_onto_apex_rows(problem_arrays, candidate, z_full)
                signed_full, scale_full = signed_jump_coefficients(problem, problem_arrays, full)
                if not in_p or np.max(signed_full) > margin * scale_full:
                    continue
                for omitted in in_p:
                    rows = [row for row in z_full if row != omitted]
                    projected = project_onto_apex_rows(problem_arrays, candidate, rows)
                    after, after_scale = signed_jump_coefficients(problem, problem_arrays, projected)
                    if after[omitted[1]] > 1e3 * margin * after_scale:
                        return {"element": element, "local_node": local_node, "perturbation": perturbation,
                                "apex_zero": [list(row) for row in rows], "omitted": list(omitted)}
    fail("negative_fixture", "no apex omission instance")


def build_negative_fixtures(documents):
    P = case_stem("vm_plate_shear")
    D = case_stem("dp_uniaxial_tension")
    G = case_stem("vm_dead_gravity_shear_layer")
    H = case_stem("vm_ground_patch_halfspace")
    A = case_stem("dp_free_hydrostatic")
    L = f"{record_stem('vm_plate_shear', LOWER_ID)}.result"
    U = f"{record_stem('vm_plate_shear', CONTROL_ID)}.result"
    J = f"{record_stem('vm_plate_shear', JUMP_ID)}.result"
    GL = f"{record_stem('vm_dead_gravity_shear_layer', LOWER_ID)}.result"
    GU = f"{record_stem('vm_dead_gravity_shear_layer', CONTROL_ID)}.result"
    DU = f"{record_stem('dp_uniaxial_tension', CONTROL_ID)}.result"
    HU = f"{record_stem('vm_ground_patch_halfspace', JUMP_ID)}.result"
    AU = f"{record_stem('dp_free_hydrostatic', JUMP_ID)}.result"
    T = case_stem("vm_uniaxial_tension")
    TJ = f"{record_stem('vm_uniaxial_tension', JUMP_ID)}.result"
    R = case_stem("vm_plate_shear_restrained")
    RU = f"{record_stem('vm_plate_shear_restrained', CONTROL_ID)}.result"

    plate, plate_doc = documents[f"{P}.problem"][1], documents[f"{P}.problem"][0]
    nodes, tets = plate["nodes"], plate["tets"]
    centre = int(np.argmin(np.linalg.norm(nodes - 0.5, axis=1)))
    footprint = int(np.nonzero((np.abs(nodes[:, 1] - 1.0) < 1e-12) & (np.abs(nodes[:, 0]) < 1e-12))[0][0])
    base_node = int(np.nonzero((np.abs(nodes[:, 1]) < 1e-12) & (np.abs(nodes[:, 0] - 1.0) < 1e-12))[0][0])
    grads, _ = tet_gradients(nodes[tets[0][:4]])
    vertex = int(np.argmax(np.abs(grads[:, 0])))
    amount = 3.0 * tau_eq() * float(np.sum(np.abs(grads))) / abs(grads[vertex, 0])

    # Amendment 1, item 9: an edge of one element whose two facets are both interior. Moving its mid-edge node along
    # grad L_a x grad L_b keeps theta = 0 at every vertex and touches no boundary row; only n . c_B = 0 breaks.
    exterior = {(int(e), int(f)) for e, f in plate["facets"]}
    normal_jump = None
    for e in range(len(tets)):
        for q, (a, b) in enumerate(EDGES):
            through = [f for f in range(4) if a in FACET_VERTICES[f] and b in FACET_VERTICES[f]]
            if all((e, f) not in exterior for f in through):
                normal_jump = {"kind": "edge_normal_jump", "element": e, "edge": q, "amount": 1e-3}
                break
        if normal_jump:
            break

    # The von Mises cube: the side y > x of the interior plane x = y translates rigidly along +y, so every coefficient
    # on that plane opens; no element of that side has a facet on x = 1 or y = 0, and +y is tangent to x = 0 and z = 0.
    t_arrays = documents[f"{T}.problem"][1]
    t_centroids = t_arrays["nodes"][t_arrays["tets"][:, :4]].mean(axis=1)
    t_side = [int(e) for e in np.nonzero(t_centroids[:, 1] > t_centroids[:, 0])[0]]

    # The interior body diagonal of the Drucker-Prager cube: no boundary row touches its mid node.
    dp_problem, dp_arrays = documents[f"{D}.problem"]
    dp_velocity = documents[DU][1]["velocity"]
    diagonal = int(np.argmin(np.linalg.norm(dp_arrays["nodes"] - 0.5, axis=1)))
    strong = _first_violating_direction(dp_problem, dp_arrays, dp_velocity, diagonal, 1e-3)

    h_problem, h_arrays = documents[f"{H}.problem"]
    xmin = [item["name"] for item in h_problem["boundary_sets"]].index("xmin")
    row = int(np.nonzero(h_arrays["facet_set"] == xmin)[0][0])
    wall_e, wall_f = (int(v) for v in h_arrays["facets"][row])
    wall_node = FACET_NODES[wall_f][0]

    a_arrays = documents[f"{A}.problem"][1]
    a_centroids = a_arrays["nodes"][a_arrays["tets"][:, :4]].mean(axis=1)
    side = [int(e) for e in np.nonzero(a_centroids[:, 1] > a_centroids[:, 0])[0]]
    opening = [-1.0 / S2, 1.0 / S2, 0.0]
    a_problem = documents[f"{A}.problem"][0]
    omission = _apex_omission(a_problem, a_arrays, documents[AU][1]["element_velocity"], side, opening)
    z_full = sorted(omission["apex_zero"] + [omission["omitted"]])

    entries = [
        ("perturbed_stress_component", "one perturbed stress component", P, L, "checker",
         {"kind": "perturb_sigma", "element": 0, "vertex": 0, "component": 0, "amount": 1e-3},
         "projection_move_too_large"),
        ("stress_scaled_1p01", "stress scaled by 1.01", P, L, "checker",
         {"kind": "scale_lower", "factor": 1.01}, "cone_violation_too_large"),
        ("facet_traction_jump", "a facet traction jump", P, L, "checker",
         {"kind": "element_offset", "element": 0, "component": 0, "amount": 1e-2}, "projection_move_too_large"),
        ("wrong_lower_value", "a Lower value that is not the multiplier of its field", P, L, "checker",
         {"kind": "scale_value", "factor": 1.01}, "projection_move_too_large"),
        ("wrong_load_normalisation", "a wrong load normalisation", P, U, "checker",
         {"kind": "scale_value", "factor": 1.1}, "value_mismatch"),
        ("swapped_bound_labels", "swapped LB and UB labels", P, L, "schema",
         {"kind": "relabel_upper"}, "result_wrong_bound_arrays"),
        ("curved_tetrahedron", "a curved tetrahedron", P, None, "schema",
         {"kind": "move_mid_node", "element": 0, "local_node": 4, "amount": 0.01}, "element_curved"),
        ("mismatched_problem_hash", "a mismatched problem hash", P, U, "schema",
         {"kind": "problem_hash"}, "hash_mismatch"),
        ("unsupported_domain_claim", "an unsupported domain_claim", P, U, "schema",
         {"kind": "certificate_claim", "claim": "halfspace_by_zero_extension"}, "unsupported_domain_claim"),
        ("unlisted_side_face", "a side face outside every declared set under a zero-extension claim", H, HU, "schema",
         {"kind": "drop_truncation_set", "name": "xmin"}, "unsupported_domain_claim"),
        ("wrong_gravity_in_equilibrium", "wrong gravity in equilibrium", G, GL, "checker",
         {"kind": "zero_gravity_lower"}, "projection_move_too_large"),
        ("isochory_one_vertex", "tr(D) != 0 at one vertex of an isochoric Upper field", P, U, "checker",
         {"kind": "perturb_velocity", "node": centre, "component": 0, "amount": 1e-3}, "projection_move_too_large"),
        ("dp_flow_cone_one_vertex", "Drucker-Prager flow-cone violation at one vertex", D, DU, "checker",
         {"kind": "perturb_velocity", "node": diagonal, "component": strong, "amount": 1e-3},
         "flow_violation_too_large"),
        ("footprint_off_rigid", "a footprint node off the rigid motion", P, U, "checker",
         {"kind": "perturb_velocity", "node": footprint, "component": 1, "amount": 1e-3}, "projection_move_too_large"),
        ("support_velocity_nonzero", "non-zero support velocity on a box target", P, U, "checker",
         {"kind": "perturb_velocity", "node": base_node, "component": 0, "amount": 1e-3}, "projection_move_too_large"),
        ("wall_velocity_zero_extension", "nonzero wall velocity under a zero-extension claim", H, HU, "checker",
         {"kind": "perturb_element_velocity", "element": wall_e, "local_node": wall_node, "component": 0, "amount": 1e-3},
         "projection_move_too_large"),
        ("normal_jump_isochoric", "non-zero normal jump at one facet node under isochoric flow", P, J, "checker",
         normal_jump, "projection_move_too_large"),
        ("vm_jump_opening", "a von Mises jump that opens, admissible only if read as the cone with tau = 0", T, TJ,
         "checker", {"kind": "side_translation", "elements": t_side, "direction": [0.0, 1.0, 0.0], "amount": 1e-3},
         "projection_move_too_large"),
        ("restrained_footing_component", "a non-zero velocity on a restrained footing component", R, RU, "checker",
         {"kind": "restrained_shear", "footing": 0, "component": 2, "axis": 1, "amount": 1e-3},
         "projection_move_too_large"),
        ("jump_between_nodes", "a jump admissible at facet nodes but not across the facet", A, AU, "checker",
         {"kind": "side_opening", "elements": side, "normal": opening, "amount": 1e-3}, "flow_violation_too_large"),
        ("wrong_gravity_work", "wrong gravity work or wrong dissipation quadrature", G, GU, "checker",
         {"kind": "value_rule", "body": "lumped", "quadrature": "vertex"}, "value_mismatch"),
        ("wrong_dissipation_quadrature", "wrong gravity work or wrong dissipation quadrature", G, GU,
         "checker", {"kind": "value_rule", "body": "exact", "quadrature": "centroid"}, "value_mismatch"),
        ("nonpositive_reference_work", "a non-positive reference work rate", P, U, "checker",
         {"kind": "negate_velocity"}, "work_nonpositive"),
        ("dissipation_charge_as_repair", "a finite dissipation charge offered as repair", D, DU, "checker",
         {"kind": "dissipation_charge", "node": diagonal, "component": strong, "amount": 1e-8, "charge": 1.0},
         "no_admissible_correction"),
        ("equality_just_above_floor", "an equality residual just above the floor", P, L, "final_check",
         {"kind": "perturb_sigma", "element": 0, "vertex": vertex, "component": 0, "amount": amount,
          "tau_multiple": 3.0}, "equality_residual_above_floor"),
        # Amendment 2 (interface 19.2.1): the declared apex face Z.
        ("apex_zero_omits_a_condition", "a declared apex face that omits a needed condition", A, AU, "checker",
         {"kind": "apex_zero_omits", "elements": side, "direction": opening, "amount": APEX_FIXTURE_DELTA, **omission},
         "interior_field_invalid"),
        # M1: Z lists a condition violated at about 1e-3 D_scale; the Z rows must move the field (soundness).
        ("apex_zero_hides_violation", "a declared apex face over a violated condition", A, AU, "checker",
         {"kind": "apex_zero_hides", "elements": side, "direction": opening, "amount": APEX_FIXTURE_DELTA,
          "element": omission["element"], "local_node": omission["local_node"],
          "perturbation": [APEX_FIXTURE_HIDE if v else 0.0 for v in omission["perturbation"]], "apex_zero": z_full},
         "projection_move_too_large"),
        # M1, positive companion: the omission fixture with the full face certifies.
        ("apex_zero_full_face", "a declared apex face that covers every needed condition", A, AU, "certify",
         {"kind": "apex_zero_full", "elements": side, "direction": opening, "amount": APEX_FIXTURE_DELTA,
          "element": omission["element"], "local_node": omission["local_node"],
          "perturbation": omission["perturbation"], "apex_zero": z_full},
         "certified_bound"),
        ("apex_zero_on_lower", "a declared apex face on a Lower bound", P, L, "schema",
         {"kind": "add_apex_zero", "apex_zero": [[0, 0]]}, "apex_zero_lower"),
        ("apex_zero_outside_hash", "a declared apex face not covered by the result hash", A, AU, "schema",
         {"kind": "add_apex_zero_unhashed", "apex_zero": [[0, 0]]}, "hash_mismatch"),
    ]
    assert plate_doc["format"] == FORMAT_PROBLEM
    return {
        "format": "ladruno-limit-negative-fixtures", "format_version": FORMAT_VERSION,
        "fixtures": [
            {"id": item[0], "harness_item": item[1], "problem": item[2], "result": item[3], "stage": item[4],
             "mutation": item[5], "expect": item[6]}
            for item in entries
        ],
    }


def build_fixture_set():
    cases = parse_cases()
    documents = {}
    for name in ("vm_uniaxial_tension", "vm_plate_shear"):
        problem, arrays = build_case_problem(cases[name], 0)
        _add_exact_set(documents, name, problem, arrays)

    # A half-space target whose whole exterior off the ground is fixed, with an admissible compact field.
    name = "vm_ground_patch_halfspace"
    case = parse_case_text(GROUND_PATCH_CASE)[name]
    problem, arrays = build_case_problem(case, 0)
    problem["domain"] = {"target": "halfspace", "truncation_sets": ["xmin", "xmax", "ymin", "ymax", "base"],
                         "symmetry_sets": [], "ground": {"axis": "z", "coordinate": 0.0}, "exterior_load_free": True}
    _with_hash(problem, arrays)
    footing = np.zeros((0, 6))
    field = ground_patch_field(arrays)
    _store_problem(documents, name, problem, arrays)
    upper = _store_result(documents, name, problem, "upper", fixture_upper_value(problem, arrays, field, footing),
                          _float_arrays(element_velocity=field, footing_velocity=footing), JUMP_ID,
                          objective=fixture_upper_objective(problem, arrays, field, footing))
    certificate = _certificate(problem, arrays, upper)
    documents[f"{record_stem(name, JUMP_ID)}.certificate"] = (certificate, None)
    documents[f"{record_stem(name, JUMP_ID)}.certificate.method"] = (_method_report(certificate), None)

    # Drucker-Prager tension: the exact field, continuous and in the broken space.
    name = "dp_uniaxial_tension"
    problem, arrays = build_case_problem(cases[name], 0)
    alpha, k = problem["regions"][0]["alpha"], problem["regions"][0]["k"]
    a, c = 1.0 / S3 + alpha, -0.5 / S3 + alpha
    x = arrays["nodes"]
    velocity = np.column_stack([x[:, 0], c / a * x[:, 1], c / a * x[:, 2]])
    value = k / a
    _store_problem(documents, name, problem, arrays)
    _store_result(documents, name, problem, "upper", value, _float_arrays(velocity=velocity, footing_velocity=footing))
    _store_result(documents, name, problem, "upper", value,
                  _float_arrays(element_velocity=velocity[arrays["tets"]], footing_velocity=footing), JUMP_ID)

    # Free cube at the apex: v = x is strictly inside the flow cone and no equality row exists.
    name = "dp_free_hydrostatic"
    problem, arrays = build_case_problem(parse_case_text(FREE_APEX_CASE)[name], 0)
    field = arrays["nodes"][arrays["tets"]]
    _store_problem(documents, name, problem, arrays)
    _store_result(documents, name, problem, "upper", fixture_upper_value(problem, arrays, field, footing),
                  _float_arrays(element_velocity=field, footing_velocity=footing), JUMP_ID,
                  objective=fixture_upper_objective(problem, arrays, field, footing))

    # Dead gravity: the exact Lower field, and an admissible quadratic Upper field whose work and
    # dissipation are sensitive to the quadrature rule.
    name = "vm_dead_gravity_shear_layer"
    problem, arrays = build_case_problem(cases[name], 0)
    gamma, height, target = 0.3, 1.0, cases[name]["target"]
    sigma = np.zeros((len(arrays["tets"]), 4, 6))
    sigma[..., 5] = target + gamma * (height - arrays["nodes"][arrays["tets"][:, :4], 2])
    z = arrays["nodes"][:, 2]
    velocity = np.column_stack([z * z - 0.5 * z, np.zeros_like(z), np.zeros_like(z)])
    _store_problem(documents, name, problem, arrays)
    _store_result(documents, name, problem, "lower", target, _float_arrays(sigma=sigma, resultant=footing))
    _store_result(documents, name, problem, "upper", fixture_upper_value(problem, arrays, velocity, footing),
                  _float_arrays(velocity=velocity, footing_velocity=footing),
                  objective=fixture_upper_objective(problem, arrays, velocity, footing))

    # Plate shear with the plate's uz restrained: the exact field has V_uz = 0 (amendment 1, item 9).
    name = "vm_plate_shear_restrained"
    case = copy.deepcopy(cases["vm_plate_shear"])
    case.update(name=name, restrained=[("plate", "uz")])
    problem, arrays = build_case_problem(case, 0)
    value, _, _, velocity, plate = _exact_fields("vm_plate_shear", arrays)
    _store_problem(documents, name, problem, arrays)
    _store_result(documents, name, problem, "upper", value, _float_arrays(velocity=velocity, footing_velocity=plate))

    documents["negative_fixtures"] = (build_negative_fixtures(documents), None)
    return documents


def build_problem_hashes():
    """problem_sha256 of every active case on every listed mesh (interface 15.1)."""
    table = {}
    for name, case in parse_cases().items():
        if case["status"] == "active" and case["boundaries"]:
            table[name] = {f"m{i}": build_case_problem(case, i)[0]["problem_sha256"] for i in range(len(case["mesh"]))}
    return table


def write_fixtures():
    DATA.mkdir(parents=True, exist_ok=True)
    for path in DATA.iterdir():
        if not path.name.endswith(".mesh.npz"):
            path.unlink()
    for stem, (doc, arrays) in build_fixture_set().items():
        (DATA / f"{stem}.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        if arrays is not None:
            np.savez(DATA / doc["arrays"], **arrays)
    (DATA / "problem_hashes.json").write_text(json.dumps(build_problem_hashes(), indent=0, sort_keys=True) + "\n",
                                              encoding="utf-8")


# ---------------------------------------------------------------------------------------------
# Tests: conventions
# ---------------------------------------------------------------------------------------------

def random_tet(rng):
    while True:
        x = rng.normal(size=(4, 3))
        if np.linalg.det(x[1:] - x[0]) > 0.2:
            return x


def tet10_coordinates(x):
    return np.vstack([x, [0.5 * (x[a] + x[b]) for a, b in EDGES]])


def gauss_triangle(order=8):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes, weights = 0.5 * (nodes + 1.0), 0.5 * weights
    for u, wu in zip(nodes, weights):
        for v, wv in zip(nodes, weights):
            l1, l2 = u, v * (1.0 - u)
            yield np.array([1.0 - l1 - l2, l1, l2]), 2.0 * wu * wv * (1.0 - u)


def gauss_tet(order=8):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes, weights = 0.5 * (nodes + 1.0), 0.5 * weights
    for (a, wa), (b, wb), (c, wc) in itertools.product(zip(nodes, weights), repeat=3):
        l1 = a
        l2 = b * (1.0 - a)
        l3 = c * (1.0 - a) * (1.0 - b)
        yield np.array([1.0 - l1 - l2 - l3, l1, l2, l3]), 6.0 * wa * wb * wc * (1.0 - a) ** 2 * (1.0 - b)


def test_facet_tables_are_outward_and_consistent():
    rng = np.random.default_rng(1)
    for _ in range(20):
        x = random_tet(rng)
        centroid = x.mean(axis=0)
        for f in range(4):
            normal, area = facet_geometry(x, f)
            assert area > 0.0
            assert normal @ (x[FACET_VERTICES[f][0]] - centroid) > 0.0
            assert f not in FACET_VERTICES[f]
    for f, nodes in enumerate(FACET_NODES):
        a, b, c = FACET_VERTICES[f]
        assert nodes[:3] == (a, b, c)
        for position, (p, q) in zip(nodes[3:], ((a, b), (b, c), (c, a))):
            assert set(EDGES[position - 4]) == {p, q}


def test_gradient_routes_agree():
    rng = np.random.default_rng(2)
    x = random_tet(rng)
    grads, volume = tet_gradients(x)
    for i in range(4):
        normal, area = facet_geometry(x, i)
        np.testing.assert_allclose(grads[i], -area / (3.0 * volume) * normal, atol=1e-12)


def test_p2_shape_functions_and_vertex_gradients():
    rng = np.random.default_rng(3)
    x = random_tet(rng)
    grads, _ = tet_gradients(x)
    nodal = np.vstack([np.eye(4), [0.5 * (np.eye(4)[a] + np.eye(4)[b]) for a, b in EDGES]])
    np.testing.assert_allclose([p2_values(b) for b in nodal], np.eye(10), atol=1e-14)
    step = 1e-6
    for m in range(4):
        analytic = p2_gradients_at_vertex(grads, m)
        for axis in range(3):
            dl = grads[:, axis] * step
            bary = np.eye(4)[m]
            numeric = (p2_values(bary + dl) - p2_values(bary - dl)) / (2.0 * step)
            np.testing.assert_allclose(analytic[:, axis], numeric, atol=1e-8)


def test_voigt_operator_gives_traction_divergence_and_strain_rate():
    rng = np.random.default_rng(4)
    sigma = rng.normal(size=6)
    n = rng.normal(size=3)
    tensor = np.array([[sigma[0], sigma[3], sigma[5]], [sigma[3], sigma[1], sigma[4]], [sigma[5], sigma[4], sigma[2]]])
    np.testing.assert_allclose(voigt_n(n) @ sigma, tensor @ n)
    x = random_tet(rng)
    grads, _ = tet_gradients(x)
    vertex_stress = rng.normal(size=(4, 6))
    divergence = sum(voigt_n(grads[i]) @ vertex_stress[i] for i in range(4))
    gradient = grads.T @ vertex_stress  # d sigma_c / d x_j, shape (3, 6)
    expected = [
        gradient[0, 0] + gradient[1, 3] + gradient[2, 5],
        gradient[0, 3] + gradient[1, 1] + gradient[2, 4],
        gradient[0, 5] + gradient[1, 4] + gradient[2, 2],
    ]
    np.testing.assert_allclose(divergence, expected)
    # A linear velocity gives its symmetric gradient, with engineering shear, at every vertex.
    grad_v = rng.normal(size=(3, 3))
    nodes = tet10_coordinates(x)
    velocity = nodes @ grad_v.T
    expected_d = [grad_v[0, 0], grad_v[1, 1], grad_v[2, 2], grad_v[0, 1] + grad_v[1, 0],
                  grad_v[1, 2] + grad_v[2, 1], grad_v[2, 0] + grad_v[0, 2]]
    for vertex in range(4):
        np.testing.assert_allclose(strain_rate_at_vertex(nodes, velocity, vertex), expected_d, atol=1e-12)


def test_rigid_velocity_has_zero_strain_rate():
    rng = np.random.default_rng(5)
    x = random_tet(rng)
    nodes = tet10_coordinates(x)
    u, w = rng.normal(size=3), rng.normal(size=3)
    velocity = u + np.cross(w, nodes)
    for vertex in range(4):
        np.testing.assert_allclose(strain_rate_at_vertex(nodes, velocity, vertex), 0.0, atol=1e-12)


def test_tangent_basis_rule():
    rng = np.random.default_rng(6)
    for _ in range(20):
        n = rng.normal(size=3)
        n /= np.linalg.norm(n)
        t1, t2 = tangent_basis(n)
        np.testing.assert_allclose([t1 @ n, t2 @ n, t1 @ t2], 0.0, atol=1e-14)
        np.testing.assert_allclose(np.cross(t1, t2), n, atol=1e-14)
    np.testing.assert_allclose(tangent_basis([0.0, 0.0, 1.0])[0], [0.0, -1.0, 0.0])


def test_cone_map_matches_invariants_and_principal_route():
    rng = np.random.default_rng(7)
    np.testing.assert_allclose(Q @ ONE, 0.0, atol=1e-15)
    for _ in range(200):
        sigma = rng.normal(size=6)
        k, alpha = rng.uniform(0.1, 2.0), rng.uniform(0.0, 0.5)
        i1, j2 = invariants(sigma)
        assert np.linalg.norm(Q @ sigma) ** 2 == pytest.approx(2.0 * j2)
        margin = S2 * (k - alpha * i1) - np.linalg.norm(Q @ sigma)
        assert margin == pytest.approx(-S2 * yield_value(sigma, k, alpha))
        assert yield_from_principal(sigma, k, alpha) == pytest.approx(yield_value(sigma, k, alpha))


def test_support_function_is_the_supremum_and_matches_conic_form():
    rng = np.random.default_rng(8)
    for alpha in (0.0, 0.2):
        k = 1.3
        for _ in range(100):
            z = rng.normal(size=5)
            r = np.linalg.norm(z) * (1.0 if alpha == 0.0 else rng.uniform(1.0, 2.0))
            d = S2 * alpha * r * ONE + Q.T @ z
            value = support_function(d, k, alpha)
            assert value == pytest.approx(S2 * k * r)
            theta, enorm = strain_split(d)
            if alpha > 0.0:
                maximiser = k / (3.0 * alpha) * ONE
            else:
                dev = np.concatenate([d[:3] - theta / 3.0, 0.5 * d[3:]])
                maximiser = S2 * k * dev / enorm + rng.normal() * ONE
            assert yield_value(maximiser, k, alpha) <= 1e-9
            assert maximiser @ d == pytest.approx(value)
            for _ in range(5):
                trial = rng.normal(size=6)
                g = yield_value(trial, k, alpha)
                if g > 0.0:
                    trial = trial * k / (g + k) if alpha == 0.0 else trial * 0.0
                assert trial @ d <= value + 1e-9
    assert support_function(np.array([1.0, 0, 0, 0, 0, 0]), 1.0, 0.0) == math.inf
    assert support_function(np.array([-1.0, -1, -1, 0, 0, 0]), 1.0, 0.2) == math.inf


def test_plane_strain_sections():
    rng = np.random.default_rng(9)
    for phi in (0.0, 10.0, 20.0, 30.0):
        c = 1.0
        k, alpha = plane_strain_match(c, phi)
        sin_phi, cos_phi = math.sin(math.radians(phi)), math.cos(math.radians(phi))
        for _ in range(20):
            m = rng.uniform(-3.0, 0.5)
            radius = c * cos_phi - m * sin_phi
            if radius <= 0.0:
                continue
            deltas = np.linspace(-10.0, 10.0, 400001)
            values = np.sqrt(radius ** 2 + deltas ** 2 / 3.0) + alpha * (3.0 * m + deltas) - k
            assert values.min() == pytest.approx(0.0, abs=1e-7)


def test_facet_work_table_and_volume_weights_are_exact():
    table = np.zeros((3, 6))
    area = 1.0
    for i in range(3):
        for j in range(3):
            table[i, j] = area / 30.0 if i == j else -area / 60.0
        for q, (a, b) in enumerate(((0, 1), (1, 2), (2, 0))):
            table[i, 3 + q] = 2.0 * area / 15.0 if i in (a, b) else area / 15.0
    numeric = np.zeros((3, 6))
    for bary, weight in gauss_triangle():
        shape = [bary[j] * (2.0 * bary[j] - 1.0) for j in range(3)]
        shape += [4.0 * bary[a] * bary[b] for a, b in ((0, 1), (1, 2), (2, 0))]
        numeric += weight * np.outer(bary, shape)
    np.testing.assert_allclose(table, numeric, atol=1e-13)
    weights = np.zeros(10)
    for bary, weight in gauss_tet():
        weights += weight * p2_values(bary)
    np.testing.assert_allclose(weights, [-0.05] * 4 + [0.2] * 6, atol=1e-13)


def test_footing_resultant_formula_and_power_are_exact():
    rng = np.random.default_rng(10)
    x = rng.normal(size=(3, 3))
    x_ref = rng.normal(size=3)
    traction = rng.normal(size=(3, 3))
    area = 0.5 * np.linalg.norm(np.cross(x[1] - x[0], x[2] - x[0]))
    force = sum(area / 3.0 * traction[a] for a in range(3))
    moment = sum(
        area * (1.0 / 6.0 if a == b else 1.0 / 12.0) * np.cross(x[a] - x_ref, traction[b])
        for a in range(3) for b in range(3)
    )
    numeric_force, numeric_moment = np.zeros(3), np.zeros(3)
    u, w = rng.normal(size=3), rng.normal(size=3)
    power = 0.0
    for bary, weight in gauss_triangle():
        point, t = bary @ x, bary @ traction
        numeric_force += weight * area * t
        numeric_moment += weight * area * np.cross(point - x_ref, t)
        power += weight * area * t @ (u + np.cross(w, point - x_ref))
    np.testing.assert_allclose(force, numeric_force)
    np.testing.assert_allclose(moment, numeric_moment)
    assert power == pytest.approx(force @ u + moment @ w)


def test_bernstein_coefficients_bound_a_quadratic_jump_over_the_facet():
    """Interface 12: vertex values do not cover a quadratic jump; Bernstein coefficients do."""
    def field(bary, nodal):
        shape = [bary[j] * (2.0 * bary[j] - 1.0) for j in range(3)]
        shape += [4.0 * bary[a] * bary[b] for a, b in ((0, 1), (1, 2), (2, 0))]
        return np.asarray(shape) @ nodal

    # The normal opening vanishes at the vertices and is negative at the edge midpoints.
    n = np.array([0.0, 0.0, 1.0])
    nodal = np.zeros((6, 3))
    nodal[3:, 2] = -1.0
    assert all(nodal[i] @ n >= 0.0 for i in range(3))
    centre = field(np.array([1.0, 1.0, 1.0]) / 3.0, nodal)
    assert centre @ n < 0.0  # a facet-vertex check would have passed
    pairs = ((0, 1), (1, 2), (2, 0))
    coefficients = [nodal[i] for i in range(3)]
    coefficients += [2.0 * nodal[3 + q] - 0.5 * (nodal[a] + nodal[b]) for q, (a, b) in enumerate(pairs)]
    assert min(c @ n for c in coefficients) < 0.0  # the Bernstein check refuses it
    rng = np.random.default_rng(11)
    nodal = rng.normal(size=(6, 3))
    coefficients = [nodal[i] for i in range(3)]
    coefficients += [2.0 * nodal[3 + q] - 0.5 * (nodal[a] + nodal[b]) for q, (a, b) in enumerate(pairs)]
    total = 0.0
    for bary, weight in gauss_triangle():
        basis = [bary[j] ** 2 for j in range(3)] + [2.0 * bary[a] * bary[b] for a, b in pairs]
        assert min(basis) >= 0.0 and sum(basis) == pytest.approx(1.0)
        np.testing.assert_allclose(sum(b * c for b, c in zip(basis, coefficients)), field(bary, nodal))
        total += weight * basis[0]
    assert total == pytest.approx(1.0 / 6.0)


def test_jump_support_uses_the_same_cone_map():
    rng = np.random.default_rng(12)
    n = rng.normal(size=3)
    n /= np.linalg.norm(n)
    t1, _ = tangent_basis(n)
    assert support_function(voigt_n(n).T @ (2.5 * t1), 1.2, 0.0, tol=1e-10) == pytest.approx(1.2 * 2.5)
    alpha = 0.2
    slip = 1.0
    opening = 3.0 * alpha * slip / math.sqrt(1.0 - 12.0 * alpha * alpha)
    jump = slip * t1 + opening * n
    d = voigt_n(n).T @ jump
    theta, enorm = strain_split(d)
    assert theta == pytest.approx(opening)
    assert theta == pytest.approx(3.0 * S2 * alpha * enorm)


def test_jump_constants_match_the_coulomb_form():
    for phi in (0.0, 10.0, 20.0, 30.0):
        c = 1.7
        k, alpha = plane_strain_match(c, phi)
        tau, root = jump_constants(alpha)
        assert tau == pytest.approx(math.tan(math.radians(phi)), abs=1e-14)
        assert k / root == pytest.approx(c), "k_J is the cohesion under the plane-strain match"
        if phi > 0.0:
            assert k / (3.0 * alpha) * tau == pytest.approx(c), "slip with o = s tan(phi) dissipates c s"
    with pytest.raises(ContractError):
        jump_constants(0.3)


def test_bernstein_map_reproduces_the_quadratic_trace():
    rng = np.random.default_rng(13)
    nodal = rng.normal(size=(6, 3))
    coefficients = bernstein_from_nodal(nodal)
    pairs = ((0, 1), (1, 2), (2, 0))
    for bary, _ in gauss_triangle(3):
        shape = [bary[j] * (2.0 * bary[j] - 1.0) for j in range(3)]
        shape += [4.0 * bary[a] * bary[b] for a, b in pairs]
        basis = [bary[j] ** 2 for j in range(3)] + [2.0 * bary[a] * bary[b] for a, b in pairs]
        np.testing.assert_allclose(np.array(shape) @ nodal, np.array(basis) @ coefficients)


def test_equality_level_carries_no_conditioning_factor():
    assert tau_eq() == 64.0 * 2.0 ** -53
    assert TOLERANCES["TOL_KAPPA_CAP"] * tau_eq() <= TOLERANCES["TOL_FLOOR_MAX"]


def test_row_backward_error_accepts_a_homogeneous_row_of_negligible_entries():
    # Componentwise, this row has residual equal to its denominator (omega = 1). Row-normwise it is rounding.
    matrix = np.array([[1.0, -1.0], [1.0, 1.0]])
    rhs = np.array([0.0, 2.0])
    x = np.array([1.0, 1.0])
    assert row_backward_error(matrix, rhs, x, np.ones(2)) == 0.0
    tiny = np.array([1e-30, 3e-30])
    row = matrix[:1]
    componentwise = abs(row @ tiny)[0] / (np.abs(row) @ np.abs(tiny))[0]
    assert componentwise == pytest.approx(0.5)
    assert row_backward_error(row, [0.0], tiny, np.ones(2)) < tau_eq()
    assert row_backward_error(row, [0.0], [1.0, 1.0 + 1e-9], np.ones(2)) > tau_eq()
    with pytest.raises(ContractError):
        row_backward_error([[1.0]], [0.0], [1.0], [0.0])


def test_level_ground_von_mises_superposition_is_blind_to_gravity():
    rng = np.random.default_rng(14)
    gamma = 1.7
    for _ in range(20):
        sigma = rng.normal(size=6)
        z = rng.uniform(-2.0, 0.0)
        shifted = sigma + gamma * z * ONE
        assert yield_value(shifted, 1.0) == pytest.approx(yield_value(sigma, 1.0))
    grads, _ = tet_gradients(random_tet(rng))
    x = random_tet(rng)
    grads, _ = tet_gradients(x)
    vertex_shift = np.array([gamma * point[2] * ONE for point in x])
    divergence = sum(voigt_n(grads[i]) @ vertex_shift[i] for i in range(4))
    np.testing.assert_allclose(divergence + np.array([0.0, 0.0, -gamma]), 0.0, atol=1e-12)
    for normal in np.eye(3):
        for tangent in tangent_basis(normal):
            assert tangent @ (voigt_n(normal) @ ONE) == pytest.approx(0.0)


# ---------------------------------------------------------------------------------------------
# Tests: hash and meshes
# ---------------------------------------------------------------------------------------------

def test_hash_treats_negative_zero_as_zero_and_has_the_published_digest():
    arrays = {"v": np.array([0.5, 0.0])}
    doc = {"b": [1, None, True], "a": "x"}
    digest = document_hash(doc, arrays, "h")
    assert digest == document_hash({"a": "x", "b": [1.0, None, True]}, {"v": np.array([0.5, -0.0])}, "h")
    assert document_hash({"z": -0.0}, {}, "h") == document_hash({"z": 0.0}, {}, "h")
    assert re.search(r"Worked digest.*?hashes to `([0-9a-f]{64})`", interface_text(), re.S).group(1) == digest


@pytest.mark.parametrize("text", ['{"a": 1e999}', '{"a": -1e999}', '{"a": NaN}', '{"a": Infinity}'])
def test_non_finite_json_is_refused(text):
    with pytest.raises(ContractError) as error:
        strict_json_loads(text)
    assert error.value.code == "json_non_finite"


def test_canonical_encoding_is_stable():
    assert canonical_encode(1) == canonical_encode(1.0)
    assert canonical_encode({"b": 1, "a": [None, True]}) == canonical_encode({"a": [None, True], "b": 1.0})
    assert canonical_encode("a") == b"s" + struct.pack("<Q", 1) + b"a"
    doc = {"x": 1, "h": "ignored"}
    arrays = {"a": np.array([1.0], dtype="<f8")}
    assert document_hash(doc, arrays, "h") == document_hash({"x": 1.0, "h": "other"}, arrays, "h")
    assert document_hash(doc, arrays, "h") != document_hash(doc, {"a": np.array([2.0], dtype="<f8")}, "h")
    with pytest.raises(ContractError):
        strict_json_loads('{"a": NaN}')


def test_kuhn_mesh_rule():
    vertices, tets = kuhn_tet4((0, 0, 0), (1, 2, 3), (2, 1, 1))
    assert len(vertices) == 3 * 2 * 2 and len(tets) == 12
    assert all(np.linalg.det(vertices[t[1:]] - vertices[t[0]]) > 0.0 for t in tets)
    volume = sum(np.linalg.det(vertices[t[1:]] - vertices[t[0]]) / 6.0 for t in tets)
    assert volume == pytest.approx(6.0)
    nodes, tet10 = tet10_from_tet4(*kuhn_tet4((0, 0, 0), (1, 1, 1), (1, 1, 1)))
    assert nodes.shape == (27, 3) and tet10.shape == (6, 10)


def test_doubling_every_count_nests_the_kuhn_mesh():
    coarse_v, coarse_t = kuhn_tet4((0, 0, 0), (1, 2, 1), (1, 2, 1))
    fine_v, fine_t = kuhn_tet4((0, 0, 0), (1, 2, 1), (2, 4, 2))
    inverse = [np.linalg.inv(coarse_v[t[1:]].T - coarse_v[t[0]][:, None]) for t in coarse_t]

    def inside(points, index):
        lam = inverse[index] @ (points - coarse_v[coarse_t[index][0]]).T
        return np.all(lam > -1e-12) and np.all(lam.sum(axis=0) < 1.0 + 1e-12)

    for tet in fine_t:
        assert any(inside(fine_v[tet], index) for index in range(len(coarse_t)))


def test_uniform_refinement_rule():
    rng = np.random.default_rng(15)
    x = random_tet(rng)
    points, children = refine_tet4(x, np.array([[0, 1, 2, 3]]))
    assert len(points) == 10 and len(children) == 8
    _, parent = tet_gradients(x)
    volumes = [tet_gradients(points[c])[1] for c in children]
    assert min(volumes) > 0.0
    np.testing.assert_allclose(volumes, parent / 8.0)
    vertices, tets = kuhn_tet4((0, 0, 0), (1, 1, 1), (1, 1, 1))
    fine_v, fine_t = refine_tet4(vertices, tets)
    interior, exterior = mesh_facets(fine_t)
    assert len(fine_t) == 48 and len(exterior) == 4 * len(mesh_facets(tets)[1])
    assert 2 * len(interior) + len(exterior) == 4 * len(fine_t)


# ---------------------------------------------------------------------------------------------
# Tests: case file
# ---------------------------------------------------------------------------------------------

ALL_CASES = parse_cases()
ACTIVE = sorted(name for name, case in ALL_CASES.items() if case["status"] == "active")
GATE_3D = [name for name in ACTIVE if ALL_CASES[name]["class"] == "upper_gate_3d"]
OPERATOR = [name for name in ACTIVE if ALL_CASES[name]["class"] == "operator_work"]
SANDWICH = [name for name in ACTIVE if ALL_CASES[name]["class"] == "sandwich"]
REFERENCE = [name for name in ACTIVE if name not in GATE_3D + OPERATOR + SANDWICH]
STRIPS = [name for name in REFERENCE if ALL_CASES[name]["class"] == "one_sided_reference"]
TODO = sorted(name for name, case in ALL_CASES.items() if case["status"] == "todo")
GRAVITY_KILLERS = ["vm_dead_gravity_shear_layer", "dp_weighted_apex_column", "vm_weighted_column_sandwich"]


def assertion_params(tokens):
    return tokens[1], parse_key_values(tokens[2:])


def first_layer(case, mesh_index, params):
    axis = int(params["axis"])
    return (case["box"][3 + axis] - case["box"][axis]) / case["mesh"][mesh_index][1][axis]


def test_case_inventory_covers_the_required_set():
    assert {"tresca_strip_vm_thin_slab", "dp_prandtl_cohesion_thin_slab", "dp_prandtl_surcharge_thin_slab"} <= set(STRIPS)
    assert {"gravity_work_rigid_motion"} == set(OPERATOR)
    assert set(GRAVITY_KILLERS) <= set(ACTIVE)
    assert {ALL_CASES[name]["phi_deg"] for name in GATE_3D} == {0.0, 10.0, 20.0, 30.0}
    assert len(GATE_3D) == 8
    assert "strip_self_weight_ngamma" in TODO


@pytest.mark.parametrize("name", REFERENCE)
def test_reference_case_values_and_assertions(name):
    case = ALL_CASES[name]
    assert case["source"] and case["checked"] and case["units"]
    if name in ("tresca_strip_vm_thin_slab", "tresca_strip_weighted_insensitive"):
        # Amendment 2, item B (OD3): continuous P2 is a control on the Tresca strips; the jump pair binds.
        assert set(case["formulation"]) == {LOWER_ID, JUMP_ID} and case["control"] == [CONTROL_ID]
    else:
        assert set(case["formulation"]) == {LOWER_ID, CONTROL_ID, JUMP_ID}
    assert case["target"] == pytest.approx(TARGET_FORMULAS[case["formula"]](case["formula_params"]), rel=1e-14)
    kinds = {tokens[0] for tokens in case["assert"]}
    assert {"lower_le_target", "upper_ge_target"} <= kinds, "both one-sided assertions are required"
    for tokens in case["assert"]:
        if tokens[0] in ("lower_le_target", "upper_ge_target", "lower_ge_target"):
            assert float(tokens[1]) == TOLERANCES["TOL_CASE_ONE_SIDED"]
    if case["class"] == "exact_on_mesh":
        assert ["gap_le", "1e-6"] in case["assert"]
        assert float("1e-6") == TOLERANCES["TOL_CASE_EXACT_GAP"]
    if case["class"] == "exact_lower_on_mesh":
        assert "lower_ge_target" in kinds and "upper_le_formula" in kinds
    if case["class"] == "one_sided_reference":
        assert "lower_ge_target" not in kinds
        limits = [tokens[1] for tokens in case["gate"] if tokens[0] == "gap_le"]
        assert len(limits) == 1 and limits[0] in TOLERANCES


@pytest.mark.parametrize("name", STRIPS)
def test_strip_cases_satisfy_the_lemma_conditions(name):
    case = ALL_CASES[name]
    phi = case.get("phi_deg", 0.0)
    width = 2.0
    x_max, depth = prandtl_extent(0.5 * width, phi)
    x0, y0, z0, x1, y1, z1 = case["box"]
    assert x0 == 0.0 and y0 == 0.0 and z1 == 0.0
    assert x1 >= max(1.5 * width, x_max), "lateral extent from the centreline"
    assert -z0 >= max(width / S2, depth), "depth"
    kinds = {item["name"]: item["kind"] for item in case["boundaries"]}
    assert kinds["far"] == "fixed" and kinds["base"] == "fixed", "fixed x and z walls"
    assert kinds["y0"] == kinds["y1"] == kinds["axis"] == "symmetry"
    assert ("strip", "ux") in case["restrained"] and ("strip", "ry") in case["restrained"]
    region = case["regions"][0]
    if region["criterion"] == "drucker_prager":
        c = 0.0 if region["k"] == 0.0 else 1.0
        k, alpha = plane_strain_match(c, phi)
        assert region["alpha"] == pytest.approx(alpha, rel=1e-14)
        assert region["k"] == pytest.approx(k, rel=1e-14, abs=1e-14)
    gates = {tokens[0] for tokens in case["gate"]}
    assert {"certified_bracket", "gap_le", "upper_nonincrease", "lower_nondecrease"} <= gates


def gate_mesh_stems(phi):
    prefix = "upper3d" if phi == 0.0 else f"upper3d_phi{int(phi):02d}"
    return [f"{prefix}_{letter}" for letter in "abcd"]


@pytest.mark.parametrize("name", GATE_3D)
def test_true_3d_gate_cases_are_complete(name):
    case = ALL_CASES[name]
    phi = case["phi_deg"]
    assert case["target"] is None
    assert JUMP_ID in case["formulation"] and case["control"] == [CONTROL_ID]
    gates = [tokens[0] for tokens in case["gate"]]
    for required in ("kinematic_status", "certified_upper", "certified_bracket", "upper_nonincrease",
                     "lower_nondecrease"):
        assert required in gates, "feasibility alone must not pass"
    assert "gap_le" not in gates, "amendment 3 (OD5): the 3D gap limits are P4 targets, not P0 gates"
    assert [int(t[2]) for t in case["gate"] if t[0] == "certified_upper"] == list(range(5)), "amendment 3: meshes 0-4"
    assert case["record"] == [["gap", "mesh", str(i)] for i in range(6)], "the gap trend is recorded, mesh 5 included"
    stems = gate_mesh_stems(phi)
    assert case["mesh"] == [("file", f"{stem}.mesh.npz", 0) for stem in stems] + [
        ("file", f"{stems[0]}.mesh.npz", 1), ("file", f"{stems[3]}.mesh.npz", 1)]
    assert case["nested"] == [(0, 4), (3, 5)]
    region = case["regions"][0]
    if phi == 0.0:
        assert region["criterion"] == "von_mises"
    else:
        k, alpha = plane_strain_match(1.0, phi)
        assert (region["k"], region["alpha"]) == (pytest.approx(k, rel=1e-14), pytest.approx(alpha, rel=1e-14))
    wall = "fixed" if "fixed" in name else "symmetry"
    kinds = {item["name"]: item["kind"] for item in case["boundaries"]}
    assert all(kinds[key] == wall for key in ("xmin", "xmax", "ymin", "ymax", "base"))


AMENDMENT_3_GATE = ([["kinematic_status", "formulation_feasible"]]
                    + [["certified_upper", "mesh", str(i)] for i in range(5)]
                    + [["certified_bracket", "mesh", str(i)] for i in range(5)]
                    + [["upper_nonincrease", "0", "4"], ["lower_nondecrease", "0", "4"]])


def test_true_3d_gate_is_meshes_0_to_4_after_amendment_3():
    """Amendment 3 (owner decision OD5): the eight 3D cases gate exactly this set, and no gate line names mesh 5."""
    expected = sorted(f"upper3d_{wall}_phi{phi:02d}" for wall in ("fixed", "roller") for phi in (0, 10, 20, 30))
    assert sorted(GATE_3D) == expected
    for name in GATE_3D:
        case = ALL_CASES[name]
        assert [list(tokens) for tokens in case["gate"]] == AMENDMENT_3_GATE, name
        assert not any("5" in tokens[1:] for tokens in case["gate"]), f"{name}: mesh 5 gates nothing in P0"
        assert len(case["mesh"]) == 6 and (3, 5) in case["nested"], f"{name}: mesh 5 stays as evidence"
        assert ["gap", "mesh", "5"] in case["record"], name


def test_three_d_gap_limits_are_p4_targets_with_unchanged_values():
    """Amendment 3: the values stay (certificate tolerance blocks must not churn); no case gates on them in P0."""
    assert TOLERANCES["MAX_GAP_3D_VM"] == 0.30 and TOLERANCES["MAX_GAP_3D_DP"] == 0.50
    used = {tokens[1] for case in ALL_CASES.values() for tokens in case["gate"] if tokens[0] == "gap_le"}
    assert used == {"MAX_GAP_STRIP", "MAX_GAP_PRANDTL"}
    text = interface_text()
    table = text[text.index("### 17.8 Named tolerances"):text.index("### 17.9")]
    for name in ("MAX_GAP_3D_VM", "MAX_GAP_3D_DP"):
        row = next(line for line in table.splitlines() if line.startswith(f"| `{name}` |"))
        assert "P4 exit target, not a P0 gate" in row, name


@pytest.mark.parametrize("name", [n for n in GATE_3D if ALL_CASES[n]["phi_deg"] > 0.0])
def test_friction_gate_box_holds_the_mechanism_scale(name):
    case = ALL_CASES[name]
    x_extent, depth = prandtl_extent(0.5, case["phi_deg"])
    x0, y0, z0, x1, y1, z1 = case["box"]
    assert x1 == -x0 and y1 == -y0 and x1 == y1 and z1 == 0.0
    assert x1 >= x_extent and -z0 >= depth
    planes = {item["name"]: item["plane"] for item in case["boundaries"]}
    assert planes["xmax"] == (0, x1) and planes["ymin"] == (1, y0) and planes["base"] == (2, z0)


def committed_gate_meshes():
    return [(name, phi) for phi in (0.0, 10.0, 20.0, 30.0) for name in gate_mesh_stems(phi)]


@pytest.mark.parametrize("file_name,phi", committed_gate_meshes())
def test_gate_meshes_are_unstructured_with_a_conformal_footprint(file_name, phi):
    with np.load(DATA / f"{file_name}.mesh.npz", allow_pickle=False) as npz:
        assert set(npz.files) == {"nodes", "tets"}
        vertices, tets = npz["nodes"], npz["tets"]
    assert vertices.dtype.str == "<f8" and tets.dtype.str == "<i8" and tets.shape[1] == 4
    assert all(np.linalg.det(vertices[t[1:]] - vertices[t[0]]) > 0.0 for t in tets)
    case = next(ALL_CASES[n] for n in GATE_3D if ALL_CASES[n]["phi_deg"] == phi)
    lo, hi = np.array(case["box"][:3]), np.array(case["box"][3:])
    assert np.all(vertices >= lo) and np.all(vertices <= hi)
    interior = np.all((vertices > lo + 1e-9) & (vertices < hi - 1e-9), axis=1)
    assert interior.sum() >= 5
    assert len(np.unique(np.round(vertices[interior, 0], 9))) > 0.5 * interior.sum(), "not a lattice"
    valence = np.bincount(tets.ravel())
    assert len(np.unique(valence[interior])) >= 3, "interior vertex valence varies"
    index = [i for i, entry in enumerate(case["mesh"]) if entry[1].startswith(file_name + ".") and entry[2] == 0][0]
    doc, arrays = build_case_problem(case, index)
    validate_problem(doc, arrays)
    names = [item["name"] for item in doc["boundary_sets"]]
    rows = np.nonzero(arrays["facet_set"] == names.index("contact"))[0]
    area = sum(facet_geometry(arrays["nodes"][arrays["tets"][e][:4]], int(f))[1] for e, f in arrays["facets"][rows])
    assert area == pytest.approx(1.0, abs=1e-12), "the footprint is exactly the conforming unit square"


def test_gate_mesh_sizes_increase_within_each_box():
    for phi in (0.0, 10.0, 20.0, 30.0):
        counts = []
        for stem in gate_mesh_stems(phi):
            with np.load(DATA / f"{stem}.mesh.npz", allow_pickle=False) as npz:
                counts.append(len(npz["tets"]))
        assert counts == sorted(counts) and len(set(counts)) == 4


def test_gate_mesh_files_match_their_published_digests():
    table = dict(re.findall(r"^\| `(upper3d[a-z0-9_]*\.mesh\.npz)` \| `([0-9a-f]{64})` \|", interface_text(), re.M))
    files = sorted(path.name for path in DATA.glob("*.mesh.npz"))
    assert sorted(table) == files and len(files) == 16
    for name, digest in table.items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == digest, name


def test_gate_refined_mesh_is_a_valid_problem():
    case = ALL_CASES["upper3d_roller_phi20"]
    assert case["mesh"][4] == ("file", "upper3d_phi20_a.mesh.npz", 1)
    assert case["mesh"][5] == ("file", "upper3d_phi20_d.mesh.npz", 1)
    doc, arrays = build_case_problem(case, 4)
    validate_problem(doc, arrays)
    with np.load(DATA / "upper3d_phi20_a.mesh.npz", allow_pickle=False) as npz:
        assert len(arrays["tets"]) == 8 * len(npz["tets"])


def test_committed_problem_hashes_match_the_cases():
    table = load_problem_hashes()
    assert set(table) == set(ACTIVE) - set(OPERATOR)
    checked = 0
    for name in table:
        case = ALL_CASES[name]
        assert set(table[name]) == {f"m{i}" for i in range(len(case["mesh"]))}
        for i, entry in enumerate(case["mesh"]):
            if entry[0] == "kuhn" and 6 * np.prod(entry[1]) > 800 or entry[0] == "file" and (entry[2] or "_a" not in entry[1]):
                continue
            assert build_case_problem(case, i)[0]["problem_sha256"] == table[name][f"m{i}"], (name, i)
            checked += 1
    assert checked >= len(table)


@pytest.mark.parametrize("name", TODO)
def test_todo_cases_carry_no_value(name):
    case = ALL_CASES[name]
    assert "target" not in case and not case["assert"] and not case["gate"] and "TODO" in case["desc"]


@pytest.mark.parametrize("name", REFERENCE + SANDWICH + GATE_3D)
def test_case_builds_a_valid_problem(name):
    case = ALL_CASES[name]
    doc, arrays = build_case_problem(case, 0)
    validate_problem(doc, arrays)
    for coarse, fine in case["nested"]:
        first, second = case["mesh"][coarse], case["mesh"][fine]
        if first[0] == "kuhn":
            assert second[1] == tuple(2 * n for n in first[1])
        else:
            assert second[1] == first[1] and second[2] == first[2] + 1


def test_rigid_motion_gravity_work_is_exact_at_l0():
    case = ALL_CASES["gravity_work_rigid_motion"]
    probe, region = case["probe"], case["regions"][0]
    b0, u, w, about = (np.array(v) for v in (region["b0"], probe["u"], probe["w"], probe["about"]))
    lo, hi = np.array(case["box"][:3]), np.array(case["box"][3:])
    volume = float(np.prod(hi - lo))
    moment = volume * (0.5 * (lo + hi) - about)
    assert case["target"] == pytest.approx(b0 @ (u * volume + np.cross(w, moment)), rel=1e-14)
    tolerance = float(case["assert"][0][1])
    for index in range(len(case["mesh"])):
        nodes, tets = tet10_from_tet4(*case_tet4(case, index))
        velocity = u + np.cross(w, nodes - about)
        work = 0.0
        for tet in tets:
            _, v = tet_gradients(nodes[tet[:4]])
            work += v * (np.array([-0.05] * 4 + [0.2] * 6) @ (velocity[tet] @ b0))
        assert work == pytest.approx(case["target"], rel=tolerance)


def band_velocity(formula, nodes, h, height):
    v = np.zeros_like(nodes)
    if formula == "shear_band":
        v[:, 0] = np.minimum(nodes[:, 1] / h, 1.0)
    elif formula == "dead_shear_band":
        v[:, 0] = np.minimum(nodes[:, 2] / h, 1.0)
    else:
        v[:, 2] = -np.minimum((height - nodes[:, 2]) / h, 1.0)
    return v


@pytest.mark.parametrize("name", ["vm_body_force_shear_block", "vm_dead_gravity_shear_layer", "dp_weighted_apex_column"])
def test_band_upper_fields_attain_their_formula(name):
    case = ALL_CASES[name]
    tokens = next(t for t in case["assert"] if t[0] == "upper_le_formula")
    formula, params = assertion_params(tokens)
    for index in range(len(case["mesh"])):
        h = first_layer(case, index, params)
        limit = band_formula(formula, params, h)
        assert limit > case["target"], "the band field cannot reach the exact value on a finite mesh"
        doc, arrays = build_case_problem(case, index)
        velocity = band_velocity(formula, arrays["nodes"], h, params["H"])
        assert fixture_upper_value(doc, arrays, velocity, np.zeros((0, 6))) == pytest.approx(limit, rel=1e-12)


def test_weighted_exact_lower_fields_are_admissible():
    # Dead-load layer: sigma_zx = lambda + gamma (H - z). Apex column: sigma = (lambda + gamma z) I.
    layer, apex = ALL_CASES["vm_dead_gravity_shear_layer"], ALL_CASES["dp_weighted_apex_column"]
    for case, stress in (
        (layer, lambda z: (layer["target"] + 0.3 * (1.0 - z)) * np.array([0, 0, 0, 0, 0, 1.0])),
        (apex, lambda z: (apex["target"] + 0.2 * z) * ONE),
    ):
        region = case["regions"][0]
        doc, arrays = build_case_problem(case, 1)
        worst = -math.inf
        for tet in arrays["tets"]:
            x = arrays["nodes"][tet[:4]]
            grads, _ = tet_gradients(x)
            sigma = np.array([stress(point[2]) for point in x])
            residual = sum(voigt_n(grads[i]) @ sigma[i] for i in range(4)) + np.array(region["b0"])
            np.testing.assert_allclose(residual, 0.0, atol=1e-12)
            worst = max(worst, max(yield_from_principal(s, region["k"], region.get("alpha", 0.0)) for s in sigma))
        assert worst == pytest.approx(0.0, abs=1e-12), "the field reaches yield and no further"


def test_column_sandwich_has_explicit_fields_and_is_not_claimed_exact():
    case = ALL_CASES["vm_weighted_column_sandwich"]
    assert case["target"] is None and case["class"] == "sandwich" and "NOT exact" in case["desc"]
    lower = next(t for t in case["assert"] if t[0] == "lower_ge_formula")
    upper = next(t for t in case["assert"] if t[0] == "upper_le_formula")
    (lower_name, p), (upper_name, _) = assertion_params(lower), assertion_params(upper)
    low, high = TARGET_FORMULAS[lower_name](p), TARGET_FORMULAS[upper_name](p)
    assert low < high
    k, gamma, height = p["k"], p["gamma"], p["H"]
    doc, arrays = build_case_problem(case, 1)
    worst = -math.inf
    for tet in arrays["tets"]:
        x = arrays["nodes"][tet[:4]]
        grads, _ = tet_gradients(x)
        sigma = np.zeros((4, 6))
        sigma[:, 2] = -low - gamma * (height - x[:, 2])
        residual = sum(voigt_n(grads[i]) @ sigma[i] for i in range(4)) + np.array([0.0, 0.0, -gamma])
        np.testing.assert_allclose(residual, 0.0, atol=1e-12)
        worst = max(worst, max(yield_value(s, k) for s in sigma))
    assert worst == pytest.approx(0.0, abs=1e-12)
    x = arrays["nodes"]
    velocity = np.column_stack([x[:, 0] / (2.0 * height), x[:, 1] / (2.0 * height), -x[:, 2] / height])
    assert fixture_upper_value(doc, arrays, velocity, np.zeros((0, 6))) == pytest.approx(high, rel=1e-12)


def test_gravity_sensitivity_flags_and_the_zero_gravity_kill():
    flags = {name: ALL_CASES[name].get("gravity_sensitive") for name in ACTIVE}
    assert all(flags[name] == "yes" for name in GRAVITY_KILLERS)
    assert flags["tresca_strip_weighted_insensitive"] == "no"
    # With gravity zeroed the Upper value is the weightless one, which breaks each Upper limit.
    layer = ALL_CASES["vm_dead_gravity_shear_layer"]
    _, p = assertion_params(next(t for t in layer["assert"] if t[0] == "upper_le_formula"))
    assert p["k"] > band_formula("dead_shear_band", p, first_layer(layer, 0, p))
    apex = ALL_CASES["dp_weighted_apex_column"]
    _, p = assertion_params(next(t for t in apex["assert"] if t[0] == "upper_le_formula"))
    assert p["k"] / (3.0 * p["alpha"]) > band_formula("apex_band", p, first_layer(apex, 0, p))
    column = ALL_CASES["vm_weighted_column_sandwich"]
    _, p = assertion_params(next(t for t in column["assert"] if t[0] == "upper_le_formula"))
    assert S3 * p["k"] > TARGET_FORMULAS["column_uniform"](p)


def test_gate_limits_are_named_and_ordered():
    assert 0.0 < TOLERANCES["MAX_GAP_STRIP"] <= TOLERANCES["MAX_GAP_PRANDTL"]
    assert TOLERANCES["MAX_GAP_PRANDTL"] <= TOLERANCES["MAX_GAP_3D_VM"] <= TOLERANCES["MAX_GAP_3D_DP"]


# ---------------------------------------------------------------------------------------------
# Tests: fixtures
# ---------------------------------------------------------------------------------------------

FIXTURE_NAMES = ("vm_uniaxial_tension", "vm_plate_shear")


def load_set(name):
    stem = case_stem(name)
    problem, arrays = load_document(DATA / f"{stem}.problem.json")
    lower, lower_arrays = load_document(DATA / f"{record_stem(name, LOWER_ID)}.result.json")
    upper, upper_arrays = load_document(DATA / f"{record_stem(name, CONTROL_ID)}.result.json")
    lower_cert, _ = load_document(DATA / f"{record_stem(name, LOWER_ID)}.certificate.json")
    upper_cert, _ = load_document(DATA / f"{record_stem(name, CONTROL_ID)}.certificate.json")
    bracket, _ = load_document(DATA / f"{bracket_stem(name, LOWER_ID, CONTROL_ID)}.bracket.json")
    feasibility, _ = load_document(DATA / f"{record_stem(name, CONTROL_ID)}.feasibility.json")
    return {
        "problem": problem, "arrays": arrays, "lower": lower, "lower_arrays": lower_arrays,
        "upper": upper, "upper_arrays": upper_arrays, "lower_cert": lower_cert,
        "upper_cert": upper_cert, "bracket": bracket, "feasibility": feasibility,
    }


def test_committed_fixtures_match_their_generator():
    for stem, (doc, arrays) in build_fixture_set().items():
        stored, stored_arrays = load_document(DATA / f"{stem}.json")
        assert stored == json.loads(json.dumps(doc)), stem
        if arrays is not None:
            assert set(stored_arrays) == set(arrays), stem
            for name, array in arrays.items():
                assert stored_arrays[name].dtype == array.dtype
                np.testing.assert_array_equal(stored_arrays[name], array)


def test_fixture_files_are_tiny():
    for path in DATA.iterdir():
        assert path.stat().st_size < 64 * 1024, path.name


@pytest.mark.parametrize("name", FIXTURE_NAMES)
def test_fixture_set_validates(name):
    s = load_set(name)
    validate_problem(s["problem"], s["arrays"])
    validate_result(s["lower"], s["lower_arrays"], s["problem"], s["arrays"])
    validate_result(s["upper"], s["upper_arrays"], s["problem"], s["arrays"])
    jump, jump_arrays = load_document(DATA / f"{record_stem(name, JUMP_ID)}.result.json")
    validate_result(jump, jump_arrays, s["problem"], s["arrays"])
    assert jump["formulation"] == JUMP_ID and "velocity" not in jump_arrays
    validate_certificate(s["lower_cert"], s["problem"], s["arrays"], s["lower"])
    validate_certificate(s["upper_cert"], s["problem"], s["arrays"], s["upper"])
    validate_bracket(s["bracket"], s["lower_cert"], s["upper_cert"])
    validate_feasibility(s["feasibility"], s["problem"])
    assert not gate_admissible(s["lower_cert"]) and not gate_admissible(s["upper_cert"])
    assert not gate_admissible(s["feasibility"]), "a schema fixture never passes the kinematic gate"
    for bound, formulation in (("lower", LOWER_ID), ("upper", CONTROL_ID)):
        method, _ = load_document(DATA / f"{record_stem(name, formulation)}.certificate.method.json")
        validate_method_report(method, s[f"{bound}_cert"], s["problem"], s[bound], s[f"{bound}_arrays"], s["arrays"])
        assert not gate_admissible(s[f"{bound}_cert"], method=method), "a schema fixture with a report is still no gate"


@pytest.mark.parametrize("name", FIXTURE_NAMES)
def test_fixture_fields_are_the_analytic_solution(name):
    s = load_set(name)
    case = ALL_CASES[name]
    k = case["regions"][0]["k"]
    nodes, tets = s["arrays"]["nodes"], s["arrays"]["tets"]
    sigma = s["lower_arrays"]["sigma"]
    velocity = s["upper_arrays"]["velocity"]
    assert s["lower"]["value"] == pytest.approx(case["target"])
    dissipation = 0.0
    for e, tet in enumerate(tets):
        _, volume = tet_gradients(nodes[tet[:4]])
        for vertex in range(4):
            assert yield_from_principal(sigma[e, vertex], k) == pytest.approx(0.0, abs=1e-12)
            d = strain_rate_at_vertex(nodes[tet], velocity[tet], vertex)
            dissipation += volume / 4.0 * support_function(d, k, 0.0)
    assert dissipation == pytest.approx(case["target"])
    assert s["upper"]["value"] == pytest.approx(dissipation)
    kinds = [item["kind"] for item in s["problem"]["boundary_sets"]]
    resultant = np.zeros(6)
    for row, ((e, f), index) in enumerate(zip(s["arrays"]["facets"], s["arrays"]["facet_set"])):
        normal, area = facet_geometry(nodes[tets[e][:4]], int(f))
        traction = voigt_n(normal) @ sigma[e, 0]
        kind = kinds[int(index)]
        if kind == "traction":
            np.testing.assert_allclose(traction, case["target"] * s["arrays"]["t1"][row, 0], atol=1e-12)
        elif kind == "free":
            np.testing.assert_allclose(traction, 0.0, atol=1e-12)
        elif kind == "symmetry":
            for tangent in tangent_basis(normal):
                assert tangent @ traction == pytest.approx(0.0, abs=1e-12)
        elif kind == "footing":
            resultant[:3] += area * traction
    if s["problem"]["footings"]:
        np.testing.assert_allclose(resultant[:3], s["lower_arrays"]["resultant"][0, :3], atol=1e-12)


# Malformed variants ---------------------------------------------------------------------------

def rehashed_problem(doc, arrays):
    doc["problem_sha256"] = document_hash(doc, arrays, "problem_sha256")


def problem_mutation(mutate):
    def run(code):
        s = load_set("vm_plate_shear")
        doc, arrays = copy.deepcopy(s["problem"]), {k: v.copy() for k, v in s["arrays"].items()}
        arrays = mutate(doc, arrays) or arrays
        rehashed_problem(doc, arrays)
        with pytest.raises(ContractError) as error:
            validate_problem(doc, arrays)
        assert error.value.code == code
    return run


def _set(path, value):
    def mutate(doc, arrays):
        target = doc
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
    return mutate


def _curve(doc, arrays):
    arrays["nodes"][arrays["tets"][0, 4]] += 0.01


def _invert(doc, arrays):
    arrays["tets"][0, [2, 3]] = arrays["tets"][0, [3, 2]]
    arrays["tets"][0, [5, 8]] = arrays["tets"][0, [8, 5]]
    arrays["tets"][0, [6, 7]] = arrays["tets"][0, [7, 6]]


def _drop_array(doc, arrays):
    del arrays["t1"]


def _integer_nodes(doc, arrays):
    arrays["nodes"] = arrays["nodes"].astype("<i8")


def _nan_node(doc, arrays):
    arrays["nodes"][0, 0] = np.nan


def _reorder_facets(doc, arrays):
    arrays["facets"] = arrays["facets"][::-1].copy()


def _missing_facet(doc, arrays):
    for name in ("facets", "facet_set", "t0", "t1"):
        arrays[name] = arrays[name][1:].copy()


def _traction_on_fixed(doc, arrays):
    index = [item["kind"] for item in doc["boundary_sets"]].index("fixed")
    arrays["t0"][np.nonzero(arrays["facet_set"] == index)[0][0], 0, 0] = 1.0


def _empty_set(doc, arrays):
    doc["boundary_sets"].append({"name": "unused", "kind": "free"})


def _bent_symmetry(doc, arrays):
    doc["boundary_sets"][2]["kind"] = "symmetry"
    doc["boundary_sets"][3]["kind"] = "symmetry"
    arrays["t1"][:] = 0.0
    arrays["facet_set"][arrays["facet_set"] == 3] = 2
    doc["boundary_sets"].pop(3)
    arrays["facet_set"][arrays["facet_set"] > 3] -= 1


def _restrained_with_load(doc, arrays):
    doc["queries"][0]["restrained"] = [[0, 0]]


def _no_reference_load(doc, arrays):
    doc["queries"][0]["F1"] = [[0.0] * 6]
    arrays["t1"][:] = 0.0


def _support_with_reference(doc, arrays):
    doc["queries"][0] = {"id": "q", "form": "support", "restrained": [], "direction": [[0, 1.0, 0, 0, 0, 0]]}


def _support_zero(doc, arrays):
    arrays["t1"][:] = 0.0
    doc["queries"][0] = {"id": "q", "form": "support", "restrained": [], "direction": [[0.0] * 6]}


def _fixed_on_direction(doc, arrays):
    doc["queries"][0] = {
        "id": "q", "form": "fixed_components", "restrained": [], "fixed": [[0, 0, 1.0]],
        "direction": [[1.0, 0, 0, 0, 0, 0]],
    }


PROBLEM_MUTATIONS = [
    ("wrong_format", _set(["format"], "other"), "format"),
    ("wrong_version", _set(["format_version"], 1), "format_version"),
    ("unknown_key", _set(["extra"], 1), "unknown_key"),
    ("voigt_order", _set(["conventions", "voigt"], ["xx", "yy", "zz", "yz", "zx", "xy"]), "conventions"),
    ("compression_positive", _set(["conventions", "sign"], "compression_positive"), "conventions"),
    ("reserved_tresca", _set(["regions", 0, "criterion"], "tresca"), "criterion_reserved"),
    ("unknown_criterion", _set(["regions", 0, "criterion"], "hoek_brown"), "criterion"),
    ("negative_strength", _set(["regions", 0, "k"], -1.0), "strength"),
    ("short_body_force", _set(["regions", 0, "b0"], [0.0, 0.0]), "body_force"),
    ("reserved_interface_law", _set(["footings", 0, "interface_law"], "coulomb"), "interface_law_reserved"),
    ("reserved_boundary", _set(["boundary_sets", 0, "kind"], "interface"), "boundary_kind_reserved"),
    ("unknown_boundary", _set(["boundary_sets", 0, "kind"], "roller"), "boundary_kind"),
    ("bad_domain", _set(["domain", "target"], "sphere"), "domain_target"),
    ("halfspace_without_walls", _set(["domain", "target"], "halfspace"), "domain_truncation"),
    ("zero_length_scale", _set(["length_scale"], 0.0), "length_scale"),
    ("curved_element", _curve, "element_curved"),
    ("inverted_element", _invert, "element_orientation"),
    ("missing_array", _drop_array, "array_names"),
    ("integer_nodes", _integer_nodes, "array_dtype"),
    ("nan_node", _nan_node, "array_non_finite"),
    ("facet_order", _reorder_facets, "facets_order"),
    ("missing_facet", _missing_facet, "facets_mismatch"),
    ("traction_on_fixed", _traction_on_fixed, "traction_on_other_kind"),
    ("empty_set", _empty_set, "boundary_set_empty"),
    ("bent_symmetry", _bent_symmetry, "symmetry_not_planar"),
    ("restrained_with_load", _restrained_with_load, "query_restrained_load"),
    ("no_reference_load", _no_reference_load, "query_no_reference_load"),
    ("support_with_reference_load", _support_with_reference, "query_support_reference_load"),
    ("support_zero_direction", _support_zero, "query_zero_direction"),
    ("fixed_component_in_direction", _fixed_on_direction, "query_fixed_direction"),
]


@pytest.mark.parametrize("label,mutate,code", PROBLEM_MUTATIONS, ids=[m[0] for m in PROBLEM_MUTATIONS])
def test_malformed_problem_is_refused(label, mutate, code):
    problem_mutation(mutate)(code)


def test_halfspace_needs_a_ground_plane():
    def mutate(doc, arrays):
        doc["domain"].update(target="halfspace", truncation_sets=["base"])
    problem_mutation(mutate)("domain_ground")


def test_declared_symmetry_set_must_be_a_symmetry_set():
    def mutate(doc, arrays):
        doc["domain"].update(target="halfspace", truncation_sets=["base"], symmetry_sets=["z0"],
                             ground={"axis": "y", "coordinate": 1.0})
    problem_mutation(mutate)("domain_symmetry")


def test_problem_hash_mismatch_is_refused():
    s = load_set("vm_plate_shear")
    s["arrays"]["t1"][0, 0, 1] += 1e-3
    with pytest.raises(ContractError) as error:
        validate_problem(s["problem"], s["arrays"])
    assert error.value.code == "hash_mismatch"


def result_case(bound="lower"):
    s = load_set("vm_plate_shear")
    doc = copy.deepcopy(s[bound])
    arrays = {k: v.copy() for k, v in s[f"{bound}_arrays"].items()}
    return s, doc, arrays


def check_result(s, doc, arrays, code):
    doc["result_sha256"] = document_hash(doc, arrays, "result_sha256")
    with pytest.raises(ContractError) as error:
        validate_result(doc, arrays, s["problem"], s["arrays"])
    assert error.value.code == code


def test_result_with_swapped_bound_label_is_refused():
    s, doc, arrays = result_case("lower")
    doc["bound"] = "upper"
    doc["formulation"] = "upper.p2c-vertex.v0"
    check_result(s, doc, arrays, "result_wrong_bound_arrays")


def test_result_for_another_problem_is_refused():
    s, doc, arrays = result_case("upper")
    doc["problem_sha256"] = "0" * 64
    check_result(s, doc, arrays, "hash_mismatch")


def test_result_for_unknown_query_is_refused():
    s, doc, arrays = result_case("lower")
    doc["query_id"] = "other"
    check_result(s, doc, arrays, "query_mismatch")


def test_time_limit_result_carries_no_number():
    s, doc, arrays = result_case("lower")
    doc["solver_status"] = "time_limit"
    check_result(s, doc, arrays, "result_value_without_solution")
    doc["value"] = doc["objective"] = None
    check_result(s, doc, arrays, "result_arrays_without_solution")
    doc["result_sha256"] = document_hash(doc, {}, "result_sha256")
    validate_result(doc, {}, s["problem"], s["arrays"])


def test_result_array_shape_and_interior_rules():
    s, doc, arrays = result_case("lower")
    arrays["sigma"] = arrays["sigma"][:-1]
    check_result(s, doc, arrays, "array_shape")
    s, doc, arrays = result_case("lower")
    arrays["interior_sigma"] = arrays["sigma"].copy()
    check_result(s, doc, arrays, "result_interior")
    doc["interior_value"] = 0.5
    doc["result_sha256"] = document_hash(doc, arrays, "result_sha256")
    validate_result(doc, arrays, s["problem"], s["arrays"])


def test_jump_result_needs_broken_velocity_and_an_applicable_criterion():
    s, doc, arrays = result_case("upper")
    doc["formulation"] = JUMP_ID
    check_result(s, doc, arrays, "result_wrong_bound_arrays")
    jump, jump_arrays = load_document(DATA / f"{record_stem('vm_plate_shear', JUMP_ID)}.result.json")
    problem = copy.deepcopy(s["problem"])
    problem["regions"][0] = {**problem["regions"][0], "criterion": "drucker_prager", "alpha": 0.3}
    jump["problem_sha256"] = problem["problem_sha256"]
    with pytest.raises(ContractError) as error:
        validate_result(jump, jump_arrays, problem, s["arrays"])
    assert error.value.code == "formulation_not_applicable"


@pytest.mark.parametrize("formulation,code", [
    ("upper.p2d-jump-bernstein.v1", "formulation_not_admitted"),
    ("upper.p2c-vertex.v9", "formulation_not_admitted"),
    ("lower.p1d-vertex.v0", "formulation_unknown"),
])
def test_upper_formulation_field_is_versioned(formulation, code):
    s, doc, arrays = result_case("upper")
    doc["formulation"] = formulation
    check_result(s, doc, arrays, code)


def test_formulation_is_not_part_of_the_problem():
    s = load_set("vm_plate_shear")
    assert "formulation" not in json.dumps(s["problem"])


def certificate_case(bound="lower"):
    s = load_set("vm_plate_shear")
    return s, copy.deepcopy(s[f"{bound}_cert"])


def check_certificate(s, doc, code, bound="lower", problem=None):
    with pytest.raises(ContractError) as error:
        validate_certificate(doc, problem or s["problem"], s["arrays"], s[bound])
    assert error.value.code == code


@pytest.mark.parametrize("key,value,code", [
    ("domain_claim", "box", "domain_claim"),
    ("domain_claim", "halfspace_by_zero_extension", "unsupported_domain_claim"),
    ("domain_claim", "halfspace_by_extension_elements", "unsupported_domain_claim"),
    ("domain_claim", "unverified", "unsupported_domain_claim"),
    ("flow_rule", "non_associated", "certificate_flow_rule"),
    ("classification", "maybe", "certificate_classification"),
    ("reasons", ["not_a_code"], "certificate_reason"),
    ("reasons", ["value_mismatch"], "certificate_reason"),
    ("result_sha256", "0" * 64, "hash_mismatch"),
    ("criteria", ["drucker_prager"], "certificate_criteria"),
    ("tolerances", {}, "tolerances"),
])
def test_malformed_certificate_is_refused(key, value, code):
    s, doc = certificate_case()
    doc[key] = value
    check_certificate(s, doc, code)


def test_estimate_certificate_needs_a_reason_and_no_value():
    s, doc = certificate_case()
    doc["classification"] = "estimate"
    check_certificate(s, doc, "certificate_reason")
    doc["reasons"] = ["no_admissible_correction"]
    check_certificate(s, doc, "certificate_value")
    doc["value_certified"] = doc["floor"] = None
    validate_certificate(doc, s["problem"], s["arrays"], s["lower"])
    assert not gate_admissible(doc)


def test_gate_refuses_reserved_checker_names_and_accepts_a_real_one():
    _, doc = certificate_case()
    assert not gate_admissible(doc, method=_method_report(doc))
    doc["checker"]["name"] = PROCEDURE_DEMO_CHECKER
    assert not gate_admissible(doc, method=_method_report(doc))
    doc["checker"]["name"] = "independent-checker"
    assert gate_admissible(doc, method=_method_report(doc))
    assert not gate_admissible(doc), "amendment 1: a gate certificate needs its method report"


@pytest.mark.parametrize("bound,key,value", [
    ("lower", "omega_eq", None), ("lower", "omega_eq", 1e-13), ("lower", "kappa_eq", None), ("lower", "kappa_eq", 2e8),
    ("lower", "tau_eq", 1e-9), ("lower", "projection_move", None), ("lower", "projection_move", 1e-5),
    ("lower", "rho_cone", None), ("lower", "rho_cone", 1e-9), ("lower", "repair_mu", None), ("lower", "repair_mu", 0.01),
    ("lower", "straightness_max", None), ("upper", "rho_flow", None), ("upper", "rho_flow", 1e-9),
    ("upper", "work_reference", None), ("upper", "work_reference", 0.0), ("upper", "work_reference", -1.0),
])
def test_certified_bound_must_show_its_numbers(bound, key, value):
    s, doc = certificate_case(bound)
    doc["checks"][key] = value
    doc["checker"]["name"] = "independent-checker"
    check_certificate(s, doc, "certificate_checks", bound)
    assert not gate_admissible(doc)


def test_certified_jump_bound_must_show_the_jump_cone():
    s, doc = certificate_case("upper")
    doc["formulation"] = JUMP_ID
    assert "rho_jump_cone" in certified_check_failures(doc, "ray")


def test_certified_floor_must_follow_its_rule():
    s, doc = certificate_case()
    doc["floor"] *= 2.0
    check_certificate(s, doc, "certificate_checks")
    s, doc = certificate_case()
    doc["checks"].update(omega_eq=1e-15, kappa_eq=1e4)
    check_certificate(s, doc, "certificate_checks")
    doc["floor"] = 1e-11 * abs(doc["value_certified"])
    validate_certificate(doc, s["problem"], s["arrays"], s["lower"])


def halfspace_problem():
    name = "vm_ground_patch_halfspace"
    return load_document(DATA / f"{case_stem(name)}.problem.json")


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
def test_domain_claim_is_per_bound(target, kind, weight, free, bound, expected):
    problem, arrays = halfspace_problem()
    walls = problem["domain"]["truncation_sets"]
    for item in problem["boundary_sets"]:
        if item["name"] in walls:
            item["kind"] = kind
    problem["regions"][0]["b0"] = [0.0, -1.0 if weight else 0.0, 0.0]
    problem["domain"]["exterior_load_free"] = free
    if target == "box":
        problem["domain"].update(target="box", truncation_sets=[], ground=None)
    assert supported_domain_claim(problem, arrays, bound) == expected


def test_zero_extension_needs_every_face_off_the_ground_accounted_for():
    problem, arrays = halfspace_problem()
    assert supported_domain_claim(problem, arrays, "upper") == "halfspace_by_zero_extension"
    problem["domain"]["truncation_sets"].remove("xmin")
    assert supported_domain_claim(problem, arrays, "upper") == "unverified", "a listed-walls rule alone is not enough"
    problem["domain"]["truncation_sets"].append("xmin")
    problem["domain"]["ground"]["coordinate"] = -1.0
    assert supported_domain_claim(problem, arrays, "upper") == "unverified", "a load off the ground plane"


def test_old_plate_field_cannot_claim_a_halfspace():
    s = load_set("vm_plate_shear")
    problem = copy.deepcopy(s["problem"])
    problem["domain"].update(target="halfspace", truncation_sets=["base"], ground={"axis": "y", "coordinate": 1.0})
    assert supported_domain_claim(problem, s["arrays"], "upper") == "unverified"


def test_bracket_claim_is_what_both_bounds_support():
    zero = "halfspace_by_zero_extension"
    assert bracket_claim(zero, zero) == zero
    assert bracket_claim("unverified", zero) == "box_stated"
    assert bracket_claim(zero, "unverified") == "box_stated"
    assert bracket_claim("box_stated", "box_stated") == "box_stated"


def bracket_set():
    s = load_set("vm_plate_shear")
    return copy.deepcopy(s["bracket"]), copy.deepcopy(s["lower_cert"]), copy.deepcopy(s["upper_cert"])


def expect_bracket(code, doc, low, high):
    with pytest.raises(ContractError) as error:
        validate_bracket(doc, low, high)
    assert error.value.code == code


def test_bracket_refuses_swapped_labels_and_estimates():
    doc, low, high = bracket_set()
    expect_bracket("bracket_labels", doc, high, low)
    doc, low, high = bracket_set()
    low["classification"] = "estimate"
    expect_bracket("bracket_uncertified", doc, low, high)


def test_bracket_names_its_formulations():
    doc, low, high = bracket_set()
    doc["upper_formulation"] = JUMP_ID
    expect_bracket("bracket_formulation", doc, low, high)


def test_bracket_reports_a_contradiction():
    doc, low, high = bracket_set()
    low["value_certified"] = doc["lower"] = 1.5
    expect_bracket("bracket_status", doc, low, high)
    doc["status"] = "contradiction"
    doc["gap"] = None
    validate_bracket(doc, low, high)


def test_bracket_order_within_floor_is_not_a_bracket():
    doc, low, high = bracket_set()
    low["value_certified"] = doc["lower"] = high["value_certified"] + 0.5 * (low["floor"] + high["floor"])
    doc["status"] = "order_within_floor"
    doc["gap"] = None
    validate_bracket(doc, low, high)


NO_WITNESS = {"witness_result_sha256": None, "witness_residual": None, "witness_projection_move": None,
              "witness_work_reference": None}


def feasibility_doc(**changes):
    s = load_set("vm_plate_shear")
    doc = copy.deepcopy(s["feasibility"])
    doc.update(changes)
    return s["problem"], doc


def test_feasibility_failure_is_a_recordable_result():
    problem, doc = feasibility_doc(status="formulation_infeasible", numerator=0.1, denominator=1.0,
                                   rho_feas=0.1, **NO_WITNESS)
    validate_feasibility(doc, problem)
    assert not gate_admissible(doc)


def test_feasibility_from_a_real_checker_is_gate_admissible():
    problem, doc = feasibility_doc(checker={"name": "independent-checker", "version": "1"})
    validate_feasibility(doc, problem)
    assert gate_admissible(doc)
    doc["status"] = "inconclusive"
    assert not gate_admissible(doc)


@pytest.mark.parametrize("changes,code", [
    ({"status": "formulation_feasible", "numerator": 0.1, "denominator": 1.0, "rho_feas": 0.1}, "feasibility_verdict"),
    (NO_WITNESS, "feasibility_verdict"),
    ({"witness_result_sha256": None, "witness_residual": None}, "feasibility_witness"),
    # Amendment 2, m2: the witness keeps its projection move within the limit and does positive reference work.
    ({"witness_projection_move": 2e-6}, "feasibility_verdict"),
    ({"witness_work_reference": 0.0}, "feasibility_verdict"),
    ({"witness_work_reference": 1e-13}, "feasibility_verdict"),  # below TOL_WORK_MIN, as 17.4 step 5
    ({"witness_work_reference": None}, "feasibility_witness"),
    ({"witness_projection_move": -1.0}, "feasibility_witness"),
    ({"witness_residual": 1e-3}, "feasibility_verdict"),
    ({"witness_residual": None}, "feasibility_witness"),
    ({"witness_result_sha256": "abc"}, "feasibility_witness"),
    ({"status": "formulation_infeasible"}, "feasibility_verdict"),
    ({"status": "inconclusive", "numerator": 1.0, "denominator": 1.0, "rho_feas": 1.0}, "feasibility_verdict"),
    ({"status": "solver_said_ok"}, "feasibility_status"),
    ({"rho_feas": 0.5}, "feasibility_ratio"),
    ({"formulation": "upper.p2d-jump-bernstein.v1"}, "formulation_not_admitted"),
    ({"denominator": 0.0}, "feasibility_number"),
    ({"nullity": 3}, "feasibility_rank"),
    ({"footing_modes": 0}, "feasibility_rank"),
])
def test_feasibility_record_cannot_hide_a_failure(changes, code):
    problem, doc = feasibility_doc(**changes)
    with pytest.raises(ContractError) as error:
        validate_feasibility(doc, problem)
    assert error.value.code == code


def test_nullity_does_not_decide_feasibility():
    problem, doc = feasibility_doc(rank=10, nullity=0)
    validate_feasibility(doc, problem)


def test_frictional_feasibility_rests_on_the_witness_alone():
    problem, doc = feasibility_doc()
    problem = copy.deepcopy(problem)
    problem["regions"][0] = {**problem["regions"][0], "criterion": "drucker_prager", "alpha": 0.1}
    doc.update(rho_feas=2e-11, numerator=2e-11, witness_residual=2e-11)
    validate_feasibility(doc, problem)
    doc.update(status="formulation_infeasible", rho_feas=1e-3, numerator=1e-3, witness_residual=1e-3)
    with pytest.raises(ContractError) as error:
        validate_feasibility(doc, problem)
    assert error.value.code == "feasibility_verdict"
    doc["status"] = "inconclusive"
    validate_feasibility(doc, problem)


def schema_program(s, bound, formulation=None):
    query = s["problem"]["queries"][0]
    formulation = formulation or (LOWER_ID if bound == "lower" else CONTROL_ID)
    n_var, n_row, n_zero, cones = expected_program_shape(s["problem"], s["arrays"], query, bound, formulation)
    assert n_row == n_zero + sum(c["dim"] * c["count"] for c in cones)
    doc = {
        "format": FORMAT_PROGRAM, "format_version": FORMAT_VERSION, "arrays": "unused.npz",
        "problem_sha256": s["problem"]["problem_sha256"], "query_id": "q", "bound": bound,
        "formulation": formulation, "objective_sign": -1 if bound == "lower" else 1,
        "n_variables": n_var, "n_rows": n_row, "n_zero": n_zero, "cones": cones,
    }
    arrays = {
        "A_indptr": np.zeros(n_var + 1, dtype="<i8"), "A_indices": np.zeros(0, dtype="<i8"),
        "A_data": np.zeros(0, dtype="<f8"), "b": np.zeros(n_row, dtype="<f8"), "c": np.zeros(n_var, dtype="<f8"),
    }
    return doc, arrays


@pytest.mark.parametrize("bound,formulation", [("lower", LOWER_ID), ("upper", CONTROL_ID), ("upper", JUMP_ID)])
def test_program_schema_and_canonical_shape(bound, formulation):
    s = load_set("vm_plate_shear")
    doc, arrays = schema_program(s, bound, formulation)
    validate_program(doc, arrays, s["problem"], s["arrays"])
    n_e = len(s["arrays"]["tets"])
    interior = len(mesh_facets(s["arrays"]["tets"])[0])
    if formulation == LOWER_ID:
        assert doc["n_variables"] == 1 + 24 * n_e
    elif formulation == CONTROL_ID:
        assert doc["n_variables"] == 3 * len(s["arrays"]["nodes"]) + 6 + 24 * n_e
    else:
        assert doc["n_variables"] == 30 * n_e + 6 + 24 * n_e + 18 * interior
        assert doc["cones"] == [{"dim": 6, "count": 4 * n_e}, {"dim": 3, "count": 6 * interior}]
    doc["n_rows"] += 1
    with pytest.raises(ContractError) as error:
        validate_program(doc, arrays, s["problem"], s["arrays"])
    assert error.value.code == "program_shape"


def test_program_refuses_a_bad_sparse_matrix_and_sign():
    s = load_set("vm_plate_shear")
    doc, arrays = schema_program(s, "upper")
    arrays["A_indptr"][-1] = 3
    with pytest.raises(ContractError) as error:
        validate_program(doc, arrays, s["problem"], s["arrays"])
    assert error.value.code == "program_csc"
    doc, arrays = schema_program(s, "upper")
    doc["objective_sign"] = -1
    with pytest.raises(ContractError) as error:
        validate_program(doc, arrays, s["problem"], s["arrays"])
    assert error.value.code == "program_objective_sign"


# ---------------------------------------------------------------------------------------------
# Tests: the documents agree with this module
# ---------------------------------------------------------------------------------------------

def interface_text():
    return INTERFACE.read_text(encoding="utf-8")


def test_tolerance_table_matches_the_interface():
    found = dict(re.findall(r"^\| `((?:TOL|FLOOR|MAX_GAP)_[A-Z0-9_]+)` \| ([0-9.e+-]+) \|", interface_text(), re.M))
    assert {name: float(value) for name, value in found.items()} == TOLERANCES


@pytest.mark.parametrize("terms", [
    sorted(REASON_CODES), sorted(DOMAIN_CLAIMS), sorted(CRITERIA | RESERVED_CRITERIA),
    sorted(BOUNDARY_KINDS | RESERVED_BOUNDARY_KINDS), sorted(INTERFACE_LAWS | RESERVED_INTERFACE_LAWS),
    sorted(QUERY_FORMS), sorted(SOLVER_STATUS_WITH_VALUE | SOLVER_STATUS_NO_VALUE),
    sorted(FEASIBILITY_STATUS), sorted(CHECK_KEYS),
    sorted(FORMULATIONS["lower"] | FORMULATIONS["upper"]),
    [FORMAT_PROBLEM, FORMAT_RESULT, FORMAT_PROGRAM, FORMAT_FEASIBILITY, FORMAT_CERTIFICATE, FORMAT_BRACKET,
     FORMAT_METHOD, FORMAT_REFUSAL],
    sorted(KAPPA_SVD_METHODS | KAPPA_EIGEN_METHODS), sorted(DEPENDENT_ROUTES), sorted(INTERIOR_SOURCES),
    sorted(REFUSAL_KINDS), sorted(REFUSAL_STAGES),
])
def test_interface_defines_every_term(terms):
    text = interface_text()
    missing = [term for term in terms if f"`{term}`" not in text]
    assert not missing


def test_interface_conventions_object_is_literal():
    block = re.search(r"### 15\.7.*?```json\n(.*?)```", interface_text(), re.S).group(1)
    assert json.loads(block) == CONVENTIONS


def test_interface_separates_owner_amendments_from_orchestrator_rulings():
    text = interface_text()
    owner = text.index("### 18.1 Accepted owner amendments")
    rulings = text.index("### 18.2 Orchestrator P0 rulings")
    assert owner < rulings
    for item in ("A1", "A2", "A3", "A4"):
        assert f"| {item} |" in text[owner:rulings]
    for item in ("OD1", "OD2", "OD3", "OD4", "OD5"):
        assert f"| {item} |" in text[owner:rulings]
    for item in ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9", "R10", "R11"):
        assert f"| {item} |" in text[rulings:]
    assert "frozen at `03f59c1a4`" in text[:600] and "amendment 1" in text[:600]
    assert "pending" not in text[text.index("### 18.3"):text.index("## 19.")].lower(), "the owner has ruled"


# Negative fixtures -----------------------------------------------------------------------------

NEGATIVE = load_document(DATA / "negative_fixtures.json")[0]["fixtures"] if (DATA / "negative_fixtures.json").exists() else []
HARNESS_ITEMS = {
    "one perturbed stress component", "stress scaled by 1.01", "a facet traction jump",
    "a wrong load normalisation", "swapped LB and UB labels", "a curved tetrahedron",
    "a mismatched problem hash", "an unsupported domain_claim",
    "tr(D) != 0 at one vertex of an isochoric Upper field", "a footprint node off the rigid motion",
    "nonzero wall velocity under a zero-extension claim",
    "a jump admissible at facet nodes but not across the facet",
    "wrong gravity work or wrong dissipation quadrature", "an equality residual just above the floor",
}
# ADR 237 section 11.5, row by row, and the fixture that discharges each row.
ADR_11_5 = {
    "wrong load normalisation": ["wrong_load_normalisation", "wrong_lower_value"],
    "swapped Lower and Upper labels": ["swapped_bound_labels"],
    "curved tetrahedron": ["curved_tetrahedron"],
    "mismatched problem hash": ["mismatched_problem_hash"],
    "unsupported domain_claim": ["unsupported_domain_claim", "unlisted_side_face"],
    "equality residual just above the threshold": ["equality_just_above_floor"],
    "one perturbed stress component": ["perturbed_stress_component"],
    "stress scaled by 1.01": ["stress_scaled_1p01"],
    "facet traction jump": ["facet_traction_jump"],
    "wrong gravity in equilibrium": ["wrong_gravity_in_equilibrium"],
    "tr(D) != 0 at one vertex of an isochoric field": ["isochory_one_vertex"],
    "Drucker-Prager flow-cone violation at one vertex": ["dp_flow_cone_one_vertex"],
    "footprint node off the rigid motion": ["footprint_off_rigid"],
    "non-zero support velocity": ["support_velocity_nonzero"],
    "non-zero wall velocity under a zero-extension claim": ["wall_velocity_zero_extension"],
    "non-zero normal jump at one facet node under isochoric flow": ["normal_jump_isochoric"],
    "jump admissible at all six facet nodes but not across the facet": ["jump_between_nodes"],
    "wrong gravity work": ["wrong_gravity_work"],
    "wrong dissipation quadrature": ["wrong_dissipation_quadrature"],
    "non-positive reference work rate": ["nonpositive_reference_work"],
    "a finite dissipation charge offered in place of restoring a violated cone": ["dissipation_charge_as_repair"],
}
RECOMPUTED_VALUE = {"perturb_velocity", "perturb_element_velocity", "edge_normal_jump", "side_opening", "negate_velocity",
                    "dissipation_charge", "side_translation", "restrained_shear", "apex_zero_omits", "apex_zero_hides",
                    "apex_zero_full"}


def apply_negative(entry):
    problem, problem_arrays = load_document(DATA / f"{entry['problem']}.problem.json")
    result, arrays = load_document(DATA / f"{entry['result']}.json") if entry["result"] else (None, {})
    m = entry["mutation"]
    kind = m["kind"]
    certificate = None
    if kind == "perturb_sigma":
        arrays["sigma"][m["element"], m["vertex"], m["component"]] += m["amount"]
    elif kind == "scale_lower":
        arrays["sigma"] *= m["factor"]
        arrays["resultant"] *= m["factor"]
        result["value"] *= m["factor"]
        result["objective"] *= m["factor"]
    elif kind == "element_offset":
        arrays["sigma"][m["element"], :, m["component"]] += m["amount"]
    elif kind == "zero_gravity_lower":
        arrays["sigma"][:] = 0.0
        arrays["sigma"][..., 5] = result["value"]
    elif kind == "scale_value":
        result["value"] *= m["factor"]
        result["objective"] *= m["factor"]
    elif kind == "relabel_upper":
        result["bound"], result["formulation"] = "upper", CONTROL_ID
    elif kind == "move_mid_node":
        problem_arrays["nodes"][problem_arrays["tets"][m["element"], m["local_node"]]] += m["amount"]
        problem["problem_sha256"] = document_hash(problem, problem_arrays, "problem_sha256")
    elif kind == "problem_hash":
        result["problem_sha256"] = "0" * 64
    elif kind == "certificate_claim":
        certificate = load_document(DATA / f"{entry['result'].replace('.result', '.certificate')}.json")[0]
        certificate["domain_claim"] = m["claim"]
    elif kind == "drop_truncation_set":
        certificate = load_document(DATA / f"{entry['result'].replace('.result', '.certificate')}.json")[0]
        problem["domain"]["truncation_sets"].remove(m["name"])
        problem["problem_sha256"] = document_hash(problem, problem_arrays, "problem_sha256")
        result["problem_sha256"] = certificate["problem_sha256"] = problem["problem_sha256"]
    elif kind == "perturb_velocity":
        arrays["velocity"][m["node"], m["component"]] += m["amount"]
    elif kind == "perturb_element_velocity":
        arrays["element_velocity"][m["element"], m["local_node"], m["component"]] += m["amount"]
    elif kind == "edge_normal_jump":
        grads, _ = tet_gradients(problem_arrays["nodes"][problem_arrays["tets"][m["element"]][:4]])
        a, b = EDGES[m["edge"]]
        direction = np.cross(grads[a], grads[b])
        arrays["element_velocity"][m["element"], 4 + m["edge"]] += m["amount"] * direction / np.linalg.norm(direction)
    elif kind == "side_opening":
        arrays["element_velocity"][np.array(m["elements"]), :4] += m["amount"] * np.array(m["normal"])
    elif kind == "side_translation":
        arrays["element_velocity"][np.array(m["elements"])] += m["amount"] * np.array(m["direction"])
    elif kind == "restrained_shear":
        arrays["velocity"][:, m["component"]] += m["amount"] * problem_arrays["nodes"][:, m["axis"]]
        arrays["footing_velocity"][m["footing"], m["component"]] += m["amount"]
    elif kind == "negate_velocity":
        arrays["velocity"] *= -1.0
        arrays["footing_velocity"] *= -1.0
    elif kind == "dissipation_charge":
        arrays["velocity"][m["node"], m["component"]] += m["amount"]
        conditions, _ = signed_flow_conditions(problem, problem_arrays, arrays["velocity"])
        arrays["dissipation"] = np.where(conditions > 0.0, m["charge"], 0.0).astype("<f8")
    elif kind == "value_rule":
        value = fixture_upper_value(problem, problem_arrays, arrays["velocity"], arrays["footing_velocity"],
                                    quadrature=m["quadrature"], body=m["body"])
        result["value"] = result["objective"] = value
    elif kind in ("apex_zero_omits", "apex_zero_hides", "apex_zero_full"):
        candidate, interior = apex_fixture_fields(problem_arrays, arrays["element_velocity"], m["elements"],
                                                  m["direction"], m["element"], m["local_node"], m["perturbation"])
        arrays["element_velocity"] = candidate
        arrays["interior_element_velocity"] = interior
        arrays["interior_footing_velocity"] = arrays["footing_velocity"].copy()
        arrays["apex_zero"] = np.array(m["apex_zero"], dtype="<i8").reshape(-1, 2)
    elif kind in ("add_apex_zero", "add_apex_zero_unhashed"):
        arrays["apex_zero"] = np.array(m["apex_zero"], dtype="<i8").reshape(-1, 2)
    else:
        fail("negative_fixture", kind)
    if kind in RECOMPUTED_VALUE:
        name = "element_velocity" if "element_velocity" in arrays else "velocity"
        result["value"] = fixture_upper_value(problem, problem_arrays, arrays[name], arrays["footing_velocity"])
        result["objective"] = fixture_upper_objective(problem, problem_arrays, arrays[name], arrays["footing_velocity"])
    if "sigma" in arrays and kind in ("perturb_sigma", "element_offset", "zero_gravity_lower"):
        arrays["resultant"] = fixture_resultant(problem, problem_arrays, arrays["sigma"]).astype("<f8")
    if result is not None and kind not in ("problem_hash", "add_apex_zero_unhashed"):
        result["result_sha256"] = document_hash(result, arrays, "result_sha256")
        if certificate is not None:
            certificate["result_sha256"] = result["result_sha256"]
    return problem, problem_arrays, result, arrays, certificate


def test_negative_fixtures_cover_the_harness_and_every_adr_row():
    items = {entry["harness_item"] for entry in NEGATIVE}
    assert HARNESS_ITEMS <= items
    ids = [entry["id"] for entry in NEGATIVE]
    assert len(set(ids)) == len(ids)
    for row, names in ADR_11_5.items():
        assert names and set(names) <= set(ids), row
    assert len(ADR_11_5) == 21
    for entry in NEGATIVE:
        assert isinstance(entry["expect"], str), "one expected code per fixture"
        if entry["stage"] == "schema":
            assert entry["expect"] in ("result_wrong_bound_arrays", "element_curved", "hash_mismatch",
                                       "unsupported_domain_claim", "apex_zero_lower")
        elif entry["stage"] == "certify":  # amendment 2, M1: a positive companion
            assert entry["expect"] == "certified_bound"
        else:
            assert entry["expect"] in REASON_CODES


@pytest.mark.parametrize("entry", [e for e in NEGATIVE if e["stage"] == "schema"], ids=lambda e: e["id"])
def test_schema_negative_fixture_is_refused(entry):
    problem, problem_arrays, result, arrays, certificate = apply_negative(entry)
    with pytest.raises(ContractError) as error:
        if entry["mutation"]["kind"] == "move_mid_node":
            validate_problem(problem, problem_arrays)
        elif certificate is not None:
            validate_certificate(certificate, problem, problem_arrays, result)
        else:
            validate_result(result, arrays, problem, problem_arrays)
    assert error.value.code == entry["expect"]


def _plate_vertex_row_error(entry, problem_arrays, arrays):
    m = entry["mutation"]
    grads, _ = tet_gradients(problem_arrays["nodes"][problem_arrays["tets"][m["element"]][:4]])
    return abs(grads[m["vertex"], 0] * m["amount"]) / float(np.sum(np.abs(grads)))


def assert_rows_untouched(problem, problem_arrays, field, footing_velocity):
    """Isochory at every vertex, the fixed and symmetry rows (set normal) and the rigid footprint hold for a fixture
    field: no projection row of interface 17.4 step 2 other than the normal jumps and restrained components."""
    nodes, tets = problem_arrays["nodes"], problem_arrays["tets"]
    per_element = element_velocities(problem_arrays, field)
    for e, tet in enumerate(tets):
        for vertex in range(4):
            theta = strain_split(strain_rate_at_vertex(nodes[tet], per_element[e], vertex))[0]
            assert abs(theta) <= 1e-12, "isochory holds at every vertex"
    normals = {}
    for (e, f), index in zip(problem_arrays["facets"], problem_arrays["facet_set"]):
        e, f, item = int(e), int(f), problem["boundary_sets"][int(index)]
        normal, _ = facet_geometry(nodes[tets[e][:4]], f)
        n_set = normals.setdefault(int(index), symmetry_set_normal(problem, problem_arrays, int(index)))
        local = per_element[e][list(FACET_NODES[f])]
        if item["kind"] == "fixed":
            np.testing.assert_array_equal(local, 0.0)
        elif item["kind"] == "symmetry":
            np.testing.assert_allclose(local @ n_set, 0.0, atol=1e-14)
        elif item["kind"] == "footing":
            x_ref = np.array(problem["footings"][item["footing"]]["x_ref"])
            u, w = footing_velocity[item["footing"], :3], footing_velocity[item["footing"], 3:]
            rigid = u + np.cross(w, nodes[tets[e][list(FACET_NODES[f])]] - x_ref)
            np.testing.assert_allclose(local, rigid, atol=1e-14)


@pytest.mark.parametrize("entry", [e for e in NEGATIVE if e["stage"] != "schema"], ids=lambda e: e["id"])
def test_checker_negative_fixture_is_well_formed_and_live(entry):
    """The named defect is present and every earlier stage of the procedure passes, so one code is expected."""
    problem, problem_arrays, result, arrays, _ = apply_negative(entry)
    validate_result(result, arrays, problem, problem_arrays)
    base, base_arrays = load_document(DATA / f"{entry['result']}.json")
    m = entry["mutation"]
    nodes, tets = problem_arrays["nodes"], problem_arrays["tets"]
    fixture = entry["id"]
    move = TOLERANCES["TOL_PROJECTION_MOVE"]
    velocity_name = "element_velocity" if "element_velocity" in arrays else "velocity"
    if result["bound"] == "upper":
        true_value = fixture_upper_value(problem, problem_arrays, arrays[velocity_name], arrays["footing_velocity"])
        if entry["expect"] == "value_mismatch":
            assert abs(result["value"] - true_value) > TOLERANCES["TOL_VALUE_MATCH"] * abs(true_value)
            assert true_value == pytest.approx(base["value"], rel=1e-12)
        else:
            assert result["value"] == pytest.approx(true_value, rel=1e-14), "the value is that of the mutated field"
    else:
        np.testing.assert_allclose(arrays["resultant"], fixture_resultant(problem, problem_arrays, arrays["sigma"]),
                                   atol=1e-12, err_msg="the stored resultant is that of the mutated field")
    if fixture == "perturbed_stress_component":
        assert m["amount"] > move * np.max(np.abs(base_arrays["sigma"]))
    elif fixture == "equality_just_above_floor":
        assert 2.9 * tau_eq() <= _plate_vertex_row_error(entry, problem_arrays, arrays) <= 3.1 * tau_eq()
        assert m["tau_multiple"] == 3.0
    elif fixture == "stress_scaled_1p01":
        k = problem["regions"][0]["k"]
        s_ref = float(np.max(np.abs(arrays["sigma"])))
        assert max(yield_from_principal(s, k) for s in arrays["sigma"].reshape(-1, 6)) / s_ref > TOLERANCES["TOL_CONE_ACCEPT"]
    elif fixture == "facet_traction_jump":
        interior, _ = mesh_facets(tets)
        worst = 0.0
        for (e1, f1), (e2, _) in interior:
            if m["element"] in (e1, e2):
                normal, _ = facet_geometry(nodes[tets[e1][:4]], f1)
                vertex1 = FACET_VERTICES[f1][0]
                vertex2 = [int(v) for v in tets[e2]].index(int(tets[e1][vertex1]))
                jump = voigt_n(normal) @ (arrays["sigma"][e1, vertex1] - arrays["sigma"][e2, vertex2])
                worst = max(worst, float(np.max(np.abs(jump))))
        assert worst > move
    elif fixture == "wrong_lower_value":
        np.testing.assert_array_equal(arrays["sigma"], base_arrays["sigma"])
        assert abs(result["value"] - base["value"]) > move * abs(base["value"])
    elif fixture == "wrong_gravity_in_equilibrium":
        body = np.array(problem["regions"][0]["b0"])
        for e, tet in enumerate(tets):
            grads, _ = tet_gradients(nodes[tet[:4]])
            divergence = sum(voigt_n(grads[i]) @ arrays["sigma"][e, i] for i in range(4))
            np.testing.assert_allclose(divergence, 0.0, atol=1e-12)
        assert np.max(np.abs(body)) > move
    elif fixture in ("wrong_load_normalisation", "wrong_gravity_work", "wrong_dissipation_quadrature"):
        pass  # covered by the value comparison above
    elif fixture == "isochory_one_vertex":
        theta = [strain_split(strain_rate_at_vertex(nodes[tet], arrays["velocity"][tet], v))[0]
                 for tet in tets if m["node"] in tet for v in range(4)]
        assert max(abs(t) for t in theta) > TOLERANCES["TOL_FLOW_ACCEPT"]
    elif fixture == "footprint_off_rigid":
        assert nodes[m["node"], 1] == 1.0
        rigid = base_arrays["footing_velocity"][0, :3]
        assert np.linalg.norm(arrays["velocity"][m["node"]] - rigid) == pytest.approx(m["amount"])
    elif fixture == "support_velocity_nonzero":
        assert nodes[m["node"], 1] == 0.0 and problem["domain"]["target"] == "box"
        assert np.linalg.norm(arrays["velocity"][m["node"]]) == pytest.approx(m["amount"])
    elif fixture == "wall_velocity_zero_extension":
        assert supported_domain_claim(problem, problem_arrays, "upper") == "halfspace_by_zero_extension"
        assert np.max(np.abs(arrays["element_velocity"][m["element"], m["local_node"]])) == pytest.approx(m["amount"])
    elif fixture == "normal_jump_isochoric":
        field = arrays["element_velocity"]
        assert_rows_untouched(problem, problem_arrays, field, arrays["footing_velocity"])
        worst = 0.0
        for (e1, f1), (e2, _) in mesh_facets(tets)[0]:
            normal, _ = facet_geometry(nodes[tets[e1][:4]], f1)
            worst = max(worst, float(np.max(np.abs(bernstein_from_nodal(facet_nodal_jump(tets, field, e1, f1, e2))
                                                   @ normal))))
        assert worst > move, "the normal jump is violated"
    elif fixture == "vm_jump_opening":
        field = arrays["element_velocity"]
        assert problem["regions"][0]["criterion"] == "von_mises"
        assert_rows_untouched(problem, problem_arrays, field, arrays["footing_velocity"])
        openings = []
        for (e1, f1), (e2, _) in mesh_facets(tets)[0]:
            normal, _ = facet_geometry(nodes[tets[e1][:4]], f1)
            opening = bernstein_from_nodal(facet_nodal_jump(tets, field, e1, f1, e2)) @ normal
            if np.max(np.abs(opening)) > 1e-12:
                openings.append(opening)
        assert openings and all(np.all(o >= 0.5 * m["amount"]) for o in openings), "every moved coefficient opens"
    elif fixture == "restrained_footing_component":
        f, c = m["footing"], m["component"]
        validate_problem(problem, problem_arrays)
        assert [f, c] in problem["queries"][0]["restrained"]
        assert arrays["footing_velocity"][f, c] == m["amount"], "the restrained component moves"
        assert_rows_untouched(problem, problem_arrays, arrays["velocity"], arrays["footing_velocity"])
    elif fixture == "nonpositive_reference_work":
        assert result["value"] == pytest.approx(-base["value"]), "W1 = -1 with unchanged dissipation"
    elif fixture == "dp_flow_cone_one_vertex":
        conditions, size = signed_flow_conditions(problem, problem_arrays, arrays["velocity"])
        assert conditions.max() / size > TOLERANCES["TOL_FLOW_ACCEPT"]
    elif fixture == "dissipation_charge_as_repair":
        conditions, size = signed_flow_conditions(problem, problem_arrays, arrays["velocity"])
        assert TOLERANCES["TOL_CONE_MARGIN"] < conditions.max() / size < TOLERANCES["TOL_FLOW_ACCEPT"]
        assert np.any(arrays["dissipation"] > 0.0) and "interior_velocity" not in arrays
    elif fixture == "jump_between_nodes":
        field = arrays["element_velocity"]
        conditions, size = signed_flow_conditions(problem, problem_arrays, field)
        assert conditions.max() < -0.1 * size, "every volume condition stays strictly inside"
        assert not any(item["kind"] in ("fixed", "symmetry", "footing") for item in problem["boundary_sets"]), "no equality row"
        rows = signed_jump_conditions(problem, problem_arrays, field)
        assert max(float(nodal.max()) for nodal, _ in rows) <= 1e-12 * m["amount"], "passes at all six nodes of every facet"
        assert max(float(coefficient.max()) for _, coefficient in rows) >= 0.9 * m["amount"], "fails across a facet"
    elif fixture == "apex_zero_omits_a_condition":
        assert_apex_omission_is_live(problem, problem_arrays, arrays, m)
    elif fixture == "apex_zero_hides_violation":
        # Amendment 2, M1: a checker that drops Z from the tests but forgets the Z rows would certify this field.
        rows = apex_zero_rows(arrays)
        field = arrays["element_velocity"]
        signed, scale = signed_jump_coefficients(problem, problem_arrays, field)
        worst = max(signed[index] for _, index in rows) / scale
        assert 1e-4 < worst < 1e-2, "Z hides a condition violated at about 1e-3 D_scale"
        others = [i for i in range(len(signed)) if (1, i) not in rows]
        assert max(signed[i] for i in others) <= 0.0, "outside Z every condition holds"
        flow, size = signed_flow_conditions(problem, problem_arrays, field)
        assert flow.max() < -0.1 * size
        projected = project_onto_apex_rows(problem_arrays, field, rows)
        move = float(np.max(np.abs(projected - field))) / float(np.max(np.abs(field)))
        assert move > TOLERANCES["TOL_PROJECTION_MOVE"], "the Z rows move the field past the limit"
    elif fixture == "apex_zero_full_face":
        assert_full_face_certifies(problem, problem_arrays, arrays, result)
    else:
        fail("negative_fixture", fixture)


def assert_full_face_certifies(problem, problem_arrays, arrays, result):
    """Amendment 2, M1 positive companion: with every apex condition in Z, the step-2 projection is small, every
    condition outside Z holds within the margin (P empty, no repair), the Z arguments are at rounding level,
    W1 is positive and the value of the projected field agrees with the candidate's within TOL_VALUE_MATCH."""
    margin = TOLERANCES["TOL_CONE_MARGIN"]
    field = arrays["element_velocity"]
    rows = apex_zero_rows(arrays)
    projected = project_onto_apex_rows(problem_arrays, field, rows)
    move = float(np.max(np.abs(projected - field))) / float(np.max(np.abs(field)))
    assert 0.0 < move <= TOLERANCES["TOL_PROJECTION_MOVE"]
    signed, scale = signed_jump_coefficients(problem, problem_arrays, projected)
    outside = [i for i in range(len(signed)) if (1, i) not in rows]
    assert max(signed[i] for i in outside) <= margin * scale, "P is empty"
    _, coefficients = facet_coefficients(problem_arrays, projected)
    rho_apex = max(float(np.linalg.norm(coefficients[i // 6, i % 6])) for _, i in rows) / scale
    assert rho_apex <= margin, "rho_apex within the cone rounding margin"
    flow, size = signed_flow_conditions(problem, problem_arrays, projected)
    assert flow.max() < -0.1 * size
    _, _, w1 = fixture_upper_terms(problem, problem_arrays, projected, arrays["footing_velocity"])
    assert w1 > 0.0
    value = fixture_upper_value(problem, problem_arrays, projected, arrays["footing_velocity"])
    assert abs(value - result["value"]) <= TOLERANCES["TOL_VALUE_MATCH"] * abs(value)


def assert_apex_omission_is_live(problem, problem_arrays, arrays, m):
    """Amendment 2: with Z the projection is small, the volume conditions stay strictly inside, and exactly the
    omitted condition is in P where the interior field (on the face) is at the apex; with Z plus it, P is empty."""
    margin, move_limit = TOLERANCES["TOL_CONE_MARGIN"], TOLERANCES["TOL_PROJECTION_MOVE"]
    field, interior = arrays["element_velocity"], arrays["interior_element_velocity"]
    rows = apex_zero_rows(arrays)
    omitted = tuple(m["omitted"])
    assert omitted not in rows and rows == sorted(rows)
    assert not any(item["kind"] in ("fixed", "symmetry", "footing") for item in problem["boundary_sets"]), \
        "the Z rows are the only projection rows"
    projected = project_onto_apex_rows(problem_arrays, field, rows)
    v_ref = float(np.max(np.abs(field)))
    assert 0.0 < float(np.max(np.abs(projected - field))) / v_ref <= move_limit
    flow, size = signed_flow_conditions(problem, problem_arrays, projected)
    assert flow.max() < -0.1 * size, "every volume condition stays strictly inside"
    signed, scale = signed_jump_coefficients(problem, problem_arrays, projected)
    assert margin < signed[omitted[1]] / scale < TOLERANCES["TOL_FLOW_ACCEPT"], "the omitted condition is in P"
    assert np.max(np.abs(apex_zero_row_values(problem_arrays, projected, rows))) <= 1e-12 * scale, "Z rows hold"
    np.testing.assert_array_equal(apex_zero_row_values(problem_arrays, interior, rows + [omitted]), 0.0)
    signed_int, scale_int = signed_jump_coefficients(problem, problem_arrays, interior)
    assert scale_int > 0.0 and signed_int.max() <= 0.0
    assert signed_int[omitted[1]] > -TOLERANCES["TOL_INTERIOR_MIN"] * scale_int, "the interior cannot clear it"
    full = project_onto_apex_rows(problem_arrays, field, sorted(rows + [omitted]))
    signed_full, scale_full = signed_jump_coefficients(problem, problem_arrays, full)
    assert signed_full.max() <= margin * scale_full, "with the omitted condition declared, P is empty"


def test_good_halfspace_and_apex_fixtures_are_admissible():
    problem, arrays = halfspace_problem()
    result, fields = load_document(DATA / f"{record_stem('vm_ground_patch_halfspace', JUMP_ID)}.result.json")
    validate_problem(problem, arrays)
    validate_result(result, fields, problem, arrays)
    field, nodes, tets = fields["element_velocity"], arrays["nodes"], arrays["tets"]
    assert result["value"] > 0.0
    for e, tet in enumerate(tets):
        for m in range(4):
            assert abs(strain_split(strain_rate_at_vertex(nodes[tet], field[e], m))[0]) < 1e-13
    for (e1, f1), (e2, _) in mesh_facets(tets)[0]:
        normal, _ = facet_geometry(nodes[tets[e1][:4]], f1)
        np.testing.assert_allclose(bernstein_from_nodal(facet_nodal_jump(tets, field, e1, f1, e2)) @ normal, 0.0, atol=1e-13)
    names = [item["name"] for item in problem["boundary_sets"]]
    for (e, f), index in zip(arrays["facets"], arrays["facet_set"]):
        if names[int(index)] != "ground":
            np.testing.assert_array_equal(field[e][list(FACET_NODES[int(f)])], 0.0)
    certificate = load_document(DATA / f"{record_stem('vm_ground_patch_halfspace', JUMP_ID)}.certificate.json")[0]
    validate_certificate(certificate, problem, arrays, result)
    assert certificate["domain_claim"] == "halfspace_by_zero_extension"
    apex, apex_arrays = load_document(DATA / f"{case_stem('dp_free_hydrostatic')}.problem.json")
    upper = load_document(DATA / f"{record_stem('dp_free_hydrostatic', JUMP_ID)}.result.json")[0]
    assert upper["value"] == pytest.approx(TARGET_FORMULAS["dp_hydrostatic"]({"k": 1.0, "alpha": 0.2}))
    validate_problem(apex, apex_arrays)


def test_correct_gravity_and_dissipation_rules_differ_from_the_wrong_ones():
    problem, arrays = load_document(DATA / f"{case_stem('vm_dead_gravity_shear_layer')}.problem.json")
    upper, fields = load_document(DATA / f"{record_stem('vm_dead_gravity_shear_layer', CONTROL_ID)}.result.json")
    velocity, footing = fields["velocity"], fields["footing_velocity"]
    good = fixture_upper_value(problem, arrays, velocity, footing)
    assert good == pytest.approx(upper["value"])
    assert good >= ALL_CASES["vm_dead_gravity_shear_layer"]["target"], "an admissible Upper field bounds from above"
    assert fixture_upper_value(problem, arrays, velocity, footing, body="lumped") != pytest.approx(good, rel=1e-6)
    assert fixture_upper_value(problem, arrays, velocity, footing, quadrature="centroid") < good


def probe_rows():
    rows = {}
    for line in PROBE_OUTPUT.read_text(encoding="ascii").splitlines():
        parts = line.split()
        if len(parts) == 10 and parts[0] in ("a", "b", "c"):
            rows[(parts[0], parts[2], parts[3])] = parts
    return rows


def test_tracked_probe_record_is_labelled_and_shows_the_contrast():
    """Reads the committed record of the probe. It does not rerun the probe."""
    if not PROBE_OUTPUT.exists():
        pytest.skip("the fork-local testbed folder is not part of this checkout")
    text = PROBE_OUTPUT.read_text(encoding="ascii")
    assert text.startswith("LINEAR-ALGEBRA PROBE, NOT A BOUND.")
    assert "NOT A BOUND" in PROBE_SCRIPT.read_text(encoding="utf-8")
    rows = probe_rows()
    assert len(rows) == 12
    for mesh, tets in (("a", 243), ("b", 463), ("c", 949)):
        control = rows[(mesh, "fixed", "continuous")]
        assert int(control[1]) == tets and float(control[9]) >= TOLERANCES["TOL_FEAS_LOCKED"]
        for walls in ("fixed", "roller"):
            assert float(rows[(mesh, walls, "jump")][9]) <= TOLERANCES["TOL_FEAS_OK"]


def test_certificate_demonstration_record_reached_the_procedure():
    """Reads the committed record of the demonstration. It proves reachability, not a bound."""
    if not DEMO_RECORD.exists():
        pytest.skip("the fork-local testbed folder is not part of this checkout")
    record = json.loads(DEMO_RECORD.read_text(encoding="ascii"))
    assert record["title"].startswith("PROCEDURE DEMONSTRATION, NOT A CERTIFIED BOUND")
    assert "NOT THE INDEPENDENT CHECKER" in DEMO_SCRIPT.read_text(encoding="utf-8")
    required = {("tresca_strip_vm_thin_slab", 0), ("upper3d_fixed_phi00", 0)}
    seen = set()
    for run in record["runs"]:
        if run["mesh_label"] == "listed":
            case = ALL_CASES[run["case"]]
            problem, arrays = build_case_problem(case, run["mesh"])
            assert problem["problem_sha256"] == run["problem_sha256"] == load_problem_hashes()[run["case"]][f"m{run['mesh']}"]
            seen.add((run["case"], run["mesh"]))
            for item in (run["lower"], run["upper"]):
                result = item["result"]
                validate_certificate(item["certificate"], problem, arrays, result)
                assert item["certificate"]["classification"] == "certified_bound"
                assert item["certificate"]["checker"]["name"] == PROCEDURE_DEMO_CHECKER
                assert not gate_admissible(item["certificate"]), "the demonstration satisfies no gate"
        else:
            for item in (run["lower"], run["upper"]):
                assert item["certificate"]["classification"] == "certified_bound"
                assert not certified_check_failures(item["certificate"], "ray")
                assert not gate_admissible(item["certificate"])
    assert required <= seen
    assert any(run["upper"]["certificate"]["checks"]["repair_mu"] > 0.0 for run in record["runs"]), "a repair was exercised"


def test_public_files_carry_no_private_material():
    forbidden = ["TI" + "Ms", "Op" + "tum", "Fa" + "ble", "C:" + "\\Users", "nmo" + "ra"]
    extra = (PROBE_SCRIPT, PROBE_OUTPUT, DEMO_SCRIPT, DEMO_RECORD, MESH_SCRIPT)
    paths = [INTERFACE, CASES, Path(__file__)] + [p for p in extra if p.exists()]
    paths += [path for path in DATA.iterdir() if path.suffix == ".json"]
    for path in paths:
        text = path.read_text(encoding="utf-8")
        for word in forbidden:
            assert word not in text, f"{word!r} in {path.name}"
    for path in DATA.glob("*.npz"):
        with np.load(path, allow_pickle=False) as npz:
            assert all(npz[name].dtype.str in ("<f8", "<i8") for name in npz.files), path.name


# Case evaluator ---------------------------------------------------------------------------------

def synthetic_pair(case, mesh, bound, formulation, value, checker="independent-checker"):
    """A result document, the certificate a checker would write for it, valid for the case's problem, and its
    method report (amendment 1)."""
    problem, arrays = case_problem_cached(case, mesh)
    result = {"result_sha256": hashlib.sha256(f"{case['name']}|{mesh}|{bound}|{formulation}|{value!r}".encode()).hexdigest(),
              "problem_sha256": problem["problem_sha256"], "query_id": "q", "bound": bound, "formulation": formulation,
              "value": value}
    certificate = _certificate(problem, arrays, result)
    certificate["checker"] = {"name": checker, "version": "1"}
    return result, certificate, _method_report(certificate)


def synthetic_records(case, low, high, uppers=None, meshes=None, checker="independent-checker", values=None):
    """(certificates, results, feasibility, methods) for every listed mesh; values maps (mesh, bound) to a value."""
    uppers = uppers or [f for f in case["formulation"] if f.startswith("upper.")]
    certificates, results, feasibility, methods = [], [], [], []
    for mesh in (range(len(case["mesh"])) if meshes is None else meshes):
        for bound, formulations, default in (("lower", [LOWER_ID], low), ("upper", uppers, high)):
            for formulation in formulations:
                value = (values or {}).get((mesh, bound), default)
                result, certificate, method = synthetic_pair(case, mesh, bound, formulation, value, checker)
                results.append(result)
                certificates.append(certificate)
                methods.append(method)
                if bound == "upper":
                    record = _feasibility(case_problem_cached(case, mesh)[0], result)
                    record["checker"] = {"name": checker, "version": "1"}
                    feasibility.append(record)
    return certificates, results, feasibility, methods


def synthetic_refusal(case, mesh, result, checker="independent-checker", kind="capacity", stage="kappa_bound"):
    """The refusal a checker would write instead of a certificate (amendment 1, interface 15.8)."""
    problem = case_problem_cached(case, mesh)[0]
    return {
        "format": FORMAT_REFUSAL, "format_version": FORMAT_VERSION, "problem_sha256": problem["problem_sha256"],
        "result_sha256": result["result_sha256"], "query_id": result["query_id"], "bound": result["bound"],
        "formulation": result["formulation"], "kind": kind, "stage": stage,
        "detail": "dense blocks of the reduced system did not fit; no lower bound on sigma_r was proved",
        "resources": {"n_rows": 1000, "n_unknowns": 2000, "nnz": 9000},
        "budget": {"memory_bytes": 6.4e10, "wall_s": 3600.0, "exceeded": "memory" if kind == "capacity" else None},
        "checker": {"name": checker, "version": "1"}, "tolerances": dict(TOLERANCES),
    }


STRIP = "tresca_strip_vm_thin_slab"


def failed(verdict):
    return {(item["kind"], item["cause"]) for item in verdict["items"] if item["status"] == "fail"}


def test_evaluator_passes_a_consistent_certified_strip_without_bracket_records():
    case = ALL_CASES[STRIP]
    verdict = evaluate_case(case, *synthetic_records(case, 5.05, 5.25))
    assert verdict["passed"] and not verdict["stop_program"] and not verdict["p0_blocking"]


def test_evaluator_reads_the_named_gap_limit():
    case = ALL_CASES[STRIP]
    verdict = evaluate_case(case, *synthetic_records(case, 5.0, 5.3))
    assert not verdict["passed"] and ("gap_le", "violated") in failed(verdict) and verdict["p0_blocking"]


def test_evaluator_refuses_estimates_and_reserved_checkers():
    case = ALL_CASES[STRIP]
    certificates, results, feasibility, methods = synthetic_records(case, 5.05, 5.25)
    for doc in certificates:
        doc.update(classification="estimate", reasons=["no_interior_field"], value_certified=None, floor=None)
    verdict = evaluate_case(case, certificates, results, feasibility, methods)
    assert not verdict["passed"] and ("lower_le_target", "estimate") in failed(verdict)
    for checker in sorted(NON_GATE_CHECKERS):
        verdict = evaluate_case(case, *synthetic_records(case, 5.05, 5.25, checker=checker))
        assert not verdict["passed"] and ("certified_bracket", "estimate") in failed(verdict), checker


def test_evaluator_needs_every_listed_mesh():
    case = ALL_CASES[STRIP]
    verdict = evaluate_case(case, *synthetic_records(case, 5.05, 5.25, meshes=[0, 1, 2]))
    assert not verdict["passed"] and ("certified_bracket", "missing") in failed(verdict)


def test_evaluator_refuses_records_of_another_case():
    case = ALL_CASES[STRIP]
    verdict = evaluate_case(case, *synthetic_records(ALL_CASES["tresca_strip_weighted_insensitive"], 5.05, 5.25))
    assert not verdict["passed"] and verdict["foreign_records"]


def test_evaluator_binds_every_formulation_line():
    case = ALL_CASES[STRIP]
    verdict = evaluate_case(case, *synthetic_records(case, 5.05, 5.25, uppers=[CONTROL_ID]))
    assert not verdict["passed"], "a control-only record set does not satisfy the jump formulation"
    assert any(item["formulation"] == JUMP_ID and item["status"] == "fail" for item in verdict["items"])


def test_evaluator_validates_every_certificate_it_uses():
    case = ALL_CASES[STRIP]
    certificates, results, feasibility, methods = synthetic_records(case, 5.05, 5.25)
    certificates[0]["checks"]["omega_eq"] = None
    certificates[1]["tolerances"] = {}
    results.pop()  # the last certificate loses its result document
    verdict = evaluate_case(case, certificates, results, feasibility, methods)
    codes = {item["code"] for item in verdict["invalid_records"]}
    assert {"certificate_checks", "tolerances", "certificate_without_candidate"} <= codes
    assert not verdict["passed"] and verdict["p0_blocking"]


def test_evaluator_validates_feasibility_and_bracket_records():
    case = ALL_CASES["upper3d_fixed_phi00"]
    certificates, results, feasibility, methods = synthetic_records(case, 4.0, 4.5, uppers=[JUMP_ID])
    assert evaluate_case(case, certificates, results, feasibility, methods)["passed"]
    broken = copy.deepcopy(feasibility)
    broken[0]["witness_result_sha256"] = None  # feasible without a witness
    verdict = evaluate_case(case, certificates, results, broken, methods)
    assert {"kind": "feasibility", "id": broken[0]["problem_sha256"], "code": "feasibility_witness"} in verdict["invalid_records"]
    assert not verdict["passed"] and ("kinematic_status", "missing") in failed(verdict)
    lower, upper = certificates[0], certificates[1]
    record = _bracket(lower, upper)
    record["gap"] = 0.0  # a bracket record is checked, never trusted
    verdict = evaluate_case(case, certificates, results, feasibility, methods, brackets=[record])
    assert [item["code"] for item in verdict["invalid_records"]] == ["bracket_gap"] and not verdict["passed"]


def test_evaluator_recomputes_status_and_gap_from_the_certificates():
    case = ALL_CASES[STRIP]
    verdict = evaluate_case(case, *synthetic_records(case, 5.05, 5.25))
    gaps = [item["detail"] for item in verdict["items"] if item["kind"] == "gap_le"]
    assert gaps and all(detail == f"gap {(5.25 - 5.05) / 5.15!r}" for detail in gaps)


def test_evaluator_gate_case_ignores_the_control_and_needs_the_jump_formulation():
    case = ALL_CASES["upper3d_fixed_phi00"]
    certificates, results, feasibility, methods = synthetic_records(case, 4.0, 4.5, uppers=[JUMP_ID])
    verdict = evaluate_case(case, certificates, results, feasibility, methods)
    assert verdict["passed"]
    assert len(verdict["records"]) == 6 and all(item["gap"] is not None for item in verdict["records"])
    control = evaluate_case(case, *synthetic_records(case, 4.0, 4.5, uppers=[CONTROL_ID]))
    assert not control["passed"]
    no_kinematic = evaluate_case(case, certificates, results, methods=methods)
    assert not no_kinematic["passed"] and ("kinematic_status", "missing") in failed(no_kinematic)


def test_evaluator_monotonicity_uses_the_nested_pairs():
    case = ALL_CASES[STRIP]
    verdict = evaluate_case(case, *synthetic_records(case, 5.05, 5.25, values={(1, "upper"): 5.3}))
    assert ("upper_nonincrease", "violated") in failed(verdict)


def test_order_within_floor_counts_as_zero_gap_only_for_exact_cases():
    case = ALL_CASES["vm_plate_shear"]
    target = case["target"]
    verdict = evaluate_case(case, *synthetic_records(case, target * (1.0 + 5e-13), target))
    assert verdict["passed"], failed(verdict)
    strip = ALL_CASES[STRIP]
    verdict = evaluate_case(strip, *synthetic_records(strip, 5.2 * (1.0 + 5e-13), 5.2))
    assert not verdict["passed"] and ("gap_le", "order_within_floor") in failed(verdict)
    assert not verdict["stop_program"]


def test_contradiction_stops_the_program():
    case = ALL_CASES[STRIP]
    verdict = evaluate_case(case, *synthetic_records(case, 5.2, 5.1))
    assert verdict["stop_program"] and not verdict["passed"] and verdict["p0_blocking"]
    assert evaluate_program([verdict])["stop_program"]


def test_column_sandwich_kills_a_zero_gravity_lower_mutant():
    case = ALL_CASES["vm_weighted_column_sandwich"]
    verdict = evaluate_case(case, *synthetic_records(case, S3, 1.45))
    assert not verdict["passed"] and ("lower_le_formula", "violated") in failed(verdict)
    assert evaluate_case(case, *synthetic_records(case, 1.36, 1.3835))["passed"]


PHI30 = "upper3d_fixed_phi30"


def gate_verdict(name, low, high, **changes):
    case = ALL_CASES[name]
    return evaluate_case(case, *synthetic_records(case, low, high, uppers=[JUMP_ID], **changes))


def test_thirty_degree_wide_gap_no_longer_gates_p0():
    """Amendment 3 (OD5): gap 0.545 > MAX_GAP_3D_DP on every mesh passes P0; the gap is only recorded."""
    verdict = gate_verdict(PHI30, 4.0, 7.0)
    assert verdict["passed"] and not verdict["claim_withheld"] and not verdict["p0_blocking"]
    assert not any(item["kind"] == "gap_le" for item in verdict["items"])
    assert [record["mesh"] for record in verdict["records"]] == list(range(6))
    assert all(record["gap"] == pytest.approx(3.0 / 5.5) for record in verdict["records"])


def test_thirty_degree_gap_branch_of_od2_still_withholds():
    """OD2's gap branch has no P0 gate line left (amendment 3); `_withholdable` keeps it, shown on a modified copy."""
    case = copy.deepcopy(ALL_CASES[PHI30])
    case["gate"].append(["gap_le", "MAX_GAP_3D_DP", "mesh", "4"])
    verdict = evaluate_case(case, *synthetic_records(case, 4.0, 7.0, uppers=[JUMP_ID]))
    assert failed(verdict) == {("gap_le", "violated")}
    assert not verdict["passed"] and verdict["claim_withheld"] and not verdict["p0_blocking"]
    program = evaluate_program([verdict])
    assert program["withheld_claims"] == [PHI30] and not program["p0_open"]


@pytest.mark.parametrize("name", ["upper3d_fixed_phi00", "upper3d_roller_phi20", PHI30])
def test_mesh_five_gates_nothing_in_p0(name):
    """Amendment 3: no certificate on mesh 5, or values there that broke the old (3, 5) pair, leave P0 passed."""
    case = ALL_CASES[name]
    certificates, results, feasibility, methods = synthetic_records(case, 4.0, 4.5, uppers=[JUMP_ID])
    sha = load_problem_hashes()[name]["m5"]
    verdict = evaluate_case(case, [doc for doc in certificates if doc["problem_sha256"] != sha],
                            [doc for doc in results if doc["problem_sha256"] != sha], feasibility, methods)
    assert verdict["passed"], failed(verdict)
    assert next(r for r in verdict["records"] if r["mesh"] == 5)["detail"] == "missing"
    loose = gate_verdict(name, 4.0, 4.5, values={(5, "lower"): 3.9, (5, "upper"): 4.6})
    assert loose["passed"], failed(loose)


def test_thirty_degree_missing_certificate_withholds_the_claim():
    case = ALL_CASES[PHI30]
    certificates, results, feasibility, methods = synthetic_records(case, 4.0, 4.5, uppers=[JUMP_ID])
    mesh2 = load_problem_hashes()[PHI30]["m2"]
    certificates = [doc for doc in certificates if doc["problem_sha256"] != mesh2]
    results = [doc for doc in results if doc["problem_sha256"] != mesh2]  # no result offered: truly missing
    verdict = evaluate_case(case, certificates, results, feasibility, methods)
    assert failed(verdict) == {("certified_upper", "missing"), ("certified_bracket", "missing")}  # amendment 3
    assert verdict["claim_withheld"] and not verdict["p0_blocking"]
    no_feasibility = [doc for doc in feasibility if doc["problem_sha256"] != mesh2]
    verdict = evaluate_case(case, certificates, results, no_feasibility, methods)
    assert not verdict["claim_withheld"] and verdict["p0_blocking"], "a missing kinematic record is not a missing certificate"


@pytest.mark.parametrize("values,kind", [
    ({(4, "upper"): 4.6}, "upper_nonincrease"),
    ({(4, "lower"): 3.9}, "lower_nondecrease"),  # amendment 3: the pair is (0, 4); (3, 5) is not checked
])
def test_thirty_degree_monotonicity_failure_keeps_p0_open(values, kind):
    verdict = gate_verdict(PHI30, 4.0, 4.5, values=values)
    assert (kind, "violated") in failed(verdict)
    assert not verdict["claim_withheld"] and verdict["p0_blocking"]
    assert evaluate_program([verdict])["p0_open"]


def test_thirty_degree_contradiction_stops_the_program():
    """A contradiction stops the program even on mesh 5, which gates nothing since amendment 3: its record finds it."""
    verdict = gate_verdict(PHI30, 4.0, 4.5, values={(5, "lower"): 4.6})
    assert verdict["stop_program"] and not verdict["claim_withheld"]
    assert evaluate_program([verdict])["stop_program"]


def test_thirty_degree_kinematic_failure_is_not_withheld():
    case = ALL_CASES[PHI30]
    certificates, results, _, methods = synthetic_records(case, 4.0, 4.5, uppers=[JUMP_ID])
    verdict = evaluate_case(case, certificates, results, methods=methods)
    assert not verdict["claim_withheld"] and verdict["p0_blocking"]


def test_twenty_degree_wide_gap_passes_p0_and_a_monotonicity_failure_keeps_it_open():
    """Amendment 3: a gap above MAX_GAP_3D_DP no longer fails a 20-degree case; a (0, 4) violation still does."""
    verdict = gate_verdict("upper3d_fixed_phi20", 4.0, 7.0)
    assert verdict["passed"] and not verdict["p0_blocking"]
    verdict = gate_verdict("upper3d_fixed_phi20", 4.0, 4.5, values={(4, "upper"): 4.6})
    assert ("upper_nonincrease", "violated") in failed(verdict)
    assert not verdict["claim_withheld"] and verdict["p0_blocking"]


def test_program_needs_every_required_case():
    verdict = gate_verdict(PHI30, 4.0, 4.5)
    assert verdict["passed"]
    program = evaluate_program([verdict], required=[PHI30, STRIP])
    assert program["p0_open"] and program["missing_cases"] == [STRIP]


# Amendment 1 -------------------------------------------------------------------------------------

def test_interface_records_the_freeze_and_amendment_1():
    text = interface_text()
    assert "**frozen at `03f59c1a4`**" in text[:600]
    log = text[text.index("### 19.1 Amendment 1"):]
    assert "`03f59c1a4`" in log
    for item in range(1, 11):
        assert f"\n| {item} |" in log, item
    assert "the 2026-11-19 review-fix round was recorded here" in text, "earlier record intact (wording: amendment 2, m7)"
    for term in ("certificate_sha256", "rows_skipped", "interior_field", "D_scale_vol", "D_scale_jump", "n_S",
                 "refused", "refusal_conflict", "method_report_missing"):
        assert f"`{term}`" in text, term


@pytest.mark.parametrize("name,mesh,label", [
    ("tresca_strip_vm_thin_slab", 0, "8x1x4"), ("vm_plate_shear", 1, "2x2x2"),
    ("upper3d_fixed_phi00", 0, "upper3d_a+0"), ("upper3d_fixed_phi00", 5, "upper3d_d+1"),
    ("upper3d_roller_phi30", 4, "upper3d_phi30_a+1"),
])
def test_case_problem_name_rule(name, mesh, label):
    """Amendment 1, item 4: the hashed problem name is <case>:<label> (interface 15.2)."""
    case = ALL_CASES[name]
    problem, _ = case_problem_cached(case, mesh)
    assert problem["name"] == f"{name}:{label}"
    assert problem["problem_sha256"] == load_problem_hashes()[name][f"m{mesh}"]


def test_upper_objective_is_c_transpose_x_and_value_is_the_field_value():
    """Amendment 1, item 2: they coincide only when W1 = 1."""
    problem, arrays = halfspace_problem()
    result, fields = load_document(DATA / f"{record_stem('vm_ground_patch_halfspace', JUMP_ID)}.result.json")
    dissipation, w0, w1 = fixture_upper_terms(problem, arrays, fields["element_velocity"], fields["footing_velocity"])
    assert result["value"] == pytest.approx((dissipation - w0) / w1, rel=1e-14)
    assert result["objective"] == pytest.approx(dissipation - w0, rel=1e-14)
    assert abs(w1 - 1.0) > 1e-3 and result["objective"] != pytest.approx(result["value"])
    lower = load_set("vm_plate_shear")["lower"]
    assert lower["value"] == -lower["objective"], "Lower: value = objective_sign * objective"


def method_case(bound="upper"):
    s = load_set("vm_plate_shear")
    certificate = copy.deepcopy(s[f"{bound}_cert"])
    certificate["checker"] = {"name": "independent-checker", "version": "1"}
    return s, certificate, _method_report(certificate)


def test_method_report_is_bound_to_its_certificate():
    s, certificate, method = method_case()
    validate_method_report(method, certificate, s["problem"], s["upper"], s["upper_arrays"], s["arrays"])
    assert gate_admissible(certificate, method=method)
    certificate["value_certified"] *= 1.0 + 1e-15
    with pytest.raises(ContractError) as error:
        validate_method_report(method, certificate)
    assert error.value.code == "method_certificate_hash"
    assert not gate_admissible(certificate, method=method)


def _kappa(**changes):
    def mutate(method, certificate):
        method["kappa"].update(changes)
    return mutate


def _proofs(proofs, skipped, rank=None):
    def mutate(method, certificate):
        method["dependent_rows"]["proofs"] = proofs
        method["rows_skipped"] = skipped
        method["dependent_rows"]["rank"] = method["dependent_rows"]["n_rows"] - skipped if rank is None else rank
    return mutate


def _interior(source, sha=None, repair=None):
    def mutate(method, certificate):
        method["interior_field"] = {"source": source, "sha256": sha}
        if repair is not None:
            certificate["checks"]["repair_mu"] = repair
            method["certificate_sha256"] = certificate_hash(certificate)
    return mutate


def _set_kappa_eq(value):
    def mutate(method, certificate):
        certificate["checks"]["kappa_eq"] = value
        certificate["floor"] = max(TOLERANCES["FLOOR_MIN_REL"], value * certificate["checks"]["omega_eq"]) * abs(
            certificate["value_certified"])
        method["certificate_sha256"] = certificate_hash(certificate)
    return mutate


METHOD_MUTATIONS = [
    ("estimated_kappa", _kappa(kappa_upper=0.5), "method_kappa_not_proved"),
    ("kappa_not_the_certificate_one", _set_kappa_eq(2.0), "method_kappa_mismatch"),
    ("zero_lower_bound", _kappa(sigma_min_lower=0.0), "method_kappa"),
    ("unknown_method", _kappa(method="lanczos_estimate"), "method_kappa"),
    ("eigen_without_residual", _kappa(method="eigen_normal_columns"), "method_kappa"),
    ("svd_with_completeness", _kappa(completeness={"route": "inertia", "shift": 1.0, "count_below": 0}),
     "method_kappa"),
    ("svd_with_residual", _kappa(residual=1e-14), "method_kappa"),
    ("skipped_without_proof", _proofs([], 2), "method_rows_skipped"),
    ("proof_count_mismatch", _proofs([{"route": "shape_forced_zeros", "count": 1, "residual_bound": None}], 2),
     "method_rows_skipped"),
    ("rank_mismatch", _proofs([{"route": "shape_forced_zeros", "count": 1, "residual_bound": None}], 1, rank=4),
     "method_rows_skipped"),
    ("unknown_route", _proofs([{"route": "threshold", "count": 1, "residual_bound": None}], 1), "method_dependent_rows"),
    ("null_space_without_bound", _proofs([{"route": "null_space_basis", "count": 1, "residual_bound": None}], 1),
     "method_dependent_rows"),
    ("null_space_residual_too_large", _proofs([{"route": "null_space_basis", "count": 1, "residual_bound": 1e-3}], 1),
     "method_dependent_rows"),
    ("repair_without_interior", _interior("none", repair=1e-7), "method_interior"),
    ("zero_interior_for_upper", _interior("zero"), "method_interior"),
    ("result_interior_without_hash", _interior("result"), "method_interior"),
    ("result_interior_not_in_result", _interior("result", "0" * 64), "method_interior"),
]


@pytest.mark.parametrize("label,mutate,code", METHOD_MUTATIONS, ids=[m[0] for m in METHOD_MUTATIONS])
def test_malformed_method_report_is_refused(label, mutate, code):
    s, certificate, method = method_case()
    mutate(method, certificate)
    with pytest.raises(ContractError) as error:
        validate_method_report(method, certificate, s["problem"], s["upper"], s["upper_arrays"], s["arrays"])
    assert error.value.code == code
    if label != "result_interior_not_in_result":  # needs the result arrays, which gate_admissible does not see
        assert not gate_admissible(certificate, method=method)


def proved_eigen_case():
    """Three proved dependent rows and an eigenvalue route with a counted inertia below the shift."""
    s, certificate, method = method_case()
    method["dependent_rows"]["proofs"] = [
        {"route": "shape_forced_zeros", "count": 1, "residual_bound": None},
        {"route": "dense_svd_blocks", "count": 1, "residual_bound": 1e-14},
        {"route": "null_space_basis", "count": 1, "residual_bound": 1e-12},
    ]
    method["rows_skipped"] = 3
    method["dependent_rows"]["rank"] = 1
    # A_hat^T A_hat has n_unknowns - rank = 9 zero eigenvalues; nine lie below the shift 1.0, and 1.0^2 <= 1.0 - 0.
    method["kappa"].update(method="eigen_normal_columns", residual=1e-15, sigma_min_lower=0.999,
                           kappa_upper=certificate["checks"]["kappa_eq"],
                           completeness={"route": "inertia", "shift": 1.0, "count_below": 9})
    method["kappa"]["sigma_max_upper"] = method["kappa"]["kappa_upper"] * method["kappa"]["sigma_min_lower"]
    return s, certificate, method


def test_unproved_eigenvalue_report_is_refused():
    """M1: a Ritz value with a residual but no completeness proof bounds nothing."""
    s, certificate, method = proved_eigen_case()
    method["kappa"]["completeness"] = None
    with pytest.raises(ContractError) as error:
        validate_method_report(method, certificate, s["problem"], s["upper"], s["upper_arrays"], s["arrays"])
    assert error.value.code == "method_kappa_completeness"
    assert not gate_admissible(certificate, method=method)


def test_proved_dependent_rows_and_an_eigenvalue_route_are_accepted():
    s, certificate, method = proved_eigen_case()
    validate_method_report(method, certificate, s["problem"], s["upper"], s["upper_arrays"], s["arrays"])
    assert gate_admissible(certificate, method=method)
    method["kappa"]["completeness"]["route"] = "full_spectrum"
    validate_method_report(method, certificate)


def _completeness(**changes):
    def mutate(method, certificate):
        method["kappa"]["completeness"].update(changes)
    return mutate


@pytest.mark.parametrize("mutate,code", [
    (_completeness(route="lanczos"), "method_kappa_completeness"),
    (_completeness(count_below=8), "method_kappa_completeness"),
    (_completeness(shift=1e-16), "method_kappa_completeness"),
    (_completeness(shift=0.5), "method_kappa_not_proved"),
], ids=["unknown_route", "count_not_the_zeros", "shift_below_residual", "lower_bound_above_shift"])
def test_eigenvalue_completeness_must_prove_the_bound(mutate, code):
    s, certificate, method = proved_eigen_case()
    mutate(method, certificate)
    with pytest.raises(ContractError) as error:
        validate_method_report(method, certificate)
    assert error.value.code == code


@pytest.mark.parametrize("kappa,sigma_min,bound", [(1.0, 1.0, 1e-9), (3e10, 5e-11, 1e-11)],
                         ids=["perturbation_above_threshold", "kept_value_not_separated"])
def test_dense_block_dependent_rows_need_a_separated_split(kappa, sigma_min, bound):
    """M1: residual_bound <= TOL_RANK sigma_max_upper and sigma_min_lower > TOL_RANK sigma_max_upper + residual_bound."""
    s, certificate, method = method_case()
    certificate["checks"]["kappa_eq"] = kappa
    method["certificate_sha256"] = certificate_hash(certificate)
    method["kappa"].update(sigma_max_upper=1.0, sigma_min_lower=sigma_min, kappa_upper=kappa)
    method["dependent_rows"]["proofs"] = [{"route": "dense_svd_blocks", "count": 1, "residual_bound": bound}]
    method["rows_skipped"] = 1
    method["dependent_rows"]["rank"] = method["dependent_rows"]["n_rows"] - 1
    with pytest.raises(ContractError) as error:
        validate_method_report(method, certificate)
    assert error.value.code == "method_dependent_rows"


def test_result_interior_field_is_identified_by_hash():
    s, certificate, method = method_case()
    arrays = {k: v.copy() for k, v in s["upper_arrays"].items()}
    arrays["interior_velocity"] = 0.5 * arrays["velocity"]
    arrays["interior_footing_velocity"] = 0.5 * arrays["footing_velocity"]
    method["interior_field"] = {"source": "result", "sha256": interior_field_sha256("upper", CONTROL_ID, arrays)}
    validate_method_report(method, certificate, s["problem"], None, arrays, s["arrays"])
    arrays["interior_velocity"][0, 0] += 1e-9
    with pytest.raises(ContractError) as error:
        validate_method_report(method, certificate, s["problem"], None, arrays, s["arrays"])
    assert error.value.code == "method_interior"


def test_zero_interior_field_needs_a_lower_bound_without_fixed_loads():
    s, certificate, method = method_case("lower")
    method["interior_field"]["source"] = "zero"
    validate_method_report(method, certificate, s["problem"], s["lower"], s["lower_arrays"], s["arrays"])
    problem = copy.deepcopy(s["problem"])
    problem["regions"][0]["b0"] = [0.0, -1.0, 0.0]
    with pytest.raises(ContractError) as error:
        validate_method_report(method, certificate, problem)
    assert error.value.code == "method_interior"


def test_bracket_of_gate_certificates_needs_both_method_reports():
    doc, low, high = bracket_set()
    for certificate in (low, high):
        certificate["checker"] = {"name": "independent-checker", "version": "1"}
    expect_bracket("method_report_missing", doc, low, high)
    with pytest.raises(ContractError) as error:
        validate_bracket(doc, low, high, _method_report(low), None)
    assert error.value.code == "method_report_missing"
    broken = _method_report(high)
    broken["rows_skipped"] = 1
    with pytest.raises(ContractError) as error:
        validate_bracket(doc, low, high, _method_report(low), broken)
    assert error.value.code == "method_rows_skipped"
    validate_bracket(doc, low, high, _method_report(low), _method_report(high))


def test_evaluator_needs_a_valid_method_report_for_every_gate_certificate():
    case = ALL_CASES[STRIP]
    certificates, results, feasibility, methods = synthetic_records(case, 5.05, 5.25)
    verdict = evaluate_case(case, certificates, results, feasibility)
    assert not verdict["passed"] and ("certified_bracket", "unproved") in failed(verdict)
    methods[0]["kappa"]["kappa_upper"] = 0.5  # an estimate below the ratio it should bound
    verdict = evaluate_case(case, certificates, results, feasibility, methods)
    assert {"kind": "method", "id": methods[0]["result_sha256"], "code": "method_kappa_not_proved"} in verdict["invalid_records"]
    assert not verdict["passed"] and verdict["p0_blocking"]


def refused_records(name, mesh, bound="lower"):
    case = ALL_CASES[name]
    certificates, results, feasibility, methods = synthetic_records(case, 4.0, 4.5, uppers=[JUMP_ID])
    sha = load_problem_hashes()[name][f"m{mesh}"]
    keep = [doc for doc in certificates if not (doc["problem_sha256"] == sha and doc["bound"] == bound)]
    gone = next(doc for doc in certificates if doc["problem_sha256"] == sha and doc["bound"] == bound)
    result = next(doc for doc in results if doc["result_sha256"] == gone["result_sha256"])
    return case, keep, results, feasibility, methods, synthetic_refusal(case, mesh, result), gone


def without_result(results, certificate):
    return [doc for doc in results if doc["result_sha256"] != certificate["result_sha256"]]


def test_refusal_record_validates_and_is_not_a_certificate():
    case, _, results, _, _, refusal, gone = refused_records(PHI30, 5)
    problem = case_problem_cached(case, 5)[0]
    result = next(doc for doc in results if doc["result_sha256"] == gone["result_sha256"])
    validate_refusal(refusal, problem, result)
    for key, value, code in (("kind", "estimate", "refusal_kind"), ("stage", "solve", "refusal_stage"),
                             ("detail", "", "refusal_detail"), ("bound", "upper", "refusal_identity"),
                             ("tolerances", {}, "tolerances"), ("value", 4.0, "unknown_key")):
        broken = copy.deepcopy(refusal)
        broken[key] = value
        with pytest.raises(ContractError) as error:
            validate_refusal(broken, problem, result)
        assert error.value.code == code, key
    broken = copy.deepcopy(refusal)
    broken["resources"]["nnz"] = -1
    with pytest.raises(ContractError) as error:
        validate_refusal(broken, problem, result)
    assert error.value.code == "refusal_resources"


def test_thirty_degree_refusal_is_not_a_missing_certificate():
    """Amendment 1, item 7: a capacity refusal on a gate mesh keeps P0 open, where a missing certificate would withhold.
    Amendment 3: mesh 4 is the finest gate mesh, so the refusal sits there."""
    case, certificates, results, feasibility, methods, refusal, gone = refused_records(PHI30, 4)
    missing = evaluate_case(case, certificates, without_result(results, gone), feasibility, methods)
    assert missing["claim_withheld"] and not missing["p0_blocking"], "with no result offered OD2 withholds"
    verdict = evaluate_case(case, certificates, results, feasibility, methods, refusals=[refusal])
    causes = failed(verdict)
    assert ("certified_bracket", "refused") in causes and ("lower_nondecrease", "refused") in causes
    assert not any(cause == "missing" for _, cause in causes)
    assert not verdict["claim_withheld"] and verdict["p0_blocking"] and not verdict["stop_program"]
    program = evaluate_program([verdict])
    assert program["p0_open"] and program["refused_cases"] == [PHI30] and not program["withheld_claims"]


def test_refusal_from_a_reserved_checker_or_beside_a_certificate_does_not_count():
    case, certificates, results, feasibility, methods, refusal, gone = refused_records(PHI30, 4)
    demo = copy.deepcopy(refusal)
    demo["checker"]["name"] = PROCEDURE_DEMO_CHECKER
    verdict = evaluate_case(case, certificates, results, feasibility, methods, refusals=[demo])
    assert ("certified_bracket", "unchecked") in failed(verdict), "a demonstration refusal is not the checker's"
    assert not verdict["claim_withheld"] and verdict["p0_blocking"]
    verdict = evaluate_case(case, certificates + [gone], results, feasibility, methods, refusals=[refusal])
    assert [item["code"] for item in verdict["invalid_records"]] == ["refusal_conflict"] and not verdict["passed"]


def test_von_mises_refusal_blocks_on_mesh_four_and_not_on_mesh_five():
    """Amendment 1, item 8, read under amendment 3: a refused jump Upper on gate mesh 4 fails its gate items; on mesh 5,
    which gates nothing in P0, it fails none and shows only in the gap record."""
    case, certificates, results, feasibility, methods, refusal, _ = refused_records("upper3d_fixed_phi00", 4, "upper")
    verdict = evaluate_case(case, certificates, results, feasibility, methods, refusals=[refusal])
    kinds = {kind for kind, cause in failed(verdict) if cause == "refused"}
    assert {"certified_upper", "certified_bracket", "upper_nonincrease"} <= kinds
    assert verdict["p0_blocking"]
    case, certificates, results, feasibility, methods, refusal, _ = refused_records("upper3d_fixed_phi00", 5, "upper")
    verdict = evaluate_case(case, certificates, results, feasibility, methods, refusals=[refusal])
    assert verdict["passed"] and not verdict["p0_blocking"], failed(verdict)
    assert next(r for r in verdict["records"] if r["mesh"] == 5)["detail"] == "refused"
    assert evaluate_program([verdict])["refused_cases"] == [], "a mesh-5 refusal is evidence, not a refused item"


def test_thirty_degree_gate_mesh_result_without_certificate_or_refusal_is_unchecked():
    """B1: an offered result that no checker answered is not a missing certificate (mesh 4 since amendment 3)."""
    case, certificates, results, feasibility, methods, _, gone = refused_records(PHI30, 4)
    verdict = evaluate_case(case, certificates, results, feasibility, methods)
    causes = failed(verdict)
    assert ("certified_bracket", "unchecked") in causes and ("lower_nondecrease", "unchecked") in causes
    assert not any(cause == "missing" for _, cause in causes)
    assert not verdict["claim_withheld"] and verdict["p0_blocking"]
    missing = evaluate_case(case, certificates, without_result(results, gone), feasibility, methods)
    assert ("certified_bracket", "missing") in failed(missing) and missing["claim_withheld"]


def test_thirty_degree_gate_certificate_without_method_report_is_unproved():
    """B1: a gate certificate whose method report is absent or invalid is unproved, never withheld."""
    case = ALL_CASES[PHI30]
    certificates, results, feasibility, methods = synthetic_records(case, 4.0, 4.5, uppers=[JUMP_ID])
    sha = load_problem_hashes()[PHI30]["m4"]
    lower = next(doc for doc in certificates if doc["problem_sha256"] == sha and doc["bound"] == "lower")
    kept = [m for m in methods if m["certificate_sha256"] != certificate_hash(lower)]
    verdict = evaluate_case(case, certificates, results, feasibility, kept)
    causes = failed(verdict)
    assert ("certified_bracket", "unproved") in causes and ("lower_nondecrease", "unproved") in causes
    assert not verdict["claim_withheld"] and verdict["p0_blocking"]


def test_thirty_degree_ill_conditioned_estimate_is_not_withheld():
    """B1: an estimate refused for conditioning keeps P0 open; another estimate at 30 degrees is withheld."""
    case = ALL_CASES[PHI30]
    certificates, results, feasibility, methods = synthetic_records(case, 4.0, 4.5, uppers=[JUMP_ID])
    sha = load_problem_hashes()[PHI30]["m4"]
    lower = next(doc for doc in certificates if doc["problem_sha256"] == sha and doc["bound"] == "lower")
    for reason, cause, withheld in (("equality_ill_conditioned", "ill_conditioned", False),
                                    ("no_interior_field", "estimate", True)):
        changed = copy.deepcopy(certificates)
        target = next(doc for doc in changed if doc["result_sha256"] == lower["result_sha256"])
        target.update(classification="estimate", reasons=[reason], value_certified=None, floor=None)
        verdict = evaluate_case(case, changed, results, feasibility, methods)
        assert ("certified_bracket", cause) in failed(verdict), reason
        assert verdict["claim_withheld"] is withheld and verdict["p0_blocking"] is not withheld, reason


def test_feasibility_without_a_footing_has_no_footing_modes():
    s = load_set("vm_uniaxial_tension")
    assert not s["problem"]["footings"] and s["feasibility"]["footing_modes"] is None
    validate_feasibility(s["feasibility"], s["problem"])
    for value in (0, 1):
        doc = copy.deepcopy(s["feasibility"])
        doc["footing_modes"] = value
        with pytest.raises(ContractError) as error:
            validate_feasibility(doc, s["problem"])
        assert error.value.code == "feasibility_rank"
    doc = copy.deepcopy(s["feasibility"])
    doc.update(rank=None, nullity=None)
    validate_feasibility(doc, s["problem"])
    plate = load_set("vm_plate_shear")
    doc = copy.deepcopy(plate["feasibility"])
    doc["footing_modes"] = None  # with a footing and a computed rank, footing_modes is computed too
    with pytest.raises(ContractError) as error:
        validate_feasibility(doc, plate["problem"])
    assert error.value.code == "feasibility_rank"
    doc.update(rank=None, nullity=None)
    validate_feasibility(doc, plate["problem"])


# Amendment 2 -------------------------------------------------------------------------------------

def test_interface_records_amendment_2():
    text = interface_text()
    assert "**amendment 2**" in text[:900]
    log = text[text.index("### 19.2 Amendment 2"):]
    for item in ("A", "m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8"):
        assert f"\n| {item} |" in log, item
    assert "The file is not frozen, so" not in text, "m7: the stale wording is gone"
    for term in ("apex_zero", "apex_zero_dual", "z_forced", "apex_zero_lower", "apex_zero_invalid", "method_apex_zero",
                 "witness_projection_move", "witness_work_reference", "budget", "refusal_budget", "result_objective",
                 "interior_field_invalid", "TOL_REPAIR_MAX"):
        assert f"`{term}`" in text, term
    for value in sorted(Z_FORCED):
        assert f"`{value}`" in log, value
    mesh5 = text.index("**Mesh-5 conditioning route")
    assert text.index("### 18.3") < mesh5 < text.index("## 19."), "m3: the capacity statement is evidence (18.3)"
    assert "K1" in text[mesh5:mesh5 + 2500]


# Amendment 3 -------------------------------------------------------------------------------------

def test_interface_records_amendment_3():
    """Owner decision OD5: mesh 5 out of the P0 3D gate, decided after results and disclosed as such."""
    text = interface_text()
    assert "**amendment 3** (section 19.3)" in text[:900]
    log = text[text.index("### 19.3 Amendment 3, 2026-11-19"):]
    for item in ("A", "B", "C", "D", "E", "F"):
        assert f"\n| {item} |" in log, item
    assert "**Disclosure: decided after results.**" in log and "not presented as pre-registered" in log
    owner = text[text.index("### 18.1"):text.index("### 18.2")]
    assert "| OD5 | 2026-11-19, after mesh-4 and mesh-5 results were seen" in owner
    quality = text[text.index("### 11.2 Quality test"):text.index("### 11.3")]
    assert "They are closed from this\nrevision: a limit does not move" in quality, "OD1 sentence kept"
    assert "**That closure did not hold for the true-3D gate.**" in quality
    evidence = text[text.index("### 18.3"):text.index("## 19.")]
    assert "| Amendment 3, 2026-11-19 |" in evidence
    for number in ("9.24394", "9.23579", "3.34843", "3.34841", "5.92866", "5.92915", "10,808", "9,475", "`17ee1b4fc`"):
        assert number in evidence, number
    mesh5 = evidence[evidence.index("**Mesh-5 conditioning route"):]
    assert "cannot pass until WP-237d" not in mesh5 and "no longer exist" in mesh5


def free_apex_jump():
    problem, problem_arrays = load_document(DATA / f"{case_stem('dp_free_hydrostatic')}.problem.json")
    result, arrays = load_document(DATA / f"{record_stem('dp_free_hydrostatic', JUMP_ID)}.result.json")
    return problem, problem_arrays, result, arrays


def with_face(result, arrays, rows, dual=None, rehash=True):
    result, arrays = copy.deepcopy(result), {k: v.copy() for k, v in arrays.items()}
    arrays["apex_zero"] = np.array(rows, dtype="<i8").reshape(-1, 2)
    if dual is not None:
        arrays["apex_zero_dual"] = np.array(dual, dtype="<f8")
    if rehash:
        result["result_sha256"] = document_hash(result, arrays, "result_sha256")
    return result, arrays


def test_declared_apex_face_validates_and_is_hashed():
    problem, problem_arrays, result, arrays = free_apex_jump()
    n_e, n_f = len(problem_arrays["tets"]), len(mesh_facets(problem_arrays["tets"])[0])
    rows = [[0, 0], [0, 4 * n_e - 1], [1, 0], [1, 6 * n_f - 1]]
    faced, faced_arrays = with_face(result, arrays, rows, dual=[[1, 1, 1, 0, 0, 0]] * 2 + [[1, 0, 0, 0, 0, 0]] * 2)
    validate_result(faced, faced_arrays, problem, problem_arrays)
    assert faced["result_sha256"] != result["result_sha256"], "Z is part of the result identity"
    stale, stale_arrays = with_face(result, arrays, rows, rehash=False)
    with pytest.raises(ContractError) as error:
        validate_result(stale, stale_arrays, problem, problem_arrays)
    assert error.value.code == "hash_mismatch"


@pytest.mark.parametrize("rows,dual", [
    ([[0, 4], [0, 0]], None),  # unsorted
    ([[0, 0], [0, 0]], None),  # duplicate
    (np.zeros((0, 2)), None),  # empty: an absent face is no array, so the identity is unique
    ([[2, 0]], None),  # unknown kind
    ([[0, 24]], None),  # flow index past 4 E
    ([[1, 36]], None),  # jump index past 6 F
    ([[0, 0]], [[1.0, 0, 0, 0, 0]]),  # dual of the wrong shape
    ([[1, 0]], [[1.0, 0, 0, 1.0, 0, 0]]),  # a jump dual has three components
], ids=["unsorted", "duplicate", "empty", "kind", "flow_index", "jump_index", "dual_shape", "jump_dual"])
def test_malformed_apex_face_is_refused(rows, dual):
    problem, problem_arrays, result, arrays = free_apex_jump()
    faced, faced_arrays = with_face(result, arrays, rows, dual)
    with pytest.raises(ContractError) as error:
        validate_result(faced, faced_arrays, problem, problem_arrays)
    assert error.value.code in ("apex_zero_invalid", "array_shape")


def test_apex_face_needs_inequality_conditions_of_its_formulation():
    plate = load_set("vm_plate_shear")
    jump, jump_arrays = load_document(DATA / f"{record_stem('vm_plate_shear', JUMP_ID)}.result.json")
    control = load_document(DATA / f"{record_stem('dp_uniaxial_tension', CONTROL_ID)}.result.json")
    tension = load_document(DATA / f"{case_stem('dp_uniaxial_tension')}.problem.json")
    for (result, arrays), (problem, problem_arrays), rows in (
            ((jump, jump_arrays), (plate["problem"], plate["arrays"]), [[0, 0]]),  # von Mises: alpha = 0
            ((jump, jump_arrays), (plate["problem"], plate["arrays"]), [[1, 0]]),  # von Mises: tau = 0
            (control, tension, [[1, 0]])):  # no jump condition in the control
        faced, faced_arrays = with_face(result, arrays, rows)
        with pytest.raises(ContractError) as error:
            validate_result(faced, faced_arrays, problem, problem_arrays)
        assert error.value.code == "apex_zero_invalid", rows
    faced, faced_arrays = with_face(*control, [[0, 0]])
    validate_result(faced, faced_arrays, *tension)  # a flow condition of the control may be declared
    dual_only = {**control[1], "apex_zero_dual": np.zeros((1, 6))}
    with pytest.raises(ContractError) as error:
        validate_result(control[0], dual_only, *tension)
    assert error.value.code == "apex_zero_invalid"


def test_apex_face_is_refused_on_a_lower():
    s = load_set("vm_plate_shear")
    faced, faced_arrays = with_face(s["lower"], s["lower_arrays"], [[0, 0]])
    with pytest.raises(ContractError) as error:
        validate_result(faced, faced_arrays, s["problem"], s["arrays"])
    assert error.value.code == "apex_zero_lower"
    _, certificate, method = method_case("lower")
    method["apex_zero"] = {"count": 1, "sha256": "0" * 64, "z_forced": "absent"}
    with pytest.raises(ContractError) as error:
        validate_method_report(method, certificate)
    assert error.value.code == "method_apex_zero"


def faced_method(z_forced="absent", dual=None):
    """A jump result with Z, its schema certificate under a gate checker name, and a matching method report."""
    problem, problem_arrays, result, arrays = free_apex_jump()
    faced, faced_arrays = with_face(result, arrays, [[0, 0], [1, 0]], dual)
    certificate = _certificate(problem, problem_arrays, faced)
    certificate["checker"] = {"name": "independent-checker", "version": "1"}
    certificate["checks"]["rho_apex"] = 1e-15  # M3
    method = _method_report(certificate)
    method["apex_zero"] = {"count": 2, "sha256": apex_zero_sha256(faced_arrays), "z_forced": z_forced}
    return problem, problem_arrays, faced, faced_arrays, certificate, method


def test_rho_apex_is_reported_and_bounded():
    """Amendment 2, M3: a certified Z-declaring result reports rho_apex <= TOL_CONE_MARGIN; without Z it is null."""
    problem, problem_arrays, faced, faced_arrays, certificate, method = faced_method()
    assert not certified_check_failures(certificate, "ray")
    for value in (2e-12, -1.0, None):
        broken = copy.deepcopy(certificate)
        broken["checks"]["rho_apex"] = value
        if value is None:
            report = copy.deepcopy(method)
            report["certificate_sha256"] = certificate_hash(broken)
            with pytest.raises(ContractError) as error:
                validate_method_report(report, broken, problem, None, faced_arrays, problem_arrays)
            assert error.value.code == "method_apex_zero"
        else:
            assert "rho_apex" in certified_check_failures(broken, "ray"), value
    _, plain, plain_method = method_case()
    plain["checks"]["rho_apex"] = 1e-15
    plain_method["certificate_sha256"] = certificate_hash(plain)
    with pytest.raises(ContractError) as error:
        validate_method_report(plain_method, plain)
    assert error.value.code == "method_apex_zero"


@pytest.mark.parametrize("rows", [[[0, 0]], [[0, 10 ** 9]], [[0, 1], [0, 0]]],
                         ids=["von_mises_element", "index_out_of_range", "unsorted"])
def test_gate_path_validates_the_declared_face(rows):
    """Amendment 2, M2: supplied result arrays go through validate_result, so a malformed Z never reaches a gate."""
    name = STRIP if rows == [[0, 0]] else PHI30  # the strip is von Mises: no flow condition can be declared
    case = ALL_CASES[name]
    uppers = None if name == STRIP else [JUMP_ID]
    low, high = (5.05, 5.25) if name == STRIP else (4.0, 4.5)
    real = real_upper_records(case, 0, high, rows=rows)
    records = swap_mesh_upper(case, synthetic_records(case, low, high, uppers=uppers), 0, real)
    verdict = evaluate_case(case, *records, result_arrays={real[0]["result_sha256"]: real[1]})
    assert {"kind": "result_arrays", "id": real[0]["result_sha256"], "code": "apex_zero_invalid"} in verdict["invalid_records"]
    assert not verdict["passed"] and verdict["p0_blocking"]


def test_method_report_names_the_declared_face():
    problem, problem_arrays, faced, faced_arrays, certificate, method = faced_method()
    validate_method_report(method, certificate, problem, faced, faced_arrays, problem_arrays)
    assert gate_admissible(certificate, method=method, result_arrays=faced_arrays)
    assert not gate_admissible(certificate, method=method), "m1: Z is verified against the result arrays only"
    dual = [[1.0, 0, 0, 0, 0, 0], [1.0, 0, 0, 0, 0, 0]]
    for change, arrays in (({"count": 1}, faced_arrays), ({"sha256": "0" * 64}, faced_arrays),
                           ({"z_forced": "maybe"}, faced_arrays), ({"z_forced": "proved"}, faced_arrays),
                           ({}, with_face(faced, faced_arrays, [[0, 0], [1, 0]], dual)[1])):
        broken = copy.deepcopy(method)
        broken["apex_zero"].update(change)
        with pytest.raises(ContractError) as error:
            validate_method_report(broken, certificate, problem, None, arrays, problem_arrays)
        assert error.value.code == "method_apex_zero", change
    silent = copy.deepcopy(method)
    silent["apex_zero"] = None
    with pytest.raises(ContractError) as error:
        validate_method_report(silent, certificate, problem, None, faced_arrays, problem_arrays)
    assert error.value.code == "method_apex_zero"


@pytest.mark.parametrize("z_forced", ["proved", "unproved"])
def test_z_forced_is_a_diagnostic_never_a_gate(z_forced):
    """A wrong Z only makes the Upper worse, never unsound, so `proved` and `unproved` are both admissible."""
    dual = [[1.0, 0, 0, 0, 0, 0], [1.0, 0, 0, 0, 0, 0]]
    problem, problem_arrays, faced, faced_arrays, certificate, method = faced_method(z_forced, dual)
    validate_method_report(method, certificate, problem, faced, faced_arrays, problem_arrays)
    assert gate_admissible(certificate, method=method, result_arrays=faced_arrays)


def real_upper_records(case, mesh, value, rows=None, z_forced="absent"):
    """A jump Upper result with arrays (an interior field, and Z when `rows` is given), its certificate under a gate
    checker name and a method report naming both by hash. The field is schema-only: the evaluator never reads it."""
    problem, problem_arrays = case_problem_cached(case, mesh)
    n_e, n_f = len(problem_arrays["tets"]), len(problem["footings"])
    arrays = _float_arrays(element_velocity=np.zeros((n_e, 10, 3)), footing_velocity=np.zeros((n_f, 6)),
                           interior_element_velocity=np.zeros((n_e, 10, 3)),
                           interior_footing_velocity=np.zeros((n_f, 6)))
    if rows is not None:
        arrays["apex_zero"] = np.array(rows, dtype="<i8")
    result = _result(case["name"], problem, "upper", value, arrays, JUMP_ID, mesh)
    certificate = _certificate(problem, problem_arrays, result)
    certificate["checker"] = {"name": "independent-checker", "version": "1"}
    certificate["checks"]["repair_mu"] = 1e-7
    if rows is not None:
        certificate["checks"]["rho_apex"] = 1e-15  # M3
    method = _method_report(certificate)
    method["interior_field"] = {"source": "result", "sha256": interior_field_sha256("upper", JUMP_ID, arrays)}
    if rows is not None:
        method["apex_zero"] = {"count": len(rows), "sha256": apex_zero_sha256(arrays), "z_forced": z_forced}
    return result, arrays, certificate, method


def swap_mesh_upper(case, records, mesh, real):
    """Replace the synthetic jump Upper of one mesh by a real record set."""
    certificates, results, feasibility, methods = records
    sha = load_problem_hashes()[case["name"]][f"m{mesh}"]
    drop = {doc["result_sha256"] for doc in certificates
            if doc["problem_sha256"] == sha and doc["formulation"] == JUMP_ID}
    result, _, certificate, method = real
    certificates = [doc for doc in certificates if doc["result_sha256"] not in drop] + [certificate]
    results = [doc for doc in results if doc["result_sha256"] not in drop] + [result]
    methods = [doc for doc in methods if doc["result_sha256"] not in drop] + [method]
    for record in feasibility:
        if record["problem_sha256"] == sha and record["witness_result_sha256"] in drop:
            record["witness_result_sha256"] = result["result_sha256"]
    return certificates, results, feasibility, methods


def test_gate_path_checks_the_interior_field_hash():
    """Amendment 2, m1: on the gate path a `result` interior field is checked against the result's arrays."""
    case = ALL_CASES[STRIP]
    real = real_upper_records(case, 0, 5.25)
    records = swap_mesh_upper(case, synthetic_records(case, 5.05, 5.25), 0, real)
    result, arrays = real[0], real[1]
    unseen = evaluate_case(case, *records)
    assert not unseen["passed"] and not unseen["invalid_records"], "an unverified report is not counted, not invalid"
    assert any(item["cause"] == "unproved" and item["mesh"] == 0 for item in unseen["items"])
    seen = evaluate_case(case, *records, result_arrays={result["result_sha256"]: arrays})
    assert seen["passed"], failed(seen)
    tampered = {k: v.copy() for k, v in arrays.items()}
    tampered["interior_element_velocity"][0, 0, 0] = 1e-9
    verdict = evaluate_case(case, *records, result_arrays={result["result_sha256"]: tampered})
    assert {"kind": "result_arrays", "id": result["result_sha256"], "code": "hash_mismatch"} in verdict["invalid_records"]
    assert not verdict["passed"]


def test_thirty_degree_face_needs_the_result_arrays_on_the_gate_path():
    case = ALL_CASES[PHI30]
    real = real_upper_records(case, 0, 4.5, rows=[[0, 0], [0, 1]])
    records = swap_mesh_upper(case, synthetic_records(case, 4.0, 4.5, uppers=[JUMP_ID]), 0, real)
    verdict = evaluate_case(case, *records)
    assert ("certified_upper", "unproved") in failed(verdict)
    assert not verdict["claim_withheld"] and verdict["p0_blocking"], "unproved is never withheld under OD2"
    verdict = evaluate_case(case, *records, result_arrays={real[0]["result_sha256"]: real[1]})
    assert verdict["passed"], failed(verdict)


def test_face_that_makes_the_equalities_unprovable_is_refused_and_keeps_p0_open():
    """Amendment 2: kappa_eq is proved on the augmented system; when no route proves it the checker refuses."""
    case = ALL_CASES[PHI30]
    real = real_upper_records(case, 0, 4.5, rows=[[0, 0], [0, 1]])
    certificates, results, feasibility, methods = swap_mesh_upper(
        case, synthetic_records(case, 4.0, 4.5, uppers=[JUMP_ID]), 0, real)
    result, certificate = real[0], real[2]
    certificates = [doc for doc in certificates if doc is not certificate]
    refusal = synthetic_refusal(case, 0, result, kind="conditioning", stage="dependent_rows")
    validate_refusal(refusal, case_problem_cached(case, 0)[0], result)
    verdict = evaluate_case(case, certificates, results, feasibility, methods, refusals=[refusal],
                            result_arrays={result["result_sha256"]: real[1]})
    assert ("certified_upper", "refused") in failed(verdict)
    assert not verdict["claim_withheld"] and verdict["p0_blocking"]
    program = evaluate_program([verdict])
    assert program["p0_open"] and program["refused_cases"] == [PHI30]


@pytest.mark.parametrize("change,code", [
    ({"budget": None}, "refusal_budget"),
    ({"budget": {"memory_bytes": 1.0, "wall_s": 1.0}}, "missing_key"),
    ({"budget": {"memory_bytes": 0.0, "wall_s": 1.0, "exceeded": "memory"}}, "refusal_budget"),
    ({"budget": {"memory_bytes": 1.0, "wall_s": 1.0, "exceeded": None}}, "refusal_budget"),
    ({"budget": {"memory_bytes": 1.0, "wall_s": 1.0, "exceeded": "patience"}}, "refusal_budget"),
    ({"kind": "conditioning"}, "refusal_budget"),
], ids=["absent", "no_exceeded", "zero_memory", "capacity_without_exceeded", "unknown_exceeded",
        "conditioning_with_exceeded"])
def test_capacity_refusal_declares_its_budget(change, code):
    """Amendment 2, m4."""
    case, _, results, _, _, refusal, gone = refused_records(PHI30, 5)
    problem = case_problem_cached(case, 5)[0]
    result = next(doc for doc in results if doc["result_sha256"] == gone["result_sha256"])
    broken = copy.deepcopy(refusal)
    broken.update(change)
    with pytest.raises(ContractError) as error:
        validate_refusal(broken, problem, result)
    assert error.value.code == code


def test_lower_objective_is_checked_and_upper_objective_is_not():
    """Amendment 2, m5."""
    s = load_set("vm_plate_shear")
    lower, arrays = copy.deepcopy(s["lower"]), s["lower_arrays"]
    lower["objective"] = -lower["value"] * (1.0 + 1e-15)
    lower["result_sha256"] = document_hash(lower, arrays, "result_sha256")
    with pytest.raises(ContractError) as error:
        validate_result(lower, arrays, s["problem"], s["arrays"])
    assert error.value.code == "result_objective"
    upper, arrays = copy.deepcopy(s["upper"]), s["upper_arrays"]
    upper["objective"] = 123.0
    upper["result_sha256"] = document_hash(upper, arrays, "result_sha256")
    validate_result(upper, arrays, s["problem"], s["arrays"])


TRESCA_STRIPS = ("tresca_strip_vm_thin_slab", "tresca_strip_weighted_insensitive")
PRANDTL_SLABS = ("dp_prandtl_cohesion_thin_slab", "dp_prandtl_surcharge_thin_slab")


@pytest.mark.parametrize("name", TRESCA_STRIPS + PRANDTL_SLABS)
def test_only_the_tresca_strips_make_continuous_p2_a_control(name):
    """Amendment 2, item B (owner decision OD3): continuous P2 is a control on the two Tresca strips; the Prandtl
    thin slabs still bind it until the owner decides them case by case."""
    case = ALL_CASES[name]
    if name in TRESCA_STRIPS:
        assert case["formulation"] == [LOWER_ID, JUMP_ID] and case["control"] == [CONTROL_ID]
    else:
        assert case["formulation"] == [LOWER_ID, CONTROL_ID, JUMP_ID] and not case["control"]


def test_strip_control_miss_does_not_fail_the_gate_and_a_jump_miss_does():
    """OD3: on the strip the P2c pair may miss MAX_GAP_STRIP (gap 0.079 as recorded in 19.2.2) without failing the
    gate; a jump-pair miss still fails it."""
    case = ALL_CASES[STRIP]
    lower, jump_upper, control_upper = 5.00525, 5.18571, 5.4168
    assert (control_upper - lower) / (0.5 * (control_upper + lower)) > TOLERANCES["MAX_GAP_STRIP"]
    certificates, results, feasibility, methods = synthetic_records(case, lower, jump_upper)
    for mesh in range(len(case["mesh"])):
        result, certificate, method = synthetic_pair(case, mesh, "upper", CONTROL_ID, control_upper)
        results.append(result)
        certificates.append(certificate)
        methods.append(method)
    verdict = evaluate_case(case, certificates, results, feasibility, methods)
    assert verdict["passed"], failed(verdict)
    assert not any(CONTROL_ID in item["formulation"] for item in verdict["items"]), "the control binds nothing"
    miss = evaluate_case(case, *synthetic_records(case, lower, 5.3))
    assert not miss["passed"] and ("gap_le", "violated") in failed(miss) and miss["p0_blocking"]


def test_symmetry_set_normal_is_the_stated_formula_bitwise_on_a_tilted_plane():
    """Amendment 2, M5: a tilted facet, the formula evaluated independently in scalar IEEE double, identical bits."""
    nodes = np.array([[0.1, 0.2, 0.3], [1.3, 0.1, 0.7], [0.2, 1.1, 0.9], [0.3, 0.4, 1.9]])
    tets = np.array([list(range(4)) + [0] * 6])
    arrays = {"nodes": nodes, "tets": tets, "facets": np.array([[0, 3]]), "facet_set": np.array([0])}
    xa, xb, xc = (nodes[v].tolist() for v in FACET_VERTICES[3])
    u = [xb[i] - xa[i] for i in range(3)]
    w = [xc[i] - xa[i] for i in range(3)]
    q = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
    norm = math.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2])
    expected = [q[i] / norm for i in range(3)]
    normal = symmetry_set_normal({}, arrays, 0)
    assert np.count_nonzero(normal) == 3, "the plane is tilted"
    assert normal.tobytes() == np.array(expected).tobytes()
    np.testing.assert_allclose(normal, facet_geometry(nodes, 3)[0], rtol=0, atol=4 * UNIT_ROUNDOFF)


def test_symmetry_set_normal_formula_is_exact():
    for name in ("tresca_strip_vm_thin_slab", "upper3d_roller_phi30"):
        problem, arrays = case_problem_cached(ALL_CASES[name], 0)
        for index, item in enumerate(problem["boundary_sets"]):
            if item["kind"] != "symmetry":
                continue
            normal = symmetry_set_normal(problem, arrays, index)
            row = int(np.nonzero(arrays["facet_set"] == index)[0][0])
            e, f = (int(v) for v in arrays["facets"][row])
            np.testing.assert_allclose(normal, facet_geometry(arrays["nodes"][arrays["tets"][e][:4]], f)[0],
                                       rtol=0, atol=4 * UNIT_ROUNDOFF)
            assert np.count_nonzero(normal) == 1 and abs(normal).max() == 1.0, "an axis-aligned set gets a unit axis"


if __name__ == "__main__":
    if "--write-fixtures" in sys.argv:
        write_fixtures()
    else:
        sys.exit(pytest.main([__file__, "-q"]))
