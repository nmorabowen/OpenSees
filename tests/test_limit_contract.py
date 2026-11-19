"""WP-237a: the limit-analysis contract, executable part.

Schema validators for the version-0 JSON + NPZ formats, the parser and mesh builder for
``tests/limit_cases.txt``, convention checks and fixtures. The mathematics is decided in
``tests/limit_interface.md``. Nothing here assembles a conic program or certifies a bound.

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

FORMAT_VERSION = 0
FORMAT_PROBLEM = "ladruno-limit-problem"
FORMAT_RESULT = "ladruno-limit-result"
FORMAT_PROGRAM = "ladruno-limit-program"
FORMAT_FEASIBILITY = "ladruno-limit-feasibility"
FORMAT_CERTIFICATE = "ladruno-limit-certificate"
FORMAT_BRACKET = "ladruno-limit-bracket"

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
    "unsupported_domain_claim",
}
CHECK_KEYS = {
    "rho_div", "rho_jump", "rho_bc", "rho_resultant", "projection_move", "omega_eq", "kappa_eq",
    "tau_eq", "rho_cone", "repair_mu", "straightness_max", "rho_flow", "rho_jump_cone",
    "work_reference",
}
FEASIBILITY_STATUS = {"formulation_feasible", "formulation_infeasible", "inconclusive"}
SCHEMA_FIXTURE_CHECKER = "schema-fixture"

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
        return b"d" + struct.pack("<d", float(obj))
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
        digest.update(b"a" + canonical_encode(name) + canonical_encode(array.dtype.str))
        digest.update(_u64(array.ndim) + b"".join(_u64(dim) for dim in array.shape))
        digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def strict_json_loads(text):
    def refuse(token):
        fail("json_non_finite", token)
    return json.loads(text, parse_constant=refuse)


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


def equality_floor(kappa):
    """The conditioning-aware final equality threshold (interface 17.2)."""
    return TOLERANCES["TOL_EQ_ROUNDING_FACTOR"] * UNIT_ROUNDOFF * max(1.0, kappa)


def jump_constants(alpha):
    """tau (dilation per unit slip) and k_J of interface section 12."""
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
    cases, current = {}, None
    for raw in Path(path).read_text(encoding="ascii").splitlines():
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
        "assert": [], "gate": [], "nested": [],
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
        else:
            fail("case_format", f"unknown key {key}")
    return case


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
        "arrays": f"{case['name']}.problem.npz", "name": f"{case['name']}:{mesh_label}",
        "units": {"length": units["length"], "force": units["force"]},
        "conventions": copy.deepcopy(CONVENTIONS),
        "domain": {"target": "box", "truncation_sets": [], "exterior_load_free": True},
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
    _keys(domain, {"target", "truncation_sets", "exterior_load_free"})
    require(domain["target"] in ("box", "halfspace"), "domain_target", repr(domain["target"]))
    require(isinstance(domain["exterior_load_free"], bool), "domain_target")
    require(set(domain["truncation_sets"]) <= set(names), "domain_truncation", "unknown set")
    if domain["target"] == "box":
        require(not domain["truncation_sets"], "domain_truncation", "box has no truncation sets")
    else:
        require(domain["truncation_sets"], "domain_truncation", "halfspace needs truncation sets")

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
    other = "upper" if bound == "lower" else "lower"
    foreign = (set(shapes[other][0]) | set(shapes[other][1]) | wrong_velocity) & set(arrays)
    require(not foreign, "result_wrong_bound_arrays", str(sorted(foreign)))
    unknown = set(arrays) - set(required) - set(allowed)
    require(not unknown, "array_names", str(sorted(unknown)))
    if status in SOLVER_STATUS_WITH_VALUE:
        _finite_number(doc["value"], "result_value")
        _finite_number(doc["objective"], "result_value")
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
    require(doc["result_sha256"] == document_hash(doc, arrays, "result_sha256"), "hash_mismatch", "result")


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
        "witness_checked", "producer", "tolerances",
    })
    require(doc["problem_sha256"] == problem["problem_sha256"], "hash_mismatch")
    query_by_id(problem, doc["query_id"])
    validate_formulation("upper", doc["formulation"], problem)
    require(doc["status"] in FEASIBILITY_STATUS, "feasibility_status", repr(doc["status"]))
    require(doc["tolerances"] == TOLERANCES, "tolerances")
    for key in ("rho_feas", "numerator", "denominator"):
        _finite_number(doc[key], "feasibility_number")
    require(doc["numerator"] >= 0.0 and doc["denominator"] > 0.0, "feasibility_number")
    rho = doc["rho_feas"]
    require(abs(rho - doc["numerator"] / doc["denominator"]) <= 1e-12 * max(1.0, abs(rho)), "feasibility_ratio")
    if doc["rank"] is not None:
        require(doc["nullity"] == doc["n_unknowns"] - doc["rank"], "feasibility_rank")
    if doc["status"] == "formulation_feasible":
        require(rho <= TOLERANCES["TOL_FEAS_OK"] and doc["witness_checked"] is True, "feasibility_verdict")
        require(doc["nullity"] is None or doc["nullity"] >= 1, "feasibility_rank")
        require(doc["footing_modes"] is None or doc["footing_modes"] >= 1, "feasibility_rank")
    elif doc["status"] == "formulation_infeasible":
        require(rho >= TOLERANCES["TOL_FEAS_LOCKED"], "feasibility_verdict")
    else:
        require(rho > TOLERANCES["TOL_FEAS_OK"] or doc["witness_checked"] is not True, "feasibility_verdict")
        require(rho < TOLERANCES["TOL_FEAS_LOCKED"], "feasibility_verdict")


def supported_domain_claim(problem, arrays, bound):
    """The strongest claim version 0 can support for a certified bound (interface 14)."""
    domain = problem["domain"]
    if domain["target"] == "box":
        return "box_stated"
    kinds = {item["name"]: item["kind"] for item in problem["boundary_sets"]}
    walls = [kinds[name] for name in domain["truncation_sets"]]
    if bound == "upper":
        return "halfspace_by_zero_extension" if all(kind == "fixed" for kind in walls) else "unverified"
    weightless = all(not any(region["b0"]) and not any(region["b1"]) for region in problem["regions"])
    if all(kind == "free" for kind in walls) and weightless and domain["exterior_load_free"]:
        return "halfspace_by_zero_extension"
    return "unverified"


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
    else:
        require(doc["reasons"], "certificate_reason", "estimate needs a reason")
        require(doc["value_certified"] is None and doc["floor"] is None, "certificate_value")


def gate_admissible(certificate):
    """Only a real checker's certified bound can satisfy a gate."""
    return certificate["classification"] == "certified_bound" and certificate["checker"]["name"] != SCHEMA_FIXTURE_CHECKER


def bracket_claim(lower_claim, upper_claim):
    if lower_claim == upper_claim == "halfspace_by_zero_extension":
        return lower_claim
    return "box_stated"


def validate_bracket(doc, lower, upper):
    require(doc.get("format") == FORMAT_BRACKET, "format", repr(doc.get("format")))
    require(doc.get("format_version") == FORMAT_VERSION, "format_version")
    _keys(doc, {
        "format", "format_version", "problem_sha256", "query_id", "lower_result_sha256",
        "upper_result_sha256", "lower", "upper", "floor_lower", "floor_upper", "status", "gap",
        "domain_claim", "lower_domain_claim", "upper_domain_claim",
    })
    require(lower["bound"] == "lower" and upper["bound"] == "upper", "bracket_labels")
    require(lower["classification"] == upper["classification"] == "certified_bound", "bracket_uncertified")
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
# Fixtures: exact fields on the smallest meshes. No solver is involved.
# ---------------------------------------------------------------------------------------------

KAPPA_PLATE_LOWER = 51.88  # sigma_max / sigma_min of the row-equilibrated Lower equalities, 6-tet plate


def facet_work_table(area):
    table = np.zeros((3, 6))
    for i in range(3):
        for j in range(3):
            table[i, j] = area / 30.0 if i == j else -area / 60.0
        for q, (a, b) in enumerate(((0, 1), (1, 2), (2, 0))):
            table[i, 3 + q] = 2.0 * area / 15.0 if i in (a, b) else area / 15.0
    return table


def fixture_upper_value(problem, arrays, velocity, footing_velocity, quadrature="vertex", body="exact"):
    """Value of a continuous P2 fixture field. Fixture tooling only; it is not the checker."""
    nodes, tets = arrays["nodes"], arrays["tets"]
    query = problem["queries"][0]
    dissipation, w0, w1 = 0.0, 0.0, 0.0
    for e, tet in enumerate(tets):
        region = problem["regions"][int(arrays["element_region"][e])]
        k, alpha = region["k"], region.get("alpha", 0.0)
        _, volume = tet_gradients(nodes[tet[:4]])
        rates = [strain_rate_at_vertex(nodes[tet], velocity[tet], m) for m in range(4)]
        if quadrature == "vertex":
            dissipation += sum(volume / 4.0 * support_function(d, k, alpha, tol=1e-9) for d in rates)
        else:
            dissipation += volume * support_function(np.mean(rates, axis=0), k, alpha, tol=1e-9)
        weights = volume * (np.array([-0.05] * 4 + [0.2] * 6) if body == "exact" else np.array([0.25] * 4 + [0.0] * 6))
        moment = weights @ velocity[tet]
        w0 += np.dot(region["b0"], moment)
        w1 += np.dot(region["b1"], moment)
    for row, (e, f) in enumerate(arrays["facets"]):
        _, area = facet_geometry(nodes[tets[e][:4]], int(f))
        table = facet_work_table(area)
        facet_velocity = velocity[tets[e][list(FACET_NODES[int(f)])]]
        w0 += float(np.sum(arrays["t0"][row] * (table @ facet_velocity)))
        w1 += float(np.sum(arrays["t1"][row] * (table @ facet_velocity)))
    w0 += float(np.sum(np.array(query["F0"]).reshape(-1, 6) * footing_velocity))
    w1 += float(np.sum(np.array(query["F1"]).reshape(-1, 6) * footing_velocity))
    return (dissipation - w0) / w1


def _result(stem, problem, bound, value, arrays, formulation=None, label=None):
    formulation = formulation or (LOWER_ID if bound == "lower" else CONTROL_ID)
    doc = {
        "format": FORMAT_RESULT, "format_version": FORMAT_VERSION,
        "arrays": f"{stem}.q.{label or bound}.result.npz", "problem_sha256": problem["problem_sha256"],
        "query_id": "q", "bound": bound, "formulation": formulation,
        "producer": {"name": "analytic-fixture", "version": "0"},
        "solver": {"name": "none", "version": "0", "settings": {}},
        "solver_status": "solved", "objective": -value if bound == "lower" else value,
        "value": value, "wall_time_s": None,
    }
    doc["result_sha256"] = document_hash(doc, arrays, "result_sha256")
    return doc


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
        "checks": {key: None for key in sorted(CHECK_KEYS)},
        "checker": {"name": SCHEMA_FIXTURE_CHECKER, "version": "0"},
        "tolerances": dict(TOLERANCES),
    }


def _bracket(lower, upper):
    low, high = lower["value_certified"], upper["value_certified"]
    return {
        "format": FORMAT_BRACKET, "format_version": FORMAT_VERSION,
        "problem_sha256": lower["problem_sha256"], "query_id": lower["query_id"],
        "lower_result_sha256": lower["result_sha256"], "upper_result_sha256": upper["result_sha256"],
        "lower": low, "upper": high, "floor_lower": lower["floor"], "floor_upper": upper["floor"],
        "status": "bracket", "gap": (high - low) / (0.5 * (high + low)),
        "domain_claim": bracket_claim(lower["domain_claim"], upper["domain_claim"]),
        "lower_domain_claim": lower["domain_claim"], "upper_domain_claim": upper["domain_claim"],
    }


def _feasibility(problem):
    return {
        "format": FORMAT_FEASIBILITY, "format_version": FORMAT_VERSION,
        "problem_sha256": problem["problem_sha256"], "query_id": "q", "formulation": CONTROL_ID,
        "status": "formulation_feasible", "rho_feas": 0.0, "numerator": 0.0, "denominator": 1.0,
        "n_unknowns": 10, "n_rows": 4, "rank": 4, "nullity": 6, "footing_modes": None,
        "method": "analytic witness; counts and denominator are schema-only", "witness_checked": True,
        "producer": {"name": SCHEMA_FIXTURE_CHECKER, "version": "0"}, "tolerances": dict(TOLERANCES),
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


def _add_exact_set(documents, name, problem, arrays):
    value, sigma, resultant, velocity, footing_velocity = _exact_fields(name, arrays)
    lower_arrays = _float_arrays(sigma=sigma, resultant=resultant)
    upper_arrays = _float_arrays(velocity=velocity, footing_velocity=footing_velocity)
    jump_arrays = _float_arrays(element_velocity=velocity[arrays["tets"]], footing_velocity=footing_velocity)
    lower = _result(name, problem, "lower", value, lower_arrays)
    upper = _result(name, problem, "upper", value, upper_arrays)
    jump = _result(name, problem, "upper", value, jump_arrays, JUMP_ID, "upper-jump")
    lower_cert = _certificate(problem, arrays, lower)
    upper_cert = _certificate(problem, arrays, upper)
    documents[f"{name}.problem"] = (problem, arrays)
    documents[f"{name}.q.lower.result"] = (lower, lower_arrays)
    documents[f"{name}.q.upper.result"] = (upper, upper_arrays)
    documents[f"{name}.q.upper-jump.result"] = (jump, jump_arrays)
    documents[f"{name}.q.lower.certificate"] = (lower_cert, None)
    documents[f"{name}.q.upper.certificate"] = (upper_cert, None)
    documents[f"{name}.q.bracket"] = (_bracket(lower_cert, upper_cert), None)
    documents[f"{name}.q.feasibility"] = (_feasibility(problem), None)


def build_negative_fixtures(documents):
    plate = documents["vm_plate_shear.problem"][1]
    nodes, tets = plate["nodes"], plate["tets"]
    centre = int(np.argmin(np.linalg.norm(nodes - 0.5, axis=1)))
    footprint = int(np.nonzero((np.abs(nodes[:, 1] - 1.0) < 1e-12) & (np.abs(nodes[:, 0]) < 1e-12))[0][0])
    wall = int(np.nonzero((np.abs(nodes[:, 1]) < 1e-12) & (np.abs(nodes[:, 0] - 1.0) < 1e-12))[0][0])
    grads, _ = tet_gradients(nodes[tets[0][:4]])
    vertex = int(np.argmax(np.abs(grads[:, 0])))
    row_scale = float(np.sum(np.abs(grads[:, 1])))
    amount = 3.0 * equality_floor(KAPPA_PLATE_LOWER) * row_scale / abs(grads[vertex, 0])
    dp_tets = documents["dp_uniaxial_tension.problem"][1]["tets"]
    interior, _ = mesh_facets(dp_tets)
    (e1, f1), (e2, _) = interior[0]
    shared = int(dp_tets[e1][FACET_VERTICES[f1][0]])
    local2 = [int(v) for v in dp_tets[e2]].index(shared)
    P, L, U = "vm_plate_shear", "vm_plate_shear.q.lower.result", "vm_plate_shear.q.upper.result"
    G = "vm_dead_gravity_shear_layer"
    H = "vm_plate_shear_halfspace"
    entries = [
        ("perturbed_stress_component", "one perturbed stress component", P, L, "checker",
         {"kind": "perturb_sigma", "element": 0, "vertex": 0, "component": 0, "amount": 1e-3},
         ["projection_move_too_large"]),
        ("stress_scaled_1p01", "stress scaled by 1.01", P, L, "checker",
         {"kind": "scale_lower", "factor": 1.01}, ["repair_too_large", "cone_violation_too_large"]),
        ("facet_traction_jump", "a facet traction jump", P, L, "checker",
         {"kind": "element_offset", "element": 0, "component": 0, "amount": 1e-2}, ["projection_move_too_large"]),
        ("wrong_load_normalisation", "a wrong load normalisation", P, U, "checker",
         {"kind": "scale_value", "factor": 1.1}, ["value_mismatch"]),
        ("swapped_bound_labels", "swapped LB and UB labels", P, L, "schema",
         {"kind": "relabel_upper"}, "result_wrong_bound_arrays"),
        ("curved_tetrahedron", "a curved tetrahedron", P, None, "schema",
         {"kind": "move_mid_node", "element": 0, "local_node": 4, "amount": 0.01}, "element_curved"),
        ("mismatched_problem_hash", "a mismatched problem hash", P, U, "schema",
         {"kind": "problem_hash"}, "hash_mismatch"),
        ("unsupported_domain_claim", "an unsupported domain_claim", P, U, "schema",
         {"kind": "certificate_claim", "claim": "halfspace_by_zero_extension"}, "unsupported_domain_claim"),
        ("isochory_one_vertex", "tr(D) != 0 at one vertex of an isochoric Upper field", P, U, "checker",
         {"kind": "perturb_velocity", "node": centre, "component": 0, "amount": 1e-3}, ["projection_move_too_large"]),
        ("footprint_off_rigid", "a footprint node off the rigid motion", P, U, "checker",
         {"kind": "perturb_velocity", "node": footprint, "component": 1, "amount": 1e-3}, ["projection_move_too_large"]),
        ("wall_velocity_zero_extension", "nonzero wall velocity under a zero-extension claim", H,
         f"{H}.q.upper.result", "checker",
         {"kind": "perturb_velocity", "node": wall, "component": 0, "amount": 1e-3}, ["projection_move_too_large"]),
        ("jump_between_nodes", "a jump admissible at facet nodes but not across the facet", "dp_uniaxial_tension",
         "dp_uniaxial_tension.q.upper-jump.result", "checker",
         {"kind": "node_jump", "element": e2, "local_node": local2, "facet": [e1, f1], "amount": 1e-3},
         ["flow_violation_too_large", "no_admissible_correction"]),
        ("wrong_gravity_work", "wrong gravity work or wrong dissipation quadrature", G, f"{G}.q.upper.result", "checker",
         {"kind": "value_rule", "body": "lumped", "quadrature": "vertex"}, ["value_mismatch"]),
        ("wrong_dissipation_quadrature", "wrong gravity work or wrong dissipation quadrature", G, f"{G}.q.upper.result",
         "checker", {"kind": "value_rule", "body": "exact", "quadrature": "centroid"}, ["value_mismatch"]),
        ("equality_just_above_floor", "an equality residual just above the floor", P, L, "final_check",
         {"kind": "perturb_sigma", "element": 0, "vertex": vertex, "component": 0, "amount": amount,
          "kappa_reference": KAPPA_PLATE_LOWER, "floor_multiple": 3.0}, ["equality_residual_above_floor"]),
    ]
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

    # The same plate as a half-space target, with the fixed base as the truncation boundary.
    problem, arrays = build_case_problem(cases["vm_plate_shear"], 0)
    name = "vm_plate_shear_halfspace"
    problem["arrays"], problem["name"] = f"{name}.problem.npz", name
    problem["domain"] = {"target": "halfspace", "truncation_sets": ["base"], "exterior_load_free": True}
    _with_hash(problem, arrays)
    value, _, _, velocity, footing_velocity = _exact_fields("vm_plate_shear", arrays)
    upper_arrays = _float_arrays(velocity=velocity, footing_velocity=footing_velocity)
    documents[f"{name}.problem"] = (problem, arrays)
    documents[f"{name}.q.upper.result"] = (_result(name, problem, "upper", value, upper_arrays), upper_arrays)

    # Drucker-Prager tension with the broken velocity space: the exact field has no jumps.
    name = "dp_uniaxial_tension"
    problem, arrays = build_case_problem(cases[name], 0)
    alpha, k = problem["regions"][0]["alpha"], problem["regions"][0]["k"]
    a, c = 1.0 / S3 + alpha, -0.5 / S3 + alpha
    x = arrays["nodes"]
    velocity = np.column_stack([x[:, 0], c / a * x[:, 1], c / a * x[:, 2]])
    value = k / a
    footing = np.zeros((0, 6))
    jump_arrays = _float_arrays(element_velocity=velocity[arrays["tets"]], footing_velocity=footing)
    documents[f"{name}.problem"] = (problem, arrays)
    documents[f"{name}.q.upper-jump.result"] = (
        _result(name, problem, "upper", value, jump_arrays, JUMP_ID, "upper-jump"), jump_arrays)

    # Dead gravity: the exact Lower field, and an admissible quadratic Upper field whose work and
    # dissipation are sensitive to the quadrature rule.
    name = "vm_dead_gravity_shear_layer"
    problem, arrays = build_case_problem(cases[name], 0)
    gamma, height, target = 0.3, 1.0, cases[name]["target"]
    sigma = np.zeros((len(arrays["tets"]), 4, 6))
    sigma[..., 5] = target + gamma * (height - arrays["nodes"][arrays["tets"][:, :4], 2])
    z = arrays["nodes"][:, 2]
    velocity = np.column_stack([z * z - 0.5 * z, np.zeros_like(z), np.zeros_like(z)])
    footing = np.zeros((0, 6))
    value = fixture_upper_value(problem, arrays, velocity, footing)
    lower_arrays = _float_arrays(sigma=sigma, resultant=footing)
    upper_arrays = _float_arrays(velocity=velocity, footing_velocity=footing)
    documents[f"{name}.problem"] = (problem, arrays)
    documents[f"{name}.q.lower.result"] = (_result(name, problem, "lower", target, lower_arrays), lower_arrays)
    documents[f"{name}.q.upper.result"] = (_result(name, problem, "upper", value, upper_arrays), upper_arrays)

    documents["negative_fixtures"] = (build_negative_fixtures(documents), None)
    return documents


def write_fixtures():
    DATA.mkdir(parents=True, exist_ok=True)
    for stem, (doc, arrays) in build_fixture_set().items():
        (DATA / f"{stem}.json").write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        if arrays is not None:
            np.savez(DATA / doc["arrays"], **arrays)


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
        k, alpha = plane_strain_match(1.0, phi)
        tau, root = jump_constants(alpha)
        assert tau == pytest.approx(math.tan(math.radians(phi)), abs=1e-14)
        if phi > 0.0:
            assert k / root * (1.0 / tau) * tau == pytest.approx(k / root)
            assert k / (3.0 * alpha) * tau == pytest.approx(1.0)  # slip dissipation is c = 1
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


def test_equality_floor_grows_with_conditioning():
    base = TOLERANCES["TOL_EQ_ROUNDING_FACTOR"] * UNIT_ROUNDOFF
    assert equality_floor(0.5) == base
    assert equality_floor(1e4) == pytest.approx(1e4 * base)
    assert 1e4 * equality_floor(1e4) < TOLERANCES["TOL_FLOOR_MAX"]


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


@pytest.mark.parametrize("name", GATE_3D)
def test_true_3d_gate_cases_are_complete(name):
    case = ALL_CASES[name]
    assert case["target"] is None
    assert JUMP_ID in case["formulation"] and case["control"] == [CONTROL_ID]
    gates = [tokens[0] for tokens in case["gate"]]
    for required in ("kinematic_status", "certified_bracket", "gap_le", "upper_nonincrease", "lower_nondecrease"):
        assert required in gates, "feasibility alone must not pass"
    limit = next(tokens[1] for tokens in case["gate"] if tokens[0] == "gap_le")
    assert limit == ("MAX_GAP_3D_VM" if case["phi_deg"] == 0.0 else "MAX_GAP_3D_DP")
    files = [entry for entry in case["mesh"] if entry[2] == 0]
    assert len(files) == 4 and all(entry[0] == "file" for entry in case["mesh"])
    assert case["nested"] == [(0, 4), (3, 5)]
    region = case["regions"][0]
    if case["phi_deg"] == 0.0:
        assert region["criterion"] == "von_mises"
    else:
        k, alpha = plane_strain_match(1.0, case["phi_deg"])
        assert (region["k"], region["alpha"]) == (pytest.approx(k, rel=1e-14), pytest.approx(alpha, rel=1e-14))
    wall = "fixed" if "fixed" in name else "symmetry"
    kinds = {item["name"]: item["kind"] for item in case["boundaries"]}
    assert all(kinds[key] == wall for key in ("xmin", "xmax", "ymin", "ymax", "base"))


@pytest.mark.parametrize("file_name,count", [("upper3d_a", 243), ("upper3d_b", 463), ("upper3d_c", 949), ("upper3d_d", 2840)])
def test_gate_meshes_are_unstructured_with_a_conformal_footprint(file_name, count):
    with np.load(DATA / f"{file_name}.mesh.npz", allow_pickle=False) as npz:
        vertices, tets = npz["nodes"], npz["tets"]
    assert vertices.dtype.str == "<f8" and tets.dtype.str == "<i8" and tets.shape == (count, 4)
    assert all(np.linalg.det(vertices[t[1:]] - vertices[t[0]]) > 0.0 for t in tets)
    lo, hi = np.array([-2.0, -2.0, -2.0]), np.array([2.0, 2.0, 0.0])
    assert np.all(vertices >= lo) and np.all(vertices <= hi)
    interior = np.all((vertices > lo + 1e-9) & (vertices < hi - 1e-9), axis=1)
    assert interior.sum() >= 5
    assert len(np.unique(np.round(vertices[interior, 0], 9))) > 0.5 * interior.sum(), "not a lattice"
    valence = np.bincount(tets.ravel())
    assert len(np.unique(valence[interior])) >= 3, "interior vertex valence varies"
    case = ALL_CASES["upper3d_fixed_phi00"]
    index = [i for i, entry in enumerate(case["mesh"]) if entry[1].startswith(file_name) and entry[2] == 0][0]
    doc, arrays = build_case_problem(case, index)
    validate_problem(doc, arrays)
    names = [item["name"] for item in doc["boundary_sets"]]
    rows = np.nonzero(arrays["facet_set"] == names.index("contact"))[0]
    area = sum(facet_geometry(arrays["nodes"][arrays["tets"][e][:4]], int(f))[1] for e, f in arrays["facets"][rows])
    assert area == pytest.approx(1.0, abs=1e-12), "the footprint is exactly the conforming unit square"


def test_gate_refined_mesh_is_a_valid_problem():
    case = ALL_CASES["upper3d_roller_phi20"]
    assert case["mesh"][4] == ("file", "upper3d_a.mesh.npz", 1)
    assert case["mesh"][5] == ("file", "upper3d_d.mesh.npz", 1)
    doc, arrays = build_case_problem(case, 4)
    validate_problem(doc, arrays)
    assert len(arrays["tets"]) == 8 * 243


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
    problem, arrays = load_document(DATA / f"{name}.problem.json")
    lower, lower_arrays = load_document(DATA / f"{name}.q.lower.result.json")
    upper, upper_arrays = load_document(DATA / f"{name}.q.upper.result.json")
    lower_cert, _ = load_document(DATA / f"{name}.q.lower.certificate.json")
    upper_cert, _ = load_document(DATA / f"{name}.q.upper.certificate.json")
    bracket, _ = load_document(DATA / f"{name}.q.bracket.json")
    feasibility, _ = load_document(DATA / f"{name}.q.feasibility.json")
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
    jump, jump_arrays = load_document(DATA / f"{name}.q.upper-jump.result.json")
    validate_result(jump, jump_arrays, s["problem"], s["arrays"])
    assert jump["formulation"] == JUMP_ID and "velocity" not in jump_arrays
    validate_certificate(s["lower_cert"], s["problem"], s["arrays"], s["lower"])
    validate_certificate(s["upper_cert"], s["problem"], s["arrays"], s["upper"])
    validate_bracket(s["bracket"], s["lower_cert"], s["upper_cert"])
    validate_feasibility(s["feasibility"], s["problem"])
    assert not gate_admissible(s["lower_cert"]) and not gate_admissible(s["upper_cert"])


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
    jump, jump_arrays = load_document(DATA / "vm_plate_shear.q.upper-jump.result.json")
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


def test_gate_refuses_schema_fixture_and_accepts_a_real_checker_name():
    _, doc = certificate_case()
    assert not gate_admissible(doc)
    doc["checker"]["name"] = "independent-checker"
    assert gate_admissible(doc)


@pytest.mark.parametrize("target,walls,kind,weight,free,bound,expected", [
    ("box", [], "fixed", False, True, "upper", "box_stated"),
    ("box", [], "free", False, True, "lower", "box_stated"),
    ("halfspace", ["wall"], "fixed", False, True, "upper", "halfspace_by_zero_extension"),
    ("halfspace", ["wall"], "symmetry", False, True, "upper", "unverified"),
    ("halfspace", ["wall"], "fixed", False, True, "lower", "unverified"),
    ("halfspace", ["wall"], "free", False, True, "lower", "halfspace_by_zero_extension"),
    ("halfspace", ["wall"], "free", True, True, "lower", "unverified"),
    ("halfspace", ["wall"], "free", False, False, "lower", "unverified"),
    ("halfspace", ["wall"], "free", True, True, "upper", "unverified"),
])
def test_domain_claim_is_per_bound(target, walls, kind, weight, free, bound, expected):
    problem = {
        "domain": {"target": target, "truncation_sets": walls, "exterior_load_free": free},
        "boundary_sets": [{"name": "wall", "kind": kind}],
        "regions": [{"b0": [0.0, -1.0 if weight else 0.0, 0.0], "b1": [0.0, 0.0, 0.0]}],
    }
    assert supported_domain_claim(problem, {}, bound) == expected


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


def feasibility_doc(**changes):
    s = load_set("vm_plate_shear")
    doc = copy.deepcopy(s["feasibility"])
    doc.update(changes)
    return s["problem"], doc


def test_feasibility_failure_is_a_recordable_result():
    problem, doc = feasibility_doc(status="formulation_infeasible", numerator=0.1, denominator=1.0,
                                   rho_feas=0.1, witness_checked=False)
    validate_feasibility(doc, problem)


@pytest.mark.parametrize("changes,code", [
    ({"status": "formulation_feasible", "numerator": 0.1, "denominator": 1.0, "rho_feas": 0.1}, "feasibility_verdict"),
    ({"status": "formulation_feasible", "witness_checked": False}, "feasibility_verdict"),
    ({"status": "formulation_infeasible"}, "feasibility_verdict"),
    ({"status": "inconclusive", "numerator": 1.0, "denominator": 1.0, "rho_feas": 1.0}, "feasibility_verdict"),
    ({"status": "solver_said_ok"}, "feasibility_status"),
    ({"rho_feas": 0.5}, "feasibility_ratio"),
    ({"formulation": "upper.p2d-jump-bernstein.v1"}, "formulation_not_admitted"),
    ({"denominator": 0.0}, "feasibility_number"),
    ({"nullity": 0, "rank": 10}, "feasibility_rank"),
    ({"nullity": 3}, "feasibility_rank"),
    ({"footing_modes": 0}, "feasibility_rank"),
])
def test_feasibility_record_cannot_hide_a_failure(changes, code):
    problem, doc = feasibility_doc(**changes)
    with pytest.raises(ContractError) as error:
        validate_feasibility(doc, problem)
    assert error.value.code == code


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
    [FORMAT_PROBLEM, FORMAT_RESULT, FORMAT_PROGRAM, FORMAT_FEASIBILITY, FORMAT_CERTIFICATE, FORMAT_BRACKET],
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
    for item in ("R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8"):
        assert f"| {item} |" in text[rulings:]
    assert "not frozen" in text[:600]


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
        certificate = load_document(DATA / f"{entry['problem']}.q.upper.certificate.json")[0]
        certificate["domain_claim"] = m["claim"]
    elif kind == "perturb_velocity":
        arrays["velocity"][m["node"], m["component"]] += m["amount"]
    elif kind == "node_jump":
        e1, f1 = m["facet"]
        normal, _ = facet_geometry(problem_arrays["nodes"][problem_arrays["tets"][e1][:4]], f1)
        arrays["element_velocity"][m["element"], m["local_node"]] += m["amount"] * normal
    elif kind == "value_rule":
        value = fixture_upper_value(problem, problem_arrays, arrays["velocity"], arrays["footing_velocity"],
                                    quadrature=m["quadrature"], body=m["body"])
        result["value"] = result["objective"] = value
    else:
        fail("negative_fixture", kind)
    if result is not None and kind != "problem_hash":
        result["result_sha256"] = document_hash(result, arrays, "result_sha256")
    return problem, problem_arrays, result, arrays, certificate


def test_negative_fixtures_cover_every_harness_item():
    assert {entry["harness_item"] for entry in NEGATIVE} == HARNESS_ITEMS
    assert len({entry["id"] for entry in NEGATIVE}) == len(NEGATIVE)
    for entry in NEGATIVE:
        if entry["stage"] == "schema":
            assert entry["expect"] in ("result_wrong_bound_arrays", "element_curved", "hash_mismatch",
                                       "unsupported_domain_claim")
        else:
            assert set(entry["expect"]) <= REASON_CODES and entry["expect"]


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


@pytest.mark.parametrize("entry", [e for e in NEGATIVE if e["stage"] != "schema"], ids=lambda e: e["id"])
def test_checker_negative_fixture_is_well_formed_and_live(entry):
    problem, problem_arrays, result, arrays, _ = apply_negative(entry)
    validate_result(result, arrays, problem, problem_arrays)
    base, base_arrays = load_document(DATA / f"{entry['result']}.json")
    m = entry["mutation"]
    nodes, tets = problem_arrays["nodes"], problem_arrays["tets"]
    kind = m["kind"]
    if kind == "perturb_sigma" and entry["stage"] == "checker":
        assert m["amount"] > TOLERANCES["TOL_PROJECTION_MOVE"] * np.max(np.abs(base_arrays["sigma"]))
    elif kind == "perturb_sigma":
        grads, _ = tet_gradients(nodes[tets[m["element"]][:4]])
        residual = abs(grads[m["vertex"], 0] * m["amount"])
        scale = residual + float(np.sum(np.abs(grads[:, 1]) * np.abs(base_arrays["sigma"][m["element"], :, 3])))
        floor = equality_floor(m["kappa_reference"])
        assert 2.9 * floor <= residual / scale <= 3.1 * floor
        assert (residual / scale) * m["kappa_reference"] < TOLERANCES["TOL_FLOOR_MAX"]
    elif kind == "scale_lower":
        k = problem["regions"][0]["k"]
        assert max(yield_from_principal(s, k) for s in arrays["sigma"].reshape(-1, 6)) > TOLERANCES["TOL_CONE_ACCEPT"]
    elif kind == "element_offset":
        interior, _ = mesh_facets(tets)
        worst = 0.0
        for (e1, f1), (e2, _) in interior:
            if m["element"] in (e1, e2):
                normal, _ = facet_geometry(nodes[tets[e1][:4]], f1)
                vertex1 = FACET_VERTICES[f1][0]
                vertex2 = [int(v) for v in tets[e2]].index(int(tets[e1][vertex1]))
                jump = voigt_n(normal) @ (arrays["sigma"][e1, vertex1] - arrays["sigma"][e2, vertex2])
                worst = max(worst, float(np.max(np.abs(jump))))
        assert worst > TOLERANCES["TOL_PROJECTION_MOVE"]
    elif kind in ("scale_value", "value_rule"):
        true_value = fixture_upper_value(problem, problem_arrays, arrays["velocity"], arrays["footing_velocity"])
        assert abs(result["value"] - true_value) > TOLERANCES["TOL_VALUE_MATCH"] * abs(true_value)
        assert true_value == pytest.approx(base["value"], rel=1e-12)
    elif kind == "perturb_velocity":
        node, changed = m["node"], arrays["velocity"] - base_arrays["velocity"]
        assert np.count_nonzero(changed) == 1 and m["amount"] > TOLERANCES["TOL_PROJECTION_MOVE"]
        x = nodes[node]
        if entry["id"] == "isochory_one_vertex":
            theta = []
            for tet in tets:
                if node in tet:
                    theta += [strain_split(strain_rate_at_vertex(nodes[tet], arrays["velocity"][tet], v))[0]
                              for v in range(4)]
            assert max(abs(t) for t in theta) > TOLERANCES["TOL_FLOW_ACCEPT"]
        elif entry["id"] == "footprint_off_rigid":
            assert x[1] == 1.0
            rigid = base_arrays["footing_velocity"][0, :3]
            assert np.linalg.norm(arrays["velocity"][node] - rigid) == pytest.approx(m["amount"])
        else:
            assert x[1] == 0.0 and problem["domain"]["target"] == "halfspace"
            assert supported_domain_claim(problem, problem_arrays, "upper") == "halfspace_by_zero_extension"
            assert np.linalg.norm(arrays["velocity"][node]) == pytest.approx(m["amount"])
    elif kind == "node_jump":
        e1, f1 = m["facet"]
        e2 = m["element"]
        tau, _ = jump_constants(problem["regions"][0]["alpha"])
        normal, _ = facet_geometry(nodes[tets[e1][:4]], f1)
        tangents = tangent_basis(normal)
        order2 = [int(v) for v in tets[e2]]
        minus = arrays["element_velocity"][e1][list(FACET_NODES[f1])]
        plus = arrays["element_velocity"][e2][[order2.index(int(tets[e1][a])) for a in FACET_NODES[f1]]]
        nodal = plus - minus

        def admissible(c):
            return c @ normal >= tau * np.hypot(c @ tangents[0], c @ tangents[1]) - 1e-12 * m["amount"]

        assert all(admissible(c) for c in nodal), "the jump passes at all six facet nodes"
        assert not all(admissible(c) for c in bernstein_from_nodal(nodal)), "and fails across the facet"
    else:
        fail("negative_fixture", kind)


def test_correct_gravity_and_dissipation_rules_differ_from_the_wrong_ones():
    problem, arrays = load_document(DATA / "vm_dead_gravity_shear_layer.problem.json")
    upper, fields = load_document(DATA / "vm_dead_gravity_shear_layer.q.upper.result.json")
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


def test_tracked_probe_output_is_labelled_and_shows_the_contrast():
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


def test_public_files_carry_no_private_material():
    forbidden = ["TI" + "Ms", "Op" + "tum", "Fa" + "ble", "C:" + "\\Users", "nmo" + "ra"]
    paths = [INTERFACE, CASES, Path(__file__)] + [p for p in (PROBE_SCRIPT, PROBE_OUTPUT) if p.exists()]
    paths += [path for path in DATA.iterdir() if path.suffix == ".json"]
    for path in paths:
        text = path.read_text(encoding="utf-8")
        for word in forbidden:
            assert word not in text, f"{word!r} in {path.name}"


if __name__ == "__main__":
    if "--write-fixtures" in sys.argv:
        write_fixtures()
    else:
        sys.exit(pytest.main([__file__, "-q"]))
