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

TOLERANCES = {
    "TOL_GEOM_STRAIGHT": 1e-9,
    "TOL_GEOM_VOLUME": 1e-10,
    "TOL_GEOM_PLANE": 1e-9,
    "TOL_EQ_ACCEPT": 1e-6,
    "TOL_EQ_FINAL": 1e-11,
    "TOL_CONE_ACCEPT": 1e-6,
    "TOL_CONE_MARGIN": 1e-12,
    "TOL_INTERIOR_MIN": 1e-6,
    "TOL_REPAIR_MAX": 1e-3,
    "TOL_FLOW_ACCEPT": 1e-6,
    "TOL_FLOW_FINAL": 1e-11,
    "TOL_WORK_MIN": 1e-12,
    "TOL_FEAS_OK": 1e-9,
    "TOL_FEAS_LOCKED": 1e-6,
    "FLOOR_MIN_REL": 1e-12,
    "TOL_CASE_ONE_SIDED": 1e-6,
    "TOL_CASE_EXACT_GAP": 1e-6,
}

CRITERIA = {"von_mises", "drucker_prager"}
RESERVED_CRITERIA = {"tresca", "mohr_coulomb"}
BOUNDARY_KINDS = {"free", "fixed", "symmetry", "traction", "footing"}
RESERVED_BOUNDARY_KINDS = {"interface"}
INTERFACE_LAWS = {"bonded"}
RESERVED_INTERFACE_LAWS = {"no_tension", "coulomb"}
QUERY_FORMS = {"ray", "fixed_components", "support"}
FORMULATIONS = {"lower": {"lower.p1d-vertex.v0"}, "upper": {"upper.p2c-vertex.v0"}}
PROPOSED_FORMULATIONS = {"upper.p2d-jump-bernstein.v0"}
SOLVER_STATUS_WITH_VALUE = {"solved", "almost_solved"}
SOLVER_STATUS_NO_VALUE = {
    "primal_infeasible", "dual_infeasible", "time_limit", "iteration_limit", "numerical_error", "not_run",
}
DOMAIN_CLAIMS = {
    "box_stated", "unverified", "halfspace_by_zero_extension", "halfspace_by_extension_elements",
}
REASON_CODES = {
    "no_candidate", "hash_mismatch", "query_mismatch", "formulation_not_admitted", "curved_element",
    "degenerate_element", "boundary_violation", "equality_residual_too_large",
    "equality_correction_failed", "cone_violation_too_large", "no_interior_field",
    "interior_field_invalid", "repair_too_large", "flow_violation_too_large",
    "no_admissible_correction", "work_nonpositive", "value_mismatch", "unsupported_domain_claim",
}
CHECK_KEYS = {
    "rho_div", "rho_jump", "rho_bc", "rho_resultant", "rho_eq_final", "rho_cone", "repair_mu",
    "straightness_max", "rho_flow", "rho_iso_final", "work_reference",
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
    "prandtl_tresca": lambda p: (2.0 + math.pi) * p["cu"],
    "prandtl_nc": lambda p: p["c"] * prandtl_nc(p["phi_deg"]),
    "prandtl_nq": lambda p: p["q0"] * prandtl_nq(p["phi_deg"]),
}


def shear_band_limit(params, first_layer):
    return params["k"] / (params["gamma"] * (params["H"] - 0.5 * first_layer))


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
        "name": raw["name"], "source": [], "checked": [], "mesh": [], "formulation": [], "regions": [],
        "footings": [], "boundaries": [], "queries": [], "F0": {}, "F1": {}, "restrained": [],
        "assert": [], "nested": [],
    }
    for key, rest in raw["lines"]:
        tokens = rest.split()
        if key in ("status", "class", "desc", "units", "on_failure"):
            case[key] = rest
        elif key in ("source", "checked", "formulation"):
            case[key].append(rest)
        elif key == "box":
            case["box"] = _vector(tokens, 6)
        elif key == "mesh":
            case["mesh"].append(("kuhn", tuple(int(token) for token in tokens)))
        elif key == "mesh_file":
            case["mesh"].append(("file", rest))
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
        elif key == "assert":
            case["assert"].append(tokens)
        else:
            fail("case_format", f"unknown key {key}")
    return case


def build_case_problem(case, mesh_index=0):
    """Build the problem document and arrays of one case on one of its meshes."""
    kind, spec = case["mesh"][mesh_index]
    lo, hi = np.array(case["box"][:3]), np.array(case["box"][3:])
    if kind == "kuhn":
        vertices, tet4 = kuhn_tet4(lo, hi, spec)
    else:
        with np.load(DATA / spec, allow_pickle=False) as npz:
            vertices, tet4 = npz["nodes"], npz["tets"]
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
    mesh_label = "x".join(map(str, spec)) if kind == "kuhn" else spec.split(".")[0]
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


def validate_formulation(bound, formulation):
    require(formulation not in PROPOSED_FORMULATIONS, "formulation_not_admitted", repr(formulation))
    require(formulation in FORMULATIONS.get(bound, ()), "formulation_unknown", repr(formulation))


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
    validate_formulation(bound, doc["formulation"])
    _keys(doc["producer"], {"name", "version"})
    _keys(doc["solver"], {"name", "version", "settings"})
    status = doc["solver_status"]
    require(status in SOLVER_STATUS_WITH_VALUE | SOLVER_STATUS_NO_VALUE, "solver_status", repr(status))

    n_e, n_n = len(problem_arrays["tets"]), len(problem_arrays["nodes"])
    n_f = len(problem["footings"])
    shapes = {
        "lower": ({"sigma": (n_e, 4, 6), "resultant": (n_f, 6)},
                  {"interior_sigma": (n_e, 4, 6), "equilibrium_multiplier": (n_e, 3)}),
        "upper": ({"velocity": (n_n, 3), "footing_velocity": (n_f, 6)},
                  {"interior_velocity": (n_n, 3), "interior_footing_velocity": (n_f, 6),
                   "strain_rate": (n_e, 4, 6), "dissipation": (n_e, 4)}),
    }
    required, allowed = shapes[bound]
    other = "upper" if bound == "lower" else "lower"
    foreign = (set(shapes[other][0]) | set(shapes[other][1])) & set(arrays)
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
    require(("interior_velocity" in arrays) == ("interior_footing_velocity" in arrays), "result_interior")
    require(doc["result_sha256"] == document_hash(doc, arrays, "result_sha256"), "hash_mismatch", "result")


def expected_program_shape(problem, arrays, query, bound):
    n_e, n_n = len(arrays["tets"]), len(arrays["nodes"])
    n_f = len(problem["footings"])
    free_components = 6 * n_f - len(query["restrained"])
    is_ray = query["form"] != "support"
    kinds = [item["kind"] for item in problem["boundary_sets"]]
    tets, facets, facet_set = arrays["tets"], arrays["facets"], arrays["facet_set"]
    if bound == "lower":
        interior, _ = mesh_facets(tets)
        rows = 3 * n_e + 9 * len(interior)
        for index in facet_set:
            rows += {"free": 9, "traction": 9, "symmetry": 6}.get(kinds[int(index)], 0)
        rows += free_components if is_ray else 0
        return (1 if is_ray else 0) + 24 * n_e, rows + 24 * n_e, rows, 4 * n_e
    node_sets = [set() for _ in kinds]
    for (e, f), index in zip(facets, facet_set):
        node_sets[int(index)].update(int(a) for a in tets[int(e)][list(FACET_NODES[int(f)])])
    fixed = set().union(*[s for s, kind in zip(node_sets, kinds) if kind == "fixed"], set())
    rows = 24 * n_e + 3 * len(fixed)
    footing_nodes = {}
    for s, item in zip(node_sets, problem["boundary_sets"]):
        if item["kind"] == "symmetry":
            rows += len(s - fixed)
        elif item["kind"] == "footing":
            footing_nodes.setdefault(item["footing"], set()).update(s)
    rows += 3 * sum(len(s) for s in footing_nodes.values()) + len(query["restrained"])
    rows += 1 if is_ray else free_components
    return 3 * n_n + 6 * n_f + 24 * n_e, rows + 24 * n_e, rows, 4 * n_e


def validate_program(doc, arrays, problem, problem_arrays):
    require(doc.get("format") == FORMAT_PROGRAM, "format", repr(doc.get("format")))
    require(doc.get("format_version") == FORMAT_VERSION, "format_version")
    _keys(doc, {
        "format", "format_version", "arrays", "problem_sha256", "query_id", "bound", "formulation",
        "objective_sign", "n_variables", "n_rows", "n_zero", "n_soc", "soc_dim",
    })
    require(doc["problem_sha256"] == problem["problem_sha256"], "hash_mismatch")
    query = query_by_id(problem, doc["query_id"])
    bound = doc["bound"]
    require(bound in ("lower", "upper"), "bound", repr(bound))
    validate_formulation(bound, doc["formulation"])
    require(doc["objective_sign"] == (-1 if bound == "lower" else 1), "program_objective_sign")
    require(doc["soc_dim"] == 6, "program_cone")
    n_var, n_row, n_zero, n_soc = expected_program_shape(problem, problem_arrays, query, bound)
    require((doc["n_variables"], doc["n_rows"], doc["n_soc"]) == (n_var, n_row, n_soc), "program_shape")
    require(doc["n_zero"] == n_row - 6 * n_soc, "program_shape", "zero rows")
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
        "numerator", "denominator", "method", "witness_checked", "producer", "tolerances",
    })
    require(doc["problem_sha256"] == problem["problem_sha256"], "hash_mismatch")
    query_by_id(problem, doc["query_id"])
    validate_formulation("upper", doc["formulation"])
    require(doc["status"] in FEASIBILITY_STATUS, "feasibility_status", repr(doc["status"]))
    require(doc["tolerances"] == TOLERANCES, "tolerances")
    for key in ("rho_feas", "numerator", "denominator"):
        _finite_number(doc[key], "feasibility_number")
    require(doc["numerator"] >= 0.0 and doc["denominator"] > 0.0, "feasibility_number")
    rho = doc["rho_feas"]
    require(abs(rho - doc["numerator"] / doc["denominator"]) <= 1e-12 * max(1.0, abs(rho)), "feasibility_ratio")
    if doc["status"] == "formulation_feasible":
        require(rho <= TOLERANCES["TOL_FEAS_OK"] and doc["witness_checked"] is True, "feasibility_verdict")
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
    validate_formulation(doc["bound"], doc["formulation"])
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

def _result(stem, problem, bound, value, arrays):
    formulation = next(iter(FORMULATIONS[bound]))
    doc = {
        "format": FORMAT_RESULT, "format_version": FORMAT_VERSION,
        "arrays": f"{stem}.q.{bound}.result.npz", "problem_sha256": problem["problem_sha256"],
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


def build_fixture_set():
    cases = parse_cases()
    documents = {}
    for name in ("vm_uniaxial_tension", "vm_plate_shear"):
        problem, arrays = build_case_problem(cases[name], 0)
        value, sigma, resultant, velocity, footing_velocity = _exact_fields(name, arrays)
        lower_arrays = {"sigma": sigma.astype("<f8"), "resultant": resultant.astype("<f8")}
        upper_arrays = {"velocity": velocity.astype("<f8"), "footing_velocity": footing_velocity.astype("<f8")}
        lower = _result(name, problem, "lower", value, lower_arrays)
        upper = _result(name, problem, "upper", value, upper_arrays)
        lower_cert = _certificate(problem, arrays, lower)
        upper_cert = _certificate(problem, arrays, upper)
        documents[f"{name}.problem"] = (problem, arrays)
        documents[f"{name}.q.lower.result"] = (lower, lower_arrays)
        documents[f"{name}.q.upper.result"] = (upper, upper_arrays)
        documents[f"{name}.q.lower.certificate"] = (lower_cert, None)
        documents[f"{name}.q.upper.certificate"] = (upper_cert, None)
        documents[f"{name}.q.bracket"] = (_bracket(lower_cert, upper_cert), None)
        documents[f"{name}.q.feasibility"] = ({
            "format": FORMAT_FEASIBILITY, "format_version": FORMAT_VERSION,
            "problem_sha256": problem["problem_sha256"], "query_id": "q",
            "formulation": "upper.p2c-vertex.v0", "status": "formulation_feasible", "rho_feas": 0.0,
            "numerator": 0.0, "denominator": 1.0, "method": "analytic witness; denominator is schema-only",
            "witness_checked": True, "producer": {"name": SCHEMA_FIXTURE_CHECKER, "version": "0"},
            "tolerances": dict(TOLERANCES),
        }, None)
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


# ---------------------------------------------------------------------------------------------
# Tests: case file
# ---------------------------------------------------------------------------------------------

ALL_CASES = parse_cases()
ACTIVE = sorted(name for name, case in ALL_CASES.items() if case["status"] == "active")
REFERENCE = [name for name in ACTIVE if ALL_CASES[name]["class"] != "upper_feasibility_3d"]
FEASIBILITY = [name for name in ACTIVE if ALL_CASES[name]["class"] == "upper_feasibility_3d"]
TODO = sorted(name for name, case in ALL_CASES.items() if case["status"] == "todo")


def test_case_inventory_covers_the_required_set():
    assert "tresca_strip_vm_thin_slab" in REFERENCE
    assert {"dp_prandtl_cohesion_thin_slab", "dp_prandtl_surcharge_thin_slab"} <= set(REFERENCE)
    assert {ALL_CASES[name]["phi_deg"] for name in FEASIBILITY} >= {0.0, 10.0, 20.0, 30.0}
    assert "strip_self_weight_ngamma" in TODO


@pytest.mark.parametrize("name", REFERENCE)
def test_reference_case_values_and_assertions(name):
    case = ALL_CASES[name]
    assert case["source"] and case["checked"] and case["units"] and case["formulation"]
    expected = TARGET_FORMULAS[case["formula"]](case["formula_params"])
    assert case["target"] == pytest.approx(expected, rel=1e-14)
    kinds = {tokens[0] for tokens in case["assert"]}
    assert {"lower_le_target", "upper_ge_target"} <= kinds, "both one-sided assertions are required"
    assert "closeness" not in kinds
    for tokens in case["assert"]:
        if tokens[0] in ("lower_le_target", "upper_ge_target", "lower_ge_target"):
            assert float(tokens[1]) == TOLERANCES["TOL_CASE_ONE_SIDED"]
    if case["class"] == "exact_on_mesh":
        assert ["gap_le", "1e-6"] in case["assert"]
        assert float("1e-6") == TOLERANCES["TOL_CASE_EXACT_GAP"]
    if case["class"] == "exact_lower_on_mesh":
        assert "lower_ge_target" in kinds and "upper_le_formula" in kinds
    if case["class"] == "one_sided_reference":
        assert "gap_le" not in kinds and "lower_ge_target" not in kinds


@pytest.mark.parametrize("name", [n for n in REFERENCE if ALL_CASES[n]["class"] == "one_sided_reference"])
def test_strip_cases_satisfy_the_box_conditions(name):
    case = ALL_CASES[name]
    phi = case.get("phi_deg", 0.0)
    x_max, depth = prandtl_extent(1.0, phi)
    x0, y0, _, x1, y1, _ = case["box"]
    assert x0 == 0.0 and y1 == 0.0
    assert x1 >= x_max and y1 - y0 >= depth, "the mechanism must fit inside the box"
    kinds = {item["name"]: item["kind"] for item in case["boundaries"]}
    assert kinds["far"] == "fixed" and kinds["base"] == "fixed", "the reference needs fixed walls"
    assert kinds["z0"] == kinds["z1"] == kinds["axis"] == "symmetry"
    assert ("strip", "ux") in case["restrained"] and ("strip", "rz") in case["restrained"]
    region = case["regions"][0]
    if region["criterion"] == "drucker_prager":
        c = 0.0 if region["k"] == 0.0 else 1.0
        k, alpha = plane_strain_match(c, phi)
        assert region["alpha"] == pytest.approx(alpha, rel=1e-14)
        assert region["k"] == pytest.approx(k, rel=1e-14, abs=1e-14)


@pytest.mark.parametrize("name", FEASIBILITY)
def test_feasibility_cases_are_true_3d_and_unstructured(name):
    case = ALL_CASES[name]
    assert case["target"] is None and case["on_failure"] == "formulation_result"
    assert case["formulation"] == ["upper.p2c-vertex.v0"]
    assert ["feasibility_status", "formulation_feasible"] in case["assert"]
    assert len(case["mesh"]) >= 3 and all(kind == "file" for kind, _ in case["mesh"])
    region = case["regions"][0]
    if case["phi_deg"] == 0.0:
        assert region["criterion"] == "von_mises"
    else:
        k, alpha = plane_strain_match(1.0, case["phi_deg"])
        assert (region["k"], region["alpha"]) == (pytest.approx(k, rel=1e-14), pytest.approx(alpha, rel=1e-14))
    lo, hi = np.array(case["box"][:3]), np.array(case["box"][3:])
    for _, file_name in case["mesh"]:
        with np.load(DATA / file_name, allow_pickle=False) as npz:
            vertices, tets = npz["nodes"], npz["tets"]
        assert vertices.dtype.str == "<f8" and tets.dtype.str == "<i8" and tets.shape[1] == 4
        assert np.all(vertices >= lo) and np.all(vertices <= hi)
        assert np.ptp(vertices[:, 0]) == np.ptp(vertices[:, 1]) == 4.0 and np.ptp(vertices[:, 2]) == 2.0
        # Not a thin slab, and not a lattice: interior vertices are off any regular grid.
        interior = np.all((vertices > lo + 1e-9) & (vertices < hi - 1e-9), axis=1)
        assert interior.sum() >= 5
        assert len(np.unique(np.round(vertices[interior, 0], 9))) > 0.5 * interior.sum()


@pytest.mark.parametrize("name", TODO)
def test_todo_cases_carry_no_value(name):
    case = ALL_CASES[name]
    assert "target" not in case and not case["assert"] and "TODO" in case["desc"]


@pytest.mark.parametrize("name", ACTIVE)
def test_case_builds_a_valid_problem(name):
    case = ALL_CASES[name]
    meshes = range(len(case["mesh"])) if case["class"] == "upper_feasibility_3d" else [0]
    for index in meshes:
        doc, arrays = build_case_problem(case, index)
        validate_problem(doc, arrays)
    for coarse, fine in case["nested"]:
        assert case["mesh"][fine][1] == tuple(2 * n for n in case["mesh"][coarse][1])


def test_shear_band_limit_is_above_the_target():
    case = ALL_CASES["vm_body_force_shear_block"]
    params = case["formula_params"]
    for _, (_, ny, _) in case["mesh"]:
        assert shear_band_limit(params, params["H"] / ny) > case["target"]


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


@pytest.mark.parametrize("formulation,code", [
    ("upper.p2d-jump-bernstein.v0", "formulation_not_admitted"),
    ("upper.p2c-vertex.v9", "formulation_unknown"),
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
    ({"formulation": "upper.p2d-jump-bernstein.v0"}, "formulation_not_admitted"),
    ({"denominator": 0.0}, "feasibility_number"),
])
def test_feasibility_record_cannot_hide_a_failure(changes, code):
    problem, doc = feasibility_doc(**changes)
    with pytest.raises(ContractError) as error:
        validate_feasibility(doc, problem)
    assert error.value.code == code


def schema_program(s, bound):
    query = s["problem"]["queries"][0]
    n_var, n_row, n_zero_soc, n_soc = expected_program_shape(s["problem"], s["arrays"], query, bound)
    assert n_zero_soc == n_row - 6 * n_soc
    doc = {
        "format": FORMAT_PROGRAM, "format_version": FORMAT_VERSION, "arrays": "unused.npz",
        "problem_sha256": s["problem"]["problem_sha256"], "query_id": "q", "bound": bound,
        "formulation": next(iter(FORMULATIONS[bound])), "objective_sign": -1 if bound == "lower" else 1,
        "n_variables": n_var, "n_rows": n_row, "n_zero": n_row - 6 * n_soc, "n_soc": n_soc, "soc_dim": 6,
    }
    arrays = {
        "A_indptr": np.zeros(n_var + 1, dtype="<i8"), "A_indices": np.zeros(0, dtype="<i8"),
        "A_data": np.zeros(0, dtype="<f8"), "b": np.zeros(n_row, dtype="<f8"), "c": np.zeros(n_var, dtype="<f8"),
    }
    return doc, arrays


@pytest.mark.parametrize("bound", ["lower", "upper"])
def test_program_schema_and_canonical_shape(bound):
    s = load_set("vm_plate_shear")
    doc, arrays = schema_program(s, bound)
    validate_program(doc, arrays, s["problem"], s["arrays"])
    n_e = len(s["arrays"]["tets"])
    if bound == "lower":
        assert doc["n_variables"] == 1 + 24 * n_e
    else:
        assert doc["n_variables"] == 3 * len(s["arrays"]["nodes"]) + 6 + 24 * n_e
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
    found = dict(re.findall(r"^\| `((?:TOL|FLOOR)_[A-Z_]+)` \| ([0-9.e+-]+) \|", interface_text(), re.M))
    assert {name: float(value) for name, value in found.items()} == TOLERANCES


@pytest.mark.parametrize("terms", [
    sorted(REASON_CODES), sorted(DOMAIN_CLAIMS), sorted(CRITERIA | RESERVED_CRITERIA),
    sorted(BOUNDARY_KINDS | RESERVED_BOUNDARY_KINDS), sorted(INTERFACE_LAWS | RESERVED_INTERFACE_LAWS),
    sorted(QUERY_FORMS), sorted(SOLVER_STATUS_WITH_VALUE | SOLVER_STATUS_NO_VALUE),
    sorted(FEASIBILITY_STATUS), sorted(CHECK_KEYS),
    sorted(FORMULATIONS["lower"] | FORMULATIONS["upper"] | PROPOSED_FORMULATIONS),
    [FORMAT_PROBLEM, FORMAT_RESULT, FORMAT_PROGRAM, FORMAT_FEASIBILITY, FORMAT_CERTIFICATE, FORMAT_BRACKET],
])
def test_interface_defines_every_term(terms):
    text = interface_text()
    missing = [term for term in terms if f"`{term}`" not in text]
    assert not missing


def test_interface_conventions_object_is_literal():
    block = re.search(r"### 15\.7.*?```json\n(.*?)```", interface_text(), re.S).group(1)
    assert json.loads(block) == CONVENTIONS


def test_interface_marks_the_design_changes_as_proposed():
    text = interface_text()
    assert "pending the owner's decision" in text
    for item in ("P-1", "P-2", "P-3", "P-4", "P-5", "P-6", "P-7"):
        assert f"| {item} |" in text
    assert "`upper.p2c-vertex.v0` | Upper | admitted **as the baseline under test**" in text


def test_public_files_carry_no_private_material():
    forbidden = ["TI" + "Ms", "Op" + "tum", "Fa" + "ble", "C:" + "\\Users", "nmo" + "ra"]
    paths = [INTERFACE, CASES, Path(__file__)]
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
