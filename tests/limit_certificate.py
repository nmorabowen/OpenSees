"""WP-237c: the independent certificate checker for limit-analysis bounds.

Written from the frozen contract alone: ``tests/limit_interface.md`` (sections 2 to 6, 11, 12, 14, 15 and 17),
``tests/limit_cases.txt`` and the documents under ``tests/data/limit/``. It reads a problem document and a
result document, rebuilds every operator from the mesh, and runs, in this order: affine projection, inequality
repair, a recheck of the final field from scratch (interface 17.3 and 17.4). It never reads solver cones,
solver tolerances, stored strain rates, stored dissipation or any claim made by the producer.

Routes chosen to differ from a conic solver's:

* gradients from the facet route ``grad L_i = -(A_i / (3 V)) n_i`` (interface 3.3);
* yield from principal stresses, and the flow cone from principal strain rates;
* the jump cone from opening and slip taken directly, with no tangent basis and no second-order cone;
* quadrature weights from the factorial formula for barycentric monomials.

Entry points: ``check`` (a certificate), ``final_check`` (the recheck alone: a verdict, never a certificate),
``bracket``, ``feasibility`` (the record of interface 15.5), ``conditioning_report`` and the ``conditioning`` command
(kappa_eq before any result exists), and the ``write_*`` and ``read_*`` helpers.

The frozen certificate has no key for the method behind ``kappa_eq``. ``check`` therefore returns a method report
with every certificate, ``write_certificate`` stores it beside the certificate, and ``gate_admissible`` and
``bracket`` take no certificate without it.

``kappa_eq`` comes from a dense singular value decomposition up to ``DENSE_ENTRIES_MAX`` entries. Beyond that it is
an upper bound proved by Sylvester counts on a factorised normal matrix, with every dependent row proved by a dense
decomposition of its row block. A system the checker cannot evaluate raises ``CapacityRefusal``: that is neither an
estimate nor a withheld claim, and it says nothing about the result.

Uses the standard library, numpy and scipy. ``tests/test_limit_certificate.py`` cross-checks every table held
here against the contract test.
"""
import hashlib
import json
import math
import struct
from pathlib import Path

import numpy as np
import scipy.linalg
import scipy.sparse as sp
import scipy.sparse.linalg as spla

CHECKER_NAME = "ladruno-limit-checker"
CHECKER_VERSION = "0.1.0"
SCHEMA_FIXTURE_CHECKER = "schema-fixture"
PROCEDURE_DEMO_CHECKER = "procedure-demo"

FORMAT_VERSION = 0
FORMAT_PROBLEM = "ladruno-limit-problem"
FORMAT_RESULT = "ladruno-limit-result"
FORMAT_FEASIBILITY = "ladruno-limit-feasibility"
FORMAT_CERTIFICATE = "ladruno-limit-certificate"
FORMAT_BRACKET = "ladruno-limit-bracket"
FORMAT_METHOD = "ladruno-limit-certificate-method"      # the checker's own record beside a certificate; not a frozen format
FORMAT_CONDITIONING = "ladruno-limit-conditioning-report"
INTERIOR_FIELDS = ("none", "zero", "supplied")

LOWER_ID = "lower.p1d-vertex.v0"
CONTROL_ID = "upper.p2c-vertex.v0"
JUMP_ID = "upper.p2d-jump-bernstein.v0"
ADMITTED = {"lower": (LOWER_ID,), "upper": (CONTROL_ID, JUMP_ID)}

# Interface 17.8, transcribed. The test file pins this table to the interface text and to the contract test.
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
UNIT_ROUNDOFF = 2.0 ** -53
TAU_EQ = TOLERANCES["TOL_EQ_ROUNDING_FACTOR"] * UNIT_ROUNDOFF

CONVENTIONS = {
    "sign": "tension_positive", "voigt": ["xx", "yy", "zz", "xy", "yz", "zx"], "stress_shear": "tensor",
    "strain_rate_shear": "engineering", "tet10": "vtk", "facets": "opposite_vertex_outward",
    "footing_resultant": "footing_on_soil",
}
CHECK_KEYS = (
    "rho_div", "rho_jump", "rho_bc", "rho_resultant", "projection_move", "omega_eq", "kappa_eq", "tau_eq",
    "rho_cone", "repair_mu", "straightness_max", "rho_flow", "rho_jump_cone", "work_reference",
)
REASON_CODES = (
    "no_candidate", "hash_mismatch", "query_mismatch", "formulation_not_admitted", "formulation_not_applicable",
    "curved_element", "degenerate_element", "boundary_violation", "projection_move_too_large",
    "equality_residual_above_floor", "equality_ill_conditioned", "cone_violation_too_large", "no_interior_field",
    "interior_field_invalid", "repair_too_large", "flow_violation_too_large", "no_admissible_correction",
    "work_nonpositive", "value_mismatch", "unsupported_domain_claim", "literal_zero_violated",
)

# Interface 3.1 and 3.2.
EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
FACET_VERTICES = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))
FACET_EDGES = ((0, 1), (1, 2), (2, 0))   # facet edge q joins facet vertices (a, b)
AXES = {"x": 0, "y": 1, "z": 2}
# Traction component c reads stress component v times normal component j (interface 2.1).
TRACTION_MAP = (((0, 0), (3, 1), (5, 2)), ((1, 1), (3, 0), (4, 2)), ((2, 2), (4, 1), (5, 0)))

# Capacity limits of this implementation. They bound cost; none is an acceptance threshold.
DENSE_ENTRIES_MAX = 3.6e7      # a dense singular value decomposition is used up to this many matrix entries
SPARSE_ROWS_MAX = 250_000      # largest normal matrix handed to the sparse direct factorisation
SHIFT_REL = 1e-12              # positive shift of the shift-invert step, relative to sigma_max^2 (interface 17.1)
PROBE_REL = 1e-6               # rows above this share of a null-space probe are examined for dependence
PROBES = 2                     # random probes of the null space of the normal matrix
NULL_BASIS_ENTRIES_MAX = 1.0e8 # largest explicit null-space basis, with its image, held in memory
FEAS_DENSE_UNKNOWNS = 4000     # rank, nullity and footing modes are reported up to this many unknowns
LSMR_COMPATIBLE = 1e-9         # relative residual below which a least-squares system counts as solved exactly


class Refusal(ValueError):
    """The documents are malformed, or a record is not fit for the use asked of it: nothing is written."""

    def __init__(self, code, detail=""):
        super().__init__(f"{code}: {detail}" if detail else code)
        self.code = code


class CapacityRefusal(Refusal):
    """The checker cannot evaluate this system: too large for its factorisations, or a conditioning it cannot prove.

    It is never an estimate and never a withheld claim. It says nothing about the result.
    """

    def __init__(self, detail=""):
        super().__init__("checker_capacity", detail)


def _require(condition, code, detail=""):
    if not condition:
        raise Refusal(code, detail)


def _facet_nodes():
    """All six P2 nodes of each local facet: three vertices, then edges (a, b), (b, c), (c, a)."""
    edge_node = {frozenset(edge): 4 + q for q, edge in enumerate(EDGES)}
    return tuple(tuple(v) + tuple(edge_node[frozenset((v[a], v[b]))] for a, b in FACET_EDGES)
                 for v in FACET_VERTICES)


FACET_NODES = _facet_nodes()


# ---------------------------------------------------------------------------------------------
# Documents: canonical hash (interface 15.1) and files
# ---------------------------------------------------------------------------------------------

def _encode(obj):
    if obj is None:
        return b"n"
    if obj is True:
        return b"t"
    if obj is False:
        return b"f"
    if isinstance(obj, (int, float)):
        return b"d" + struct.pack("<d", float(obj) + 0.0)
    if isinstance(obj, str):
        raw = obj.encode("utf-8")
        return b"s" + struct.pack("<Q", len(raw)) + raw
    if isinstance(obj, (list, tuple)):
        return b"l" + struct.pack("<Q", len(obj)) + b"".join(_encode(item) for item in obj)
    if isinstance(obj, dict):
        keys = sorted(obj, key=lambda key: key.encode("utf-8"))
        return b"m" + struct.pack("<Q", len(keys)) + b"".join(_encode(key) + _encode(obj[key]) for key in keys)
    raise Refusal("json_type", repr(type(obj)))


def document_hash(doc, arrays, hash_key):
    digest = hashlib.sha256()
    digest.update(_encode({key: value for key, value in doc.items() if key != hash_key}))
    for name in sorted(arrays):
        array = np.asarray(arrays[name])
        if array.dtype.kind == "f":
            array = array + 0.0
        digest.update(b"a" + _encode(name) + _encode(array.dtype.str))
        digest.update(struct.pack("<Q", array.ndim) + b"".join(struct.pack("<Q", dim) for dim in array.shape))
        digest.update(np.ascontiguousarray(array).tobytes(order="C"))
    return digest.hexdigest()


def load_document(path):
    """A JSON document and the arrays of its NPZ. Non-finite numbers are refused (interface 15)."""
    def refuse(token):
        raise Refusal("json_non_finite", token)

    def finite(token):
        value = float(token)
        _require(math.isfinite(value), "json_non_finite", token)
        return value

    def finite_int(token):
        _require(math.isfinite(float(token)), "json_non_finite", token)
        return int(token)

    path = Path(path)
    doc = json.loads(path.read_text(encoding="utf-8"), parse_constant=refuse, parse_float=finite, parse_int=finite_int)
    arrays = {}
    if doc.get("arrays"):
        with np.load(path.parent / doc["arrays"], allow_pickle=False) as npz:
            arrays = {name: npz[name] for name in npz.files}
    return doc, arrays


def write_document(path, doc):
    Path(path).write_text(json.dumps(doc, indent=1, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def record_stem(result):
    """``<case>.m<m>.<query>.<formulation-id>`` from the result's array file name (interface 15)."""
    name = result["arrays"]
    _require(isinstance(name, str) and name.endswith(".result.npz"), "result_name", repr(name))
    return name[: -len(".result.npz")]


def _number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


# ---------------------------------------------------------------------------------------------
# Mesh, regions and boundary, rebuilt from the problem arrays
# ---------------------------------------------------------------------------------------------

def monomial_integral(powers, measure):
    """Integral of a barycentric monomial over a simplex of the given measure (length, area or volume)."""
    dimension = len(powers) - 1
    top = math.factorial(dimension)
    for power in powers:
        top *= math.factorial(power)
    return measure * top / math.factorial(sum(powers) + dimension)


def volume_weights():
    """int N_a dV / V for the ten P2 shape functions, N_i = 2 L_i^2 - L_i and N_ab = 4 L_a L_b."""
    vertex = 2.0 * monomial_integral((2, 0, 0, 0), 1.0) - monomial_integral((1, 0, 0, 0), 1.0)
    edge = 4.0 * monomial_integral((1, 1, 0, 0), 1.0)
    return np.array([vertex] * 4 + [edge] * 6)


def facet_work_table():
    """int L_i N_j dA / A: linear traction basis i against the six P2 facet functions j."""
    table = np.zeros((3, 6))
    for i in range(3):
        for j in range(3):
            square = [0, 0, 0]
            square[j] += 2
            square[i] += 1
            linear = [0, 0, 0]
            linear[j] += 1
            linear[i] += 1
            table[i, j] = 2.0 * monomial_integral(square, 1.0) - monomial_integral(linear, 1.0)
        for q, (a, b) in enumerate(FACET_EDGES):
            powers = [0, 0, 0]
            powers[a] += 1
            powers[b] += 1
            powers[i] += 1
            table[i, 3 + q] = 4.0 * monomial_integral(powers, 1.0)
    return table


def traction_moment_weights():
    """int L_a L_b dA / A (interface 6): 1/6 on the diagonal, 1/12 off it."""
    table = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            powers = [0, 0, 0]
            powers[a] += 1
            powers[b] += 1
            table[a, b] = monomial_integral(powers, 1.0)
    return table


def bernstein_map():
    """Bernstein coefficients of a quadratic from its six nodal values (interface 12.1)."""
    out = np.eye(6)
    for q, (a, b) in enumerate(FACET_EDGES):
        out[3 + q, 3 + q] = 2.0
        out[3 + q, a] = out[3 + q, b] = -0.5
    return out


def p2_vertex_gradient_table():
    """T[m, a, i]: grad N_a at vertex m is sum_i T[m, a, i] grad L_i (interface 3.3)."""
    table = np.zeros((4, 10, 4))
    for m in range(4):
        for i in range(4):
            table[m, i, i] = 4.0 * (1.0 if i == m else 0.0) - 1.0
        for q, (a, b) in enumerate(EDGES):
            table[m, 4 + q, b] += 4.0 * (1.0 if a == m else 0.0)
            table[m, 4 + q, a] += 4.0 * (1.0 if b == m else 0.0)
    return table


def traction_matrix(normal):
    """N(n) of interface 2.1 for normals of shape (..., 3): shape (..., 3, 6)."""
    normal = np.asarray(normal, dtype=float)
    out = np.zeros(normal.shape[:-1] + (3, 6))
    for c, pairs in enumerate(TRACTION_MAP):
        for v, j in pairs:
            out[..., c, v] = normal[..., j]
    return out


def tangent_pair(normal):
    """T1, T2 of interface 5: m is the first index of the smallest |n_m|."""
    normal = np.asarray(normal, dtype=float)
    m = int(np.argmin(np.abs(normal)))
    axis = np.zeros(3)
    axis[m] = 1.0
    t1 = np.cross(axis, normal)
    t1 = t1 / np.linalg.norm(t1)
    return t1, np.cross(normal, t1)


def skew(vector):
    x, y, z = vector
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


class Model:
    """Geometry, strength, loads and boundary of one problem, from its arrays alone."""

    def __init__(self, problem, arrays):
        _require(problem.get("format") == FORMAT_PROBLEM and problem.get("format_version") == FORMAT_VERSION, "format")
        _require(problem.get("conventions") == CONVENTIONS, "conventions")
        for name, kind in (("nodes", "f"), ("tets", "i"), ("element_region", "i"), ("facets", "i"), ("facet_set", "i"),
                           ("t0", "f"), ("t1", "f")):
            _require(name in arrays, "missing_array", name)
            _require(arrays[name].dtype.str == ("<f8" if kind == "f" else "<i8"), "array_dtype", name)
        self.problem = problem
        self.nodes = np.asarray(arrays["nodes"], dtype=float)
        self.tets = np.asarray(arrays["tets"], dtype=np.int64)
        _require(self.nodes.ndim == 2 and self.nodes.shape[1] == 3 and np.all(np.isfinite(self.nodes)), "array_shape", "nodes")
        _require(self.tets.ndim == 2 and self.tets.shape[1] == 10 and len(self.tets) > 0, "array_shape", "tets")
        _require(self.tets.min() >= 0 and self.tets.max() < len(self.nodes), "tet_index")
        self.n_elements, self.n_nodes = len(self.tets), len(self.nodes)
        self.ell = problem["length_scale"]
        _require(_number(self.ell) and self.ell > 0.0, "length_scale")
        self._geometry()
        self._regions(arrays)
        self._facets(arrays)
        self._boundary(arrays)
        self._domain()

    def _domain(self):
        """The domain declaration of interface 14. An unknown or wrongly typed set is refused by name."""
        domain, sets = self.problem.get("domain"), self.problem["boundary_sets"]
        _require(isinstance(domain, dict) and domain.get("target") in ("box", "halfspace"), "domain_target")
        _require(isinstance(domain.get("exterior_load_free"), bool), "domain_target", "exterior_load_free")
        names = [item.get("name") for item in sets]
        _require(len(set(names)) == len(names), "boundary_set_name", "duplicate")
        walls, mirrors = domain.get("truncation_sets"), domain.get("symmetry_sets")
        _require(isinstance(walls, list) and all(name in names for name in walls), "domain_truncation", "unknown truncation set")
        _require(isinstance(mirrors, list) and all(name in names for name in mirrors), "domain_symmetry", "unknown symmetry set")
        _require(all(self.set_kind[names.index(name)] == "symmetry" for name in mirrors), "domain_symmetry",
                 "a declared symmetry set is not of kind symmetry")
        if domain["target"] == "box":
            _require(not walls and not mirrors and domain.get("ground") is None, "domain_truncation",
                     "a box has no truncation sets, symmetry sets or ground plane")
        else:
            ground = domain.get("ground")
            _require(bool(walls), "domain_truncation", "a half-space needs truncation sets")
            _require(isinstance(ground, dict) and ground.get("axis") in AXES and _number(ground.get("coordinate")), "domain_ground")

    # geometry -------------------------------------------------------------------------------

    def _geometry(self):
        x = self.nodes[self.tets[:, :4]]
        fv = np.array(FACET_VERTICES)
        a, b, c = x[:, fv[:, 0]], x[:, fv[:, 1]], x[:, fv[:, 2]]
        cross = np.cross(b - a, c - a)                     # twice the facet area times the outward normal
        twice = np.linalg.norm(cross, axis=2)
        edges = np.array([[i, j] for i in range(4) for j in range(i + 1, 4)])
        lengths = np.linalg.norm(x[:, edges[:, 0]] - x[:, edges[:, 1]], axis=2)
        self.x = x
        self.h = lengths.max(axis=1)
        self.volume = np.einsum("ei,ei->e", x[:, 1] - x[:, 0], np.cross(x[:, 2] - x[:, 0], x[:, 3] - x[:, 0])) / 6.0
        self.degenerate = bool(np.any(~(self.volume > 0.0)) or np.any(self.volume < TOLERANCES["TOL_GEOM_VOLUME"] * self.h ** 3)
                               or np.any(twice == 0.0))
        if self.degenerate:
            self.straightness = None
            return
        self.area = 0.5 * twice
        self.normal = cross / twice[:, :, None]
        self.grad = -cross / (6.0 * self.volume[:, None, None])      # grad L_i = -(A_i / (3 V)) n_i
        self.dn = np.einsum("mai,eic->emac", p2_vertex_gradient_table(), self.grad)   # grad N_a at vertex m
        edge = np.array(EDGES)
        mid = 0.5 * (x[:, edge[:, 0]] + x[:, edge[:, 1]])
        offset = np.linalg.norm(self.nodes[self.tets[:, 4:]] - mid, axis=2)
        self.straightness = float(np.max(offset / self.h[:, None]))
        self.x10 = np.concatenate([x, mid], axis=1)       # affine geometry: mid-edge nodes at chord midpoints
        self.facet_h = np.stack([np.max(np.stack([np.linalg.norm(x[:, v[i]] - x[:, v[j]], axis=1)
                                                  for i, j in FACET_EDGES]), axis=0) for v in FACET_VERTICES], axis=1)
        self.l_ref = float(np.linalg.norm(self.nodes.max(axis=0) - self.nodes.min(axis=0)))

    # strength and body force -----------------------------------------------------------------

    def _regions(self, arrays):
        regions = self.problem["regions"]
        _require(isinstance(regions, list) and regions, "regions")
        k, alpha, b0, b1 = [], [], [], []
        for region in regions:
            criterion = region.get("criterion")
            _require(criterion in ("von_mises", "drucker_prager"), "criterion", repr(criterion))
            _require(_number(region.get("k")), "strength")
            a = region.get("alpha", 0.0) if criterion == "drucker_prager" else 0.0
            if criterion == "von_mises":
                _require("alpha" not in region and region["k"] > 0.0, "strength", "von_mises")
            else:
                _require(_number(a) and 0.0 < a < 1.0 / math.sqrt(3.0) and region["k"] >= 0.0, "strength", "drucker_prager")
            for key in ("b0", "b1"):
                _require(isinstance(region.get(key), list) and len(region[key]) == 3 and all(_number(v) for v in region[key]),
                         "body_force")
            k.append(float(region["k"]))
            alpha.append(float(a))
            b0.append(region["b0"])
            b1.append(region["b1"])
        index = np.asarray(arrays["element_region"])
        _require(index.shape == (self.n_elements,) and index.min() >= 0 and index.max() < len(regions), "element_region")
        self.k, self.alpha = np.array(k)[index], np.array(alpha)[index]
        self.b0, self.b1 = np.array(b0, dtype=float)[index], np.array(b1, dtype=float)[index]
        self.region_k, self.region_alpha = np.array(k), np.array(alpha)
        self.region_b0, self.region_b1 = np.array(b0, dtype=float), np.array(b1, dtype=float)

    # facets ---------------------------------------------------------------------------------

    def _facets(self, arrays):
        tets, n_e = self.tets, self.n_elements
        fv = np.array(FACET_VERTICES)
        triples = np.sort(tets[:, fv], axis=2).reshape(-1, 3)
        order = np.lexsort((triples[:, 2], triples[:, 1], triples[:, 0]))
        ranked = triples[order]
        same = np.all(ranked[1:] == ranked[:-1], axis=1)
        _require(not np.any(same[1:] & same[:-1]), "mesh_nonconforming", "facet with more than two owners")
        first = np.nonzero(same)[0]
        one, two = order[first], order[first + 1]
        low, high = np.minimum(one, two), np.maximum(one, two)
        rank = np.argsort(low)
        low, high = low[rank], high[rank]
        self.int_minus = np.stack([low // 4, low % 4], axis=1)
        self.int_plus = np.stack([high // 4, high % 4], axis=1)
        paired = np.zeros(4 * n_e, dtype=bool)
        paired[low] = paired[high] = True
        exterior = np.nonzero(~paired)[0]
        listed = np.asarray(arrays["facets"])
        _require(listed.ndim == 2 and listed.shape[1] == 2, "array_shape", "facets")
        _require(np.array_equal(listed[:, 0] * 4 + listed[:, 1], exterior), "facets_mismatch",
                 "the facet list is not the sorted list of exterior facets")
        self.ext = listed.astype(np.int64)

        # Every mid-edge node is shared by the elements of its edge, so jumps can be matched by node id.
        edge = np.array(EDGES)
        pairs = np.sort(tets[:, edge], axis=2).reshape(-1, 2)
        mids = tets[:, 4:].reshape(-1)
        key = pairs[:, 0] * (self.n_nodes + 1) + pairs[:, 1]
        rank = np.argsort(key, kind="stable")
        key, mids = key[rank], mids[rank]
        _require(not np.any((key[1:] == key[:-1]) & (mids[1:] != mids[:-1])), "mesh_nonconforming", "mid-edge node")
        _require(len(np.unique(key)) == len(np.unique(mids)), "mesh_nonconforming", "mid-edge node reused")

        # Local nodes of e+ that match the six facet nodes of e-, in the facet-node order of e-.
        fn = np.array(FACET_NODES)
        e1, f1, e2 = self.int_minus[:, 0], self.int_minus[:, 1], self.int_plus[:, 0]
        self.int_local_minus = fn[f1]
        ids = tets[e1[:, None], self.int_local_minus]
        match = tets[e2][:, None, :] == ids[:, :, None]
        _require(bool(np.all(match.sum(axis=2) == 1)), "mesh_nonconforming", "facet nodes do not match")
        self.int_local_plus = np.argmax(match, axis=2)

    def _boundary(self, arrays):
        problem = self.problem
        sets = problem["boundary_sets"]
        _require(isinstance(sets, list) and sets, "boundary_sets")
        footings = problem["footings"]
        for footing in footings:
            _require(footing.get("interface_law") == "bonded", "interface_law", repr(footing.get("interface_law")))
        kinds = []
        for item in sets:
            _require(item.get("kind") in ("free", "fixed", "symmetry", "traction", "footing"), "boundary_kind", repr(item.get("kind")))
            if item["kind"] == "footing":
                _require(isinstance(item.get("footing"), int) and 0 <= item["footing"] < len(footings), "boundary_footing")
            kinds.append(item["kind"])
        self.set_kind = kinds
        self.facet_set = np.asarray(arrays["facet_set"]).astype(np.int64)
        n_b = len(self.ext)
        _require(self.facet_set.shape == (n_b,) and self.facet_set.min() >= 0 and self.facet_set.max() < len(sets), "facet_set")
        _require(set(self.facet_set.tolist()) == set(range(len(sets))), "boundary_set_empty")
        self.facet_kind = np.array(kinds)[self.facet_set]
        self.t0, self.t1 = np.asarray(arrays["t0"], dtype=float), np.asarray(arrays["t1"], dtype=float)
        _require(self.t0.shape == (n_b, 3, 3) and self.t1.shape == (n_b, 3, 3), "array_shape", "t0, t1")
        _require(np.all(np.isfinite(self.t0)) and np.all(np.isfinite(self.t1)), "array_non_finite")
        other = self.facet_kind != "traction"
        _require(not np.any(self.t0[other]) and not np.any(self.t1[other]), "traction_on_other_kind")
        self.n_footings = len(footings)
        self.x_ref = np.array([f["x_ref"] for f in footings], dtype=float).reshape(self.n_footings, 3)
        self.facet_footing = np.array([sets[s].get("footing", -1) if kinds[s] == "footing" else -1 for s in self.facet_set],
                                      dtype=np.int64)
        _require(set(self.facet_footing[self.facet_footing >= 0].tolist()) == set(range(self.n_footings)), "footing_without_set")
        if self.degenerate:
            return
        e, f = self.ext[:, 0], self.ext[:, 1]
        self.ext_normal, self.ext_area = self.normal[e, f], self.area[e, f]
        self.footing_area = np.array([self.ext_area[self.facet_footing == q].sum() for q in range(self.n_footings)])
        # Symmetry sets are planar: one normal and one plane (interface 5). The set normal is that of its first facet.
        self.set_normal = {}
        tol = TOLERANCES["TOL_GEOM_PLANE"]
        for s, kind in enumerate(kinds):
            if kind != "symmetry":
                continue
            rows = np.nonzero(self.facet_set == s)[0]
            reference = self.ext_normal[rows[0]]
            origin = self.x[e[rows[0]], FACET_VERTICES[f[rows[0]]][0]]
            points = self.x[e[rows][:, None], np.array(FACET_VERTICES)[f[rows]]]
            _require(np.all(np.linalg.norm(self.ext_normal[rows] - reference, axis=1) <= tol), "symmetry_not_planar", sets[s]["name"])
            _require(np.all(np.abs((points - origin) @ reference) <= tol * self.l_ref), "symmetry_not_planar", sets[s]["name"])
            self.set_normal[s] = reference


class Query:
    def __init__(self, model, query):
        n_f = model.n_footings
        self.id, self.form = query["id"], query.get("form")
        _require(self.form in ("ray", "fixed_components", "support"), "query_form", repr(self.form))
        self.is_ray = self.form != "support"
        self.restrained = np.zeros((n_f, 6), dtype=bool)
        for f, c in query["restrained"]:
            self.restrained[f, c] = True

        def table(key):
            return np.array(query[key], dtype=float).reshape(n_f, 6)

        self.f0, self.f1, self.direction = np.zeros((n_f, 6)), np.zeros((n_f, 6)), np.zeros((n_f, 6))
        if self.form == "ray":
            self.f0, self.f1 = table("F0"), table("F1")
        elif self.form == "fixed_components":
            self.f1 = table("direction")
            for f, c, value in query["fixed"]:
                _require(self.f1[f, c] == 0.0, "query_fixed_direction")
                self.f0[f, c] = value
        else:
            self.direction = table("direction")
            _require(not np.any(model.t1) and not np.any(model.region_b1), "query_support_reference_load")
        _require(not np.any(self.f0[self.restrained]) and not np.any(self.f1[self.restrained])
                 and not np.any(self.direction[self.restrained]), "query_restrained_load")


def supported_domain_claim(model, bound):
    """The strongest claim version 0 supports for one bound (interface 14)."""
    domain = model.problem["domain"]
    if domain["target"] == "box":
        return "box_stated"
    names = [item["name"] for item in model.problem["boundary_sets"]]
    axis, coordinate = AXES[domain["ground"]["axis"]], domain["ground"]["coordinate"]
    tol = TOLERANCES["TOL_GEOM_PLANE"] * model.l_ref
    e, f = model.ext[:, 0], model.ext[:, 1]
    points = model.x[e[:, None], np.array(FACET_VERTICES)[f]]
    on_ground = np.all(np.abs(points[:, :, axis] - coordinate) <= tol, axis=1)
    declared = set(domain["truncation_sets"]) | set(domain["symmetry_sets"])
    for row in np.nonzero(~on_ground)[0]:
        s = int(model.facet_set[row])
        if model.set_kind[s] in ("traction", "footing") or names[s] not in declared:
            return "unverified"       # a loaded or undeclared face off the ground plane: no coverage
    walls = [model.set_kind[names.index(name)] for name in domain["truncation_sets"]]
    if bound == "upper":
        return "halfspace_by_zero_extension" if all(kind == "fixed" for kind in walls) else "unverified"
    weightless = not np.any(model.region_b0) and not np.any(model.region_b1)
    if all(kind == "free" for kind in walls) and weightless and domain["exterior_load_free"] is True:
        return "halfspace_by_zero_extension"
    return "unverified"


# ---------------------------------------------------------------------------------------------
# Affine systems
# ---------------------------------------------------------------------------------------------

class _Rows:
    """Row blocks of A x = g in coordinate form, with a kind per row."""

    def __init__(self, n_unknowns):
        self.n, self.m = n_unknowns, 0
        self.rows, self.cols, self.values, self.rhs, self.kinds = [], [], [], [], []

    def block(self, kind, rhs, entries):
        rhs = np.asarray(rhs, dtype=float).reshape(-1)
        for local, cols, values in entries:
            local, cols, values = np.broadcast_arrays(np.asarray(local), np.asarray(cols), np.asarray(values, dtype=float))
            self.rows.append(local.reshape(-1) + self.m)
            self.cols.append(cols.reshape(-1))
            self.values.append(values.reshape(-1))
        self.rhs.append(rhs)
        self.kinds.extend([kind] * len(rhs))
        self.m += len(rhs)

    def build(self):
        if self.rows:
            rows, cols, values = (np.concatenate(item) for item in (self.rows, self.cols, self.values))
        else:
            rows, cols, values = np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64), np.zeros(0)
        matrix = sp.coo_matrix((values, (rows, cols)), shape=(self.m, self.n)).tocsr()
        matrix.sum_duplicates()
        matrix.eliminate_zeros()
        rhs = np.concatenate(self.rhs) if self.rhs else np.zeros(0)
        return matrix, rhs, np.array(self.kinds, dtype=object)


class LowerSystem:
    """Rows of interface 9 and 17.3: equilibrium, facet continuity, boundary tractions, footing resultants."""

    def __init__(self, model, query):
        self.model, self.query = model, query
        n_e = model.n_elements
        self.offset = 1 if query.is_ray else 0
        self.n = self.offset + 24 * n_e
        off, rows = self.offset, _Rows(self.n)
        e = np.arange(n_e)

        entries = []
        for c, pairs in enumerate(TRACTION_MAP):
            for v, j in pairs:
                for i in range(4):
                    entries.append((3 * e + c, off + 24 * e + 6 * i + v, model.grad[:, i, j]))
            if off:
                entries.append((3 * e + c, np.zeros(n_e, dtype=np.int64), model.b1[:, c]))
        rows.block("equilibrium", -model.b0.reshape(-1), entries)

        n_i = len(model.int_minus)
        if n_i:
            e1, f1, e2 = model.int_minus[:, 0], model.int_minus[:, 1], model.int_plus[:, 0]
            normal, facet = model.normal[e1, f1], np.arange(n_i)
            entries = []
            for c, pairs in enumerate(TRACTION_MAP):
                for v, j in pairs:
                    for slot in range(3):
                        row = 9 * facet + 3 * slot + c
                        entries.append((row, off + 24 * e1 + 6 * model.int_local_minus[:, slot] + v, normal[:, j]))
                        entries.append((row, off + 24 * e2 + 6 * model.int_local_plus[:, slot] + v, -normal[:, j]))
            rows.block("continuity", np.zeros(9 * n_i), entries)

        fv = np.array(FACET_VERTICES)
        for kind in ("free", "traction"):
            pick = np.nonzero(model.facet_kind == kind)[0]
            if len(pick) == 0:
                continue
            eb, fb, normal, facet = model.ext[pick, 0], model.ext[pick, 1], model.ext_normal[pick], np.arange(len(pick))
            entries = []
            for c, pairs in enumerate(TRACTION_MAP):
                for slot in range(3):
                    row = 9 * facet + 3 * slot + c
                    for v, j in pairs:
                        entries.append((row, off + 24 * eb + 6 * fv[fb, slot] + v, normal[:, j]))
                    if off and kind == "traction":
                        entries.append((row, np.zeros(len(pick), dtype=np.int64), -model.t1[pick, slot, c]))
            rows.block("boundary", model.t0[pick].reshape(-1) if kind == "traction" else np.zeros(9 * len(pick)), entries)

        pick = np.nonzero(model.facet_kind == "symmetry")[0]
        if len(pick):
            eb, fb, normal, facet = model.ext[pick, 0], model.ext[pick, 1], model.ext_normal[pick], np.arange(len(pick))
            tangents = np.array([tangent_pair(n) for n in normal])                    # (n, 2, 3)
            coefficient = np.einsum("ntc,ncv->ntv", tangents, traction_matrix(normal))  # (n, 2, 6)
            entries = []
            for slot in range(3):
                for t in range(2):
                    for v in range(6):
                        entries.append((6 * facet + 2 * slot + t, off + 24 * eb + 6 * fv[fb, slot] + v, coefficient[:, t, v]))
            rows.block("boundary", np.zeros(6 * len(pick)), entries)

        # Footing resultants (interface 6), all six components of every footing, as an operator on x.
        n_f = model.n_footings
        res_rows, res_cols, res_values = [], [], []
        weights = traction_moment_weights()
        for row in np.nonzero(model.facet_footing >= 0)[0]:
            q, eb, fb = int(model.facet_footing[row]), int(model.ext[row, 0]), int(model.ext[row, 1])
            area, nmat = model.ext_area[row], traction_matrix(model.ext_normal[row])
            arm = model.x[eb, list(FACET_VERTICES[fb])] - model.x_ref[q]
            lever = area * weights @ arm                    # lever[b] = A sum_a w_ab (x_a - x_ref)
            for b, vertex in enumerate(FACET_VERTICES[fb]):
                block = np.vstack([area * monomial_integral((1, 0, 0), 1.0) * nmat, skew(lever[b]) @ nmat])
                for c in range(6):
                    for v in range(6):
                        res_rows.append(6 * q + c)
                        res_cols.append(off + 24 * eb + 6 * vertex + v)
                        res_values.append(block[c, v])
        self.resultant = sp.coo_matrix((res_values, (res_rows, res_cols)), shape=(6 * n_f, self.n)).tocsr()
        if query.is_ray and n_f:
            free = np.nonzero(~query.restrained.reshape(-1))[0]
            block = self.resultant[free].tocoo()
            entries = [(block.row, block.col, block.data)]
            if np.any(query.f1.reshape(-1)[free]):
                entries.append((np.arange(len(free)), np.zeros(len(free), dtype=np.int64), -query.f1.reshape(-1)[free]))
            rows.block("resultant", query.f0.reshape(-1)[free], entries)
        self.matrix, self.rhs, self.kinds = rows.build()

    def pack(self, value, sigma):
        x = np.zeros(self.n)
        if self.offset:
            x[0] = value
        x[self.offset:] = np.asarray(sigma, dtype=float).reshape(-1)
        return x

    def sigma(self, x):
        return x[self.offset:].reshape(self.model.n_elements, 4, 6)

    def scales(self, s_ref, value):
        d = np.full(self.n, s_ref)
        if self.offset:
            d[0] = abs(value) if value != 0.0 else 1.0
        return d


class UpperSystem:
    """Rows of interface 17.4 step 2 for the continuous control or the broken jump formulation."""

    def __init__(self, model, query, formulation):
        self.model, self.query, self.broken = model, query, formulation == JUMP_ID
        n_e, n_f = model.n_elements, model.n_footings
        self.nv = 30 * n_e if self.broken else 3 * model.n_nodes
        self.n = self.nv + 6 * n_f
        local = np.arange(10)[None, :]
        element = np.arange(n_e)[:, None]
        self.key = (10 * element + local) if self.broken else model.tets          # (E, 10): velocity node of (e, a)
        self.n_keys = 10 * n_e if self.broken else model.n_nodes
        self.key_x = np.zeros((self.n_keys, 3))
        self.key_x[self.key.reshape(-1)] = model.x10.reshape(-1, 3)
        rows = _Rows(self.n)
        fn = np.array(FACET_NODES)
        facet_keys = self.key[model.ext[:, 0][:, None], fn[model.ext[:, 1]]]          # (Fb, 6)
        comp = np.arange(3)

        fixed = np.unique(facet_keys[model.facet_kind == "fixed"])
        if len(fixed):
            index = np.arange(3 * len(fixed))
            rows.block("fixed", np.zeros(3 * len(fixed)), [(index, (3 * fixed[:, None] + comp).reshape(-1), 1.0)])
        for s, kind in enumerate(model.set_kind):
            if kind != "symmetry":
                continue
            keys = np.setdiff1d(np.unique(facet_keys[model.facet_set == s]), fixed)
            if len(keys):
                index = np.repeat(np.arange(len(keys)), 3)
                rows.block("symmetry", np.zeros(len(keys)),
                           [(index, (3 * keys[:, None] + comp).reshape(-1), np.tile(model.set_normal[s], len(keys)))])
        self.footing_keys = []
        for q in range(n_f):
            keys = np.unique(facet_keys[model.facet_footing == q])
            self.footing_keys.append(keys)
            arm = self.key_x[keys] - model.x_ref[q]
            entries = []
            for c in range(3):
                row = 3 * np.arange(len(keys)) + c
                entries.append((row, 3 * keys + c, 1.0))
                entries.append((row, np.full(len(keys), self.nv + 6 * q + c), -1.0))
                for j in range(3):       # -(w x r)_c = (r x w)_c = sum_j skew(r)[c, j] w_j
                    if j != c:
                        entries.append((row, np.full(len(keys), self.nv + 6 * q + 3 + j),
                                        np.array([skew(r)[c, j] for r in arm])))
            rows.block("footing", np.zeros(3 * len(keys)), entries)
        restrained = np.nonzero(query.restrained.reshape(-1))[0]
        if len(restrained):
            rows.block("restrained", np.zeros(len(restrained)), [(np.arange(len(restrained)), self.nv + restrained, 1.0)])
        if not query.is_ray:
            free = np.nonzero(~query.restrained.reshape(-1))[0]
            rows.block("support", query.direction.reshape(-1)[free], [(np.arange(len(free)), self.nv + free, 1.0)])

        isochoric = np.nonzero(model.alpha == 0.0)[0]
        if len(isochoric):
            entries = []
            for m in range(4):
                for a in range(10):
                    for c in range(3):
                        entries.append((4 * np.arange(len(isochoric)) + m, 3 * self.key[isochoric, a] + c,
                                        model.dn[isochoric, m, a, c]))
            rows.block("isochory", np.zeros(4 * len(isochoric)), entries)

        self.bernstein = bernstein_map()
        self._normal_jump_rows(rows)
        self.matrix, self.rhs, self.kinds = rows.build()

    def _normal_jump_rows(self, rows):
        """n . c_B = 0 on all six Bernstein coefficients of every interior facet with tau = 0 (interface 12.3)."""
        model = self.model
        if not self.broken or len(model.int_minus) == 0:
            return
        e1, f1, e2 = model.int_minus[:, 0], model.int_minus[:, 1], model.int_plus[:, 0]
        tight = np.nonzero(model.alpha[e1] == 0.0)[0]
        if len(tight) == 0:
            return
        normal = model.normal[e1[tight], f1[tight]]
        minus = self.key[e1[tight][:, None], model.int_local_minus[tight]]
        plus = self.key[e2[tight][:, None], model.int_local_plus[tight]]
        entries = []
        for b in range(6):
            for p in np.nonzero(self.bernstein[b])[0]:
                for c in range(3):
                    row = 6 * np.arange(len(tight)) + b
                    entries.append((row, 3 * plus[:, p] + c, self.bernstein[b, p] * normal[:, c]))
                    entries.append((row, 3 * minus[:, p] + c, -self.bernstein[b, p] * normal[:, c]))
        rows.block("normal_jump", np.zeros(6 * len(tight)), entries)

    def pack(self, velocity, footing_velocity):
        x = np.zeros(self.n)
        x[: self.nv] = np.asarray(velocity, dtype=float).reshape(-1)
        x[self.nv:] = np.asarray(footing_velocity, dtype=float).reshape(-1)
        return x

    def element_field(self, x):
        return x[: self.nv].reshape(-1, 3)[self.key]

    def footing(self, x):
        return x[self.nv:].reshape(self.model.n_footings, 6)

    def scales(self, x):
        footing = self.footing(x)
        v_ref = max(float(np.max(np.abs(x[: self.nv]), initial=0.0)), float(np.max(np.abs(footing[:, :3]), initial=0.0)))
        if v_ref == 0.0:
            v_ref = 1.0
        d = np.full(self.n, v_ref)
        rotation = max(v_ref / self.model.ell, float(np.max(np.abs(footing[:, 3:]), initial=0.0)))
        d[self.nv:].reshape(-1, 6)[:, 3:] = rotation
        return d


def literal_zeros(matrix, rhs):
    """Unknowns fixed to zero by homogeneous rows with one remaining entry, to a fixed point (interface 17.1)."""
    pattern = matrix.copy()
    pattern.data = np.ones_like(pattern.data)
    homogeneous = (rhs == 0.0)
    zero = np.zeros(matrix.shape[1], dtype=bool)
    while True:
        alive = (~zero).astype(float)
        count = pattern @ alive
        single = np.nonzero(homogeneous & (count == 1.0))[0]
        if len(single) == 0:
            return zero
        column = pattern[single] @ (alive * np.arange(matrix.shape[1]))     # the one surviving column of each row
        new = np.unique(np.rint(column).astype(np.int64))
        if np.all(zero[new]):
            return zero
        zero[new] = True


def backward_error(matrix, rhs, x, scale):
    """omega of interface 17.1 over all rows. A row with zero denominator must have zero residual."""
    if matrix.shape[0] == 0:
        return 0.0
    residual = np.abs(matrix @ x - rhs)
    denominator = abs(matrix) @ scale + np.abs(rhs)
    if np.any((denominator == 0.0) & (residual != 0.0)):
        return math.inf
    active = denominator > 0.0
    return float(np.max(residual[active] / denominator[active], initial=0.0))


class Conditioning:
    """The reduced, scaled system A_hat of interface 17.1, its conditioning, and its minimum-norm solver."""

    def __init__(self, matrix, rhs, scale, zero, dense_entries_max=DENSE_ENTRIES_MAX, sparse_rows_max=SPARSE_ROWS_MAX,
                 block_entries_max=DENSE_ENTRIES_MAX):
        self.sparse_rows_max, self.block_entries_max = sparse_rows_max, block_entries_max
        keep_cols = np.nonzero(~zero)[0]
        reduced = (matrix[:, keep_cols] @ sp.diags(scale[keep_cols])).tocsr()
        row_max = abs(reduced).max(axis=1).toarray().reshape(-1) if reduced.shape[1] else np.zeros(reduced.shape[0])
        self.rows = np.nonzero(row_max > 0.0)[0]
        used = np.unique(reduced[self.rows].indices) if len(self.rows) else np.zeros(0, dtype=np.int64)
        self.cols = keep_cols[used]
        self.row_max = row_max[self.rows]
        self.a_hat = (sp.diags(1.0 / self.row_max) @ reduced[self.rows][:, used]).tocsr() if len(self.rows) else None
        self.scale = scale[self.cols]
        self.n_rows, self.n_cols = len(self.rows), len(self.cols)
        self.kappa, self.rank, self.sigma_max, self.sigma_min = None, None, None, None
        self.ritz_residual, self.skipped, self.refusal, self.proof, self.kappa_estimate = None, None, None, None, None
        self._dense, self._lu, self._normal_shift = None, None, None
        if self.n_rows == 0:
            self.method, self.kappa, self.rank, self.skipped = "no equality rows remain: kappa_eq is 1 by convention", 1.0, 0, 0
            self.proof, self.kappa_estimate = "there is no equality system", 1.0
        elif self.n_rows * self.n_cols <= dense_entries_max:
            self._factor_dense()
        else:
            self._factor_sparse()

    def _factor_dense(self):
        u, s, vt = scipy.linalg.svd(self.a_hat.toarray(), full_matrices=False)
        keep = s > TOLERANCES["TOL_RANK"] * s[0]
        self._dense = (u[:, keep], s[keep], vt[keep])
        self.rank, self.sigma_max, self.sigma_min = int(keep.sum()), float(s[0]), float(s[keep][-1])
        self.kappa = self.kappa_estimate = self.sigma_max / self.sigma_min
        self.skipped = self.n_rows - self.rank
        self.method = "dense singular value decomposition"
        self.proof = "every singular value is computed"

    def _refuse(self, detail):
        self.method, self.refusal, self.kappa = "unavailable", detail, None

    def _factor_sparse(self):
        """Shift-invert Lanczos on a normal matrix of A_hat with positive shifts (interface 17.1).

        kappa_eq is reported only as a bound this route proves, and the route refuses when it cannot prove one:

        * a first shift just above zero counts the eigenvalues below it (Sylvester, diagonal pivots only);
        * that count is never read as dependent rows on its own. It is accepted when it equals the zero eigenvalues
          forced by the shape of A_hat, or when a dense singular value decomposition of the rows involved, block by
          block, finds exactly that many singular values at or below TOL_RANK sigma_max. A row that is only nearly
          dependent fails this, and no kappa_eq is given;
        * a second shift with the same count shows that no eigenvalue lies between the two, so the smallest
          non-zero singular value is at least its square root;
        * sigma_max is bounded above by a norm bound, not by a Ritz value.

        A_hat A_hat^T and A_hat^T A_hat have the same non-zero eigenvalues. The interface names the first; the second
        is tried first when it is the sparser one, as for footing rows that all share the six footing unknowns.
        """
        counts = np.bincount(self.a_hat.indices, minlength=self.n_cols).astype(float)
        fill_rows, fill_cols = float(np.sum(counts ** 2)), float(np.sum(np.diff(self.a_hat.indptr).astype(float) ** 2))
        for by_rows in ((True,) if fill_rows <= fill_cols else (False, True)):
            try:
                if self._normal_route(by_rows):
                    return
            except (RuntimeError, MemoryError, spla.ArpackError) as error:
                self._refuse(f"sparse factorisation of {self.n_rows} rows by {self.n_cols} unknowns failed: {type(error).__name__}")
                if by_rows:
                    return

    def _normal_route(self, by_rows):
        """One normal matrix. Returns True when it settles the matter, False to hand over to the row matrix."""
        size, other = (self.n_rows, self.n_cols) if by_rows else (self.n_cols, self.n_rows)
        if size > self.sparse_rows_max:
            self._refuse(f"{self.n_rows} rows and {self.n_cols} unknowns exceed the sparse factorisation limit of "
                         f"{self.sparse_rows_max}")
            return by_rows
        normal = (self.a_hat @ self.a_hat.T if by_rows else self.a_hat.T @ self.a_hat).tocsc()
        identity = sp.identity(size, format="csc")

        def factor(shift):
            lu = spla.splu((normal - shift * identity).tocsc(), permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                           options={"SymmetricMode": True})
            # With diagonal pivots only, the negative pivots are the eigenvalues below the shift.
            count = int(np.sum(lu.U.diagonal() < 0.0)) if np.array_equal(lu.perm_r, lu.perm_c) else None
            return lu, count, spla.LinearOperator((size, size), matvec=lu.solve, dtype=float)

        top = float(spla.eigsh(normal, k=1, which="LA", return_eigenvectors=False, tol=1e-10)[0])     # from below
        bound = min(float(abs(normal).sum(axis=1).max()),
                    float(abs(self.a_hat).sum(axis=0).max()) * float(abs(self.a_hat).sum(axis=1).max()))
        low = SHIFT_REL * top
        lu, count, inverse = factor(low)
        if count is None:
            self._refuse("the factorisation pivoted off the diagonal, so the eigenvalues below the shift could not be counted")
            return True
        if count == max(0, size - other):
            proof = f"the {count} eigenvalues below the first shift are the zeros forced by the shape of A_hat"
        elif not by_rows:
            return False
        else:
            proved, blocks, largest = self._dependent_rows(lu, low, math.sqrt(top))
            if proved < 0:
                # One block too large to decompose: prove the count with an explicit basis of the null space instead.
                worst = self._null_space_bound(lu, low, count)
                if worst is None or not worst <= TOLERANCES["TOL_RANK"] * math.sqrt(top):
                    self._refuse(f"{count} eigenvalues of the normal matrix lie below the first shift; a row block of {largest} "
                                 "rows is too large for a dense decomposition, and an explicit null-space basis "
                                 + ("is too large as well" if worst is None else f"reaches only {worst / math.sqrt(top):.1e} sigma_max")
                                 + ": the dependent rows are not proved")
                    return True
                proved = count
                blocks, largest = 0, 0
            if proved != count:
                self._refuse(f"{count} eigenvalues of the normal matrix lie below the first shift, and the row blocks prove "
                             f"{proved} dependent rows at TOL_RANK: a singular value may lie between TOL_RANK sigma_max and "
                             f"{math.sqrt(SHIFT_REL):.0e} sigma_max, so no kappa_eq is proved")
                return True
            if blocks:
                proof = (f"{count} dependent rows proved by dense singular value decompositions of {blocks} row blocks, "
                         f"the largest of {largest} rows")
            else:
                proof = (f"{count} dependent rows proved by an orthonormal basis Y of {count} vectors with "
                         f"||A_hat^T y|| <= {worst / math.sqrt(top):.1e} sigma_max ||y|| for every y in its span")
        rough = float(spla.eigsh(normal, k=1, sigma=low, which="LA", OPinv=inverse, return_eigenvectors=False, tol=1e-6)[0])
        del lu, inverse                    # one factorisation in memory at a time
        shift, value = None, None
        for share in (0.9, 0.7, 0.5, 0.25, 0.1, 0.01, 0.001):       # the rough value may sit above the eigenvalue
            trial = share * rough
            if trial <= 2.0 * low:
                break
            lu_two, count_two, inverse_two = factor(trial)
            if count_two == count:
                shift = trial
                try:
                    pair = spla.eigsh(normal, k=1, sigma=shift, which="LA", OPinv=inverse_two, tol=1e-10)
                    z, value = pair[1][:, 0], float(pair[0][0])
                    self.ritz_residual = float(np.linalg.norm(normal @ z - value * z) / (abs(value) * np.linalg.norm(z)))
                except spla.ArpackError:
                    value = None
                break
            del lu_two, inverse_two        # an eigenvalue lies below this shift: look lower
        if shift is None:
            self._refuse(f"no second shift above {low:.3e} separates the smallest non-zero eigenvalue (first estimate {rough:.3e})")
            return True
        if value is None or not value >= shift:
            self._refuse("the Ritz pair above the second shift did not converge, so its residual cannot be reported")
            return True
        self._normal_shift = (normal, low) if by_rows else None
        self.refusal = None
        self.sigma_max = math.sqrt(top)
        self.sigma_min = math.sqrt(value) if value is not None and value > 0.0 else None
        self.kappa_estimate = self.sigma_max / self.sigma_min if self.sigma_min else None
        self.kappa = math.sqrt(bound / shift)
        self.rank = size - count
        self.skipped = self.n_rows - self.rank
        self.proof = proof
        self.method = (f"shift-invert Lanczos on {'A_hat A_hat^T' if by_rows else 'A_hat^T A_hat'} with Sylvester counts at two "
                       "shifts; kappa_eq is the proved upper bound sqrt(norm bound on lambda_max / second shift)")
        return True

    def _null_space_bound(self, lu, low, count):
        """The largest ||A_hat^T y|| over unit vectors of an explicit ``count``-dimensional subspace, or None.

        Subspace iteration with -low (M - low I)^-1 gives an orthonormal Y. If ||A_hat^T Y||_2 is at most
        TOL_RANK sigma_max, then A_hat has at least ``count`` singular values at or below TOL_RANK sigma_max
        (Courant-Fischer), and the Sylvester count says there are no more below the shift.
        """
        if (self.n_rows + self.n_cols) * count > NULL_BASIS_ENTRIES_MAX:
            return None
        basis = np.random.default_rng(1).standard_normal((self.n_rows, count))
        for _ in range(3):
            basis, triangle = np.linalg.qr(-low * lu.solve(basis))
            diagonal = np.abs(np.diag(triangle))
            if not diagonal.min() > 1e-8 * diagonal.max():
                return None                 # the iterate lost rank: fewer independent null vectors than counted
        image = self.a_hat.T @ basis
        return math.sqrt(max(0.0, float(scipy.linalg.eigvalsh(image.T @ image, subset_by_index=[count - 1, count - 1])[0])))

    def _dependent_rows(self, lu, low, sigma_max):
        """Dependent rows proved by dense decompositions of the row blocks that carry the null space.

        -low (M - low I)^-1 is, to rounding, the projector onto the null space of M = A_hat A_hat^T. Applied to random
        vectors it marks the rows that null vectors touch. Marked rows that share an unknown form a block; blocks share
        no unknown, so the null space splits over them. Each block is decomposed densely, and only singular values at
        or below TOL_RANK sigma_max count. Missing a row can lose a proof; it cannot create one.
        """
        from scipy.sparse.csgraph import connected_components
        rng = np.random.default_rng(0)
        marked = np.zeros(self.n_rows, dtype=bool)
        for _ in range(PROBES):
            probe = rng.standard_normal(self.n_rows)
            for _ in range(3):
                probe = -low * lu.solve(probe)
            size = float(np.max(np.abs(probe)))
            if size > 0.0:
                marked |= np.abs(probe) > PROBE_REL * size
        rows = np.nonzero(marked)[0]
        if len(rows) == 0:
            return 0, 0, 0
        sub = self.a_hat[rows].tocsr()
        pattern = sub.copy()
        pattern.data = np.ones_like(pattern.data)
        blocks, label = connected_components(pattern @ pattern.T, directed=False)
        order = np.argsort(label, kind="stable")
        ends = np.cumsum(np.bincount(label, minlength=blocks))
        proved, largest, start = 0, 0, 0
        for stop in ends:
            members = order[start:stop]
            start = stop
            piece = sub[members]
            cols = np.unique(piece.indices)
            if len(members) * len(cols) > self.block_entries_max:
                return -1, blocks, len(members)        # a block too large to decompose: nothing is proved
            values = scipy.linalg.svdvals(piece[:, cols].toarray())
            proved += len(members) - int(np.sum(values > TOLERANCES["TOL_RANK"] * sigma_max))
            largest = max(largest, len(members))
        return proved, blocks, largest

    def minimum_norm(self, residual_hat):
        """The minimum-norm correction of the scaled unknowns for one scaled residual."""
        if self._dense is not None:
            u, s, vt = self._dense
            return vt.T @ ((u.T @ residual_hat) / s)
        if self._normal_shift is not None:
            if self._lu is None:           # factorised again on first use: only one factorisation is held at a time
                normal, low = self._normal_shift
                self._lu = spla.splu((normal - low * sp.identity(normal.shape[0], format="csc")).tocsc(),
                                     permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0, options={"SymmetricMode": True})
            return self.a_hat.T @ self._lu.solve(residual_hat)
        return spla.lsmr(self.a_hat, residual_hat, atol=1e-15, btol=1e-15, conlim=1e14, maxiter=20 * max(self.a_hat.shape))[0]


CAPACITY_OPTIONS = ("dense_entries_max", "sparse_rows_max", "block_entries_max", "feasibility_dense_unknowns")


def _capacity(options):
    """The capacity options a caller may set. They bound cost and change no acceptance rule. No other option exists."""
    unknown = sorted(set(options) - set(CAPACITY_OPTIONS))
    if unknown:
        raise TypeError(f"unknown option {unknown}: only capacity limits can be set")
    return {"dense_entries_max": options.get("dense_entries_max", DENSE_ENTRIES_MAX),
            "sparse_rows_max": options.get("sparse_rows_max", SPARSE_ROWS_MAX),
            "block_entries_max": options.get("block_entries_max", DENSE_ENTRIES_MAX)}


def _conditioning(system, scale, zero, options):
    """The conditioning of one system, or a capacity refusal: a kappa_eq that is not proved is never used."""
    item = Conditioning(system.matrix, system.rhs, scale, zero, **_capacity(options))
    if item.kappa is None:
        raise CapacityRefusal(item.refusal)
    return item


def project(matrix, rhs, scale, x, zero, conditioning):
    """Literal zeros, then the scaled minimum-norm correction with iterative refinement (interface 17.3 step 2)."""
    start = x.copy()
    x = x.copy()
    x[zero] = 0.0
    omega = backward_error(matrix, rhs, x, scale)
    sweeps = 0
    if conditioning.n_rows:
        for _ in range(5):
            residual = (rhs - matrix @ x)[conditioning.rows] / conditioning.row_max
            trial = x.copy()
            trial[conditioning.cols] += conditioning.scale * conditioning.minimum_norm(residual)
            trial_omega = backward_error(matrix, rhs, trial, scale)
            if not trial_omega < omega:
                break
            x, omega, sweeps = trial, trial_omega, sweeps + 1
    move = float(np.max(np.abs(x - start) / scale, initial=0.0))
    return x, move, omega, sweeps


# ---------------------------------------------------------------------------------------------
# Yield and flow, by principal values
# ---------------------------------------------------------------------------------------------

def principal_yield(sigma, k, alpha):
    """sqrt(J2) + alpha I1 - k from the principal stresses. sigma (..., 6); k, alpha broadcast."""
    s = np.asarray(sigma, dtype=float)
    tensor = np.empty(s.shape[:-1] + (3, 3))
    tensor[..., 0, 0], tensor[..., 1, 1], tensor[..., 2, 2] = s[..., 0], s[..., 1], s[..., 2]
    tensor[..., 0, 1] = tensor[..., 1, 0] = s[..., 3]
    tensor[..., 1, 2] = tensor[..., 2, 1] = s[..., 4]
    tensor[..., 0, 2] = tensor[..., 2, 0] = s[..., 5]
    p = np.linalg.eigvalsh(tensor)
    j2 = ((p[..., 0] - p[..., 1]) ** 2 + (p[..., 1] - p[..., 2]) ** 2 + (p[..., 2] - p[..., 0]) ** 2) / 6.0
    return np.sqrt(j2) + alpha * p.sum(axis=-1) - k


def principal_rates(model, field):
    """theta, |e| and |D|_F at the four vertices of every element, from the principal strain rates."""
    gradient = np.einsum("eac,emad->emcd", field, model.dn)
    p = np.linalg.eigvalsh(0.5 * (gradient + np.swapaxes(gradient, 2, 3)))
    theta = p.sum(axis=2)
    deviator = p - theta[:, :, None] / 3.0
    return theta, np.sqrt(np.sum(deviator ** 2, axis=2)), np.sqrt(np.sum(p ** 2, axis=2))


def jump_coefficients(model, system, field):
    """Bernstein coefficients c_B (F, 6, 3) of v+ - v- on every interior facet, in the facet-node order of e-."""
    e1, e2 = model.int_minus[:, 0], model.int_plus[:, 0]
    nodal = field[e2[:, None], model.int_local_plus] - field[e1[:, None], model.int_local_minus]
    return np.einsum("bp,fpc->fbc", system.bernstein, nodal)


def jump_parts(model, coefficients):
    """Opening o = n . c and slip s = |c - o n| of every coefficient."""
    normal = model.normal[model.int_minus[:, 0], model.int_minus[:, 1]]
    opening = np.einsum("fbc,fc->fb", coefficients, normal)
    slip = np.linalg.norm(coefficients - opening[:, :, None] * normal[:, None, :], axis=2)
    return opening, slip


class UpperState:
    """Everything the Upper procedure evaluates on one field."""

    def __init__(self, model, system, query, x):
        self.field, self.footing = system.element_field(x), system.footing(x)
        self.theta, self.enorm, self.dnorm = principal_rates(model, self.field)
        frictional = model.alpha > 0.0
        self.flow = np.where(frictional[:, None], 3.0 * math.sqrt(2.0) * model.alpha[:, None] * self.enorm - self.theta, -math.inf)
        self.d_scale = float(np.max(self.dnorm, initial=0.0))
        self.jump = np.zeros((0, 6))
        self.coefficients = None
        if system.broken and len(model.int_minus):
            alpha = model.alpha[model.int_minus[:, 0]]
            _require(bool(np.all(12.0 * alpha ** 2 < 1.0)), "formulation_not_applicable")
            self.coefficients = jump_coefficients(model, system, self.field)
            self.opening, self.slip = jump_parts(model, self.coefficients)
            self.tau = 3.0 * alpha / np.sqrt(1.0 - 12.0 * alpha ** 2)
            self.jump = np.where((alpha > 0.0)[:, None], self.tau[:, None] * self.slip - self.opening, -math.inf)
            self.d_scale = max(self.d_scale, float(np.max(np.linalg.norm(self.coefficients, axis=2), initial=0.0)))
        self.signed = np.concatenate([self.flow.reshape(-1), self.jump.reshape(-1)])
        scale = self.d_scale if self.d_scale > 0.0 else 1.0
        self.rho_flow = max(0.0, float(np.max(self.flow, initial=-math.inf))) / scale
        self.rho_jump = max(0.0, float(np.max(self.jump, initial=-math.inf))) / scale

        # Dissipation by the value formulas of interface 4.3 and 12.2, with no test of the domain.
        alpha, k = model.alpha[:, None], model.k[:, None]
        safe = np.where(alpha > 0.0, alpha, 1.0)
        density = np.where(alpha > 0.0, k * self.theta / (3.0 * safe), math.sqrt(2.0) * k * self.enorm)
        self.dissipation = float(np.sum(0.25 * model.volume[:, None] * density))
        if self.coefficients is not None:
            e1, f1 = model.int_minus[:, 0], model.int_minus[:, 1]
            alpha, k = model.alpha[e1][:, None], model.k[e1][:, None]
            safe = np.where(alpha > 0.0, alpha, 1.0)
            density = np.where(alpha > 0.0, k * self.opening / (3.0 * safe), k * self.slip)
            self.dissipation += float(np.sum(model.area[e1, f1][:, None] / 6.0 * density))
        self.w0, _ = work_rate(model, self.field, self.footing, model.b0, model.t0, query.f0 if query.is_ray else None)
        self.w1, self.w_scale = work_rate(model, self.field, self.footing, model.b1, model.t1, query.f1)
        if query.is_ray:
            self.value = (self.dissipation - self.w0) / self.w1 if self.w1 != 0.0 else None
        else:
            self.value = self.dissipation - self.w0


def work_rate(model, field, footing, body, traction, load):
    """A work rate of interface 10 and the sum of the absolute values of its terms."""
    terms = [np.einsum("ec,ec->e", body, np.einsum("a,eac->ec", volume_weights(), field)) * model.volume]
    pick = np.nonzero(model.facet_kind == "traction")[0]
    if len(pick):
        e, f = model.ext[pick, 0], model.ext[pick, 1]
        facet_field = field[e[:, None], np.array(FACET_NODES)[f]]
        terms.append(model.ext_area[pick] * np.einsum("fic,ij,fjc->f", traction[pick], facet_work_table(), facet_field))
    if load is not None:
        terms.append((load * footing).reshape(-1))
    flat = np.concatenate([np.asarray(t, dtype=float).reshape(-1) for t in terms])
    return float(flat.sum()), float(np.abs(flat).sum())


# ---------------------------------------------------------------------------------------------
# The certificate procedure
# ---------------------------------------------------------------------------------------------

class Outcome:
    """What one run of the checker gives.

    ``check`` gives a certificate in the frozen format and its method report. ``final_check`` gives a verdict only:
    ``passed`` is True or False, an estimate certificate carries the code when it is False, and there is never a
    certified bound.
    """

    def __init__(self, certificate, method, report, fields=None, passed=None):
        self.certificate, self.method, self.report, self.fields, self.passed = certificate, method, report, fields or {}, passed

    @property
    def certified(self):
        return self.certificate is not None and self.certificate["classification"] == "certified_bound"

    @property
    def reasons(self):
        return [] if self.certificate is None else list(self.certificate["reasons"])


def method_report(certificate, report):
    """The record kept beside a certificate: how kappa_eq was obtained and which interior field repaired the bound.

    The frozen certificate has no key for these. The record is bound to its certificate by a hash of the whole
    certificate, so a method cannot be moved to another certificate.
    """
    def whole(key):
        value = report.get(key)
        return None if value is None else int(value)

    return {
        "format": FORMAT_METHOD, "format_version": FORMAT_VERSION,
        "problem_sha256": certificate["problem_sha256"], "result_sha256": certificate["result_sha256"],
        "certificate_sha256": document_hash(certificate, {}, ""), "checker": dict(certificate["checker"]),
        "query_form": report.get("query_form"),
        "kappa_eq": certificate["checks"]["kappa_eq"], "kappa_estimate": _clean(report.get("kappa_estimate")),
        "kappa_method": report.get("kappa_method"), "kappa_proof": report.get("kappa_proof"),
        "rows": whole("rows"), "unknowns": whole("unknowns"), "literal_zeros": whole("literal_zeros"),
        "reduced_rows": whole("reduced_rows"), "reduced_unknowns": whole("reduced_unknowns"), "rank": whole("rank"),
        "dependent_rows_skipped": whole("dependent_rows_skipped"), "ritz_residual": _clean(report.get("ritz_residual")),
        "sigma_max": _clean(report.get("sigma_max")), "sigma_min": _clean(report.get("sigma_min")),
        "refinement_sweeps": whole("refinement_sweeps"), "interior_field": report.get("interior_field", "none"),
    }


def method_failures(certificate, method):
    """Why a method report does not support its certificate at a gate. Empty means it does."""
    if not isinstance(method, dict) or method.get("format") != FORMAT_METHOD or method.get("format_version") != FORMAT_VERSION:
        return ["certificate_method_missing"]
    failures = []
    if method.get("certificate_sha256") != document_hash(certificate, {}, ""):
        failures.append("certificate_sha256")
    for key in ("problem_sha256", "result_sha256", "checker"):
        if method.get(key) != certificate.get(key):
            failures.append(key)
    if method.get("query_form") not in ("ray", "fixed_components", "support"):
        failures.append("query_form")
    kappa_method = method.get("kappa_method")
    if not isinstance(kappa_method, str) or not kappa_method or kappa_method == "unavailable":
        failures.append("kappa_method")
    if not isinstance(method.get("kappa_proof"), str) or not method.get("kappa_proof"):
        failures.append("kappa_proof")
    if method.get("kappa_eq") != certificate["checks"].get("kappa_eq"):
        failures.append("kappa_eq")
    counts = [method.get(key) for key in ("reduced_rows", "rank", "dependent_rows_skipped")]
    if not all(isinstance(item, int) and not isinstance(item, bool) and item >= 0 for item in counts) or counts[0] - counts[1] != counts[2]:
        failures.append("dependent_rows_skipped")
    if "ritz_residual" not in method or (method["ritz_residual"] is not None and not _number(method["ritz_residual"])):
        failures.append("ritz_residual")
    elif isinstance(kappa_method, str) and kappa_method.startswith("shift-invert") and method["ritz_residual"] is None:
        failures.append("ritz_residual")
    interior = method.get("interior_field")
    repaired = _number(certificate["checks"].get("repair_mu")) and certificate["checks"]["repair_mu"] > 0.0
    if interior not in INTERIOR_FIELDS or (repaired and interior == "none"):
        failures.append("interior_field")
    return failures


def _clean(value):
    if value is None:
        return None
    value = float(value)
    return value if math.isfinite(value) else None


class _Run:
    """State shared by the Lower and Upper procedures: identity, checks, reasons."""

    def __init__(self, problem, result, model, bound, requested_claim, options):
        self.problem, self.result, self.model, self.bound = problem, result, model, bound
        self.requested_claim, self.options = requested_claim, options
        self.checks = {key: None for key in CHECK_KEYS}
        self.checks["tau_eq"] = TAU_EQ
        self.report = {"checker": CHECKER_NAME, "kappa_method": None}
        self.value_certified = None
        self.fields = {}
        self.final_only = False

    def passed(self):
        """A final check that found nothing: a verdict, with no certificate."""
        return Outcome(None, None, self.report, self.fields, passed=True)

    def stop(self, code):
        _require(code in REASON_CODES, "reason_code", code)
        return self.finish([code])

    def finish(self, reasons):
        model, result, checks = self.model, self.result, self.checks
        candidate = result.get("value")
        certified = not reasons
        _require(not (certified and self.final_only), "checker_internal", "a final check never certifies")
        claim = "box_stated"
        if model is not None and not model.degenerate:
            claim = supported_domain_claim(model, self.bound)
            if certified and self.requested_claim is not None and self.requested_claim != claim:
                reasons, certified = ["unsupported_domain_claim"], False
        floor = None
        if certified:
            kappa, omega = checks["kappa_eq"], checks["omega_eq"]
            floor = max(TOLERANCES["FLOOR_MIN_REL"], kappa * omega) * abs(self.value_certified)
        # A producer that is this checker's own code cannot be certified by it for a gate.
        own = isinstance(result.get("producer"), dict) and result["producer"].get("name") in (CHECKER_NAME, PROCEDURE_DEMO_CHECKER)
        name = PROCEDURE_DEMO_CHECKER if own else CHECKER_NAME
        document = {
            "format": FORMAT_CERTIFICATE, "format_version": FORMAT_VERSION,
            "problem_sha256": self.problem.get("problem_sha256"), "result_sha256": result.get("result_sha256"),
            "query_id": result.get("query_id"), "bound": result.get("bound"), "formulation": result.get("formulation"),
            "classification": "certified_bound" if certified else "estimate", "reasons": list(reasons),
            "value_candidate": candidate if _number(candidate) else None,
            "value_certified": _clean(self.value_certified) if certified else None, "floor": _clean(floor),
            "domain_claim": claim,
            "criteria": sorted({region.get("criterion") for region in self.problem.get("regions", [])}),
            "flow_rule": "associated",
            "interface_laws": sorted({item.get("interface_law") for item in self.problem.get("footings", [])}),
            "checks": {key: _clean(checks[key]) for key in CHECK_KEYS},
            "checker": {"name": name, "version": CHECKER_VERSION},
            "tolerances": dict(TOLERANCES),
        }
        method = method_report(document, self.report)
        if certified:
            failures = certified_check_failures(document, self.report.get("query_form")) + method_failures(document, method)
            _require(not failures, "checker_internal", f"a certified bound must show {failures}")
        return Outcome(document, method, self.report, self.fields, passed=False if self.final_only else None)

    def conditioning(self, system, scale, zero):
        item = _conditioning(system, scale, zero, self.options)
        self.checks["kappa_eq"] = item.kappa
        self.report.update(kappa_method=item.method, kappa_proof=item.proof, kappa_estimate=item.kappa_estimate, rank=item.rank,
                           reduced_rows=item.n_rows, reduced_unknowns=item.n_cols, rows=int(system.matrix.shape[0]),
                           unknowns=int(system.n), literal_zeros=int(zero.sum()), sigma_max=item.sigma_max,
                           sigma_min=item.sigma_min, ritz_residual=item.ritz_residual, dependent_rows_skipped=item.skipped)
        return item

    def recheck_equalities(self, system, x, scale, zero, conditioning):
        """Interface 17.3 step 6 and 17.4 step 7, the equality part. Returns a reason code or None."""
        if np.any(x[zero] != 0.0):
            return "literal_zero_violated"
        omega = backward_error(system.matrix, system.rhs, x, scale)
        self.checks["omega_eq"] = omega
        if not omega <= TAU_EQ:
            return "equality_residual_above_floor"
        kappa = conditioning.kappa
        if kappa is None:
            raise CapacityRefusal("kappa_eq is not available")
        self.checks["kappa_eq"] = kappa
        if not kappa <= TOLERANCES["TOL_KAPPA_CAP"] or not kappa * omega <= TOLERANCES["TOL_FLOOR_MAX"]:
            return "equality_ill_conditioned"
        return None


def certified_check_failures(document, query_form):
    """The numbers a certified_bound must show (interface 17.7), read back from the finished document."""
    checks, t = document["checks"], TOLERANCES
    need = ["omega_eq", "kappa_eq", "tau_eq", "projection_move", "repair_mu", "straightness_max"]
    if document["bound"] == "lower":
        need.append("rho_cone")
    else:
        need.append("rho_flow")
        if document["formulation"] == JUMP_ID:
            need.append("rho_jump_cone")
        if query_form in ("ray", "fixed_components"):
            need.append("work_reference")
    missing = [key for key in need if not _number(checks.get(key))]
    if missing or not _number(document["value_certified"]) or not _number(document["floor"]):
        return missing or ["value_certified"]
    omega, kappa = checks["omega_eq"], checks["kappa_eq"]
    tests = {
        "omega_eq": 0.0 <= omega <= checks["tau_eq"], "tau_eq": checks["tau_eq"] == TAU_EQ,
        "kappa_eq": 1.0 <= kappa <= t["TOL_KAPPA_CAP"] and kappa * omega <= t["TOL_FLOOR_MAX"],
        "projection_move": 0.0 <= checks["projection_move"] <= t["TOL_PROJECTION_MOVE"],
        "repair_mu": 0.0 <= checks["repair_mu"] <= t["TOL_REPAIR_MAX"],
        "straightness_max": checks["straightness_max"] <= t["TOL_GEOM_STRAIGHT"],
    }
    if document["bound"] == "lower":
        tests["rho_cone"] = checks["rho_cone"] <= 0.0
    else:
        tests["rho_flow"] = checks["rho_flow"] <= t["TOL_CONE_MARGIN"]
        if document["formulation"] == JUMP_ID:
            tests["rho_jump_cone"] = checks["rho_jump_cone"] <= t["TOL_CONE_MARGIN"]
        if "work_reference" in need:
            tests["work_reference"] = checks["work_reference"] > 0.0
    failures = [key for key, passed in tests.items() if not passed]
    expected = max(t["FLOOR_MIN_REL"], kappa * omega) * abs(document["value_certified"])
    if not math.isclose(document["floor"], expected, rel_tol=1e-9, abs_tol=1e-300):
        failures.append("floor")
    return failures


def gate_admissible(record, method=None):
    """Only this checker's certified bound with its method report, or its feasible verdict with a witness, can
    satisfy a gate. A certificate without the report that says how kappa_eq was obtained is rejected."""
    if record["checker"]["name"] != CHECKER_NAME:
        return False
    if record.get("format") == FORMAT_FEASIBILITY:
        return record["status"] == "formulation_feasible" and record["witness_result_sha256"] is not None
    if record["classification"] != "certified_bound" or method_failures(record, method):
        return False
    return not certified_check_failures(record, method["query_form"])


def _identity(problem, arrays, result, result_arrays):
    """Step 1, the identity part. Returns (reason or None, model, query). A malformed document is refused."""
    try:
        return _identify(problem, arrays, result, result_arrays)
    except (KeyError, IndexError, TypeError, AttributeError) as error:
        raise Refusal("malformed_document", repr(error)) from error


def _identify(problem, arrays, result, result_arrays):
    _require(result.get("format") == FORMAT_RESULT and result.get("format_version") == FORMAT_VERSION, "format")
    bound = result.get("bound")
    _require(bound in ("lower", "upper"), "bound", repr(bound))
    model = Model(problem, arrays)
    if problem.get("problem_sha256") != document_hash(problem, arrays, "problem_sha256"):
        return "hash_mismatch", model, None
    if result.get("problem_sha256") != problem["problem_sha256"]:
        return "hash_mismatch", model, None
    if result.get("result_sha256") != document_hash(result, result_arrays, "result_sha256"):
        return "hash_mismatch", model, None
    queries = [query for query in problem["queries"] if query.get("id") == result.get("query_id")]
    if len(queries) != 1:
        return "query_mismatch", model, None
    formulation = result.get("formulation")
    if formulation not in ADMITTED[bound]:
        return "formulation_not_admitted", model, None      # interface 8: any other id, prefixed or not
    if formulation == JUMP_ID and np.any(12.0 * model.region_alpha ** 2 >= 1.0):
        return "formulation_not_applicable", model, None
    if model.degenerate:
        return "degenerate_element", model, None
    if model.straightness > TOLERANCES["TOL_GEOM_STRAIGHT"]:
        return "curved_element", model, None
    query = Query(model, queries[0])
    _result_arrays(model, result, result_arrays)
    if result.get("solver_status") not in ("solved", "almost_solved") or not _number(result.get("value")):
        return "no_candidate", model, query
    return None, model, query


def _result_arrays(model, result, arrays):
    """Interface 15.3: a result holds the arrays of its own bound and formulation, and no others."""
    n_e, n_n, n_f = model.n_elements, model.n_nodes, model.n_footings
    broken = result["formulation"] == JUMP_ID
    velocity, other = ("element_velocity", "velocity") if broken else ("velocity", "element_velocity")
    shape = (n_e, 10, 3) if broken else (n_n, 3)
    lower = {"sigma": (n_e, 4, 6), "resultant": (n_f, 6), "interior_sigma": (n_e, 4, 6), "equilibrium_multiplier": (n_e, 3)}
    upper = {velocity: shape, "footing_velocity": (n_f, 6), f"interior_{velocity}": shape, "interior_footing_velocity": (n_f, 6),
             "strain_rate": (n_e, 4, 6), "dissipation": (n_e, 4)}
    own, foreign = (lower, upper) if result["bound"] == "lower" else (upper, lower)
    wrong = (set(foreign) | {other, f"interior_{other}"}) & set(arrays)
    _require(not wrong, "result_wrong_bound_arrays", str(sorted(wrong)))
    _require(set(arrays) <= set(own), "array_names", str(sorted(set(arrays) - set(own))))
    for name in arrays:
        array = np.asarray(arrays[name])
        _require(array.dtype.str == "<f8" and array.shape == own[name] and bool(np.all(np.isfinite(array))), "array_shape", name)
    _require(("interior_sigma" in arrays) == ("interior_value" in result), "result_interior")
    _require((f"interior_{velocity}" in arrays) == ("interior_footing_velocity" in arrays), "result_interior")


def _field(arrays, name, shape):
    _require(name in arrays, "missing_array", name)
    array = np.asarray(arrays[name])
    _require(array.dtype.str == "<f8" and array.shape == tuple(shape) and bool(np.all(np.isfinite(array))), "array_shape", name)
    return array.astype(float)


def check(problem, arrays, result, result_arrays, requested_claim=None, **options):
    """Run interface 17.3 or 17.4 on one result. Returns an Outcome whose certificate is in the frozen format.

    ``requested_claim`` names a domain claim the caller wants: one the problem does not support gives an estimate.
    A system beyond the checker's capacity raises ``CapacityRefusal``: no certificate of either kind is written.
    """
    return _run(problem, arrays, result, result_arrays, requested_claim, False, options)


def final_check(problem, arrays, result, result_arrays, **options):
    """The recheck of interface 17.3 step 6 or 17.4 step 7 on the stored field, with no projection and no repair.

    It never certifies. ``passed`` says whether the stored field meets every final condition; when it does not, the
    Outcome carries an estimate with the one code.
    """
    return _run(problem, arrays, result, result_arrays, None, True, options)


def _run(problem, arrays, result, result_arrays, requested_claim, final_only, options):
    _capacity(options)
    reason, model, query = _identity(problem, arrays, result, result_arrays)
    run = _Run(problem, result, model, result["bound"], requested_claim, options)
    run.final_only = final_only
    if model is not None and not model.degenerate:
        run.checks["straightness_max"] = model.straightness
    if reason is not None:
        return run.stop(reason)
    run.report["query_form"] = query.form
    if result["bound"] == "lower":
        return _check_lower(run, model, query, result, result_arrays, final_only)
    return _check_upper(run, model, query, result, result_arrays, final_only)


def _reference_stress(model, query, value, sigma, resultant):
    """S_ref of interface 17.7. Body forces and footing forces enter by their Euclidean norms."""
    load = abs(value) if query.is_ray else 0.0
    terms = [float(np.max(model.region_k)), float(np.max(np.abs(model.t0), initial=0.0)),
             load * float(np.max(np.abs(model.t1), initial=0.0)),
             model.l_ref * float(np.max(np.linalg.norm(model.region_b0, axis=1) + load * np.linalg.norm(model.region_b1, axis=1)))]
    if sigma is not None:
        terms.append(float(np.max(np.abs(sigma), initial=0.0)))
    for q in range(model.n_footings):
        force = np.linalg.norm(query.f0[q, :3]) + load * np.linalg.norm(query.f1[q, :3])
        if resultant is not None:
            force = max(force, float(np.linalg.norm(resultant[q, :3])))
        terms.append(force / model.footing_area[q])
    return max(terms) or 1.0


def _lower_value(system, query, x):
    if query.is_ray:
        return float(x[0])
    return float(query.direction.reshape(-1) @ (system.resultant @ x))


def _check_lower(run, model, query, result, arrays, final_check):
    t, checks = TOLERANCES, run.checks
    n_e, n_f = model.n_elements, model.n_footings
    sigma = _field(arrays, "sigma", (n_e, 4, 6))
    stored = _field(arrays, "resultant", (n_f, 6))
    value = float(result["value"])
    system = LowerSystem(model, query)
    x = system.pack(value, sigma)

    # Step 1: the value and the resultant are those of the stored field.
    resultant = (system.resultant @ x).reshape(n_f, 6)
    s_ref = _reference_stress(model, query, value, sigma, resultant)
    run.report["s_ref"] = s_ref
    if n_f:
        size = s_ref * model.footing_area[:, None] * np.array([1.0] * 3 + [model.ell] * 3)
        checks["rho_resultant"] = float(np.max(np.abs(resultant - stored) / size))
        if np.any(np.abs(resultant - stored) > t["TOL_VALUE_MATCH"] * size):
            return run.stop("value_mismatch")
    if not query.is_ray:
        recomputed = _lower_value(system, query, x)
        if abs(recomputed - value) > t["TOL_VALUE_MATCH"] * abs(value):
            return run.stop("value_mismatch")
    residual = np.abs(system.matrix @ x - system.rhs)

    def worst(kind, weight=None):
        rows = system.kinds == kind
        if not np.any(rows):
            return None
        return float(np.max(residual[rows] * (1.0 if weight is None else weight)))

    checks["rho_div"] = worst("equilibrium", np.repeat(model.h, 3)) / s_ref
    for key, kind in (("rho_jump", "continuity"), ("rho_bc", "boundary")):
        checks[key] = None if worst(kind) is None else worst(kind) / s_ref

    scale = system.scales(s_ref, value)
    zero = literal_zeros(system.matrix, system.rhs)
    conditioning = run.conditioning(system, scale, zero)

    def cone(vector):
        return principal_yield(system.sigma(vector), model.k[:, None], model.alpha[:, None])

    if final_check:
        checks["projection_move"], checks["repair_mu"] = 0.0, 0.0
    else:
        # Step 2: exact affine projection.
        x, move, omega, run.report["refinement_sweeps"] = project(system.matrix, system.rhs, scale, x, zero, conditioning)
        checks["projection_move"], checks["omega_eq"] = move, omega
        if move > t["TOL_PROJECTION_MOVE"]:
            return run.stop("projection_move_too_large")
        # Step 3: yield from principal stresses.
        g = cone(x)
        checks["rho_cone"] = float(np.max(g)) / s_ref
        if checks["rho_cone"] > t["TOL_CONE_ACCEPT"]:
            return run.stop("cone_violation_too_large")
        # Steps 4 and 5: admissible repair with a strictly admissible field.
        margin = t["TOL_CONE_MARGIN"] * s_ref
        checks["repair_mu"] = 0.0
        if np.any(g > -margin):
            interior = _lower_interior(run, model, query, system, zero, conditioning, result, arrays)
            if isinstance(interior, str):
                return run.stop(interior)
            g_int = cone(interior)
            if not np.all(g_int <= -t["TOL_INTERIOR_MIN"] * s_ref):
                return run.stop("interior_field_invalid")
            need = g > -margin
            mu = max(0.0, float(np.max((g[need] + margin) / (g[need] - g_int[need]))))
            checks["repair_mu"] = mu
            if mu > t["TOL_REPAIR_MAX"]:
                return run.stop("repair_too_large")
            x = (1.0 - mu) * x + mu * interior

    # Step 6: the final field, from scratch.
    final_value = _lower_value(system, query, x)
    final_scale = system.scales(_reference_stress(model, query, final_value, system.sigma(x),
                                                  (system.resultant @ x).reshape(n_f, 6)), final_value)
    if float(np.max(np.abs(final_scale / scale - 1.0))) > 1e-6:
        conditioning = run.conditioning(system, final_scale, zero)
    reason = run.recheck_equalities(system, x, final_scale, zero, conditioning)
    if reason:
        return run.stop(reason)
    checks["rho_cone"] = float(np.max(cone(x))) / float(final_scale[-1])
    if checks["rho_cone"] > 0.0:
        return run.stop("cone_violation_too_large")
    run.value_certified = final_value
    run.fields = {"sigma": system.sigma(x).copy(), "resultant": (system.resultant @ x).reshape(n_f, 6)}
    return run.passed() if final_check else run.finish([])


def _lower_interior(run, model, query, system, zero, conditioning, result, arrays):
    """The strictly admissible field of interface 17.3 step 5, projected, or a reason code."""
    t = TOLERANCES
    homogeneous = not np.any(model.region_b0) and not np.any(model.t0) and not np.any(query.f0)
    if homogeneous and np.all(model.region_k > 0.0):
        run.report["interior_field"] = "zero"
        return np.zeros(system.n)
    if "interior_sigma" not in arrays:
        return "no_interior_field"
    run.report["interior_field"] = "supplied"
    sigma = np.asarray(arrays["interior_sigma"], dtype=float)
    value = result.get("interior_value")
    if sigma.shape != (model.n_elements, 4, 6) or not np.all(np.isfinite(sigma)) or not _number(value):
        return "interior_field_invalid"
    x = system.pack(value, sigma)
    if not query.is_ray and abs(_lower_value(system, query, x) - value) > t["TOL_VALUE_MATCH"] * max(abs(value), 1e-300):
        return "interior_field_invalid"
    s_ref = _reference_stress(model, query, value, sigma, (system.resultant @ x).reshape(-1, 6))
    scale = system.scales(s_ref, value)
    own = _conditioning(system, scale, zero, run.options)
    x, move, omega, _ = project(system.matrix, system.rhs, scale, x, zero, own)
    if move > t["TOL_PROJECTION_MOVE"] or not omega <= TAU_EQ:
        return "interior_field_invalid"
    return x


def _check_upper(run, model, query, result, arrays, final_check):
    t, checks = TOLERANCES, run.checks
    n_e, n_f = model.n_elements, model.n_footings
    system = UpperSystem(model, query, result["formulation"])
    name = "element_velocity" if system.broken else "velocity"
    shape = (n_e, 10, 3) if system.broken else (model.n_nodes, 3)
    x = system.pack(_field(arrays, name, shape), _field(arrays, "footing_velocity", (n_f, 6)))
    value = float(result["value"])

    # Step 1: the value is that of the stored field, by the formulas alone.
    state = UpperState(model, system, query, x)
    if state.value is None:
        checks["work_reference"] = state.w1
        return run.stop("work_nonpositive")
    if abs(state.value - value) > t["TOL_VALUE_MATCH"] * abs(state.value):
        return run.stop("value_mismatch")
    scale = system.scales(x)
    run.report["v_ref"] = float(scale[0]) if system.nv else None
    residual = np.abs(system.matrix @ x - system.rhs)
    velocity_rows = np.isin(system.kinds, ("fixed", "symmetry", "footing", "restrained", "support"))
    checks["rho_bc"] = float(np.max(residual[velocity_rows], initial=0.0)) / float(scale[0]) if np.any(velocity_rows) else None

    zero = literal_zeros(system.matrix, system.rhs)
    conditioning = run.conditioning(system, scale, zero)
    if final_check:
        checks["projection_move"], checks["repair_mu"] = 0.0, 0.0
    else:
        # Step 2: exact affine projection.
        x, move, omega, run.report["refinement_sweeps"] = project(system.matrix, system.rhs, scale, x, zero, conditioning)
        checks["projection_move"], checks["omega_eq"] = move, omega
        if move > t["TOL_PROJECTION_MOVE"]:
            return run.stop("projection_move_too_large")
        # Step 3: the signed conditions with the checker's own operators.
        state = UpperState(model, system, query, x)
        checks["rho_flow"], checks["rho_jump_cone"] = state.rho_flow, (state.rho_jump if system.broken else None)
        if state.rho_flow > t["TOL_FLOW_ACCEPT"] or state.rho_jump > t["TOL_FLOW_ACCEPT"]:
            return run.stop("flow_violation_too_large")
        # Step 4: repair only by an admissible construction.
        checks["repair_mu"] = 0.0
        repair = state.signed > t["TOL_CONE_MARGIN"] * state.d_scale
        if np.any(repair):
            interior = _upper_interior(run, model, query, system, zero, arrays, name, shape, repair)
            if isinstance(interior, str):
                return run.stop(interior)
            x_int, signed_int = interior
            mu = float(np.max(state.signed[repair] / (state.signed[repair] - signed_int[repair])))
            checks["repair_mu"] = mu
            if mu > t["TOL_REPAIR_MAX"]:
                return run.stop("repair_too_large")
            x = (1.0 - mu) * x + mu * x_int
        state = UpperState(model, system, query, x)

    # Step 5: reference work.
    if query.is_ray:
        checks["work_reference"] = state.w1
        if not (state.w1 > 0.0 and state.w1 >= t["TOL_WORK_MIN"] * state.w_scale):
            return run.stop("work_nonpositive")
    # Steps 6 and 7: the value and the final field, from scratch.
    final_scale = system.scales(x)
    if float(np.max(np.abs(final_scale / scale - 1.0))) > 1e-6:
        conditioning = run.conditioning(system, final_scale, zero)
    reason = run.recheck_equalities(system, x, final_scale, zero, conditioning)
    if reason:
        return run.stop(reason)
    checks["rho_flow"], checks["rho_jump_cone"] = state.rho_flow, (state.rho_jump if system.broken else None)
    if state.rho_flow > t["TOL_CONE_MARGIN"] or state.rho_jump > t["TOL_CONE_MARGIN"]:
        return run.stop("flow_violation_too_large")
    run.value_certified = state.value
    run.report.update(d_scale=state.d_scale, dissipation=state.dissipation, work_fixed=state.w0)
    run.fields = {name: (state.field if system.broken else x[: system.nv].reshape(-1, 3)).copy(), "footing_velocity": state.footing.copy()}
    return run.passed() if final_check else run.finish([])


def _upper_interior(run, model, query, system, zero, arrays, name, shape, repair):
    """The relative-interior field of interface 17.4 step 4, projected and checked, or a reason code."""
    t = TOLERANCES
    if f"interior_{name}" not in arrays or "interior_footing_velocity" not in arrays:
        return "no_admissible_correction"
    velocity = np.asarray(arrays[f"interior_{name}"], dtype=float)
    footing = np.asarray(arrays["interior_footing_velocity"], dtype=float)
    if velocity.shape != tuple(shape) or footing.shape != (model.n_footings, 6):
        return "interior_field_invalid"
    if not (np.all(np.isfinite(velocity)) and np.all(np.isfinite(footing))):
        return "interior_field_invalid"
    x = system.pack(velocity, footing)
    scale = system.scales(x)
    own = _conditioning(system, scale, zero, run.options)
    x, move, omega, _ = project(system.matrix, system.rhs, scale, x, zero, own)
    if move > t["TOL_PROJECTION_MOVE"] or not omega <= TAU_EQ:
        return "interior_field_invalid"
    state = UpperState(model, system, query, x)
    if not state.d_scale > 0.0:
        return "interior_field_invalid"
    if np.any(state.signed > t["TOL_CONE_MARGIN"] * state.d_scale):
        return "interior_field_invalid"
    if np.any(state.signed[repair] > -t["TOL_INTERIOR_MIN"] * state.d_scale):
        return "interior_field_invalid"
    run.report["interior_field"] = "supplied"
    return x, state.signed


# ---------------------------------------------------------------------------------------------
# Brackets (interface 15.6 and 17.9)
# ---------------------------------------------------------------------------------------------

def bracket(lower, upper, lower_method=None, upper_method=None):
    """The bracket record of two gate-admissible bounds of one problem and query.

    Each certificate comes with its method report: an estimate, a certificate under a reserved checker name, or one
    without its report forms no bracket.
    """
    _require(lower["bound"] == "lower" and upper["bound"] == "upper", "bracket_labels")
    _require(lower["classification"] == upper["classification"] == "certified_bound", "bracket_uncertified")
    _require(gate_admissible(lower, lower_method) and gate_admissible(upper, upper_method), "bracket_not_gate_admissible")
    _require(lower["problem_sha256"] == upper["problem_sha256"], "hash_mismatch")
    _require(lower["query_id"] == upper["query_id"], "query_mismatch")
    low, high = lower["value_certified"], upper["value_certified"]
    floors = lower["floor"] + upper["floor"]
    status = "bracket" if low <= high else ("order_within_floor" if low <= high + floors else "contradiction")
    both = lower["domain_claim"] == upper["domain_claim"] == "halfspace_by_zero_extension"
    return {
        "format": FORMAT_BRACKET, "format_version": FORMAT_VERSION, "problem_sha256": lower["problem_sha256"],
        "query_id": lower["query_id"], "lower_formulation": lower["formulation"], "upper_formulation": upper["formulation"],
        "lower_result_sha256": lower["result_sha256"], "upper_result_sha256": upper["result_sha256"],
        "lower": low, "upper": high, "floor_lower": lower["floor"], "floor_upper": upper["floor"], "status": status,
        "gap": (high - low) / (0.5 * (high + low)) if status == "bracket" and low > 0.0 else None,
        "domain_claim": "halfspace_by_zero_extension" if both else "box_stated",
        "lower_domain_claim": lower["domain_claim"], "upper_domain_claim": upper["domain_claim"],
    }


def write_certificate(directory, result, outcome):
    """Write ``<stem>.certificate.json`` (frozen format) and ``<stem>.certificate.method.json`` beside it."""
    _require(outcome.certificate is not None and outcome.method is not None, "no_certificate", "a final check writes nothing")
    path = Path(directory) / f"{record_stem(result)}.certificate.json"
    write_document(path, outcome.certificate)
    write_document(path.with_name(f"{record_stem(result)}.certificate.method.json"), outcome.method)
    return path


def read_certificate(path):
    """A certificate and its method report, or None for the report when the file beside it is missing."""
    path = Path(path)
    certificate, _ = load_document(path)
    beside = path.with_name(path.name[: -len(".json")] + ".method.json")
    return certificate, (load_document(beside)[0] if beside.exists() else None)


def read_gate_certificate(path):
    """A certificate for a gate: refused unless it is certified, real, and comes with a matching method report."""
    certificate, method = read_certificate(path)
    _require(method is not None, "certificate_method_missing", Path(path).name)
    failures = method_failures(certificate, method)
    _require(not failures, "certificate_method_invalid", str(failures))
    _require(gate_admissible(certificate, method), "certificate_not_gate_admissible", Path(path).name)
    return certificate, method


def write_feasibility(directory, result, record):
    path = Path(directory) / f"{record_stem(result)}.feasibility.json"
    write_document(path, record)
    return path


def write_bracket(directory, lower_result, upper_result, record):
    low, high = record_stem(lower_result), record_stem(upper_result)
    head = low[: -len(lower_result["formulation"])]
    _require(low.endswith(lower_result["formulation"]) and high == head + upper_result["formulation"], "bracket_names")
    path = Path(directory) / f"{head}{lower_result['formulation']}__{upper_result['formulation']}.bracket.json"
    write_document(path, record)
    return path


# ---------------------------------------------------------------------------------------------
# Kinematic feasibility (interface 11.1 and 15.5)
# ---------------------------------------------------------------------------------------------

def kinematic_residual(model, system, state):
    """Numerator and denominator of rho(v, V) of interface 11.1 at one field."""
    w_vertex = 0.25 * model.volume[:, None]
    frictional = (model.alpha > 0.0)[:, None]
    c_vertex = np.where(frictional, np.maximum(0.0, np.where(frictional, state.flow, 0.0)), np.abs(state.theta))
    top = float(np.sum(w_vertex * c_vertex ** 2))
    bottom = float(np.sum(w_vertex * state.dnorm ** 2))
    if state.coefficients is not None:
        e1, f1 = model.int_minus[:, 0], model.int_minus[:, 1]
        w_facet = (model.area[e1, f1] / (6.0 * model.facet_h[e1, f1]))[:, None]
        loose = (model.alpha[e1] > 0.0)[:, None]
        c_facet = np.where(loose, np.maximum(0.0, np.where(loose, state.jump, 0.0)), np.abs(state.opening))
        top += float(np.sum(w_facet * c_facet ** 2))
        bottom += float(np.sum(w_facet * np.sum(state.coefficients ** 2, axis=2)))
    return math.sqrt(top), math.sqrt(bottom)


def _weighted_rows(model, system):
    """Sparse C and G with ||C x||^2 and ||G x||^2 the squared numerator and denominator at alpha = 0."""
    n_e = model.n_elements
    root_v = np.sqrt(0.25 * model.volume)
    numerator, denominator = _Rows(system.n), _Rows(system.n)
    e = np.arange(n_e)
    entries = []
    for m in range(4):
        for a in range(10):
            for c in range(3):
                entries.append((4 * e + m, 3 * system.key[:, a] + c, root_v * model.dn[:, m, a, c]))
    numerator.block("theta", np.zeros(4 * n_e), entries)
    # |D|_F^2 = sum_c D_cc^2 + 2 sum_{c<d} D_cd^2, with D the symmetric part of the velocity gradient.
    entries, row = [], 0
    for m in range(4):
        for c in range(3):
            for d in range(c, 3):
                index = 24 * e + row
                for a in range(10):
                    if c == d:
                        entries.append((index, 3 * system.key[:, a] + c, root_v * model.dn[:, m, a, c]))
                    else:
                        half = root_v * math.sqrt(2.0) * 0.5
                        entries.append((index, 3 * system.key[:, a] + c, half * model.dn[:, m, a, d]))
                        entries.append((index, 3 * system.key[:, a] + d, half * model.dn[:, m, a, c]))
                row += 1
    denominator.block("rate", np.zeros(24 * n_e), entries)
    if system.broken and len(model.int_minus):
        e1, f1, e2 = model.int_minus[:, 0], model.int_minus[:, 1], model.int_plus[:, 0]
        n_i = len(e1)
        root_f = np.sqrt(model.area[e1, f1] / (6.0 * model.facet_h[e1, f1]))
        normal = model.normal[e1, f1]
        minus, plus = system.key[e1[:, None], model.int_local_minus], system.key[e2[:, None], model.int_local_plus]
        top, bottom = [], []
        for b in range(6):
            for p in np.nonzero(system.bernstein[b])[0]:
                weight = root_f * system.bernstein[b, p]
                for c in range(3):
                    top.append((6 * np.arange(n_i) + b, 3 * plus[:, p] + c, weight * normal[:, c]))
                    top.append((6 * np.arange(n_i) + b, 3 * minus[:, p] + c, -weight * normal[:, c]))
                    bottom.append((18 * np.arange(n_i) + 3 * b + c, 3 * plus[:, p] + c, weight))
                    bottom.append((18 * np.arange(n_i) + 3 * b + c, 3 * minus[:, p] + c, -weight))
        numerator.block("opening", np.zeros(6 * n_i), top)
        denominator.block("jump", np.zeros(18 * n_i), bottom)
    return numerator.build()[0], denominator.build()[0]


def _admissible_basis(system, query):
    """x = x_p + Z y for every (v, V) meeting the boundary and footing rows of A, without the normalisation.

    Three exact steps keep it sparse: literal zeros leave; a row holding an unknown that no other row holds gives that
    unknown in terms of the rest (every footing node follows the footing); the few rows left are solved densely.
    """
    rows = np.isin(system.kinds, ("fixed", "symmetry", "footing", "restrained", "support"))
    matrix, rhs = system.matrix[np.nonzero(rows)[0]].tocsr(), system.rhs[rows]
    n = system.n
    zero = literal_zeros(matrix, rhs)
    matrix = (matrix @ sp.diags((~zero).astype(float))).tocsr()
    matrix.eliminate_zeros()
    live = np.nonzero(np.diff(matrix.indptr) > 0)[0]
    _require(not np.any(rhs[np.diff(matrix.indptr) == 0]), "boundary_rows_inconsistent")
    matrix, rhs = matrix[live], rhs[live]

    # Rows with a private unknown: that unknown is a slave of the others in its row.
    count = np.bincount(matrix.indices, minlength=n)
    coo = matrix.tocoo()
    private = count[coo.col] == 1
    slave_of = np.full(matrix.shape[0], -1, dtype=np.int64)
    pivot = np.zeros(matrix.shape[0])
    for row, col, value in zip(coo.row[private][::-1], coo.col[private][::-1], coo.data[private][::-1]):
        slave_of[row], pivot[row] = col, value          # the first private entry of each row wins
    with_slave = np.nonzero(slave_of >= 0)[0]
    is_slave = np.zeros(n, dtype=bool)
    is_slave[slave_of[with_slave]] = True
    master = ~zero & ~is_slave
    master_index = np.cumsum(master) - 1
    n_m = int(master.sum())

    # The rows left over constrain masters only: a particular solution and a null space, densely.
    rest = np.nonzero(slave_of < 0)[0]
    particular_m = np.zeros(n_m)
    blocks = [sp.identity(n_m, format="csr")]
    if len(rest):
        block = matrix[rest]
        touched = np.unique(block.indices)
        u, sv, vt = scipy.linalg.svd(block[:, touched].toarray(), full_matrices=True)
        rank = int(np.sum(sv > TOLERANCES["TOL_RANK"] * sv[0])) if len(sv) else 0
        local = vt[:rank].T @ ((u[:, :rank].T @ rhs[rest]) / sv[:rank])
        _require(float(np.max(np.abs(block[:, touched] @ local - rhs[rest]), initial=0.0))
                 <= 1e-9 * max(1.0, float(np.max(np.abs(rhs[rest]), initial=0.0))), "boundary_rows_inconsistent")
        particular_m[master_index[touched]] = local
        null = vt[rank:].T
        null[np.abs(null) < 1e-14] = 0.0
        untouched = np.setdiff1d(np.arange(n_m), master_index[touched])
        keep = sp.coo_matrix((np.ones(len(untouched)), (untouched, np.arange(len(untouched)))), shape=(n_m, len(untouched)))
        span = sp.coo_matrix((null.reshape(-1), (np.repeat(master_index[touched], null.shape[1]),
                                                 np.tile(np.arange(null.shape[1]), len(touched)))), shape=(n_m, null.shape[1]))
        blocks = [sp.hstack([keep, span]).tocsr()]
    masters = blocks[0]

    # Lift to all unknowns: masters in place, slaves from their rows, literal zeros at zero.
    lift = sp.coo_matrix((np.ones(n_m), (np.nonzero(master)[0], np.arange(n_m))), shape=(n, n_m)).tocsr()
    offset = np.zeros(n)
    if len(with_slave):
        slave_rows = (sp.diags(-1.0 / pivot[with_slave]) @ matrix[with_slave]).tocsr()
        on_masters = slave_rows @ lift                                   # the slave's own entry drops out here
        place = sp.coo_matrix((np.ones(len(with_slave)), (slave_of[with_slave], np.arange(len(with_slave)))),
                              shape=(n, len(with_slave))).tocsr()
        lift = lift + place @ on_masters
        offset[slave_of[with_slave]] = rhs[with_slave] / pivot[with_slave]
    particular = lift @ particular_m + offset
    basis = (lift @ masters).tocsr()
    check_rows = system.matrix[np.nonzero(rows)[0]]
    _require(float(np.max(np.abs(check_rows @ particular - system.rhs[rows]), initial=0.0)) <= 1e-9 * max(1.0, float(np.max(np.abs(system.rhs[rows]), initial=0.0)))
             and float(abs(check_rows @ basis).max() if basis.shape[1] and check_rows.shape[0] else 0.0) <= 1e-9, "boundary_rows_inconsistent")
    return particular, basis


def _lsmr(matrix, rhs):
    """LSMR with its stop status. A run that did not converge is a capacity refusal, never a number."""
    out = spla.lsmr(matrix, rhs, atol=1e-14, btol=1e-14, conlim=1e14, maxiter=50 * max(matrix.shape))
    status = {"istop": int(out[1]), "iterations": int(out[2]), "residual": float(out[3]), "normal_residual": float(out[4])}
    if status["istop"] not in (0, 1, 2, 4, 5):
        raise CapacityRefusal(f"LSMR stopped with status {status['istop']} after {status['iterations']} iterations on a "
                              f"{matrix.shape[0]} by {matrix.shape[1]} system")
    return out[0], status


def _least_seminorm(operator, particular, basis, reference):
    """min ||operator x|| over x = particular + basis y, with reference^T x = 1 when a reference is given.

    Returns ``upper``, the value at an explicit admissible point, and ``lower``, a bound from below: the minimum of a
    converged least-squares problem, or the dual value 1 / ||z|| when reference = operator^T z holds to
    LSMR_COMPATIBLE, and 0 otherwise. The LSMR stop statuses and residuals are returned with them.
    """
    reduced = (operator @ basis).tocsr()
    if reference is None:
        offset = operator @ particular
        y, status = _lsmr(reduced, -offset)
        value = float(np.linalg.norm(offset + reduced @ y))
        return {"upper": value, "lower": value, "lsmr": [status]}
    w = basis.T @ reference
    size = float(np.linalg.norm(w))
    _require(size > 0.0, "work_nonpositive", "no admissible field does reference work")
    # w = reduced^T z + y_n, with y_n the part of w in the null space of the operator: a free mechanism.
    z, first = _lsmr(reduced.T.tocsr(), w)
    mechanism = w - reduced.T @ z
    first["relative_residual"] = float(np.linalg.norm(mechanism)) / size
    back, second = _lsmr(reduced, z)
    candidates = []
    if float(w @ mechanism) > 1e-12 * size * size:
        candidates.append(float(np.linalg.norm(reduced @ mechanism)) / float(w @ mechanism))
    if float(w @ back) > 0.0:
        candidates.append(float(np.linalg.norm(reduced @ back)) / float(w @ back))
    _require(bool(candidates), "feasibility_minimum")
    upper = min(candidates)
    # With w = reduced^T z, w^T y = 1 gives 1 <= ||z|| ||reduced y||.
    size_z = float(np.linalg.norm(z))
    lower = 1.0 / size_z if first["relative_residual"] <= LSMR_COMPATIBLE and size_z > 0.0 else 0.0
    return {"upper": upper, "lower": min(lower, upper), "lsmr": [first, second]}


def feasibility(problem, arrays, witness, witness_arrays, **options):
    """The feasibility record of interface 15.5 for the formulation of one witness result.

    Returns (record, report). The record is in the frozen format. Readings this implementation had to choose,
    each reported to the contract owner:

    * the witness residual is evaluated on the field the certificate procedure ends with, and only when that
      procedure certifies the result; the stored field's own residual goes into the report;
    * ``rank`` is that of the equality rows E of interface 11.1, without the query normalisation rows, with each
      row divided by its largest entry and rotations scaled by the length scale;
    * ``rho_feas`` is not below the true ratio: its numerator is an explicit admissible point and its denominator a
      bound from below. ``formulation_infeasible`` is given only when the opposite pair, a numerator from below
      over a denominator at an explicit point, is locked too. LSMR stop statuses and residuals go into ``method``;
    * ``footing_modes`` is ``null`` for a problem without footings, where the count does not apply.
    """
    _capacity(options)
    reason, model, query = _identity(problem, arrays, witness, witness_arrays)
    _require(witness.get("bound") == "upper", "feasibility_bound")
    _require(reason in (None, "no_candidate"), reason or "")
    formulation = witness["formulation"]
    system = UpperSystem(model, query, formulation)
    frictional = bool(np.any(model.alpha > 0.0))
    t = TOLERANCES
    report = {"query_form": query.form}

    # The witness is the stored field after the certificate procedure's projection and repair (17.4 steps 2 to 4).
    # It counts only when that procedure certifies it, so a certified Upper is a witness and nothing weaker is.
    witness_rho, witness_parts, stored_parts, sha = None, None, None, None
    if reason is None:
        name = "element_velocity" if system.broken else "velocity"
        outcome = check(problem, arrays, witness, witness_arrays, **options)
        report["witness_reasons"] = outcome.reasons
        stored = system.pack(witness_arrays[name], witness_arrays["footing_velocity"])
        stored_parts = kinematic_residual(model, system, UpperState(model, system, query, stored))
        report["stored_field_residual"] = stored_parts[0] / stored_parts[1] if stored_parts[1] > 0.0 else None
        if outcome.certified:
            x = system.pack(outcome.fields[name], outcome.fields["footing_velocity"])
            state = UpperState(model, system, query, x)
            witness_parts = kinematic_residual(model, system, state)
            witness_work = state.w1 if query.is_ray else 1.0
            if witness_parts[1] > 0.0:
                witness_rho, sha = witness_parts[0] / witness_parts[1], witness["result_sha256"]

    rank = nullity = modes = None
    equality = system.matrix[np.nonzero(system.kinds != "support")[0]]       # E of 11.1: the normalisation is not in it
    n_unknowns, n_rows = int(system.n), int(equality.shape[0])
    if n_unknowns <= options.get("feasibility_dense_unknowns", FEAS_DENSE_UNKNOWNS):
        scale = np.ones(system.n)
        scale[system.nv:].reshape(-1, 6)[:, 3:] = 1.0 / model.ell
        scaled = (equality @ sp.diags(scale)).toarray()
        norms = np.max(np.abs(scaled), axis=1) if n_rows else np.zeros(0)
        scaled = scaled[norms > 0.0] / norms[norms > 0.0][:, None]
        if scaled.shape[0]:
            _, s, vt = scipy.linalg.svd(scaled, full_matrices=True)
            rank = int(np.sum(s > t["TOL_RANK"] * s[0]))
            kernel = vt[rank:].T
        else:
            rank, kernel = 0, np.eye(system.n)
        nullity = n_unknowns - rank
        if model.n_footings:
            part = kernel[system.nv:]
            values = scipy.linalg.svdvals(part) if part.size else np.zeros(0)
            modes = int(np.sum(values > t["TOL_RANK"] * max(1.0, float(values[0])))) if len(values) else 0

    if frictional:
        method = "frictional: rho_feas is the residual of the witness after the checker's projection and repair"
        if witness_rho is not None:
            numerator, denominator = witness_parts
        else:
            _require(stored_parts is not None, "no_candidate", "no witness field to evaluate")
            numerator, denominator = stored_parts
    else:
        top, bottom = _weighted_rows(model, system)
        particular, basis = _admissible_basis(system, query)
        reference = _work_functional(model, system, query) if query.is_ray else None
        above, below = _least_seminorm(top, particular, basis, reference), _least_seminorm(bottom, particular, basis, reference)
        report["least_squares"] = {"numerator": above, "denominator": below}
        numerator, denominator = above["upper"], below["lower"]
        if witness_rho is not None:
            # The witness, scaled to the query normalisation, is a point of the admissible set: the minimum is not above it.
            numerator = min(numerator, witness_parts[0] / witness_work)
        if not denominator > 0.0:
            raise CapacityRefusal("the least-squares minimum of the denominator of rho_feas has no bound from below")
        locked_ratio = above["lower"] / below["upper"]
        runs = above["lsmr"] + below["lsmr"]
        method = ("two linear least-squares problems over the admissible set by LSMR: numerator at an explicit point, "
                  "denominator bounded below by duality; stop statuses " + ",".join(str(item["istop"]) for item in runs)
                  + "; largest relative residual of the dual systems "
                  + f"{max(item.get('relative_residual', 0.0) for item in runs):.1e}")
    _require(denominator > 0.0, "feasibility_denominator", "a mechanism without strain rate or jump")
    rho = numerator / denominator
    passes = witness_rho is not None and witness_rho <= t["TOL_FEAS_OK"]
    if rho <= t["TOL_FEAS_OK"] and passes and (modes is None or modes >= 1):
        status = "formulation_feasible"
    elif not frictional and rho >= t["TOL_FEAS_LOCKED"]:
        if not locked_ratio >= t["TOL_FEAS_LOCKED"]:
            raise CapacityRefusal(f"rho_feas is {rho:.3e} at explicit points, but the least-squares minima bound it below only "
                                  f"by {locked_ratio:.3e}: the locked verdict is not established")
        status = "formulation_infeasible"
    else:
        status = "inconclusive"
    own = isinstance(witness.get("producer"), dict) and witness["producer"].get("name") in (CHECKER_NAME, PROCEDURE_DEMO_CHECKER)
    record = {
        "format": FORMAT_FEASIBILITY, "format_version": FORMAT_VERSION, "problem_sha256": problem["problem_sha256"],
        "query_id": query.id, "formulation": formulation, "status": status, "rho_feas": rho, "numerator": numerator,
        "denominator": denominator, "n_unknowns": n_unknowns, "n_rows": n_rows, "rank": rank, "nullity": nullity,
        "footing_modes": modes, "method": method, "witness_result_sha256": sha, "witness_residual": witness_rho,
        "checker": {"name": PROCEDURE_DEMO_CHECKER if own else CHECKER_NAME, "version": CHECKER_VERSION},
        "tolerances": dict(TOLERANCES),
    }
    return record, report


def _work_functional(model, system, query):
    """W1 of interface 10 as a vector over the unknowns (v, V)."""
    out = np.zeros(system.n)
    weights = volume_weights()
    for c in range(3):
        np.add.at(out, (3 * system.key + c).reshape(-1), (model.volume[:, None] * weights[None, :] * model.b1[:, c][:, None]).reshape(-1))
    pick = np.nonzero(model.facet_kind == "traction")[0]
    if len(pick):
        e, f = model.ext[pick, 0], model.ext[pick, 1]
        keys = system.key[e[:, None], np.array(FACET_NODES)[f]]
        nodal = model.ext_area[pick][:, None, None] * np.einsum("fic,ij->fjc", model.t1[pick], facet_work_table())
        for c in range(3):
            np.add.at(out, (3 * keys + c).reshape(-1), nodal[:, :, c].reshape(-1))
    out[system.nv:] += query.f1.reshape(-1)
    return out


def conditioning_report(problem, arrays, query_id, formulation, **options):
    """kappa_eq of one formulation on one problem at nominal scales, before any result exists.

    Nominal scales: every stress scale and the multiplier scale are 1; every velocity scale is 1 and the rotation
    scale is 1 / length_scale. A certificate recomputes kappa_eq with the scales of its own candidate.
    """
    model = Model(problem, arrays)
    _require(not model.degenerate, "degenerate_element")
    query = Query(model, [item for item in problem["queries"] if item["id"] == query_id][0])
    if formulation == LOWER_ID:
        system = LowerSystem(model, query)
        scale = system.scales(1.0, 1.0)
    else:
        _require(formulation in ADMITTED["upper"], "formulation_not_admitted", repr(formulation))
        system = UpperSystem(model, query, formulation)
        scale = system.scales(np.zeros(system.n))
    zero = literal_zeros(system.matrix, system.rhs)
    item = Conditioning(system.matrix, system.rhs, scale, zero, **_capacity(options))
    return {"formulation": formulation, "elements": model.n_elements, "unknowns": int(system.n), "rows": int(system.matrix.shape[0]),
            "literal_zeros": int(zero.sum()), "reduced_rows": item.n_rows, "reduced_unknowns": item.n_cols, "rank": item.rank,
            "kappa_eq": item.kappa, "kappa_estimate": item.kappa_estimate, "sigma_max": item.sigma_max, "sigma_min": item.sigma_min,
            "method": item.method, "proof": item.proof, "ritz_residual": item.ritz_residual,
            "dependent_rows_skipped": item.skipped, "capacity_refusal": item.refusal}


# ---------------------------------------------------------------------------------------------
# Command: kappa_eq on the gate meshes (interface 17.1)
# ---------------------------------------------------------------------------------------------

GATE_CLASSES = ("upper_gate_3d", "one_sided_reference")


def conditioning_table(case_names=None, meshes=None, out=None, echo=None, **options):
    """kappa_eq at nominal scales for every binding formulation of the gate cases, on every listed mesh.

    The cases are built by the contract's own case builder. A system beyond the checker's capacity is recorded with
    its ``capacity_refusal`` and no number. Returns the document, and writes it to ``out`` after every system.
    """
    import datetime
    import time
    import test_limit_contract as contract      # the frozen case parser and builder; nothing else is taken from it

    cases = contract.parse_cases()
    if case_names is None:
        case_names = [name for name, case in cases.items() if case.get("status") == "active" and case.get("class") in GATE_CLASSES]
    document = {
        "format": FORMAT_CONDITIONING, "format_version": FORMAT_VERSION,
        "checker": {"name": CHECKER_NAME, "version": CHECKER_VERSION}, "date": datetime.date.today().isoformat(),
        "scales": "nominal: stress and multiplier scales 1, velocity scale 1, rotation scale 1 / length_scale",
        "capacity": {key: float(value) for key, value in _capacity(options).items()}, "rows": [],
    }
    for name in case_names:
        _require(name in cases, "case_unknown", name)
        case = cases[name]
        for mesh in (range(len(case["mesh"])) if meshes is None else meshes):
            problem, arrays = contract.build_case_problem(case, mesh)
            for formulation in case["formulation"]:
                started = time.time()
                try:
                    row = conditioning_report(problem, arrays, case["queries"][0]["id"], formulation, **options)
                except MemoryError:
                    row = {"formulation": formulation, "kappa_eq": None, "capacity_refusal": "out of memory while assembling the system"}
                row.update(case=name, mesh=int(mesh), problem_sha256=problem["problem_sha256"], seconds=round(time.time() - started, 1))
                document["rows"].append({key: row[key] for key in sorted(row)})
                if echo:
                    echo(f"{name} m{mesh} {formulation}: kappa_eq {row.get('kappa_eq')} rank {row.get('rank')} "
                         f"{row.get('capacity_refusal') or row.get('proof')} ({row['seconds']} s)")
                if out:
                    write_document(out, document)
    return document


def main(argv=None):
    """``python tests/limit_certificate.py conditioning --out FILE [--cases a,b] [--meshes 0,1]``"""
    import argparse
    import sys

    parser = argparse.ArgumentParser(description="Limit-analysis certificate checker (WP-237c).")
    commands = parser.add_subparsers(dest="command", required=True)
    table = commands.add_parser("conditioning", help="kappa_eq at nominal scales on the gate meshes, with capacity refusals")
    table.add_argument("--out", required=True, help="JSON file to write")
    table.add_argument("--cases", help="comma-separated case names; default: the eight true-3D cases and the strips")
    table.add_argument("--meshes", help="comma-separated mesh indices; default: every listed mesh")
    table.add_argument("--dense-entries-max", type=float, default=DENSE_ENTRIES_MAX)
    table.add_argument("--sparse-rows-max", type=float, default=SPARSE_ROWS_MAX)
    args = parser.parse_args(argv)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    document = conditioning_table(
        case_names=args.cases.split(",") if args.cases else None,
        meshes=[int(item) for item in args.meshes.split(",")] if args.meshes else None, out=args.out,
        echo=lambda line: print(line, flush=True), dense_entries_max=args.dense_entries_max, sparse_rows_max=args.sparse_rows_max)
    refused = sum(1 for row in document["rows"] if row["kappa_eq"] is None)
    print(f"{len(document['rows'])} systems, {refused} capacity refusals, written to {args.out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
