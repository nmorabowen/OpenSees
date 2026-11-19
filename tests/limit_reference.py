"""WP-237b: the Python reference assembler for limit analysis, version 0.

Written from ``tests/limit_interface.md``, ``tests/limit_cases.txt`` and the executable contract
``tests/test_limit_contract.py`` (case parser, problem builder, hash and validators) alone. It covers:

* the geometry of straight tet10 meshes and the frozen conventions (interface 2, 3);
* the canonical conic programs of the three admitted formulations (interface 9, 10, 12, 15.4):
  ``lower.p1d-vertex.v0``, the Upper control ``upper.p2c-vertex.v0`` and the P0 Upper candidate
  ``upper.p2d-jump-bernstein.v0``;
* von Mises and Drucker-Prager regions, fixed and reference body forces, linear facet tractions, bonded
  rigid footings with restrained components (interface 4 to 6);
* the ``ray``, ``fixed_components`` and ``support`` queries (interface 7);
* export of problems, programs and results in the frozen formats (interface 15);
* an optional open-solver path through Clarabel, imported lazily, so assembly needs only numpy and scipy;
  tests that solve are gated on Clarabel being installed. It also builds the optional interior fields of
  interface 15.3: a strictly admissible Lower field, and a relative-interior Upper field for every problem
  with ``alpha > 0`` (interface 17.4 step 4);
* the declared apex face ``Z`` of amendment 2 (interface 19.2.1): :func:`find_apex_zero` finds the Upper
  conditions forced onto their cone's apex by the cones together with the rows, the result carries them
  as the hashed array ``apex_zero``, and the Upper interior field is built on that face.

Every problem is checked with the contract's ``validate_problem`` before assembly and export, and refusals
carry the contract validator's error codes.

Every value this module produces is an **estimate**. Certification belongs to the independent checker
(interface 17); nothing here writes a certificate, a method report, a feasibility record (with its witness
keys), a refusal (with its budget) or a bracket. Amendment 2 asks of a producer only the result arrays
``apex_zero`` and, optionally, ``apex_zero_dual``; this module does not compute the dual.

Program identity. The exported matrix is canonical CSC: duplicate entries summed, explicit zeros removed,
row indices sorted within each column. Operators use the four vertex coordinates only; mid-edge nodes
sit at ``0.5 * (x_a + x_b)`` of their element's vertices (interface 3.1). Symmetry rows use the set's
normal ``n_S`` in both bounds (interface 5, amendment 1, item 3), evaluated by the contract's
``symmetry_set_normal``, the exact formula and order of amendment 2 (m6).

Result values (interface 15.3, amendment 1, item 2). ``objective`` is ``c^T x``. A Lower ``value`` is ``objective_sign * c^T x``: ``lambda``, or
``d^T R(sigma)`` for ``support``. An Upper ``value`` is the theorem's value of the stored field, ``(D - W0) /
W1`` or ``D - W0`` with the formulas of interface 4.3 and 12.2 (17.4 steps 1 and 6), the terms of a
declared face charged as ``max(0, pi)`` (19.2.1 item 8). Up to equality residuals it does not exceed the
objective, and differs from it by the solver's cone slack.

Command line::

    python tests/limit_reference.py <case> <mesh> <formulation> [<formulation> ...] --out DIR
                                    [--solve] [--interior] [--no-apex-face]
"""
import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import test_limit_contract as contract  # noqa: E402  (the frozen contract: parser, builder, hash, validators)

PRODUCER = {"name": "ladruno-limit-reference", "version": "0.1"}

LOWER_ID = "lower.p1d-vertex.v0"
CONTROL_ID = "upper.p2c-vertex.v0"
JUMP_ID = "upper.p2d-jump-bernstein.v0"
ADMITTED = {LOWER_ID: "lower", CONTROL_ID: "upper", JUMP_ID: "upper"}

# Open-solver budget. Not a certificate tolerance: Clarabel 0.11.1 runs at its defaults except for this
# wall-clock stop, so a run cannot hang a battery. A run ending on it gives solver_status time_limit and no
# number (interface 15.3). Uncalibrated: 237d owns solver settings.
CLARABEL_TIME_LIMIT_S = 900.0

S2, S3, S6 = math.sqrt(2.0), math.sqrt(3.0), math.sqrt(6.0)

# Interface 4.2.
Q = np.array([
    [1 / S2, -1 / S2, 0, 0, 0, 0],
    [-1 / S6, -1 / S6, 2 / S6, 0, 0, 0],
    [0, 0, 0, S2, 0, 0],
    [0, 0, 0, 0, S2, 0],
    [0, 0, 0, 0, 0, S2],
])
ONE = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])

# Interface 3.1 and 3.2.
EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
FACET_VERTICES = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))
FACET_NODES = ((1, 2, 3, 5, 9, 8), (0, 3, 2, 7, 9, 6), (0, 1, 3, 4, 8, 7), (0, 2, 1, 6, 5, 4))
FACET_EDGES = ((0, 1), (1, 2), (2, 0))

# Interface 10: int N dV over a straight tet, per unit volume.
BODY_WEIGHTS = np.array([-1.0 / 20.0] * 4 + [1.0 / 5.0] * 6)


def facet_work_table():
    """``int L_i N_j dA`` per unit area: traction vertex i (3), facet P2 node j (6). Interface 10."""
    table = np.zeros((3, 6))
    for i in range(3):
        for j in range(3):
            table[i, j] = 1.0 / 30.0 if i == j else -1.0 / 60.0
        for q, (a, b) in enumerate(FACET_EDGES):
            table[i, 3 + q] = 2.0 / 15.0 if i in (a, b) else 1.0 / 15.0
    return table


FACET_WORK = facet_work_table()


def bernstein_matrix():
    """Rows: Bernstein coefficient B; columns: nodal value k, facet-node order. Interface 12.1."""
    m = np.zeros((6, 6))
    for a in range(3):
        m[a, a] = 1.0
    for q, (a, b) in enumerate(FACET_EDGES):
        m[3 + q, 3 + q] = 2.0
        m[3 + q, a] = -0.5
        m[3 + q, b] = -0.5
    return m


BERNSTEIN = bernstein_matrix()

# Resultant weights w_ab of interface 6.
RESULTANT_W = np.full((3, 3), 1.0 / 12.0) + np.eye(3) * (1.0 / 6.0 - 1.0 / 12.0)


def fail(code, detail=""):
    contract.fail(code, detail)


# ---------------------------------------------------------------------------------------------------
# Voigt operators (interface 2.1)
# ---------------------------------------------------------------------------------------------------

def voigt_n(n):
    """``N(n)``, shape (..., 3, 6), for vectors n of shape (..., 3)."""
    n = np.asarray(n, dtype=float)
    out = np.zeros(n.shape[:-1] + (3, 6))
    out[..., 0, 0] = n[..., 0]
    out[..., 1, 1] = n[..., 1]
    out[..., 2, 2] = n[..., 2]
    out[..., 0, 3] = n[..., 1]
    out[..., 1, 3] = n[..., 0]
    out[..., 1, 4] = n[..., 2]
    out[..., 2, 4] = n[..., 1]
    out[..., 0, 5] = n[..., 2]
    out[..., 2, 5] = n[..., 0]
    return out


def skew(r):
    """``[r]_x`` with ``[r]_x w = r x w``, shape (..., 3, 3)."""
    r = np.asarray(r, dtype=float)
    out = np.zeros(r.shape[:-1] + (3, 3))
    out[..., 0, 1], out[..., 0, 2] = -r[..., 2], r[..., 1]
    out[..., 1, 0], out[..., 1, 2] = r[..., 2], -r[..., 0]
    out[..., 2, 0], out[..., 2, 1] = -r[..., 1], r[..., 0]
    return out


def tangent_basis(n):
    """Interface 5: ``m`` = first index of the smallest ``|n_m|``; ``T1 = e_m x n / |.|``, ``T2 = n x T1``."""
    n = np.asarray(n, dtype=float)
    m = int(np.argmin(np.abs(n)))
    e = np.zeros(3)
    e[m] = 1.0
    t1 = np.cross(e, n)
    t1 = t1 / np.linalg.norm(t1)
    return t1, np.cross(n, t1)


def tangent_bases(normals):
    """``[T1, T2]`` of :func:`tangent_basis` for many normals, shape (n, 2, 3)."""
    normals = np.asarray(normals, dtype=float).reshape(-1, 3)
    e = np.zeros_like(normals)
    e[np.arange(len(normals)), np.argmin(np.abs(normals), axis=1)] = 1.0
    t1 = np.cross(e, normals)
    t1 = t1 / np.linalg.norm(t1, axis=1)[:, None]
    return np.stack([t1, np.cross(normals, t1)], axis=1)


def strain_split(d):
    """``theta = tr D`` and the Frobenius norm of the deviator, engineering shear in ``d``. Vectorised."""
    d = np.asarray(d, dtype=float)
    theta = d[..., 0] + d[..., 1] + d[..., 2]
    dev = d[..., :3] - theta[..., None] / 3.0
    enorm = np.sqrt(np.sum(dev * dev, axis=-1) + 0.5 * np.sum(d[..., 3:] * d[..., 3:], axis=-1))
    return theta, enorm


def support_value(d, k, alpha):
    """Interface 4.3 value formula, without testing the domain (as 17.4 step 1)."""
    theta, enorm = strain_split(d)
    k, alpha = np.broadcast_arrays(np.asarray(k, float), np.asarray(alpha, float))
    out = S2 * k * enorm
    frictional = alpha > 0.0
    safe = np.where(frictional, alpha, 1.0)
    return np.where(frictional, k * theta / (3.0 * safe), out)


def jump_constants(alpha):
    """``tau``, ``k_J / k`` of interface 12.2. Refuses ``alpha >= 1/sqrt(12)``."""
    alpha = np.asarray(alpha, dtype=float)
    if np.any(12.0 * alpha * alpha >= 1.0):
        fail("formulation_not_applicable", "alpha >= 1/sqrt(12) admits only a zero jump")
    root = np.sqrt(1.0 - 12.0 * alpha * alpha)
    return 3.0 * alpha / root, 1.0 / root


# ---------------------------------------------------------------------------------------------------
# Geometry (interface 3)
# ---------------------------------------------------------------------------------------------------

class Geometry:
    """Element, facet and region data of one problem, from the four vertex coordinates."""

    def __init__(self, problem, arrays):
        nodes = np.asarray(arrays["nodes"], dtype=float)
        tets = np.asarray(arrays["tets"], dtype=np.int64)
        self.problem, self.arrays = problem, arrays
        self.nodes, self.tets = nodes, tets
        self.E, self.N = len(tets), len(nodes)
        xv = nodes[tets[:, :4]]
        jac = np.transpose(xv[:, 1:] - xv[:, :1], (0, 2, 1))
        det = np.linalg.det(jac)
        edge = np.stack([np.linalg.norm(xv[:, a] - xv[:, b], axis=1) for a, b in EDGES], axis=1)
        self.h = edge.max(axis=1)
        if np.any(det <= 0.0):
            fail("element_orientation", f"element {int(np.argmin(det))}")
        tol = contract.TOLERANCES
        if np.any(det / 6.0 < tol["TOL_GEOM_VOLUME"] * self.h ** 3):
            fail("element_degenerate", "volume below TOL_GEOM_VOLUME h^3")
        inv = np.linalg.inv(jac)
        grad = np.empty((self.E, 4, 3))
        grad[:, 1:] = inv
        grad[:, 0] = -inv.sum(axis=1)
        self.xv, self.grad, self.volume = xv, grad, det / 6.0
        xn = np.empty((self.E, 10, 3))
        xn[:, :4] = xv
        for q, (a, b) in enumerate(EDGES):
            xn[:, 4 + q] = 0.5 * (xv[:, a] + xv[:, b])
        self.xn = xn
        offset = np.linalg.norm(nodes[tets[:, 4:]] - xn[:, 4:], axis=2).max(axis=1)
        self.straightness = offset / self.h
        if np.any(self.straightness > tol["TOL_GEOM_STRAIGHT"]):
            fail("element_curved", f"element {int(np.argmax(self.straightness))}")

        regions = problem["regions"]
        region = np.asarray(arrays["element_region"], dtype=np.int64)
        for item in regions:
            if item["criterion"] not in ("von_mises", "drucker_prager"):
                fail("criterion_reserved" if item["criterion"] in ("tresca", "mohr_coulomb") else "criterion",
                     repr(item["criterion"]))
        self.k = np.array([float(item["k"]) for item in regions])[region]
        self.alpha = np.array([float(item.get("alpha", 0.0)) for item in regions])[region]
        self.b0 = np.array([item["b0"] for item in regions], dtype=float).reshape(-1, 3)[region]
        self.b1 = np.array([item["b1"] for item in regions], dtype=float).reshape(-1, 3)[region]
        for item in problem["footings"]:
            if item["interface_law"] != "bonded":
                fail("interface_law_reserved", repr(item["interface_law"]))

        # Facet geometry of every (e, f): outward unit normal and area.
        fv = np.asarray(FACET_VERTICES)
        xa, xb, xc = xv[:, fv[:, 0]], xv[:, fv[:, 1]], xv[:, fv[:, 2]]
        cross = np.cross(xb - xa, xc - xa)
        area = 0.5 * np.linalg.norm(cross, axis=2)
        self.facet_normal = cross / (2.0 * area[..., None])
        self.facet_area = area
        self._facets()
        self._boundary()

    def _facets(self):
        tets = self.tets
        keys = np.sort(tets[:, np.asarray(FACET_VERTICES)].reshape(-1, 3), axis=1)
        order = np.lexsort((keys[:, 2], keys[:, 1], keys[:, 0]))
        sk = keys[order]
        same = np.all(sk[1:] == sk[:-1], axis=1)
        if np.any(same[1:] & same[:-1]):
            fail("mesh_nonconforming", "facet with more than two owners")
        first = np.nonzero(same)[0]
        a, b = order[first], order[first + 1]
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        sort = np.argsort(lo, kind="stable")
        lo, hi = lo[sort], hi[sort]
        paired = np.zeros(len(keys), dtype=bool)
        paired[lo] = paired[hi] = True
        self.exterior_ids = np.nonzero(~paired)[0]
        self.e_minus, self.f_minus = lo // 4, lo % 4
        self.e_plus, self.f_plus = hi // 4, hi % 4
        self.F = len(lo)
        local = np.asarray(FACET_NODES)[self.f_minus]                       # (F, 6) in e-
        glob = tets[self.e_minus[:, None], local]
        match = tets[self.e_plus][:, None, :] == glob[:, :, None]           # (F, 6, 10)
        if not np.all(match.sum(axis=2) == 1):
            fail("mesh_nonconforming", "interior facet nodes do not match")
        self.minus_nodes = local
        self.plus_nodes = np.argmax(match, axis=2)
        self.interior_normal = self.facet_normal[self.e_minus, self.f_minus]
        self.interior_area = self.facet_area[self.e_minus, self.f_minus]
        longest = np.zeros(self.F)
        xf = self.xn[self.e_minus[:, None], local[:, :3]]
        for a, b in FACET_EDGES:
            longest = np.maximum(longest, np.linalg.norm(xf[:, a] - xf[:, b], axis=1))
        self.interior_longest_edge = longest

    def _boundary(self):
        problem, arrays = self.problem, self.arrays
        facets = np.asarray(arrays["facets"], dtype=np.int64).reshape(-1, 2)
        listed = facets[:, 0] * 4 + facets[:, 1]
        if len(listed) != len(self.exterior_ids) or np.any(listed != self.exterior_ids):
            fail("facets_mismatch", "facets is not the sorted exterior facet list")
        self.facets = facets
        self.facet_set = np.asarray(arrays["facet_set"], dtype=np.int64)
        self.sets = problem["boundary_sets"]
        self.kinds = np.array([item["kind"] for item in self.sets])
        for kind in self.kinds:
            if kind not in contract.BOUNDARY_KINDS:
                fail("boundary_kind_reserved" if kind == "interface" else "boundary_kind", repr(kind))
        self.boundary_normal = self.facet_normal[facets[:, 0], facets[:, 1]]
        self.boundary_area = self.facet_area[facets[:, 0], facets[:, 1]]
        # n_S of interface 5 by the exact formula of amendment 2 (m6): the contract's own evaluation order, so
        # independent assemblies write the same symmetry rows bit for bit, also on a tilted plane.
        self.set_normal = {}
        for index in np.unique(self.facet_set):
            index = int(index)
            if self.kinds[index] == "symmetry":
                self.set_normal[index] = contract.symmetry_set_normal(problem, arrays, index)

    # Strain rate operator (interface 3.3, 10) ------------------------------------------------------

    def p2_vertex_gradients(self):
        """``grad N_a(x_m)``, shape (E, 4 vertices m, 10 nodes a, 3)."""
        g = self.grad
        out = np.zeros((self.E, 4, 10, 3))
        for m in range(4):
            for i in range(4):
                out[:, m, i] = (3.0 if i == m else -1.0) * g[:, i]
            for q, (a, b) in enumerate(EDGES):
                if a == m:
                    out[:, m, 4 + q] = 4.0 * g[:, b]
                elif b == m:
                    out[:, m, 4 + q] = 4.0 * g[:, a]
        return out

    def strain_operator(self):
        """``B_ei``, shape (E, 4, 6, 30), on the element's nodal velocity ordered (node, component)."""
        nt = np.swapaxes(voigt_n(self.p2_vertex_gradients()), -1, -2)       # (E, 4, 10, 6, 3)
        return np.transpose(nt, (0, 1, 3, 2, 4)).reshape(self.E, 4, 6, 30)

    def strain_rates(self, element_velocity):
        """``d_ei`` of an element velocity field (E, 10, 3)."""
        v = np.asarray(element_velocity, dtype=float).reshape(self.E, 30)
        return np.einsum("eiab,eb->eia", self.strain_operator(), v)

    # Boundary node sets (interface 15.4) -----------------------------------------------------------

    def set_element_nodes(self, index):
        """Sorted unique (e, a) pairs of the six P2 nodes of every facet in boundary set ``index``."""
        rows = np.nonzero(self.facet_set == index)[0]
        e = np.repeat(self.facets[rows, 0], 6)
        a = np.asarray(FACET_NODES)[self.facets[rows, 1]].reshape(-1)
        pairs = np.unique(np.stack([e, a], axis=1), axis=0) if len(rows) else np.zeros((0, 2), np.int64)
        return pairs

    def set_global_nodes(self, index):
        pairs = self.set_element_nodes(index)
        return np.unique(self.tets[pairs[:, 0], pairs[:, 1]]) if len(pairs) else np.zeros(0, np.int64)

    def global_node_coordinates(self):
        """Affine coordinates of every global node (mid-edge nodes from their element's vertices)."""
        out = np.zeros((self.N, 3))
        out[self.tets.reshape(-1)] = self.xn.reshape(-1, 3)
        return out


# ---------------------------------------------------------------------------------------------------
# Operators shared by the bounds
# ---------------------------------------------------------------------------------------------------

def footing_resultant_matrix(geo):
    """``R = Rmat @ sigma``, shape (6 n_f, 24 E), sigma ordered (e, i, component). Interface 6."""
    n_f = len(geo.problem["footings"])
    rows, cols, vals = [], [], []
    for row, (e, f) in enumerate(geo.facets):
        item = geo.sets[int(geo.facet_set[row])]
        if item["kind"] != "footing":
            continue
        fi = int(item["footing"])
        x_ref = np.asarray(geo.problem["footings"][fi]["x_ref"], dtype=float)
        verts = FACET_VERTICES[int(f)]
        area = geo.boundary_area[row]
        nmat = voigt_n(geo.boundary_normal[row])                            # (3, 6)
        r = geo.xv[e, list(verts)] - x_ref                                  # (3, 3)
        for jb, b in enumerate(verts):
            moment = area * sum(RESULTANT_W[ja, jb] * skew(r[ja]) for ja in range(3))
            block = np.vstack([area / 3.0 * nmat, moment @ nmat])           # (6, 6)
            rr, cc = np.meshgrid(6 * fi + np.arange(6), 24 * e + 6 * b + np.arange(6), indexing="ij")
            rows.append(rr.ravel())
            cols.append(cc.ravel())
            vals.append(block.ravel())
    if not rows:
        return sp.csr_matrix((6 * n_f, 24 * geo.E))
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                         shape=(6 * n_f, 24 * geo.E))


def element_node_work(geo, which):
    """Coefficients of ``sum_e b_e . int v dV + sum_facets int t . v dA`` per element node, (E, 10, 3)."""
    b = geo.b0 if which == 0 else geo.b1
    out = (geo.volume[:, None] * BODY_WEIGHTS[None, :])[:, :, None] * b[:, None, :]
    t = np.asarray(geo.arrays["t0" if which == 0 else "t1"], dtype=float)
    for row, (e, f) in enumerate(geo.facets):
        if not np.any(t[row]):
            continue
        out[e, list(FACET_NODES[int(f)])] += geo.boundary_area[row] * (FACET_WORK.T @ t[row])
    return out


def query_tables(problem, query):
    """``(form, F0, F1, direction, restrained mask)`` with tables of shape (n_f, 6). Interface 7."""
    n_f = len(problem["footings"])
    form = query["form"]
    mask = np.zeros((n_f, 6), dtype=bool)
    for f, c in query["restrained"]:
        mask[f, c] = True
    zero = np.zeros((n_f, 6))
    if form == "ray":
        return form, np.asarray(query["F0"], float).reshape(n_f, 6), np.asarray(query["F1"], float).reshape(n_f, 6), \
            None, mask
    if form == "fixed_components":
        f0 = zero.copy()
        for f, c, value in query["fixed"]:
            f0[f, c] = value
        return form, f0, np.asarray(query["direction"], float).reshape(n_f, 6), None, mask
    if form == "support":
        return form, zero, zero.copy(), np.asarray(query["direction"], float).reshape(n_f, 6), mask
    fail("query_form", repr(form))


class _Triplets:
    def __init__(self):
        self.rows, self.cols, self.vals = [], [], []

    def add(self, rows, cols, vals):
        rows, cols, vals = np.broadcast_arrays(np.asarray(rows), np.asarray(cols), np.asarray(vals, dtype=float))
        self.rows.append(rows.ravel())
        self.cols.append(cols.ravel())
        self.vals.append(vals.ravel())

    def matrix(self, n_rows, n_cols):
        if not self.rows:
            return sp.csc_matrix((n_rows, n_cols))
        a = sp.coo_matrix((np.concatenate(self.vals), (np.concatenate(self.rows), np.concatenate(self.cols))),
                          shape=(n_rows, n_cols)).tocsc()
        a.sum_duplicates()
        a.eliminate_zeros()
        a.sort_indices()
        return a


class Program:
    """A canonical program ``min c^T x, A x + s = b, s in K`` (interface 15.4) and its layout."""

    def __init__(self, problem, query, formulation, bound, matrix, b, c, n_zero, cones, layout):
        self.problem, self.query, self.formulation, self.bound = problem, query, formulation, bound
        self.A, self.b, self.c = matrix, b, c
        self.n_zero, self.cones, self.layout = n_zero, cones, layout
        self.objective_sign = -1 if bound == "lower" else 1

    def document(self, arrays_name):
        return {
            "format": contract.FORMAT_PROGRAM, "format_version": contract.FORMAT_VERSION,
            "arrays": arrays_name, "problem_sha256": self.problem["problem_sha256"],
            "query_id": self.query["id"], "bound": self.bound, "formulation": self.formulation,
            "objective_sign": self.objective_sign, "n_variables": int(self.A.shape[1]),
            "n_rows": int(self.A.shape[0]), "n_zero": int(self.n_zero), "cones": [dict(c) for c in self.cones],
        }

    def arrays(self):
        return {
            "A_indptr": self.A.indptr.astype("<i8"), "A_indices": self.A.indices.astype("<i8"),
            "A_data": self.A.data.astype("<f8"), "b": np.asarray(self.b, dtype="<f8"),
            "c": np.asarray(self.c, dtype="<f8"),
        }

    def value(self, x):
        return float(self.objective_sign * (self.c @ x))

    def zero_residual(self, x):
        """``A x - b`` on the zero cone rows."""
        return (self.A @ x - self.b)[: self.n_zero]

    def cone_slacks(self, x):
        """``s = b - A x`` split into the cones, as a list of (count, dim) arrays."""
        s = (self.b - self.A @ x)[self.n_zero:]
        out, start = [], 0
        for cone in self.cones:
            size = cone["dim"] * cone["count"]
            out.append(s[start:start + size].reshape(cone["count"], cone["dim"]))
            start += size
        return out


def query_by_id(problem, query_id):
    return contract.query_by_id(problem, query_id)


def check_formulation(formulation, problem):
    """The contract's refusal codes: ``formulation_not_admitted``, ``formulation_unknown``, ``..._not_applicable``."""
    bound = ADMITTED.get(formulation)
    if bound is None:
        prefix = formulation.split(".", 1)[0] if isinstance(formulation, str) else None
        bound = prefix if prefix in ("lower", "upper") else "lower"
    contract.validate_formulation(bound, formulation, problem)


def check_query(problem, query):
    """Refuse query content that the programs would otherwise ignore (beyond ``validate_problem``)."""
    if query["form"] != "fixed_components":
        return
    restrained = {tuple(item) for item in query["restrained"]}
    seen = set()
    for f, c, _ in query["fixed"]:
        if (f, c) in seen:
            fail("query_fixed", f"duplicate fixed component {[f, c]}")
        seen.add((f, c))
        if (f, c) in restrained:
            fail("query_restrained_load", f"fixed component {[f, c]} is restrained and has no load row")


def validate(problem, arrays):
    """The contract's problem validator, then the query rules of :func:`check_query`."""
    contract.validate_problem(problem, arrays)
    for query in problem["queries"]:
        check_query(problem, query)


def assemble(problem, arrays, formulation, query_id="q", geometry=None, validated=False):
    """The canonical program of ``formulation`` for one query (interface 15.4).

    The problem is validated first unless the caller has just done so (``validated=True``).
    """
    if not validated:
        validate(problem, arrays)
    check_formulation(formulation, problem)
    geo = geometry or Geometry(problem, arrays)
    query = query_by_id(problem, query_id)
    if formulation == LOWER_ID:
        return _assemble_lower(geo, query)
    return _assemble_upper(geo, query, broken=formulation == JUMP_ID)


# ---------------------------------------------------------------------------------------------------
# Lower, lower.p1d-vertex.v0 (interface 9, 15.4)
# ---------------------------------------------------------------------------------------------------

def _assemble_lower(geo, query):
    problem = geo.problem
    form, f0, f1, direction, mask = query_tables(problem, query)
    ray = form != "support"
    n_f = len(problem["footings"])
    off = 1 if ray else 0
    E = geo.E
    n_var = off + 24 * E

    def sig(e, i):
        return off + 24 * np.asarray(e)[..., None] + 6 * np.asarray(i)[..., None] + np.arange(6)

    t = _Triplets()
    rhs = []
    row = 0
    # 1. Equilibrium: [b1_e | sum_i N(grad L_i)] = -b0_e.
    ndiv = voigt_n(geo.grad)                                                 # (E, 4, 3, 6)
    e_idx = np.arange(E)
    for i in range(4):
        rr = row + 3 * e_idx[:, None, None] + np.arange(3)[None, :, None]
        t.add(rr, sig(e_idx, np.full(E, i))[:, None, :], ndiv[:, i])
    if ray:
        t.add(row + 3 * e_idx[:, None] + np.arange(3), 0, geo.b1)
    rhs.append(-geo.b0.reshape(-1))
    row += 3 * E
    # 2. Interior facets: N(n) sigma_- - N(n) sigma_+ = 0 at the three vertices of e-'s facet.
    nmat = voigt_n(geo.interior_normal)                                      # (F, 3, 6)
    f_idx = np.arange(geo.F)
    for k in range(3):
        rr = row + 9 * f_idx[:, None, None] + 3 * k + np.arange(3)[None, :, None]
        t.add(rr, sig(geo.e_minus, geo.minus_nodes[:, k])[:, None, :], nmat)
        t.add(rr, sig(geo.e_plus, geo.plus_nodes[:, k])[:, None, :], -nmat)
    rhs.append(np.zeros(9 * geo.F))
    row += 9 * geo.F
    # 3. Exterior facets in facets order.
    t0 = np.asarray(geo.arrays["t0"], dtype=float)
    t1 = np.asarray(geo.arrays["t1"], dtype=float)
    for r, (e, f) in enumerate(geo.facets):
        index = int(geo.facet_set[r])
        kind = geo.kinds[index]
        verts = FACET_VERTICES[int(f)]
        if kind in ("free", "traction"):
            n = voigt_n(geo.boundary_normal[r])
            for k, a in enumerate(verts):
                t.add(row + np.arange(3)[:, None], sig(e, a)[None, :], n)
                if ray:
                    t.add(row + np.arange(3), 0, -t1[r, k])
                rhs.append(t0[r, k])
                row += 3
        elif kind == "symmetry":
            normal = geo.set_normal[index]
            t1v, t2v = tangent_basis(normal)
            tn = np.vstack([t1v, t2v]) @ voigt_n(normal)                     # (2, 6)
            for a in verts:
                t.add(row + np.arange(2)[:, None], sig(e, a)[None, :], tn)
                rhs.append(np.zeros(2))
                row += 2
    # 4. Ray forms: R - lambda F1 = F0 per footing and unrestrained component.
    rmat = footing_resultant_matrix(geo).tocsr()
    if ray:
        for fi in range(n_f):
            for c in range(6):
                if mask[fi, c]:
                    continue
                line = rmat.getrow(6 * fi + c).tocoo()
                t.add(row, off + line.col, line.data)
                t.add(row, 0, -f1[fi, c])
                rhs.append([f0[fi, c]])
                row += 1
    n_zero = row
    # 5. Cones per (e, i): [sqrt2 alpha one^T ; -Q] sigma_ei + s = (sqrt2 k, 0).
    cone_block = np.zeros((E, 6, 6))
    cone_block[:, 0, :] = S2 * geo.alpha[:, None] * ONE[None, :]
    cone_block[:, 1:, :] = -Q[None]
    for i in range(4):
        rr = row + 24 * e_idx[:, None, None] + 6 * i + np.arange(6)[None, :, None]
        t.add(rr, sig(e_idx, np.full(E, i))[:, None, :], cone_block)
    cone_rhs = np.zeros((E, 4, 6))
    cone_rhs[:, :, 0] = S2 * geo.k[:, None]
    rhs.append(cone_rhs.reshape(-1))
    row += 24 * E
    b = np.concatenate([np.atleast_1d(np.asarray(item, dtype=float)).reshape(-1) for item in rhs])
    assert len(b) == row
    c = np.zeros(n_var)
    if ray:
        c[0] = -1.0
    else:
        c[off:] = -(direction.reshape(-1) @ rmat)
    layout = {"lambda": 0 if ray else None, "sigma": off, "resultant": rmat, "geometry": geo, "form": form,
              "direction": direction}
    return Program(problem, query, LOWER_ID, "lower", t.matrix(row, n_var), b, c, n_zero,
                   [{"dim": 6, "count": 4 * E}], layout)


# ---------------------------------------------------------------------------------------------------
# Upper: control upper.p2c-vertex.v0 and jump upper.p2d-jump-bernstein.v0 (interface 10, 12, 15.4)
# ---------------------------------------------------------------------------------------------------

def _assemble_upper(geo, query, broken):
    problem = geo.problem
    form, f0, f1, direction, mask = query_tables(problem, query)
    ray = form != "support"
    n_f = len(problem["footings"])
    E, F = geo.E, geo.F
    n_v = 30 * E if broken else 3 * geo.N
    v_off, V_off = 0, n_v
    rz_off = V_off + 6 * n_f
    jump_off = rz_off + 24 * E
    n_var = jump_off + (18 * F if broken else 0)
    tets = geo.tets
    e_idx = np.arange(E)

    def vel(e, a):
        """Column indices of the three velocity components at element node (e, a)."""
        e, a = np.asarray(e), np.asarray(a)
        base = 30 * e + 3 * a if broken else 3 * tets[e, a]
        return np.asarray(base)[..., None] + np.arange(3)

    def rz(e, i):
        return rz_off + 24 * np.asarray(e)[..., None] + 6 * np.asarray(i)[..., None] + np.arange(6)

    t = _Triplets()
    rhs = []
    row = 0
    # 1. Flow: B_ei v - sqrt2 alpha r one - Q^T z = 0.
    bop = geo.strain_operator()                                              # (E, 4, 6, 30)
    elem_cols = vel(e_idx[:, None], np.arange(10)[None, :]).reshape(E, 30)
    flow_rz = np.zeros((E, 6, 6))
    flow_rz[:, :, 0] = -S2 * geo.alpha[:, None] * ONE[None, :]
    flow_rz[:, :, 1:] = -Q.T[None]
    for i in range(4):
        rr = row + 24 * e_idx[:, None, None] + 6 * i + np.arange(6)[None, :, None]
        t.add(rr, elem_cols[:, None, :], bop[:, i])
        t.add(rr, rz(e_idx, np.full(E, i))[:, None, :], flow_rz)
    rhs.append(np.zeros(24 * E))
    row += 24 * E
    # 2. Jumps (jump formulation): n.c_B - tau rho = 0, T1.c_B - y1 = 0, T2.c_B - y2 = 0.
    if broken and F:
        tau, kj = jump_constants(geo.alpha[geo.e_minus])
        basis = np.concatenate([geo.interior_normal[:, None, :], tangent_bases(geo.interior_normal)], axis=1)
        cm = vel(geo.e_minus[:, None], geo.minus_nodes)                       # (F, 6 nodes, 3)
        cp = vel(geo.e_plus[:, None], geo.plus_nodes)
        f_idx = np.arange(F)
        # rows (F, B, q); columns (node k, component c); coefficient BERNSTEIN[B, k] * basis[F, q, c]
        rr = row + 18 * f_idx[:, None, None, None, None] + 3 * np.arange(6)[None, :, None, None, None] \
            + np.arange(3)[None, None, :, None, None]
        coeff = BERNSTEIN[None, :, None, :, None] * basis[:, None, :, None, :]  # (F, 6, 3, 6, 3)
        t.add(rr, cp[:, None, None, :, :], coeff)
        t.add(rr, cm[:, None, None, :, :], -coeff)
        slack = np.empty((F, 6, 3))
        slack[:, :, 0] = -tau[:, None]
        slack[:, :, 1:] = -1.0
        t.add(row + np.arange(18 * F), jump_off + np.arange(18 * F), slack.reshape(-1))
        row += 18 * F
        rhs.append(np.zeros(18 * F))
    # 3.-5. Boundary rows.
    kinds = geo.kinds
    fixed_sets = [i for i, kind in enumerate(kinds) if kind == "fixed"]
    if broken:
        def nodes_of(index):
            return [tuple(p) for p in geo.set_element_nodes(index)]
    else:
        def nodes_of(index):
            return [int(a) for a in geo.set_global_nodes(index)]
    fixed = sorted(set().union(*[set(nodes_of(i)) for i in fixed_sets])) if fixed_sets else []
    fixed_set = set(fixed)

    def cols_of(node):
        if broken:
            return vel(node[0], node[1])
        return 3 * node + np.arange(3)

    def coords_of(node):
        if broken:
            return geo.xn[node[0], node[1]]
        return xg[node]

    xg = None if broken else geo.global_node_coordinates()
    for node in fixed:
        t.add(row + np.arange(3), cols_of(node), 1.0)
        rhs.append(np.zeros(3))
        row += 3
    for index, kind in enumerate(kinds):
        if kind != "symmetry":
            continue
        normal = geo.set_normal[index]
        for node in nodes_of(index):
            if node in fixed_set:
                continue
            t.add(row, cols_of(node), normal)
            rhs.append([0.0])
            row += 1
    for fi in range(n_f):
        members = set()
        for index, item in enumerate(geo.sets):
            if item["kind"] == "footing" and int(item["footing"]) == fi:
                members |= set(nodes_of(index))
        x_ref = np.asarray(problem["footings"][fi]["x_ref"], dtype=float)
        for node in sorted(members):
            r = coords_of(node) - x_ref
            rr = row + np.arange(3)
            t.add(rr, cols_of(node), 1.0)
            t.add(rr, V_off + 6 * fi + np.arange(3), -1.0)
            # -(w x r) = r x w = [r]_x w
            t.add(rr[:, None], (V_off + 6 * fi + 3 + np.arange(3))[None, :], skew(r))
            rhs.append(np.zeros(3))
            row += 3
    # 6. Restrained components, then the query normalisation.
    for fi, comp in query["restrained"]:
        t.add(row, V_off + 6 * fi + comp, 1.0)
        rhs.append([0.0])
        row += 1
    w0, w1 = work_vectors(geo, f0, f1, broken, include_f0=ray)
    work_row = row if ray else None
    if ray:
        nz = np.nonzero(w1)[0]
        t.add(row, nz, w1[nz])
        rhs.append([1.0])
        row += 1
    else:
        for fi in range(n_f):
            for comp in range(6):
                if mask[fi, comp]:
                    continue
                t.add(row, V_off + 6 * fi + comp, 1.0)
                rhs.append([direction[fi, comp]])
                row += 1
    n_zero = row
    # 7. Cones: -(r, z) + s = 0, then -(rho, y1, y2) + s = 0.
    t.add(row + np.arange(24 * E), rz_off + np.arange(24 * E), -1.0)
    rhs.append(np.zeros(24 * E))
    row += 24 * E
    cones = [{"dim": 6, "count": 4 * E}]
    if broken:
        t.add(row + np.arange(18 * F), jump_off + np.arange(18 * F), -1.0)
        rhs.append(np.zeros(18 * F))
        row += 18 * F
        cones.append({"dim": 3, "count": 6 * F})
    b = np.concatenate([np.atleast_1d(np.asarray(item, dtype=float)).reshape(-1) for item in rhs]) if rhs \
        else np.zeros(0)
    assert len(b) == row
    c = np.zeros(n_var)
    c[rz_off + 24 * e_idx[:, None] + 6 * np.arange(4)[None, :]] = (geo.volume / 4.0 * S2 * geo.k)[:, None]
    if broken and F:
        tau, kj = jump_constants(geo.alpha[geo.e_minus])
        c[jump_off + 18 * np.arange(F)[:, None] + 3 * np.arange(6)[None, :]] = \
            (geo.interior_area / 6.0 * geo.k[geo.e_minus] * kj)[:, None]
    c[: rz_off] -= w0
    layout = {"broken": broken, "v": v_off, "V": V_off, "rz": rz_off, "jump": jump_off if broken else None,
              "geometry": geo, "w0": w0, "w1": w1, "form": form, "direction": direction, "work_row": work_row}
    formulation = JUMP_ID if broken else CONTROL_ID
    return Program(problem, query, formulation, "upper", t.matrix(row, n_var), b, c, n_zero, cones, layout)


def work_vectors(geo, f0, f1, broken, include_f0=True):
    """``W0`` and ``W1`` as coefficient vectors over ``(v, V)`` (interface 10)."""
    n_f = f0.shape[0]
    n_v = 30 * geo.E if broken else 3 * geo.N
    out = []
    for which, table in ((0, f0), (1, f1)):
        per_node = element_node_work(geo, which)
        vec = np.zeros(n_v + 6 * n_f)
        if broken:
            vec[:n_v] = per_node.reshape(-1)
        else:
            np.add.at(vec, (3 * geo.tets[:, :, None] + np.arange(3)).reshape(-1), per_node.reshape(-1))
        if which == 1 or include_f0:
            vec[n_v:] = table.reshape(-1)
        out.append(vec)
    return out[0], out[1]


# ---------------------------------------------------------------------------------------------------
# Field values, independent of the program (producer-side check of the theorem's value)
# ---------------------------------------------------------------------------------------------------

def upper_field_value(problem, arrays, query_id, velocity, footing_velocity, geometry=None, apex_zero=None):
    """``(U, D, W0, W1)`` of an Upper field by the formulas of interface 4.3, 12.2 and 17.4 step 6.

    ``velocity`` is continuous ``(N, 3)`` or broken ``(E, 10, 3)``. The domain is not tested. The terms of
    the conditions of a declared apex face ``apex_zero`` (interface 19.2.1 item 8) are charged as
    ``max(0, pi)``: their argument is zero up to rounding, and a negative charge is never credited.
    """
    geo = geometry or Geometry(problem, arrays)
    query = query_by_id(problem, query_id)
    form, f0, f1, direction, mask = query_tables(problem, query)
    velocity = np.asarray(velocity, dtype=float)
    field = velocity if velocity.ndim == 3 else velocity[geo.tets]
    flow_z, jump_z = apex_zero_masks(geo, apex_zero)
    d = geo.strain_rates(field)
    pi = support_value(d, geo.k[:, None], geo.alpha[:, None])
    pi = np.where(flow_z, np.maximum(pi, 0.0), pi)
    dissipation = float(np.sum(geo.volume[:, None] / 4.0 * pi))
    if velocity.ndim == 3 and geo.F:
        coeff = jump_coefficients(geo, field)
        n = geo.interior_normal[:, None, :]
        opening = np.sum(coeff * n, axis=2)
        slip = np.linalg.norm(coeff - opening[..., None] * n, axis=2)
        k, alpha = geo.k[geo.e_minus][:, None], geo.alpha[geo.e_minus][:, None]
        density = np.where(alpha > 0.0, k * opening / (3.0 * np.where(alpha > 0, alpha, 1.0)), k * slip)
        density = np.where(jump_z, np.maximum(density, 0.0), density)
        dissipation += float(np.sum(geo.interior_area[:, None] / 6.0 * density))
    ray = form != "support"
    w = []
    for which, table in ((0, f0), (1, f1)):
        per_node = element_node_work(geo, which)
        value = float(np.sum(per_node * field))
        if which == 1 or ray:
            value += float(np.sum(table * np.asarray(footing_velocity, float).reshape(-1, 6)))
        w.append(value)
    w0, w1 = w
    upper = (dissipation - w0) / w1 if ray else dissipation - w0
    return upper, dissipation, w0, w1


def jump_coefficients(geo, element_velocity):
    """Bernstein coefficients ``c_FB`` of ``v+ - v-`` on every interior facet, shape (F, 6, 3)."""
    v = np.asarray(element_velocity, dtype=float)
    plus = v[geo.e_plus[:, None], geo.plus_nodes]
    minus = v[geo.e_minus[:, None], geo.minus_nodes]
    return np.einsum("bk,fkc->fbc", BERNSTEIN, plus - minus)


def lower_resultant(problem, arrays, sigma, geometry=None):
    geo = geometry or Geometry(problem, arrays)
    return (footing_resultant_matrix(geo) @ np.asarray(sigma, dtype=float).reshape(-1)).reshape(-1, 6)


def yield_margin(problem, arrays, sigma, geometry=None):
    """``sqrt(J2) + alpha I1 - k`` at every element vertex, by invariants. Shape (E, 4)."""
    geo = geometry or Geometry(problem, arrays)
    s = np.asarray(sigma, dtype=float)
    i1 = s[..., 0] + s[..., 1] + s[..., 2]
    dev = s[..., :3] - i1[..., None] / 3.0
    j2 = 0.5 * np.sum(dev * dev, axis=-1) + np.sum(s[..., 3:] ** 2, axis=-1)
    return np.sqrt(j2) + geo.alpha[:, None] * i1 - geo.k[:, None]


# ---------------------------------------------------------------------------------------------------
# Files (interface 15)
# ---------------------------------------------------------------------------------------------------

def problem_stem(problem):
    name = problem.get("arrays") or ""
    suffix = ".problem.npz"
    return name[: -len(suffix)] if name.endswith(suffix) else problem["name"].replace(":", "_")


def record_stem(problem, query_id, formulation):
    return f"{problem_stem(problem)}.{query_id}.{formulation}"


def _json_ready(obj):
    if isinstance(obj, dict):
        return {str(key): _json_ready(value) for key, value in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_ready(item) for item in obj]
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, float) and not math.isfinite(obj):
        fail("json_non_finite", repr(obj))
    return obj


def write_document(directory, stem, doc, arrays=None):
    """``<stem>.json`` and, when given, the NPZ named by ``doc['arrays']``, beside it."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    doc = _json_ready(doc)
    path = directory / f"{stem}.json"
    path.write_text(json.dumps(doc, indent=1, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    if arrays:
        np.savez(directory / doc["arrays"], **arrays)
    return path


def export_problem(directory, problem, arrays):
    return write_document(directory, f"{problem_stem(problem)}.problem", problem, arrays)


def export_program(directory, program):
    """Validate the problem and the program document against the contract, then write both files."""
    geo = program.layout["geometry"]
    validate(program.problem, geo.arrays)
    stem = record_stem(program.problem, program.query["id"], program.formulation)
    doc, arrays = program.document(f"{stem}.program.npz"), program.arrays()
    contract.validate_program(doc, arrays, program.problem, geo.arrays)
    return write_document(directory, f"{stem}.program", doc, arrays)


def result_document(program, solver_status, x=None, solver=None, wall_time_s=None, interior=None, apex_zero=None,
                    apex_zero_dual=None):
    """The result document and arrays of a solved program (interface 15.3). Values are estimates.

    A ray-form Upper field whose reference work ``W1`` is not positive gives no number: it is recorded as
    ``numerical_error`` with null ``value`` and ``objective`` and no arrays (interface 13, ``W1 > 0``).

    ``apex_zero`` is the declared apex face ``Z`` of interface 19.2.1 (from :func:`find_apex_zero`), written as
    the result array ``apex_zero`` before hashing, so ``result_sha256`` covers it; ``None`` or an empty face
    writes no array (an empty array is refused). ``apex_zero_dual`` is the optional ``(n_Z, 6)`` dual of item 11.
    A face is refused on a Lower (``apex_zero_lower``) and dropped from a result without a field (item 1).
    """
    problem, query = program.problem, program.query
    stem = record_stem(problem, query["id"], program.formulation)
    has_value = solver_status in contract.SOLVER_STATUS_WITH_VALUE and x is not None
    face = None if apex_zero is None or len(apex_zero) == 0 else apex_zero_array(apex_zero)
    if face is None and apex_zero_dual is not None:
        fail("apex_zero_invalid", "a dual without a declared face")
    if face is not None and program.bound == "lower":
        fail("apex_zero_lower", "a Lower result declares no apex face")
    value, arrays = None, {}
    if has_value:
        x = np.asarray(x, dtype=float)
        arrays = extract_fields(program, x)
        if program.bound == "lower":
            value = program.value(x)
        else:
            # Amendment 1, item 2: the value of the stored field (17.4 step 6), not c^T x. An interior-point
            # solution keeps r slightly above |z|, which would put the objective a few 1e-6 above the
            # field's own value and fail TOL_VALUE_MATCH for no defect.
            velocity = arrays["element_velocity"] if program.layout["broken"] else arrays["velocity"]
            geo = program.layout["geometry"]
            value, _, _, w1 = upper_field_value(problem, geo.arrays, query["id"], velocity, arrays["footing_velocity"],
                                                geo, apex_zero=face)
            if query["form"] != "support" and not w1 > 0.0:
                solver_status, has_value, value = "numerical_error", False, None
    arrays = arrays if has_value else {}
    doc = {
        "format": contract.FORMAT_RESULT, "format_version": contract.FORMAT_VERSION,
        "arrays": f"{stem}.result.npz" if has_value else None,
        "problem_sha256": problem["problem_sha256"], "query_id": query["id"], "bound": program.bound,
        "formulation": program.formulation, "producer": dict(PRODUCER),
        "solver": solver or {"name": "none", "version": "0", "settings": {}},
        "solver_status": solver_status, "objective": float(program.c @ x) if has_value else None, "value": value,
        "wall_time_s": None if wall_time_s is None else float(wall_time_s),
    }
    if has_value:
        if interior is not None:
            if program.bound == "lower":
                arrays["interior_sigma"] = np.asarray(interior[0], dtype="<f8")
                doc["interior_value"] = float(interior[1])
            else:
                name = "interior_element_velocity" if program.layout["broken"] else "interior_velocity"
                arrays[name] = np.asarray(interior[0], dtype="<f8")
                arrays["interior_footing_velocity"] = np.asarray(interior[1], dtype="<f8")
        if face is not None:
            arrays["apex_zero"] = face
            if apex_zero_dual is not None:
                dual = np.asarray(apex_zero_dual, dtype="<f8")
                if dual.shape != (len(face), 6):
                    fail("apex_zero_invalid", f"dual shape {dual.shape}")
                arrays["apex_zero_dual"] = dual
            contract.validate_apex_zero(doc, arrays, problem, program.layout["geometry"].arrays)
    doc["result_sha256"] = contract.document_hash(doc, arrays, "result_sha256")
    return doc, arrays


def extract_fields(program, x):
    """Required result arrays (interface 15.3) from a program vector."""
    layout = program.layout
    geo = layout["geometry"]
    n_f = len(program.problem["footings"])
    if program.bound == "lower":
        sigma = x[layout["sigma"]: layout["sigma"] + 24 * geo.E].reshape(geo.E, 4, 6)
        resultant = (layout["resultant"] @ sigma.reshape(-1)).reshape(n_f, 6)
        return {"sigma": sigma.astype("<f8"), "resultant": resultant.astype("<f8")}
    footing = x[layout["V"]: layout["V"] + 6 * n_f].reshape(n_f, 6).astype("<f8")
    if layout["broken"]:
        return {"element_velocity": x[: 30 * geo.E].reshape(geo.E, 10, 3).astype("<f8"), "footing_velocity": footing}
    return {"velocity": x[: 3 * geo.N].reshape(geo.N, 3).astype("<f8"), "footing_velocity": footing}


def export_result(directory, program, *args, **kwargs):
    """Write the result; a no-value status also removes a stale ``.result.npz`` of an earlier run."""
    doc, arrays = result_document(program, *args, **kwargs)
    stem = record_stem(program.problem, program.query["id"], program.formulation)
    stale = Path(directory) / f"{stem}.result.npz"
    if not arrays and stale.exists():
        stale.unlink()
    write_document(directory, f"{stem}.result", doc, arrays)
    return doc, arrays


# ---------------------------------------------------------------------------------------------------
# Optional open solver (Zone B only): Clarabel, imported lazily
# ---------------------------------------------------------------------------------------------------

CLARABEL_STATUS = {
    "Solved": "solved", "AlmostSolved": "almost_solved",
    "PrimalInfeasible": "primal_infeasible", "DualInfeasible": "dual_infeasible",
    "MaxIterations": "iteration_limit", "MaxTime": "time_limit",
    "NumericalError": "numerical_error", "InsufficientProgress": "numerical_error",
    # An almost-certificate of infeasibility is not recorded as infeasibility.
    "AlmostPrimalInfeasible": "numerical_error", "AlmostDualInfeasible": "numerical_error",
}


def _clarabel_cones(clarabel, n_zero, cones, extra_nonneg=0):
    out = [clarabel.ZeroConeT(int(n_zero))] if n_zero else []
    for cone in cones:
        out.extend(clarabel.SecondOrderConeT(int(cone["dim"])) for _ in range(int(cone["count"])))
    if extra_nonneg:
        out.append(clarabel.NonnegativeConeT(int(extra_nonneg)))
    return out


def _clarabel_solve(matrix, b, c, cones, time_limit):
    import clarabel
    n = matrix.shape[1]
    settings = clarabel.DefaultSettings()
    settings.verbose = False
    settings.time_limit = float(time_limit)
    solver = clarabel.DefaultSolver(sp.csc_matrix((n, n)), np.asarray(c, float), sp.csc_matrix(matrix),
                                    np.asarray(b, float), cones, settings)
    start = time.perf_counter()
    solution = solver.solve()
    wall = time.perf_counter() - start
    status = CLARABEL_STATUS.get(str(solution.status).split(".")[-1], "numerical_error")
    record = {
        "name": "clarabel", "version": str(clarabel.__version__),
        "settings": {"defaults": True, "time_limit": float(time_limit), "verbose": False,
                     "tol_gap_abs": float(settings.tol_gap_abs), "tol_gap_rel": float(settings.tol_gap_rel),
                     "tol_feas": float(settings.tol_feas), "max_iter": int(settings.max_iter),
                     "status_raw": str(solution.status).split(".")[-1]},
    }
    return status, np.asarray(solution.x, dtype=float), wall, record


def solve(program, time_limit=CLARABEL_TIME_LIMIT_S):
    """Solve an exported program with Clarabel. Returns ``(solver_status, x, wall_time_s, solver_record)``."""
    import clarabel
    cones = _clarabel_cones(clarabel, program.n_zero, program.cones)
    return _clarabel_solve(program.A, program.b, program.c, cones, time_limit)


def lower_interior_field(program, time_limit=CLARABEL_TIME_LIMIT_S):
    """A strictly admissible Lower field for the checker's repair (interface 17.3 step 4), or ``None``.

    Maximises a uniform yield clearance ``t`` (``g <= -t`` at every vertex) over the same equalities, with
    ``lambda`` free, capped at ``t <= S0`` where ``S0`` is the largest strength or fixed-load stress scale.
    Returns ``(sigma_int, value_int, t)``. Optional, never trusted (interface 15.3).
    """
    import clarabel
    layout, geo = program.layout, program.layout["geometry"]
    n_rows, n_var = program.A.shape
    first_rows = program.n_zero + 6 * np.arange(4 * geo.E)
    column = sp.csc_matrix((np.full(len(first_rows), S2), (first_rows, np.zeros(len(first_rows), dtype=int))),
                           shape=(n_rows, 1))
    cap_row = sp.csc_matrix(([1.0], ([0], [n_var])), shape=(1, n_var + 1))
    matrix = sp.vstack([sp.hstack([program.A, column]), cap_row]).tocsc()
    t0 = np.abs(np.asarray(geo.arrays["t0"], float))
    l_ref = float(np.linalg.norm(geo.nodes.max(axis=0) - geo.nodes.min(axis=0)))
    cap = max(float(geo.k.max()), float(t0.max(initial=0.0)), l_ref * float(np.abs(geo.b0).max(initial=0.0)))
    query = program.query
    if query["form"] in ("ray", "fixed_components"):
        f0 = query_tables(program.problem, query)[1]
        areas = _footing_areas(geo)
        for fi in range(len(areas)):
            cap = max(cap, float(np.abs(f0[fi]).max(initial=0.0)) / areas[fi])
    if cap <= 0.0:
        return None
    b = np.concatenate([program.b, [cap]])
    c = np.zeros(n_var + 1)
    c[-1] = -1.0
    cones = _clarabel_cones(clarabel, program.n_zero, program.cones, extra_nonneg=1)
    status, x, _, _ = _clarabel_solve(matrix, b, c, cones, time_limit)
    if status not in contract.SOLVER_STATUS_WITH_VALUE or x[-1] <= 0.0:
        return None
    xs = x[:-1]
    sigma = xs[layout["sigma"]: layout["sigma"] + 24 * geo.E].reshape(geo.E, 4, 6)
    if layout["lambda"] is not None:
        value = float(xs[layout["lambda"]])
    else:
        value = float(layout["direction"].reshape(-1) @ (layout["resultant"] @ sigma.reshape(-1)))
    return sigma, value, float(x[-1])


def upper_conditions(geo, velocity):
    """Signed Upper conditions of interface 17.4 step 3 for a continuous (N, 3) or broken (E, 10, 3) field.

    Returns ``(c_vertex (E, 4), c_jump (F, 6), D_scale_vol, D_scale_jump)``: ``3 sqrt2 alpha |e| - theta`` at
    vertices with ``alpha > 0`` and ``tau |c_t| - n . c_B`` at coefficients with ``tau > 0``, ``-inf``
    elsewhere. The scales are per kind (amendment 1, item 10): ``D_scale_vol`` the largest ``|D_ei|_F``
    (strain-rate units), ``D_scale_jump`` the largest ``|c_B|`` (velocity units). A continuous field has no
    jump conditions and ``D_scale_jump = 0``. The producer's own evaluation; the checker recomputes it.
    """
    velocity = np.asarray(velocity, dtype=float)
    broken = velocity.ndim == 3
    field = velocity if broken else velocity[geo.tets]
    d = geo.strain_rates(field)
    theta, enorm = strain_split(d)
    frobenius = np.sqrt(np.sum(d[..., :3] ** 2, axis=-1) + 0.5 * np.sum(d[..., 3:] ** 2, axis=-1))
    c_vertex = np.where(geo.alpha[:, None] > 0.0, 3.0 * S2 * geo.alpha[:, None] * enorm - theta, -np.inf)
    scale_vol, scale_jump = float(frobenius.max(initial=0.0)), 0.0
    c_jump = np.full((geo.F, 6), -np.inf)
    if broken and geo.F:
        coeff = jump_coefficients(geo, field)
        n = geo.interior_normal[:, None, :]
        opening = np.sum(coeff * n, axis=2)
        slip = np.linalg.norm(coeff - opening[..., None] * n, axis=2)
        tau, _ = jump_constants(geo.alpha[geo.e_minus])
        c_jump = np.where(tau[:, None] > 0.0, tau[:, None] * slip - opening, -np.inf)
        scale_jump = float(np.linalg.norm(coeff, axis=2).max(initial=0.0))
    return c_vertex, c_jump, scale_vol, scale_jump


def repair_preview(geo, velocity, interior_velocity, apex_zero=None):
    """What the repair of interface 17.4 step 4 would need, evaluated WITHOUT the checker's projection.

    Each kind of condition (volume, jump) uses its own scale (amendment 1, item 10). ``P`` is the set of
    conditions above ``TOL_CONE_MARGIN`` times their kind's ``D_scale``. Returns the size of ``P``, the
    interior field's worst clearance on ``P`` and worst condition overall, each relative to its kind's
    ``D_scale_int``, ``mu`` and whether every requirement of step 4 holds. A preview, never a certificate.
    The conditions of a declared apex face ``apex_zero`` leave ``P`` and the clearance tests (19.2.1 item 7);
    the scales still run over every condition of their kind (item 6).
    """
    tol = contract.TOLERANCES
    cv, cj, vol, jump = upper_conditions(geo, velocity)
    iv, ij, vol_int, jump_int = upper_conditions(geo, interior_velocity)
    flow_z, jump_z = apex_zero_masks(geo, apex_zero)
    cand, inter, scale, scale_int = [], [], [], []
    for c, i, s, s_int, z in ((cv, iv, vol, vol_int, flow_z), (cj, ij, jump, jump_int, jump_z)):
        live = np.isfinite(c) & ~z
        cand.append(c[live])
        inter.append(i[live])
        scale.append(np.full(live.sum(), s))
        scale_int.append(np.full(live.sum(), s_int))
    cand, inter, scale, scale_int = (np.concatenate(item) for item in (cand, inter, scale, scale_int))
    in_p = cand > tol["TOL_CONE_MARGIN"] * scale
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(scale_int > 0.0, inter / scale_int, np.where(inter <= 0.0, 0.0, np.inf))
    mu = float(np.max(cand[in_p] / (cand[in_p] - inter[in_p]), initial=0.0)) if np.all(inter[in_p] < cand[in_p])         else math.inf
    ready = (bool(np.all(scale_int[in_p] > 0.0)) and bool(np.all(inter <= tol["TOL_CONE_MARGIN"] * scale_int))
             and bool(np.all(rel[in_p] <= -tol["TOL_INTERIOR_MIN"])) and mu <= tol["TOL_REPAIR_MAX"])
    return {
        "n_conditions": int(len(cand)), "n_repaired": int(in_p.sum()), "n_apex_zero": int(flow_z.sum() + jump_z.sum()),
        "d_scale_vol_int": vol_int, "d_scale_jump_int": jump_int,
        "worst_clearance_on_P": float(np.max(rel[in_p], initial=-math.inf)),
        "worst_condition": float(np.max(rel, initial=-math.inf)),
        "mu": mu, "ready": ready,
    }


# ---------------------------------------------------------------------------------------------------
# Declared apex face Z (interface 19.2.1, amendment 2 item A)
# ---------------------------------------------------------------------------------------------------

def apex_zero_array(rows):
    """``Z`` as the result array of 19.2.1 item 1: ``<i8`` ``(n_Z, 2)`` rows ``(kind, index)``, sorted, unique."""
    rows = np.asarray(rows, dtype=np.int64).reshape(-1, 2)
    unique = np.unique(rows, axis=0)                     # lexicographic order of (kind, index)
    if len(unique) != len(rows):
        fail("apex_zero_invalid", "rows are not unique")
    return np.ascontiguousarray(unique, dtype="<i8")


def apex_zero_masks(geo, apex_zero):
    """Boolean masks ``(E, 4)`` of the flow entries and ``(F, 6)`` of the jump entries of ``Z``."""
    flow, jump = np.zeros((geo.E, 4), dtype=bool), np.zeros((geo.F, 6), dtype=bool)
    if apex_zero is not None and len(apex_zero):
        rows = np.asarray(apex_zero, dtype=np.int64).reshape(-1, 2)
        size = np.where(rows[:, 0] == 0, flow.size, np.where(rows[:, 0] == 1, jump.size, 0))
        if np.any((rows[:, 1] < 0) | (rows[:, 1] >= size)):
            fail("apex_zero_invalid", "kind or index out of range")
        flow.reshape(-1)[rows[rows[:, 0] == 0, 1]] = True
        jump.reshape(-1)[rows[rows[:, 0] == 1, 1]] = True
    return flow, jump


def upper_inequality_conditions(geo, broken):
    """Every Upper inequality condition of 17.4 step 3, in the order and encoding of 19.2.1 item 1.

    ``key`` ``(n, 2)``: ``(0, 4 e + i)`` for the flow condition at a vertex with ``alpha_e > 0``, then
    ``(1, 6 F + B)`` for the jump condition at a coefficient with ``tau > 0`` (jump formulation only); the
    rows come out sorted. ``cone_row`` is the first row of the condition's cone inside the program's cone
    block, ``head`` the coefficient that turns a clearance ``s`` into a shift of the cone head
    (``(r - s / (3 sqrt2 alpha), z)`` in ``SOC6``, ``(rho - s / tau, y)`` in ``SOC3``), and ``natural`` the
    natural size of the clearance per unit velocity: ``1 / h_e`` for a strain rate, ``1`` for a jump.
    """
    E = geo.E
    elements, vertices = np.nonzero(np.repeat(geo.alpha[:, None] > 0.0, 4, axis=1))
    flat = 4 * elements + vertices
    key = [np.stack([np.zeros_like(flat), flat], axis=1)]
    cone_row = [6 * flat]
    head = [1.0 / (3.0 * S2 * geo.alpha[elements])]
    natural = [1.0 / geo.h[elements]]
    if broken and geo.F:
        tau, _ = jump_constants(geo.alpha[geo.e_minus])
        facets = np.nonzero(tau > 0.0)[0]
        f_idx, b_idx = np.repeat(facets, 6), np.tile(np.arange(6), len(facets))
        flat = 6 * f_idx + b_idx
        key.append(np.stack([np.ones_like(flat), flat], axis=1))
        cone_row.append(24 * E + 3 * flat)
        head.append(1.0 / tau[f_idx])
        natural.append(np.ones(len(flat)))
    return {"key": np.concatenate(key).astype(np.int64), "cone_row": np.concatenate(cone_row),
            "head": np.concatenate(head), "natural": np.concatenate(natural)}


def apex_zero_cone_rows(program, apex_zero):
    """Exact homogeneous rows that put every ``Z`` condition at its cone's apex, over the program variables.

    Flow entry ``(e, i)``: ``(r_ei, z_ei) = 0`` (six rows); with the flow row ``B_ei v = sqrt2 alpha r one +
    Q^T z`` this is ``d_ei = 0``. Jump entry ``(F, B)``: ``(rho, y1, y2)_FB = 0`` (three rows); with the
    jump rows ``[n, T1, T2] c_B = (tau rho, y1, y2)`` and ``tau > 0`` this is ``c_B = 0``. The same face as
    the velocity rows of 19.2.1 item 4 (:func:`face_rows`); here in the variables the solver sees.
    """
    layout = program.layout
    rows = apex_zero_array(apex_zero) if apex_zero is not None and len(apex_zero) else np.zeros((0, 2), np.int64)
    cols = []
    for kind, index in rows:
        if kind == 0:
            cols.append(layout["rz"] + 6 * index + np.arange(6))
        else:
            if not layout["broken"]:
                fail("apex_zero_invalid", "jump entries exist in the jump formulation only")
            cols.append(layout["jump"] + 3 * index + np.arange(3))
    cols = np.concatenate(cols) if cols else np.zeros(0, np.int64)
    return sp.csr_matrix((np.ones(len(cols)), (np.arange(len(cols)), cols)), shape=(len(cols), program.A.shape[1]))


def face_rows(program, apex_zero=None):
    """The projection rows of 17.4 step 2 over ``(v, V)``, augmented by ``Z`` (19.2.1 item 4): ``(M, g)``.

    Fixed, symmetry, footing and restrained rows and, for ``support``, ``V = d``, taken from the program
    without the work row (R9); isochory ``one^T B_ei v = 0`` where ``alpha = 0``; ``n . c_B = 0`` where
    ``tau = 0``; then six rows ``B_ei v`` per flow entry and three rows ``[n, T1, T2] c_B`` per jump entry
    of ``Z``, in ``Z`` order. Producer tooling for :func:`face_projection`; the checker builds its own.
    """
    layout, geo = program.layout, program.layout["geometry"]
    E, F, broken = geo.E, geo.F, layout["broken"]
    a = program.A.tocsr()[:, : layout["rz"]]
    jump_rows = 18 * F if broken else 0
    boundary = np.arange(24 * E + jump_rows, program.n_zero)
    if layout["work_row"] is not None:
        boundary = boundary[boundary != layout["work_row"]]
    blocks, rhs = [a[boundary]], [program.b[boundary]]
    vm = np.nonzero(geo.alpha <= 0.0)[0]
    if len(vm):
        flow = (24 * vm[:, None, None] + 6 * np.arange(4)[None, :, None] + np.arange(3)[None, None, :]).reshape(-1, 3)
        summing = sp.csr_matrix((np.ones(flow.size), (np.repeat(np.arange(len(flow)), 3), flow.reshape(-1))),
                                shape=(len(flow), program.A.shape[0]))
        blocks.append(summing @ program.A.tocsr()[:, : layout["rz"]])
        rhs.append(np.zeros(len(flow)))
    if broken and F:
        tau, _ = jump_constants(geo.alpha[geo.e_minus])
        facets = np.nonzero(tau <= 0.0)[0]
        normal_rows = (24 * E + 18 * facets[:, None] + 3 * np.arange(6)[None, :]).reshape(-1)
        blocks.append(a[normal_rows])
        rhs.append(np.zeros(len(normal_rows)))
    if apex_zero is not None and len(apex_zero):
        z_rows = []
        for kind, index in apex_zero_array(apex_zero):
            if kind == 0:
                z_rows.append(24 * (index // 4) + 6 * (index % 4) + np.arange(6))
            else:
                z_rows.append(24 * E + 3 * index + np.arange(3))
        z_rows = np.concatenate(z_rows)
        blocks.append(a[z_rows])
        rhs.append(np.zeros(len(z_rows)))
    return sp.vstack(blocks).tocsr(), np.concatenate(rhs)


def face_projection(program, velocity, footing_velocity, apex_zero=None):
    """Minimum-norm correction of a field onto :func:`face_rows` (unscaled least squares, LSQR).

    Returns ``(velocity, footing_velocity, move)`` with ``move`` the largest correction relative to the
    field's largest entry. Producer tooling: the checker's projection (17.1) uses its own scales.
    """
    from scipy.sparse.linalg import lsqr
    layout = program.layout
    n_f = len(program.problem["footings"])
    x = np.concatenate([np.asarray(velocity, float).reshape(-1), np.asarray(footing_velocity, float).reshape(-1)])
    matrix, g = face_rows(program, apex_zero)
    correction = lsqr(matrix, g - matrix @ x, atol=1e-15, btol=1e-15, iter_lim=50 * matrix.shape[1])[0]
    y = x + correction
    move = float(np.abs(correction).max(initial=0.0) / max(np.abs(x).max(initial=0.0), np.finfo(float).tiny))
    n_v = layout["V"]
    return y[:n_v].reshape(np.shape(velocity)), y[n_v: n_v + 6 * n_f].reshape(n_f, 6), move


class _ClearanceProgram:
    """The capped clearance program shared by :func:`find_apex_zero` and :func:`upper_interior_field`.

    Over the projection rows of 17.4 step 2 (every zero row except the work normalisation, R9), plus the
    exact rows of a declared face ``Z``, one bounded clearance ``s`` per inequality condition outside ``Z``
    and not forced to zero by the rows alone (:func:`row_forced_conditions`) enters its cone, so the
    condition reads ``c~ <= -s`` (the head shift of
    :func:`upper_inequality_conditions`). Velocities are boxed by ``v_box`` (and the footing rotations by
    ``v_box / length_scale``) unless ``V = d`` fixes them; a clearance's natural bound is ``natural * v_box``.
    :meth:`solve` maximises the weighted sum of ``s / s_max`` with ``held <= s <= cap``.
    """

    def __init__(self, program, apex_zero=None):
        layout, geo = program.layout, program.layout["geometry"]
        broken, n_f = layout["broken"], len(program.problem["footings"])
        self.program, self.n_var = program, program.A.shape[1]
        conditions = upper_inequality_conditions(geo, broken)
        face = apex_zero_array(apex_zero) if apex_zero is not None and len(apex_zero) else np.zeros((0, 2), np.int64)
        on_face = np.zeros(len(conditions["key"]), dtype=bool)
        if len(face):
            lookup = {tuple(row): j for j, row in enumerate(conditions["key"].tolist())}
            for row in face.tolist():
                if tuple(row) not in lookup:
                    fail("apex_zero_invalid", f"{row} is not an inequality condition (alpha > 0, tau > 0)")
                on_face[lookup[tuple(row)]] = True
        # Conditions the rows force to zero cannot be raised by any field: they get no clearance.
        by_rows = row_forced_conditions(program, conditions["key"]) & ~on_face
        self.row_forced_key = conditions["key"][by_rows]
        raised = ~on_face & ~by_rows
        self.key = conditions["key"][raised]
        n_s = len(self.key)
        self.n_s = n_s
        a_csr = program.A.tocsr()
        keep = np.ones(program.n_zero, dtype=bool)
        if layout["work_row"] is not None:
            keep[layout["work_row"]] = False
        a_eq = sp.vstack([a_csr[: program.n_zero][keep], apex_zero_cone_rows(program, face)]).tocsr()
        b_eq = np.concatenate([program.b[: program.n_zero][keep], np.zeros(a_eq.shape[0] - int(keep.sum()))])
        a_cone, b_cone = a_csr[program.n_zero:], program.b[program.n_zero:]
        ell = float(program.problem["length_scale"])
        v_box = 1.0
        if layout["form"] == "support":
            d = layout["direction"]
            v_box = max(1.0, float(np.abs(d[:, :3]).max(initial=0.0)), ell * float(np.abs(d[:, 3:]).max(initial=0.0)))
        self.v_box = v_box
        self.s_max = conditions["natural"][raised] * v_box
        slack_block = sp.csr_matrix((conditions["head"][raised], (conditions["cone_row"][raised], np.arange(n_s))),
                                    shape=(a_cone.shape[0], n_s))
        n_v = 30 * geo.E if broken else 3 * geo.N
        boxed, bounds = list(range(n_v)), [v_box] * n_v
        if layout["form"] != "support":
            for fi in range(n_f):
                boxed += [layout["V"] + 6 * fi + c for c in range(6)]
                bounds += [v_box] * 3 + [v_box / ell] * 3
        boxed, self.bounds = np.asarray(boxed), np.asarray(bounds)
        n_box = len(boxed)
        box = sp.csr_matrix((np.concatenate([np.ones(n_box), -np.ones(n_box)]),
                             (np.arange(2 * n_box), np.concatenate([boxed, boxed]))), shape=(2 * n_box, self.n_var + n_s))
        upper_s = sp.csr_matrix((np.ones(n_s), (np.arange(n_s), self.n_var + np.arange(n_s))),
                                shape=(n_s, self.n_var + n_s))
        self.matrix = sp.vstack([
            sp.hstack([a_eq, sp.csr_matrix((a_eq.shape[0], n_s))]),
            sp.hstack([a_cone, slack_block]),
            box, upper_s, -upper_s,
        ]).tocsc()
        self.b_eq, self.b_cone = b_eq, b_cone
        self.n_eq, self.n_box = a_eq.shape[0], n_box

    def solve(self, cap, held, weight, time_limit):
        import clarabel
        cones = _clarabel_cones(clarabel, self.n_eq, self.program.cones, extra_nonneg=2 * self.n_box + 2 * self.n_s)
        b = np.concatenate([self.b_eq, self.b_cone, self.bounds, self.bounds, cap, -held])
        c = np.zeros(self.n_var + self.n_s)
        c[self.n_var:] = -weight / self.s_max
        status, x, _, _ = _clarabel_solve(self.matrix, b, c, cones, time_limit)
        return status, x

    def field(self, x):
        fields = extract_fields(self.program, x[: self.n_var])
        velocity = fields["element_velocity"] if self.program.layout["broken"] else fields["velocity"]
        return velocity, fields["footing_velocity"]


# Apex-face finder (19.2.1; the evidence of 18.3). First the conditions that the step-2 rows force to zero on
# their own are set aside (row_forced_conditions): the checker already keeps them out of P, and Z rows for
# them would only be dependent rows of the augmented system (19.2.1 item 5). Then one solve of the clearance
# program, every other clearance capped at Z_CAP_FRACTION of its natural bound and the plain sum of
# s / s_max maximised. The admissible set without the work row is a convex cone (an affine set for support
# or fixed components) and every condition is sublinear, so the average of the one-condition maximisers
# raises every raisable condition at once: at the optimum a raisable condition sits at its cap. A condition
# left below Z_SLACK_FRACTION of its cap could not be raised by any admissible field; it is in Z when its
# FULL cone argument (|d_ei|_F, or |c_B|) is at most Z_APEX_ARGUMENT_REL of its kind's scale at that field
# (the 17.4 step 3 scales of the same field), so it sits at the apex and not on a ray of the cone boundary.
# A stuck condition off the apex is counted in info["stuck_off_apex"] and is not declared: Z cannot express
# it. So Z holds exactly the conditions forced onto the apex by the cones together with the rows.
#
# Calibration (producer settings, not certificate tolerances; a wrong Z only costs bound quality, never
# soundness, 19.2.1 item 10). upper3d_fixed_phi30 mesh 0 (323 tets), jump Upper, Clarabel 0.11.1 at its
# defaults, 2026-11-19, scratch probe of this finder (the same run the slow test
# test_apex_face_on_fixed_phi30_mesh_0 repeats):
# - the row probe found 230 row-forced conditions at <= 3.6e-16 of their kind's scale against >= 2.3e-3
#   for every other condition: Z_ROW_FORCED_REL = 1e-10 sits in a gap of twelve decades. Two seeded probes
#   guard against a probe that is small at one condition by chance;
# - of the other 4,278 clearances, 3,973 reached > 0.99 of their cap and 305 stayed below 7.7e-11 of their
#   natural bound, none in between: Z_SLACK_FRACTION (1e-6 of the cap, 1e-9 of the natural bound) sits in a
#   gap of more than six decades;
# - the stuck conditions' full arguments were at most 2.4e-9 of their kind's scale (the rows hold only to
#   the solver's tolerances), the smallest raised condition's argument 2.5e-3 of it: Z_APEX_ARGUMENT_REL =
#   1e-6 sits between, more than two decades from each. (The 18.3 diagnosis eliminated the rows exactly and
#   saw about 1e-11; this program keeps them as solver rows.)
# Z = 305 conditions (26 flow, 279 jump), the count of the 18.3 diagnosis.
Z_CAP_FRACTION = 1e-3
Z_SLACK_FRACTION = 1e-6
Z_APEX_ARGUMENT_REL = 1e-6
Z_ROW_FORCED_REL = 1e-10
Z_ROW_PROBES = 2
Z_ROW_PROBE_SEED = 237


def condition_arguments(geo, velocity):
    """Full cone arguments of every condition and the kind scales of 17.4 step 3.

    Returns ``(flow (E, 4), jump (F, 6), D_scale_vol, D_scale_jump)``: ``|d_ei|_F`` (tensor Frobenius norm,
    engineering shear halved) and ``|c_B|``. A continuous field has no jump arguments.
    """
    velocity = np.asarray(velocity, dtype=float)
    broken = velocity.ndim == 3
    field = velocity if broken else velocity[geo.tets]
    d = geo.strain_rates(field)
    flow = np.sqrt(np.sum(d[..., :3] ** 2, axis=-1) + 0.5 * np.sum(d[..., 3:] ** 2, axis=-1))
    jump = np.zeros((geo.F, 6))
    if broken and geo.F:
        jump = np.linalg.norm(jump_coefficients(geo, field), axis=2)
    return flow, jump, float(flow.max(initial=0.0)), float(jump.max(initial=0.0))


def rho_apex(geo, velocity, apex_zero):
    """``rho_apex`` of 19.2.1 item 6 at a field, or ``None`` without a face. Producer-side preview."""
    if apex_zero is None or not len(apex_zero):
        return None
    flow, jump, vol, jmp = condition_arguments(geo, velocity)
    flow_z, jump_z = apex_zero_masks(geo, apex_zero)
    parts = []
    if flow_z.any():
        parts.append(float(flow[flow_z].max()) / vol if vol > 0.0 else 0.0)
    if jump_z.any():
        parts.append(float(jump[jump_z].max()) / jmp if jmp > 0.0 else 0.0)
    return max(parts)


def find_apex_zero(program, time_limit=CLARABEL_TIME_LIMIT_S):
    """The conditions of an Upper program forced onto their cone's apex: the face ``Z`` of 19.2.1.

    Method and constants above :data:`Z_CAP_FRACTION`. Returns ``(apex_zero, info)``: ``apex_zero`` the
    ``<i8`` ``(n_Z, 2)`` array of :func:`apex_zero_array`, or ``None`` when no condition is forced (or the
    problem has no inequality condition). ``info`` counts the conditions, the row-forced ones set aside, the
    stuck ones and those stuck off the apex, with the relative arguments that bracket Z_APEX_ARGUMENT_REL.
    Refuses a Lower program (``apex_zero_lower``) and an unsolved clearance program (``no_interior_field``).
    """
    if program.bound != "upper":
        fail("apex_zero_lower", "a Lower program has no apex face")
    geo = program.layout["geometry"]
    clear = _ClearanceProgram(program)
    info = {"n_conditions": clear.n_s + len(clear.row_forced_key), "n_row_forced": len(clear.row_forced_key),
            "n_raisable_checked": clear.n_s, "n_stuck": 0, "n_apex_zero": 0, "stuck_off_apex": 0,
            "v_box": clear.v_box}
    if clear.n_s == 0:
        return None, info
    cap = Z_CAP_FRACTION * clear.s_max
    status, x = clear.solve(cap, np.zeros(clear.n_s), np.ones(clear.n_s), time_limit)
    info["status"] = status
    if status not in contract.SOLVER_STATUS_WITH_VALUE:
        fail("no_interior_field", f"apex-face clearance program not solved ({status})")
    s = x[clear.n_var:]
    stuck = s < Z_SLACK_FRACTION * cap
    at_cap = s > 0.99 * cap
    velocity, _ = clear.field(x)
    flow, jump, vol, jmp = condition_arguments(geo, velocity)
    key = clear.key
    flow_rows, jump_rows = key[:, 0] == 0, key[:, 0] == 1
    relative = np.zeros(len(key))
    if vol > 0.0:
        relative[flow_rows] = flow.reshape(-1)[key[flow_rows, 1]] / vol
    if jmp > 0.0:
        relative[jump_rows] = jump.reshape(-1)[key[jump_rows, 1]] / jmp
    apex = stuck & (relative <= Z_APEX_ARGUMENT_REL)
    off = stuck & ~apex

    def extreme(fn, values):
        return float(fn(values)) if len(values) else None

    info.update(n_at_cap=int(at_cap.sum()), n_partial=int(np.sum(~at_cap & ~stuck)), n_stuck=int(stuck.sum()),
                n_apex_zero=int(apex.sum()), stuck_off_apex=int(off.sum()),
                max_clearance_stuck=extreme(np.max, s[stuck] / clear.s_max[stuck]),
                max_argument_in_z=extreme(np.max, relative[apex]),
                min_argument_raised=extreme(np.min, relative[~stuck]),
                min_argument_stuck_off_apex=extreme(np.min, relative[off]))
    if not apex.any():
        return None, info
    return apex_zero_array(key[apex]), info


def row_forced_conditions(program, key, probes=Z_ROW_PROBES, seed=Z_ROW_PROBE_SEED):
    """Which conditions ``key`` the projection rows of 17.4 step 2 force to zero on their own.

    Such a condition has a zero argument on every field satisfying the homogeneous rows of
    :func:`face_rows` (without ``Z``), so the checker already keeps it out of ``P`` (17.4 step 4) and a
    ``Z`` row for it would only add dependent rows to the augmented system (19.2.1 item 5). Probe: a
    seeded Gaussian vector minus its minimum-norm LSQR correction onto those rows, ``probes`` times; a
    condition is row-forced when its full argument stays at or below ``Z_ROW_FORCED_REL`` of its kind's
    scale in every probe. Returns a boolean mask over ``key``.
    """
    from scipy.sparse.linalg import lsqr
    layout, geo = program.layout, program.layout["geometry"]
    matrix, _ = face_rows(program, None)
    rng = np.random.default_rng(seed)
    forced = np.ones(len(key), dtype=bool)
    flat = key[:, 1]
    flow_rows, jump_rows = key[:, 0] == 0, key[:, 0] == 1
    for _ in range(probes):
        r = rng.standard_normal(matrix.shape[1])
        x = r - lsqr(matrix, matrix @ r, atol=1e-16, btol=1e-16, iter_lim=20 * matrix.shape[1])[0]
        n_v = layout["V"]
        velocity = x[:n_v].reshape(geo.E, 10, 3) if layout["broken"] else x[:n_v].reshape(geo.N, 3)
        flow, jump, vol, jmp = condition_arguments(geo, velocity)
        relative = np.zeros(len(key))           # a kind whose scale is zero is zero everywhere: forced
        if vol > 0.0:
            relative[flow_rows] = flow.reshape(-1)[flat[flow_rows]] / vol
        if jmp > 0.0:
            relative[jump_rows] = jump.reshape(-1)[flat[jump_rows]] / jmp
        forced &= relative <= Z_ROW_FORCED_REL
    return forced


# Relative-interior construction of 17.4 step 4. Each slack is capped at a small fraction of its natural
# bound, so the capped sum is reachable jointly (the average of the one-condition maximisers clears every
# raisable condition) and an LP vertex does not starve the others. A slack that did not reach half its
# cap is re-raised next round with a cap ten times smaller while the reached ones are held at half their
# value. Producer settings, uncalibrated; the checker judges the field.
INTERIOR_CAP_FRACTION = 1e-3
INTERIOR_REACHED_FRACTION = 0.5
INTERIOR_HOLD_FRACTION = 0.5
INTERIOR_CAP_SHRINK = 0.1
INTERIOR_ROUNDS = 4


def upper_interior_field(program, time_limit=CLARABEL_TIME_LIMIT_S, weights=None, apex_zero=None):
    """A relative-interior Upper field for the checker's repair (interface 17.4 step 4).

    It satisfies the projection rows of 17.4 step 2: every zero row of the program except the work
    normalisation (R9), so fixed, roller, footing, restrained, ``V = d`` for ``support``, the flow rows
    (isochory where ``alpha = 0``) and the jump rows (``n . c_B = 0`` where ``tau = 0``). With a declared
    face ``apex_zero`` (19.2.1 item 7) it also satisfies the exact ``Z`` rows (:func:`apex_zero_cone_rows`),
    so it lies on that face, and only the conditions outside ``Z`` are raised. The clearance program is
    :class:`_ClearanceProgram`; its slacks are capped at ``INTERIOR_CAP_FRACTION`` of their natural bound.
    The sum of scaled slacks is maximised; slacks that did not reach their cap are re-raised with smaller
    caps while the others are held (the producer route of 17.4 step 4) until none moves. Conditions the rows
    force to zero stay at zero; on a problem with an implicit face (18.3) they are the ones
    :func:`find_apex_zero` declares. ``weights`` (one per raised slack, flow first, in
    :func:`upper_inequality_conditions` order without ``Z``) is a diagnostic override of the objective.

    Returns ``(velocity, footing_velocity, info)`` with ``velocity`` shaped as the result's field. Refuses
    with ``no_interior_field`` when the problem has no inequality condition (every ``alpha = 0``).
    """
    geo = program.layout["geometry"]
    if program.bound != "upper":
        fail("no_interior_field", "upper_interior_field needs an Upper program")
    if not np.any(geo.alpha > 0.0):
        fail("no_interior_field", "every alpha = 0: no inequality condition exists, so no repair is needed")
    clear = _ClearanceProgram(program, apex_zero)
    n_s, s_max = clear.n_s, clear.s_max
    weight = np.ones(n_s) if weights is None else np.asarray(weights, dtype=float)
    cap = INTERIOR_CAP_FRACTION * s_max
    held = np.zeros(n_s)
    active = weight > 0.0
    best, rounds, status = None, 0, "not_run"
    for rounds in range(1, INTERIOR_ROUNDS + 1):
        status, x = clear.solve(cap, held, np.where(active, weight, 0.0), time_limit)
        if status not in contract.SOLVER_STATUS_WITH_VALUE:
            break
        best = x
        s = x[clear.n_var:]
        reached = s >= INTERIOR_REACHED_FRACTION * cap
        held = np.where(reached | ~active, INTERIOR_HOLD_FRACTION * s, 0.0)
        if not np.any(active & reached) and rounds > 1:
            break
        active = active & ~reached
        if not np.any(active):
            break
        cap = np.where(active, INTERIOR_CAP_SHRINK * cap, cap)
    if best is None:
        fail("no_interior_field", f"interior program not solved ({status})")
    velocity, footing = clear.field(best)
    s = best[clear.n_var:]
    info = {"rounds": rounds, "n_conditions": n_s, "n_unraised": int(np.sum(s < 1e-6 * s_max)), "v_box": clear.v_box,
            "n_apex_zero": 0 if apex_zero is None else int(len(apex_zero))}
    return velocity, footing, info


def _footing_areas(geo):
    areas = np.zeros(len(geo.problem["footings"]))
    for row in range(len(geo.facets)):
        item = geo.sets[int(geo.facet_set[row])]
        if item["kind"] == "footing":
            areas[int(item["footing"])] += geo.boundary_area[row]
    return areas


# ---------------------------------------------------------------------------------------------------
# Cases
# ---------------------------------------------------------------------------------------------------

def case_problem(name, mesh_index=0):
    """Problem document and arrays of a case on one mesh, through the contract's builder."""
    cases = contract.parse_cases()
    return contract.build_case_problem(cases[name], mesh_index)


def run_case(name, mesh_index, formulations, directory=None, solve_it=False, interior=False,
             time_limit=CLARABEL_TIME_LIMIT_S, apex_face=True):
    """Build, validate, assemble, optionally solve and export one case. Returns ``{formulation: (program, result)}``.

    Upper solves of a problem with any ``alpha > 0`` always carry an interior field
    (:func:`upper_interior_field`). With ``apex_face`` (the default) they first look for the conditions
    forced onto the apex (:func:`find_apex_zero`), declare them as ``apex_zero`` and build the interior
    field on that face (interface 19.2.1). ``interior=True`` adds the Lower interior field; it refuses an
    Upper formulation on a problem with every ``alpha = 0``, which has no inequality to repair.
    """
    problem, arrays = case_problem(name, mesh_index)
    validate(problem, arrays)
    geo = Geometry(problem, arrays)
    frictional = bool(np.any(geo.alpha > 0.0))
    for formulation in formulations:
        check_formulation(formulation, problem)
        if interior and ADMITTED[formulation] == "upper" and not frictional:
            fail("no_interior_field", f"{formulation}: every alpha = 0, so no Upper interior field exists")
    if directory is not None:
        export_problem(directory, problem, arrays)
    out = {}
    for formulation in formulations:
        program = assemble(problem, arrays, formulation, "q", geometry=geo, validated=True)
        if directory is not None:
            export_program(directory, program)
        result = None
        if solve_it:
            status, x, wall, record = solve(program, time_limit)
            extra, face = None, None
            solved_ok = status in contract.SOLVER_STATUS_WITH_VALUE
            if interior and program.bound == "lower" and solved_ok:
                found = lower_interior_field(program, time_limit)
                extra = None if found is None else found[:2]
            if program.bound == "upper" and frictional and solved_ok:
                if apex_face:
                    face = find_apex_zero(program, time_limit)[0]
                extra = upper_interior_field(program, time_limit, apex_zero=face)[:2]
            result = (result_document(program, status, x, record, wall, extra, face) if directory is None
                      else export_result(directory, program, status, x, record, wall, extra, face))
        out[formulation] = (program, result)
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("case")
    parser.add_argument("mesh", type=int)
    parser.add_argument("formulations", nargs="+")
    parser.add_argument("--out", required=True)
    parser.add_argument("--solve", action="store_true")
    parser.add_argument("--interior", action="store_true")
    parser.add_argument("--time-limit", type=float, default=CLARABEL_TIME_LIMIT_S)
    parser.add_argument("--no-apex-face", action="store_true", help="do not look for or declare an apex face Z")
    args = parser.parse_args(argv)
    runs = run_case(args.case, args.mesh, args.formulations, args.out, args.solve, args.interior, args.time_limit,
                    not args.no_apex_face)
    for formulation, (program, result) in runs.items():
        line = f"{args.case} m{args.mesh} {formulation}: n_var={program.A.shape[1]} n_rows={program.A.shape[0]}"
        if result is not None:
            doc = result[0]
            line += f" status={doc['solver_status']} estimate={doc['value']} wall={doc['wall_time_s']:.2f}s"
        print(line)


if __name__ == "__main__":
    main()
