"""WP-237a: can the section 17 certificate procedure certify a solver-produced field?

PROCEDURE DEMONSTRATION, NOT THE REFERENCE AND NOT THE INDEPENDENT CHECKER.

The authors of WP-237b (reference), WP-237c (checker), WP-237e (external reference) and WP-237f (C++) must not
read this file: they implement from tests/limit_interface.md. It exists so that the contract's acceptance rule is
known to be reachable before the contract is frozen. It assembles the programs of interface 15.4 for the two
formulations it needs, solves them with Clarabel (Apache-2.0), and runs interface 17.3 and 17.4 on the result.

Producer and checker here share one author and one set of operators, so a pass shows only that the procedure is
reachable. The certificates carry the reserved checker name `procedure-demo`, which no gate accepts, and every
value stays an ESTIMATE for the limit-analysis program.

Usage: python certificate_demo.py [--write] [case ...]
Needs numpy, scipy and clarabel. --write refreshes certificate_demo_record.json beside this file.
The dense factorisations limit it to a few hundred tetrahedra; the method named for large meshes is interface 17.2.
"""
import copy
import json
import sys
import time
from pathlib import Path

import clarabel
import numpy as np
import scipy
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "tests"))
import test_limit_contract as T  # noqa: E402

TOL = T.TOLERANCES
PAIRS = ((0, 1), (1, 2), (2, 0))
BERNSTEIN = np.eye(6)
for q, (a, b) in enumerate(PAIRS):
    BERNSTEIN[3 + q, 3 + q] = 2.0
    BERNSTEIN[3 + q, a] = BERNSTEIN[3 + q, b] = -0.5
SKEW = ((1, 2), (2, 0), (0, 1))


class Rows:
    def __init__(self):
        self.i, self.j, self.v, self.b, self.n = [], [], [], [], 0

    def add(self, cols, vals, rhs=0.0):
        for c, v in zip(cols, vals):
            if v != 0.0:
                self.i.append(self.n), self.j.append(int(c)), self.v.append(float(v))
        self.b.append(float(rhs))
        self.n += 1

    def matrix(self, n_cols):
        return sp.csr_matrix((self.v, (self.i, self.j)), shape=(self.n, n_cols)), np.array(self.b)


class Model:
    """Geometry of one problem, built from the problem document and arrays alone."""

    def __init__(self, problem, arrays):
        self.problem, self.arrays = problem, arrays
        self.query = problem["queries"][0]
        self.nodes, self.tets = arrays["nodes"], arrays["tets"]
        self.n_e, self.n_f = len(self.tets), len(problem["footings"])
        self.grads, self.vol = np.zeros((self.n_e, 4, 3)), np.zeros(self.n_e)
        for e, tet in enumerate(self.tets):
            self.grads[e], self.vol[e] = T.tet_gradients(self.nodes[tet[:4]])
        self.interior, _ = T.mesh_facets(self.tets)
        self.region = [problem["regions"][int(r)] for r in arrays["element_region"]]
        self.k = np.array([r["k"] for r in self.region])
        self.alpha = np.array([r.get("alpha", 0.0) for r in self.region])
        self.F0 = np.array(self.query["F0"], dtype=float).reshape(self.n_f, 6)
        self.F1 = np.array(self.query["F1"], dtype=float).reshape(self.n_f, 6)
        self.restrained = {tuple(item) for item in self.query["restrained"]}
        self.kinds = [item["kind"] for item in problem["boundary_sets"]]
        self.ell = problem["length_scale"]

    def facet(self, e, f):
        return T.facet_geometry(self.nodes[self.tets[e][:4]], int(f))

    def boundary(self):
        for row, ((e, f), index) in enumerate(zip(self.arrays["facets"], self.arrays["facet_set"])):
            yield row, int(e), int(f), int(index)

    def s_ref_common(self, value):
        a = self.arrays
        out = max(float(self.k.max()), float(np.abs(a["t0"]).max(initial=0.0)),
                  abs(value) * float(np.abs(a["t1"]).max(initial=0.0)))
        span = float(np.linalg.norm(self.nodes.max(axis=0) - self.nodes.min(axis=0)))
        for r in self.problem["regions"]:
            out = max(out, span * (np.abs(r["b0"]).max() + abs(value) * np.abs(r["b1"]).max()))
        areas = np.zeros(self.n_f)
        for _, e, f, index in self.boundary():
            item = self.problem["boundary_sets"][index]
            if item["kind"] == "footing":
                areas[item["footing"]] += self.facet(e, f)[1]
        for f in range(self.n_f):
            force = np.abs(self.F0[f] + value * self.F1[f]).max()
            out = max(out, force / areas[f])
        return out


# ---------------------------------------------------------------------------------------------
# Interface 17.1 and 17.2: reduction, projection, backward error
# ---------------------------------------------------------------------------------------------

def literal_zero_columns(A, g):
    """Columns forced to zero by homogeneous rows with one remaining entry, to a fixed point."""
    A = A.tocsr()
    csc = A.tocsc()
    alive = np.ones(A.shape[1], dtype=bool)
    count = np.diff(A.indptr).copy()
    stack = [i for i in np.nonzero((count == 1) & (g == 0.0))[0]]
    while stack:
        row = stack.pop()
        if count[row] != 1:
            continue
        cols = A.indices[A.indptr[row]:A.indptr[row + 1]]
        col = cols[alive[cols]][0]
        alive[col] = False
        for other in csc.indices[csc.indptr[col]:csc.indptr[col + 1]]:
            count[other] -= 1
            if count[other] == 1 and g[other] == 0.0:
                stack.append(other)
    return ~alive


def backward_error(A, g, x, scale):
    denominator = np.abs(A) @ scale + np.abs(g)
    residual = np.abs(A @ x - g)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(denominator > 0.0, residual / denominator, np.where(residual > 0.0, np.inf, 0.0))
    return float(ratio.max(initial=0.0))


def project(A, g, x, scale, sweeps=4):
    """Literal zeros, then the minimum-norm correction in scaled unknowns with iterative refinement."""
    A = A.tocsr()
    zero = literal_zero_columns(A, g)
    out = x.copy()
    out[zero] = 0.0
    free = np.nonzero(~zero)[0]
    reduced = A[:, free]
    keep = np.nonzero(np.diff(reduced.indptr) > 0)[0]
    reduced, rhs = reduced[keep], g[keep]
    column = reduced @ sp.diags(scale[free])
    row_scale = np.abs(column).max(axis=1).toarray().ravel()
    hat = (sp.diags(1.0 / row_scale) @ column).toarray()
    left, sing, right = np.linalg.svd(hat, full_matrices=False)
    rank = int(np.sum(sing > TOL["TOL_RANK"] * sing[0]))
    kappa = float(sing[0] / sing[rank - 1])
    left, sing_kept, right = left[:, :rank], sing[:rank], right[:rank]
    for _ in range(sweeps):
        residual = (reduced @ out[free] - rhs) / row_scale
        out[free] -= scale[free] * (right.T @ ((left.T @ residual) / sing_kept))
    return out, {"kappa": kappa, "rank": rank, "rows": hat.shape[0], "cols": hat.shape[1],
                 "literal_zero": int(zero.sum()), "zero_mask": zero}


# ---------------------------------------------------------------------------------------------
# Lower: interface 9, 15.4 and 17.3
# ---------------------------------------------------------------------------------------------

def sv(e, i):
    return 1 + 24 * e + 6 * i


def lower_rows(m):
    rows = Rows()
    resultant = sp.lil_matrix((6 * m.n_f, 1 + 24 * m.n_e))
    for e in range(m.n_e):
        region = m.region[e]
        block = np.hstack([T.voigt_n(m.grads[e, i]) for i in range(4)])
        for c in range(3):
            rows.add([0] + list(range(sv(e, 0), sv(e, 0) + 24)), [region["b1"][c]] + list(block[c]), -region["b0"][c])
    for (e1, f1), (e2, _) in m.interior:
        n, _ = m.facet(e1, f1)
        N = T.voigt_n(n)
        local2 = [int(v) for v in m.tets[e2]]
        for a in T.FACET_VERTICES[f1]:
            b = local2.index(int(m.tets[e1][a]))
            for c in range(3):
                rows.add(list(range(sv(e1, a), sv(e1, a) + 6)) + list(range(sv(e2, b), sv(e2, b) + 6)),
                         list(N[c]) + list(-N[c]))
    for row, e, f, index in m.boundary():
        kind = m.kinds[index]
        n, area = m.facet(e, f)
        N = T.voigt_n(n)
        vertices = T.FACET_VERTICES[f]
        if kind in ("free", "traction"):
            for j, a in enumerate(vertices):
                for c in range(3):
                    rows.add([0] + list(range(sv(e, a), sv(e, a) + 6)), [-m.arrays["t1"][row, j, c]] + list(N[c]),
                             m.arrays["t0"][row, j, c])
        elif kind == "symmetry":
            for tangent in T.tangent_basis(n):
                for a in vertices:
                    rows.add(list(range(sv(e, a), sv(e, a) + 6)), list(tangent @ N))
        elif kind == "footing":
            f_index = m.problem["boundary_sets"][index]["footing"]
            x_ref = np.array(m.problem["footings"][f_index]["x_ref"])
            x = m.nodes[m.tets[e][list(vertices)]]
            for j, b in enumerate(vertices):
                column = sv(e, b)
                for c in range(6):
                    for d in range(3):
                        resultant[6 * f_index + d, column + c] += area / 3.0 * N[d, c]
                arm = sum((area / 6.0 if p == j else area / 12.0) * (x[p] - x_ref) for p in range(3))
                moment = np.cross(arm, N.T).T  # column-wise arm x (N sigma)
                for c in range(6):
                    for d in range(3):
                        resultant[6 * f_index + 3 + d, column + c] += moment[d, c]
    resultant = resultant.tocsr()
    for f in range(m.n_f):
        for c in range(6):
            if (f, c) in m.restrained:
                continue
            line = resultant[6 * f + c]
            rows.add([0] + list(line.indices), [-m.F1[f, c]] + list(line.data), m.F0[f, c])
    A, g = rows.matrix(1 + 24 * m.n_e)
    return A, g, resultant


def solve_conic(P_c, A_eq, g_eq, A_cone, b_cone, cone_dims):
    A = sp.vstack([A_eq, A_cone]).tocsc()
    b = np.concatenate([g_eq, b_cone])
    cones = [clarabel.ZeroConeT(A_eq.shape[0])] + [clarabel.SecondOrderConeT(d) for d in cone_dims]
    settings = clarabel.DefaultSettings()
    settings.verbose = False
    n = A.shape[1]
    solution = clarabel.DefaultSolver(sp.csc_matrix((n, n)), P_c, A, b, cones, settings).solve()
    return np.array(solution.x), str(solution.status)


def solve_lower(m, A, g):
    rows = Rows()
    n = A.shape[1]
    for e in range(m.n_e):
        alpha, k = m.alpha[e], m.k[e]
        for i in range(4):
            column = list(range(sv(e, i), sv(e, i) + 6))
            rows.add(column, list(T.S2 * alpha * T.ONE), T.S2 * k)
            for q in range(5):
                rows.add(column, list(-T.Q[q]))
    G, h = rows.matrix(n)
    c = np.zeros(n)
    c[0] = -1.0
    return solve_conic(c, A, g, G, h, [6] * (4 * m.n_e))


def yield_values(m, sigma):
    s = sigma.reshape(-1, 6)
    tensor = np.zeros((len(s), 3, 3))
    for (p, q), c in {(0, 0): 0, (1, 1): 1, (2, 2): 2, (0, 1): 3, (1, 2): 4, (0, 2): 5}.items():
        tensor[:, p, q] = tensor[:, q, p] = s[:, c]
    w = np.linalg.eigvalsh(tensor)
    j2 = ((w[:, 0] - w[:, 1]) ** 2 + (w[:, 1] - w[:, 2]) ** 2 + (w[:, 2] - w[:, 0]) ** 2) / 6.0
    return np.sqrt(j2) + np.repeat(m.alpha, 4) * w.sum(axis=1) - np.repeat(m.k, 4)


def certify_lower(m, A, g, x_c, resultant):
    """Interface 17.3 for a ray query whose fixed loads are all zero (zero stress is then interior)."""
    value = float(x_c[0])
    s_ref = max(m.s_ref_common(value), float(np.abs(x_c[1:]).max()))
    scale = np.full(A.shape[1], s_ref)
    scale[0] = abs(value) if value != 0.0 else 1.0
    x_p, info = project(A, g, x_c, scale)
    move = float(np.max(np.abs(x_p - x_c) / scale))
    if move > TOL["TOL_PROJECTION_MOVE"]:
        return None, "projection_move_too_large", {"projection_move": move}
    g_y = yield_values(m, x_p[1:])
    rho_before = float(g_y.max() / s_ref)
    if rho_before > TOL["TOL_CONE_ACCEPT"]:
        return None, "cone_violation_too_large", {"rho_before": rho_before}
    margin = TOL["TOL_CONE_MARGIN"] * s_ref
    mu = 0.0
    x_f = x_p.copy()
    if np.any(g_y > -margin):
        if np.any(g != 0.0) or np.any(m.k <= 0.0):
            return None, "no_interior_field", {"rho_before": rho_before}  # this demonstration supplies no interior_sigma
        interior = -np.repeat(m.k, 4)
        mu = float(np.max((g_y + margin) / (g_y - interior)))
        if mu > TOL["TOL_REPAIR_MAX"]:
            return None, "repair_too_large", {"repair_mu": mu}
        x_f = (1.0 - mu) * x_p
    g_final = yield_values(m, x_f[1:])
    omega = backward_error(A, g, x_f, scale)
    kappa = info["kappa"]
    tau = T.tau_eq()
    metrics = {
        "omega_eq": omega, "tau_eq": tau, "kappa_eq": kappa, "projection_move": move,
        "rho_cone": float(g_final.max() / s_ref), "repair_mu": mu, "rho_before_repair": rho_before,
        "rank": info["rank"], "rows": info["rows"], "cols": info["cols"], "literal_zero": info["literal_zero"],
    }
    if np.any(x_f[info["zero_mask"]] != 0.0):
        return None, "literal_zero_violated", metrics
    if omega > tau:
        return None, "equality_residual_above_floor", metrics
    if kappa > TOL["TOL_KAPPA_CAP"] or kappa * omega > TOL["TOL_FLOOR_MAX"]:
        return None, "equality_ill_conditioned", metrics
    if g_final.max() > 0.0:
        return None, "cone_violation_too_large", metrics
    return float(x_f[0]), None, metrics


# ---------------------------------------------------------------------------------------------
# Upper, jump formulation: interface 12, 15.4 and 17.4
# ---------------------------------------------------------------------------------------------

class Jump:
    def __init__(self, m):
        self.m = m
        self.n_aff = 30 * m.n_e + 6 * m.n_f
        self.B = np.zeros((m.n_e, 4, 6, 30))
        for e in range(m.n_e):
            for i in range(4):
                dn = T.p2_gradients_at_vertex(m.grads[e], i)
                self.B[e, i] = np.hstack([T.voigt_n(dn[a]).T for a in range(10)])
        self.facets = []
        for (e1, f1), (e2, _) in m.interior:
            n, area = m.facet(e1, f1)
            t1, t2 = T.tangent_basis(n)
            local1 = list(T.FACET_NODES[f1])
            order = [int(v) for v in m.tets[e2]]
            local2 = [order.index(int(m.tets[e1][a])) for a in local1]
            region = m.region[e1]
            tau, root = T.jump_constants(region.get("alpha", 0.0))
            self.facets.append((e1, e2, local1, local2, n, t1, t2, area, tau, region["k"] / root))
        self.n_jump = len(self.facets)
        self.J = self._jump_operator()
        self.flow = self._flow_operator()
        self.W0, self.W1 = self._work()

    def dof(self, e, a, c):
        return 30 * e + 3 * a + c

    def footing_dof(self, f, c):
        return 30 * self.m.n_e + 6 * f + c

    def _flow_operator(self):
        i, j, v = [], [], []
        for e in range(self.m.n_e):
            for vertex in range(4):
                for c in range(6):
                    for q in np.nonzero(self.B[e, vertex, c])[0]:
                        i.append(24 * e + 6 * vertex + c), j.append(30 * e + int(q)), v.append(self.B[e, vertex, c, q])
        return sp.csr_matrix((v, (i, j)), shape=(24 * self.m.n_e, self.n_aff))

    def _jump_operator(self):
        """c_B components in the global frame: 18 rows per facet."""
        i, j, v = [], [], []
        for index, (e1, e2, local1, local2, *_rest) in enumerate(self.facets):
            for coefficient in range(6):
                for q in range(6):
                    w = BERNSTEIN[coefficient, q]
                    if w == 0.0:
                        continue
                    for c in range(3):
                        row = 18 * index + 3 * coefficient + c
                        i += [row, row]
                        j += [self.dof(e2, local2[q], c), self.dof(e1, local1[q], c)]
                        v += [w, -w]
        return sp.csr_matrix((v, (i, j)), shape=(18 * self.n_jump, self.n_aff))

    def _work(self):
        m = self.m
        W0, W1 = np.zeros(self.n_aff), np.zeros(self.n_aff)
        weights = np.array([-0.05] * 4 + [0.2] * 6)
        for e in range(m.n_e):
            region = m.region[e]
            for a in range(10):
                for c in range(3):
                    W0[self.dof(e, a, c)] += m.vol[e] * weights[a] * region["b0"][c]
                    W1[self.dof(e, a, c)] += m.vol[e] * weights[a] * region["b1"][c]
        for row, e, f, _ in m.boundary():
            table = T.facet_work_table(m.facet(e, f)[1])
            nodes = T.FACET_NODES[f]
            for i in range(3):
                for q in range(6):
                    for c in range(3):
                        W0[self.dof(e, nodes[q], c)] += table[i, q] * m.arrays["t0"][row, i, c]
                        W1[self.dof(e, nodes[q], c)] += table[i, q] * m.arrays["t1"][row, i, c]
        for f in range(m.n_f):
            for c in range(6):
                W0[self.footing_dof(f, c)] += m.F0[f, c]
                W1[self.footing_dof(f, c)] += m.F1[f, c]
        return W0, W1

    def boundary_rows(self, rows):
        m = self.m
        sets = [set() for _ in m.kinds]
        normals = {}
        for _, e, f, index in m.boundary():
            sets[index].update((e, a) for a in T.FACET_NODES[f])
            normals.setdefault(index, m.facet(e, f)[0])
        fixed = set().union(*[s for s, kind in zip(sets, m.kinds) if kind == "fixed"], set())
        for e, a in sorted(fixed):
            for c in range(3):
                rows.add([self.dof(e, a, c)], [1.0])
        for index, kind in enumerate(m.kinds):
            if kind == "symmetry":
                n = normals[index]
                for e, a in sorted(sets[index] - fixed):
                    rows.add([self.dof(e, a, c) for c in range(3)], list(n))
        for index, item in enumerate(m.problem["boundary_sets"]):
            if item["kind"] != "footing":
                continue
            f = item["footing"]
            x_ref = np.array(m.problem["footings"][f]["x_ref"])
            for e, a in sorted(sets[index]):
                r = m.nodes[m.tets[e][a]] - x_ref
                for c in range(3):
                    p, q = SKEW[c]
                    # (w x r)_c = w_p r_q - w_q r_p
                    rows.add([self.dof(e, a, c), self.footing_dof(f, c), self.footing_dof(f, 3 + p),
                              self.footing_dof(f, 3 + q)], [1.0, -1.0, -r[q], r[p]])
        for f, c in sorted(m.restrained):
            rows.add([self.footing_dof(f, c)], [1.0])

    def certificate_rows(self):
        """Interface 17.4 step 2. The work normalisation is not a row."""
        m = self.m
        rows = Rows()
        self.boundary_rows(rows)
        for e in range(m.n_e):
            if m.alpha[e] == 0.0:
                for i in range(4):
                    theta = T.ONE @ self.B[e, i]
                    rows.add([self.dof(e, 0, 0) + q for q in range(30)], list(theta))
        J = self.J.tocsr()
        for index, facet in enumerate(self.facets):
            n, tau = facet[4], facet[8]
            if tau == 0.0:
                for coefficient in range(6):
                    start = 18 * index + 3 * coefficient
                    line = n[0] * J[start] + n[1] * J[start + 1] + n[2] * J[start + 2]
                    rows.add(line.indices, line.data)
        return rows.matrix(self.n_aff)

    def interior_field(self, rounds=5):
        """A point of the relative interior: raise the slacks still at zero while holding the others."""
        weights, floors, x, status, previous = None, None, None, None, None
        for _ in range(rounds):
            x, slack, status = self.solve(interior=True, weights=weights, floors=floors)
            zero = slack < 1e-7
            self.interior_zero_slacks = int(zero.sum())
            if weights is not None and np.array_equal(zero, previous):
                break
            previous = zero
            weights = zero.astype(float)
            floors = np.where(zero, 0.0, 0.5 * slack)
            if not zero.any():
                break
        return x, status

    def solve(self, interior=False, weights=None, floors=None):
        """The Upper program, or with interior=True a field inside the cones, one bounded slack per condition."""
        m = self.m
        n_aff, n_r = self.n_aff, 24 * m.n_e
        n_base = n_aff + n_r + 18 * self.n_jump
        n_slack = (4 * m.n_e + 6 * self.n_jump) if interior else 0
        n = n_base + n_slack
        def r_col(e, v):
            return n_aff + 24 * e + 6 * v

        def y_col(f, b):
            return n_aff + n_r + 18 * f + 3 * b

        def s_vol(e, v):
            return n_base + 4 * e + v

        def s_jump(f, b):
            return n_base + 4 * m.n_e + 6 * f + b
        rows = Rows()
        flow = self.flow.tocsr()
        for e in range(m.n_e):
            for v in range(4):
                for c in range(6):
                    line = flow[24 * e + 6 * v + c]
                    coeff = np.concatenate([[-T.S2 * m.alpha[e] * T.ONE[c]], -T.Q[:, c]])
                    rows.add(list(line.indices) + [r_col(e, v) + q for q in range(6)], list(line.data) + list(coeff))
        J = self.J.tocsr()
        for index, facet in enumerate(self.facets):
            normal, t1, t2, tau = facet[4], facet[5], facet[6], facet[8]
            for b in range(6):
                start = 18 * index + 3 * b
                comps = [J[start], J[start + 1], J[start + 2]]
                for vec, col, weight in ((normal, 0, tau), (t1, 1, 1.0), (t2, 2, 1.0)):
                    line = vec[0] * comps[0] + vec[1] * comps[1] + vec[2] * comps[2]
                    cols, vals = list(line.indices) + [y_col(index, b) + col], list(line.data) + [-weight]
                    if interior and col == 0:
                        cols.append(s_jump(index, b)), vals.append(-1.0)  # n.c = tau rho + slack
                    rows.add(cols, vals)
        self.boundary_rows(rows)
        if not interior:
            nz = np.nonzero(self.W1)[0]
            rows.add(nz, self.W1[nz], 1.0)
        A, g = rows.matrix(n)
        cone = Rows()
        dims = []
        for e in range(m.n_e):
            for v in range(4):
                for q in range(6):
                    cols, vals = [r_col(e, v) + q], [-1.0]
                    if interior and q == 0:
                        cols.append(s_vol(e, v)), vals.append(1.0)  # |z| <= r - slack
                    cone.add(cols, vals)
                dims.append(6)
        for index in range(self.n_jump):
            for b in range(6):
                for q in range(3):
                    cone.add([y_col(index, b) + q], [-1.0])
                dims.append(3)
        G, h = cone.matrix(n)
        c = np.zeros(n)
        if interior:
            c[n_base:] = -1.0 if weights is None else -weights
            box = Rows()
            for index, col in enumerate(range(n_base, n)):
                box.add([col], [-1.0], 0.0 if floors is None else -floors[index])
                box.add([col], [1.0], 1.0)
            for e in range(m.n_e):
                for v in range(4):
                    box.add([r_col(e, v)], [1.0], 2.0)
            B, d = box.matrix(n)
            cones = [clarabel.ZeroConeT(A.shape[0])] + [clarabel.SecondOrderConeT(q) for q in dims]
            cones.append(clarabel.NonnegativeConeT(B.shape[0]))
            settings = clarabel.DefaultSettings()
            settings.verbose = False
            solution = clarabel.DefaultSolver(sp.csc_matrix((n, n)), c, sp.vstack([A, G, B]).tocsc(),
                                              np.concatenate([g, h, d]), cones, settings).solve()
            solved = np.array(solution.x)
            return solved[:n_aff], solved[n_base:], str(solution.status)
        c[:n_aff] = -self.W0
        for e in range(m.n_e):
            for v in range(4):
                c[r_col(e, v)] = m.vol[e] / 4.0 * T.S2 * m.k[e]
        for index, facet in enumerate(self.facets):
            for b in range(6):
                c[y_col(index, b)] = facet[7] / 6.0 * facet[9]
        x, status = solve_conic(c, A, g, G, h, dims)
        return x[:n_aff], status

    def conditions(self, x):
        """Signed conditions of interface 17.4 and the two dissipation sums."""
        m = self.m
        d = (self.flow @ x).reshape(-1, 4, 6)
        theta = d[..., 0] + d[..., 1] + d[..., 2]
        dev = d[..., :3] - theta[..., None] / 3.0
        enorm = np.sqrt((dev * dev).sum(axis=-1) + 0.5 * (d[..., 3:] ** 2).sum(axis=-1))
        alpha = m.alpha[:, None]
        k = m.k[:, None]
        c_vol = np.where(alpha > 0.0, 3.0 * T.S2 * alpha * enorm - theta, -np.inf)
        with np.errstate(divide="ignore", invalid="ignore"):
            pi_vol = np.where(alpha > 0.0, k * theta / (3.0 * np.where(alpha > 0.0, alpha, 1.0)), T.S2 * k * enorm)
        D = float(np.sum(m.vol[:, None] / 4.0 * pi_vol))
        coefficient = (self.J @ x).reshape(-1, 6, 3)
        c_jump, size = [], [np.sqrt((d * d).sum(axis=-1)).max(initial=0.0)]
        for index, facet in enumerate(self.facets):
            n, t1, t2, area, tau, k_j = facet[4], facet[5], facet[6], facet[7], facet[8], facet[9]
            cb = coefficient[index]
            opening, s1, s2 = cb @ n, cb @ t1, cb @ t2
            slip = np.hypot(s1, s2)
            size.append(float(np.sqrt((cb * cb).sum(axis=1)).max()))
            if tau > 0.0:
                c_jump.append(tau * slip - opening)
                D += float(area / 6.0 * np.sum(k_j / tau * opening))
            else:
                D += float(area / 6.0 * k_j * np.sum(slip))
        c_jump = np.concatenate(c_jump) if c_jump else np.zeros(0)
        return c_vol, c_jump, D, float(max(size))

    def scale_of(self, x):
        m = self.m
        translations = x[30 * m.n_e:][np.arange(6 * m.n_f) % 6 < 3]
        v_ref = float(np.abs(np.concatenate([x[:30 * m.n_e], translations])).max())
        scale = np.full(self.n_aff, v_ref)
        for f in range(m.n_f):
            scale[self.footing_dof(f, 3):self.footing_dof(f, 6)] = v_ref / m.ell
        return scale

    def certify(self, A, g, x_c, x_int):
        scale = self.scale_of(x_c)
        value_candidate = (self.conditions(x_c)[2] - self.W0 @ x_c) / (self.W1 @ x_c)
        x_p, info = project(A, g, x_c, scale)
        move = float(np.max(np.abs(x_p - x_c) / scale))
        if move > TOL["TOL_PROJECTION_MOVE"]:
            return None, "projection_move_too_large", {"projection_move": move}
        c_vol, c_jump, _, size = self.conditions(x_p)
        c_all = np.concatenate([c_vol[np.isfinite(c_vol)], c_jump])
        rho_before = float(np.max(c_all, initial=0.0) / size)
        if rho_before > TOL["TOL_FLOW_ACCEPT"]:
            return None, "flow_violation_too_large", {"rho_before": rho_before}
        margin = TOL["TOL_CONE_MARGIN"] * size
        mu = 0.0
        x_f = x_p
        interior_margin = None
        active = c_all > margin
        if np.any(active):
            if x_int is None:
                return None, "no_admissible_correction", {"rho_before": rho_before}
            x_ip, _ = project(A, g, x_int, self.scale_of(x_int))
            ci_vol, ci_jump, _, size_int = self.conditions(x_ip)
            ci_all = np.concatenate([ci_vol[np.isfinite(ci_vol)], ci_jump])
            if not size_int > 0.0:
                return None, "interior_field_invalid", {"interior_size": size_int}
            interior_margin = float(np.max(ci_all[active]) / size_int)
            if ci_all.max() > TOL["TOL_CONE_MARGIN"] * size_int or interior_margin > -TOL["TOL_INTERIOR_MIN"]:
                return None, "interior_field_invalid", {"interior_margin": interior_margin, "active": int(active.sum())}
            mu = float(np.max(c_all[active] / (c_all[active] - ci_all[active])))
            if mu > TOL["TOL_REPAIR_MAX"]:
                return None, "repair_too_large", {"repair_mu": mu}
            x_f = (1.0 - mu) * x_p + mu * x_ip
        c_vol, c_jump, D, size = self.conditions(x_f)
        c_all = np.concatenate([c_vol[np.isfinite(c_vol)], c_jump])
        W0, W1 = float(self.W0 @ x_f), float(self.W1 @ x_f)
        w_scale = float(np.abs(self.W1) @ np.abs(x_f))
        omega = backward_error(A, g, x_f, scale)
        tau = T.tau_eq()
        kappa = info["kappa"]
        n_vol = int(np.isfinite(c_vol).sum())
        metrics = {
            "omega_eq": omega, "tau_eq": tau, "kappa_eq": kappa, "projection_move": move,
            "rho_flow": float(np.max(c_all[:n_vol], initial=0.0) / size),
            "rho_jump_cone": float(np.max(c_all[n_vol:], initial=0.0) / size), "repair_mu": mu,
            "work_reference": W1, "rho_before_repair": rho_before, "interior_margin": interior_margin,
            "repaired_conditions": int(active.sum()), "value_candidate_recomputed": float(value_candidate),
            "rank": info["rank"], "rows": info["rows"], "cols": info["cols"], "literal_zero": info["literal_zero"],
        }
        if np.any(x_f[info["zero_mask"]] != 0.0):
            return None, "literal_zero_violated", metrics
        if omega > tau:
            return None, "equality_residual_above_floor", metrics
        if kappa > TOL["TOL_KAPPA_CAP"] or kappa * omega > TOL["TOL_FLOOR_MAX"]:
            return None, "equality_ill_conditioned", metrics
        if c_all.size and c_all.max() > TOL["TOL_CONE_MARGIN"] * size:
            return None, "flow_violation_too_large", metrics
        if W1 < TOL["TOL_WORK_MIN"] * w_scale or W1 <= 0.0:
            return None, "work_nonpositive", metrics
        return (D - W0) / W1, None, metrics, x_f


# ---------------------------------------------------------------------------------------------
# Runs
# ---------------------------------------------------------------------------------------------

# (case, mesh index, or None with an explicit coarse Kuhn mesh that is not one of the case's listed meshes)
RUNS = [
    ("tresca_strip_vm_thin_slab", 0, None),
    ("upper3d_fixed_phi00", 0, None),
    ("dp_prandtl_cohesion_thin_slab", None, (8, 1, 4)),
]
STATUS = {"Solved": "solved", "AlmostSolved": "almost_solved"}


def straightness(m):
    worst = 0.0
    for tet in m.tets:
        x = m.nodes[tet]
        h = max(np.linalg.norm(x[a] - x[b]) for a in range(4) for b in range(a))
        for index, (a, b) in enumerate(T.EDGES):
            worst = max(worst, float(np.linalg.norm(x[4 + index] - 0.5 * (x[a] + x[b]))) / h)
    return worst


def documents(name, problem, arrays, m, bound, formulation, status, candidate, fields, value, reason, metrics):
    result = T._result(name, problem, bound, candidate, fields, formulation)
    result["producer"] = {"name": "certificate_demo", "version": "0"}
    result["solver"] = {"name": "clarabel", "version": clarabel.__version__, "settings": {"defaults": True}}
    result["solver_status"] = STATUS.get(status, "numerical_error")
    result.pop("result_sha256")
    result["result_sha256"] = T.document_hash(result, fields, "result_sha256")
    certificate = T._certificate(problem, arrays, result)
    certificate["checker"] = {"name": T.PROCEDURE_DEMO_CHECKER, "version": "0"}
    checks = {key: None for key in sorted(T.CHECK_KEYS)}
    checks.update({key: item for key, item in metrics.items() if key in T.CHECK_KEYS})
    checks["straightness_max"] = straightness(m)
    certificate["checks"] = checks
    if reason is None:
        kappa, omega = metrics["kappa_eq"], metrics["omega_eq"]
        certificate["value_certified"] = value
        certificate["floor"] = max(TOL["FLOOR_MIN_REL"], kappa * omega) * abs(value)
    else:
        certificate.update(classification="estimate", reasons=[reason], value_certified=None, floor=None)
    return result, certificate


def run_case(name, mesh, coarse):
    case = copy.deepcopy(T.ALL_CASES[name])
    label = "listed"
    if coarse is not None:
        case["mesh"], mesh, label = [("kuhn", coarse)], 0, "coarse-demo"
    problem, arrays = T.build_case_problem(case, mesh)
    m = Model(problem, arrays)
    out = {"case": name, "mesh": mesh, "mesh_label": label, "kuhn": list(coarse) if coarse else None, "tets": m.n_e,
           "problem_sha256": problem["problem_sha256"]}

    A, g, resultant = lower_rows(m)
    x, status = solve_lower(m, A, g)
    x = x[:A.shape[1]]
    value, reason, metrics = certify_lower(m, A, g, x, resultant)
    sigma = x[1:].reshape(m.n_e, 4, 6)
    table = np.zeros((m.n_f, 6))
    if m.n_f:
        table = (resultant @ x).reshape(m.n_f, 6)
    fields = T._float_arrays(sigma=sigma, resultant=table)
    result, certificate = documents(name, problem, arrays, m, "lower", T.LOWER_ID, status, float(x[0]), fields, value,
                                    reason, metrics)
    out["lower"] = {"result": result, "certificate": certificate,
                    "diagnostics": {k: v for k, v in metrics.items() if k not in T.CHECK_KEYS}}

    jump = Jump(m)
    A, g = jump.certificate_rows()
    x, status = jump.solve()
    interior, interior_status = None, None
    if m.alpha.max() > 0.0:
        interior, interior_status = jump.interior_field()
    outcome = jump.certify(A, g, x, interior)
    value, reason, metrics = outcome[:3]
    candidate = metrics.get("value_candidate_recomputed", float((jump.conditions(x)[2] - jump.W0 @ x) / (jump.W1 @ x)))
    fields = T._float_arrays(element_velocity=x[:30 * m.n_e].reshape(m.n_e, 10, 3),
                             footing_velocity=x[30 * m.n_e:].reshape(m.n_f, 6))
    result, certificate = documents(name, problem, arrays, m, "upper", T.JUMP_ID, status, candidate, fields, value,
                                    reason, metrics)
    out["upper"] = {"result": result, "certificate": certificate, "interior_solver_status": interior_status,
                    "diagnostics": {k: v for k, v in metrics.items() if k not in T.CHECK_KEYS}}
    return out


def main():
    arguments = sys.argv[1:]
    write = "--write" in arguments
    names = [item for item in arguments if item != "--write"]
    record = {
        "title": "PROCEDURE DEMONSTRATION, NOT A CERTIFIED BOUND: interface 17.3 and 17.4 on solver fields. "
                 "Producer and checker share one author, so every value stays an estimate for the program.",
        "date": "2026-11-19", "driver": "Ladruno_files/testbed/limit_analysis/certificate_demo.py",
        "solver": {"name": "clarabel", "version": clarabel.__version__, "settings": "defaults"},
        "numpy": np.__version__, "scipy": scipy.__version__, "runs": [],
    }
    for name, mesh, coarse in RUNS:
        if names and name not in names:
            continue
        start = time.time()
        entry = run_case(name, mesh, coarse)
        record["runs"].append(entry)
        print(name, entry["mesh_label"], "tets", entry["tets"], f"{time.time() - start:.0f}s", flush=True)
        for bound in ("lower", "upper"):
            item = entry[bound]
            doc = item["certificate"]
            shown = {k: float(f"{v:.3g}") for k, v in doc["checks"].items() if isinstance(v, float)}
            print("  ", bound, item["result"]["solver_status"], "candidate", item["result"]["value"], doc["classification"],
                  doc["value_certified"], doc["reasons"], shown, flush=True)
    if write:
        path = HERE / "certificate_demo_record.json"
        path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="ascii", newline="\n")
    return record


if __name__ == "__main__":
    main()
