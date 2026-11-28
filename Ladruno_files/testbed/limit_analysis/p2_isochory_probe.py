"""Linear-algebra probe of isochoric velocity spaces on unstructured tetrahedra (WP-237a).

THIS IS A LINEAR-ALGEBRA PROBE, NOT A BOUND. It solves no conic program and produces no collapse load.
It asks one question: can a rigid square footing move while tr(D) = 0 holds at every tetrahedron vertex?

Two velocity spaces from tests/limit_interface.md are compared on the committed meshes
tests/data/limit/upper3d_{a,b,c}.mesh.npz, with fixed walls and with roller walls:

  continuous   upper.p2c-vertex.v0, the control
  jump         upper.p2d-jump-bernstein.v0, with zero normal jump on all facet Bernstein coefficients

Reported per run: unknowns, equality rows, numerical rank, nullity, how many footing motion components
survive in the kernel, and rho_feas. rho_feas is |constraint residual|_w / |strain measure|_w, each
minimised by sparse least squares with the footing's vertical velocity prescribed. See section 11 of the
interface. A rank that is not computed is written "-": the dense factorisation is too large there.

Run from the repository root; it writes nothing unless --write is given:

    python Ladruno_files/testbed/limit_analysis/p2_isochory_probe.py --write

Needs numpy and scipy only.
"""
import argparse
import platform
import sys
from pathlib import Path

import numpy as np
import scipy
import scipy.sparse as sp
import scipy.sparse.linalg as sl

ROOT = Path(__file__).resolve().parents[3]
MESHES = ROOT / "tests" / "data" / "limit"
OUTPUT = Path(__file__).with_name("p2_isochory_probe_output.txt")

TOL_RANK = 1e-10
TOL_PLANE = 1e-12
RANK_LIMIT_CONTINUOUS = 10000
RANK_LIMIT_JUMP = 8000
EDGES = ((0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3))
FACET_VERTICES = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))
FACET_NODES = ((1, 2, 3, 5, 9, 8), (0, 3, 2, 7, 9, 6), (0, 1, 3, 4, 8, 7), (0, 2, 1, 6, 5, 4))
BERNSTEIN = np.array([
    [1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0],
    [-0.5, -0.5, 0, 2, 0, 0], [0, -0.5, -0.5, 0, 2, 0], [-0.5, 0, -0.5, 0, 0, 2],
])
SHEAR_WEIGHT = np.array([1.0, 1.0, 1.0, 2.0 ** -0.5, 2.0 ** -0.5, 2.0 ** -0.5])


def voigt_n(n):
    return np.array([
        [n[0], 0, 0, n[1], 0, n[2]],
        [0, n[1], 0, n[0], n[2], 0],
        [0, 0, n[2], 0, n[1], n[0]],
    ])


def load_mesh(name):
    with np.load(MESHES / f"upper3d_{name}.mesh.npz", allow_pickle=False) as data:
        vertices, tets = data["nodes"], data["tets"]
    edges = sorted({tuple(sorted((int(t[a]), int(t[b])))) for t in tets for a, b in EDGES})
    index = {edge: len(vertices) + i for i, edge in enumerate(edges)}
    mids = np.array([0.5 * (vertices[a] + vertices[b]) for a, b in edges])
    rows = [[int(v) for v in t] + [index[tuple(sorted((int(t[a]), int(t[b]))))] for a, b in EDGES] for t in tets]
    return np.vstack([vertices, mids]), np.array(rows)


class Mesh:
    def __init__(self, name):
        self.name = name
        self.x, self.t = load_mesh(name)
        self.n_e = len(self.t)
        self.grad, self.vol = [], []
        for tet in self.t:
            p = self.x[tet[:4]]
            jac = (p[1:] - p[0]).T
            inverse = np.linalg.inv(jac)
            self.grad.append(np.vstack([-inverse.sum(axis=0), inverse]))
            self.vol.append(np.linalg.det(jac) / 6.0)
        owners = {}
        for e, tet in enumerate(self.t):
            for f, local in enumerate(FACET_VERTICES):
                owners.setdefault(tuple(sorted(int(tet[a]) for a in local)), []).append((e, f))
        self.interior = sorted(tuple(sorted(v)) for v in owners.values() if len(v) == 2)
        self.exterior = sorted(v[0] for v in owners.values() if len(v) == 1)

    def facet(self, e, f):
        p = self.x[self.t[e][list(FACET_VERTICES[f])]]
        cross = np.cross(p[1] - p[0], p[2] - p[0])
        area = 0.5 * np.linalg.norm(cross)
        return p, cross / (2.0 * area), area

    def vertex_gradients(self, e, m):
        g = self.grad[e]
        out = np.zeros((10, 3))
        for i in range(4):
            out[i] = (3.0 if i == m else -1.0) * g[i]
        for q, (a, b) in enumerate(EDGES):
            if a == m:
                out[4 + q] = 4.0 * g[b]
            elif b == m:
                out[4 + q] = 4.0 * g[a]
        return out


def classify(p, walls):
    """Boundary kind of one exterior facet; walls are every face except the ground z = 0."""
    if np.all(np.abs(p[:, 2]) < TOL_PLANE):
        inside = np.all(np.abs(p[:, 0]) <= 0.5 + TOL_PLANE) and np.all(np.abs(p[:, 1]) <= 0.5 + TOL_PLANE)
        return "footing" if inside else "free"
    return walls


class Space:
    """Velocity unknowns after removing fixed, roller and rigid-footing degrees of freedom exactly."""

    def __init__(self, mesh, walls, broken):
        self.mesh, self.broken = mesh, broken
        self.state = {}
        for e, f in mesh.exterior:
            p, normal, _ = mesh.facet(e, f)
            kind = classify(p, walls)
            for a in FACET_NODES[f]:
                key = (e, a) if broken else int(mesh.t[e][a])
                entry = self.state.setdefault(key, {"zero": set(), "footing": False})
                if kind == "fixed":
                    entry["zero"].update((0, 1, 2))
                elif kind == "symmetry":
                    entry["zero"].add(int(np.argmax(np.abs(normal))))
                elif kind == "footing":
                    entry["footing"] = True
        keys = ((e, a) for e in range(mesh.n_e) for a in range(10)) if broken else range(len(mesh.x))
        self.columns, count = {}, 0
        for key in keys:
            entry = self.state.get(key)
            if entry is not None and entry["footing"]:
                continue
            zero = entry["zero"] if entry else set()
            for c in range(3):
                if c not in zero:
                    self.columns[(key, c)] = count
                    count += 1
        self.n_free = count
        self.n = count + 6  # footing motion is appended: u, then w

    def terms(self, e, a, c):
        key = (e, a) if self.broken else int(self.mesh.t[e][a])
        entry = self.state.get(key)
        if entry is not None and entry["footing"]:
            r = self.mesh.x[self.mesh.t[e][a]]
            rotation = np.array([[0, r[2], -r[1]], [-r[2], 0, r[0]], [r[1], -r[0], 0]])
            out = [(self.n_free + c, 1.0)]
            out += [(self.n_free + 3 + q, rotation[c, q]) for q in range(3) if rotation[c, q] != 0.0]
            return out
        column = self.columns.get((key, c))
        return [] if column is None else [(column, 1.0)]


class Rows:
    def __init__(self):
        self.i, self.j, self.v, self.count = [], [], [], 0

    def add(self, terms):
        for column, value in terms:
            self.i.append(self.count)
            self.j.append(column)
            self.v.append(value)
        self.count += 1

    def matrix(self, n):
        return sp.csr_matrix((self.v, (self.i, self.j)), shape=(self.count, n))


def build(space):
    """Return constraint rows C and strain-measure rows S, both already carrying sqrt(weights)."""
    mesh, constraint, strain = space.mesh, Rows(), Rows()
    for e in range(mesh.n_e):
        weight = np.sqrt(mesh.vol[e] / 4.0)
        for m in range(4):
            gradients = mesh.vertex_gradients(e, m)
            theta = []
            for a in range(10):
                for c in range(3):
                    if gradients[a, c] != 0.0:
                        theta += [(j, weight * gradients[a, c] * v) for j, v in space.terms(e, a, c)]
            constraint.add(theta)
            for q in range(6):
                row = []
                for a in range(10):
                    operator = voigt_n(gradients[a]).T[q]
                    for c in range(3):
                        if operator[c] != 0.0:
                            row += [(j, weight * SHEAR_WEIGHT[q] * operator[c] * v) for j, v in space.terms(e, a, c)]
                strain.add(row)
    if space.broken:
        for (e1, f1), (e2, _) in mesh.interior:
            p, normal, area = mesh.facet(e1, f1)
            size = max(np.linalg.norm(p[i] - p[(i + 1) % 3]) for i in range(3))
            weight = np.sqrt(area / 6.0 / size)
            nodes = [int(mesh.t[e1][a]) for a in FACET_NODES[f1]]
            local1 = list(FACET_NODES[f1])
            order2 = [int(v) for v in mesh.t[e2]]
            local2 = [order2.index(node) for node in nodes]
            axis = int(np.argmin(np.abs(normal)))
            helper = np.zeros(3)
            helper[axis] = 1.0
            t1 = np.cross(helper, normal)
            t1 /= np.linalg.norm(t1)
            t2 = np.cross(normal, t1)
            for coefficient in BERNSTEIN:
                for direction, targets in ((normal, (constraint, strain)), (t1, (strain,)), (t2, (strain,))):
                    row = []
                    for q in range(6):
                        if coefficient[q] == 0.0:
                            continue
                        for c in range(3):
                            scale = weight * coefficient[q] * direction[c]
                            if scale != 0.0:
                                row += [(j, scale * v) for j, v in space.terms(e2, local2[q], c)]
                                row += [(j, -scale * v) for j, v in space.terms(e1, local1[q], c)]
                    for target in targets:
                        target.add(row)
    return constraint.matrix(space.n), strain.matrix(space.n)


def residual_with_prescribed_vz(matrix, n_free):
    """min |A u| with footing velocity u_z = -1 (the unit downward motion)."""
    column = n_free + 2
    keep = np.array([j for j in range(matrix.shape[1]) if j != column])
    right = np.asarray(matrix[:, column].todense()).ravel()
    solution = sl.lsmr(matrix[:, keep], right, atol=1e-14, btol=1e-14, conlim=0.0, maxiter=400000)[0]
    return float(np.linalg.norm(matrix[:, keep] @ solution - right))


def rank_and_mobility(matrix, n_free, need_kernel):
    dense = matrix.toarray()
    if not need_kernel:
        s = np.linalg.svd(dense, compute_uv=False)
        return int(np.sum(s > TOL_RANK * s[0])), None
    _, s, vt = np.linalg.svd(dense, full_matrices=True)
    rank = int(np.sum(s > TOL_RANK * s[0]))
    kernel = vt[rank:]
    if kernel.size == 0:
        return rank, 0
    part = kernel[:, n_free:]
    singular = np.linalg.svd(part, compute_uv=False)
    return rank, int(np.sum(singular > 1e-8))


def run_case(mesh, walls, broken):
    space = Space(mesh, walls, broken)
    constraint, strain = build(space)
    numerator = residual_with_prescribed_vz(constraint, space.n_free)
    denominator = residual_with_prescribed_vz(strain, space.n_free)
    rank, mobile = None, None
    if not broken and constraint.shape[1] <= RANK_LIMIT_CONTINUOUS:
        rank, mobile = rank_and_mobility(constraint, space.n_free, True)
    elif broken and constraint.shape[0] <= RANK_LIMIT_JUMP:
        rank, mobile = rank_and_mobility(constraint, space.n_free, False)
    return {
        "unknowns": space.n, "rows": constraint.shape[0], "rank": rank, "mobile": mobile,
        "rho": numerator / denominator,
    }


def text(value):
    return "-" if value is None else str(value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="write the output file next to this script")
    args = parser.parse_args()
    lines = [
        "LINEAR-ALGEBRA PROBE, NOT A BOUND. No conic program was solved and no collapse load is implied.",
        "Question: can the rigid footing move with tr(D) = 0 at every tet vertex (and, for 'jump', zero normal",
        "jump on every facet Bernstein coefficient)? rho_feas near 0 means yes; a value far above 0 means no.",
        f"numpy {np.__version__}, scipy {scipy.__version__}, python {platform.python_version()}",
        f"rank threshold {TOL_RANK:g} of the largest singular value; '-' means not computed at this size",
        "",
        "mesh tets walls    space        unknowns   rows   rank nullity footing_modes rho_feas",
    ]
    for name in "abc":
        mesh = Mesh(name)
        for walls in ("fixed", "symmetry"):
            for broken, label in ((False, "continuous"), (True, "jump")):
                out = run_case(mesh, walls, broken)
                nullity = None if out["rank"] is None else out["unknowns"] - out["rank"]
                lines.append(
                    f"{name:4s} {mesh.n_e:4d} {'roller' if walls == 'symmetry' else 'fixed':8s} {label:12s}"
                    f"{out['unknowns']:8d} {out['rows']:6d} {text(out['rank']):>6s} {text(nullity):>7s}"
                    f" {text(out['mobile']):>13s} {out['rho']:8.1e}"
                )
                print(lines[-1], flush=True)
    report = "\n".join(lines) + "\n"
    if args.write:
        OUTPUT.write_text(report, encoding="ascii", newline="\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
