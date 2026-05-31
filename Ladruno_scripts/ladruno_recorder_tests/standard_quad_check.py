"""Standard-rule QUADRATURE + GLOBAL_GP_COORDS gate -- checker (Step B).

Self-contained h5py inspector (does NOT use ladruno_format, so it does not depend
on the Step C reader changes). For each sq_*.ladruno produced by
standard_quad_model.py it:

  1. prints TOPOLOGY/FAMILY/ORDER/NDIR/NUM_GP + GP_PARAM/GP_WEIGHT/GLOBAL_GP_COORDS;
  2. asserts NDIR/NUM_GP/GP_PARAM match the verified per-element ordering tables
     (quad CCW, tri centroid bary, brick i,j,k lexicographic);
  3. WRITE-TIME ROUND-TRIP ORACLE -- reconstructs x(GP_PARAM[k]) = Sum N_i(xi_k) X_i
     from the file's own MODEL/NODES + CONNECTIVITY using an INDEPENDENT Python
     reimplementation of the element shape functions, and compares to
     GLOBAL_GP_COORDS[k] (<= 1e-12). This catches any GP-ordering / basis bug.

Run with the venv python (has h5py/numpy):
    python standard_quad_check.py <out_dir>
"""

from __future__ import annotations

import os
import sys

import h5py
import numpy as np

TOL = 1e-12
A = 1.0 / np.sqrt(3.0)            # 2-pt GL abscissa
# expected GP_PARAM tables (must match getStandardQuadrature, element GP order)
EXPECT_GP = {
    "quad": np.array([[-A, -A], [A, -A], [A, A], [-A, A]]),                 # CCW
    "tri": np.array([[1.0 / 3.0, 1.0 / 3.0]]),                              # centroid
    "hex": np.array([[-A, -A, -A], [-A, -A, A], [-A, A, -A], [-A, A, A],    # i,j,k
                     [A, -A, -A], [A, -A, A], [A, A, -A], [A, A, A]]),
}


# --- independent shape-function reimplementation (mirrors computeGlobalGP) --- #
def shape(topology: str, num_nodes: int, xi: np.ndarray) -> np.ndarray:
    if topology == "line" and num_nodes >= 2:
        s = xi[0]
        return np.array([0.5 * (1 - s), 0.5 * (1 + s)])
    if topology == "tri" and num_nodes >= 3:
        x, e = xi[0], xi[1]
        return np.array([x, e, 1 - x - e])
    if topology == "quad" and num_nodes == 4:
        s, t = xi[0], xi[1]
        return 0.25 * np.array([(1 - s) * (1 - t), (1 + s) * (1 - t),
                                (1 + s) * (1 + t), (1 - s) * (1 + t)])
    if topology == "quad" and num_nodes == 9:
        s, t = xi[0], xi[1]
        return np.array([
            (1 - s) * (1 - t) * s * t / 4, -(1 + s) * (1 - t) * s * t / 4,
            (1 + s) * (1 + t) * s * t / 4, -(1 - s) * (1 + t) * s * t / 4,
            -(1 - s * s) * (1 - t) * t / 2, (1 + s) * (1 - t * t) * s / 2,
            (1 - s * s) * (1 + t) * t / 2, -(1 - s) * (1 - t * t) * s / 2,
            (1 - s * s) * (1 - t * t)])
    if topology == "tet" and num_nodes >= 4:
        r, s, u = xi[0], xi[1], xi[2]
        return np.array([r, s, u, 1 - r - s - u])
    if topology == "hex" and num_nodes == 8:
        r, s, u = xi[0], xi[1], xi[2]
        sgn = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                        [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
        return 0.125 * np.prod(1 + sgn * np.array([r, s, u]), axis=1)
    raise ValueError(f"no Python basis for topology={topology} num_nodes={num_nodes}")


def _attr(obj, key):
    v = obj.attrs[key]
    a = np.atleast_1d(v)
    return a[0] if a.size == 1 else a


def check_file(path: str, expect_topo: str) -> int:
    problems = 0
    print(f"\n=== {os.path.basename(path)} (expect topology={expect_topo}) ===")
    with h5py.File(path, "r") as f:
        stages = [k for k in f if k.startswith("MODEL_STAGE")]
        stage = stages[0]
        nid = f[f"{stage}/MODEL/NODES/ID"][...].ravel()
        ncoord = f[f"{stage}/MODEL/NODES/COORDINATES"][...]
        ndim = ncoord.shape[1]
        coord_of = {int(t): ncoord[i] for i, t in enumerate(nid)}

        base = f[f"{stage}/MODEL/ELEMENTS"]
        for name in base:
            grp = base[name]
            topo = _attr(grp, "TOPOLOGY")
            topo = topo.decode() if isinstance(topo, bytes) else str(topo)
            fam = _attr(grp, "FAMILY")
            fam = fam.decode() if isinstance(fam, bytes) else str(fam)
            order = tuple(int(x) for x in np.atleast_1d(_attr(grp, "ORDER")))
            has_ndir = "NDIR" in grp.attrs
            ndir = int(_attr(grp, "NDIR")) if has_ndir else None
            num_gp = int(_attr(grp, "NUM_GP")) if "NUM_GP" in grp.attrs else None
            conn = grp["CONNECTIVITY"][...]
            print(f"  group {name}: TOPOLOGY={topo} FAMILY={fam} ORDER={order} "
                  f"NDIR={ndir} NUM_GP={num_gp} nElem={conn.shape[0]}")

            if "QUADRATURE" not in grp:
                print("    [skip] no QUADRATURE group")
                continue
            gp_param = grp["QUADRATURE/GP_PARAM"][...]
            gp_weight = grp["QUADRATURE/GP_WEIGHT"][...] if "QUADRATURE/GP_WEIGHT" in grp else None
            print(f"    GP_PARAM shape={gp_param.shape}\n{gp_param}")
            print(f"    GP_WEIGHT={None if gp_weight is None else gp_weight.ravel()}")

            # (a) NDIR present + equals GP_PARAM cols
            if not has_ndir:
                print("    [FAIL] NDIR attribute missing"); problems += 1
            elif ndir != gp_param.shape[1]:
                print(f"    [FAIL] NDIR={ndir} != GP_PARAM cols {gp_param.shape[1]}"); problems += 1
            if num_gp != gp_param.shape[0]:
                print(f"    [FAIL] NUM_GP={num_gp} != GP_PARAM rows {gp_param.shape[0]}"); problems += 1

            # (b) GP_PARAM matches the verified ordering table
            if topo in EXPECT_GP:
                exp = EXPECT_GP[topo]
                if gp_param.shape != exp.shape or not np.allclose(gp_param, exp, atol=1e-9):
                    print(f"    [FAIL] GP_PARAM != expected {topo} ordering\n      expected\n{exp}")
                    problems += 1
                else:
                    print(f"    [OK] GP_PARAM matches verified {topo} ordering")

            # (c) GLOBAL_GP_COORDS landed + round-trip oracle
            if "GLOBAL_GP_COORDS" not in grp:
                print("    [FAIL] GLOBAL_GP_COORDS missing"); problems += 1
                continue
            ggp = grp["GLOBAL_GP_COORDS"][...]
            n_elem = conn.shape[0]
            if ggp.shape != (n_elem, num_gp * ndim):
                print(f"    [FAIL] GLOBAL_GP_COORDS shape {ggp.shape} "
                      f"!= ({n_elem}, {num_gp * ndim})"); problems += 1
                continue
            print(f"    GLOBAL_GP_COORDS shape={ggp.shape} (=[nElem x nGP*ndim])")
            ggp3 = ggp.reshape(n_elem, num_gp, ndim)

            max_err = 0.0
            num_nodes = conn.shape[1] - 1
            for e in range(n_elem):
                node_tags = conn[e, 1:]
                X = np.array([coord_of[int(t)][:ndim] for t in node_tags])
                for k in range(num_gp):
                    N = shape(topo, num_nodes, gp_param[k])
                    x_rec = N @ X
                    err = np.max(np.abs(x_rec - ggp3[e, k]))
                    max_err = max(max_err, err)
            status = "OK" if max_err <= TOL else "FAIL"
            if status == "FAIL":
                problems += 1
            print(f"    [{status}] round-trip oracle x(GP_PARAM) vs GLOBAL_GP_COORDS "
                  f"max_err={max_err:.3e} (tol {TOL:.0e})")
            print(f"    GLOBAL_GP_COORDS[0]=\n{ggp3[0]}")
    return problems


def main() -> int:
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    total = 0
    total += check_file(os.path.join(out, "sq_quad.ladruno"), "quad")
    total += check_file(os.path.join(out, "sq_tri.ladruno"), "tri")
    total += check_file(os.path.join(out, "sq_brick.ladruno"), "hex")
    print("\n" + ("=" * 60))
    if total == 0:
        print("STANDARD_QUAD_CHECK: ALL PASS")
        return 0
    print(f"STANDARD_QUAD_CHECK: {total} PROBLEM(S)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
