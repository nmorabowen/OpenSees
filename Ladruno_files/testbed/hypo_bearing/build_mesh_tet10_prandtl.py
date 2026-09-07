"""ADR-95 P3: a structured tet10 mesh of the SAME plane-strain Prandtl-Reissner
strip domain as `h20_prandtl.strip_mesh` (NOT the ADR-79 square-footing box
`build_mesh_tet10.py` builds -- `bearing_mesh_tet10.npz` is a different problem:
a 10x10x8 m 3-D box under a 2x2 m SQUARE footing with PDMY sand, no plane-strain
slab.  Verified by inspection before writing this file: its node coordinate
ranges are x,y in [-10,10], z in [-8,0], footing 49 nodes -- the ADR-79 mesh,
not usable here without silently changing the problem (a square footing has no
q0*Nq oracle).  So this script builds a NEW mesh, one element thick in y (the
plane-strain slab, thickness = h20_prandtl.THICK) and graded in x/z exactly on
h20_prandtl.strip_mesh's own block boundaries (XLIM=30, ZBOT=-20, B_FOOT=2,
h0=1.0 -> 9 graded + 2 uniform + 9 graded in x, 7 graded + 3 uniform in z --
200 hex-shaped cells, matching the H20/H8 h0=1.0 legs' cell count), each split
into 6 structured tets by gmsh (transfinite volume, NO recombine).

Node order convention: TenNodeTetrahedron order (v1 v2 v3 v4, e12 e23 e13 e14
e34 e24), which BezierTet10 was built to match (both source files say so) --
gmsh TET10 (type 11) edge slots 4..9 are (0,1)(1,2)(0,2)(0,3)(2,3)(1,3), the
SAME order, so gmsh connectivity maps with NO permutation, verified below by
the same midpoint assertion `build_mesh_tet10.py` uses.

Run:  C:\\Users\\nmb\\venv\\opensees_env\\Scripts\\python.exe build_mesh_tet10_prandtl.py
Writes bearing_mesh_tet10_prandtl.npz: nodes(N,3), tets(M,10) int64,
is_vertex(N,) bool, sets: top/bottom/xface/footing (node indices, 0-based),
h0, b_foot, thick, xlim, zbot -- the deck parameters this mesh was built at.
"""
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "bearing_mesh_tet10_prandtl.npz")

# --- identical to h20_prandtl.py's problem constants ------------------------
B_FOOT = 2.0
THICK = 0.5
XLIM, ZBOT = 30.0, -20.0
R_GRADE = 1.35
H0 = 1.0                         # 2 quadratic-tet elements across B (matches
                                  # the H20/H8 h0=1.0 legs' cell count, 200)
TOL = 1e-6

GMSH_TET10_EDGES = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]


def n_graded(length, h0, r=R_GRADE):
    return max(1, int(np.ceil(np.log(1.0 + length * (r - 1.0) / h0) / np.log(r))))


def growth(length, n, h_first):
    if abs(h_first * n - length) < 1e-12:
        return 1.0

    def total(r):
        return (h_first * n if abs(r - 1.0) < 1e-12
                else h_first * (r ** n - 1.0) / (r - 1.0))

    lo, hi = 1.0, 4.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if total(mid) < length:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def spec(blocks):
    out = {}
    for lo, hi, n, fine in blocks:
        r = 1.0 if fine is None else growth(hi - lo, n, H0)
        out[(lo, hi)] = (n + 1, r, fine)
    return out


def main():
    from apeGmsh import apeGmsh
    import gmsh

    no = n_graded(XLIM - 1.0, H0)          # outer x block cell count
    nf = int(round(B_FOOT / H0))           # inner x block (under+beside footing)
    nz_out = n_graded(-ZBOT - 3.0, H0)     # outer z block
    nz_in = int(round(3.0 / H0))           # inner z block
    print(f"[mesh] blocks: x outer {no} (x2) + inner {nf} = "
          f"{2 * no + nf}; z outer {nz_out} + inner {nz_in} = {nz_out + nz_in}; "
          f"y 1 (thickness {THICK}) -> {(2 * no + nf) * (nz_out + nz_in)} "
          f"hex-shaped cells")

    XBLOCKS = [(-XLIM, -1.0, no, -1.0), (-1.0, 1.0, nf, None), (1.0, XLIM, no, 1.0)]
    YBLOCKS = [(0.0, THICK, 1, None)]
    ZBLOCKS = [(ZBOT, -3.0, nz_out, -3.0), (-3.0, 0.0, nz_in, None)]

    sx, sy, sz = spec(XBLOCKS), spec(YBLOCKS), spec(ZBLOCKS)

    with apeGmsh(model_name="adr95_prandtl_tet10") as g:
        for i, (x0, x1, *_) in enumerate(XBLOCKS):
            for j, (y0, y1, *_) in enumerate(YBLOCKS):
                for k, (z0, z1, *_) in enumerate(ZBLOCKS):
                    g.model.geometry.add_box(x0, y0, z0, x1 - x0, y1 - y0,
                                             z1 - z0, label=f"b{i}{j}{k}")
        vols = [t for (d, t) in gmsh.model.getEntities(3)]
        g.model.boolean.fragment([(3, vols[0])], [(3, t) for t in vols[1:]])
        gmsh.model.occ.synchronize()

        for _, ct in gmsh.model.getEntities(1):
            t0, t1 = gmsh.model.getParametrizationBounds(1, ct)
            p0 = np.array(gmsh.model.getValue(1, ct, [t0[0]]))
            p1 = np.array(gmsh.model.getValue(1, ct, [t1[0]]))
            d = p1 - p0
            axis = int(np.argmax(np.abs(d)))
            if np.abs(d).sum() - abs(d[axis]) > TOL:
                raise RuntimeError(f"curve {ct} is not axis-aligned: {d}")
            s = (sx, sy, sz)[axis]
            lo, hi = min(p0[axis], p1[axis]), max(p0[axis], p1[axis])
            key = next((k for k in s
                        if abs(k[0] - lo) < TOL and abs(k[1] - hi) < TOL), None)
            if key is None:
                raise RuntimeError(f"curve {ct} spans {lo}..{hi} axis {axis}")
            nn, r, fine = s[key]
            coef = (1.0 if r == 1.0
                    else (r if abs(p0[axis] - fine) < TOL else 1.0 / r))
            g.mesh.structured.set_transfinite_curve(ct, nn, coef=coef)

        for _, st in gmsh.model.getEntities(2):
            g.mesh.structured.set_transfinite_surface(st)
        for _, vt in gmsh.model.getEntities(3):
            g.mesh.structured.set_transfinite_volume(vt)

        g.mesh.generation.generate(dim=3)
        gmsh.model.mesh.setOrder(2)

        ntags, coords = gmsh.model.mesh.getNodes()[:2]
        coords = np.asarray(coords).reshape(-1, 3)
        order = np.argsort(ntags)
        gid = np.asarray(ntags)[order]
        xyz = coords[order]
        remap = {int(t): i for i, t in enumerate(gid)}

        tets = []
        for et, _, enodes in zip(*gmsh.model.mesh.getElements(dim=3)):
            if et != 11:
                raise RuntimeError(f"expected tet10 (type 11), got type {et}")
            conn = np.asarray(enodes).reshape(-1, 10)
            tets.append(np.vectorize(remap.get)(conn))
        tets = np.vstack(tets)

    # ---- straight-side + edge-order verification ---------------------------
    bad = 0
    for a, b in enumerate(GMSH_TET10_EDGES):
        mid_xyz = xyz[tets[:, 4 + a]]
        want = 0.5 * (xyz[tets[:, b[0]]] + xyz[tets[:, b[1]]])
        err = np.abs(mid_xyz - want).max()
        if err > 1e-9:
            bad += 1
            print(f"[mesh] EDGE SLOT {4 + a} does not match midpoint of "
                  f"{b}: max err {err:.3e}")
    if bad:
        raise RuntimeError("gmsh tet10 edge ordering does not match "
                           "TenNodeTetrahedron/BezierTet10's convention, or "
                           "sides are curved")
    print(f"[mesh] tet10 edge slots match {GMSH_TET10_EDGES} to 1e-9 -- "
          f"straight-sided, TenNodeTetrahedron/BezierTet10 ordering confirmed")

    # ---- positive Jacobian on the vertex tet --------------------------------
    p = xyz[tets[:, :4]]
    vol6 = np.einsum("ij,ij->i", p[:, 3] - p[:, 0],
                     np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]))
    flip = vol6 < 0
    if flip.any():
        print(f"[mesh] flipping {flip.sum()} inverted tets")
        tets[flip][:, [1, 2]] = tets[flip][:, [2, 1]]
        t = tets[flip]
        t[:, [4, 6]] = t[:, [6, 4]]
        t[:, [8, 9]] = t[:, [9, 8]]
        tets[flip] = t
        p = xyz[tets[:, :4]]
        vol6 = np.einsum("ij,ij->i", p[:, 3] - p[:, 0],
                         np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]))
        assert (vol6 > 0).all(), "tet flip failed"
    vol = vol6.sum() / 6.0
    exact = 2 * XLIM * THICK * (-ZBOT)
    print(f"[mesh] volume {vol:.6f} m3 against exact {exact:.6f} "
          f"(rel {abs(vol - exact) / exact:.2e})")
    assert abs(vol - exact) / exact < 1e-8, "mesh volume does not match the domain"

    is_vertex = np.zeros(len(xyz), dtype=bool)
    is_vertex[np.unique(tets[:, :4])] = True
    nmid = int((~is_vertex).sum())
    print(f"[mesh] {len(xyz)} nodes ({int(is_vertex.sum())} vertex + {nmid} "
          f"mid-edge), {len(tets)} tet10, {3 * len(xyz)} DOF (ndf=3)")

    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    sets = {
        "top": np.where(np.isclose(z, 0.0, atol=TOL))[0],
        "bottom": np.where(np.isclose(z, ZBOT, atol=TOL))[0],
        "xface": np.where(np.isclose(np.abs(x), XLIM, atol=TOL))[0],
    }
    half = 0.5 * B_FOOT
    sets["footing"] = np.where(np.isclose(z, 0.0, atol=TOL)
                               & (np.abs(x) <= half + TOL))[0]
    for k, v in sets.items():
        print(f"[mesh] set {k}: {len(v)} nodes")

    np.savez_compressed(
        OUT, nodes=xyz, tets=tets, is_vertex=is_vertex,
        h0=H0, b_foot=B_FOOT, thick=THICK, xlim=XLIM, zbot=ZBOT,
        **{f"set_{k}": v for k, v in sets.items()})
    print(f"[mesh] -> {OUT}")


if __name__ == "__main__":
    main()
