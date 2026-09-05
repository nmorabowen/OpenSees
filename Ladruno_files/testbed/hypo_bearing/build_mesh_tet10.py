"""Graded Bezier-TET10 mesh of the SAME bearing domain as build_mesh.py, for
the locking-free coupled legs.

WHY A SECOND MESH. `-formulation bbar` on a LINEAR element buys one of the two
ingredients note 15 of the vault identifies; reaching the unlocked floor needs
QUADRATIC interpolation as well. `LadrunoUP` has exactly one quadratic 3D
provider — the 10-node Bezier tet, Taylor-Hood (quadratic u, linear vertex p) —
so an unlocked coupled run needs a tet mesh. `-geom hypo` is H8-only and plays
no part here; the lane is `-geom corot`.

WHY COARSER CELLS. Taylor-Hood tet10 is DOF-expensive: vertices carry 4 DOF and
every mid-edge node carries 3, and a structured hex cell splits into 6 tets with
roughly 7 new edges. Splitting the H8 mesh 1:1 would cost ~73k DOF against its
13.9k — 5x, on elements that are also dearer per DOF. So the mesh is matched on
DOF rather than on cell count, exactly as note 15 matched its tet meshes, with
COARSER cells carrying QUADRATIC interpolation instead of finer linear ones.

WHAT THAT COSTS, STATED UP FRONT. The H8 mesh puts 4 linear elements across the
footing; this one puts 3 quadratic. Absolute capacities are therefore NOT
comparable with the H8 legs. What these legs measure is the bbar-vs-std
DIFFERENCE on a quadratic element, and both legs of a pair share this mesh, so
the ratio is clean even where the absolute value is not.

Run:  C:\\Users\\nmb\\venv\\opensees_env\\Scripts\\python.exe build_mesh_tet10.py
Writes bearing_mesh_tet10.npz (nodes, tets[:,10], is_vertex, tributary, sets).
"""
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "bearing_mesh_tet10.npz")

B_FOOT = 2.0
XLIM, ZBOT = 10.0, -8.0
H0 = 2.0 / 3.0                 # 3 quadratic elements across the footing

# (lo, hi, n_elements, coordinate at which the elements are FINE; None=uniform)
XBLOCKS = [(-XLIM, -1.0, 3, -1.0), (-1.0, 1.0, 3, None), (1.0, XLIM, 3, 1.0)]
YBLOCKS = XBLOCKS
ZBLOCKS = [(ZBOT, -3.0, 3, -3.0), (-3.0, 0.0, 4, None)]
TOL = 1e-6

# gmsh TET10 (type 11) edge slots 4..9 are (0,1),(1,2),(0,2),(0,3),(2,3),(1,3)
# — the SAME order LadrunoUP expects for its 6 mid-edge tags (verified against
# tests/test_ladruno_up_element_th.py::_btet10_single). So gmsh connectivity
# maps onto (vt + mt) with no permutation. Asserted below, not assumed.
GMSH_TET10_EDGES = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]


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

    sx, sy, sz = spec(XBLOCKS), spec(YBLOCKS), spec(ZBLOCKS)

    with apeGmsh(model_name="adr79_bearing_tet10") as g:
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

        # transfinite surfaces/volumes but NO recombine -> structured TETS
        for _, st in gmsh.model.getEntities(2):
            g.mesh.structured.set_transfinite_surface(st)
        for _, vt in gmsh.model.getEntities(3):
            g.mesh.structured.set_transfinite_volume(vt)

        g.mesh.generation.generate(dim=3)
        gmsh.model.mesh.setOrder(2)             # tet4 -> tet10
        # Straight sides: LadrunoUP's Bezier tet10 guards against curved edges,
        # and every boundary here is planar, so setOrder must have placed each
        # mid node at its edge midpoint. Verified below.

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

    # ---- straight-side + edge-order verification (both are assertions) ------
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
        raise RuntimeError("gmsh tet10 edge ordering does not match LadrunoUP's "
                           "(vt+mt) convention, or sides are curved")
    print(f"[mesh] tet10 edge slots match {GMSH_TET10_EDGES} to 1e-9 — "
          f"straight-sided, LadrunoUP ordering confirmed")

    # ---- positive Jacobian on the vertex tet --------------------------------
    p = xyz[tets[:, :4]]
    vol6 = np.einsum("ij,ij->i", p[:, 3] - p[:, 0],
                     np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]))
    flip = vol6 < 0
    if flip.any():
        print(f"[mesh] flipping {flip.sum()} inverted tets")
        # swap vertices 1<->2 and the edge slots that follow them
        tets[flip][:, [1, 2]] = tets[flip][:, [2, 1]]
        t = tets[flip]
        t[:, [4, 6]] = t[:, [6, 4]]      # (0,1)<->(0,2)
        t[:, [8, 9]] = t[:, [9, 8]]      # (2,3)<->(1,3)
        tets[flip] = t
        p = xyz[tets[:, :4]]
        vol6 = np.einsum("ij,ij->i", p[:, 3] - p[:, 0],
                         np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]))
        assert (vol6 > 0).all(), "tet flip failed"
    vol = vol6.sum() / 6.0
    exact = (2 * XLIM) ** 2 * (-ZBOT)
    print(f"[mesh] volume {vol:.6f} m3 against exact {exact:.6f} "
          f"(rel {abs(vol - exact) / exact:.2e})")

    # ---- vertex vs mid-edge nodes (ndf 4 vs ndf 3) --------------------------
    is_vertex = np.zeros(len(xyz), dtype=bool)
    is_vertex[np.unique(tets[:, :4])] = True
    nmid = int((~is_vertex).sum())
    ndof = 4 * int(is_vertex.sum()) + 3 * nmid
    print(f"[mesh] {len(xyz)} nodes ({int(is_vertex.sum())} vertex + {nmid} "
          f"mid-edge), {len(tets)} tet10, {ndof} u-p DOF")

    # ---- node sets ----------------------------------------------------------
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    sets = {
        "xface": np.where(np.isclose(np.abs(x), XLIM, atol=TOL))[0],
        "yface": np.where(np.isclose(np.abs(y), XLIM, atol=TOL))[0],
        "bottom": np.where(np.isclose(z, ZBOT, atol=TOL))[0],
        "top": np.where(np.isclose(z, 0.0, atol=TOL))[0],
    }
    half = 0.5 * B_FOOT
    sets["footing"] = np.where(np.isclose(z, 0.0, atol=TOL)
                               & (np.abs(x) <= half + TOL)
                               & (np.abs(y) <= half + TOL))[0]

    # ---- surface tributary: LUMPED, deliberately, not consistent ------------
    # The CONSISTENT T6 rule for uniform pressure is zero at the vertices and
    # A/3 at each midside node, and it is computed below as `tributary_consist`.
    # It is NOT what the runner loads with, because it leaves all 100 surface
    # VERTEX nodes carrying exactly zero — and scoping finding 2 records that
    # this surcharge exists to REGULARISE the zero-confinement surface DOFs as
    # much as to be a physical overburden. With the consistent rule the tet10
    # legs cannot take their first plastic step at ANY increment: subdividing to
    # dt = 0.06 s (13 halvings) still leaves a force residual of order 250,
    # which is finding 2's failure and not a step-size one. The vertices are
    # also the only nodes carrying a pressure DOF, so leaving them unloaded is
    # doubly bad here.
    #
    # So the loaded weights are LUMPED — area/6 to each of a face's six nodes,
    # which sums to the exact area and leaves no surface DOF unloaded. The cost
    # is that a uniform pressure is no longer reproduced exactly in the weak
    # sense; against a 10 kPa surcharge on a 600-3300 kPa signal that is the
    # sub-one-percent bias the campaign already records for the H8 legs.
    trib_consist = np.zeros(len(xyz))
    trib = np.zeros(len(xyz))
    tri_faces = [(0, 1, 2, 4, 5, 6), (0, 1, 3, 4, 9, 7),
                 (1, 2, 3, 5, 8, 9), (0, 2, 3, 6, 8, 7)]
    ntop = 0
    for tet in tets:
        for f in tri_faces:
            n = [tet[q] for q in f]
            if not np.all(np.isclose(z[n], 0.0, atol=TOL)):
                continue
            a, b, c = xyz[n[0]], xyz[n[1]], xyz[n[2]]
            area = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
            for m in n[3:]:
                trib_consist[m] += area / 3.0
            for m in n:
                trib[m] += area / 6.0
            ntop += 1
    print(f"[mesh] {ntop} top faces, tributary sum {trib.sum():.6f} m2 "
          f"against exact {(2 * XLIM) ** 2:.1f} "
          f"(rel {abs(trib.sum() - (2 * XLIM) ** 2) / (2 * XLIM) ** 2:.2e})")
    print(f"[mesh] footing tributary {trib[sets['footing']].sum():.6f} m2 "
          f"against B^2 = {B_FOOT ** 2:.1f}")

    nv_top = [n for n in sets["top"] if is_vertex[n]]
    print(f"[mesh] lumped tributary puts {trib[nv_top].sum():.3f} m2 on the "
          f"{len(nv_top)} surface VERTEX nodes (consistent rule gives "
          f"{trib_consist[nv_top].sum():.3f} — every one unloaded)")
    np.savez_compressed(OUT, nodes=xyz, tets=tets, is_vertex=is_vertex,
                        tributary=trib, tributary_consist=trib_consist,
                        **{f"set_{k}": v for k, v in sets.items()})
    print(f"[mesh] -> {OUT}")


if __name__ == "__main__":
    main()
