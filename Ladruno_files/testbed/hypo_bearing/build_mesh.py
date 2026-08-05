"""ADR-79 bearing campaign — graded all-hex benchmark mesh (apeGmsh).

Scoping finding 3 of this testbed: every uniform-1 m-hex box small enough for
a quick run is an OEDOMETER — roller sides at 0.5-1.5 B clearance confine the
mechanism, PDMY's G proportional to sqrt(p') self-reinforces, and q runs to
tens of MPa with no bearing meaning. A credible bearing verdict needs >= 4-5 B
of side clearance and ~3-4 B of depth, which at a uniform 0.5 m footing-scale
resolution would be ~1e5 elements. The answer is a GRADED mesh.

Domain (B = 2 m square footing centred at the origin on the top face):

    x, y in [-10, +10] m   -> 9 m of clearance each side = 4.5 B
    z    in [ -8,   0] m   -> 4 B deep, top face z = 0 is the ground surface

Structured tensor-product grading, built as a 3 x 3 x 2 block decomposition
fragmented into one conformal solid, every block transfinite + recombined so
the result is pure H8 (LadrunoUP's hypo lane is H8-only):

    x: [-10,-1] 6 el (graded, fine at -1) | [-1,1] 4 el @ 0.5 m | [1,10] 6 el
    y: same as x
    z: [-8,-3]  5 el (graded, fine at -3) | [-3,0] 6 el @ 0.5 m

=> 16 x 16 x 11 = 2816 hexes, 3468 nodes, 13 872 u-p DOF. The footing is the
central 4 x 4-element patch of the top face (5 x 5 nodes, exactly 2 x 2 m), so
B is resolved by 4 elements -- coarse for an absolute capacity, adequate for
the CROSS-METHOD comparison this campaign is actually about.

Run:  C:\\Users\\nmb\\venv\\opensees_env\\Scripts\\python.exe build_mesh.py
Writes bearing_mesh.npz next to this script (nodes, hexes, node sets).
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "bearing_mesh.npz")

B_FOOT = 2.0
XLIM, ZBOT = 10.0, -8.0
H0 = 0.5                       # target element size under and near the footing

# (lo, hi, n_elements, coordinate at which the elements are FINE; None = uniform)
XBLOCKS = [(-XLIM, -1.0, 6, -1.0), (-1.0, 1.0, 4, None), (1.0, XLIM, 6, 1.0)]
YBLOCKS = XBLOCKS
ZBLOCKS = [(ZBOT, -3.0, 5, -3.0), (-3.0, 0.0, 6, None)]
TOL = 1e-6


def growth(length, n, h_first):
    """Progression ratio r with sum_{i<n} h_first*r^i == length (bisection)."""
    if abs(h_first * n - length) < 1e-12:
        return 1.0

    def total(r):
        return h_first * n if abs(r - 1.0) < 1e-12 else h_first * (r ** n - 1.0) / (r - 1.0)

    lo, hi = 1.0, 4.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if total(mid) < length:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def spec(blocks):
    """(lo, hi) -> (n_nodes, ratio, fine_coord) for the transfinite curves."""
    out = {}
    for lo, hi, n, fine in blocks:
        r = 1.0 if fine is None else growth(hi - lo, n, H0)
        out[(lo, hi)] = (n + 1, r, fine)
    return out


def main():
    from apeGmsh import apeGmsh
    import gmsh

    sx, sy, sz = spec(XBLOCKS), spec(YBLOCKS), spec(ZBLOCKS)
    for ax, s in (("x", sx), ("y", sy), ("z", sz)):
        for (lo, hi), (nn, r, fine) in s.items():
            print(f"  {ax} [{lo:+.1f},{hi:+.1f}] {nn - 1} el  r={r:.4f}  "
                  f"h_max={H0 * r ** (nn - 2):.3f} m")

    with apeGmsh(model_name="adr79_bearing") as g:
        for i, (x0, x1, *_ ) in enumerate(XBLOCKS):
            for j, (y0, y1, *_) in enumerate(YBLOCKS):
                for k, (z0, z1, *_) in enumerate(ZBLOCKS):
                    g.model.geometry.add_box(x0, y0, z0, x1 - x0, y1 - y0,
                                             z1 - z0, label=f"b{i}{j}{k}")
        vols = [t for (d, t) in gmsh.model.getEntities(3)]
        g.model.boolean.fragment([(3, vols[0])], [(3, t) for t in vols[1:]])
        gmsh.model.occ.synchronize()

        # every axis-aligned curve gets the transfinite spec of the (axis,
        # interval) it spans; orientation is read from the curve's own
        # parametrisation so "fine end" always lands where we asked.
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
                raise RuntimeError(f"curve {ct} spans {lo}..{hi} on axis "
                                   f"{axis} — no block matches")
            nn, r, fine = s[key]
            if r == 1.0:
                coef = 1.0
            else:
                # progression grows from the curve's START; invert when the
                # start is the coarse end.
                coef = r if abs(p0[axis] - fine) < TOL else 1.0 / r
            g.mesh.structured.set_transfinite_curve(ct, nn, coef=coef)

        for _, st in gmsh.model.getEntities(2):
            g.mesh.structured.set_transfinite_surface(st)
            g.mesh.structured.set_recombine(st, dim=2)
        for _, vt in gmsh.model.getEntities(3):
            g.mesh.structured.set_transfinite_volume(vt)

        g.mesh.generation.generate(dim=3)
        fem = g.mesh.queries.get_fem_data(dim=3)

        ntags, coords = gmsh.model.mesh.getNodes()[:2]
        coords = np.asarray(coords).reshape(-1, 3)
        order = np.argsort(ntags)
        gid = np.asarray(ntags)[order]
        xyz = coords[order]
        remap = {int(t): i for i, t in enumerate(gid)}

        hexes = []
        for et, _, enodes in zip(*gmsh.model.mesh.getElements(dim=3)):
            if et != 5:
                raise RuntimeError(f"non-hex element type {et} in the mesh")
            conn = np.asarray(enodes).reshape(-1, 8)
            hexes.append(np.vectorize(remap.get)(conn))
        hexes = np.vstack(hexes)
        print(f"[mesh] {len(xyz)} nodes, {len(hexes)} hexes  ({fem.info})")

    # positive-Jacobian check (OpenSees bricks need it); gmsh hex8 ordering is
    # already bottom-CCW / top-CCW, so this should be a no-op — but verify.
    p = xyz[hexes]
    vol = np.einsum("ij,ij->i", p[:, 6] - p[:, 0],
                    np.cross(p[:, 1] - p[:, 0], p[:, 3] - p[:, 0]))
    flip = vol < 0
    if flip.any():
        print(f"[mesh] flipping {flip.sum()} inverted hexes")
        hexes[flip] = hexes[flip][:, [4, 5, 6, 7, 0, 1, 2, 3]]
        p = xyz[hexes]
        vol = np.einsum("ij,ij->i", p[:, 6] - p[:, 0],
                        np.cross(p[:, 1] - p[:, 0], p[:, 3] - p[:, 0]))
    assert (vol > 0).all(), "inverted hexes survive the flip"

    top = np.abs(xyz[:, 2]) < TOL
    inside = (np.abs(xyz[:, 0]) <= B_FOOT / 2 + TOL) & (np.abs(xyz[:, 1]) <= B_FOOT / 2 + TOL)
    sets = dict(
        top=np.where(top)[0],
        bottom=np.where(np.abs(xyz[:, 2] - ZBOT) < TOL)[0],
        xface=np.where(np.abs(np.abs(xyz[:, 0]) - XLIM) < TOL)[0],
        yface=np.where(np.abs(np.abs(xyz[:, 1]) - XLIM) < TOL)[0],
        footing=np.where(top & inside)[0],
    )
    for k, v in sets.items():
        print(f"[sets] {k}: {len(v)} nodes")
    assert len(sets["footing"]) == 25, sets["footing"]

    # Surcharge tributary areas over the WHOLE top face, from the top quad
    # faces of the surface hexes. Whole-face coverage is not a stylistic
    # choice: an earlier build of this mesh excluded the footing patch (to
    # keep q0 a purely adjacent overburden and avoid biasing the measured
    # reaction) and every leg then reproduced scoping finding 2 exactly —
    # unbounded first Newton iterate (corot reported inverted elements) then a
    # KrylovNewton stall, on push step 1. The zero-confinement DOFs that need
    # regularizing are the ones UNDER the footing. The reaction bias this
    # reintroduces is exactly q0 * (footing tributary area) and the runner
    # subtracts it in closed form.
    HEX_FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
                 (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    trib = np.zeros(len(xyz))
    for h in hexes:
        for f in HEX_FACES:
            face = h[list(f)]
            q = xyz[face]
            if np.any(np.abs(q[:, 2]) > TOL):
                continue
            a = 0.5 * abs(np.cross(q[2, :2] - q[0, :2], q[3, :2] - q[1, :2]))
            trib[face] += a / 4.0
    # The footing-node share EXCEEDS B^2: perimeter nodes of the patch also
    # collect tributary from the faces just outside it. That is the ordinary
    # nodal-load ambiguity at a patch boundary, and it is fine — the runner's
    # reaction correction uses this exact sum (it recovers the internal force
    # at the constrained nodes identically), not an idealized q0*B^2.
    fshare = float(trib[sets["footing"]].sum())
    print(f"[sets] surcharge area {trib.sum():.2f} m^2 "
          f"(top face = {(2 * XLIM) ** 2:.2f}); footing-node share "
          f"{fshare:.3f} m^2 (B^2 = {B_FOOT ** 2:.2f})")
    assert abs(trib.sum() - (2 * XLIM) ** 2) < 1e-6
    assert B_FOOT ** 2 < fshare < 2.5 * B_FOOT ** 2

    np.savez_compressed(OUT, nodes=xyz, hexes=hexes.astype(np.int32),
                        tributary=trib,
                        **{f"set_{k}": v.astype(np.int32) for k, v in sets.items()})
    print(f"[out] {OUT}  ({os.path.getsize(OUT) / 1024:.0f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
