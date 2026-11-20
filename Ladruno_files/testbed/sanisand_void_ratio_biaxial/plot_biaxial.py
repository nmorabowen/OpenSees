"""WP-238 demo: plot the runs of biaxial.py.

    python plot_biaxial.py out  [case ...]      -> out/bands.png, out/curves.png, out/summary.md

Per case: the void-ratio field, and the maximum shear strain gamma_max of each element
(from the nodal displacements, at the centroid, RELATIVE to the stage-0 end is not removed:
stage 0 is 1e-3 at most) at the last snapshot, on the deformed mesh. gamma_max =
sqrt((eps_xx - eps_yy)^2 + gamma_xy^2).
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                     # noqa: E402
from matplotlib.collections import PolyCollection                   # noqa: E402


def gamma_max(xy, conn, u):
    """Centroid strain of a bilinear quad from its nodal displacements."""
    X = xy[conn - 1]                       # (ne, 4, 2)
    U = u[conn - 1]
    dN = 0.25 * np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]], float)   # dN/dxi, dN/deta at 0,0
    J = np.einsum("ak,eaj->ekj", dN, X)    # J[k, j] = d x_j / d xi_k
    dNx = np.einsum("ejk,ak->eaj", np.linalg.inv(J), dN)
    G = np.einsum("eaj,eai->eij", dNx, U)  # G[i, j] = d u_i / d x_j
    exx, eyy, gxy = G[:, 0, 0], G[:, 1, 1], G[:, 0, 1] + G[:, 1, 0]
    return np.sqrt((exx - eyy) ** 2 + gxy ** 2), exx + eyy


def band_stats(g, area):
    """How concentrated the shear is: the area fraction that carries half of the total
    shear (0.5 = uniform), and the max / mean ratio."""
    o = np.argsort(g)[::-1]
    c = np.cumsum(g[o] * area[o])
    k = int(np.searchsorted(c, 0.5 * c[-1]))
    return float(np.sum(area[o][:k + 1]) / np.sum(area)), float(g.max() / g.mean())


def main():
    root = sys.argv[1]
    cases = sys.argv[2:] or sorted(d for d in os.listdir(root)
                                   if os.path.isfile(os.path.join(root, d, "run.npz")))
    R = {c: np.load(os.path.join(root, c, "run.npz"), allow_pickle=True) for c in cases}
    n = len(cases)

    # ---- fields and bands ----
    fig, ax = plt.subplots(2, n, figsize=(2.3 * n + 1.2, 8.6), squeeze=False)
    gtop = max(np.percentile(gamma_max(r["xy"], r["conn"], r["snaps"][-1])[0], 99.5)
               for r in R.values())
    emid = 0.6944
    espan = max(float(np.abs(r["e"] - emid).max()) for r in R.values()) or 0.01
    rows = []
    for k, c in enumerate(cases):
        r = R[c]
        xy, conn, u = r["xy"], r["conn"], r["snaps"][-1]
        g, ev = gamma_max(xy, conn, u)
        X0 = xy[conn - 1]
        area = 0.5 * np.abs((X0[:, 2, 0] - X0[:, 0, 0]) * (X0[:, 3, 1] - X0[:, 1, 1])
                            - (X0[:, 3, 0] - X0[:, 1, 0]) * (X0[:, 2, 1] - X0[:, 0, 1]))
        frac, ratio = band_stats(g, area)
        cv = r["curve"]
        if "holds" in r.files:                 # dynamic relaxation: the HELD points are static
            h = r["holds"]
            cv = np.column_stack([h[:, 0], h[:, 2]])
        ip = int(np.argmax(cv[:, 1]))
        rows.append((c, int(r["nx"]), int(r["ny"]), float(r["e"].std()), cv[ip, 1], cv[ip, 0],
                     cv[-1, 1], cv[-1, 0], frac, ratio, str(r["status"])))
        pc = PolyCollection(X0, array=r["e"], cmap="RdBu_r", edgecolors="none")
        pc.set_clim(emid - espan, emid + espan)
        ax[0, k].add_collection(pc)
        ax[0, k].set_title(f"{c}\ninitial void ratio", fontsize=9)
        pd = PolyCollection((xy + u)[conn - 1], array=g, cmap="magma_r", edgecolors="none")
        pd.set_clim(0.0, gtop)
        ax[1, k].add_collection(pd)
        ax[1, k].set_title(f"$\\gamma_{{max}}$ at $\\varepsilon_a$ = {100*r['snap_at'][-1]:.1f} %",
                           fontsize=9)
        for a in ax[:, k]:
            a.set_xlim(-0.03, float(r["W"]) + 0.03)
            a.set_ylim(-0.005, float(r["H"]) + 0.005)
            a.set_aspect("equal")
            a.set_xticks([])
            a.set_yticks([])
    fig.colorbar(pc, ax=ax[0, :], shrink=0.8, label="e")
    fig.colorbar(pd, ax=ax[1, :], shrink=0.8, label=r"$\gamma_{max}$")
    fig.savefig(os.path.join(root, "bands.png"), dpi=150, bbox_inches="tight")

    # ---- curves ----
    fig2, a2 = plt.subplots(figsize=(7.2, 4.6))
    for c in cases:
        cv = R[c]["curve"]
        ln, = a2.plot(100 * cv[:, 0], cv[:, 1], label=c, lw=1.2,
                      ls="--" if "homog" in c else "-")
        if "holds" in R[c].files:              # ramped line, held (static) points as markers
            h = R[c]["holds"]
            ok = h[:, 4] > 0
            a2.plot(100 * h[ok, 0], h[ok, 2], "o", ms=3.5, color=ln.get_color())
            a2.plot(100 * h[~ok, 0], h[~ok, 2], "x", ms=5, color=ln.get_color())
    a2.set_xlabel("axial strain (%)")
    a2.set_ylabel("axial stress (kPa)")
    a2.grid(alpha=0.3)
    a2.legend(fontsize=8)
    fig2.savefig(os.path.join(root, "curves.png"), dpi=150, bbox_inches="tight")

    with open(os.path.join(root, "summary.md"), "w", encoding="utf-8") as f:
        f.write("| case | mesh | std(e) | peak sigma_a (kPa) | at eps_a (%) | end sigma_a (kPa) | "
                "end eps_a (%) | area carrying half the shear | gamma max/mean | status |\n")
        f.write("|---|---|---|---|---|---|---|---|---|---|\n")
        for c, nx, ny, sd, sp, ep, se, ee, fr, ra, st in rows:
            f.write(f"| {c} | {nx}x{ny} | {sd:.4f} | {sp:.1f} | {100*ep:.2f} | {se:.1f} | "
                    f"{100*ee:.2f} | {fr:.3f} | {ra:.1f} | {st} |\n")
    print(open(os.path.join(root, "summary.md"), encoding="utf-8").read())


if __name__ == "__main__":
    main()
