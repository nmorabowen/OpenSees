"""WP-238 demo: one run of biaxial.py at its LAST saved state, three panels: the void-ratio
field, the total gamma_max, and the gamma_max of the last displacement increment (between
the last two snapshots), each on the deformed mesh.

    python plot_stall.py out/<case> [png]
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                     # noqa: E402
from matplotlib.collections import PolyCollection                   # noqa: E402

from plot_biaxial import gamma_max                                  # noqa: E402


def main():
    d = sys.argv[1]
    r = np.load(os.path.join(d, "run.npz"), allow_pickle=True)
    xy, conn, S, at = r["xy"], r["conn"], r["snaps"], r["snap_at"]
    g, _ = gamma_max(xy, conn, S[-1])
    k0 = len(S) - 2
    dg, _ = gamma_max(xy, conn, S[-1] - S[k0])
    de = at[-1] - at[k0]
    X0, X1 = xy[conn - 1], (xy + S[-1])[conn - 1]
    fig, ax = plt.subplots(1, 3, figsize=(9.6, 5.6))
    panels = [(X0, r["e"], "RdBu_r", "initial void ratio", None),
              (X1, g, "magma_r", f"$\\gamma_{{max}}$ at $\\varepsilon_a$ = {100*at[-1]:.2f} %", 0.0),
              (X1, dg / de, "magma_r",
               f"$\\Delta\\gamma_{{max}}/\\Delta\\varepsilon_a$, "
               f"{100*at[k0]:.2f} to {100*at[-1]:.2f} %", 0.0)]
    for a, (X, v, cm, title, lo) in zip(ax, panels):
        pc = PolyCollection(X, array=v, cmap=cm, edgecolors="none")
        if lo is not None:
            pc.set_clim(lo, float(np.percentile(v, 99.5)))
        a.add_collection(pc)
        a.set_title(title, fontsize=9)
        a.set_xlim(-0.01, float(r["W"]) + 0.01)
        a.set_ylim(-0.005, float(r["H"]) + 0.005)
        a.set_aspect("equal")
        a.set_xticks([])
        a.set_yticks([])
        fig.colorbar(pc, ax=a, shrink=0.75)
    fig.suptitle(f"{os.path.basename(d)}: {r['status']}", fontsize=9)
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(d, "stall.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(out, "| gamma max/mean total", round(float(g.max() / g.mean()), 2),
          "| incremental", round(float(dg.max() / dg.mean()), 2),
          "| uniform would be 1.0")


if __name__ == "__main__":
    main()
