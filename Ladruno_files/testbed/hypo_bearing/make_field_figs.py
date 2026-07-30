"""Deformation and stress figures from fields_hypo.npz, for TIMs vault note 17.

Three figures, each making one point:

  hypo_deform.png   the deformed section and the surface profile. A forming
                    Prandtl mechanism HEAVES beside the footing; pure punching
                    does not. This is the picture that distinguishes them.
  hypo_stress.png   mean effective stress p' beside the settled state, which
                    shows whether the load is being carried by a mechanism or
                    by a growing confined bulb under the footing.
  hypo_mobilised.png  sin(phi_mob) against the MEASURED failure value of 0.770
                    (element_probe.py) rather than sin(33 deg) -- where the soil
                    is actually at its strength.

Vault visual contract: Archivo Narrow (resolution asserted), cividis for
continuous fields, no red-vs-green.

Run:  LADRUNO_OPENSEES_QUIET=1 python make_field_figs.py
"""
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.collections import PolyCollection

userdir = os.path.join(os.environ.get("LOCALAPPDATA", ""),
                       "Microsoft", "Windows", "Fonts")
if os.path.isdir(userdir):
    for fn in os.listdir(userdir):
        if "archivo" in fn.lower() and fn.lower().endswith((".ttf", ".otf")):
            fm.fontManager.addfont(os.path.join(userdir, fn))
FONT = ("Archivo Narrow" if "Archivo Narrow" in
        {f.name for f in fm.fontManager.ttflist} else "DejaVu Sans")
print("font resolved:", FONT)

CB = dict(blue="#0072B2", orange="#E69F00", sky="#56B4E9",
          vermillion="#D55E00", grey="#8C8C8C", black="#000000")
mpl.rcParams.update({
    "font.family": FONT, "font.size": 12,
    "axes.grid": True, "grid.color": "#DDDDDD", "grid.linewidth": 0.7,
    "axes.edgecolor": "#444444", "axes.linewidth": 0.9,
    "axes.titlesize": 12, "axes.labelsize": 12,
    "legend.frameon": False, "figure.dpi": 130, "savefig.dpi": 200,
    "savefig.bbox": "tight",
})

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = r"C:\Users\nmb\Dropbox\obsidian\Ladruño\TIMs\attachments"
B = 2.0
MAG = 8.0                 # displacement magnification for the deformed section


def slab(nodes, hexes):
    """Elements in the y-slab just below y = 0, as (x, z) quads for plotting."""
    keep, quads = [], []
    for e, h in enumerate(hexes):
        q = nodes[[int(i) for i in h]]
        if q[:, 1].min() < -0.5 or not (q[:, 1].min() < 1e-6 < q[:, 1].max() + 1e-6):
            continue
        keep.append(e)
        x0, x1 = q[:, 0].min(), q[:, 0].max()
        z0, z1 = q[:, 2].min(), q[:, 2].max()
        quads.append([(x0, z0), (x1, z0), (x1, z1), (x0, z1)])
    return np.array(keep), quads


def deformed_quads(nodes, hexes, keep, u):
    """Same slab but with nodal displacement (magnified) applied."""
    out = []
    for e in keep:
        h = [int(i) for i in hexes[e]]
        q = nodes[h] + MAG * u[h, :3]
        x0, x1 = q[:, 0].min(), q[:, 0].max()
        z0, z1 = q[:, 2].min(), q[:, 2].max()
        out.append([(x0, z0), (x1, z0), (x1, z1), (x0, z1)])
    return out


def main():
    f = os.path.join(HERE, "fields_hypo.npz")
    if not os.path.exists(f):
        raise SystemExit("fields_hypo.npz missing — run field_run.py first")
    z = np.load(f)
    nodes, hexes = z["nodes"], z["hexes"]
    cps = sorted(float(k.split("_")[1]) for k in z.files if k.startswith("u_")
                 and k != "u_grav")
    if not cps:
        raise SystemExit("no checkpoints in fields_hypo.npz yet")
    last = cps[-1]
    tag = f"{last:.4f}"
    print(f"checkpoints available: {cps}; plotting s/B = {last}")
    u = z[f"u_{tag}"]
    sig = z[f"sig_{tag}"]
    sphi = z[f"sphi_{tag}"]
    sin_fail = float(z["sin_fail"])
    keep, quads = slab(nodes, hexes)

    # ---------------------------------------------------------- deformation
    fig, ax = plt.subplots(1, 2, figsize=(13.8, 4.6))
    ax[0].add_collection(PolyCollection(quads, facecolors="none",
                                        edgecolors="#C8C8C8", linewidths=0.5))
    ax[0].add_collection(PolyCollection(
        deformed_quads(nodes, hexes, keep, u), facecolors="none",
        edgecolors=CB["blue"], linewidths=0.7))
    ax[0].plot([-B / 2, B / 2], [MAG * u[:, 2].min() * 0 + 0, 0], lw=0)
    ax[0].set_xlim(-6.5, 6.5); ax[0].set_ylim(-6.0, 1.2)
    ax[0].set_aspect("equal")
    ax[0].set_xlabel("x [m]"); ax[0].set_ylabel("z [m]")
    ax[0].set_title(f"Deformed section at $s/B$ = {last:.0%}, "
                    f"displacement x{MAG:g}\ngrey = undeformed, "
                    f"blue = deformed", fontsize=11)

    top = np.abs(nodes[:, 2]) < 1e-6
    band = top & (np.abs(nodes[:, 1]) < 1e-6)
    o = np.argsort(nodes[band, 0])
    xs = nodes[band, 0][o]
    uz = u[band, 2][o] * 1000.0
    ax[1].plot(xs, uz, "-o", ms=4, color=CB["blue"], lw=1.8)
    ax[1].axvspan(-B / 2, B / 2, color=CB["orange"], alpha=0.25)
    ax[1].annotate("footing", xy=(0, uz.min() * 0.5), ha="center", fontsize=10)
    ax[1].axhline(0.0, color=CB["black"], ls=":", lw=0.9)
    ax[1].set_xlim(-10.3, 10.3)
    ax[1].set_xlabel("x [m]"); ax[1].set_ylabel("surface $u_z$ [mm]")
    heave = uz[np.abs(xs) > B / 2].max()
    ax[1].set_title(f"Surface profile — max heave beside the footing "
                    f"{heave:+.2f} mm\na Prandtl mechanism heaves; punching "
                    f"does not", fontsize=11)
    fig.tight_layout()
    p = os.path.join(OUT, "hypo_deform.png")
    fig.savefig(p); plt.close(fig)
    print("wrote", p, f"(max heave {heave:+.3f} mm)")

    # -------------------------------------------------------------- stresses
    pmean = -(sig[:, 0] + sig[:, 1] + sig[:, 2]) / 3.0
    fig, ax = plt.subplots(1, 2, figsize=(13.8, 4.6))
    for a, v, ttl, cb_lbl in (
            (ax[0], pmean[keep], "Mean effective stress $p'$ — the load is "
             "carried by a\nconfined bulb under the footing, not by a "
             "slip surface", "$p'$ [kPa]"),
            (ax[1], z[f"J_{tag}"][keep], "Volume ratio $J$ from `hypoState` — "
             "below 1 is compaction\nthe near-footing soil is being crushed, "
             "not sheared aside", "$J$ [-]")):
        pc = PolyCollection(quads, array=v, cmap="cividis",
                            edgecolors="#00000018", linewidths=0.2)
        a.add_collection(pc)
        a.set_xlim(-6.5, 6.5); a.set_ylim(-6.0, 0.4)
        a.set_aspect("equal"); a.set_axisbelow(True)
        a.set_xlabel("x [m]"); a.set_ylabel("z [m]")
        a.set_title(ttl, fontsize=11)
        plt.colorbar(pc, ax=a, shrink=0.85, label=cb_lbl)
    fig.tight_layout()
    p = os.path.join(OUT, "hypo_stress.png")
    fig.savefig(p); plt.close(fig)
    print("wrote", p)

    # ------------------------------------------------------ mobilised friction
    fig, ax = plt.subplots(1, 2, figsize=(13.8, 4.6))
    v = sphi[keep] / sin_fail
    pc = PolyCollection(quads, array=v, cmap="cividis", vmin=0.0, vmax=1.0,
                        edgecolors="#00000018", linewidths=0.2)
    ax[0].add_collection(pc)
    ax[0].set_xlim(-6.5, 6.5); ax[0].set_ylim(-6.0, 0.4)
    ax[0].set_aspect("equal")
    ax[0].set_xlabel("x [m]"); ax[0].set_ylabel("z [m]")
    ax[0].set_title("Mobilised fraction of the MEASURED strength\n"
                    r"$\sin\varphi_\mathrm{mob}/0.770$ — 1.0 is at failure",
                    fontsize=11)
    plt.colorbar(pc, ax=ax[0], shrink=0.85, label="fraction of strength")

    good = np.isfinite(sphi)
    ax[1].hist(sphi[good] / sin_fail, bins=40, range=(0, 1.2),
               color=CB["blue"], edgecolor="white")
    ax[1].axvline(1.0, color=CB["vermillion"], ls="--", lw=2.0)
    ax[1].annotate("measured failure\n(0.770)", xy=(1.0, 0), xytext=(0.62, 0.62),
                   textcoords="axes fraction", fontsize=10,
                   arrowprops=dict(arrowstyle="->", lw=0.9))
    frac = float(np.mean(sphi[good] >= 0.95 * sin_fail))
    ax[1].set_xlabel(r"$\sin\varphi_\mathrm{mob}/0.770$")
    ax[1].set_ylabel("elements")
    ax[1].set_title(f"Only {frac:.1%} of elements are within 5 % of failure\n"
                    f"the mechanism is nowhere near fully mobilised",
                    fontsize=11)
    fig.tight_layout()
    p = os.path.join(OUT, "hypo_mobilised.png")
    fig.savefig(p); plt.close(fig)
    print("wrote", p, f"(fraction within 5% of failure: {frac:.3f})")

    print(f"\nq at this checkpoint = {float(z[f'q_{tag}']):.1f} kPa")


if __name__ == "__main__":
    main()
