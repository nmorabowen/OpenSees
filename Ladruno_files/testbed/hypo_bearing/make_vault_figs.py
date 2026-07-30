"""Figures for TIMs vault note 16 — the geometric-kinematics study.

Obeys the vault's visual contract ([[02 — Style, Figures and Colour]]):
Archivo Narrow registered explicitly from the per-user font directory (the
matplotlib cache misses it otherwise and falls back SILENTLY, so the resolved
family is printed), Okabe-Ito palette, no red-versus-green encoding, and every
colour distinction backed by a second channel (line style + direct label) so
the figures survive greyscale.

Produces, into the vault attachments directory:
  hypo_mesh.png      the graded benchmark domain and a section through it
  hypo_curves.png    the four-rung backbones, absolute and Vesic-normalized
  hypo_ladder.png    each rung against the linear base — the decomposition
  hypo_tangent.png   dq/ds, the limit-point diagnostic

Run:  LADRUNO_OPENSEES_QUIET=1 python make_vault_figs.py
Needs the apeGmsh env (matplotlib); reads only committed CSV + npz artifacts.
"""
import csv
import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

userdir = os.path.join(os.environ.get("LOCALAPPDATA", ""),
                       "Microsoft", "Windows", "Fonts")
if os.path.isdir(userdir):
    for fn in os.listdir(userdir):
        if "archivo" in fn.lower() and fn.lower().endswith((".ttf", ".otf")):
            fm.fontManager.addfont(os.path.join(userdir, fn))
_have = {f.name for f in fm.fontManager.ttflist}
FONT = "Archivo Narrow" if "Archivo Narrow" in _have else "DejaVu Sans"
print("font resolved:", FONT)

CB = dict(blue="#0072B2", orange="#E69F00", sky="#56B4E9",
          vermillion="#D55E00", green="#009E73", yellow="#F0E442",
          purple="#CC79A7", grey="#8C8C8C", black="#000000")

mpl.rcParams.update({
    "font.family": FONT, "font.size": 12,
    "axes.grid": True, "grid.color": "#DDDDDD", "grid.linewidth": 0.7,
    "axes.edgecolor": "#444444", "axes.linewidth": 0.9,
    "axes.titlesize": 13, "axes.labelsize": 12,
    "legend.frameon": False, "figure.dpi": 130, "savefig.dpi": 200,
    "savefig.bbox": "tight",
})

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = r"C:\Users\nmb\Dropbox\obsidian\Ladruño\TIMs\attachments"
B_FOOT, Q_VESIC = 2.0, 637.5

# grey = baseline/de-emphasised, orange = the rung that goes furthest.
# Style is the second channel: base dash-dot, corot solid, hypo solid,
# hypo+kc dashed on top of hypo so coincidence is visible rather than hidden.
LEGS = [
    ("linear",  "linear (base)",      dict(color=CB["grey"],       ls="-.", lw=2.0)),
    ("corot",   "corot",              dict(color=CB["blue"],       ls="-",  lw=2.0)),
    ("hypo",    "hypo",               dict(color=CB["orange"],     ls="-",  lw=2.2)),
    ("hypo_kc", "hypo + kozenyCarman", dict(color=CB["vermillion"], ls="--", lw=1.4)),
]


def load(leg):
    p = os.path.join(HERE, f"backbone_{leg}.csv")
    if not os.path.exists(p):
        return None
    rows = []
    with open(p) as f:
        for r in list(csv.reader(f))[1:]:
            try:
                rows.append([float(v) for v in r])
            except ValueError:
                continue
    if not rows:
        return None
    a = np.array(rows)
    return dict(s=a[:, 0], q=a[:, 3], qn=a[:, 4])


def smooth_tangent(d, w=25):
    k = np.gradient(d["q"], d["s"])
    if len(k) <= 2 * w:
        return d["s"], k
    ker = np.ones(w) / w
    return d["s"][w - 1:], np.convolve(k, ker, mode="valid")


# ---------------------------------------------------------------- mesh figure
def fig_mesh(data):
    z = np.load(os.path.join(HERE, "bearing_mesh.npz"))
    xyz, hexes = z["nodes"], z["hexes"]
    FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    # exterior faces = those belonging to exactly one hex
    seen = {}
    for h in hexes:
        for f in FACES:
            key = tuple(sorted(int(h[i]) for i in f))
            seen[key] = seen.get(key, 0) + 1
    ext = []
    foot = []
    for h in hexes:
        for f in FACES:
            idx = [int(h[i]) for i in f]
            if seen[tuple(sorted(idx))] != 1:
                continue
            q = xyz[idx]
            if np.all(np.abs(q[:, 2]) < 1e-6) and \
               abs(q[:, 0].mean()) < B_FOOT / 2 and abs(q[:, 1].mean()) < B_FOOT / 2:
                foot.append(q)
            else:
                ext.append(q)

    fig = plt.figure(figsize=(13.5, 5.2))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    ax.add_collection3d(Poly3DCollection(
        ext, facecolor="#F2F2F2", edgecolor="#9A9A9A", linewidths=0.30, alpha=1.0))
    ax.add_collection3d(Poly3DCollection(
        foot, facecolor=CB["orange"], edgecolor=CB["black"], linewidths=0.6))
    ax.set_xlim(-10, 10); ax.set_ylim(-10, 10); ax.set_zlim(-8, 0)
    ax.set_box_aspect((20, 20, 8))
    ax.view_init(elev=26, azim=-56)
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_zlabel("z [m]")
    ax.set_title("2816 LadrunoUP H8 — 0.5 m under the footing, graded outward\n"
                 "footing patch B = 2 m (orange), 4.5 B clearance, 4 B depth",
                 fontsize=11)
    ax.grid(False)

    # section at y = 0: draw the element edges of the slab straddling y = 0
    ax2 = fig.add_subplot(1, 2, 2)
    tol = 1e-6
    for h in hexes:
        q = xyz[[int(i) for i in h]]
        if not (q[:, 1].min() < tol and q[:, 1].max() > -tol):
            continue
        if q[:, 1].min() < -0.5:          # keep only the slab just below y=0
            continue
        xs, zs = q[:, 0], q[:, 2]
        x0, x1, z0, z1 = xs.min(), xs.max(), zs.min(), zs.max()
        ax2.add_patch(plt.Rectangle((x0, z0), x1 - x0, z1 - z0,
                                    facecolor="#F2F2F2", edgecolor="#9A9A9A",
                                    linewidth=0.5))
    ax2.plot([-B_FOOT / 2, B_FOOT / 2], [0, 0], color=CB["orange"], lw=4,
             solid_capstyle="butt", label="driven footing patch")
    ax2.annotate("B = 2 m, 4 elements", xy=(0, 0), xytext=(3.2, -1.2),
                 fontsize=10, color=CB["black"],
                 arrowprops=dict(arrowstyle="->", color=CB["black"], lw=0.9))
    ax2.set_xlim(-10.3, 10.3); ax2.set_ylim(-8.3, 0.9)
    ax2.set_aspect("equal")
    ax2.set_xlabel("x [m]"); ax2.set_ylabel("z [m]")
    ax2.set_title("Section at y = 0 — the grading the campaign depends on\n"
                  "0.5 m near-field, r ≈ 1.44 in plan and 1.35 in depth",
                  fontsize=11)
    ax2.legend(loc="lower right", fontsize=10)
    fig.tight_layout()
    p = os.path.join(OUT, "hypo_mesh.png")
    fig.savefig(p); plt.close(fig)
    print("wrote", p)


# -------------------------------------------------------------- curve figures
def fig_curves(data):
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 4.8))
    for leg, lbl, st in LEGS:
        d = data.get(leg)
        if d is None:
            continue
        x = d["s"] / B_FOOT * 100
        ax[0].plot(x, d["q"], label=lbl, **st)
        ax[1].plot(x, d["qn"], label=lbl, **st)
        ax[0].plot(x[-1], d["q"][-1], "o", ms=6, color=st["color"], zorder=5)
    ax[1].axhline(1.0, color=CB["black"], ls=":", lw=1.0)
    ax[1].annotate("Vesic capacity", xy=(0.4, 1.0), xytext=(0.6, 1.45),
                   fontsize=10, arrowprops=dict(arrowstyle="->", lw=0.8))
    ax[0].set_ylabel("footing pressure $q$ [kPa]")
    ax[1].set_ylabel("$q\\,/\\,q_\\mathrm{Vesic}$")
    ax[0].set_title("No rung develops a limit point\n"
                    "markers = last converged point", fontsize=11)
    ax[1].set_title("All four cross Vesic at $s/B \\approx 1\\,\\%$ and keep climbing",
                    fontsize=11)
    for a in ax:
        a.set_xlabel("settlement $s/B$ [%]")
        a.legend(loc="upper left", fontsize=10)
    fig.tight_layout()
    p = os.path.join(OUT, "hypo_curves.png")
    fig.savefig(p); plt.close(fig)
    print("wrote", p)


def fig_ladder(data):
    base = data.get("linear")
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    if base is None:
        print("no base leg — skipping ladder figure")
        plt.close(fig); return
    smax = base["s"][-1]
    grid = np.linspace(base["s"][20], smax, 200)
    qb = np.interp(grid, base["s"], base["q"])
    for leg, lbl, st in LEGS:
        if leg == "linear":
            continue
        d = data.get(leg)
        if d is None:
            continue
        ax.plot(grid / B_FOOT * 100,
                (np.interp(grid, d["s"], d["q"]) / qb - 1.0) * 100,
                label=f"{lbl} vs base", **st)
    ax.axhline(0.0, color=CB["grey"], ls="-.", lw=2.0)
    ax.annotate("linear base", xy=(3.5, 0), xytext=(3.6, 0.8), fontsize=10,
                color=CB["black"])
    ax.set_xlabel("settlement $s/B$ [%]")
    ax.set_ylabel("change in $q$ against the linear base [%]")
    ax.set_title("The rungs pull in opposite directions, and both are small\n"
                 "rotation softens; genuine large strain stiffens", fontsize=11)
    ax.legend(loc="upper left", fontsize=10)
    fig.tight_layout()
    p = os.path.join(OUT, "hypo_ladder.png")
    fig.savefig(p); plt.close(fig)
    print("wrote", p)


def fig_tangent(data):
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    for leg, lbl, st in LEGS:
        d = data.get(leg)
        if d is None:
            continue
        s, k = smooth_tangent(d)
        ax.plot(s / B_FOOT * 100, k / 1000.0, label=lbl, **st)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("settlement $s/B$ [%]")
    ax.set_ylabel("tangent $\\mathrm{d}q/\\mathrm{d}s$ [MPa/m]")
    ax.set_title("A limit point would bring this to zero; none does\n"
                 "hypo retains the most tangent of the three", fontsize=11)
    ax.legend(loc="upper right", fontsize=10)
    fig.tight_layout()
    p = os.path.join(OUT, "hypo_tangent.png")
    fig.savefig(p); plt.close(fig)
    print("wrote", p)


def fig_probe():
    """The single-element deviatoric probe: does the MATERIAL plateau?"""
    f = os.path.join(HERE, "element_probe.npz")
    if not os.path.exists(f):
        print("no element_probe.npz - skipping probe figure")
        return
    z = np.load(f)
    keys = sorted(z.files, key=lambda k: float(k.split("_")[1]))
    cols = [CB["blue"], CB["sky"], CB["orange"], CB["vermillion"], CB["black"]]
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 4.6))
    for c, k in zip(cols, keys):
        sv = float(k.split("_")[1])
        r = z[k]
        lbl = rf"$\sigma_v$ = {sv:g} kPa"
        ax[0].plot(r[:, 0], r[:, 1] / sv, color=c, lw=1.8, label=lbl)
        ax[1].plot(r[:, 0], r[:, 4], color=c, lw=1.8, label=lbl)
    ax[0].set_ylabel(r"$\tau / \sigma_v$")
    ax[0].set_title("The material DOES reach perfect plasticity\n"
                    "flat to 1.0000 of peak, and the ratio is "
                    "confinement-independent", fontsize=11)
    sin33 = np.sin(np.radians(33.0))
    ax[1].axhline(sin33, color=CB["grey"], ls="-.", lw=2.0)
    ax[1].annotate(r"$\sin\varphi$ at the input $\varphi = 33\degree$",
                   xy=(0.30, sin33), xytext=(0.20, 0.40), fontsize=10)
    ax[1].set_ylabel(r"$\sin\varphi_\mathrm{mob}$")
    ax[1].set_ylim(0, 0.95)
    ax[1].set_title("...but it mobilises 50.4 deg, not 33 deg\n"
                    "PDMY's cone has no Lode dependence and is stronger "
                    "in shear", fontsize=11)
    for a in ax:
        a.set_xlabel(r"imposed shear strain $\gamma$")
        a.legend(fontsize=9, loc="lower right")
    fig.tight_layout()
    out = os.path.join(OUT, "hypo_probe.png")
    fig.savefig(out); plt.close(fig)
    print("wrote", out)


def main():
    os.makedirs(OUT, exist_ok=True)
    data = {leg: load(leg) for leg, _, _ in LEGS}
    for leg, _, _ in LEGS:
        d = data[leg]
        print(f"  {leg:8s}: "
              + ("absent" if d is None
                 else f"{len(d['s'])} steps to s/B={d['s'][-1] / B_FOOT:.4f}, "
                      f"q_end={d['q'][-1]:.0f} kPa = {d['qn'][-1]:.2f} x Vesic"))
    fig_mesh(data)
    fig_curves(data)
    fig_ladder(data)
    fig_tangent(data)
    fig_probe()


if __name__ == "__main__":
    main()
