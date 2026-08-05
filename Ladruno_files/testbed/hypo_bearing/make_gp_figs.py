"""Gauss-point field figures: stresses and equivalent PLASTIC strain.

Reads a `.ladruno` written by viewer_run.py with `-E stress strain porePressure`
and produces, for TIMs vault note 17:

  hypo_gp_stress.png   mean effective stress p' and deviatoric stress q
  hypo_gp_plastic.png  equivalent plastic deviatoric strain, and the total
                       deviatoric strain it is derived from

WHY THE PLASTIC STRAIN IS DERIVED, AND WHAT THAT COSTS
------------------------------------------------------
`PressureDependMultiYield` exposes NO plastic strain: its setResponse offers
only stress / strain / tangent / backbone (PressureDependMultiYield.cpp
:1306-1326). So the plastic part is recovered by removing the elastic part
from the total strain, using PDMY's own pressure-dependent moduli

    G(p') = G_r (p'/p_ref)^d ,      K(p') = B_r (p'/p_ref)^d

    e_e = s / (2G)            (deviatoric)
    eps_v,e = p' / K          (volumetric)
    eps_p = eps - e_e - eps_v,e/3 * I

and then the equivalent plastic deviatoric strain

    eps_q^p = sqrt(2/3 * e_p : e_p).

This is an APPROXIMATION and its assumption should be read with the figure:
PDMY's elastic moduli are functions of the current p', so the elastic strain is
not a state function of the current stress alone along a path where p' changed,
and the recovery is exact only to the extent that G and K are evaluated at the
state where the elastic strain was accumulated. The TOTAL deviatoric strain is
plotted alongside precisely so the reader can see how much of it the
correction removed; where the two panels agree, the strain is essentially all
plastic and the approximation does not matter.

Everything is referenced to the FIRST recorded frame, which is the gravity +
surcharge state, so the fields show what the footing push did rather than the
datum it started from.

Run:  LADRUNO_OPENSEES_QUIET=1 python make_gp_figs.py [file.ladruno]
"""
import os
import sys

import numpy as np
import h5py
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
# PDMY, identical to bearing_backbone.py
G_REF, B_REF, P_REF, D_EXP = 5.5e4, 1.5e5, 80.0, 0.5
BUCKET = "33017-LadrunoUP[0:0:0]"


def read(path):
    """Nodes, hexes and the per-GP stress/strain of the first and last frame."""
    with h5py.File(path, "r") as f:
        stage = [k for k in f if k.startswith("MODEL_STAGE")][0]
        M = f[stage]["MODEL"]
        coords = np.array(M["NODES"]["COORDINATES"])
        ids = np.array(M["NODES"]["ID"]).ravel()
        order = np.argsort(ids)
        nodes = coords[order]
        nid = ids[order]
        # The element group carries a single CONNECTIVITY dataset of shape
        # (nEl, 1 + nodesPerElement): column 0 is the ELEMENT tag and the rest
        # are its node tags. There is no separate ID dataset here (unlike
        # MODEL/NODES, which does have one).
        emap = {}
        for grp in M["ELEMENTS"]:
            g = M["ELEMENTS"][grp]
            if "CONNECTIVITY" not in g:
                continue
            raw = np.array(g["CONNECTIVITY"])
            emap[grp] = (raw[:, 0], raw[:, 1:])
        # Each result bucket holds ONE DATA dataset of shape
        # (nSteps, nElements, nGP * nComp) — not a group of per-step datasets.
        R = f[stage]["RESULTS"]["ON_ELEMENTS"]
        out = {}
        for nm in ("stress", "strain"):
            if nm not in R:
                continue
            bucket = R[nm][next(iter(R[nm]))]
            d = bucket["DATA"]
            out[nm] = (np.array(d[0]), np.array(d[-1]))
            out[nm + "_nsteps"] = d.shape[0]
            out[nm + "_eid"] = np.array(bucket["ID"]).ravel()
        eid, conn = next(iter(emap.values()))
    # map node ids -> row index
    idx = {int(v): i for i, v in enumerate(nid)}
    hexes = np.vectorize(idx.get)(conn)
    return nodes, hexes, eid, out


def invariants(gp):
    """(nEl, nGP, 6) Voigt -> p' (compression +) and q = sqrt(3 J2)."""
    xx, yy, zz, xy, yz, zx = (gp[..., i] for i in range(6))
    I1 = xx + yy + zz
    p = -I1 / 3.0
    sxx, syy, szz = xx + p, yy + p, zz + p
    J2 = 0.5 * (sxx ** 2 + syy ** 2 + szz ** 2) + xy ** 2 + yz ** 2 + zx ** 2
    return p, np.sqrt(3.0 * J2)


def eq_plastic(sig, eps):
    """Equivalent plastic deviatoric strain per GP; see the module docstring.

    Voigt strain convention: components 3..5 are ENGINEERING shears (gamma),
    so the tensor shears are gamma/2 -- getting this wrong inflates the
    deviatoric measure by 2x on the shear terms.
    """
    p, _ = invariants(sig)
    pp = np.clip(p, 1.0, None)                     # moduli undefined at p'<=0
    G = G_REF * (pp / P_REF) ** D_EXP
    K = B_REF * (pp / P_REF) ** D_EXP
    # deviatoric stress
    s = sig.copy()
    for i in range(3):
        s[..., i] = sig[..., i] + p
    # elastic strain: deviatoric from s/2G, volumetric from p'/K
    ee = np.zeros_like(eps)
    for i in range(3):
        ee[..., i] = s[..., i] / (2.0 * G) - (p / K) / 3.0
    for i in range(3, 6):
        ee[..., i] = s[..., i] / G                 # engineering shear = tau/G
    ep = eps - ee
    # equivalent deviatoric measure of a Voigt strain (engineering shears)
    def eqdev(e):
        m = (e[..., 0] + e[..., 1] + e[..., 2]) / 3.0
        d = [e[..., i] - m for i in range(3)]
        g2 = sum((e[..., i] / 2.0) ** 2 for i in range(3, 6))
        return np.sqrt(2.0 / 3.0 * (sum(x ** 2 for x in d) + 2.0 * g2))
    return eqdev(ep), eqdev(eps)


def slab(nodes, hexes):
    """ONE element layer straddling y = 0, as (x, z) quads.

    The straddle test matters: a filter that merely keeps y > -0.5 retains
    every layer from the symmetry plane outward, they all project onto the same
    (x, z) rectangle, and whichever is drawn last — a far-field element at low
    strain — paints over the bulb. The result looks like a uniformly zero
    field, which is exactly what a first version produced.
    """
    keep, quads = [], []
    for e, h in enumerate(hexes):
        q = nodes[[int(i) for i in h]]
        if q[:, 1].min() < -0.5 or not (q[:, 1].min() < 1e-6 < q[:, 1].max() + 1e-6):
            continue
        keep.append(e)
        quads.append([(q[:, 0].min(), q[:, 2].min()), (q[:, 0].max(), q[:, 2].min()),
                      (q[:, 0].max(), q[:, 2].max()), (q[:, 0].min(), q[:, 2].max())])
    return np.array(keep), quads


def panel(ax, quads, v, title, label, cmap="cividis", vmax=None):
    """One field panel, zoomed to the near field and colour-clipped.

    Both matter. The action is confined to a bulb a couple of B across, so a
    full-domain view renders it as a dot; and these fields are extremely
    peaked — plastic strain reaches 4.4 % under the footing while most of the
    mesh is at ~0 — so a linear scale to the true maximum paints everything
    else flat dark. The colour range is therefore clipped to a high percentile
    and the clipping is stated on the colourbar.
    """
    pc = PolyCollection(quads, array=v, cmap=cmap,
                        edgecolors="#00000018", linewidths=0.2)
    if vmax is not None:
        pc.set_clim(0.0, vmax)
    ax.add_collection(pc)
    ax.set_xlim(-4.5, 4.5); ax.set_ylim(-3.5, 0.4)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]"); ax.set_ylabel("z [m]")
    ax.set_title(title, fontsize=11)
    plt.colorbar(pc, ax=ax, shrink=0.85, label=label)


def main():
    path = os.path.join(HERE, sys.argv[1] if len(sys.argv) > 1
                        else "bearing_gp.ladruno")
    if not os.path.exists(path):
        raise SystemExit(f"{path} missing — run viewer_run.py with "
                         "'-E stress strain' first")
    nodes, hexes, eid, R = read(path)
    if "strain" not in R:
        raise SystemExit("this .ladruno has no `strain` channel — re-run "
                         "viewer_run.py (it records stress+strain now)")
    print(f"{len(nodes)} nodes, {len(hexes)} hex, "
          f"{R['stress_nsteps']} frames")

    s0, s1 = (a.reshape(len(hexes), -1, 6) for a in R["stress"])
    e0, e1 = (a.reshape(len(hexes), -1, 6) for a in R["strain"])
    keep, quads = slab(nodes, hexes)

    p1, q1 = invariants(s1)
    epq1, etq1 = eq_plastic(s1, e1)
    epq0, etq0 = eq_plastic(s0, e0)
    # reference to the gravity frame: the mechanism, not the datum
    dep = np.clip(epq1 - epq0, 0.0, None).mean(axis=1)
    det = np.clip(etq1 - etq0, 0.0, None).mean(axis=1)

    fig, ax = plt.subplots(1, 2, figsize=(13.8, 4.6))
    panel(ax[0], quads, p1.mean(axis=1)[keep],
          "Mean effective stress $p'$ at the last frame\nthe load is carried by "
          "a confined bulb, not a slip surface", "$p'$ [kPa]")
    panel(ax[1], quads, q1.mean(axis=1)[keep],
          "Deviatoric stress $q=\\sqrt{3J_2}$\nhighest under the footing edge, "
          "as a punching mode predicts", "$q$ [kPa]")
    fig.tight_layout()
    pth = os.path.join(OUT, "hypo_gp_stress.png")
    fig.savefig(pth); plt.close(fig)
    print("wrote", pth)

    fig, ax = plt.subplots(1, 2, figsize=(13.8, 4.6))
    panel(ax[0], quads, dep[keep] * 100.0,
          "Equivalent PLASTIC deviatoric strain (DERIVED)\nelastic part removed "
          "with PDMY's own $G(p'),K(p')$ — see note", "$\\varepsilon_q^p$ [%]")
    panel(ax[1], quads, det[keep] * 100.0,
          "Total equivalent deviatoric strain, for comparison\nwhere the two "
          "agree the strain is essentially all plastic", "$\\varepsilon_q$ [%]")
    fig.tight_layout()
    pth = os.path.join(OUT, "hypo_gp_plastic.png")
    fig.savefig(pth); plt.close(fig)
    print("wrote", pth)
    print(f"max eps_q^p = {dep.max() * 100:.3f} %, "
          f"max eps_q = {det.max() * 100:.3f} % "
          f"(plastic fraction {dep.max() / max(det.max(), 1e-12):.2f})")
    print(f"p' range {p1.min():.1f} .. {p1.max():.1f} kPa, "
          f"q range {q1.min():.1f} .. {q1.max():.1f} kPa")


if __name__ == "__main__":
    main()
