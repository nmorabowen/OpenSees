"""WP-150 R2 analyser (memo section 9.1 steps 3-6): mesh study of the WP-138 footing legs.

    python r2_analysis.py LABEL=RUNDIR [LABEL=RUNDIR ...] [--family b4,b8,b16] [--sb 0.005,0.01,...] [--png out.png]

RUNDIR is a footing_ab.py output directory (steps.csv + ckpt/field_*.npz).
Reports:
  1. q at matched s/B for every leg. For the mesh family (coarse -> fine) it gives the successive differences, their
     ratio (< 1 = contracting) and, when contracting, a Richardson-type estimate of the h -> 0 limit.
  2. Band geometry from the last two saved checkpoints of each leg (incremental gamma_max = sqrt((exx-eyy)^2+gxy^2)):
     for depth rows under the RIGHT footing edge,
       - x_band: the x of the peak incremental gamma_max in 0 < x < 3 m;
       - FWHM: the full width at half maximum of p = (gamma_inc - row median)^+ around the peak (linear crossings).
         It is the robust width, quoted as FWHM/h: ~1 = a one-element band;
       - w2: the ADR-90 section 7.3 width, sqrt(12 Var) with Var = sum p[(x - xbar)^2 + c^2/12] / sum p, in a
         +-3 GP-column window (c = the local GP column spacing). It is h-bounded below and sensitive to the shoulder;
       - the band's inclination from vertical, by a least-squares fit of x_band(y) over -2.25 < y < -0.2 m.
         Negative = leaning INWARD (toward the footing centre, the general-shear wedge); ~0 = vertical, the
         mesh-aligned punching seen at b8/b16; positive = outward.
A leg whose band paths follow its own mesh lines (vertical at b8/b16, leaning with the columns on shear:15, erratic
on jitter) is mesh-biased. A path that is the same across variants is physical.
"""
import argparse, csv, glob, os, re, sys
import numpy as np

B = 1.5


def curve(d):
    rows = list(csv.DictReader(open(os.path.join(d, "steps.csv"))))
    return (np.array([float(r["s_over_B"]) for r in rows]), np.array([float(r["q_kPa"]) for r in rows]))


def ckpts(d):
    fs = glob.glob(os.path.join(d, "ckpt", "field_step*.npz"))
    fs.sort(key=lambda f: int(re.search(r"step(\d+)", f).group(1)))
    last = os.path.join(d, "ckpt", "field_last_converged.npz")
    if os.path.exists(last):
        fs.append(last)
    return fs


def gmax(e):
    return np.sqrt((e[:, 0] - e[:, 1]) ** 2 + e[:, 2] ** 2)


def band_geometry(f0, f1):
    a, b = np.load(f0, allow_pickle=True), np.load(f1, allow_pickle=True)
    dg = gmax(b["eps"] - a["eps"])
    x, y = b["gx"], b["gy"]
    fine = (y > -2.3) & (y < 0) & (x > 0) & (x < 3.0)
    ys = np.unique(np.round(y[fine], 3))
    dy = np.median(np.diff(ys)) if len(ys) > 1 else 0.1
    xs = np.unique(np.round(x[fine], 3))
    cx = np.median(np.diff(xs)) if len(xs) > 1 else 0.1
    out = []
    for y0 in np.arange(-0.2, -2.3, -0.25):
        dmin = float(np.min(np.abs(ys - y0)))
        m = fine & (np.abs(y - y0) <= dmin + 0.25 * dy + 1e-9)     # the nearest GP row (a band of rows on jitter)
        if m.sum() < 5:
            continue
        xx, gg = x[m], dg[m]
        o = np.argsort(xx); xx, gg = xx[o], gg[o]
        k = int(np.argmax(gg))
        xb = float(xx[k])
        base = float(np.median(gg))
        p_all = np.clip(gg - base, 0.0, None)
        # w2 (ADR-90 7.3) in a +-3 GP-column window: tail-robust against the broad shoulder
        w = np.abs(xx - xb) <= 3.0 * cx + 1e-9
        p = p_all[w]
        if p.sum() <= 0:
            continue
        xbar = float(np.sum(p * xx[w]) / p.sum())
        var = float(np.sum(p * ((xx[w] - xbar) ** 2 + cx ** 2 / 12.0)) / p.sum())
        # FWHM of the background-removed profile around the peak (linear crossings)
        half = 0.5 * p_all[k]
        il = k
        while il > 0 and p_all[il - 1] >= half:
            il -= 1
        ir = k
        while ir < len(xx) - 1 and p_all[ir + 1] >= half:
            ir += 1
        def cross(i0, i1):
            y0_, y1_ = p_all[i0], p_all[i1]
            return xx[i0] + (half - y0_) * (xx[i1] - xx[i0]) / (y1_ - y0_) if y1_ != y0_ else xx[i0]
        xl = cross(il - 1, il) if il > 0 else xx[il]
        xr = cross(ir + 1, ir) if ir < len(xx) - 1 else xx[ir]
        out.append(dict(y=float(y0), x_band=xb, w2=float(np.sqrt(12 * var)), fwhm=float(xr - xl), c=float(cx)))
    fit = [o for o in out if -2.25 <= o["y"] <= -0.2]
    incl = float("nan")
    if len(fit) >= 3:
        sl = np.polyfit([o["y"] for o in fit], [o["x_band"] for o in fit], 1)[0]   # dx/dy, y negative downward
        incl = float(np.degrees(np.arctan(-sl)))    # x decreasing with depth (inward) -> negative
    return out, incl, float(b["s_over_B"]), float(a["s_over_B"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("legs", nargs="+")
    ap.add_argument("--family", default="", help="comma list of labels, coarse to fine (e.g. b4,b8,b16)")
    ap.add_argument("--sb", default="0.002,0.005,0.01,0.02,0.03,0.05,0.075,0.10,0.15")
    ap.add_argument("--h", default="", help="label=h pairs for w2/h, e.g. b4=0.375,b8=0.1875")
    ap.add_argument("--png", default="")
    a = ap.parse_args()
    legs = dict(l.split("=", 1) for l in a.legs)
    hmap = dict((k, float(v)) for k, v in (p.split("=") for p in a.h.split(",") if p)) if a.h else {}
    grid = [float(v) for v in a.sb.split(",")]
    C = {k: curve(v) for k, v in legs.items()}

    print("## q (kPa) at matched s/B\n")
    print("| leg | reached s/B | " + " | ".join(f"{s:g}" for s in grid) + " |")
    print("|---|---|" + "---|" * len(grid))
    for k, (sb, q) in C.items():
        cells = [f"{np.interp(s, sb, q):.1f}" if s <= sb[-1] else "—" for s in grid]
        print(f"| {k} | {sb[-1]:.4f} | " + " | ".join(cells) + " |")
    fam = [f for f in a.family.split(",") if f]
    if len(fam) >= 2:
        print(f"\n## mesh family {' → '.join(fam)}: successive differences (%) and contraction\n")
        print("| s/B | " + " | ".join(f"{fam[i]}→{fam[i+1]} %" for i in range(len(fam) - 1))
              + (" | ratio | Richardson h→0 kPa |" if len(fam) == 3 else " |"))
        print("|---|" + "---|" * (len(fam) - 1) + ("---|---|" if len(fam) == 3 else ""))
        for s in grid:
            if any(s > C[f][0][-1] for f in fam):
                continue
            qs = [float(np.interp(s, *C[f])) for f in fam]
            d = [qs[i + 1] - qs[i] for i in range(len(qs) - 1)]
            cells = [f"{100 * d[i] / qs[i]:+.2f}" for i in range(len(d))]
            tail = ""
            if len(fam) == 3:
                r = d[1] / d[0] if d[0] != 0 else float("nan")
                rich = qs[2] + d[1] * r / (1 - r) if 0 < r < 1 else float("nan")
                tail = f" | {r:.2f} | {rich:.1f}" if np.isfinite(r) else " | — | —"
            print(f"| {s:g} | " + " | ".join(cells) + tail + " |")

    print("\n## band geometry under the right footing edge (last two saved checkpoints of each leg)\n")
    print("| leg | interval s/B | inclination from vertical ° (− = inward) | x_band at y ≈ −0.2 / −1.2 / −2.2 m | FWHM median (m) | FWHM/h | w2 median (m, ±3 GP cols) |")
    print("|---|---|---|---|---|---|---|")
    paths = {}
    for k, d in legs.items():
        fs = ckpts(d)
        if len(fs) < 2:
            print(f"| {k} | no checkpoint pair | | | | |")
            continue
        out, incl, s1, s0 = band_geometry(fs[-2], fs[-1])
        paths[k] = out
        pick = {round(o["y"], 2): o["x_band"] for o in out}
        xb = " / ".join(f"{pick.get(v, float('nan')):.2f}" for v in (-0.2, -1.2, -2.2))
        w2 = np.median([o["w2"] for o in out]) if out else float("nan")
        fw = np.median([o["fwhm"] for o in out]) if out else float("nan")
        h = hmap.get(k, float("nan"))
        print(f"| {k} | {s0:.4f}→{s1:.4f} | {incl:+.1f} | {xb} | {fw:.3f} | {fw / h if np.isfinite(h) else float('nan'):.2f} | {w2:.3f} |")
    if a.png and paths:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
        for k, (sb, q) in C.items():
            ax[0].plot(sb, q, label=k)
        ax[0].set_xlabel("s/B"); ax[0].set_ylabel("q (kPa)"); ax[0].grid(alpha=.3); ax[0].legend(fontsize=7)
        for k, out in paths.items():
            ax[1].plot([o["x_band"] for o in out], [o["y"] for o in out], "o-", ms=3, label=k)
        ax[1].axvline(0.75, color="k", lw=.5); ax[1].set_xlabel("x of peak incremental γ (m)"); ax[1].set_ylabel("y (m)")
        ax[1].set_title("band path, right edge (footing edge at x = 0.75)", fontsize=9); ax[1].grid(alpha=.3); ax[1].legend(fontsize=7)
        fig.tight_layout(); fig.savefig(a.png, dpi=130)
        print(f"\nwrote {a.png}")


if __name__ == "__main__":
    main()
