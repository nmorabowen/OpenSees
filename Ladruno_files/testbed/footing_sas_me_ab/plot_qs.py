"""WP-138: q-s plot (and cost plot) from runs/<leg>/steps.csv.

    python plot_qs.py --out qs.png leg1[:label] leg2[:label] ...
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def load(leg):
    p = os.path.join(HERE, "runs", leg, "steps.csv")
    rows = list(csv.DictReader(open(p, newline="")))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(HERE, "qs.png"))
    ap.add_argument("--xmax", type=float, default=None)
    ap.add_argument("legs", nargs="+")
    a = ap.parse_args()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))
    for spec in a.legs:
        leg, _, lab = spec.partition(":")
        rows = load(leg)
        sb = [0.0] + [float(r["s_over_B"]) for r in rows]
        q = [12.2667] + [float(r["q_kPa"]) for r in rows]
        wall = [0.0] + [float(r["wall_total_s"]) / 3600.0 for r in rows]
        l, = ax1.plot(sb, q, lw=1.6, label=lab or leg)
        ax1.plot(sb[-1], q[-1], "x", color=l.get_color(), ms=9)
        ax2.plot(sb, wall, lw=1.6, color=l.get_color(), label=lab or leg)
    ax1.set_xlabel("s / B")
    ax1.set_ylabel("q = footing contact force / B  [kPa]")
    ax1.set_title("load-settlement (x = where the run stopped)")
    ax2.set_xlabel("s / B")
    ax2.set_ylabel("push wall clock [h]")
    ax2.set_title("cost")
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3)
        if a.xmax:
            ax.set_xlim(0, a.xmax)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
