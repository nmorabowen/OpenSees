"""WP-138 final q-s figure: every Esmeralda arm to its end + the DP control.

    python plot_qs_final.py            -> qs_final.png (full) and qs_final_zoom.png (s/B <= 0.06)

Circles = where each run stopped (all SANISAND arms: MODE = FLOOR); triangles = the
first converged step whose refusal census contains loadingNonPosH (wall_table.py).
Right panel: push wall clock (each Esmeralda arm alone on its node, 1 thread).
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

from wall_table import load       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
LEGS = [("E_A", "E_A ModifiedEuler (IS 1)", "C3"),
        ("E_B", "E_B SAS-ME (IS 129), TanType 0, TolR 1e-4", "C0"),
        ("E_D", "E_D SAS-ME, TolR 1e-3", "C2"),
        ("E_C2", "E_C2 SAS-ME, TanType 1, Krylov tol x1", "C4"),
        ("E_B16", "E_B16 SAS-ME, B/16 mesh", "C1"),
        ("ctrl_dp38", "control: DruckerPrager 38 deg (local)", "k")]


def main():
    for xmax, name in ((None, "qs_final.png"), (0.06, "qs_final_zoom.png")):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.2))
        for leg, lab, c in LEGS:
            S, sm, codes, first = load(leg)
            sb = [0.0] + [float(r["s_over_B"]) for r in S]
            q = [12.2667] + [float(r["q_kPa"]) for r in S]
            wall = [0.0] + [float(r["wall_total_s"]) / 3600.0 for r in S]
            ls = "--" if leg == "ctrl_dp38" else "-"
            ax1.plot(sb, q, ls, color=c, lw=1.5, label=lab)
            ax1.plot(sb[-1], q[-1], "o", color=c, ms=5)
            if first:
                ax1.plot(sb[first], q[first], "v", color=c, ms=8, mfc="none", mew=1.5)
            ax2.plot(sb, wall, ls, color=c, lw=1.5, label=lab)
        ax1.axvspan(0.026, 0.041, color="grey", alpha=0.15, label="TIMs' ModifiedEuler wall, s/B 0.026-0.041")
        ax1.set_xlabel("s / B")
        ax1.set_ylabel("q = footing force / B  [kPa]")
        ax1.set_title("q-s to the end of every arm (o = FLOOR stop, v = first loadingNonPosH)")
        ax2.set_xlabel("s / B")
        ax2.set_ylabel("push wall clock [h]")
        ax2.set_title("cost (Esmeralda build 7936ed6e0, 1 thread per arm)")
        for ax in (ax1, ax2):
            ax.grid(alpha=0.3)
            if xmax:
                ax.set_xlim(0, xmax)
            ax.legend(fontsize=7)
        if xmax:
            ax1.set_ylim(0, 1100)
            ax2.set_ylim(0, 13)
        fig.tight_layout()
        fig.savefig(os.path.join(HERE, name), dpi=130)
        print("wrote", name)


if __name__ == "__main__":
    main()
