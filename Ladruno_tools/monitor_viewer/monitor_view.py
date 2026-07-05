"""
monitor_view.py  —  CLI plotter for a Ladruno Analysis-monitor sink.

Static:
    python monitor_view.py sink.h5                 # show all channels vs TIME
    python monitor_view.py sink.h5 --save out.png  # headless -> PNG
    python monitor_view.py sink.h5 --x step        # x-axis = STEP index
    python monitor_view.py sink.h5 --channels node2.disp.dof1,node5.disp.dof1

Live tail (auto-refresh while a run is streaming):
    python monitor_view.py sink.h5 --watch         # redraw every --interval s
                                                   # until the writer stops growing

The live web dashboard (monitor_server.py) is the richer option; this is the
zero-web-stack, drop-into-a-notebook plotter.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

from monitor_reader import MonitorReader


def _select(cols, wanted):
    if not wanted:
        return list(range(len(cols)))
    idx = []
    for w in wanted:
        if w in cols:
            idx.append(cols.index(w))
        else:
            print(f"  warning: channel '{w}' not in sink; skipped", file=sys.stderr)
    return idx or list(range(len(cols)))


def main(argv=None):
    ap = argparse.ArgumentParser(description="Plot a Ladruno monitor sink (.h5).")
    ap.add_argument("sink", help="path to the monitor SWMR-HDF5 sink")
    ap.add_argument("--watch", action="store_true",
                    help="live tail: redraw as new frames arrive")
    ap.add_argument("--interval", type=float, default=0.5,
                    help="poll/redraw interval in seconds (--watch)")
    ap.add_argument("--x", choices=["time", "step"], default="time",
                    help="x-axis: pseudo-time (default) or step index")
    ap.add_argument("--channels", default="",
                    help="comma-separated channel labels (default: all)")
    ap.add_argument("--save", default="",
                    help="write a PNG instead of showing a window")
    args = ap.parse_args(argv)

    import matplotlib
    if args.save and not args.watch:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    reader = MonitorReader(args.sink)
    m = reader.meta()
    cols = m["columns"]
    sel = _select(cols, [c for c in args.channels.split(",") if c])
    ylabel = f"value ({m['units']})" if m["units"] else "value"
    xlabel = "pseudo-time" if args.x == "time" else "step index"

    def load():
        d = reader.frames_since(0)
        x = d["t"] if args.x == "time" else d["step"]
        return x, d["rows"], d["n"]

    x, rows, n = load()

    fig, ax = plt.subplots(figsize=(9, 4.4), dpi=110)
    lines = {}
    for j in sel:
        (ln,) = ax.plot(x, [r[j] for r in rows], lw=1.5, label=cols[j])
        lines[j] = ln
    ax.axhline(0, color="0.75", lw=0.8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right", fontsize=9)

    def title(n, live=False):
        tag = "LIVE " if live else ""
        ax.set_title(f"{tag}monitor — {os.path.basename(args.sink)}  ({n} frames)")

    if not args.watch:
        title(n)
        fig.tight_layout()
        if args.save:
            fig.savefig(args.save)
            print("wrote", args.save)
        else:
            plt.show()
        return 0

    # --watch: poll until the writer stops growing for a few intervals.
    plt.ion()
    title(n, live=True)
    fig.tight_layout()
    plt.show(block=False)
    stagnant = 0
    while stagnant < 6:
        plt.pause(args.interval)
        if not plt.fignum_exists(fig.number):
            break
        x, rows, n2 = load()
        for j in sel:
            lines[j].set_data(x, [r[j] for r in rows])
        ax.relim()
        ax.autoscale_view()
        # Only count stagnation once at least one frame exists, so launching the
        # tail before the writer starts waits for it rather than giving up.
        stagnant = (stagnant + 1) if (n2 == n and n2 > 0) else 0
        n = n2
        title(n, live=(stagnant == 0))
        fig.canvas.draw_idle()
    title(n, live=False)
    print(f"done: {n} frames (writer idle)")
    if plt.fignum_exists(fig.number):
        plt.ioff()
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
