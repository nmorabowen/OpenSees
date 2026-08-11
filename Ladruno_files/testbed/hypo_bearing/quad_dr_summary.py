"""quad_dr_summary -- roll up the DR legs and overlay them on the STATIC legs.

Two things the per-leg log cannot do:

1.  **The leg table with termination modes**, so no number is ever quoted
    without the mode that produced it (notes 81/82's standing rule).

2.  **The point-by-point agreement between the tangent-free path and the
    Newton path.**  Note 82's static legs are in-tree as CSVs on the SAME mesh,
    element, material and loading, so the DR curve can be compared against the
    static curve at MATCHED settlement rather than only at the headline.  That
    is the validation that actually matters: a headline that agrees while the
    curve does not would mean the two paths found different states and met at
    one point by luck.  The static curve is interpolated onto the DR
    settlements over the overlap only, and both the RMS and the worst
    disagreement are printed.

Run:  py -3.12 quad_dr_summary.py
"""
import glob
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
Q_EXACT = 138.907          # q0 * N_q, phi_ps from the Chen & Han match

# note 82 sec.5's static legs on the same mesh family, for the overlay
STATIC = {"h8bbar_h5": ("qpd_h8bbar_h5.csv", "FLOOR/BUDGET: capacity 0.9977"),
          "h20uri_h5": ("qpd_h20uri_h5.csv", "FLOOR: allowance 0.6846"),
          "h20std_h5": ("qpd_h20std_h5.csv", "FLOOR: allowance 0.4587")}


def read_dr(path):
    a = np.genfromtxt(path, delimiter=",", names=True)
    if a.ndim == 0:
        a = a.reshape(1)
    return a


def read_static(name):
    p = os.path.join(HERE, STATIC[name][0])
    if not os.path.exists(p):
        return None
    a = np.genfromtxt(p, delimiter=",", names=True)
    return a


def tail_stats(s, q):
    msk = s >= 0.9 * s[-1]
    t_last = np.polyfit(s[msk], q[msk], 1)[0] if msk.sum() > 2 else float("nan")
    n0 = max(4, len(s) // 50)
    t_init = np.polyfit(s[:n0], q[:n0], 1)[0]
    return t_last, 100.0 * t_last / t_init


def main():
    rows = []
    for p in sorted(glob.glob(os.path.join(HERE, "dr_*.csv"))):
        tag = os.path.basename(p)[3:-4]
        try:
            a = read_dr(p)
        except Exception as exc:
            print(f"  (skipped {tag}: {exc})")
            continue
        if a.size < 3:
            continue
        s, q = np.atleast_1d(a["s_m"]), np.atleast_1d(a["q_kPa"])
        it = np.atleast_1d(a["dr_iters"])
        st = np.atleast_1d(a["settled"])
        t_last, t_frac = tail_stats(s, q)
        rows.append(dict(tag=tag, n=len(s), qmax=q.max(),
                         ratio=q.max() / Q_EXACT, s_end=s[-1] / 2.0,
                         t_frac=t_frac, it_max=it.max(), it_tot=it.sum(),
                         unsettled=int((st == 0).sum())))

    hdr = (f"{'DR leg':>30} {'inc':>5} {'q_max':>8} {'q/exact':>8} "
           f"{'s_end/B':>8} {'tail %':>9} {'maxDRit':>8} {'unsettled':>9}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['tag']:>30} {r['n']:>5} {r['qmax']:8.2f} {r['ratio']:8.4f} "
              f"{r['s_end']:8.4f} {r['t_frac']:9.3f} {int(r['it_max']):8d} "
              f"{r['unsettled']:9d}")

    print("\nDR vs the STATIC Newton path, at MATCHED settlement "
          "(static curve interpolated onto the DR settlements, overlap only):")
    hdr2 = (f"{'DR leg':>30} {'vs static':>14} {'overlap s/B':>16} {'pts':>5} "
            f"{'RMS %':>8} {'worst %':>8}")
    print(hdr2)
    print("-" * len(hdr2))
    for r in rows:
        for name in STATIC:
            if not r["tag"].startswith("h8bbar_h5") and \
               not r["tag"].startswith("h20uri_h5") and \
               not r["tag"].startswith("h20std_h5"):
                continue
            if not r["tag"].startswith(name):
                continue
            sa = read_static(name)
            if sa is None:
                continue
            a = read_dr(os.path.join(HERE, f"dr_{r['tag']}.csv"))
            s, q = np.atleast_1d(a["s_m"]), np.atleast_1d(a["q_kPa"])
            ss, qs = sa["s_m"], sa["q_kPa"]
            lo, hi = max(s.min(), ss.min()), min(s.max(), ss.max())
            m = (s >= lo) & (s <= hi)
            if m.sum() < 2:
                print(f"{r['tag']:>30} {name:>14} {'(no overlap)':>16}")
                continue
            qi = np.interp(s[m], ss, qs)
            d = 100.0 * (q[m] - qi) / qi
            print(f"{r['tag']:>30} {name:>14} "
                  f"{lo / 2:.5f}-{hi / 2:.5f} {int(m.sum()):>5} "
                  f"{np.sqrt((d ** 2).mean()):8.3f} "
                  f"{d[np.abs(d).argmax()]:+8.3f}")


if __name__ == "__main__":
    main()
