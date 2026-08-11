"""Re-derive the verdict of every `h20_prandtl.py` leg from its CSV alone.

The runner prints its verdict once, at the end.  This reads the CSVs back, so
a leg that is still running (or that a wall-clock cap ended) can be scored on
exactly the same rules -- the D16 2 %-drop truncation and the
tail-slope-under-2 %-of-initial plateau test -- without re-running anything.

Run:  py -3.12 h20_summary.py [glob]
"""
import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DROP_TRUNCATE = 0.02
B_FOOT = 2.0
# q0 * N_q(phi_ps) at phi_txc = 20 deg, q0 = 10 kPa: alpha = 0.148583084,
# phi_ps = 27.4701614 deg (Chen & Han plane-strain match), N_q = 13.8907010.
Q_EXACT = 138.90700970835823


def score(path):
    a = np.genfromtxt(path, delimiter=",", names=True)
    if a.size < 5:
        return None
    s, q = np.atleast_1d(a["s_m"]), np.atleast_1d(a["q_kPa"])
    ds = np.atleast_1d(a["ds_mm"])
    rel = np.atleast_1d(a["relaxed"])
    wall = float(np.atleast_1d(a["wall_s"])[-1])
    n_all = len(s)
    peak = np.maximum.accumulate(q)
    bad = np.where(q < (1.0 - DROP_TRUNCATE) * peak)[0]
    ntr = 0
    if len(bad):
        ntr = len(q) - bad[0]
        s, q = s[:bad[0]], q[:bad[0]]
    msk = s >= 0.9 * s[-1]
    t_last = np.polyfit(s[msk], q[msk], 1)[0] if msk.sum() > 2 else float("nan")
    n0 = max(4, len(s) // 50)
    t_init = np.polyfit(s[:n0], q[:n0], 1)[0]
    return dict(name=os.path.basename(path)[4:-4], n=n_all, ntr=ntr,
                qmax=float(q.max()), ratio=float(q.max()) / Q_EXACT,
                s_end=float(s[-1]) / B_FOOT,
                s_at_max=float(s[int(q.argmax())]) / B_FOOT,
                tail=100 * t_last / t_init,
                plateau=bool(abs(t_last) < 0.02 * abs(t_init)),
                wall=wall, ds_min=float(ds.min()), nrel=int(rel.sum()))


def main():
    pat = sys.argv[1] if len(sys.argv) > 1 else "h20_*.csv"
    rows = [r for r in (score(p) for p in sorted(glob.glob(os.path.join(HERE, pat))))
            if r]
    hdr = (f"{'leg':>28} {'steps':>6} {'q_max':>8} {'/exact':>7} "
           f"{'s@max/B':>8} {'s_end/B':>8} {'tail %':>8} {'plateau':>8} "
           f"{'trunc':>6} {'relax':>6} {'ds_min':>8} {'wall_s':>8}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['name']:>28} {r['n']:6d} {r['qmax']:8.2f} {r['ratio']:7.4f} "
              f"{r['s_at_max']:8.4f} {r['s_end']:8.4f} {r['tail']:8.2f} "
              f"{'yes' if r['plateau'] else 'NO':>8} {r['ntr']:6d} "
              f"{r['nrel']:6d} {r['ds_min']:8.4f} {r['wall']:8.0f}")
    print(f"\noracle q0*N_q = {Q_EXACT:.4f} kPa (exact; phi_txc = 20 deg, "
          f"Chen & Han plane-strain match, weightless)")


if __name__ == "__main__":
    main()
