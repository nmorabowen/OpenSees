"""WP-138: per-point substep distribution between consecutive checkpoints.

Reads runs/<leg>/ckpt/field_*.npz (cumulative per-point census `stats`,
column 2 = substeps in both substepStats and sasStats) and, for each interval
between two checkpoints, reports: total substeps, the share carried by the
RING under two definitions --
  ringP  : points with p' < 10 kPa at the END of the interval;
  ringG  : the TOP element row (Gauss point y > -B/8) with 0.5 <= |x|/B <= 2;
-- the top-1 / top-10 / top-100 point shares, and a log-binned histogram of
substeps per point. Also the per-STEP top-1 share from steps.csv
(sub_step_maxpt / sub_step_total).

    python census_dist.py LEG [LEG ...]   (writes runs/<leg>/census_dist.json)
"""
import csv
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
B = 1.5


def main():
    for leg in sys.argv[1:]:
        rd = os.path.join(HERE, "runs", leg)
        files = sorted(glob.glob(os.path.join(rd, "ckpt", "field_step*.npz")))
        fl = os.path.join(rd, "ckpt", "field_last_converged.npz")
        if os.path.exists(fl):
            files.append(fl)
        out = []
        prev = None
        for f in files:
            z = np.load(f)
            st = z["stats"]
            if st.ndim != 2:
                continue
            cur = dict(step=int(z["step"]), sB=float(z["s_over_B"]), sub=st[:, 2].copy(),
                       p=z["p"], gx=z["gx"], gy=z["gy"])
            if prev is not None and cur["step"] > prev["step"]:
                d = cur["sub"] - prev["sub"]
                tot = float(d.sum())
                if tot <= 0:
                    prev = cur
                    continue
                srt = np.sort(d)[::-1]
                ringP = cur["p"] < 10.0
                ax = np.abs(cur["gx"]) / B
                ringG = (cur["gy"] > -B / 8) & (ax >= 0.5) & (ax <= 2.0)
                edges = [0, 1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e9]
                hist = np.histogram(d, bins=edges)[0].tolist()
                out.append(dict(
                    from_step=prev["step"], to_step=cur["step"],
                    sB_from=prev["sB"], sB_to=cur["sB"], total=tot,
                    per_step_mean=tot / (cur["step"] - prev["step"]),
                    ringP_n=int(ringP.sum()), ringP_share=float(d[ringP].sum() / tot),
                    ringG_n=int(ringG.sum()), ringG_share=float(d[ringG].sum() / tot),
                    top1=float(srt[0] / tot), top10=float(srt[:10].sum() / tot),
                    top100=float(srt[:100].sum() / tot),
                    median_pt=float(np.median(d)), p99_pt=float(np.percentile(d, 99)),
                    hist_edges=edges, hist=hist))
            prev = cur
        # per-step top-1 share from steps.csv
        sp = os.path.join(rd, "steps.csv")
        per_step = []
        if os.path.exists(sp):
            for r in csv.DictReader(open(sp, newline="")):
                t = float(r["sub_step_total"] or 0)
                if t > 0:
                    per_step.append((int(r["step"]), float(r["s_over_B"]),
                                     float(r["sub_step_maxpt"]) / t))
        res = dict(leg=leg, intervals=out,
                   per_step_top1_max=max((x[2] for x in per_step), default=None),
                   per_step_top1_median=(float(np.median([x[2] for x in per_step]))
                                         if per_step else None))
        json.dump(res, open(os.path.join(rd, "census_dist.json"), "w"), indent=1)
        print(f"== {leg}  per-step top-1 share: median {res['per_step_top1_median']}, "
              f"max {res['per_step_top1_max']}")
        print("steps      s/B          total/step  ringP(n)       ringG(n)       top1    top10   top100  hist[0,1,10,1e2,1e3,1e4,1e5,1e6)")
        for o in out:
            print(f"{o['from_step']:3d}-{o['to_step']:3d}  {o['sB_from']:.4f}-{o['sB_to']:.4f}  "
                  f"{o['per_step_mean']:11.3e}  {o['ringP_share']:.3f}({o['ringP_n']:4d})  "
                  f"{o['ringG_share']:.3f}({o['ringG_n']:3d})  {o['top1']:.4f}  "
                  f"{o['top10']:.4f}  {o['top100']:.4f}  {o['hist']}")


if __name__ == "__main__":
    main()
