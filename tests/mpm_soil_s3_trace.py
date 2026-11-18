"""WP-232k S3 diagnosis: the HS strip footing with a per-step trace (where the soil under the
footing sits relative to the grid), to explain the r3 fbar collapse at s = 0.096 (s/h = 0.77).

Per accepted push step: q, the ladder rung that accepted it, and, for the particles under the
footing (x < B/2 = 1): the top particle centre y, and how many particle centres are still in
the top cell row (y > -h). Same deck and ladder as ladruno_mpm_soil_s_harness.s3_run.

    python tests/mpm_soil_s3_trace.py --r 2 --al fbar --s-end 0.2 --out <json>
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ladruno_mpm_soil_s_harness as S  # noqa: E402


def under_footing(ops, h):
    x = S._l(ops.mpm("particleData", "x"))
    ys = [x[2 * i + 1] for i in range(len(x) // 2) if x[2 * i] < S.B_HALF]
    return {"top_y": max(ys), "top_y_over_h": max(ys) / h, "in_top_row": sum(1 for y in ys if y > -h)}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--r", type=int, default=2)
    ap.add_argument("--al", default="fbar")
    ap.add_argument("--s-end", type=float, default=0.2)
    ap.add_argument("--out")
    a = ap.parse_args(argv)
    from _testbed import ops
    ref = S.load(f"s3_fe_hs_n{4 * 2 ** (a.r - 1)}.json")
    h = 0.5 / 2 ** (a.r - 1)
    w0 = time.time()
    S.s3_build(ops, ref["decks"]["hs"]["cmd"], a.r, a.al)
    st = {"h": 1.0, "failed": 0, "relaxed": 0, "halvings": 0}
    out = {"r": a.r, "antilocking": a.al, "h": h, "build": ops.ladrunoBuild(), "rows": [],
           "stopped": None}
    for k in range(1, 13):
        if k == 11:
            ops.mpm("stage", 1, 1)
        if not S._ladder_step(ops, float(k), st):
            out["stopped"] = f"t = {k - 1} (staging)"
            break
    out["before_push"] = under_footing(ops, h)
    n_push = int(round(a.s_end / 0.002))
    for k in range(1, n_push + 1 if out["stopped"] is None else 0):
        f0, r0, h0 = st["failed"], st["relaxed"], st["halvings"]
        if not S._ladder_step(ops, 12.0 + k, st):
            out["stopped"] = f"s = {0.002 * (k - 1):.4f}"
            break
        q = -S._l(ops.mpm("condForce", "footing"))[1] / S.B_HALF
        row = {"s": round(0.002 * k, 6), "s_over_h": round(0.002 * k / h, 4), "q": q,
               "failed": st["failed"] - f0, "relaxed": st["relaxed"] - r0,
               "halvings": st["halvings"] - h0, **under_footing(ops, h)}
        out["rows"].append(row)
        print(f"s {row['s']:.3f} s/h {row['s_over_h']:.3f} q {q:9.2f}  top_y/h {row['top_y_over_h']:+.3f}"
              f"  in_top_row {row['in_top_row']:3d}  failed {row['failed']} relaxed {row['relaxed']}"
              f" halvings {row['halvings']}", flush=True)
    out["wall_s"] = time.time() - w0
    if a.out:
        with open(a.out, "w", encoding="utf-8") as f:
            json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
