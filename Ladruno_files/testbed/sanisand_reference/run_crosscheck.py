"""WP-134 validation (3): reference vs the fork's C++ on benign states.

    python Ladruno_files/testbed/sanisand_reference/run_crosscheck.py [--deltas 1e-5 1e-4 1e-3]

(any CPython with numpy + scipy; the C++ side runs in a 3.12 subprocess, see
Ladruno_scripts/sanisand_reference/cxx.py).  Writes out/crosscheck.json and
out/crosscheck.md."""
import argparse
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "Ladruno_scripts"))

import numpy as np  # noqa: E402

from sanisand_reference import CAMPAIGN  # noqa: E402
from sanisand_reference import cxx  # noqa: E402
from sanisand_reference.crosscheck import (benign_states, job_of, probes,  # noqa: E402
                                           rel_incr_diff, run_reference, uw_variants)
from sanisand_reference.model import t2v  # noqa: E402

OUT = os.path.join(HERE, "out")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deltas", type=float, nargs="*", default=[1e-5, 1e-4, 1e-3])
    ap.add_argument("--rtol", type=float, default=1e-10)
    ap.add_argument("--protos", nargs="*", default=["ME8", "ME", "RK45"])
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    P = CAMPAIGN
    cases, labels = [], []
    for sname, st in benign_states(P):
        for d in args.deltas:
            for pname, de in probes(d).items():
                cases.append((st, de))
                labels.append(dict(state=sname, probe=pname, delta=d,
                                   start=dict(sigma=t2v(st.sigma).tolist(),
                                              alpha=t2v(st.alpha).tolist(),
                                              z=t2v(st.z).tolist(), e=st.e)))
    variants = uw_variants()
    ref = run_reference(cases, P, variants, rtol=args.rtol)
    jobs = [job_of(st, de, pr) for pr in args.protos for (st, de) in cases]
    craw = cxx.run_jobs(jobs, P.as_opensees())
    c = {pr: craw[k * len(cases):(k + 1) * len(cases)] for k, pr in enumerate(args.protos)}
    rows = []
    for i, lab in enumerate(labels):
        row = dict(lab)
        s0 = lab["start"]
        for pr in args.protos:
            cx = c[pr][i]
            row[pr] = dict(rc=cx.get("rc"), substeps=cx.get("substeps"),
                           forced=cx.get("forced"))
            for vn in variants:
                rs = ref[vn][i]
                ds, dsa = rel_incr_diff(s0, rs, cx, "sigma")
                da, daa = rel_incr_diff(s0, rs, cx, "alpha")
                row[pr][vn] = dict(dsig_rel=ds, dsig_abs=dsa, dalpha_rel=da,
                                   dalpha_abs=daa)
        for vn in variants:
            row.setdefault("ref", {})[vn] = dict(status=ref[vn][i]["status"],
                                                 segments=ref[vn][i]["segments"],
                                                 reseats=ref[vn][i]["reseats"],
                                                 f_end=ref[vn][i]["f_end"],
                                                 max_rho_b=ref[vn][i]["max_rho_b"],
                                                 negh=ref[vn][i]["negh"])
        # effect of each UW toggle alone: || ref(uw_me) - ref(uw_me with k = paper) ||
        for vn in variants:
            if vn.startswith("uw_me-"):
                es, _ = rel_incr_diff(s0, ref[vn][i], ref["uw_me"][i], "sigma")
                row.setdefault("toggle_effect", {})[vn[6:]] = es
        rows.append(row)
    json.dump(dict(rows=rows, cxx=c, ref=ref), open(os.path.join(OUT, "crosscheck.json"), "w"),
              indent=1, default=float)
    # --- markdown summary ---------------------------------------------------
    L = []
    L.append("| state | probe | delta | ME8 substeps | paper vs ME8 dsig | uw_me vs ME8 dsig | uw_me vs ME8 dalpha | uw_me UW-rule h<0 | largest single-toggle effect (from uw_me) | uw vs RK45 dsig | uw_me vs ME(campaign) dsig |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        te = r.get("toggle_effect", {})
        big = max(te.items(), key=lambda kv: kv[1] if math.isfinite(kv[1]) else -1) if te else ("-", float("nan"))
        L.append(f"| {r['state']} | {r['probe']} | {r['delta']:.0e} | {r['ME8']['substeps']:.0f} | "
                 f"{r['ME8']['paper']['dsig_rel']:.2e} | {r['ME8']['uw_me']['dsig_rel']:.2e} | "
                 f"{r['ME8']['uw_me']['dalpha_rel']:.2e} | {r['ref']['uw_me']['negh']} | {big[0]} {big[1]:.1e} | "
                 f"{r.get('RK45', {}).get('uw', {}).get('dsig_rel', float('nan')):.2e} | "
                 f"{r.get('ME', {}).get('uw_me', {}).get('dsig_rel', float('nan')):.2e} |")
    open(os.path.join(OUT, "crosscheck.md"), "w").write("\n".join(L) + "\n")
    # console digest
    for d in args.deltas:
        sub = [r for r in rows if r["delta"] == d]
        for vn, pr in [("paper", "ME8"), ("uw_me", "ME8"), ("uw", "RK45"), ("uw_me", "ME")]:
            if pr not in sub[0]:
                continue
            vals = [r[pr][vn]["dsig_rel"] for r in sub]
            vala = [r[pr][vn]["dalpha_rel"] for r in sub]
            print(f"delta {d:.0e} {vn:>6} vs {pr:>4}: dsig_rel max {max(vals):.2e} median {np.median(vals):.2e}; "
                  f"dalpha_rel max {max(vala):.2e} median {np.median(vala):.2e}")
    print("wrote", os.path.join(OUT, "crosscheck.md"))


if __name__ == "__main__":
    main()
