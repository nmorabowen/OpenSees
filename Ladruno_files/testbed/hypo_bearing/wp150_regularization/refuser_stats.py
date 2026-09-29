"""Refusing GPs (floor_refusers_in_band.csv) -> their committed H decomposition and
cumulative SAS counters at the last converged checkpoint, plus the band-wide
statistics of the h-singularity exposure: for post-peak points (b:n < 0) the
factor by which x = (alpha - alpha_in):n must shrink for H <= 0."""
import csv, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from h_decomp import decompose

AN = sys.argv[1]
SAS = ["updates", "elastic", "substeps", "accepted", "rejectedErr", "rejectedLowP",
       "rejectedNonPosH", "rejectedDrift", "rejectedAlpha", "elasticStages",
       "driftCorrections", "alphaInReseats", "hBrackets", "alphaProjected",
       "intersectFail", "refusals", "refStartF", "refStartAlpha", "refStartOther",
       "refDTmin", "refNonPosH", "refLowP", "refDrift", "refAlpha", "refCap",
       "maxSubstepsOneUpdate", "lastSubsteps", "lastRefuseCode", "maxAlphaRatio",
       "lastAlphaRatio", "lastF", "entryOverKappa", "rejectedReversal"]
rows = list(csv.DictReader(open(os.path.join(AN, "tables", "floor_refusers_in_band.csv"))))
for leg in ("E_B", "E_B16"):
    f = os.path.join(AN, "ck", leg, "ckpt", "field_last_converged.npz")
    R, sb = decompose(f)
    d = np.load(f, allow_pickle=True)
    stats = d["stats"]
    print(f"\n=== {leg} last converged s/B {sb:.4f}")
    for r in rows:
        if r["leg"] != leg:
            continue
        x, y = float(r["x"]), float(r["y"])
        k = int(np.argmin((d["gx"] - x) ** 2 + (d["gy"] - y) ** 2))
        q = R[k]
        print(f" refuser ele {r['element']} gp {r['gp']} ({x:+.3f},{y:+.3f}) {r['codes']}  -> k {k}")
        print(f"   p {q[3]:.1f} psi {q[4]:+.3f} x=(a-ain):n {q[5]:+.3e} h {q[6]:.2e} b:n {q[7]:+.3e} "
              f"Kp/T2 {q[8]/q[9]:+.3f} T3/T2 {q[10]/q[9]:+.3f} H/T2 {q[11]/q[9]:+.3f} D {q[12]:+.3f} z:n {q[13]:+.2f}")
        s = stats[k]
        print("   " + "  ".join(f"{n}={s[i]:.0f}" for i, n in enumerate(SAS)
                              if n in ("updates", "substeps", "rejectedNonPosH", "alphaInReseats",
                                       "hBrackets", "rejectedReversal", "refusals", "refNonPosH",
                                       "rejectedAlpha", "maxAlphaRatio")))
    # band-wide exposure: post-peak points and how close x is to the H=0 critical x
    bn, x, Kp, T2, T3 = R[:, 7], R[:, 5], R[:, 8], R[:, 9], R[:, 10]
    post = (bn < 0) & (x > 0)
    # H = Kp(x) + T2 + T3, Kp = c/x  -> x_crit = -c/(T2+T3), c = Kp*x
    c = Kp * x
    xcrit = np.where(post, -c / (T2 + T3), np.nan)
    ratio = xcrit / x
    print(f" post-peak GPs (b:n<0): {int(post.sum())}; x_crit/x  max {np.nanmax(ratio) if post.any() else float('nan'):.3f}"
          f"  (>1 would already be H<=0); count x_crit/x > 0.1: {int(np.nansum(ratio > 0.1))}")
    tot = stats.sum(axis=0)
    print("  domain totals: " + "  ".join(f"{n}={tot[i]:.0f}" for i, n in enumerate(SAS)
                                          if n in ("alphaInReseats", "hBrackets", "rejectedNonPosH",
                                                   "rejectedReversal", "refNonPosH", "refusals")))
