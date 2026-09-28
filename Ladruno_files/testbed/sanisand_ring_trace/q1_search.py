"""WP-128 Q1 (finding B): from benign K0 states at p0 = 2 / 5 / 20 kPa, drive
a single point along plane-strain strain paths toward eta/M^b -> 1 with
reversals (so alpha_in flips), at increment sizes 1e-6 / 1e-5 / 1e-4, and
record whether a COMMITTED alpha leaves the bounding surface
(sqrt(3/2)||alpha|| / alpha^b_theta > 1, the model's own alpha^b at the
point's Lode angle).  For each path: the first increment where that ratio
exceeds 1.05, and the census of that increment.

C++ chain (IntScheme 1, campaign flags).  Output: out/q1_search.json +
out/q1_search.txt.
"""
import json
import math
import sys

import _boot as B
import drive as D

DIRS = {
    # compression-positive, engineering shear, plane strain
    "passive": [1.0, -1.0, 0.0, 0.0, 0.0, 0.0],     # xx compress, yy extend (isochoric)
    "active": [-1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    "shear": [0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    "extShear": [0.3, -1.0, 0.0, 0.5, 0.0, 0.0],   # dilating: p falls
    "vertUnload": [0.0, -1.0, 0.0, 0.0, 0.0, 0.0],  # footing-edge heave: p falls, eta rises
}


def cycles(d, delta, n_fwd, n_back, ncyc):
    incs = []
    for c in range(ncyc):
        incs += [[delta * x for x in d]] * n_fwd
        incs += [[-delta * x for x in d]] * n_back
    return incs


def main(quick=False):
    B.define_prototypes()
    res = []
    lines = []
    for p0 in (2.0, 5.0, 20.0):
        for dname, d in DIRS.items():
            for delta in (1e-6, 1e-5, 1e-4):
                # total strain per leg ~ 2e-3 (enough to reach M^b from K0 at these p)
                n = int(round(2e-3 / delta))
                n = min(n, 2000)
                incs = cycles(d, delta, n, n // 2, 3)
                st = D.k0_state(p0)
                h = D.run(st, incs, "cpp", tag=B.TAG_ME)
                ratios = [x["alpha_over_b"] for x in h if math.isfinite(x["alpha_over_b"])]
                first = next((x for x in h if x["alpha_over_b"] > 1.05), None)
                worst = max(h, key=lambda x: x["alpha_over_b"] if math.isfinite(x["alpha_over_b"]) else -1)
                tot = {k: sum(x[k] for x in h) for k in ("substeps", "forced", "clamp",
                                                          "abandon", "cap", "pnReset",
                                                          "entryPmin")}
                nfail = sum(1 for x in h if x["rc"] != 0)
                fpos = sum(1 for x in h if x["rc"] == 0 and x["f"] > 1e-7 * max(1, x["p"]) and x["f"] > 1e-7)
                rec = dict(p0=p0, dir=dname, delta=delta, n=len(h),
                           max_alpha_over_b=max(ratios) if ratios else float("nan"),
                           worst_k=worst["k"], worst_eta=worst["eta"], worst_p=worst["p"],
                           worst_eta_alpha=worst["eta_alpha"],
                           first_k=first["k"] if first else None,
                           first_info={k: first[k] for k in ("path", "substeps", "forced",
                                                             "clamp", "abandon", "pnReset",
                                                             "entryPmin", "p", "eta",
                                                             "eta_alpha", "alpha_over_b", "f")}
                           if first else None,
                           totals=tot, refused=nfail, committed_f_pos=fpos,
                           min_p=min(x["p"] for x in h))
                res.append(rec)
                s = (f"p0={p0:>4} {dname:<10} d={delta:.0e} n={len(h):>5} "
                     f"max a/ab={rec['max_alpha_over_b']:.3f} (k={worst['k']}, p={worst['p']:.3g}, "
                     f"eta={worst['eta']:.3g}) first>1.05 k={rec['first_k']} "
                     f"minp={rec['min_p']:.3g} sub={tot['substeps']} forced={tot['forced']} "
                     f"clamp={tot['clamp']} aband={tot['abandon']} pnR={tot['pnReset']} "
                     f"cap={tot['cap']} refused={nfail} f>0 commits={fpos}")
                print(s, flush=True)
                lines.append(s)
    with open(f"{B.OUT}/q1_search.json", "w") as fh:
        json.dump(res, fh, indent=1)
    with open(f"{B.OUT}/q1_search.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
