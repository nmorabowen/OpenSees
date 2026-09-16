"""WP-106 / ADR-93 II.1 -- arm-against-arm, step for step.

`analyze_tx.py` reports each leg's own peak, which is only comparable when two
legs reach the same strain; the whole point of the floor is that they often do
NOT (the unfloored arm stalls first). This reports what IS comparable:

  * REACH  -- how many steps of the identical prescribed path each arm completed
  * COST   -- ModifiedEuler substeps over the steps BOTH arms completed, and the
              ratio (>1 = the floor is cheaper)
  * CURVE  -- the committed-stress difference over those same steps, median and
              max, with the step the max sits on (a single non-converged step
              can own the max and says nothing about the curve)
  * STRENGTH -- mobilised eta and eta/M^b at the last COMMON step

    python3.12 compare_arms.py <base.csv> <arm.csv> [<arm.csv> ...]
"""
from __future__ import annotations

import csv
import math
import statistics
import sys

MC, C_RATIO, NB = 1.3309, 0.71, 3.5
ONE3 = 1.0 / 3.0


def read(path):
    meta, body = {}, []
    for ln in open(path, encoding="utf-8"):
        if ln.startswith("#"):
            k, _, v = ln[1:].partition(":")
            meta[k.strip()] = v.strip()
        else:
            body.append(ln)
    rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(body)]
    return meta, rows


def inv(r):
    """Material convention (compression positive): eleResponse is the element's."""
    s = [-r[f"sig{i}"] for i in range(6)]
    p = ONE3 * (s[0] + s[1] + s[2])
    d = [s[0] - p, s[1] - p, s[2] - p, s[3], s[4], s[5]]
    n2 = sum(x * x for x in d[:3]) + 2.0 * sum(x * x for x in d[3:])
    sn = math.sqrt(max(n2, 0.0))
    q = math.sqrt(1.5) * sn
    if sn < 1e-14:
        cos3t = 1.0
    else:
        a, b, c_, dd, e, f = [x / sn for x in d]
        M = [[a, dd, f], [dd, b, e], [f, e, c_]]
        M2 = [[sum(M[i][k] * M[k][j] for k in range(3)) for j in range(3)]
              for i in range(3)]
        cos3t = max(-1.0, min(1.0, math.sqrt(6.0)
                              * sum(M2[i][k] * M[k][i] for i in range(3)
                                    for k in range(3))))
    g = 2.0 * C_RATIO / ((1.0 + C_RATIO) - (1.0 - C_RATIO) * cos3t)
    psi = r.get("psi", float("nan"))
    Mb = MC * g * math.exp(-NB * psi) if psi == psi else float("nan")
    return p, q, (q / p if p > 1e-12 else float("nan")), Mb


def nrm(v):
    return math.sqrt(sum(x * x for x in v))


def main():
    paths = sys.argv[1:]
    if len(paths) < 2:
        raise SystemExit(__doc__)
    bmeta, brows = read(paths[0])
    print(f"base: {paths[0]}  pRe={bmeta.get('pRe')}  p0={bmeta.get('p0')}  "
          f"reach={len(brows) - 1} steps  build={bmeta.get('build', '')[:9]}")
    head = ["arm", "pRe", "reach", "common", "substeps base", "substeps arm",
            "cost ratio", "median d(sig)", "max d(sig)", "@step",
            "eta base", "eta arm", "eta/Mb base", "eta/Mb arm"]
    print("| " + " | ".join(head) + " |")
    print("|" + "|".join("---" for _ in head) + "|")
    for path in paths[1:]:
        meta, rows = read(path)
        n = min(len(rows), len(brows))
        sb = sum(brows[i].get("substeps_me", 0.0) for i in range(n))
        sa = sum(rows[i].get("substeps_me", 0.0) for i in range(n))
        ds, worst, wstep = [], 0.0, -1
        for i in range(n):
            u = [brows[i][f"sig{j}"] for j in range(6)]
            v = [rows[i][f"sig{j}"] for j in range(6)]
            nu = nrm(u)
            rel = (nrm([x - y for x, y in zip(u, v)]) / nu) if nu > 1e-14 else 0.0
            ds.append(rel)
            if rel > worst:
                worst, wstep = rel, i
        _, _, eb, mbb = inv(brows[n - 1])
        _, _, ea, mba = inv(rows[n - 1])
        print("| " + " | ".join([
            path.split("/")[-1], f"{float(meta.get('pRe', 0)):g}",
            str(len(rows) - 1), str(n - 1), f"{sb:.0f}", f"{sa:.0f}",
            f"{sb / max(sa, 1):.3f}x",
            f"{statistics.median(ds):.2e}", f"{worst:.2e}", str(wstep),
            f"{eb:.5f}", f"{ea:.5f}",
            f"{eb / mbb:.5f}", f"{ea / mba:.5f}"]) + " |")


if __name__ == "__main__":
    main()
