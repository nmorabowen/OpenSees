"""WP-106 / ADR-93 II.1 -- read `probe_pre.py` dumps and report the II.1 numbers.

    python3.12 analyze_tx.py data/tx_p20_pre0.csv data/tx_p20_pre1.csv

Per file: peak `eta = q/p`, `M^b` at the peak step, `eta/M^b` at the peak, the
substep census, and (against the FIRST file) the largest relative difference in
committed stress over the common steps.

`M^b = Mc * g(theta, c) * exp(-nb * psi)` (Dafalias & Manzari 2004; the code's own
`GetStateDependent`), with `g = 2c / ((1 + c) - (1 - c) cos3theta)` and `psi` read
from the material's own F4 diagnostic response -- so the report cannot drift from
what the binary computed.
"""
from __future__ import annotations

import csv
import math
import sys

MC, C_RATIO, NB = 1.3309, 0.71, 3.5
ONE3 = 1.0 / 3.0


def read(path):
    meta, rows = {}, []
    with open(path, encoding="utf-8") as fh:
        lines = fh.readlines()
    body = []
    for ln in lines:
        if ln.startswith("#"):
            k, _, v = ln[1:].partition(":")
            meta[k.strip()] = v.strip()
        else:
            body.append(ln)
    for r in csv.DictReader(body):
        rows.append({k: float(v) for k, v in r.items() if v not in ("", None)})
    return meta, rows


def invariants(sig):
    """sig = (s11,s22,s33,s12,s23,s13), compression-positive as the material stores it."""
    p = ONE3 * (sig[0] + sig[1] + sig[2])
    s = [sig[0] - p, sig[1] - p, sig[2] - p, sig[3], sig[4], sig[5]]
    # ||s|| with the contravariant doubling on the shear slots
    n2 = s[0] ** 2 + s[1] ** 2 + s[2] ** 2 + 2.0 * (s[3] ** 2 + s[4] ** 2 + s[5] ** 2)
    sn = math.sqrt(max(n2, 0.0))
    q = math.sqrt(1.5) * sn
    # cos(3 theta) from the deviator's third invariant, the material's own form
    if sn < 1e-14:
        cos3t = 1.0
    else:
        nn = [x / sn for x in s]
        # tr(n^3), done explicitly on the symmetric 3x3 the Voigt slots name
        a, b, c_, d, e, f = nn
        M = [[a, d, f], [d, b, e], [f, e, c_]]
        M2 = [[sum(M[i][k] * M[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
        t3 = sum(M2[i][k] * M[k][i] for i in range(3) for k in range(3))
        cos3t = math.sqrt(6.0) * t3
        cos3t = max(-1.0, min(1.0, cos3t))
    return p, q, sn, cos3t


def g_lode(cos3t):
    return 2.0 * C_RATIO / ((1.0 + C_RATIO) - (1.0 - C_RATIO) * cos3t)


def summarise(path, base=None):
    meta, rows = read(path)
    best = None
    subs = 0.0
    caps = 0
    for r in rows:
        sig = [r[f"sig{i}"] for i in range(6)]
        p, q, sn, cos3t = invariants(sig)
        if p <= 1e-12:
            continue
        eta = q / p
        psi = r.get("psi", float("nan"))
        Mb = MC * g_lode(cos3t) * math.exp(-NB * psi) if psi == psi else float("nan")
        if best is None or eta > best["eta"]:
            best = dict(step=int(r["step"]), eta=eta, p=p, q=q, Mb=Mb, psi=psi,
                        ratio=(eta / Mb if Mb == Mb and Mb > 0 else float("nan")))
        s = r.get("substeps_me", float("nan"))
        if s == s:
            subs += s
        ch = r.get("substeps_capHit", 0.0)
        if ch == ch and ch:
            caps += 1
    out = dict(path=path, build=meta.get("build"), pRe=meta.get("pRe"),
               p0=meta.get("p0"), steps=len(rows) - 1, substeps_total=subs,
               cap_hits=caps, **{f"peak_{k}": v for k, v in best.items()})
    if base is not None:
        _, brows = read(base)
        n = min(len(rows), len(brows))
        worst, worst_i = 0.0, -1
        for i in range(n):
            a = [rows[i][f"sig{j}"] for j in range(6)]
            b = [brows[i][f"sig{j}"] for j in range(6)]
            nb_ = math.sqrt(sum(x * x for x in b))
            d = math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))
            rel = d / nb_ if nb_ > 1e-14 else d
            if rel > worst:
                worst, worst_i = rel, i
        out["max_rel_dsigma_vs_base"] = worst
        out["max_rel_at_step"] = worst_i
    return out


def main():
    paths = sys.argv[1:]
    if not paths:
        raise SystemExit(__doc__)
    base = paths[0]
    for i, p in enumerate(paths):
        s = summarise(p, base=None if i == 0 else base)
        print(f"--- {p}")
        for k, v in s.items():
            if k == "path":
                continue
            print(f"    {k:28s} {v}")
    if len(paths) == 2:
        a = summarise(paths[0])
        b = summarise(paths[1], base=paths[0])
        for key in ("peak_eta", "peak_Mb", "peak_ratio"):
            va, vb = a[key], b[key]
            if va == va and vb == vb and va != 0:
                print(f"    delta {key:22s} {100.0 * (vb - va) / va:+.4f} %")


if __name__ == "__main__":
    main()
