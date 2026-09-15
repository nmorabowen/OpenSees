"""ADR-92 F10 — reduce the leg JSONs and the census CSVs to the note's tables.

No engine import: this runs on any box against the committed artefacts.

    python3.12 f10_summary.py out
"""
from __future__ import annotations

import csv
import json
import os
import sys

import numpy as np

ORDER = ["A", "B", "C", "D", "E", "F1", "F2", "F3", "G", "H", "I", "J", "K", "L", "M"]
CENSUS = [f"{L}_ds{ds}" for L in ("A", "B", "C", "H")
          for ds in ("2e-05", "4e-05", "8e-05", "0.0002")]


def leg_table(d):
    print("| leg | s/B reached | mode | steps | subdiv | wall s | refusals "
          "total / ctl / comp / sign | floor fallbacks | f=0 guard | trial f=0 |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for leg in ORDER:
        p = os.path.join(d, f"f10_{leg}.json")
        if not os.path.exists(p):
            continue
        m = json.load(open(p))
        r, g = m["refusals"], m["guards"]
        print(f"| {leg} | {m['s_over_B']:.4f} | {m['mode']} | {m['steps']} | "
              f"{m['nsub']} | {m['wall']:.0f} | {r['total']} / {r['control']} / "
              f"{r['companion']} / {r['sign']} | {g[0]} | {g[1]} | {g[4]} |")


def census_table(d):
    print("\n| census | GPs | p' min/med/max kPa | dilatant at rest | "
          "err>0.05 | err>0.1 | err med | err max | >0.05 med p' | "
          ">0.05 med eta/M^d | >0.05 med |z| | >2B from edge |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for leg in CENSUS:
        p = os.path.join(d, f"f10_census_{leg}.csv")
        j = os.path.join(d, f"f10_census_{leg}.json")
        if not (os.path.exists(p) and os.path.exists(j)):
            continue
        meta = json.load(open(j))
        rows = list(csv.DictReader(open(p)))
        e = np.array([float(r["err"]) for r in rows])
        pp = np.array([float(r["p"]) for r in rows])
        em = np.array([float(r["eta_over_Md"]) for r in rows])
        zz = np.array([float(r["z"]) for r in rows])
        dd = np.array([float(r["dx"]) for r in rows])
        o = e > 0.05
        print(f"| {leg} | {len(e)} | {pp.min():.2f} / {np.median(pp):.2f} / "
              f"{pp.max():.1f} | {100*meta['frac_dilatant_at_rest']:.1f} % | "
              f"{int(o.sum())} ({100*o.mean():.1f} %) | "
              f"{int((e>0.1).sum())} ({100*(e>0.1).mean():.1f} %) | "
              f"{np.median(e):.4f} | {e.max():.3f} | "
              f"{np.median(pp[o]):.2f} | {np.median(em[o]):.3f} | "
              f"{-np.median(zz[o]):.2f} | {100*(dd[o] > 3.0).mean():.1f} % |"
              if o.any() else
              f"| {leg} | {len(e)} | {pp.min():.2f} / {np.median(pp):.2f} / "
              f"{pp.max():.1f} | {100*meta['frac_dilatant_at_rest']:.1f} % | 0 | "
              f"0 | {np.median(e):.4f} | {e.max():.3f} | -- | -- | -- | -- |")


def err_vs_state(d, leg):
    """Is the error a function of p' (candidate 2) or of eta/M^d (candidate 1)?"""
    p = os.path.join(d, f"f10_census_{leg}.csv")
    if not os.path.exists(p):
        return
    rows = list(csv.DictReader(open(p)))
    e = np.array([float(r["err"]) for r in rows])
    pp = np.array([float(r["p"]) for r in rows])
    em = np.array([float(r["eta_over_Md"]) for r in rows])
    print(f"\ncensus {leg}: error binned by p' and by eta/M^d")
    print("| p' bin kPa | n | median err | frac > 0.05 |")
    print("|---|---|---|---|")
    qs = np.quantile(pp, np.linspace(0, 1, 6))
    for a, b in zip(qs[:-1], qs[1:]):
        m = (pp >= a) & (pp <= b)
        if m.sum():
            print(f"| {a:.2f}-{b:.2f} | {int(m.sum())} | {np.median(e[m]):.4f} "
                  f"| {100*(e[m] > 0.05).mean():.1f} % |")
    print("| eta/M^d bin | n | median err | frac > 0.05 |")
    print("|---|---|---|---|")
    qs = np.quantile(em, np.linspace(0, 1, 6))
    for a, b in zip(qs[:-1], qs[1:]):
        m = (em >= a) & (em <= b)
        if m.sum():
            print(f"| {a:.3f}-{b:.3f} | {int(m.sum())} | {np.median(e[m]):.4f} "
                  f"| {100*(e[m] > 0.05).mean():.1f} % |")


def overlay(d, ref="G", arms=("B", "D", "F1", "F2", "F3", "I", "J", "K", "M")):
    """ADR-92 section 8's reporting condition: an -implex curve is only
    comparable to the implicit twin OVER THE OVERLAP.  `ref` is the implicit leg.
    """
    def load(leg):
        p = os.path.join(d, f"f10_{leg}.csv")
        if not os.path.exists(p):
            return None
        rows = list(csv.DictReader(open(p)))
        if not rows:
            return None
        x = np.array([float(r["s_over_B"]) for r in rows])
        y = np.array([float(r["q_kPa"]) for r in rows])
        # the adaptive controller can log two rows at the same settlement after a
        # refused-and-retried step; np.interp needs a monotone abscissa.
        o = np.argsort(x, kind="stable")
        return x[o], y[o]
    g = load(ref)
    if g is None:
        return
    print(f"\noverlay against the IMPLICIT leg {ref} "
          f"(overlap 0 < s/B <= {g[0].max():.5f}):")
    print("| arm | overlap s/B | n pts | mean |dev| % | max |dev| % |")
    print("|---|---|---|---|---|")
    for leg in arms:
        a = load(leg)
        if a is None:
            continue
        hi = min(g[0].max(), a[0].max())
        m = (g[0] > 0.2 * hi) & (g[0] <= hi)     # skip the first fifth: both
        if m.sum() < 3:                          # curves are still elastic there
            continue
        qi = np.interp(g[0][m], a[0], a[1])
        dev = 100.0 * (qi - g[1][m]) / g[1][m]
        print(f"| {leg} | {0.2*hi:.5f}-{hi:.5f} | {int(m.sum())} | "
              f"{np.abs(dev).mean():.2f} | {np.abs(dev).max():.2f} |")


def probe_table(d):
    import glob
    rows = []
    for p in glob.glob(os.path.join(d, "f10_probe_*.json")):
        rows.append(json.load(open(p)))
    if not rows:
        return
    rows.sort(key=lambda r: -r["ds"])
    print("")
    print(f"step-size refinement probe at s/B = {rows[0]['s_walked']:.5f} "
          f"(reached on a refusal-free constant-ds path):")
    print("| ds (m) | median err | max err | GPs > 0.05 | max err / previous |")
    print("|---|---|---|---|---|")
    prev = None
    for r in rows:
        ratio = "--" if prev is None else f"{prev/r['max']:.2f}"
        print(f"| {r['ds']:g} | {r['median']:.3e} | {r['max']:.3e} | "
              f"{r['n_over_005']} | {ratio} |")
        prev = r["max"]


if __name__ == "__main__":
    d = sys.argv[1] if len(sys.argv) > 1 else "out"
    leg_table(d)
    census_table(d)
    overlay(d)
    probe_table(d)
    for leg in ("B_ds4e-05", "A_ds4e-05", "C_ds4e-05", "H_ds4e-05"):
        err_vs_state(d, leg)


