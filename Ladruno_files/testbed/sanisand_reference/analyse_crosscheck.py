"""Classify the cross-check cases by the REFERENCE's path structure and report
agreement per class (reads out/crosscheck.json, writes out/crosscheck_classes.md).

Classes (from the uw_me reference, i.e. the ModifiedEuler-faithful continuum):
  A  monotonic plastic   one plastic segment from t = 0 (kink restarts allowed), status ok
  B  through the cone    an elastic segment first (unloading through the m-cone,
                         then plastic on the other side), no h < 0, status ok
  C  UW-rule h < 0       B, but the UW alpha_in rule leaves (alpha - alpha_in):n < 0
                         at the plastic onset (mechanism G) -- the C++'s own rule
  D  no admissible rate  the reference stops (H_nonpositive etc.) under the UW rule
"""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")


def klass(ref):
    r = ref["uw_me"]
    if r["status"] != "ok":
        return "D"
    modes = [m for m, _ in r["segments"]]
    if r["negh"]:
        return "C"
    if modes[0] == "elastic":
        return "B"
    return "A"


def main():
    d = json.load(open(os.path.join(OUT, "crosscheck.json")))
    rows = d["rows"]
    L = ["| delta | class | n | uw_me vs ME8 dsig (median / max) | uw_me vs ME8 dalpha (median / max) | uw vs RK45 dsig (median / max) | paper vs ME8 dsig (median) | ME8 substeps (median) |",
         "|---|---|---|---|---|---|---|---|"]
    for delta in sorted({r["delta"] for r in rows}):
        for k in "ABCD":
            sub = [r for r in rows if r["delta"] == delta and klass(r["ref"]) == k]
            if not sub:
                continue
            a = [r["ME8"]["uw_me"]["dsig_rel"] for r in sub]
            b = [r["ME8"]["uw_me"]["dalpha_rel"] for r in sub]
            b = [x for x in b if np.isfinite(x)] or [float("nan")]
            c = [r["RK45"]["uw"]["dsig_rel"] for r in sub]
            p = [r["ME8"]["paper"]["dsig_rel"] for r in sub]
            s = [r["ME8"]["substeps"] for r in sub]
            L.append(f"| {delta:.0e} | {k} | {len(sub)} | {np.median(a):.1e} / {max(a):.1e} | "
                     f"{np.median(b):.1e} / {max(b):.1e} | {np.median(c):.1e} / {max(c):.1e} | "
                     f"{np.median(p):.1e} | {np.median(s):.0f} |")
    # toggle attribution on class A: effect of each UW addition alone
    L += ["", "Per-toggle effect on class A (median over cases of ||Δσ(uw_me) − Δσ(uw_me with that addition at the paper value)|| / ||Δσ||):", "",
          "| toggle | δ=1e-5 | δ=1e-4 | δ=1e-3 |", "|---|---|---|---|"]
    names = sorted({k for r in rows for k in r.get("toggle_effect", {})})
    for nm in names:
        vals = []
        for delta in sorted({r["delta"] for r in rows}):
            sub = [r["toggle_effect"][nm] for r in rows
                   if r["delta"] == delta and klass(r["ref"]) == "A"]
            vals.append(f"{np.median(sub):.1e}" if sub else "-")
        L.append(f"| {nm} | " + " | ".join(vals) + " |")
    txt = "\n".join(L) + "\n"
    open(os.path.join(OUT, "crosscheck_classes.md"), "w", encoding="utf-8").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
