"""WP-138: one summary row per (leg, replay CSV): which binaries / oracle
preset produced it, s/B and step of the committed state, and the error of each
C++ integrator against the oracle.

ERROR NORM (the one every WP-138 figure uses): e = ||sigma_cxx - sigma_oracle|| / p'_in,
||.|| the contravariant tensor norm sqrt(s_ij s_ij) over the 6 Voigt stress
components (shear counted twice), sigma_oracle from WP-134 sanisand_reference
preset `uw_model` (Radau, rtol 1e-10), p'_in the committed p' of the point.
Only rows the oracle integrated (status ok) count.

    python replay_summary.py LEG [LEG ...]   -> markdown table on stdout
"""
import csv
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def nrm(v):
    return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2 + 2 * (v[3]**2 + v[4]**2 + v[5]**2))


def stats(xs):
    xs = sorted(xs)
    if not xs:
        return "-", "-", "-"
    return (f"{xs[len(xs)//2]:.2e}", f"{xs[int(0.9*(len(xs)-1))]:.2e}", f"{xs[-1]:.2e}")


print("| leg | dataset | step (s/B) | n | oracle | ME build | ME median / p90 / max | SAS build | SAS median / p90 / max | SAS refusals |")
print("|---|---|---|---|---|---|---|---|---|---|")
for leg in sys.argv[1:]:
    rd = os.path.join(HERE, "runs", leg)
    for c in sorted(glob.glob(os.path.join(rd, "replay", "*.csv"))):
        base = os.path.splitext(os.path.basename(c))[0]
        od = os.path.join(rd, "replay_out")
        op = os.path.join(od, f"{base}.oracle.json")
        if not os.path.exists(op):
            continue
        rows = {(int(r["element"]), int(r["gp"])): r for r in csv.DictReader(open(c, newline=""))}
        orc = json.load(open(op))
        O = {(d["element"], d["gp"]): d for d in orc["rows"]}
        out = []
        for lab in ("ME", "SAS"):
            p = os.path.join(od, f"{base}.{lab}.json")
            if not os.path.exists(p):
                out += ["-", "-", "-"]
                continue
            J = json.load(open(p))
            errs, ref = [], 0
            for d in J["rows"]:
                k = (d["element"], d["gp"])
                o = O.get(k)
                if d["rc"] != 0:
                    ref += 1
                if o and o.get("status") == "ok" and d["rc"] == 0:
                    pin = float(rows[k]["p_kPa"])
                    errs.append(nrm([x - y for x, y in zip(d["sigma"], o["sigma"])]) / pin)
            out += [J["build"][:9], " / ".join(stats(errs)), str(ref)]
        r0 = next(iter(rows.values()))
        nok = sum(1 for d in orc["rows"] if d.get("status") == "ok")
        print(f"| {leg} | {base} | {r0['step']} ({float(r0['s_over_B']):.5f}) | {len(rows)} | "
              f"{orc.get('preset')} ok {nok} | {out[0]} | {out[1]} | {out[3]} | {out[4]} | {out[5]} |")
