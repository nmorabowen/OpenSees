"""Roll the four ADR-86b T8 legs into one table, with the GATE U baseline."""
import csv
import glob
import json
import os
import sys

root = sys.argv[1]

rows = []
for d in sorted(glob.glob(os.path.join(root, "*"))):
    if not os.path.isdir(d):
        continue
    js = sorted(glob.glob(os.path.join(d, "a2_*.json")))
    if not js:
        print("NO JSON in", d)
        continue
    r = json.load(open(js[0]))
    curve = sorted(glob.glob(os.path.join(d, "a2_*_curve.csv")))
    worst = 0.0
    prev = 0.0
    n_rows = 0
    s_last = q_last = 0.0
    if curve:
        with open(curve[0]) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                break
            rd = csv.reader(fh)
            for rec in rd:
                if not rec or rec[0] == "s_m":
                    continue
                t = float(rec[6])
                worst = max(worst, t - prev)
                prev = t
                n_rows += 1
                s_last, q_last = float(rec[1]), float(rec[2])
    rows.append(dict(
        leg=os.path.basename(d), cap=r.get("max_substeps"), mode=r["mode"],
        verdict=r["verdict"], partial=r.get("partial", False),
        s_end=s_last or r["s_end_over_B"], q_end=q_last,
        q_u=r["q_u"], steps=n_rows or r["steps"], nsub=r["nsub"],
        nfail=r["nfail"], nrelax=r["nrelax"], wall=r["wall_s"],
        worst_step_s=worst,
        cps={c["cp_target"]: (c["s_over_B"], c["q_foot"]) for c in r["checkpoints"]}))

hdr = ("leg", "cap", "mode", "s/B end", "q end", "steps", "nsub", "nfail",
       "nrelax", "wall s", "worst step s")
print("%-24s %6s %-8s %9s %9s %6s %5s %6s %6s %8s %12s" % hdr)
for r in rows:
    print("%-24s %6s %-8s %9.5f %9.2f %6d %5d %6d %6d %8.1f %12.1f"
          % (r["leg"], r["cap"], r["mode"], r["s_end"], r["q_end"], r["steps"],
             r["nsub"], r["nfail"], r["nrelax"], r["wall"], r["worst_step_s"]))
print()
for r in rows:
    print(r["leg"], "->", r["mode"], "|", r["verdict"])
print()
print("checkpoints (s/B target -> actual s/B, q_foot kPa)")
for r in rows:
    print(" ", r["leg"], {k: (round(v[0], 5), round(v[1], 3))
                          for k, v in sorted(r["cps"].items())})
