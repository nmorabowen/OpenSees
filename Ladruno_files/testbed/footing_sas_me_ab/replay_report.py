"""WP-138: join a replay CSV with its C++ replay JSON(s) and its oracle JSON
into one per-point table (markdown to stdout, CSV beside the first JSON).

    python replay_report.py --csv R.csv --oracle O.json --cxx A.json [--cxx B.json ...]
        [--label A --label B] [--out table.csv]

Columns per point: p', rho_alpha in; oracle status, rho_alpha at exit, f at exit;
per C++ binary: rc (0 = integrated, else refused), substeps, cap hit,
rho_alpha out, f after, |sigma_cxx - sigma_oracle| / p'.
"""
import argparse
import csv
import json
import math


def nrm(v):
    return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2 + 2 * (v[3]**2 + v[4]**2 + v[5]**2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--oracle")
    ap.add_argument("--cxx", action="append", default=[])
    ap.add_argument("--label", action="append", default=[])
    ap.add_argument("--out")
    a = ap.parse_args()
    labels = a.label or [f"C{i}" for i in range(len(a.cxx))]
    rows = list(csv.DictReader(open(a.csv, newline="")))
    key = lambda d: (int(d["element"]), int(d["gp"]))
    orc = {key(d): d for d in json.load(open(a.oracle))["rows"]} if a.oracle else {}
    cx = [{key(d): d for d in json.load(open(p))["rows"]} for p in a.cxx]
    head = ["ele", "gp", "sel", "p_in", "rho_in", "or_status", "or_rho_end", "or_f_end"]
    for L in labels:
        head += [f"{L}_rc", f"{L}_sub", f"{L}_cap", f"{L}_rho_out", f"{L}_f_after",
                 f"{L}_dsig_or_rel"]
    table = []
    for r in rows:
        k = key(r)
        o = orc.get(k, {})
        line = [k[0], k[1], r["select"], float(r["p_kPa"]), float(r["rho_alpha"]),
                o.get("status", "-"), o.get("rho_alpha_end", float("nan")),
                o.get("f_end", float("nan"))]
        for c in cx:
            d = c.get(k)
            if d is None:
                line += ["-"] * 6
                continue
            ds = float("nan")
            if o.get("sigma") and o.get("status") == "ok":
                ds = nrm([x - y for x, y in zip(d["sigma"], o["sigma"])]) / max(float(r["p_kPa"]), 1e-9)
            if "sas" in d:
                sub, cap = d["sas_substeps"], d["sas_last_refuse_code"]
            else:
                sub, cap = int(d["stats"].get("substeps", 0)), int(d["stats"].get("capHits", 0))
            line += [d["rc"], sub, cap, d["rho_alpha_out"], d["f_after"], ds]
        table.append(line)
    fmt = lambda x: (f"{x:.4g}" if isinstance(x, float) else str(x))
    print("| " + " | ".join(head) + " |")
    print("|" + "---|" * len(head))
    for line in table:
        print("| " + " | ".join(fmt(x) for x in line) + " |")
    if a.out:
        with open(a.out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(head)
            w.writerows(table)
    # summary
    n = len(table)
    ok_or = sum(1 for t in table if t[5] == "ok")
    print(f"\nrows {n}; oracle ok {ok_or}")
    for i, L in enumerate(labels):
        base = 8 + 6 * i
        rc = [t[base] for t in table if t[base] != "-"]
        ds = [t[base + 5] for t in table if t[base] != "-" and isinstance(t[base + 5], float)
              and not math.isnan(t[base + 5])]
        print(f"{L}: refused {sum(1 for x in rc if x != 0)}/{len(rc)}; "
              f"max |dsig - oracle|/p' {max(ds) if ds else float('nan'):.3e}; "
              f"median {sorted(ds)[len(ds)//2] if ds else float('nan'):.3e}")


if __name__ == "__main__":
    main()
