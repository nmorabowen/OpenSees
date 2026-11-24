"""WP-253: the MPM throughput regression check for later MPM work packages.

A small fixed set of five rows, 3 repeats each, against a committed baseline (`baseline_nmb.json`,
recorded by the orchestrator on nmb with `--record`). A row REGRESSES when its median T falls by
more than max(10 %, 2 x the baseline row's spread). Exit code 0 (no regression), 1 (a regression
or a failed row), 2 (cannot run: build mismatch, no baseline). Wall-clock: run it on an otherwise
idle nmb, never in Zone-A and never beside a build.

    python3.12 Ladruno_scripts/mpm_bench/regress.py                      # check against the baseline
    python3.12 Ladruno_scripts/mpm_bench/regress.py --record             # (re)write the baseline
    python3.12 Ladruno_scripts/mpm_bench/regress.py --pyd X --exe Y      # another build
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bench as Bn          # noqa: E402
import benchlib as B        # noqa: E402
import summarize as S       # noqa: E402

C = B.C
BASELINE = Path(__file__).resolve().parent / "baseline_nmb.json"
FLOOR = 0.10            # a row regresses past max(FLOOR, 2 x baseline spread)
SPREAD_FACTOR = 2.0
REPEATS = 3

ROWS = [   # (mode, dim, basis, mat, N, threads)
    ("F", 3, "gimp", "elastic", 100_000, 1),
    ("F", 3, "gimp", "elastic", 100_000, 8),
    ("F", 3, "gimp", "dp", 100_000, 1),
    ("F", 2, "bspline2", "elastic", 100_000, 8),
    ("A", 3, "gimp", "elastic", 10_000, 1),
]


def configs():
    return [dict(zip(("mode", "dim", "basis", "mat", "N", "threads"), r), profile="coarse")
            for r in ROWS]


def threshold(base_spread):
    """The fraction of the baseline T below which a row regresses."""
    return max(FLOOR, SPREAD_FACTOR * (base_spread or 0.0))


def judge(base_rows, agg):
    """[(key_str, base T, now T, change, allowed, verdict)] and whether anything regressed."""
    out, bad = [], False
    for r in ROWS:
        ks = B.key_str(r)
        info = agg.get(r)
        b = base_rows.get(ks)
        now = info["T"]["median"] if info else None
        if now is None:
            out.append((ks, b and b["T_median"], None, None, None, "FAILED"))
            bad = True
            continue
        if not b:
            out.append((ks, None, now, None, None, "no baseline"))
            continue
        change = now / b["T_median"] - 1.0
        allowed = threshold(b.get("spread"))
        v = "REGRESSION" if change < -allowed else "ok"
        bad |= v == "REGRESSION"
        out.append((ks, b["T_median"], now, change, allowed, v))
    return out, bad


def main():
    ap = Bn.build_parser()
    ap.add_argument("--record", action="store_true", help="write the baseline instead of checking")
    ap.add_argument("--baseline", default=str(BASELINE))
    args = C.parse_args(ap)
    if args.child:
        return Bn.child(C.load_cfg(args.child))
    base_path = Path(args.baseline)
    if not args.record and not base_path.exists():
        print(f"no baseline at {base_path}; run with --record on the idle reference machine")
        return 2
    args.repeats = args.repeats or REPEATS
    od = C.out_dir(args, "regress")
    try:
        doc = Bn.new_doc(args, "nmb", od)
    except SystemExit as e:
        print(e)
        return 2
    Bn.run_configs(configs(), args, od, REPEATS, doc, "regress")
    agg = S.aggregate([doc["rows"]])
    C.write_json(od / "regress_results.json", doc)
    if args.record:
        rec = {"recorded": doc["conditions"]["started"], "builds": doc["builds"], "machine": doc["machine"],
               "repeats": args.repeats, "rows": {B.key_str(r): {
                   "T_median": agg[r]["T"]["median"], "spread": agg[r]["T"]["spread"],
                   "N": agg[r].get("N"), "n_ok": agg[r]["n_ok"]} for r in ROWS if r in agg}}
        C.write_json(base_path, rec)
        print(f"baseline written to {base_path}")
        return 0 if len(rec["rows"]) == len(ROWS) else 1
    base = json.loads(base_path.read_text(encoding="utf-8"))
    res, bad = judge(base["rows"], agg)
    print(f"\nbaseline build {base['builds']['exe']}  now {doc['builds']['exe']}\n")
    print(f"{'row':44s} {'baseline T':>12s} {'now T':>12s} {'change':>8s} {'allowed':>8s}  verdict")
    for ks, bt, nt, ch, al, v in res:
        f = lambda x: "-" if x is None else f"{x:.4g}"                  # noqa: E731
        p = lambda x: "-" if x is None else f"{100 * x:+.1f}%"          # noqa: E731
        print(f"{ks:44s} {f(bt):>12s} {f(nt):>12s} {p(ch):>8s} {('-' if al is None else f'-{100 * al:.0f}%'):>8s}  {v}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
