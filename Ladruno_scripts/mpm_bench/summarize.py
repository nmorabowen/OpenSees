"""WP-253: aggregate bench.py result files into the section-1 tables.

    python3.12 summarize.py results/bench_nmb_*/bench_results.json [more.json ...] \
        [--baseline base/bench_results.json] [--md out.md]

* a row's value is the MEDIAN over its coarse repeats; spread = (max - min) / median;
* a row with spread > 10 % and fewer than 5 repeats is FLAGGED for the 5-repeat re-run
  (`bench.py --rerun-from <file>`); when several files carry the same row, the one with the most
  repeats is used (the 5-repeat median replaces the 3-repeat one; both stay in the input files);
* two medians differ only if their ratio differs from 1 by more than the sum of their spreads
  (section 4); otherwise the difference is reported as "within noise";
* S(n) = T(n) / T(1) in the same file; G(n) = T_after(n) / T_baseline(1) with `--baseline`.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchlib as B   # noqa: E402

C = B.C
SPREAD_RERUN = 0.10
RERUN_REPEATS = 5
SCOPE_ORDER = ("mpm.f.dtcrit", "mpm.f.admissible", "mpm.f.support", "mpm.f.p2g", "mpm.f.grid",
               "mpm.f.g2p", "mpm.f.project2", "mpm.f.update", "mpm.f.commit", "commit",
               "mpm.advance", "other")


def rkey(r):
    return (r["mode"], r["dim"], r["basis"], r["mat"], r["N_nominal"], r["threads"])


def within_noise(ma, sa, mb, sb):
    """Section 4: the medians differ only if |ma / mb - 1| exceeds the sum of the spreads."""
    if not ma or not mb:
        return True
    return abs(ma / mb - 1.0) <= (sa or 0.0) + (sb or 0.0)


def _med(vals):
    return C.median_spread(vals)["median"]


def aggregate(row_lists):
    """{key: info} from one row list per input file (in order)."""
    groups, deeps = {}, {}
    for rows in row_lists:
        per, dper = {}, {}
        for r in rows:
            (per if r["profile"] == "coarse" else dper).setdefault(rkey(r), []).append(r)
        for k, rs in per.items():
            groups.setdefault(k, []).append(rs)
        for k, rs in dper.items():
            deeps.setdefault(k, []).append(rs)
    out = {}
    for k, gs in groups.items():
        best = None
        for g in gs:                                  # most OK repeats wins; ties go to the later file
            n_ok = sum(1 for r in g if r["ok"])
            if best is None or n_ok >= best[0]:
                best = (n_ok, g)
        n_ok, g = best
        ok = [r for r in g if r["ok"]]
        ts = C.median_spread([r["T"] for r in ok])
        info = {"key": k, "T": ts, "n_ok": n_ok, "n_runs": len(g),
                "failures": [{"rep": r["repeat"], "reason": r["reason"]} for r in g if not r["ok"]],
                "flag_rerun": bool(ts["spread"] is not None and ts["spread"] > SPREAD_RERUN
                                   and n_ok < RERUN_REPEATS),
                "groups_seen": [(len(x), sum(1 for r in x if r["ok"])) for x in gs]}
        if ok:
            info["N"] = ok[0]["N"]
            info["ns_per_particle_step"] = _med([r["ns_per_particle_step"] for r in ok])
            info["builds"] = sorted({(r["ladrunoBuild_exe"], r["ladrunoBuild_pyd"]) for r in ok})
            info["machine"] = ok[0]["machine"]
            info["mkl"] = ok[0].get("mkl")
            info["n_timed"] = ok[0]["n_timed"]
            ph = {}
            for r in ok:
                for s, v in ((r.get("phases") or {}).get("share") or {}).items():
                    ph.setdefault(s, []).append(v)
            info["phase_share"] = {s: _med(v) for s, v in ph.items()}
            info["scopes_present"] = bool((ok[0].get("phases") or {}).get("scopes_present"))
            if ok[0]["mode"] == "A":
                info["newton_k_mean"] = _med([r["newton_k_mean"] for r in ok])
                info["T_per_newton_iter"] = _med([r["T_per_newton_iter"] for r in ok
                                                  if r.get("T_per_newton_iter")])
            digs = {json.dumps(r["digest"], sort_keys=True) for r in ok if r.get("digest")}
            info["digest"] = sorted(digs)
        dk = [r for g2 in deeps.get(k, []) for r in g2 if r["ok"] and r.get("material_share") is not None]
        info["material_share"] = dk[-1]["material_share"] if dk else None
        out[k] = info
    return out


def flagged_keys(agg):
    return sorted(k for k, v in agg.items() if v["flag_rerun"])


def speedups(agg):
    """{key: S} with S = T(n) / T(1) of the same (mode, dim, basis, mat, N), plus the noise verdict."""
    out = {}
    for k, v in agg.items():
        k1 = k[:5] + (1,)
        if k1 in agg and k[5] != 1 and agg[k1]["T"]["median"] and v["T"]["median"]:
            b = agg[k1]["T"]
            out[k] = {"S": v["T"]["median"] / b["median"],
                      "within_noise": within_noise(v["T"]["median"], v["T"]["spread"],
                                                   b["median"], b["spread"])}
    return out


def compare(base, after):
    """Per key present in both: T ratio after / baseline at the same thread count, with the noise
    verdict, and the gain G(n) = T_after(n) / T_baseline(1)."""
    out = {}
    for k, a in after.items():
        b = base.get(k)
        if not b or not a["T"]["median"] or not b["T"]["median"]:
            continue
        r = a["T"]["median"] / b["T"]["median"]
        b1 = base.get(k[:5] + (1,))
        out[k] = {"ratio": r, "within_noise": within_noise(a["T"]["median"], a["T"]["spread"],
                                                           b["T"]["median"], b["T"]["spread"]),
                  "G": (a["T"]["median"] / b1["T"]["median"]) if b1 and b1["T"]["median"] else None}
    return out


# --------------------------------------------------------------------------- #
#  markdown                                                                   #
# --------------------------------------------------------------------------- #
def _cell(v):
    if v is None or v["T"]["median"] is None:
        return "fail" if v else "-"
    s = v["T"]["spread"]
    t = v["T"]["median"] / 1e6
    return f"{t:.3g} ({100 * s:.0f}%){'*' if v['flag_rerun'] else ''}"


def tables(agg, base=None):
    lines = []
    keys = sorted(agg, key=lambda k: (k[0], -k[1], k[4], k[2], k[3], k[5]))
    sp = speedups(agg)
    cmpd = compare(base, agg) if base else {}
    for mode in ("F", "A"):
        for dim in (3, 2):
            ks = [k for k in keys if k[0] == mode and k[1] == dim]
            if not ks:
                continue
            threads = sorted({k[5] for k in ks})
            lines += ["", f"### Mode {mode}, {dim}D: T [Mparticle-steps/s] (spread), `*` = spread > 10 %, "
                          "rerun owed", "",
                      "| basis | material | N (actual) | " + " | ".join(f"t={t}" for t in threads)
                      + f" | S({threads[-1]}) | ns/p-step (t=1) |"
                      + (" k (t=1) |" if mode == "A" else ""),
                      "|---|---|---|" + "---|" * (len(threads) + 2 + (mode == "A"))]
            rows = sorted({k[:5] for k in ks}, key=lambda r: (r[4], r[2], r[3]))
            for r in rows:
                cells = [_cell(agg.get(r + (t,))) for t in threads]
                a1 = agg.get(r + (1,)) or next((agg[r + (t,)] for t in threads if r + (t,) in agg), None)
                s = sp.get(r + (threads[-1],))
                sc = "-" if not s else f"{s['S']:.2f}{' (within noise)' if s['within_noise'] else ''}"
                nsp = f"{a1['ns_per_particle_step']:.0f}" if a1 and a1.get("ns_per_particle_step") else "-"
                extra = f" {a1.get('newton_k_mean', 0):.1f} |" if mode == "A" and a1 else (" - |" if mode == "A" else "")
                lines.append(f"| {r[2]} | {r[3]} | {r[4]:g} ({a1['N'] if a1 else '-'}) | "
                             + " | ".join(cells) + f" | {sc} | {nsp} |{extra}")
    # phase split, mode F, threads = 1
    ph = [k for k in keys if k[0] == "F" and k[5] == 1 and agg[k].get("phase_share")]
    if ph:
        lines += ["", "### Mode F phase split at 1 thread: share of the `step` wall [%]", ""]
        names = [s for s in SCOPE_ORDER if any(s in agg[k]["phase_share"] for k in ph)]
        lines += ["| row | " + " | ".join(names) + " | material share |", "|---|" + "---|" * (len(names) + 1)]
        for k in ph:
            v = agg[k]
            sh = v["phase_share"]
            ms = v.get("material_share")
            lines.append(f"| {k[1]}D {k[2]} {k[3]} N={k[4]:g} | "
                         + " | ".join(f"{100 * sh[s]:.1f}" if s in sh else "-" for s in names)
                         + f" | {'-' if ms is None else f'{100 * ms:.1f}'} |")
    if cmpd:
        lines += ["", "### After against baseline (same thread count), and G(n) = T_after(n) / T_baseline(1)", "",
                  "| row | threads | T ratio | verdict | G(n) |", "|---|---|---|---|---|"]
        for k, v in sorted(cmpd.items(), key=lambda kv: (kv[0][0], -kv[0][1], kv[0][4], kv[0][2], kv[0][3], kv[0][5])):
            verdict = "within noise" if v["within_noise"] else ("faster" if v["ratio"] > 1 else "slower")
            g = "-" if v["G"] is None else f"{v['G']:.2f}"
            lines.append(f"| {k[0]} {k[1]}D {k[2]} {k[3]} N={k[4]:g} | {k[5]} | {v['ratio']:.2f} | {verdict} | {g} |")
    det = [(k, v["digest"]) for k, v in agg.items() if v.get("digest")
           and any(k2[:5] == k[:5] and k2[5] != 1 for k2 in agg)]
    if det:
        by = {}
        for k, d in det:
            by.setdefault(k[:5], {})[k[5]] = d
        lines += ["", "### Determinism (particle-array digests against 1 thread)", "",
                  "| row | threads | identical to t=1 |", "|---|---|---|"]
        for r, m in sorted(by.items()):
            for t, d in sorted(m.items()):
                same = "-" if t == 1 or 1 not in m else ("yes" if d == m[1] else "NO")
                lines.append(f"| {r[0]} {r[1]}D {r[2]} {r[3]} N={r[4]:g} | {t} | {same} |")
    meta = {v["machine"]["cpu_model"] + f" | {v['machine']['logical_cores']} logical | " + str(v["machine"]["os"])
            + f" | MKL {v.get('mkl')} | builds {v['builds']}" for v in agg.values() if v.get("machine")}
    lines += ["", "Provenance: " + "; ".join(sorted(meta))] if meta else []
    fl = flagged_keys(agg)
    if fl:
        lines += ["", f"{len(fl)} row(s) with spread > 10 % owed a 5-repeat re-run: "
                  + ", ".join(B.key_str(k) for k in fl)]
    fails = [(k, v["failures"]) for k, v in agg.items() if v["failures"]]
    if fails:
        lines += ["", "Failures (data, not dropped): " + "; ".join(
            f"{B.key_str(k)}: {f['reason']}" for k, fs in fails for f in fs)]
    return "\n".join(lines) + "\n"


def load(paths):
    return [json.loads(Path(p).read_text(encoding="utf-8"))["rows"] for p in paths]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("files", nargs="+")
    ap.add_argument("--baseline", nargs="*", default=None, help="result file(s) of the baseline build")
    ap.add_argument("--md", default=None)
    args = ap.parse_args()
    agg = aggregate(load(args.files))
    base = aggregate(load(args.baseline)) if args.baseline else None
    md = tables(agg, base)
    if args.md:
        Path(args.md).write_text(md, encoding="utf-8")
    sys.stdout.buffer.write(md.encode("utf-8"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
