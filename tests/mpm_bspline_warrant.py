"""WP-232i warrant run: every MMS row of tests/mpm_bspline_cases.txt on its full ladder.

The orchestrator launches it on the 232i build (dist/bin/opensees.pyd on PYTHONPATH, the
tests/_testbed bootstrap) and commits its JSON:

    python3.12 tests/mpm_bspline_warrant.py                      # every row -> warrant.json
    python3.12 tests/mpm_bspline_warrant.py --rows G2 G3         # a subset (merged into the file)
    python3.12 tests/mpm_bspline_warrant.py --resume             # only rows not yet in the file
    python3.12 tests/mpm_bspline_warrant.py --mutant ZERO        # G6 on a ZERO mutant build
    python3.12 tests/mpm_bspline_warrant.py --mutant SCALE       # G6 on a SCALE mutant build

Output: tests/data/mpm/bspline/warrant.json (or g6_<mode>.json), one entry per row with
ladrunoBuild, dt, the validity check, e / e(T/4) / max_t e(t) / crossings / N / wall per rung,
p and the verdict against the row's pre-registered gate. A miss is recorded as a miss; no band
is widened here or anywhere.

The 1/h = 128 rung (the brief: "add 128 if a run costs under 10 min"): for each gated row the
1/h = 64 rung is timed at its own dt first; 128 is added when 8x that (4x particles, 2x steps;
quasi-static: 4x) is under 600 s, and the whole ladder then runs at the 128 rung's dt (the dt rule
is fixed by the finest rung). R1 rows keep their own ladders.

Estimated runtime on the nmb box (from the g++ core and per-step Python costs): about 45 min for
every row, plus about 5 min per mutant.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_bspline_harness as B  # noqa: E402

OUT = Path(__file__).resolve().parent / "data" / "mpm" / "bspline"
LADDER_ROWS = ["G2", "G3", "G4", "G5-G2", "G5-G3", "G5-G4", "R1-gimp-G2", "R1-gimp-G3",
               "R1-linear-G2", "R1-linear-G3", "R1-pic-G2", "R1-pic-G3", "R1-3x3-G2", "R1-3d",
               "R1-modeA-transient"]


def _ops():
    from _testbed import ops
    return ops


def _build(ops):
    try:
        return str(ops.ladrunoBuild())
    except Exception as exc:  # pragma: no cover
        return f"unknown ({exc})"


def _mutation(ops):
    try:
        return str(ops.ladrunoMutation())
    except Exception:
        return "none/unknown"


def ladder_row(ops, row, allow128=True, log=B.flushprint):
    mode = row.get("mode", "explicit")
    rungs = B.as_list(row["ladder"])
    decision = None
    if allow128 and row.get("status") == "gate" and rungs[-1] == 64:
        dt64 = float(row["dt_rule"]) * B.critical_dt(row, 64) if mode == "explicit" else None
        probe = B.mms(ops, row, 64, dt64, record=False)
        factor = 8.0 if mode == "explicit" else 4.0
        est = factor * probe["wall"]
        decision = {"probe_64_wall_s": probe["wall"], "estimated_128_wall_s": est,
                    "add_128": est < 600.0}
        log(f"{row['name']}: 1/h=64 probe {probe['wall']:.1f}s -> 128 estimated {est:.0f}s "
            f"({'added' if est < 600 else 'not added'})")
        if est < 600.0:
            rungs = rungs + [128]
    res = B.ladder(ops, row, rungs=rungs, log=log)
    ok, why = B.verdict(row, res) if row.get("status") == "gate" else (None, "reported")
    res.update({"status": row.get("status"), "gate": row.get("gate"), "verdict": ok,
                "verdict_text": why, "rung_128": decision, "mode": mode, "basis": row["basis"],
                "fmax": row["fmax"], "ppc": row["ppc"], "flip": row.get("flip")})
    if row.get("issue_below") is not None and res["p"] is not None:
        res["file_issue"] = res["p"] < float(row["issue_below"])
    log(f"{row['name']}: e={['%.3e' % e for e in res['e']]} p={res['p']} -> {why}")
    return res


def column_row(ops, row, log=B.flushprint):
    out = {}
    for case in B.as_list(row["cases"]):
        for basis in B.as_list(row["bases"]):
            for k in B.as_list(row["refine"]):
                key = f"{case}|{basis}|r{k}"
                try:
                    out[key] = B.column(ops, case, basis, int(k))
                except Exception as exc:
                    out[key] = {"error": str(exc)}
                log(f"R1-column {key}: {out[key]}")
    return {"status": "report", "gate": row["gate"], "runs": out}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", nargs="*")
    ap.add_argument("--mutant", choices=["ZERO", "SCALE"])
    ap.add_argument("--no-128", action="store_true")
    ap.add_argument("--resume", action="store_true",
                    help="skip the rows warrant.json already holds")
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args(argv)
    ops = _ops()
    rows = B.parse_rows()
    out_dir = Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if a.mutant:
        names = a.rows or ["G2", "G3"]
        path = out_dir / f"g6_{a.mutant.lower()}.json"
        data = {"ladrunoBuild": _build(ops), "ladrunoMutation": _mutation(ops),
                "mutant": a.mutant, "rows": {}}
        for n in names:
            row = dict(rows[n])
            r = B.ladder(ops, row, validity=False)
            ok, why = B.verdict(row, r)
            r.update({"verdict": ok, "verdict_text": why})
            data["rows"][n] = r
        data["red"] = [n for n, r in data["rows"].items() if not r["verdict"]]
        data["G6_killed"] = bool(data["red"])
    else:
        path = out_dir / "warrant.json"
        data = json.loads(path.read_text()) if path.exists() else {"rows": {}}
        data["ladrunoBuild"] = _build(ops)
        data["ladrunoMutation"] = _mutation(ops)
        names = a.rows or (LADDER_ROWS + ["R1-column"])
        if a.resume:
            names = [n for n in names if n not in data["rows"]]
            B.flushprint("resume: running", names)
        for n in names:
            row = rows[n]
            if n == "R1-column":
                data["rows"][n] = column_row(ops, row)
            else:
                data["rows"][n] = ladder_row(ops, row, allow128=not a.no_128)
            data["rows"][n]["ladrunoBuild"] = data["ladrunoBuild"]
            path.write_text(json.dumps(data, indent=1, default=float))   # after every row
    data["wall_s"] = time.time() - t0
    path.write_text(json.dumps(data, indent=1, default=float))
    B.flushprint(f"wrote {path} in {data['wall_s'] / 60:.1f} min")
    gated = {n: r.get("verdict") for n, r in data["rows"].items() if r.get("status") == "gate"}
    B.flushprint("gated rows:", gated)


if __name__ == "__main__":
    main()
