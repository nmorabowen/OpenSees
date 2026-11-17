"""WP-232j warrant run: the long rows of tests/mpm_antilocking_cases.txt (L3, L4, L5).

The orchestrator launches it on the 232j build (dist/bin/opensees.pyd on PYTHONPATH through the
tests/_testbed bootstrap) and commits the JSON it writes under tests/data/mpm/antilocking/:

    python3.12 tests/mpm_antilocking_warrant.py --cases footing_mildcone_r1_al footing_mildcone_r2_al
    python3.12 tests/mpm_antilocking_warrant.py --l4                  # footing_mildcone_r1..r3_al
    python3.12 tests/mpm_antilocking_warrant.py --l3                  # footing_j2_r1 only (freeze F2)
    python3.12 tests/mpm_antilocking_warrant.py --l3 --admitted       # + r2, r3 once L3 is admitted
    python3.12 tests/mpm_antilocking_warrant.py --l5                  # MMS G2, G4 with anti-locking on
    python3.12 tests/mpm_antilocking_warrant.py --off <case>          # the same case with it off (R rows)

Each footing record: ladrunoBuild, the setting, q(s) per step, every observable with its verdict,
failed attempts, relaxed steps, halvings (K6), averaged / fallback particle-steps and the fallback
share in the plastic zone x <= 3B, y >= -2B (freeze F1), wall time. A miss is recorded as a miss;
no band is widened here or anywhere. F2: L3 at r2 / r3 refuses to run without --admitted.

Estimated runtime (the nmb box, from 232f's footing runs): r1 a few minutes, r2 about 20 min,
r3 1-3 h each; L5 about 1 min.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_antilocking_harness as A  # noqa: E402

OUT = Path(__file__).resolve().parent / "data" / "mpm" / "antilocking"
L4 = ["footing_mildcone_r1_al", "footing_mildcone_r2_al", "footing_mildcone_r3_al"]
L3 = ["footing_j2_r1", "footing_j2_r2", "footing_j2_r3"]


def _ops():
    from _testbed import ops
    return ops


def run_footing(ops, name, off=False):
    c = A.case_by_name(name)
    if off:
        c["antilocking"] = "off"
    r = A.run_ladder(c, ops, log=A.flushprint)
    rec = A.run_record(c, r)
    return rec


def run_mms(ops, name):
    import ladruno_mpm_bspline_harness as B
    c = A.case_by_name(name)
    row = B.parse_rows()[c["row"]]
    orig = B.build_mms

    def build_on(o, rw, ncell):
        X = orig(o, rw, ncell)
        o.mpm("antiLocking", c["antilocking"])
        return X

    B.build_mms = build_on
    try:
        res = B.ladder(ops, row, log=A.flushprint)
    finally:
        B.build_mms = orig
    ok, why = B.verdict(row, res)
    res.update({"case": name, "row": row["name"], "antilocking": c["antilocking"],
                "ladrunoBuild": A.ladruno_build(ops), "verdict": ok, "verdict_text": why})
    return res


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="*", default=[])
    ap.add_argument("--l3", action="store_true")
    ap.add_argument("--l4", action="store_true")
    ap.add_argument("--l5", action="store_true")
    ap.add_argument("--admitted", action="store_true",
                    help="L3 is admitted (freeze F2): run footing_j2_r2 / r3 too")
    ap.add_argument("--off", nargs="*", default=[], help="footing cases to run with anti-locking off")
    a = ap.parse_args(argv)
    names = list(a.cases)
    if a.l4:
        names += L4
    if a.l3:
        names += L3 if a.admitted else L3[:1]
    for n in names:
        if n in L3[1:] and not a.admitted:
            sys.exit(f"{n}: L3 at r2 / r3 waits for L3's admission (freeze F2); pass --admitted")
    ops = _ops()
    OUT.mkdir(parents=True, exist_ok=True)
    for n in names:
        rec = run_footing(ops, n)
        (OUT / f"{n}.json").write_text(json.dumps(rec, indent=1), encoding="utf-8")
    for n in a.off:
        rec = run_footing(ops, n, off=True)
        (OUT / f"{n}_off.json").write_text(json.dumps(rec, indent=1), encoding="utf-8")
    if a.l5:
        for n in ("mms_G2_al", "mms_G4_al"):
            res = run_mms(ops, n)
            (OUT / f"{n}.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")


if __name__ == "__main__":
    main()
