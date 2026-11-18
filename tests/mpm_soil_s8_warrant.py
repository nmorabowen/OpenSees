"""WP-232k S8: the 232j warrant figures re-run bit for bit (L1 and L4 r1).

    python3.12 tests/mpm_soil_s8_warrant.py run --out <dir>          # with the build on PYTHONPATH
    python3.12 tests/mpm_soil_s8_warrant.py compare <dirA> <dirB>    # two builds' dumps
    python3.12 tests/mpm_soil_s8_warrant.py record <dir>             # L4 r1 against the 232j record

`run` writes <dir>/l1_affine.json (every L1 affine case of tests/mpm_antilocking_cases.txt through
`mpm debug affineCheck`, the outputs as printed doubles) and <dir>/l4_r1.json (footing_mildcone_r1_al
through the 232j warrant's own runner, A.run_ladder + A.run_record), each with ladrunoBuild. It
never writes under tests/data. `compare` reports equality value by value (== on floats, so bit for
bit) and the largest relative difference. `record` compares l4_r1.json with
tests/data/mpm/antilocking/footing_mildcone_r1_al.json (232j's warrant record, build 4df5764b7 on
esmeralda, GCC): the record is GCC, so it is bit for bit only against a GCC build.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))
REPO = Path(__file__).resolve().parents[1]
RECORD = REPO / "tests" / "data" / "mpm" / "antilocking" / "footing_mildcone_r1_al.json"


def _flat(x, path=""):
    if isinstance(x, dict):
        for k in sorted(x):
            yield from _flat(x[k], f"{path}.{k}")
    elif isinstance(x, (list, tuple)):
        for i, v in enumerate(x):
            yield from _flat(v, f"{path}[{i}]")
    elif isinstance(x, (int, float)) and not isinstance(x, bool):
        yield path, float(x)


def run(out):
    import ladruno_mpm_antilocking_harness as A
    import mpm_antilocking_warrant as Wt
    from _testbed import ops
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    build = A.ladruno_build(ops)
    l1 = {"ladrunoBuild": build, "cases": {}}
    for c in A.parse_cases():
        if c["kind"] == "affine":
            r = A.run_affine(c, ops)
            l1["cases"][c["name"]] = {k: v.get("out") if isinstance(v, dict) else v for k, v in r.items()}
    (out / "l1_affine.json").write_text(json.dumps(l1, indent=1, default=float), encoding="utf-8")
    rec = Wt.run_footing(ops, "footing_mildcone_r1_al")
    (out / "l4_r1.json").write_text(json.dumps(rec, indent=1, default=float), encoding="utf-8")
    print("wrote", out / "l1_affine.json", out / "l4_r1.json", "build", build)


def _compare(a, b, skip=("ladrunoBuild", "wall_s")):
    fa = {p: v for p, v in _flat(a) if not any(s in p for s in skip)}
    fb = {p: v for p, v in _flat(b) if not any(s in p for s in skip)}
    keys = sorted(set(fa) & set(fb))
    diff = [k for k in keys if fa[k] != fb[k]]
    rel = max((abs(fa[k] - fb[k]) / max(abs(fb[k]), 1e-300) for k in diff), default=0.0)
    only = sorted(set(fa) ^ set(fb))
    return {"values": len(keys), "differ": len(diff), "first": diff[:5], "max_rel": rel,
            "only_one_side": only[:5], "bitwise": not diff and not only}


def main(argv=None):
    argv = argv if argv is not None else sys.argv[1:]
    if argv[0] == "run":
        run(argv[argv.index("--out") + 1])
    elif argv[0] == "compare":
        a, b = Path(argv[1]), Path(argv[2])
        for f in ("l1_affine.json", "l4_r1.json"):
            ja, jb = json.loads((a / f).read_text()), json.loads((b / f).read_text())
            print(f, ja.get("ladrunoBuild"), "vs", jb.get("ladrunoBuild"), _compare(ja, jb))
    elif argv[0] == "record":
        ja = json.loads((Path(argv[1]) / "l4_r1.json").read_text())
        jb = json.loads(RECORD.read_text())
        print("l4_r1", ja.get("ladrunoBuild"), "vs record", jb.get("ladrunoBuild"), _compare(ja, jb))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
