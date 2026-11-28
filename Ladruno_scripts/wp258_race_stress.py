"""WP-258 race stress: repeat the soil determinism decks at 8 threads and compare digests.

A data race shows as a digest that differs on SOME runs. One pass of the pytest gates
(tests/test_wp258_soil_threads.py, tests/test_ladruno_mpm_wp258_soil_threads.py) runs each
thread count once; this script runs selected configurations R times at the top thread count
and checks every run against the 1-thread digest. It is the brief's fallback when a TSan build
is not available: GCC's libgomp is not TSan-instrumented, so OpenMP barriers would report as
races without an LLVM / Archer build.

    set PYTHONPATH=<worktree>\\dist\\bin
    python3.12 Ladruno_scripts/wp258_race_stress.py --repeats 10 --threads 8 [--out stress.json]
"""
import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TESTS = ROOT / "tests"
sys.path.insert(0, str(TESTS))
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

import test_wp258_soil_threads as FE                      # noqa: E402
import test_ladruno_mpm_wp258_soil_threads as MPM         # noqa: E402

FE_CASES = [("sanisand_is2", 3), ("sanisand_implex", 3), ("sanisand_is129", 2), ("hs_small", 3),
            ("dp", 2), ("manzari_is1", 3)]
MPM_CASES = ["sanisand_is2_3ds", "sanisand_implex_3d", "hs_implex_3d", "dp_3d", "sanisand_is1_2d"]


def _digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()


def fe_digest(mat, dim, n):
    r, _ = FE._run(mat, dim, n)
    return _digest([r["gravity_rc"], r["steps"]])


def mpm_digest(name, n):
    r = MPM._run(MPM._case(name), n)
    return _digest(r["data"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repeats", type=int, default=10)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    rows = []
    jobs = [("FE", "%s-%dd" % c, lambda n, c=c: fe_digest(c[0], c[1], n)) for c in FE_CASES]
    jobs += [("MPM", c, lambda n, c=c: mpm_digest(c, n)) for c in MPM_CASES]
    bad = 0
    for kind, name, fn in jobs:
        t0 = time.time()
        ref = fn(1)
        got = [fn(a.threads) for _ in range(a.repeats)]
        mism = sum(1 for g in got if g != ref)
        bad += mism
        rows.append(dict(kind=kind, case=name, threads=a.threads, repeats=a.repeats,
                         mismatches=mism, ref=ref[:16], wall_s=round(time.time() - t0, 1)))
        print("%-4s %-22s %d x %d threads: %d mismatch(es)  (%.0f s)"
              % (kind, name, a.repeats, a.threads, mism, time.time() - t0), flush=True)
    if a.out:
        Path(a.out).write_text(json.dumps(rows, indent=1), encoding="utf-8")
    print("TOTAL mismatches:", bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
